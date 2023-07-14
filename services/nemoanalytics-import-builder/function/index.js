'use strict';

const groupBy = require('lodash.groupby');
const map = require('lodash.map');
const pluralize = require('pluralize');
const {Storage} = require('@google-cloud/storage');
const uniq = require('lodash.uniq');
const uuid = require('uuid/v4');
const {validate} = require('jsonschema');

const storage = new Storage();
const {FUNCTION_NAME, FUNCTION_REGION, GCP_PROJECT, gcp_bucket, secret_token} = process.env;

/**
 * Bucket which should host exported JSON payloads
 * as a midpoint for NeMO Analytics' pull-based import API.
 */
const stagingBucket = storage.bucket(gcp_bucket);

/** Options to use on every object-write to the staging bucket. */
const stagingOptions = { resumable: false, contentType: 'application/json' };

/** JSON schema for this function. */
const schema = {
    "$schema": "http://json-schema.org/draft-04/schema",
    "id": `https://${FUNCTION_REGION}-${GCP_PROJECT}.cloudfunctions.net/${FUNCTION_NAME}`,
    "definitions": {
        "simpleAttribute": {
            "type": ["string", "number", "boolean", "null"]
        },
        "repeatedAttribute": {
            "type": "array",
            "items": { "$ref": "#/definitions/simpleAttribute" },
            "default": []
        },
        "entityAttribute": {
            "anyOf": [
                { "$ref": "#/definitions/simpleAttribute" },
                { "$ref": "#/definitions/repeatedAttribute" }
            ]
        },
        "entityId": {
            "type": "string",
            // NOTE: This pattern was pulled from the internals of Rawls,
            // it might change without notice in the future.
            "pattern": "^[A-Za-z0-9\\-_]+$"
        },
        "entity": {
            "type": "object",
            "properties": {
                "name": { "$ref": "#/definitions/entityId" },
                "entityType": { "type": "string" },
                "attributes": {
                    "type": "object",
                    "additionalProperties": { "$ref": "#/definitions/entityAttribute" }
                }
            },
            "required": ["name", "entityType"]
        },
        "secretToken": {
            "type": "string"
        },
        "entities": {
            "type": "array",
            "items": { "$ref": "#/definitions/entity" }
        }
    },
    "type": "object",
    "properties": {
        "entities": { "$ref": "#/definitions/entities" },
        "cohortName": { "$ref": "#/definitions/entityId" },
        "secretToken": { "ref": "#/definitions/secretToken"}
    },
    "required": ["entities", "cohortName", "secretToken"]
};

/**
 * Extract and return the unique name/type pairs from a collection
 * of entities.
 *
 * @param {Object[]} entities - JSON entities to group into sets.
 * @param {string} entities[].name - Unique ID of an entity.
 * @param {string} entities[].entityType - Type of an entity.
 */
const collectUniqueIds = (entities) => {
    return uniq(map(entities, ({entityType, name}) => {
        return {entityType, name}
    }));
}

/**
 * Build a collection of "set" entities for NeMO Analytics encapsulating
 * all the entities in an export, with one set per entity type.
 *
 * @param {Object[]} entities - JSON entities to group into sets.
 * @param {string} entities[].name - Unique ID of an entity.
 * @param {string} entities[].entityType - Type of an entity.
 * @param {string} cohortName - Name to use for this group of entities.
 */
const collectEntitySets = (entities, cohortName) => {
    const groups = groupBy(entities, 'entityType');
    return map(groups, (entityGroup, entityType) => {
        const items = map(entityGroup, ({name}) => {
            return { entityType, entityName: name };
        });

        return {
            name: cohortName,
            entityType: `${entityType}_set`,
            // NOTE: The schema for these attributes was reverse-engineered
            // from reading the Rawls code + trial-and-error; I wouldn't be
            // surprised if it changed out from under us.
            attributes: {
                [pluralize(entityType)]: {
                    itemsType: 'EntityReference',
                    items
                }
            }
        };
    });
}

/**
 * Write entity JSON to a staging location in GCS, returning a Promise
 * which (on success) will complete with a signed URL pointing to the
 * written object.
 *
 * @param {Object[]} entities - JSON entities to write into GCS.
 */
const writeStagingFile = (entities) => {
    const stagingObject = stagingBucket.file(`${uuid()}.json`);
    return stagingObject.save(JSON.stringify(entities), stagingOptions)
        .then(() => {
            // Build a signed URL valid for 1 hour.
            const urlExpirationMs = Date.now() + 3600000;
            const signingConfig = {
                action: 'read',
                expires: urlExpirationMs
            };
            return stagingObject.getSignedUrl(signingConfig)
                .then(([url]) => url);
        });
}

/**
 * HTTP entry point; converts a JSON push from a data explorer
 * into a pullable URL for NeMO Analytics.
 *
 * @param {Object} req - https://expressjs.com/en/api.html#req
 * @param {Object} res - https://expressjs.com/en/api.html#res
 */
exports.buildNemoanalyticsImport = (req, res) => {
    if (req.method === 'GET') {
        res.json(schema);
    } else if (req.method !== 'POST') {
        res.status(405).send('');
    } else if (req.get('content-type') !== 'application/json') {
        res.status(415).send('');
    } else if (!validate(req.body, schema).valid) {
        res.status(422).json({
            message: 'Request body does not match expected schema',
            schema
        })
    } else {
        const { entities, cohortName, secretToken } = req.body;
        if (secretToken !== secret_token) [
            res.status(401).json({ message: "Unauthorized."})
        ]
        const uniqueIds = collectUniqueIds(entities);

        if (uniqueIds.length !== entities.length) {
            // Fail fast here
            res.status(400).json({ message: 'Duplicate entity IDs found in payload' });
        } else {
            const entitySets = collectEntitySets(uniqueIds, cohortName);
            writeStagingFile(entities.concat(entitySets))
                .then((url) => res.status(201).json({ url }))
                .catch((err) => {
                    console.error(err);
                    res.status(500).json({ message: 'Failed to build NeMO Analytics import bundle' });
                });
        }
    }
};
