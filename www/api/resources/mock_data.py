from flask import request
from flask_restful import Resource
from pathlib import Path
import os, sys

import mysql.connector
from mysql.connector import errorcode

# Parse gEAR config
# https://stackoverflow.com/a/35904211/1368079
this = sys.modules[__name__]
from gear.serverconfig import ServerConfig
this.servercfg = ServerConfig().parse()

def connect_to_assets_db():
    try:
        cnx = mysql.connector.connect(user=this.servercfg['database']['user'], password=this.servercfg['database']['password'],
                                        host=this.servercfg['database']['host'], database="nemo_assets",
                                        buffered=True, use_pure=True)
        return cnx
    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            print("Something is wrong with your user name or password", file=sys.stderr)
        elif err.errno == errorcode.ER_BAD_DB_ERROR:
            print("Database does not exist", file=sys.stderr)
        else:
            print(err, file=sys.stderr)
        raise

class MockIdentifier(Resource):
    """
    Monkeypatched code to mock identifier metadata
    while the NeMO Archive API is being developed
    """
    def get(self, identifier):

        # using "file.nemo_id" from nemo_assets db. Only 20 FASTQ files are associated with collections,
        # so I mocked using those IDs to ensure I can grab db data
        (before, sep, nemo_id) = identifier.rpartition("der-")

        response = {"self":request.path, "success":0, "message":"", "metadata":{}}

        try:
            connection = connect_to_assets_db()
        except Exception as e:
            response["message"] = str(e)
            return response


        """ DATASET LEVEL
        * Identifier    // collection.col_id or (single) file.nemo_id (without the )
        * contact_name  // contributor.name
        * contact_email // contributor.email
        * contact_institute // contributor.organization
        * title // collection.name
        * summary   // collection.description
        * dataset_type  // file.data_type_id -> data_type.data_type
        * reference_annot_id    // Will have to assume it does not exist
        repository_accession
        * sample_taxid  // taxonony.name -> subject_taxonomy.taxonomy_id -> subject.id -> subject_assoc_sample.sample_id -> sample.id
        pubmed_id   // publication.pubmed_id -> collection.id -> file_in_collection.file_id (optional?)
        keywords    // keywords.keyword -> collection_has_keywords.collection_id ->
        device
        library_cons_method
        RNA_extraction_method
        microarray_hybridization_protocol
        * assay // assay.name
        * normalization_method  // not provided - check h5ad properties
        * log_transformation    // not provided - check h5ad properties
        * primary_analysis_completed    // not provided - check h5ad structure.
        (cell_type)
        * (tissue_type) // sample.specimen_type_id -> specimen_type.name (? unsure)
        * (file_format) // file.file_format_id -> file_format.name
        """

        """ SAMPLE LEVEL
        * sample_id // sample.id
        tissue_label
        tissue_ontology_name
        * tissue_ontology   // anatomy.name -> sample_assoc_anatomy.sample_id -> sample.id
        cell_label
        cell_ontology_name
        cell_ontology
        * treatment // subject_observations.value + subject_observations.obs_var_id -> obs_vars.name = "medication"
        treatment_label
        DRmethod_#
        DRmethod_#
        cell_clusters
        * age_value // subject_observations.value + subject_observations.obs_var_id -> obs_vars.name = "age"
        * age_unit  // subject_observations.unit + subject_observations.obs_var_id -> obs_vars.name = "age"
        age_order?
        * sex_assigned_at_birth // subject_observations.value + subject_observations.obs_var_id -> obs_vars.name= "sex" (assumed at birth)
        """

        # NOTE: Only raw files are directly associated with collections. Easiest way to link a counts file to a collection
        # is to join to the sample and then join to the collection.

        dataset_metadata = {}
        #dataset_metadata["identifier"] = get_identifier(connection, nemo_id) # should be already got
        dataset_metadata["contact_name"] = get_contact_name(connection, nemo_id)
        dataset_metadata["contact_email"] = get_contact_email(connection, nemo_id)
        dataset_metadata["contact_institute"] = get_contact_institute(connection, nemo_id)
        dataset_metadata["title"] = get_title(connection, nemo_id)
        dataset_metadata["summary"] = get_summary(connection, nemo_id)
        #dataset_metadata["dataset_type"] = get_dataset_type(connection, nemo_id)    # got in neo4j
        dataset_metadata["reference_annot_id"] = None
        #dataset_metadata["sample_taxid"] = get_taxid(connection, nemo_id)   # got in neo4j (tax name)
        #dataset_metadata["assay"] = get_assay(connection, nemo_id)
        dataset_metadata["normalization_method"] = None
        dataset_metadata["log_transformation"] = "raw"
        dataset_metadata["primary_analysis_completed"] = False
        dataset_metadata["tissue_type"] = get_tissue_type(connection, nemo_id)
        #dataset_metadata["file_format"] = get_file_format(connection, nemo_id)   # should be already got

        sample_metadata = {}
        #sample_metadata["sample_id"] = get_sample_id(connection, nemo_id)    # should be already got
        sample_metadata["tissue_ontology"] = get_tissue_ontology(connection, nemo_id)
        sample_metadata["treatment"] = get_treatment(connection, nemo_id)
        sample_metadata["age_value"] = get_age_value(connection, nemo_id)
        sample_metadata["age_unit"] = get_age_unit(connection, nemo_id)
        sample_metadata["sex_assigned_at_birth"] = get_sex_assigned_at_birth(connection, nemo_id)

        response["metadata"]["dataset"] = dataset_metadata
        response["metadata"]["sample"] = sample_metadata

        response["success"] = 1

        connection.close()

        return response

def get_contact_name(conn, nemo_id):
    cursor = conn.cursor()
    query = """
        SELECT con.name from contributor con
        JOIN collection_has_contributor chc on con.id = chc.contributor_id
        JOIN file_in_collection fic on chc.collection_id = fic.file_id
        JOIN file f on fic.file_id = f.id
        WHERE f.nemo_id = %s
        """
    #cursor.execute(query, (nemo_id, ))

    # collection_has_contributor is not populated.... retrieve random contributor
    query = "SELECT name from contributor where id = 44"
    cursor.execute(query)

    try:
        (result,) = cursor.fetchone()
        return result
    except:
        return None
    finally:
        cursor.close()

def get_contact_email(conn, nemo_id):
    cursor = conn.cursor()
    query = """
        SELECT con.email from contributor con
        JOIN collection_has_contributor chc on con.id = chc.contributor_id
        JOIN file_in_collection fic on chc.collection_id = fic.collection_id
        JOIN file f on fic.file_id = f.id
        WHERE f.nemo_id = %s
        """
    #cursor.execute(query, (nemo_id, ))

    # collection_has_contributor is not populated.... retrieve random contributor
    query = "SELECT email from contributor where id = 44"
    cursor.execute(query)

    try:
        (result,) = cursor.fetchone()
        return result
    except:
        return None
    finally:
        cursor.close()

def get_contact_institute(conn, nemo_id):
    cursor = conn.cursor()
    query = """
        SELECT con.organization from contributor con
        JOIN collection_has_contributor chc on con.id = chc.contributor_id
        JOIN file_in_collection fic on chc.collection_id = fic.collection_id
        JOIN file f on fic.file_id = f.id
        WHERE f.nemo_id = %s
        """
    #cursor.execute(query, (nemo_id, ))

    # collection_has_contributor is not populated.... retrieve random contributor
    query = "SELECT organization from contributor where id = 44"
    cursor.execute(query)

    try:
        (result,) = cursor.fetchone()
        return result
    except:
        return None
    finally:
        cursor.close()

def get_title(conn, nemo_id):
    cursor = conn.cursor()
    query = """
        SELECT c.name from collection c
        JOIN file_in_collection fic on c.id = fic.collection_id
        JOIN file f on fic.file_id = f.id
        WHERE f.nemo_id = %s
        """
    cursor.execute(query, (nemo_id, ))

    try:
        (result,) = cursor.fetchone()
        return result
    except:
        return "Placeholder title for {}".format(nemo_id)
    finally:
        cursor.close()


def get_summary(conn, nemo_id):
    cursor = conn.cursor()
    query = """
        SELECT c.description from collection c
        JOIN file_in_collection fic on c.id = fic.collection_id
        JOIN file f on fic.file_id = f.id
        WHERE f.nemo_id = %s
        """
    cursor.execute(query, (nemo_id, ))

    try:
        (result,) = cursor.fetchone()
        return result
    except:
        return "Description not available in the mock database. This is a placeholder"
    finally:
        cursor.close()

def get_dataset_type(conn, nemo_id):
    pass

def get_taxid(conn, nemo_id):
    pass

def get_assay(conn, nemo_id):
    cursor = conn.cursor()
    query = """
        SELECT a.name from assay a
        JOIN sample s on a.id = s.assay_id
        JOIN sample_assoc_file saf on s.id = saf.sample_id
        JOIN file f on saf.file_id = f.id
        WHERE f.nemo_id = %s
        """
    cursor.execute(query, (nemo_id, ))

    try:
        (result,) = cursor.fetchone()
        return result
    except:
        return None
    finally:
        cursor.close()

def get_tissue_type(conn, nemo_id):
    cursor = conn.cursor()
    query = """
        SELECT st.name from specimen_type st
        JOIN sample s on st.id = s.specimen_type_id
        JOIN sample_assoc_file saf on s.id = saf.sample_id
        JOIN file f on saf.file_id = f.id
        WHERE f.nemo_id = %s
        """
    cursor.execute(query, (nemo_id, ))

    try:
        (result,) = cursor.fetchone()
        return result
    except:
        return None
    finally:
        cursor.close()

def get_sample_id(conn, nemo_id):
    pass

def get_tissue_ontology(conn, nemo_id):
    cursor = conn.cursor()
    query = """
        SELECT a.name from anatomy a
        JOIN sample_assoc_anatomy saa on a.id = saa.anatomy_id
        JOIN sample s on saa.sample_id = s.id
        JOIN sample_assoc_file saf on s.id = saf.sample_id
        JOIN file f on saf.file_id = f.id
        WHERE f.nemo_id = %s
        """
    cursor.execute(query, (nemo_id, ))

    try:
        (result,) = cursor.fetchone()
        return result
    except:
        return None
    finally:
        cursor.close()

def get_treatment(conn, nemo_id):
    cursor = conn.cursor()
    query = """
        SELECT sbo.value from subject_observations sbo
        JOIN obs_vars ov on sbo.obs_vars_id = ov.id
        JOIN sample_assoc_subject sasb on sasb.subject_id = sbo.subject_id
        JOIN sample_assoc_file saf on saf.sample_id = sasb.sample_id
        JOIN file f on saf.file_id = f.id
        WHERE ov.var_name = "medication"
        AND f.nemo_id = %s
        """
    cursor.execute(query, (nemo_id, ))

    try:
        (result,) = cursor.fetchone()
        return result
    except:
        return None
    finally:
        cursor.close()

def get_age_value(conn, nemo_id):
    cursor = conn.cursor()
    query = """
        SELECT sbo.value from subject_observations sbo
        JOIN obs_vars ov on sbo.obs_vars_id = ov.id
        JOIN sample_assoc_subject sasb on sasb.subject_id = sbo.subject_id
        JOIN sample_assoc_file saf on saf.sample_id = sasb.sample_id
        JOIN file f on saf.file_id = f.id
        WHERE ov.var_name = "age"
        AND f.nemo_id = %s
        """
    cursor.execute(query, (nemo_id, ))

    try:
        (result,) = cursor.fetchone()
        return result
    except:
        return None
    finally:
        cursor.close()

def get_age_unit(conn, nemo_id):
    cursor = conn.cursor()
    query = """
        SELECT sbo.unit from subject_observations sbo
        JOIN obs_vars ov on sbo.obs_vars_id = ov.id
        JOIN sample_assoc_subject sasb on sasb.subject_id = sbo.subject_id
        JOIN sample_assoc_file saf on saf.sample_id = sasb.sample_id
        JOIN file f on saf.file_id = f.id
        WHERE ov.var_name = "age"
        AND f.nemo_id = %s
        """
    cursor.execute(query, (nemo_id, ))

    try:
        (result,) = cursor.fetchone()
        return result
    except:
        return None
    finally:
        cursor.close()

def get_sex_assigned_at_birth(conn, nemo_id):
    cursor = conn.cursor()
    query = """
        SELECT sbo.value from subject_observations sbo
        JOIN obs_vars ov on sbo.obs_vars_id = ov.id
        JOIN sample_assoc_subject sasb on sasb.subject_id = sbo.subject_id
        JOIN sample_assoc_file saf on saf.sample_id = sasb.sample_id
        JOIN file f on saf.file_id = f.id
        WHERE ov.var_name = "sex"
        AND f.nemo_id = %s
        """
    cursor.execute(query, (nemo_id, ))

    try:
        (result,) = cursor.fetchone()
        return result
    except:
        return None
    finally:
        cursor.close()