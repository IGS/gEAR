# NOTES

* For requester pays:
  * `archive-file-transfer@nemo-analytics.iam.gserviceaccount.com` does not have serviceusage.services.use access to the Google Cloud project
* To ensure Google has credentials to work with:
  * For the service account to use, download JSON keys
  * Put keys in VM in <gEAR_root>/www
  * Add key filepath to gear.ini->nemoarchive_import->credentials_json

## Setting up the Cloud Function service

Please see the file at file://../../services/nemoanalytics-import-builder/README.md

## API ROUTES

An OpenAPI document was created at file://../api/openapi.yml to view a description of the available routes, their inputs and outputs.

TODO: Including a server-based link instead of the YAML doc

## FLOW

The nemoarchive import page has a basic layout. At the top of the page there is some descriptive text as well as the submission ID. Beside the submission ID is a share button that, when clicked, will copy a "view-only" URL to the clipboard.

Below this is a table of datasets in the submission. For each dataset, the following columns are present:

* nemo: identifier name (with permalink to nemoarchive assets page for this identifier)
* Pulled to VM status
* Validate Metadata status
* Write to H5AD status
* Make tSNE status
* Messages, for displaying any error messages that occur in the import process
* Select metadata group (if observation metadata exists) for creating a new tSNE display with added group annotation
  * This only becomes available once the "Write H5AD" step is complete
* Permalink to view dataset only in gene expression display page.
  * This only becomes available once the "Write H5AD" step is complete

Below this is a few fields that become visible once a single dataset import is complete. This includes:

* Input to rename the submission layout label (by default named "Submission <submission_id>")
  * If the user viewing this submission page is not the same user who created this submission, this field will be disabled.
* Input to select a gene to search for when performing a gene expression search on this submission "profile"
* Button to redirect to the gene expression search on this submission "profile"
  * It is worth noting that if a "metadata group" is selected in the table or "gene" input is provided, this will quickly create a new display config with those inserted, and save that as a new default display for the user.
* Button to subscribe to email updates. The "send_updates" field for a submisison in the database will be updated to 1, and the submission user will be notified of an email update when the submission finishes (successful or not).

If a submission is currently being imported, the UI will poll for status updates from the server for each dataset and update their statuses in the table accordingly. Polling will keep happening every 10 seconds until every dataset is not in "pending" or "loading" state.

### Route 1 - View-only mode

`/nemoarchive_import/import.html?submission_id=<submission_id>`

This is intended for a user to check in on an existing submission and see how it is going. This route will not trigger the submission to start.

### Route 2 - Import mode

`/nemoarchive_import/import.html?url=<url>`

This route is intended to kick off the import of a new submmission. This is the general flow:

1. **Fetch JSON from "url" parameter.** This will include the submission ID that we can create a submission permalink from.
2. **Load submission.**
    1. Get the submission from the API (`GET /api/submissions/<submission_id>`) if it is already present. This is not the intended way to use the "import" mode, but I can see situations where a user may refresh or come back to this particular URL
        1. This GET command is expected to abort with a 404 status if the submission does not exist and that is OK
    2. If the submission does not already exist, then create a new one (`POST /api/submissions`)
        1. We also save the submission as a new layout in the database with a generic label. The layout_id is saved as a DOM data attribute for redirecting to the gene expression profile view.
3. **Load datasets for submission.** Find all dataset IDs in the "file" entity parts of the JSON. These are UUIDs (though this may change in the future, so we can rely on identifiers as a backup if need be).
    1. Datasets can be associated with multiple submissions. For this reason, we attempt to find the dataset already (`GET /api/submissions/<submission_id>/datasets/<dataset_id>`). The "submission_id" is merely a route placeholder so the user cannot go to a submission_dataset directly.
    2. If the dataset does not exist at all, we create a new submission_dataset (`POST /api/submissions/<submission_id>/datasets`) in the database.
        1. The submission_dataset ID will also create a minimal dataset table entry in the database with "pending" status. This dataset table entry will be filled out more as the import goes along.
        2. The NeMO identifier will become the dataset share_id.
4. **Associate dataset with submission.** This happens regardless if we retrieved an existing submission or created a new one. A submission_member database entry is created (`PUT /api/submissions/<submission_id>/datasets/<dataset_id>`) which associates the submission to the dataset.
    1. At this point, the dataset is added as a layout_member to the submission layout as well.
5. **Create submission table rows per dataset.**
    1. If any datasets could not be added to the submission, an error message to indicate as such will show in the "messages" column.
6. **Start the import process for each dataset.** This is triggered by (`POST /api/submissions/<submission_id>` with action="import").
7. **(optional) Receive email update when finished.**

## IMPORT PROCESS

The import process can either run inside the Apache process or a message can be sent to RabbitMQ to trigger the process (preferred).  Each of the dataset API calls will be run using async/await (python asyncio and aiohttp modules) so that the CPU is not inactive while IO stuff is happening.

### Pull GCP files to VM

This step identifies each of the component files that make up a particular "file" dataset. The Google Cloud Platform (GCP) bucket path is retrieved and each file is pulled from GCP to the `<gEAR_root>/www/uploads/files` area under a directory named after the dataset UUID.

### Validate Metadata

This step builds a JSON structure out of the file and sample metadata. There is some minor sanitation that happens. This step also parses the genes in the dataset and uses the organism ID to determine the best Ensembl release for this dataset and add to metadata (since there is zero chance this is in the metadata). The JSON structure is written out to the area the pulled dataset files are at, and is converted to a Pandas dataframe which is passed to a gear.Metadata object and validated. Upon validation, the metadata is written to the dataset's "dataset" table entry.

**NOTE**: Currently the NeMO Archive API to retrive metadata is not developed yet. So I pulled a very early build of the new nemo_assets database, set up some queries, made some mock data, and set that to be retrieved from (`GET /api/mock_identifier/<string:identifier>`). Basically when I export from the NeMO Archive portal, I grab the JSON from the staged URL, and modify it to have identifiers from this `nemo_assets.file.nemo_id` field, and save the JSON locally. The JS is also hardcoded currently to fetch from this hardcoded JSON (for demonstration purposes), and the `POST /api/submissions/<submission_id>` with action="import" route fetches metadata from the mock API route instead of the assets.nemoarchive.org API call.  Currently some metadata is used from the imported JSON (and hard-coded) and others is from the API call, but I anticipate that eventually I can just export the identifier and the files and all the sample metadata can come from the API call.

### Write to H5AD

This step looks at the uploaded dataset files (plus JSON file) and uses them in a gear.DataArchive object. This object will now convert the dataset files to H5AD inside `<gEAR_root>/www/datasets`. If the filetype is already H5AD (no conversion) or is not deemed to have Ensembl IDs present for the adata.var index, those are added using the gear database. The JSON is also copied to `<gEAR_root>/www/datasets` and the database is updated

### Make tSNE

This step will run through the basic Seurat clustering tutorial (https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html) with default arguments since we do not know anything about the data. This runs through the tSNE calculation steps to get a modified AnnData object. This is saved as a public analysis, and a mock "tsne_static" plot_config JSON is created. This plot config is saved as a dataset display and set as the user's default. If the user chooses a sample metadata category or a gene in the UI later on, those are incorporated into the plot config and a new display and default are saved (this can be applicable for a third-party who is viewing this submission using the "view-only" mode too.)

### Post-import

After the results from each dataset come in, we determine if the import for that dataset succeeded or not. If the user requests an email notification after finishing we send one. The contents will change depending if all, some, or none of the datasets were successfully imported. Any dataset import failures should also log a message into the submission_dataset table in the database, which can be retrieved when polling the dataset status.

## Future things

* Right now this only works with file-based datasets.  We will have to modify the processes some when we have sample-based datasets, which could consist of multiple files
* Also, currently we treat all datasets as unrestricted for now, since public data is the only kind of data we can access for the import process. There are some variables to store datasets and submissions as "restricted" but we currently have no way of handling them (i.e. authentication)
