# Retrieving Dataset Information from the command-line

## Terminology

* datasetId - The long-form UUID of a dataset.  It is the primary-key
* shareId- The shorter UID for a dataset, meant for sharing purposes, and for permalinked short-hand URLs

## Retrieving dataset metadata

### Retrieving all public datasets

This command will only retrieve dataset information from public datasets.

`curl "https://cancer.umgear.org/cgi/get_h5ad_dataset_list.cgi"`

The returned response is a JSON object with the following keys:

* public
* user
* shared_with_user

Each key contains an array of datasets. If you are not logged in, only the "public" array will be populated.  Each item in the array will be one dataset's metadata.  Of note is the "id" which is the defined datasetId above.

### Retrieving all datasets as a logged-in user

If you know your sessionId from having logged in, you can pass that in to also retrieve all datasets you have access to.

`curl "https://cancer.umgear.org/cgi/get_h5ad_dataset_list.cgi?session_id=${sessionId}"`

The "user" and "shared_with_user" arrays may now have datasets, depending on if you own any datasets or have access to other datasets.

### Retrieving a specific dataset using the dataset share ID

If you know the shareId for the dataset you want to retrieve, you can pass that as a query parameter
`curl "https://cancer.umgear.org/cgi/get_h5ad_dataset_list.cgi?share_id=${shareId}"`

This dataset will appear in the "shared_with_user" array in the JSON response

### Retrieving a specific dataset using the dataset primary ID

If you have the datasetId, you can retrieve that dataset information directly.

`curl "https://cancer.umgear.org/cgi/get_dataset_info.cgi?dataset_id=<${datasetId}"`\

The response is a single JSON object with key-value pairs for that dataset's metadata

appending `&include_shape=1` to the URL will also retrieve the number of genes and observations in the dataset.  Do note this will add time to the request, as the dataset has to be loaded briefly.


## Downloading dataset information

For these commands, changing the value of the "type" query parameter retrieves different information.

To download the H5AD file:

`curl -v "https://cancer.umgear.org/cgi/download_source_file.cgi?type=h5ad&share_id=${shareId}"`

To download dataset metadata:

`curl -v "https://cancer.umgear.org/cgi/download_source_file.cgi?type=metadata&share_id=${shareId}"`

