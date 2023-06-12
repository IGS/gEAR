# NOTES

* For requester pays:
  * archive-file-transfer@nemo-analytics.iam.gserviceaccount.com does not have serviceusage.services.use access to the Google Cloud project
* To ensure Google has credentials to work with:
  * For the service account to use, download JSON keys
  * Put keys in VM
  * Add key filepath to gear.ini->nemoarchive_import->credentials_json