# NeMO Analytics Import Builder

Cloud Function for staging JSON payloads in GCS for use by NeMO Analytics' pull-based import API.

Code is heavily based off of https://github.com/BICCN/terra-import-builder

I elected to build the cloud function in the online editor, and specified my GCP bucket as an environmental variable when creating the function. In addition, I created separate dev and prod buckets and functions with different CORS policies, rather than have the function create a temporary bucket during the deployment step

The service account used to execute the function must have "SignBlob" granted

## Adding cors to a bucket

```bash
gsutil -u <project_id> cors set cors-config.json gs://<bucket_name>
```

cors-config.json will be the CORS policy to apply to the bucket. See example CORS json files at https://cloud.google.com/storage/docs/cors-configurations#gsutil

## Adding a new environment (for terra, not nemo import)

To set up a new deployment env for the import-builder, do the following in the GCP Cloud Console:

1. Create a new project in GCP to host the env
2. Enable APIs in the new project:
   * IAM
   * Cloud Function
3. Create a Service Account to run the import-builder
   * Give it the "Service Account Token Creator" IAM role
   * Generate and download a JSON key for it
4. Create a GCS bucket in the project to host import bundles
   * Grant "Storage Object Creator" permissions on the bucket to the account created in 3
   * Set a lifecycle policy on the bucket to delete files after 1 day

Then, from your terminal do:

```bash
$ vault write secret/dsde/monster/${env}/biccn/nemoanalytics-import-builder/env \
    project=${the-new-project} \
    bucket=${the-new-bucket}

$ vault write secret/dsde/monster/${env}/biccn/nemoanalytics-import-builer/service-account.json \
    @${path-to-the-downloaded-json-keyfile}
```

Finally, add the env name to the `VALID_ENVS` constant declared at the top of `deploy.sh`.
