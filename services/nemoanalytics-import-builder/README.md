# NeMO Analytics Import Builder

Cloud Function for staging JSON payloads in GCS for use by NeMO Analytics' pull-based import API.

Code is heavily based off of https://github.com/BICCN/terra-import-builder

I elected to build the cloud function in the online editor, and specified my GCP bucket as a service environmental variable when creating the function. In addition, I created separate static dev and prod buckets as well as functions with different dev and prod CORS policies, rather than have the function create a temporary bucket during the deployment step.

Both dev and prod cloud function services also have a "secret_token" environment variable set that ensures only NeMO Archive can post to the function.

The service account used to execute the function must have "SignBlob" granted

## Adding cors to a bucket

```bash
gsutil -u <project_id> cors set cors-config.json gs://<bucket_name>
```

cors-config.json will be the CORS policy to apply to the bucket. See example CORS json files at https://cloud.google.com/storage/docs/cors-configurations#gsutil

## Adding a new environment

To set up a new deployment env for the import-builder, do the following in the GCP Cloud Console:

1. Create a new Cloud Run project project in GCP to host the env
2. Enable APIs in the new project:
   * IAM
   * Cloud Function
3. Create a Service Account to run the import-builder
   * As mentioned above, I elected to build the cloud function in the online editor. The source code is located in the Cloud Functions service -> <project> -> Source tab
   * The deploy.sh script is used to
   * Give it the "Service Account Token Creator" IAM role
   * Generate and download a JSON key for it. This service key will have the credential information and be stored locally to authenticate into Google Cloud services.
   * I also give the
   * NOTE: This service will also appear in the list of Cloud Functions for that console.
4. Create a GCS bucket in the project to host import bundles
   * Grant "Storage Object Creator" permissions on the bucket to the account created in 3
   * Set a lifecycle policy on the bucket to delete files after 1 day

The `deploy.sh` script is not really needed since we build the service in the online editor, but I am keeping around since it was present when I forked from the terra-import-builder. It could be useful in the future.