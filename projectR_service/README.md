# ProjectR Container service

This directory is meant to house the code to run projectR inside of a container. This container can then be deployed to a Google Cloud service, such as Google Cloud Run

Details - https://cloud.google.com/run/docs/quickstarts/build-and-deploy/deploy-python-service

## Deploying via the command line

If the service is not deployed to a Google cloud project, cd to this directory (projectR_service) and run `gcloud run deploy`.

Use the defaults for the source code location and for the service name.

For the region, I used us-east4, but I do not think it will matter too much, as long as it is one of the "us-east" regions.

If you get a message that the Artifact Registry is not enabled, enable that and follow the prompts.

Same goes for being prompted to allow unauthenticated invocations

Next up, the Cloud Build service starts building the Docker container image.  Given that Bioconductor needs to be installed from the Dockerfile, this may take a while.