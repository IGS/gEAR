# ProjectR Container service

This directory is meant to house the code to run projectR inside of a container. This container can then be deployed to a Google Cloud service, such as Google Cloud Run

Details - https://cloud.google.com/run/docs/quickstarts/build-and-deploy/deploy-python-service

NOTE: Cloud Build times out when attempting to build this image. Need to build locally. For my purposes, I named my image "projectr_service".

Also if building on an M1 Mac, you need to add platform information to make it build with Linux architecture.  Otherwise it will not be able to run in the Google Cloud services.  My docker command to build was `docker build --platform linux/amd64 --no-cache -t projectr_service .`

## Pushing Docker images to Google Artifact Registry

Google Artifact Registry has recently replaced Google Container Registry for the de facto location to store Docker container images along with other artifact types. When selecting the Artifact Registry from the Google Cloud web console, you will see two tabs: Artifact Registry and Container Registry. Container Registry is just a holdover from the older service, so we will focus on Artifact Registry

We need to have our Docker configuration file modified to support gcloud credentials when pushing images. Since we set our Cloud Run region to us-east4, we need to run `gcloud auth configure-docker us-east4-docker.pkg.dev`, which modifies `~/.docker/config.json`

To push our created image, we need to give it a specific tag after building. Follow the instructions at https://cloud.google.com/artifact-registry/docs/docker/pushing-and-pulling#pushing for pushing the image. It has instructions for both tagging the image with `docker tag` and pushing with `docker push`.  Upon pushing you should be able to see the image inside of the Artifact Registry within the repository directory.

## Deploying via the command line

Navigate to https://console.cloud.google.com/run?enableapi=true&_ga=2.240046915.1376948655.1666123756-1669168203.1659107311 which takes you to the Cloud Run service creation tool.

1. Select "Deploy one revision from an existing container image". Select the container image URL by clicking SELECT in the box, clicking the Artifact Registry tab, choosing the correct image, then SELECT.
2. For "CPU Allocation and Pricing", leave as "CPU is only allocated during request processing".
3. For "Autoscaling", leave as is.
4. For Ingress, leave as "Allow all traffic"
5. For Authentication, select "Allow unauthorized invocations". Since the gEAR API is part of a public codebase, we do not want to expose authentication credentials, so it is necessary to set it up so anything can call this service.  This should only be the projectR API call though, in reality.
  a. Currently this can only be changed if the person setting up the service has permissions to change IAM policy.
  b. I currently do not so I am leaving as is and providing auth credentials to test locally.
6. Click "Container, Connections, Security" dropdown
7. Change memory to 32Gb.  Unfortunately this means we need to set 8 CPUs as well. Mostly projectR runs should never hit 32Gb but we have had some larger runs that would.
8. Set the request timeout to 1200 seconds. The longest request I've seen so far was about 17-18 minutes, but that is more the rarity than the norm.  Most should only take a 1-3 minutes.
9. Set concurrency to 1. This is because of the odd projectR job that takes a huge amount of memory and we do not wish to have jobs collide on the same instance.
10. Ready to go, click "Create"
