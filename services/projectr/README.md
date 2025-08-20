# ProjectR Container service

This directory is meant to house the code to run projectR inside of a container. This container can then be deployed to a Google Cloud service, such as Google Cloud Run

Details - https://cloud.google.com/run/docs/quickstarts/build-and-deploy/deploy-python-service

NOTE: Cloud Build times out when attempting to build this image. Need to build locally. For my purposes, I named my image "projectr_service".

Also if building on an M1 Mac, you need to add platform information to make it build with Linux architecture.  Otherwise it will not be able to run in the Google Cloud services.  My docker command to build was `docker build --platform linux/amd64 --no-cache -t projectr_service .`

## Pushing Docker images to Google Artifact Registry

Google Artifact Registry has recently replaced Google Container Registry for the de facto location to store Docker container images along with other artifact types. When selecting the Artifact Registry from the Google Cloud web console, you will see two tabs: Artifact Registry and Container Registry. Container Registry is just a holdover from the older service, so we will focus on Artifact Registry

To create a repository for storing images, on the Artifact Registry console select Create Repository and give your service a name. I named it "cloud-run-source-deploy" since I intend on potentially hosting multiple Cloud Run Services there, not just projectR.  Also set the region to us-east1.

We need to have our Docker configuration file modified to support gcloud credentials when pushing images. Since we set our Cloud Run region to us-east1, we need to run `gcloud auth configure-docker us-east1-docker.pkg.dev`, which modifies `~/.docker/config.json`

To push our created image, we need to give it a specific tag after building. Follow the instructions at https://cloud.google.com/artifact-registry/docs/docker/pushing-and-pulling#pushing for pushing the image. It has instructions for both tagging the image with `docker tag` and pushing with `docker push`.  Upon pushing you should be able to see the image inside of the Artifact Registry within the repository directory.

## Deploying the Cloud Run service

Navigate to https://console.cloud.google.com/run?enableapi=true&_ga=2.240046915.1376948655.1666123756-1669168203.1659107311 which takes you to the Cloud Run service creation tool.

1. Select "Deploy one revision from an existing container image". Select the container image URL by clicking SELECT in the box, clicking the Artifact Registry tab, choosing the correct image, then SELECT.
2. For "CPU Allocation and Pricing", leave as "CPU is only allocated during request processing".
3. For "Autoscaling", leave as is.
4. For Ingress, leave as "Allow all traffic"
5. For Authentication, select "Allow unauthorized invocations". Since the gEAR API is part of a public codebase, we do not want to expose authentication credentials, so it is necessary to set it up so anything can call this service.  This should only be the projectR API call though, in reality.
  a. Currently this can only be changed if the person setting up the service has permissions to change IAM policy.
6. Click "Container, Connections, Security" dropdown.  Each of those are their own separate tab
7. Change memory to 16Gb.  Unfortunately this means we need to set 4 CPUs as well. Most projectR runs should never hit 16Gb due to chunking the file, but it's possible.
8. Set the request timeout to 1200 seconds. The longest request I've seen so far was about 17-18 minutes, but that is more the rarity than the norm.  Most should only take a 1-3 minutes.
9. Set concurrency to 1. This is because of the odd projectR job that takes a huge amount of memory and we do not wish to have jobs collide on the same instance.
10. Click "Security" tab and change the Service account to "gear" instead of the default Compute Engine one.
11. Ready to go, click "Create"

## Getting the service URL

After the projectR service is created, you will be redirected to a dashboard with various tabs both for the service and the current revision. Beside the "projectR_service" name above the tabs, there is a URL that can be copied.  Paste this URL in the gear.ini file under the "projectR_service" section under the "hostname" param.

## Deploying a revised service

You can use Docker to build and push revised images up to the Artifact Registry. In the "projectr-service" dashboard in Cloud Run, there is a button near the top called "Edit and Deploy New Revision" that will let you take deploy the latest image for the service. All of the previous set configurations are saved, but editable if you want, and there is an extra checkbox to reroute all traffic to this new revision after launching.  You do not need to change your URL for the API endpoint in the gear.ini file. The previous versions of the service will still show in the dashboard and can be deleted if you want.
