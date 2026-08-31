# Clone template VM to other project areas

In the gear project area, there is an "instance-template" VM that serves as a golden VM that can be used to spawn off new instances of gEAR. It is vital to ensure this VM is updated with the latest code, packages, etc before starting this process.  This image does not have any mountable disks on it.

To clone a VM across GCP projects, the best practice for a "golden VM" pattern is to create a **Custom Image** in the source project, share IAM permissions on that image, and then launch the new instance in the target project.

## Using a Custom Image (Recommended for Golden VMs)

Creating a Custom Image captures the boot disk state and allows you to spin up multiple identical instances across projects.

### 1. Stop the Source Instance

To ensure data consistency on the boot disk:

```bash
gcloud compute instances stop instance-template \
    --project=SOURCE_PROJECT_ID \
    --zone=SOURCE_ZONE

```

### 2. Create a Custom Image from the Boot Disk

Create an image using the boot disk of `instance-template`:

```bash
gcloud compute images create golden-vm-image \
    --source-disk=instance-template \
    --source-disk-zone=SOURCE_ZONE \
    --project=SOURCE_PROJECT_ID

```

### 3. Share the Image with the Target Project / User

Grant the **Compute Image User** role (`roles/compute.imageUser`) on the source project to either your account or the Compute Engine service account of the target project:

```bash
# Share with a specific user:
gcloud projects add-iam-policy-binding SOURCE_PROJECT_ID \
    --member="user:your-email@example.com" \
    --role="roles/compute.imageUser"

# OR share with the Target Project's Default Compute Service Account:
# (Format: TARGET_PROJECT_NUMBER-compute@developer.gserviceaccount.com)
gcloud projects add-iam-policy-binding SOURCE_PROJECT_ID \
    --member="serviceAccount:TARGET_PROJECT_NUMBER-compute@developer.gserviceaccount.com" \
    --role="roles/compute.imageUser"

```

### 4. Create the New VM in the Target Project

Reference the image in the source project using its full resource path (`projects/SOURCE_PROJECT_ID/global/images/IMAGE_NAME`):

```bash
gcloud compute instances create cloned-vm-instance \
    --project=TARGET_PROJECT_ID \
    --zone=TARGET_ZONE \
    --machine-type=e2-standard-2 \
    --image=projects/SOURCE_PROJECT_ID/global/images/golden-vm-image

```
