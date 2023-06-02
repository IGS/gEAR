#!/opt/bin/python3

from flask import request, abort
from flask_restful import Resource
import os,sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb


# Parse gEAR config
# https://stackoverflow.com/a/35904211/1368079
this = sys.modules[__name__]
from gear.serverconfig import ServerConfig
this.servercfg = ServerConfig().parse()

class Submission(Resource):
    """Requests to deal with a single submission"""

    def get(self, submission_id):
        """Retrieve a particular submission based on ID."""
        result = {"self": request.path, "success":False, "datasets":{}, "layout_share_id":None}

        session_id = request.cookies.get('session_id')
        # Must have a gEAR account to upload datasets
        user = geardb.get_user_from_session_id(session_id)

        submission = geardb.get_submission_by_id(submission_id)
        if not submission:
            return result

        layout = submission.get_layout_info()
        if layout:
            result["layout_share_id"] = layout.share_id

        sd = geardb.SubmissionDatasetCollection()
        submission.datasets = sd.get_by_submission_id(submission_id)

        for s_dataset in submission.datasets:
            dataset_id = s_dataset.dataset_id
            result["datasets"][dataset_id] = {"href":f"/api/submission/{submission_id}/datasets/{dataset_id}"}
        result["success"] = True
        return result

    def post(self, submission_id):
        """Start the import process for this submission."""
        session_id = request.cookies.get('session_id')
        # Must have a gEAR account to upload datasets
        user = geardb.get_user_from_session_id(session_id)

        req = request.get_json()
        metadata = req.get("metadata")
        action = req.get("action")

        submission = geardb.get_submission_by_id(submission_id)
        submission_datasets = geardb.SubmissionDatasetCollection()
        submission_datasets.get_by_submission_id(submission_id=submission_id)

        if action == "import":
            for s_dataset in submission_datasets:
                pass

            if submission.email_updates:
                pass
                # send email about submission

        abort(400)



        """
        const fileEntities = getFileEntities(jsonData);
        const sampleEntities = getSampleEntities(jsonData);

        // Initialize db information about this submission
        // Check db if datasets were loaded in previous submissions
        // (initialize_new_submission.cgi)
        try {
            await initializeNewSubmission(fileEntities, submissionId);
        } catch (error) {
            alert("Something went wrong with saving this submission to database. Please contact gEAR support.");
            console.error(error);
        }

        await Promise.allSettled(fileEntities.map(async (entity) => {
            const { attributes: fileAttributes } = entity;
            // Merge in sample attributes
            const sampleAttributes = getSampleAttributesByName(sampleEntities, fileAttributes.sample_id);
            const allAttributes = { dataset: fileAttributes, sample: sampleAttributes };
            allAttributes.name = entity.name;

            // Create JSON from metadata
            const fileMetadata = await getFileMetadata(allAttributes.dataset);

            // Merge all API metadata and metadata pulled from the portal
            // TODO: Pull sample just for specific dataset
            const allMetadata = {
                dataset: { ...allAttributes.dataset, ...fileMetadata.dataset },
                sample: { ...allAttributes.sample, ...fileMetadata.sample }
            };

            const datasetId = allMetadata.dataset.id;

            // If import has already finished, let's just focus on the display
            if (isImportComplete(datasetId) && await getDefaultDisplay(datasetId)) {
                return true;
            }

            const params = { "metadata": allMetadata};
            const response = await fetch(`/api/submission_dataset/${datasetId}`, {
                method: "POST",
                headers: {"Content-Type": "application/json"},
                body: JSON.stringify(params)
                });
            if (!response?.ok) {
                throw new Error(postResponse.statusText);
            }
            const jsonRes = await response.json();
            if (!jsonRes.success) {
                printTblMessage(datasetId, jsonRes.message, "danger");
                partialSuccess.classList.remove("is-hidden");
                return false
            }

            // Save the newly created display as a default.
            const plotType = "tsne_static";
            const label = "nemoanalytics import default plot";
            const displayId = await saveNewDisplay(datasetId, jsonRes.plot_config, plotType, label);
            await saveDefaultDisplay(datasetId, displayId);

            // enable select-obs dropdown
            if (["MEX"].includes(allMetadata.dataset.filetype)) {
                const parent = document.getElementById(`${datasetId}-obslevels`);
                parent.textContent = "Categories not found";
            } else {
                populateObsDropdown(datasetId);
            }
            // show permalink
            showSubmissionPermalink(datasetId);

            // After the first dataset is finished, we can View Datasets if desired
            return true
        }));
        """

class SubmissionEmail(Resource):
    def put(self, submission_id):
        """Update email updates."""
        if not submission_id:
            abort(404)  # "Submissiom ID not found for this request."

        result = {"success": False, "message":""}

        try:
            submission = geardb.get_submission_by_id(submission_id)
            submission.save_change(attribute="email_updates", value=1)
            result["success"] = True
        except Exception as e:
            result["message"] = str(e)

        return result

class Submissions(Resource):
    """Requests to deal with multiple submissions, including creating new ones"""

    def post(self):
        """Create a new empty submission in the database."""
        # Must have a gEAR account to create submissions
        session_id = request.cookies.get('gear_session_id')
        user = geardb.get_user_from_session_id(session_id)

        req = request.get_json()
        submission_id = req.get("submission_id")
        is_restricted = req.get("is_restricted")

        result = {"self": request.path, "success":False, "datasets":{}, "layout_share_id":None, "message":""}

        submission = geardb.get_submission_by_id(submission_id)
        if submission:
            abort(400, message="Submission already exists for this ID.")

        try:
            conn = geardb.Connection()
            cursor = conn.get_cursor()

            qry = """
                INSERT INTO submission (id, user_id, is_finished, is_restricted, date_added)
                VALUES (%s, %s, 0, %s, NOW())
            """

            cursor.execute(qry, (submission_id, user.id, is_restricted))
            cursor.close()
            conn.commit()
            conn.close()
        except:
            abort(400, message="Could not save new submission to database")

        result["href"] = result["self"] + f"/${submission_id}" # This can be GET-able

        result["success"] = True
        result["message"] = "Submission database entries have been populated"
        return result

