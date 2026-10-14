'use strict';

import { closeModal, createToast, getCurrentUser, initCommonUI, logErrorInConsole, openModal, registerPageSpecificLoginUIUpdates } from "./common.v2.js";

const FILE_UPLOAD_LIMIT = 10_485_760; // 10 MB
const SECURITY_CHECK_ANSWER = "28";   // 7 x 4

/**
 * Required form fields and the message shown when each one fails validation.
 * @type {Array<{id: string, message: string, isValid?: (input: HTMLInputElement) => boolean}>}
 */
const REQUIRED_FIELDS = [
  { id: "submitter-fullname", message: "Please enter your name." },
  { id: "submitter-email", message: "Please enter a valid email address.", isValid: (input) => !input.validity.typeMismatch },
  { id: "comment-title", message: "Please enter a title." },
  { id: "comment", message: "Please enter your question or comment." },
  { id: "super-impressive-security-check", message: "Please enter the correct answer (7 x 4).", isValid: (input) => input.value.trim() === SECURITY_CHECK_ANSWER },
];

/** Name of the screenshot file after a successful upload to contact_screenshots/, or null if none. */
let screenshot = null;

/** Tom Select instance for the tags multi-select. */
let tagSelect = null;

const uploadForm = document.getElementById("upload-form");
const submitBtn = document.getElementById("btn-submit-comment");
const actualSubmitBtn = document.getElementById("actual-submit");
const confirmModal = document.getElementById("confirm-modal");
const successModal = document.getElementById("success-modal");
const privateCheck = document.getElementById("private-check");
const alertContainer = document.getElementById("alert-c");
const alertErrorElt = document.getElementById("alert-error");

/**
 * Creates the Tom Select widget for the tags field.
 * Created before the tag list is fetched so users can still type their own tags if that request fails.
 * @returns {TomSelect} The Tom Select instance bound to #comment-tag.
 */
const initTagSelect = () => {
  const select = new TomSelect("#comment-tag", {
    plugins: ["remove_button"],
    persist: false,
    create: true,
    selectOnTab: true,
    placeholder: "Examples: mouse, sox2, upload issue",
  });

  // Tom Select replaces the <select> with its own input, so re-link the help text to it
  select.control_input.setAttribute("aria-describedby", "comment-tag-help");
  return select;
}

/**
 * Fetches existing tags from the database and adds them as autocomplete options in the tags field.
 * Failure is not fatal; the error is only logged to the console.
 * @returns {Promise<void>}
 */
const loadTagOptions = async () => {
  try {
    const { data } = await axios.get("./cgi/get_tag_list.cgi");

    if (data?.success !== 1) {
      throw new Error("Retrieving tags was not successful. Please check the server logs for more information.");
    }

    tagSelect.addOptions(data.tags.map(({ label }) => ({ value: label, text: label })));
  } catch (error) {
    // Not having populated tags is not a fatal error, but log it in the console for debugging
    console.error("Unable to generate tag list");
    logErrorInConsole(error);
  }
}

/**
 * Shows or hides the error notification above the form.
 * @param {string|null} message - Error detail to show, or null to hide the notification.
 */
const setAlertError = (message) => {
  alertErrorElt.textContent = message ? `(Error: ${message})` : "";
  alertContainer.classList.toggle("is-hidden", !message);
}

/**
 * Submits the contact form to create a GitHub issue.
 * Shows the success modal and redirects home on success, or the error notification on failure.
 * Assumes the form has already passed validateRequiredFields().
 * @returns {Promise<void>}
 */
const submitContactMessage = async () => {
  setAlertError(null);

  const formData = new FormData(uploadForm);

  // The screenshot was already uploaded on file selection, so send its filename instead of the file.
  formData.delete("files[]");
  if (screenshot) {
    formData.set("screenshot", screenshot);
  }

  // The CGI expects the string "true" or "false" (an unchecked checkbox would otherwise be omitted)
  formData.set("private_check", String(privateCheck.checked));

  // Tom Select keeps the original <select> in sync, so each tag is a separate "comment_tag" entry.
  // The CGI only reads the first value, so send one comma-separated string.
  formData.set("comment_tag", formData.getAll("comment_tag").join(", "));

  submitBtn.classList.add("is-loading");
  submitBtn.disabled = true;

  try {
    // Axios sets the multipart Content-Type (with boundary) automatically for FormData
    const { data } = await axios.post("./cgi/create_github_issue.cgi", formData);

    if (data?.success === 1) {
      openModal(successModal);
      // redirect to the home page. assign() (not replace()) keeps the contact page in the browser history,
      // so "Back" returns here instead of skipping past it.
      setTimeout(() => {
        window.location.assign("index.html");
      }, 2000);
      return;   // leave the submit button disabled so the issue is not submitted twice
    }

    setAlertError(`Error during creation: ${data?.error || "unknown error"}`);
  } catch (error) {
    logErrorInConsole(error);
    setAlertError(error.message);
  }

  alertContainer.scrollIntoView({ behavior: "smooth", block: "nearest" });
  submitBtn.classList.remove("is-loading");
  submitBtn.disabled = false;
};

/**
 * Marks a form field as invalid and shows an error message under it.
 * The message is linked with aria-describedby so screen readers announce it with the field.
 * @param {HTMLInputElement|HTMLTextAreaElement} input - The invalid field.
 * @param {string} message - The error message to display.
 */
const showFieldError = (input, message) => {
  const errorId = `${input.id}-error`;
  let errorElt = document.getElementById(errorId);

  if (!errorElt) {
    errorElt = document.createElement("p");
    errorElt.id = errorId;
    errorElt.className = "help has-text-danger-dark";
    input.closest(".field").appendChild(errorElt);
  }
  errorElt.textContent = message;

  input.classList.add("is-danger");
  input.setAttribute("aria-invalid", "true");

  const describedBy = (input.getAttribute("aria-describedby") || "").split(" ").filter(Boolean);
  if (!describedBy.includes(errorId)) {
    input.setAttribute("aria-describedby", [errorId, ...describedBy].join(" "));
  }
}

/**
 * Removes the error state and message from a form field, if present.
 * @param {HTMLInputElement|HTMLTextAreaElement} input - The field to clear.
 */
const clearFieldError = (input) => {
  const errorId = `${input.id}-error`;
  document.getElementById(errorId)?.remove();

  input.classList.remove("is-danger");
  input.removeAttribute("aria-invalid");

  const describedBy = (input.getAttribute("aria-describedby") || "").split(" ").filter((id) => id && id !== errorId);
  if (describedBy.length) {
    input.setAttribute("aria-describedby", describedBy.join(" "));
  } else {
    input.removeAttribute("aria-describedby");
  }
}

/**
 * Checks that all required fields are filled in (and valid), showing an error on each field that is not.
 * Focuses the first invalid field so keyboard and screen reader users land on it.
 * @returns {boolean} True if every required field is valid.
 */
const validateRequiredFields = () => {
  let firstInvalid = null;

  for (const field of REQUIRED_FIELDS) {
    const input = document.getElementById(field.id);
    const isValid = input.value.trim() !== "" && (field.isValid?.(input) ?? true);

    if (isValid) {
      clearFieldError(input);
      continue;
    }

    showFieldError(input, field.message);
    firstInvalid ??= input;
  }

  firstInvalid?.focus();
  return firstInvalid === null;
}

/**
 * Pre-fills the name and email fields from the logged-in user's account, if any.
 */
const prefillUserFields = () => {
  const currentUser = getCurrentUser();
  if (currentUser?.session_id) {
    document.getElementById("submitter-fullname").value = currentUser.user_name;
    document.getElementById("submitter-email").value = currentUser.email;
  }
}

/**
 * Returns the page to a blank form after a successful submission.
 * Needed when the browser restores the page from its back/forward cache, which keeps the
 * previous state (open success modal, disabled submit button, filled-in fields).
 */
const resetContactForm = () => {
  uploadForm.reset();
  tagSelect.clear();
  prefillUserFields();

  screenshot = null;
  document.getElementById("screenshot-file-name").textContent = "No file chosen";
  document.getElementById("upload-status").classList.add("is-hidden");

  for (const field of REQUIRED_FIELDS) {
    clearFieldError(document.getElementById(field.id));
  }
  setAlertError(null);
  closeModal(successModal);

  submitBtn.classList.remove("is-loading");
  submitBtn.disabled = false;
}

/**
 * Page-specific UI updates run by common.v2.js once the login state is known.
 * Pre-fills the name and email fields for logged-in users and loads the tag options.
 * @returns {Promise<void>}
 */
const handlePageSpecificLoginUIUpdates = async () => {
  document.getElementById('page-header-label').textContent = 'Contact Us!';

  prefillUserFields();

  // generate list of tags from database
  await loadTagOptions();
}

tagSelect = initTagSelect();
registerPageSpecificLoginUIUpdates(handlePageSpecificLoginUIUpdates);

// Pre-initialize some stuff
await initCommonUI();

// Handle screenshot upload
document.getElementById("screenshot").addEventListener("change", async (event) => {
  const file = event.target.files[0];
  const fileNameElt = document.getElementById("screenshot-file-name");
  const uploadStatusElt = document.getElementById("upload-status");

  // Any previously uploaded screenshot is replaced by this selection
  screenshot = null;
  fileNameElt.textContent = file ? file.name : "No file chosen";
  uploadStatusElt.classList.add("is-hidden");
  uploadStatusElt.classList.remove("has-text-success", "has-text-danger-dark");

  if (!file) return;

  /**
   * Shows a failed-upload message and clears the file input.
   * @param {string} message - Text for the status line under the file input.
   */
  const showUploadFailure = (message) => {
    uploadStatusElt.textContent = message;
    uploadStatusElt.classList.remove("is-hidden", "has-text-success");
    uploadStatusElt.classList.add("has-text-danger-dark");
    event.target.value = "";
    fileNameElt.textContent = "No file chosen";
  }

  // The "accept" attribute is only a hint to the file picker, so check the type here too
  if (!(file.type.startsWith("image/") || file.type === "application/pdf")) {
    showUploadFailure(`${file.name} is not an image or PDF file. Please choose a different file.`);
    return;
  }

  if (file.size > FILE_UPLOAD_LIMIT) {
    showUploadFailure(`${file.name} is larger than the ${FILE_UPLOAD_LIMIT / 1_048_576} MB limit. Please choose a smaller file.`);
    return;
  }

  uploadStatusElt.textContent = "Uploading...";
  uploadStatusElt.classList.remove("is-hidden");

  // The upload handler (blueimp jQuery-File-Upload PHP) expects multipart form data in the "files[]" field
  const uploadData = new FormData();
  uploadData.append("files[]", file);

  try {
    const { data } = await axios.post("contact_screenshots/", uploadData);
    const uploaded = data?.files?.[0];

    if (!uploaded || uploaded.error) {
      throw new Error(uploaded?.error || "No file information returned by the upload handler");
    }

    // Collect name of screenshot. Should be just one file
    screenshot = uploaded.name;
    uploadStatusElt.textContent = `${screenshot} uploaded successfully!`;
    uploadStatusElt.classList.add("has-text-success");
  } catch (error) {
    logErrorInConsole(error);
    showUploadFailure(`${file.name} failed to upload. Please try again.`);
  }
});

document.getElementById("btn-submit-comment").addEventListener("click", async (event) => {
  // Prevent the default form submission behavior
  event.preventDefault();

  // Validate before asking for confirmation, so the user does not confirm a form that cannot be sent
  if (!validateRequiredFields()) {
    createToast("Some required fields were left blank or are invalid. Please fix them and try again.", "is-danger");
    return;
  }

  // If user confirmed private data in the form, skip the modal
  if (privateCheck.checked) {
    await submitContactMessage();
    return;
  }

  openModal(confirmModal);
  // Move focus into the dialog; "Cancel" is the safe default since "OK" posts publicly
  confirmModal.querySelector(".modal-card-foot .js-modal-close").focus();
});

// Modal "OK" was clicked
actualSubmitBtn.addEventListener("click", async (event) => {
  event.preventDefault();

  closeModal(confirmModal);
  await submitContactMessage();
});

// Hide (rather than remove) the error notification when dismissed, so it can be shown again on a later failure.
// This listens in the capture phase so it runs before the common.v2.js ".notification .delete" handler,
// which removes the notification from the page entirely.
alertContainer.addEventListener("click", (event) => {
  if (!event.target.closest(".delete")) return;
  event.stopPropagation();
  setAlertError(null);
}, true);

// Returning with "Back" after a successful submission may restore this page exactly as it was left
// (back/forward cache), so start over with a blank form instead of the "Success!" modal.
window.addEventListener("pageshow", (event) => {
  if (event.persisted && successModal.classList.contains("is-active")) {
    resetContactForm();
  }
});

// Remove a field's error message once the user starts correcting it.
// ("input" rather than "focus", since validateRequiredFields() focuses the first invalid field.)
for (const field of REQUIRED_FIELDS) {
  const input = document.getElementById(field.id);
  input.addEventListener("input", () => clearFieldError(input));
}
