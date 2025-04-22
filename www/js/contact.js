'use strict';

const FILE_UPLOAD_LIMIT = 10_485_760; // 10 MB
let screenshot = null;

let tokenizedTags = [];
let enteredTags = [];
const tagField = document.getElementById("tag-field");
const addTag = document.getElementById("add-tag");
const tagInput = document.getElementById("tag-input");

window.onload=function() {
    // Set the page title
    document.getElementById('page-header-label').textContent = 'Contact us';

    // generate list of tags from database
    get_tag_list();

    document.getElementById("btn_submit_comment").addEventListener("click", (_e) => {
        // If user confirmed private data in the form, skip the modal
        if (document.getElementById("private_check").checked) {
            document.getElementById("actual_submit").click();
        } else {
            
        }
    });

    //remove warning label on focus
    document.querySelectorAll("input, textarea").forEach((element) => {
        element.addEventListener("focus", () => {
            const warningEl = document.getElementById(`label-warning-${element.id}`);
            if (warningEl) {
                warningEl.remove();
            }
        });
    });

    //TODO prepopulate email field if user is logged in
};

addTag.addEventListener("click", () => {
    if (tagInput.value !== "") {
      createTag(tagInput.value);
    }
    tagInput.value = "";
});
  
tagInput.addEventListener("keyup", (event) => {
    if ((event.keyCode === 13) && (tagInput.value !== "")) {
      createTag(tagInput.value);
      tagInput.value = "";
    }
});

//checks required fields for input.
function check_required_fields() {
    let pass = true;

    // Helper function to add warning labels
    function addWarningLabel(inputId, labelId, message) {
        const label = document.getElementById(labelId);
        if (!document.getElementById(`label-warning-${inputId}`)) {
            const warningSpan = document.createElement("span");
            warningSpan.id = `label-warning-${inputId}`;
            warningSpan.className = "label label-warning";
            warningSpan.textContent = message;
            label.appendChild(warningSpan);
        }
        document.getElementById(inputId).classList.add("highlight");
    }

    // Helper function to highlight fields
    function highlightField(inputId) {
        const field = document.getElementById(inputId);
        field.classList.add("highlight");
        setTimeout(() => field.classList.remove("highlight"), 1000);
    }

    // Test whether required fields were filled in
    if (!document.getElementById("submitter_firstname").value) {
        addWarningLabel("submitter_firstname", "label_submitter_firstname", "Oops. This is required.");
        highlightField("submitter_firstname");
        pass = false;
    }
    if (!document.getElementById("submitter_lastname").value) {
        addWarningLabel("submitter_lastname", "label_submitter_lastname", "Oops. This is required.");
        highlightField("submitter_lastname");
        pass = false;
    }
    if (!document.getElementById("submitter_email").value) {
        addWarningLabel("submitter_email", "label_submitter_email", "Oops. This is required.");
        highlightField("submitter_email");
        pass = false;
    }
    if (!document.getElementById("comment_title").value) {
        addWarningLabel("comment_title", "label_comment_title", "Oops. This is required.");
        highlightField("comment_title");
        pass = false;
    }
    if (!document.getElementById("comment").value) {
        addWarningLabel("comment", "label_comment", "Oops. This is required.");
        highlightField("comment");
        pass = false;
    }
    if (document.getElementById("super_impressive_security_check").value != 28) {
        addWarningLabel("super_impressive_security_check", "label_super_impressive_security_check", "Oops. This is required.");
        highlightField("super_impressive_security_check");
        return false;
    }

    return pass;
}

/* Credit (adapted from)
    https://codepen.io/yanniskatsaros/pen/ZEQByPj
*/
function createTag(message) {
    const controlDiv = document.createElement("div");
    controlDiv.classList.add("control");
    
    const tags = document.createElement("div");
    tags.classList.add("tags", "has-addons");
    
    const tagContent = document.createElement("a");
    tagContent.classList.add("tag", "is-link");
    tagContent.innerText = message;
    
    const tagDelete = document.createElement("a");
    tagDelete.classList.add("tag", "is-delete");
    tagDelete.addEventListener("click", (event) => {
      tagField.removeChild(controlDiv);
    });
    
    // finally nest all the tags together
    tags.appendChild(tagContent);
    tags.appendChild(tagDelete);
    controlDiv.appendChild(tags);
    tagField.appendChild(controlDiv);

    enteredTags.push(message);
    
    /* Each adds this HTML to the container:
    <div class="tags has-addons">
        <a class="tag is-link">${message}</a>
        <a class="tag is-delete"></a>
    </div>
    */
}

// Generate a list of existing tags from database
function get_tag_list() {

    fetch("./cgi/get_tag_list.cgi")
      .then(response => response.json())
      .then(data => {
        if (data.success == 1) {
          // put all tags into a list
          tokenizedTags = data.tags.map(tag => tag.label);

          // Create dropdown of tags for autocomplete
          /*
          const commentTag = document.getElementById("comment_tag");
          commentTag.setAttribute("placeholder", "Examples: mouse sox2 upload issue");
          commentTag.addEventListener("input", (event) => {
            const inputValue = event.target.value;
            const suggestions = tokenized_tags.filter(tag => tag.includes(inputValue));
            // Logic to display suggestions can be added here
          });
          */
        } else {
          console.error("Handle a failed report from the CGI");
        }
      })
      .catch(() => {
        console.error("Something broke. Unable to generate tag list");
      });
  } //end get_tag_list()

const handlePageSpecificLoginUIUpdates = async (event) => {
    // Nothing to do here at the moment
}