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