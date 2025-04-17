'use strict';

const FILE_UPLOAD_LIMIT = 10_485_760; // 10 MB
let screenshot = null;

window.onload=function() {
    // Set the page title
    document.getElementById('page-header-label').textContent = 'Contact us';

    // generate list of tags from database
    get_tag_list();

    //TODO prepopulate email field if user is logged in
};

// Generate a list of existing tags from database
// Autocompletes in TAG input as twitter-like tokens
function get_tag_list() {
    $.get("./cgi/get_tag_list.cgi", (data) => {
      if (data.success == 1) {
        // put all tags into a list
        const tokenized_tags = [];
        for (i = 0; i < data.tags.length; i++) {
          tokenized_tags.push(data.tags[i].label);
        }
  
        // Create dropdown of tags for autocomplete
        $("#comment_tag").select2({
          allowClear: true,
          data: tokenized_tags,
          placeholder: "Examples: mouse sox2 upload issue",
          tags: true,
          tokenSeparators: [",", " "]
        })
  
      } else {
        console.error("Handle a failed report from the CGI");
      }
    }).fail(() => {
      console.error("Something broke. Unable to generate tag list");
    });
  } //end get_tag_list()

const handlePageSpecificLoginUIUpdates = async (event) => {
    // Nothing to do here at the moment
}