const FILE_UPLOAD_LIMIT = 10_485_760; // 10 MB

let screenshot = null;

window.onload = () => {
  // generate list of tags from database
  get_tag_list();

  //TODO prepopulate email field if user is logged in
};

/*jslint unparam: true */
/*global window, $ */
$(function () {
  "use strict";

  // upload file to server
  const url = "contact_screenshots/";
  $("#screenshot").fileupload({
    url,
    dataType: "json",
    maxFileSize: FILE_UPLOAD_LIMIT,
    done(_e, data) {
      console.error("file successfully uploaded");
      screenshot = null;
      if (data.result.files.length == 0) {
        console.error("File Upload Error");
        return;
      }
      // Collect name of screenshot. Should be just one file
      screenshot = data.result.files[0].name;
      $("#upload_status").show();
      $("#upload_status").text(`${screenshot} uploaded successfully!`);
      $("#upload_status").addClass("text-success");
      $("#upload_status").removeClass("text-danger");
    },
    fail(_e, data) {
      console.error("file upload failure.");
      console.error(data.errorThrown);
      $("#upload_status").show();
      $("#upload_status").text(`${screenshot} failed to upload.`);
      $("#upload_status").addClass("text-danger");
      $("#upload_status").removeClass("text-success");
    },
  });
});

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

$("#btn_submit_comment").click((_e) => {
  // If user confirmed private data in the form, skip the modal
  if ($('#private_check').is(':checked')) {
    $("#actual_submit").click();
  } else {
    $('#confirmModal').modal("show");
  }

});

// Modal "OK" was clicked
$("#actual_submit").click((e) => {
  $('#confirmModal').modal("hide");

  //remove any existing warnings
  $(".label-warning").remove();
  e.preventDefault();

  // check if required fields were filled in
  const pass = check_required_fields();

  if (pass == true) {
    //submit the upload form
    const formData = new FormData($("#upload_form")[0]);
    // Replace the binary screenshot with the filename, since the file is uploaded.
    formData.delete("files[]");
    formData.append("screenshot", screenshot);
    formData.append(
      "private_check",
      $('#private_check').is(':checked')
    );

    // since comment tags come from "select2" now,
    // need to format this attribute to only pass a single array instead of multiple scalars
    formData.delete("comment_tag");
    comment_tags = $('#comment_tag').select2('data').map((elem) => elem.id);
    formData.append("comment_tag", comment_tags);

    $.ajax({
      url: "./cgi/create_github_issue.cgi",
      type: "POST",
      data: formData,
      dataType: "json",
      contentType: false, // contentType and processData need to be false when using FormData objects and jQuery AJAX
      processData: false,
      success(data, _textStatus, _jqXHR) {
        if (data.success == 1) {
          //activate a 'Thank you' modal
          $("#successModal").modal("show");
          // redirect to the data manager page
          setTimeout(() => {
            window.location.replace("index.html");
          }, 2000);
          return;
        }
        const msg =
        "Something went wrong. Please try again, but if this continues please contact jorvis@gmail.com for help.";
        const error = `Error during creation: ${data.error}`;
        display_error_bar(msg, error);
      },
      error(jqXHR, _textStatus, errorThrown) {
        const msg =
        "Something went wrong. Please try again, but if this continues please contact jorvis@gmail.com for help.";
        const error = `${jqXHR.status} ${errorThrown}`;
        display_error_bar(msg, error);
      },
    }); //end ajax
  } else {
    console.error("Some required fields were left blank.");
  }
}); //end #btn_submit_commit

//remove warning label on focus
$("input, textarea").focus(function () {
  const warning_el = `#label-warning-${$(this).attr("id")}`;
  $(warning_el).remove();
});

//checks required fields for input.
function check_required_fields() {
  let pass = true;

  //test whether required fields were filled in
  if (!$("#submitter_firstname").val()) {
    $("#label_submitter_firstname").append(
      ' <span id="label-warning-submitter_firstname" class="label label-warning">Oops. This is required.</span>'
    );
    $("#submitter_firstname").effect("highlight", "slow");
    pass = false;
  }
  if (!$("#submitter_lastname").val()) {
    $("#label_submitter_lastname").append(
      ' <span id="label-warning-submitter_lastname" class="label label-warning">Oops. This is required.</span>'
    );
    $("#submitter_lastname").effect("highlight", "slow");
    pass = false;
  }
  if (!$("#submitter_email").val()) {
    $("#label_submitter_email").append(
      ' <span id="label-warning-submitter_email" class="label label-warning">Oops. This is required.</span>'
    );
    $("#submitter_email").effect("highlight", "slow");
    pass = false;
  }

  if (!$("#comment_title").val()) {
    $("#label_comment_title").append(
      ' <span id="label-warning-comment_title" class="label label-warning">Oops. This is required.</span>'
    );
    $("#comment_title").effect("highlight", "slow");
    pass = false;
  }
  if (!$("#comment").val()) {
    $("#label_comment").append(
      ' <span id="label-warning-comment" class="label label-warning">Oops. This is required.</span>'
    );
    $("#comment").effect("highlight", "slow");
    pass = false;
  }
  if ($("#super_impressive_security_check").val() != 28) {
    $("#label_super_impressive_security_check").append(
      ' <span id="label-warning-super_impressive_security_check" class="label label-warning">Oops. This is required.</span>'
    );
    $("#super_impressive_security_check").effect("highlight", "slow");
    return false;
  }

  return pass;
}

// closed alert bar
$(".alert-container").on("click", "button.close-alert", () => {
  $(".alert-container").hide();
});

// error should be html message for user. Example: error = '<p>You cannot do that.</p>'
function display_error_bar(msg, error) {
  $(".alert-container")
    .html(
      '<div class="alert alert-danger alert-dismissible" role="alert">' +
        '<button type="button" class="close close-alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
        '<p class="alert-message">' +
        "<strong>Oops. </strong> " +
        msg +
        "</p>" +
        '<p style="text-align: right;">(<em>Error: ' +
        error +
        "</em>)</p>" +
        "</div>"
    )
    .show();
}
