
window.onload=function() {
    // check if the user is already logged in
    check_for_login();

    // generate list of tags from database
    get_tag_list();

	//TODO prepopulate email field if user is logged in

};

// Generate a list of existing tags from database
// Autocompletes in TAG input as twitter-like tokens
function get_tag_list() {
  	$.get('./cgi/get_tag_list.cgi', function(data) {
    		if (data['success'] == 1) {
    		    // put all tags into a list
            var tokenized_tags = [];
            for (i=0; i < data['tags'].length; i++) {
                tokenized_tags.push(data['tags'][i]['label']);
            }

            //initialize tokenfield for tags
      			$('#comment_tag').tokenfield({
      			    autocomplete: {
      			        source: tokenized_tags
      			    },
                delimiter: [',', ' ']
			       });
    		} else {
    		    console.log("Handle a failed report from the CGI");
    		}
  	})
  	.fail(function() {
  		  console.log("Something broke. Unable to generate tag list");
  	});
} //end get_tag_list()

$('#btn_submit_comment').click(function(e) {
    //remove any existing warnings
	  $('.label-warning').remove();

	  e.preventDefault();

  	// check if required fields were filled in
  	var pass = check_required_fields();

  	if ( pass == true ) {
    		//submit the upload form
    		var formData = $('#upload_form').serializeArray();
    		$.ajax({
      			url: './cgi/create_github_issue.cgi',
      			type: "POST",
      			data: formData,
      			dataType: "json",
      			success: function(data, textStatus, jqXHR) {
        				if (data['success'] == 1) {
        					//activate a 'Thank you' modal
                            $('#successModal').modal('show');
        				    // redirect to the data manager page
        				    setTimeout(function() {
        				    	window.location.replace("index.html");
        				    }, 2000);

        				} else {
        				    console.log("Error during creation: " + data["error"]);
        				}
      			},
      			error: function(jqXHR, textStatus, errorThrown) {
        				var msg = 'Something went wrong. Please try again, but if this continues please contact jorvis@gmail.com for help.';
        				var error = jqXHR.status + ' ' + errorThrown;
        				display_error_bar(msg, error);
      			}
            }); //end ajax
  	} else {
  		  console.log('Some required fields were left blank.');
  	}
}); //end #btn_submit_commit

//remove warning label on focus
$('input, textarea').focus(function() {
  	var warning_el = '#label-warning-' + $(this).attr('id');
  	$(warning_el).remove();
});


//checks required fields for input.
function check_required_fields() {
    var pass = true;

  	//test whether required fields were filled in
  	if ( !$('#submitter_firstname').val() ) {
    		$('#label_submitter_firstname').append(' <span id="label-warning-submitter_firstname" class="label label-warning">Oops. This is required.</span>');
    		$('#submitter_firstname').effect('highlight', 'slow');
    		pass=false;
  	}
  	if ( !$('#submitter_lastname').val() ) {
    		$('#label_submitter_lastname').append(' <span id="label-warning-submitter_lastname" class="label label-warning">Oops. This is required.</span>');
    		$('#submitter_lastname').effect('highlight', 'slow');
    		pass=false;
  	}
  	if ( !$('#submitter_email').val() ) {
    		$('#label_submitter_email').append(' <span id="label-warning-submitter_email" class="label label-warning">Oops. This is required.</span>');
    		$('#submitter_email').effect('highlight', 'slow');
    		pass=false;
    }

  	if ( !$('#comment_title').val() ) {
    		$('#label_comment_title').append(' <span id="label-warning-comment_title" class="label label-warning">Oops. This is required.</span>');
    		$('#comment_title').effect('highlight', 'slow');
    		pass=false;
  	}
  	if ( !$('#comment').val() ) {
    		$('#label_comment').append(' <span id="label-warning-comment" class="label label-warning">Oops. This is required.</span>');
    		$('#comment').effect('highlight', 'slow');
    		pass=false;
  	}
  	if ( $('#super_impressive_security_check').val() != 28) {
    		$('#label_super_impressive_security_check').append(' <span id="label-warning-super_impressive_security_check" class="label label-warning">Oops. This is required.</span>');
    		$('#super_impressive_security_check').effect('highlight', 'slow');
    		pass=false;
  	}

  	return pass;
} //end check_required_fields

// closed alert bar
$( ".alert-container" ).on("click", "button.close-alert", function() {
	//console.log('CLOSING ALERT!');
	$( ".alert-container" ).hide();
});

// error should be html message for user. Example: error = '<p>You cannot do that.</p>'
function display_error_bar(msg, error) {
    $('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
      '<button type="button" class="close close-alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
      '<p class="alert-message">' +
      '<strong>Oops. </strong> ' +
      msg +
      '</p>' +
      '<p style="text-align: right;">(<em>Error: ' + error + '</em>)</p>' +
      '</div>').show();
}
