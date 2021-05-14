/*
 This script relies on the source having also included the
 common.js within this project
*/

var dataset_uid = null;
var share_uid = null;
var stepper = null;

var UPLOAD_STEP = 0; //tracks which step user is on

window.onload=function() {
    // check if the user is already logged in
    check_for_login();

    // Generate the UID that will be used for this submission
    dataset_uid = guid('long');
    share_uid = guid('short');

    stepper = new Stepper($('.bs-stepper')[0]);

    // Set the curator link href for when the user gets that far
    var curator_url = `./dataset_curator.html#/dataset/${dataset_uid}/displays`;
    $(".curator_link").attr("href", curator_url);
    $("#goto_curator").click(function(){
        window.location.replace(curator_url);
    });

};

/*jslint unparam: true */
/*global window, $ */
$(function () {
    'use strict';

    var url = 'uploads/';
    $('#metadata_upload').fileupload({
        url: url,
        dataType: 'json',
        done: function (e, data) {
            var filename = null
            $.each(data.result.files, function (index, file) {
                filename = file.name
            });

            validate_metadata(filename);

        },
        progressall: function (e, data) {
            $('#metadata_upload_c').hide('slide', {direction: 'up'}, 400, function(){
                $('#metadata_uploading').show();
                $('#metadata_processing_c').show('slide', {direction: 'down'}, 400);
            });

        }
    }).prop('disabled', !$.support.fileInput)
        .parent().addClass($.support.fileInput ? undefined : 'disabled');

    $('#expression_upload').fileupload({
        url: url,
        dataType: 'json',
        done: function (e, data) {
            var filename = null
            $.each(data.result.files, function (index, file) {
                filename = file.name
            });

            validate_expression(filename);
            // $('#expression_processing_c').hide();
        },
        progress: function (e, data) {
            // disable controls while loading file
            $('#goback_spreadsheet_processing, #continue_spreadsheet_processing').prop('disabled', true);

            // hide file input; show spinner
            $('#expression_upload_c').hide('slide', {direction: 'up'}, 400, function() {
                $('#expression_uploading').show();
                $('#expression_processing_c').show('slide', {direction: 'down'}, 400);
            });
        }
    }).prop('disabled', !$.support.fileInput)
        .parent().addClass($.support.fileInput ? undefined : 'disabled');
});


function validate_metadata(filename) {
    $.ajax({
        url: "./cgi/validate_metadata.cgi",
        type: "POST",
        data: {"filename": filename, "session_id": CURRENT_USER.session_id},
        dataType: "json",
        success: function(data, textStatus, jqXHR){
            if (data['success'] == 1) {
                $('#uploaded_metadata_name').val(filename);

                var disable_continue = false;

                // Build the table containing metadata and validation messages
                metadata_html = '<thead><tr><th>field</th><th>value</th><th>message</th></tr></thead>';
                metadata_html += '<tbody>';

                $.each(JSON.parse(data['metadata']), function(index, row){
                    var row_styling = '';
                    var value = row['value'];
                    var message = row['message'];
                    var is_required = row['is_required'];

                    if (!value) {
                        value = '';
                        if (is_required) {
                            row_styling = 'danger';
                            disable_continue = true;
                        }
                    }
                    
                    if (message) {
                        row_styling = 'danger';
                        disable_continue = true;
                    } else {
                        message = ''
                    }

                    metadata_html += '<tr class="'+row_styling+'">'+
                        '<td class="field">'+index+'</td>'+
                        '<td class="value" data-field="'+index+'">'+ value +'</td>'+
                        '<td class="message" data-field="'+index+'">'+ message+'</td>'+
                        '</tr>';
                });
                metadata_html += '</tbody>'
                $('#metadata_preview').html(metadata_html);

                $('#metadata_processing_c').hide('slide', {direction: 'up'}, 400, function(){
                    $('#metadata_preview_c').show('slide', {direction: 'up'}, 400);
                    $('#continue_c, #goback_c').show();
                    $('#upload_controls_c').show('slide', {direction: 'up'}, 400);
                });

                // Disable continue btn if required fields are empty
                if (disable_continue == true) {
                    $('#goback_spreadsheet_processing').prop('disabled', false);
                    $('#continue_spreadsheet_processing').prop('disabled', true);
                } else {
                    $('#goback_spreadsheet_processing, #continue_spreadsheet_processing').prop('disabled', false);
                }
                UPLOAD_STEP += 1;
            } else {
                alert_user(data['message']);

                //reset if error
                $('#metadata_processing_c').hide('slide', {direction: 'down'}, 400, function() {
                    $('#metadata_upload_c, #metadata_upload_controller').show('slide', {direction: 'up'}, 400);
                });
            }
        },
        error: function (jqXHR, textStatus, errorThrown) {
            console.log("jqXHR: ", jqXHR);
            console.log("textStatus: ", textStatus);
            console.log("errorThrown: ", errorThrown);
            var msg = "Unable to open file for validation. Please try uploading the file again.";
            alert_user(msg);

            //reset if error
            $('#metadata_processing_c').hide('slide', {direction: 'down'}, 400, function() {
                $('#metadata_upload_c, #metadata_upload_controller').show('slide', {direction: 'up'}, 400);
            });
        }
    });
}

function validate_expression(filename) {
    $.ajax({
        url: "./cgi/validate_expression.cgi", //TODO write this file
        type: "POST",
        data: {"filename": filename, "session_id": CURRENT_USER.session_id},
        dataType: "json",
        success: function(data, textStatus, jqXHR) {
            if (data['success'] == 1) {
                $('#uploaded_expression_name').val(filename);

                // Display shapes of user's data
                var shape_html = '';
                shape_html += '<p style="padding-left:20px;">X: ' + data['shape']['X']['cols'] + ' genes by ' +  data['shape']['X']['rows'] + ' observations</p>';
                shape_html += '<p style="padding-left:20px;">obs: ' + data['shape']['obs']['cols'] + ' observations by ' + data['shape']['obs']['rows'] + ' characteristic(s)</p>';
                shape_html += '<p style="padding-left:20px;">var: ' + data['shape']['var']['cols'] + ' genes by ' + data['shape']['var']['rows'] + ' characteristic(s)</p>';

                $('#expression_preview_shape').html(shape_html);
                $('#goback_spreadsheet_processing, #continue_spreadsheet_processing').prop('disabled', false);

                $('#expression_processing_c').hide('slide', {direction: 'up'}, 400, function(){
                    $('#expression_preview_c').show('slide', {direction: 'up'}, 400);
                });

                UPLOAD_STEP += 1;

                $('#continue_c').hide();
                $('#submit_upload_c').show();
            } else {
                alert_user(data['message']);

                // reset if error
                $('#expression_processing_c').hide('slide', {direction: 'down'}, 400, function() {
                    $('#expression_upload_c, #expression_upload_controller').show('slide', {direction: 'up'}, 400);
                });

                // this isn't working either
                $('#expression_processing_c').attr('style','display:none !important');

                $('.hide_after_error').hide();

                if ($('#expression_processing_c').length ) {
                    console.log("Selector worked");
                } else {
                    console.log("Selector didn't work");
                }
            }
        },
        error: function (jqXHR, textStatus, errorThrown) {
            console.log("jqXHR: ", jqXHR);
            console.log("textStatus: ", textStatus);
            console.log("errorThrown: ", errorThrown);

            var msg = "Unable to open file for validation. Please try uploading the file again.";
            alert_user(msg);

            //reset if error
            $('#expression_processing_c').hide('slide', {direction: 'down'}, 400, function() {
                $('#expression_upload_c, #expression_upload_controller').show('slide', {direction: 'up'}, 400);
            });
        }
    });
}

$(document).on('click', '.dataset_cc', function(){
    $('.dataset_cc').removeClass('selected');
    $(this).addClass('selected');

    // Store user's default_plot_type
    $('#default_plot_type').val( $(this).data('plot-type') );

    // Enable continut btn
    $('#continue_spreadsheet_processing').prop('disabled', false);
});


$("#goback_spreadsheet_processing").click(function(){
    UPLOAD_STEP -= 1;

    // Hide the preview table (then empty it). Show the metadata file input
    if (UPLOAD_STEP <= 0) {
        $('#metadata_preview_c').hide('slide', {direction: 'up'}, 400, function(){
            $('#metadata_preview').empty();
            $('#metadata_upload_c, #metadata_upload_controller').show('slide', {direction: 'up'}, 400);
        });

        $('#continue_spreadsheet_processing').prop('disabled', true);
        $('#upload_controls_c').hide();

        // Reset to zero (for the click-happy users)
        UPLOAD_STEP = 0;
    }

    // Hide expression file input. Show metadata table
    if (UPLOAD_STEP == 1 ) {
        $('#continue_spreadsheet_processing').prop('disabled', false);
        stepper.previous();
    }

    // Hide expression preview. Show expression file input
    if (UPLOAD_STEP == 2 ) {
        $('#continue_spreadsheet_processing').prop('disabled', true);

        $("#expression_preview_c").hide('slide', {direction: 'up'}, 400, function(){
            //reset this if was hidden by MEX upload. If user uploads a different format, this needs to be showing
            $("#expression_upload_c").show('slide', {direction: 'up'}, 400);
        });
    }

    // Hide schematic file input. Show expression preview
    if (UPLOAD_STEP == 3 ) {
        // $('#continue_spreadsheet_processing').prop('disabled', false);
        $('#submit_upload_c').hide();
        $('#continue_c').show();
    }

});

$("#continue_spreadsheet_processing").click(function(){
    UPLOAD_STEP += 1;
    stepper.next();

    if (UPLOAD_STEP == 2) {
        //Disable continue so user will upload expression file
        $('#continue_spreadsheet_processing').prop('disabled', true);

        // Enable 'go back' to metadata_preview
        $('#goback_spreadsheet_processing').prop('disabled', false);
    }
});

$("#submit_upload").click(function(e) {
	// remove any warnings, if present
	$('.label-warning').remove();

    // disable the submit button and put the UID into the form
    $("#submit_upload").text("Processing");
    $("#submit_upload").attr("disabled", true);
    $('#dataset_uid').val(dataset_uid);
    $('#share_uid').val(share_uid);
    $('#session_id').val(CURRENT_USER.session_id);

    e.preventDefault();

	// hide the submit button and display 'please wait' spinner
	$('#upload_controls_c').hide();
    $('#expression_preview_c').hide();
    $('#expression_uploading').text('Finalizing dataset upload ... please wait ...');
    $('#expression_processing_c').show('slide', {direction: 'down'}, 400);

    var formData = $("#upload_form").serializeArray();
    $.ajax({
        url: './cgi/load_dataset_finalize.cgi',
        type: "POST",
        data : formData,
        dataType:"json",
        success: function(data, textStatus, jqXHR) {
            if (data['success'] == 1) {
                setTimeout(function(){
                    stepper.next();
                }, 2000);
            } else {
                var msg = "An upload error has occurred. Please try again.";
                alert_user(msg);

                //reset submit button area
                $("#submit_upload").text("Submit");
                $("#submit_upload").attr("disabled", false);
                $('#submit_upload_c').show();
            }

        },
        error: function (jqXHR, textStatus, errorThrown) {
            console.log('jqXHR: ', jqXHR);
            console.log('textStatus: ', textStatus);
            console.log('errorThrown: ', errorThrown);

            var msg = "Unable to finish dataset upload. Please try again.";
            alert_user(msg);
        }
    }); //end ajax
});

function alert_user(message){
    $('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
      '<button type="button" class="close close-alert" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
      '<p class="alert-message"><strong>Oops! </strong> ' + message + '</p></div>').show();
}


/**
 * Generates a GUID string.
 * @returns {String} The generated GUID.
 * @example af8a8416-6e18-a307-bd9c-f2c947bbb3aa
 * @author Slavik Meltser (slavik@meltser.info).
 * @link http://slavik.meltser.info/?p=142
 */
function guid(uid_length) {
    function _p8(s) {
        var p = (Math.random().toString(16)+"000000000").substr(2,8);
        return s ? "-" + p.substr(0,4) + "-" + p.substr(4,4) : p ;
    }
    if (uid_length == 'long') {
        return _p8() + _p8(true) + _p8(true) + _p8();
    }
    if (uid_length == 'short') {
        return _p8();
    }
}
