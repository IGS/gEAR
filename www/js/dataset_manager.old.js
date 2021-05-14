// this helps track datasets by layout so we don't have to keep searching the database
layout_datasets = {}
current_primary_layout = null;
current_primary_label = null;
var share_id = null;
var count = 0;
var listed_datasets = [];
var listed_owners = [];

window.onload=function() {
    // check if the user is already logged in
    check_for_login();

    // keep the add layout button disabled unless a full new name has been entered
    $(document).on('keyup', '#new_layout_name', function() {
        console.log("Called");
        if ($('#new_layout_name').val().length == 0) {
            $('#confirm_layout_add').prop('disabled', true);
        } else {
            $('#confirm_layout_add').prop('disabled', false);
        }
    });
};  // end window onload

function validate_share_id(share_id, scope) {
    $.ajax({
        url : './cgi/validate_share_id.cgi',
        type: "POST",
        data : { 'session_id': session_id, 'share_id': share_id, 'scope': scope },
        dataType:"json",
        success: function(data, textStatus, jqXHR) {
            if ( data['success'] == 1 ) {
                // // toggle panel visibility
                $('#main_panel').hide();

                if (scope== 'dataset') {
                    $('#dataset_sharing_box').show();
                } else {
                    $('#profile_sharing_box').show();
                }
            } else {
                $('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
                    '<button type="button" class="close close-alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
                    '<p class="alert-message"><strong>Oops! </strong> ' + data["error"] + '</p></div>').show();
            }
        },
        error: function (jqXHR, textStatus, errorThrown) {
            console.log('textStatus= ', textStatus);
            console.log('errorThrown= ', errorThrown);
            display_error_bar(jqXHR.status + ' ' + errorThrown.name);
        }
    });
}

$('#add_shared_profile').click(function(e) {
// function add_to_profile(is_shared, dataset_id) {
    // is_shared = 0 if dataset has NOT been shared with user (is owned by user or is a public dataset)
    // is_shared = 1 if dataset has been shared with user
    // dataset_id = shared_id IF SHARED
    session_id = Cookies.get('gear_session_id');
    $.ajax({
        url : './cgi/add_shared_profile.cgi',
        type: "POST",
        data : { 'session_id': session_id, 'layout_share_id': layout_share_id },
        dataType:"json",
        success: function(data, textStatus, jqXHR) {
            if ( data['success'] == 1 ) {
                //changes URL and reloads page so shared profile gets loaded
                window.location.href = "./dataset_manager.html";
                //check_layout_control_buttons();
                //load_dataset_list_view();
            } else {
                $('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
                    '<button type="button" class="close close-alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
                    '<p class="alert-message"><strong>Oops! </strong> ' + data["error"] + '</p></div>').show();

                $('#profile_sharing_box').hide();
                $('#main_panel').show();
            }

            // $('#profile_sharing_box').hide();
            // $('#main_panel').show();
        },
        error: function (jqXHR, textStatus, errorThrown) {
            console.log('textStatus= ', textStatus);
            console.log('errorThrown= ', errorThrown);
            display_error_bar(jqXHR.status + ' ' + errorThrown.name);
        }
    }); //end ajax
});//end of #add_shared_profile

// Public/User owned dataset. Add it to profile only
$(document).on('click', 'button.add2profile', function() {
    is_shared = 0;
    dataset_id = $(this).attr('value');
    add_to_profile(is_shared, dataset_id);
});

// Shared dataset. Add it to profile & the user's shared_datasets
$(document).on('click', 'button#add_share_to_profile', function() {
    $('.alert-container').hide();

    //Warn user that datasets cannot be added to 'Site Default'
    if ( $('#selected_layout').val() == '0' || $('#selected_layout').val() == 'Site default' ) {
        var msg = 'It looks like you were trying to add a dataset to the <strong>Site Default</strong> profile.<br />' +
                  'Please <strong>choose another</strong> or <strong>add a new</strong> profile.';

        $('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
          '<button type="button" class="close close-alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
          '<p class="alert-message"><strong>Oops! </strong> ' + msg + '</p></div>').show();
    } else {
        is_shared = 1;
        add_to_profile(is_shared, share_id);
    }
});

function add_to_profile(is_shared, dataset_id) {
    // is_shared = 0 if dataset has NOT been shared with user (is owned by user or is a public dataset)
    // is_shared = 1 if dataset has been shared with user
    // dataset_id = shared_id IF SHARED
    session_id = Cookies.get('gear_session_id');

    $.ajax({
        url : './cgi/add_dataset_to_layout.cgi',
        type: "POST",
        data : { 'session_id': session_id, 'layout_id': $('#selected_layout').val(), 'dataset_id': dataset_id, is_shared: is_shared },
        dataType:"json",
        success: function(data, textStatus, jqXHR) {
            if ( data['success'] ) {
                check_layout_control_buttons();
                load_dataset_list_view();

            }

            if ( data['error'] ) {
                $('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
                    '<button type="button" class="close close-alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
                    '<p class="alert-message"><strong>Oops! </strong> ' + data["error"] + '</p></div>').show();
            }

            share_id = null;
            $('#dataset_sharing_box').hide();
            $('#main_panel').show();
        },
        error: function (jqXHR, textStatus, errorThrown) {
            console.log('textStatus= ', textStatus);
            console.log('errorThrown= ', errorThrown);
            display_error_bar(jqXHR.status + ' ' + errorThrown.name);
        }
    }); //end ajax for .add_share_to_profile
}

$(document).on('click', 'button#btn_save_layout', function() {
    layout_id =  $('#selected_layout').val()
    session_id = Cookies.get('gear_session_id');
    dataset_id_string = "";
    dataset_widths_string = "";

    $.each( $(".sortable_tile"), function(i, item) {
        dataset_id_string += $(item).data('id') + ',';
        dataset_widths_string += $(item).css('width').replace('px', '') + ',';
    });

    $.ajax({
        url : './cgi/save_layout_arrangement.cgi',
        type: "POST",
        data : { 'session_id': session_id, 'layout_id': layout_id, 'dataset_ids': dataset_id_string, 'dataset_widths': dataset_widths_string },
        dataType:"json",
        success: function(data, textStatus, jqXHR) {
            if ( data['success'] ) {
                // This could be handled better.
                location.reload();
                //check_layout_control_buttons();
                //load_dataset_list_view();
            }

            if ( data['error'] ) {
                $('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
                    '<button type="button" class="close close-alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
                    '<p class="alert-message"><strong>Oops! </strong> ' + data["error"] + '</p></div>').show();
            }
        },
        error: function (jqXHR, textStatus, errorThrown) {
            //alert("Failure!  Status: (" + textStatus + ") Error: (" + errorThrown + ")");
            console.log('textStatus= ', textStatus);
            console.log('errorThrown= ', errorThrown);
            display_error_bar(jqXHR.status + ' ' + errorThrown.name);
        }
    }); //end ajax for .add_share_to_profile
});

$( "#btn_arrangement_view" ).click(function(e) {
    $("#dataset_list_c").hide();
    $("#dataset_arrangement_c").show();

	// disable 'site default' profile so users cannot modify the arrangement
	$('option[value="0"]').prop('disabled', true);

});

if ( $('#selected_layout').val() == '0' || $('#selected_layout').val() == 'Site default' ) {
	$('#btn_arrangement_view').prop('disabled', true);
	$('#btn_arrangement_view').prop("title", "Unable to modify this profile");
	$('#btn_arrangement_view').tooltip();
}

$( "#btn_list_view" ).click(function(e) {
    $("#dataset_arrangement_c").hide();
    $("#dataset_list_c").show();

    // enable 'site default' so users can view it and set as default
	$('option[value="0"]').prop('disabled', false);
});


//Search through or List all 'Public Datasets'
$( "#btn_search_others_datasets, #btn_list_others_datasets" ).click(function(e) {
    var search_scope = 'others';
    var search_terms = $('#dataset_search_string_others').val();
    search_datasets(search_scope, search_terms);
});
// Search 'Public datasets' on 'ENTER' kepress
$('input#dataset_search_string_others').keypress(function(e) {
    if (e.keyCode == 13) {
        $('#btn_search_others_datasets').click();
    }
});

//Search through or List all 'Your datasets'
$( "#btn_search_users_datasets, #btn_list_users_datasets" ).click(function(e) {
    var search_scope = 'self';
    var search_terms = $('#dataset_search_string_users').val();
    search_datasets(search_scope, search_terms);
});
// Search "Your datasets" on 'ENTER' kepress
$('input#dataset_search_string_users').keypress(function(e) {
    if(e.keyCode == 13) {
        $('#btn_search_users_datasets').click();
    }
});

//List all datasets 'Shared with me'
$('#btn_list_shared_with_me').click(function(e) {
    var search_scope = 'shared';
    var search_terms = '';
    search_datasets(search_scope, search_terms);
});


function search_datasets(search_scope, search_terms) {
    //searches the datasets that are available to the user
    //search_scope = 'self', 'others', or 'shared'
    //search_terms = '' - list all within the scope
    //               'coch' - list datasets containing word fragment in title or long description
    $.ajax({
        url : './cgi/get_dataset_list.cgi',
        type: "POST",
        data : { 'session_id': session_id, 'scope': search_scope, 'search_terms': search_terms },
        dataType:"json",
        success: function(data, textStatus, jqXHR) {
            data_len = data['datasets'].length;

            for (i = 0; i < data_len; i++) {
                // is the profile a domain one?
                if ( $('#selected_layout').find(':selected').data('is-domain') == 1 ) {
                    data['datasets'][i]['add2profile_enabled'] = 'disabled';
                    data['datasets'][i]['removefromprofile_enabled'] = 'disabled';
                } else {
    		        // prevent adding a dataset if its already in the profile
    		        if ( $.inArray(data['datasets'][i]['dataset_id'], listed_datasets) > -1 ) {
                        data['datasets'][i]['add2profile_enabled'] = 'disabled';
                        data['datasets'][i]['removefromprofile_enabled'] = '';
    		        } else {
                        data['datasets'][i]['add2profile_enabled'] = '';
                        data['datasets'][i]['removefromprofile_enabled'] = 'disabled';
                    }
                }

                // makes fields editable if user owns dataset
                if ( data['datasets'][i]['user_name'] == CURRENT_USER.user_name ) {
                    data['datasets'][i]['user_owns'] = true;
                } else {
                    data['datasets'][i]['user_owns'] = false;
                }

                // control elements if the user owns the dataset
                if ( data['datasets'][i]['user_owns'] == true ) {
                    data['datasets'][i]['changeaccess_enabled'] = '';
                    data['datasets'][i]['delete_enabled'] = '';
                    data['datasets'][i]['shared_enabled'] = '';
                } else {
                    data['datasets'][i]['changeaccess_enabled'] = 'disabled';
                    data['datasets'][i]['delete_enabled'] = 'disabled';
                    data['datasets'][i]['share_enabled'] = 'disabled';
                }

                // make date_added nicer looking
                data['datasets'][i]['date_formatted'] = new Date(data['datasets'][i]['date_added']);

                // some fields should never have whitespace
                if (data['datasets'][i]['pubmed_id']) {
                    data['datasets'][i]['pubmed_id'] = data['datasets'][i]['pubmed_id'].trim();
                }

                if (data['datasets'][i]['geo_id']) {
                    data['datasets'][i]['geo_id'] = data['datasets'][i]['geo_id'].trim();
                }
            }//end for

            // For the list view
            var listViewTmpl = $.templates("#dataset_list_view_tmpl");
            var listViewHtml = listViewTmpl.render(data['datasets']);
            $("#dataset_list_c").html(listViewHtml);

            // For the rearrangement view
            var arrangementViewTmpl = $.templates("#dataset_arrangement_view_tmpl");
            var arrangementViewHtml = arrangementViewTmpl.render(data['datasets']);
            $("#dataset_arrangement").html(arrangementViewHtml);

            if (data_len.length == 0) {
                $('#main_panel_header').text('No results found');
            } else {
                $('#main_panel_header').text('Search results (' + data_len + ')');
            }

            update_dataset();

            // Apply any checked filters
            $('.filterby').trigger('change');
        },
        error: function (jqXHR, textStatus, errorThrown) {
            //alert("Failure!  Status: (" + textStatus + ") Error: (" + errorThrown + ")");
	          console.log('textStatus= ', textStatus);
	          console.log('errorThrown= ', errorThrown);
            display_error_bar(jqXHR.status + ' ' + errorThrown.name);
        }
    }); //end ajax for search
}// end function search_datasets

$( "#btn_set_primary_layout" ).click(function(e) {
    $.ajax({
        url : './cgi/set_primary_layout.cgi',
        type: "POST",
        data : { 'session_id': session_id, 'layout_id': $('#selected_layout').val() },
        dataType:"json",
        success: function(data, textStatus, jqXHR) {
            // Set label as cookie so it is preserved across pages
            current_primary_label = data['label'];
            Cookies.set('gear_default_domain', current_primary_label);

            current_primary_layout = $('#selected_layout').val();
            $('#btn_set_primary_layout').prop('disabled', true);
        },
        error: function (jqXHR, textStatus, errorThrown) {
            //alert("Failure!  Status: (" + textStatus + ") Error: (" + errorThrown + ")");
		    console.log('textStatus= ', textStatus);
		    console.log('errorThrown= ', errorThrown);
	        display_error_bar(jqXHR.status + ' ' + errorThrown.name);
        }
    }); //end ajax for search
});

$('#selected_layout').on('change', function(){
    check_layout_control_buttons();
    load_dataset_list_view();
});

// setup alert-container to give 'hide' function to alert button
// alert doesn't exist at window load therefore it needs to inherit from the parent div
$( ".alert-container" ).on("click", "button.close-alert", function() {
	 $( ".alert-container" ).hide();
});

function check_layout_control_buttons() {
    if ($('#selected_layout').val() == current_primary_layout) {
        $('#btn_set_primary_layout').prop('disabled', true);
    } else {
        $('#btn_set_primary_layout').prop('disabled', false);
    }

    // if selected profile == site default --> disable arrangement button
    if ( $('#selected_layout').val() == '0' || $('#selected_layout').val() == 'Site default' ) {
		$('#btn_arrangement_view').prop('disabled', true);
        $('#btn_delete_layout').prop('disabled', true);
    }
    // if any other profile
    else {
		$('#btn_arrangement_view').prop('disabled', false);
        $('#btn_delete_layout').prop('disabled', false);
    }
} //end check_layout_control_buttons()

function display_error_bar(msg) {
    $('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
      '<button type="button" class="close close-alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
      '<p class="alert-message">' +
      '<strong>Fail. </strong> Sorry, something went wrong.  Please contact us with this message if you need help.' +
      '</p>' +
      '<p style="text-align: center;">(<em>Error: ' + msg + '</em>)</p>' +
      '</div>').show();
}

function load_dataset_list_view() {
    session_id = Cookies.get('gear_session_id');
    //NOTE caused alert-container to hide before use sees message(unsupported browser). Removing for now.
    // $('.alert-container').hide();

    $.ajax({
        url : './cgi/get_dataset_list.cgi',
        type: "POST",
        data : { 'session_id': session_id, 'layout_id': $('#selected_layout').val(),
                 'exclude_pending': 0 },
        dataType:"json",
        success: function(data, textStatus, jqXHR) {
            // TODO: Most of this is duplicated in search_datasets()
            listed_datasets = []; // reset list because new profile was selected
            listed_owners = [];

            data_len = data['datasets'].length;
            for (i = 0; i < data_len; i++) {
                // is the profile a domain one?
                if ( $('#selected_layout').find(':selected').data('is-domain') == 1 ) {
                    data['datasets'][i]['add2profile_enabled'] = 'disabled';
                    data['datasets'][i]['removefromprofile_enabled'] = 'disabled';
                } else {
                    // In this function, we're always seeing only the datasets in the selected profile
                    data['datasets'][i]['add2profile_enabled'] = 'disabled';
                    data['datasets'][i]['removefromprofile_enabled'] = '';
                }

                // makes fields editable if user owns dataset
                if ( data['datasets'][i]['user_name'] == CURRENT_USER.user_name ) {
                	data['datasets'][i]['user_owns'] = true;
                } else {
                	data['datasets'][i]['user_owns'] = false;
                }

                // control elements if the user owns the dataset
                if ( data['datasets'][i]['user_owns'] == true ) {
                    data['datasets'][i]['changeaccess_enabled'] = '';
                    data['datasets'][i]['delete_enabled'] = '';
                    data['datasets'][i]['shared_enabled'] = '';
                } else {
                    data['datasets'][i]['changeaccess_enabled'] = 'disabled';
                    data['datasets'][i]['delete_enabled'] = 'disabled';
                    data['datasets'][i]['share_enabled'] = 'disabled';
                }

                // set grid width styling for arrangement view
                if (data['datasets'][i]['grid_width'] == 12) {
                    data['datasets'][i]['style_width'] = '588px';
                } else if (data['datasets'][i]['grid_width'] == 8) {
                    data['datasets'][i]['style_width'] = '392px';
                } else { // grid_width == 4
                    data['datasets'][i]['style_width'] = '196px';
                }

                //make date_added nicer looking
                data['datasets'][i]['date_formatted'] = new Date(data['datasets'][i]['date_added']);

                //get datasets that are in the layout
                listed_datasets.push(data['datasets'][i]['dataset_id']);

                //get owners of each dataset
                listed_owners.push({ 'owner': data['datasets'][i]['user_name'], 'access': data['datasets'][i]['access'] } );
            }

            // For the list view
            var listViewTmpl = $.templates("#dataset_list_view_tmpl");
            var listViewHtml = listViewTmpl.render(data['datasets']);
            $("#dataset_list_c").html(listViewHtml);

            // For the rearrangement view
            var arrangementViewTmpl = $.templates("#dataset_arrangement_view_tmpl");
            var arrangementViewHtml = arrangementViewTmpl.render(data['datasets']);
            $("#dataset_arrangement").html(arrangementViewHtml);

            $("#dataset_arrangement").sortable({
                items: "> div",
                placeholder: "ui-state-highlight"
            });
            $("#dataset_arrangement").disableSelection();

            $('.sortable_tile').resizable({
                 handles: 'e',
                 maxWidth: 600,
                 grid: 200
            });
            $('.ui-resizable-handle').css('width', '9px');

            if (data_len == 0) {
                $('#main_panel_header').text('This profile doesn\'t have any datasets yet.  Use the search options to add some.');
            } else {
			    if ($('#selected_layout').val() == '0' || $('#selected_layout').val() == 'Site default') {
                    $('#main_panel_header').text('Profile has ' + data_len + ' datasets. It cannot be modified.');
		        } else {
                	$('#main_panel_header').text('Profile has ' + data_len + ' datasets.  Use the search options to add more.');
                }
            }

            // Shut some things down if this is a domain
            if ($('#selected_layout').find(':selected').data('is-domain') == 1) {
                $(".removefromprofile").hide();
                $(".add2profile").hide();
                $('#btn_arrangement_view').prop('disabled', true);
            } else {
                $(".removefromprofile").show();
                $(".add2profile").show();
                $('#btn_arrangement_view').prop('disabled', false);
            }

            update_dataset();
            update_profile_share();
            set_schematics_maxheight();

            // Apply any checked filters
            $('.filterby').trigger('change');
        },
        error: function (jqXHR, textStatus, errorThrown) {
            //alert("Failure!  Status: (" + textStatus + ") Error: (" + errorThrown + ")");
		    console.log('textStatus= ', textStatus);
		    console.log('errorThrown= ', errorThrown);
	        display_error_bar(jqXHR.status + ' ' + errorThrown.name);
        }
    });
}

function load_user_layouts() {
    session_id = Cookies.get('gear_session_id');

    $.ajax({
        url : './cgi/get_user_layouts.cgi',
        type: "POST",
        data : { 'session_id': session_id },
        dataType:"json",
        success: function(data, textStatus, jqXHR) {
            // clear layout selector
            $('#selected_layout').empty();

            var listTmpl = $.templates("#layout_list_tmpl");
            var listHtml = listTmpl.render(data['layouts']);
            $("#selected_layout").html(listHtml);

            // TODO: actually calculate this
            $('#selected_layout').val(0);
            current_primary_layout = 0;

            check_layout_control_buttons();
            load_dataset_list_view();
        },
        error: function (jqXHR, textStatus, errorThrown) {
	          //console.log('textStatus= ', textStatus);
	          //console.log('errorThrown= ', errorThrown);
            display_error_bar(jqXHR.status + ' ' + errorThrown.name);
        }
    });
}

function process_controls() {
  $('#btn_delete_layout').popover({
		animation: true,
		trigger: 'click',
		title: "Delete profile",
		content: "<p>Are you sure you want to delete this profile?</p>" +
		"<div class='btn-toolbar' style='width:250px'>" +
		"<button class='btn btn-default btn-danger confirm_layout_delete' data-dismiss='popover'>Delete</button>" +
		"<button class='btn btn-default cancel_delete' data-dismiss='popover' value='cancel_delete'>Cancel</button>" +
		"</div>",
		html: true,
		placement: 'auto',
    container: 'body'
	}).removeClass('actions');

  $('#btn_add_layout').popover({
		animation: true,
		trigger: 'click',
		title: "Add profile",
		content: "<p>Enter a name for a profile to add</p>" +
		"<div class='btn-toolbar' style='width:250px'>" +
        "<input type='text' name='new_layout_name' id='new_layout_name' class='form-control' placeholder='Profile name' maxlength='30'>" +
		"<button id='confirm_layout_add' class='btn btn-default btn-primary confirm_layout_add' data-dismiss='popover' disabled>Add</button>" +
		"<button class='btn btn-default cancel_add' data-dismiss='popover' value='cancel_add'>Cancel</button>" +
		"</div>",
		html: true,
		placement: 'auto',
    container: 'body'
	}).removeClass('actions');

  $("#btn_share_layout").popover({
		  animation: true,
	    trigger: 'click',
	    title: "Share profile",
	    html: true,
      placement: 'right',
      container: 'body'
	}).removeClass('actions')
    .on('show.bs.popover', function(e) {
      // Helped to insert the share_url more reliably: http://stackoverflow.com/a/25885326/2900840
      var el = $(e.target);

      // Help to close other popovers: http://stackoverflow.com/a/34320956/2900840
      el.data("bs.popover")._activeTrigger.click = false;
      $('#btn_share_layout').not(el).popover('hide');

      //add warning to popover if layout contains private datasets
      var share_warning = '';
      for ( var i=0, len=listed_owners.length; i<len; i++ ) {
          if ( listed_owners[i]['access'] == 'Private' ) {
              share_warning = '<div class="alert alert-warning" role="alert">' +
                              '<b>Caution!</b> Sharing this profile will also share private datasets.' +
                              '</div>';
              break; //only need to find one
          }
      }

      //build share_url and popover html content
      var current_url = window.location.href;
      var current_page = current_url.lastIndexOf("dataset_manager.html");
      var layout_permalink = current_url.substring(0, current_page) + 'p?l=' + $('#selected_layout').find(':selected').data('share-id');
      var html =  share_warning +
              "<p>Copy the URL below to share with others</p>" +
          		"<div class='btn-toolbar'>" +
              "<input type='text' name='profile_share_url' id='profile_share_url' class='form-control' value='" + layout_permalink + "'>" +
          		"<button class='btn btn-default cancel_share pull-right' data-dismiss='popover' value='cancel_share'>Done</button>" +
          		"</div>";

      //insert content into share popover
      el.attr('data-content', html);
  });

    $(document).on('click', '.confirm_layout_delete', function() {
        $('#btn_delete_layout').popover('hide');

        session_id = Cookies.get('gear_session_id');
        $.ajax({
            url : './cgi/remove_layout.cgi',
            type: "POST",
            data : { 'session_id': session_id, 'layout_id': $('#selected_layout').val() },
            dataType:"json",
            success: function(data, textStatus, jqXHR) {
                if ( data['success'] == 1 ) {
                    old_value = $('#selected_layout').val();

                    // loop through layouts and delete this one:
                    delete layout_datasets[old_value];

                    // Set the default layout to be the current one
                    $('#selected_layout').val('0');

                    // Delete the other value
                    $('option[value="' + old_value + '"]', $('#selected_layout')).remove();

                    check_layout_control_buttons();
                    load_dataset_list_view();
                }
                // if user does NOT own dataset
                else {
                    console.log(data);
                    $('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
                        '<button type="button" class="close close-alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
                        '<p class="alert-message"><strong>Oops! </strong> ' + data["error"] + '</p></div>').show();
                }
            },
            error: function (jqXHR, textStatus, errorThrown) {
                //alert("Failure!  Status: (" + textStatus + ") Error: (" + errorThrown + ")");
                console.log('textStatus= ', textStatus);
                console.log('errorThrown= ', errorThrown);
                display_error_bar(jqXHR.status + ' ' + errorThrown.name);
            }
        }); //end ajax for .confirm_delete
    });

    $(document).on('click', '.confirm_layout_add', function() {
        $('#btn_add_layout').popover('hide');

        session_id = Cookies.get('gear_session_id');
        $.ajax({
            url : './cgi/add_layout.cgi',
            type: "POST",
            data : { 'session_id': session_id, 'layout_name': $('#new_layout_name').val() },
            dataType:"json",
            success: function(data, textStatus, jqXHR) {
                if ( data['layout_id'] ) {

                    // Add it to the select box
                    $('#selected_layout').append('<option value="' + data['layout_id'] + '">' + data['layout_label'] + '</option>');

                    layout_datasets[data['layout_id']] = [];

                    // Set the default layout to be the new one
                    $('#selected_layout').val(data['layout_id']);

                    check_layout_control_buttons();
                    load_dataset_list_view();

                }
                // if user does NOT own dataset
                if ( data['error'] ) {
                    $('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
                        '<button type="button" class="close close-alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
                        '<p class="alert-message"><strong>Oops! </strong> ' + data["error"] + '</p></div>').show();
                }
            },
            error: function (jqXHR, textStatus, errorThrown) {
                //alert("Failure!  Status: (" + textStatus + ") Error: (" + errorThrown + ")");
                console.log('textStatus= ', textStatus);
                console.log('errorThrown= ', errorThrown);
                display_error_bar(jqXHR.status + ' ' + errorThrown.name);
            }
        }); //end ajax for .confirm_delete
    });

    $("#dataset_search_string_users").keyup(function(){
        if($(this).val().length > 0) {
            // Enable Search button
            $('#btn_list_users_datasets').hide();
            $('#btn_search_users_datasets').show();

            // Disable Share with me button
            $('#btn_list_shared_with_me').prop('disabled', true);
        } else {
            // Enable the List All button
            $('#btn_search_users_datasets').hide();
            $('#btn_list_users_datasets').show();
            $('#btn_list_shared_with_me').prop('disabled', false);
        }
    });

    $("#dataset_search_string_others").keyup(function(){
        if($(this).val().length > 0) {
            // Enable Search button
            $('#btn_list_others_datasets').hide();
            $('#btn_search_others_datasets').show();
        } else {
            // Enable the List All button
            $('#btn_search_others_datasets').hide();
            $('#btn_list_others_datasets').show();
        }
    });
}

/**
 * Generates a GUID string.
 * modified to inject and mask a known id
 * @returns {String} The generated GUID.
 * @example af8a8416-6e18-a307-bd9c-f2c947bbb3aa
 * @author Slavik Meltser (slavik@meltser.info).
 * @link http://slavik.meltser.info/?p=142
 */
function guid(id) {
    id = Number(id) + 85;
    function _p8(s) {
        var p = (Math.random().toString(16)+"000000000").substr(2,8);
        return s ? "-" + p.substr(0,4) + "-" + p.substr(4,4) : p ;
    }
    return _p8() + '-' + id + _p8(true) + _p8();
}

// initiate the manage_sharing modal
$(document).on('click', '.manage_sharing', function(){
    var dataset_id = $(this).attr('data-dataset-id');
    $(this).popover('hide');

    // reset the modal
    $('#table_c').show();
    $('#table_no_shares').hide(); //hide if previously showing
    $('#shared_users_table_c').empty(); //empty if previously populated

    $.ajax({
        url : './cgi/get_shared_users_list.cgi',
        type: "POST",
        data : { 'session_id': session_id, 'dataset_id': dataset_id },
        dataType:"json",
        success: function(data, textStatus, jqXHR) {
            if ( data['success'] == 1 ) {
                var list_len = data['user_list'].length;
                for (i = 0; i < list_len; i++) {
                    //change button to 'reshare'
                    if (data['user_list'][i]['is_allowed'] == 0) {
                        data['user_list'][i]['is_allowed'] = false;
                    } else {
                        data['user_list'][i]['is_allowed'] = true;
                    }
                }

                if (list_len > 0) {
                    //build table of user_list
                    var manageSharingViewTmpl = $.templates("#manage_sharing_view_tmpl");
                    var manageSharingViewHtml = manageSharingViewTmpl.render(data['user_list']);
                    $("#shared_users_table_c").html(manageSharingViewHtml);

                    //show modal
                    $('#shared_users_modal').modal('show');
                } else {
                  // hide table and show 'dataset not shared' message
                  $('#table_c').hide();
                  $('#table_no_shares').show();
                  $('#shared_users_modal').modal('show');
                }
            }
            // error within CGI script
            else {
                $('#shared_users_modal').modal('hide');
                $('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
                    '<button type="button" class="close close-alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
                    '<p class="alert-message"><strong>Oops! </strong> ' + data["error"] + '</p></div>').show();
            }
        },
        error: function (jqXHR, textStatus, errorThrown) {
            $('#shared_users_modal').modal('hide');

            console.log('textStatus= ', textStatus);
            console.log('errorThrown= ', errorThrown);
            display_error_bar(jqXHR.status + ' ' + errorThrown.name);
        }
    }); //end ajax
});

//user is unsharing dataset
$(document).on('click', '.revoke_share', function() {
    var dataset_id = $(this).attr('data-dataset_id');
    var id = $(this).attr('id');
    var user_id= id.substring(id.lastIndexOf('_') + 1);
    var to_share = 0; //false - user is unsharing dataset

    manage_dataset_shares(user_id, dataset_id, to_share);
}); //end .revoke_share

//user is resharing dataset
$(document).on('click', '.reallow_share', function() {
    var dataset_id = $(this).attr('data-dataset_id');
    var id = $(this).attr('id');
    var user_id= id.substring(id.lastIndexOf('_') + 1);
    var to_share = 1; //true - user is resharing dataset

    manage_dataset_shares(user_id, dataset_id, to_share);
}); //end .reallow_share

// handles both unsharing & resharing events
function manage_dataset_shares(user_id, dataset_id, to_share){
    //fade user row until ajax success or error
    $('#row_' + user_id).css('opacity', '0.5');

    $.ajax({
        url : './cgi/manage_dataset_shares.cgi',
        type: "POST",
        data : { 'session_id': session_id, 'dataset_id': dataset_id, 'user_id': user_id, 'to_share': to_share },
        dataType:"json",
        success: function(data, textStatus, jqXHR) {
            if ( data['success'] == 1 ) {
                $('#row_' + user_id).css('opacity', '1.0');
                for (i = 0; i < data['user_list'].length; i++) {
                    //change button to 'reshare'
                    if (data['user_list'][i]['is_allowed'] == 0) {
                        $('button#revoke_share_' + user_id).replaceWith('<button type="button" id="reallow_share_' +
                            data["user_list"][i]["user_id"] + '" ' +
                            'data-dataset_id="' + data["user_list"][i]["dataset_id"] + '" ' +
                            'class="btn btn-success reallow_share" autocomplete="off">Reshare</button>');
                    }
                    //change button to 'revoke'
                    else {
                        $('button#reallow_share_' + user_id).replaceWith('<button type="button" id="revoke_share_' +
                            data["user_list"][i]["user_id"] + '" ' +
                            'data-dataset_id="' + data["user_list"][i]["dataset_id"] + '" ' +
                            'class="btn btn-danger revoke_share" autocomplete="off">Unshare</button>');
                    }
                }//end for loop
            }
            // error within CGI script
            else {
                $('#row_' + user_id).css('opacity', '1.0');
                $('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
                    '<button type="button" class="close close-alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
                    '<p class="alert-message"><strong>Oops! </strong> ' + data["error"] + '</p></div>').show();
            }
        },
        error: function (jqXHR, textStatus, errorThrown) {
            $('#row_' + user_id).css('opacity', '1.0');

            console.log('textStatus= ', textStatus);
            console.log('errorThrown= ', errorThrown);
            display_error_bar(jqXHR.status + ' ' + errorThrown.name);
        }
    }); //end ajax
}


// For dataset buttons: Change_access, view, share, delete
function update_dataset() {
	$('.change_access[data-toggle="popover"]').popover({
		  animation: true,
	    trigger: 'click',
	    title: "Change Access",
	    content: "<p>Change access to dataset.</p>" +
		    "<div class='btn-toolbar'>" +
	      "<button class='btn btn-default btn-primary private' data-dismiss='popover' value='private'>Private</button>" +
	      "<button class='btn btn-default btn-primary public' data-dismiss='popover' value='public'>Public</button>" +
	      "<button class='btn btn-default cancel_change_access' data-dismiss='popover' value='cancel_change_access'>Cancel</button>" +
        "</div>",
      html: true,
	    placement: 'auto',
      viewport: {selector: 'body', padding: 0}
	}).removeClass('actions');

	$('.delete[data-toggle="popover"]').popover({
		  animation: true,
		  trigger: 'click',
		  title: "Confirm Delete",
		  content: "<p>Are you sure you want to delete this dataset?</p>" +
		    "<div class='btn-toolbar'>" +
	      "<button class='btn btn-default btn-danger confirm_delete' data-dismiss='popover'>Delete</button>" +
        "<button class='btn btn-default cancel_delete' data-dismiss='popover' value='cancel_delete'>Cancel</button>" +
        "</div>",
		  html: true,
		  placement: 'auto',
		  viewport: {selector: 'body', padding: 0}
	}).removeClass('actions');

  $('.share[data-toggle="popover"]').popover({
		  animation: true,
	    trigger: 'click',
	    title: "Share URL",
	    html: true,
      placement: 'auto'
	}).removeClass('actions')
    .on('show.bs.popover', function(e) {
      // Helped to insert the share_url more reliably: http://stackoverflow.com/a/25885326/2900840
      var el = $(e.target);
      var dataset_id = $(el).attr('data-dataset-id');
      // Help to close other popovers: http://stackoverflow.com/a/34320956/2900840
      el.data("bs.popover")._activeTrigger.click = false;
      $('.share').not(el).popover('hide');

      //build share_url and popover html content
      var dataset_share_id = $(el).attr('value');
      dataset_share_url = window.location.href + '?share_id=' + dataset_share_id;
      var html = "<p>Copy the URL below to share with others</p>" +
          		"<div class='btn-toolbar'>" +
              "<input type='text' name='dataset_share_url' id='dataset_share_url' class='form-control' value='" + dataset_share_url + "'>" +
              "<button class='btn btn-primary manage_sharing' data-dataset-id='" + dataset_id + "'data-dismiss='popover' value='manage_sharing'>Manage Sharing</button>" +
          		"<button class='btn btn-default cancel_share pull-right' data-dismiss='popover' value='cancel_share'>Done</button>" +
          		"</div>";

      //insert content into share popover
      el.attr('data-content', html);
  });

  $('.permalink[data-toggle="popover"]').popover({
		  animation: true,
	    trigger: 'click',
	    title: "Share URL",
	    html: true,
      placement: 'auto'
	}).removeClass('actions')
    .on('show.bs.popover', function(e) {
      // Helped to insert the share_url more reliably: http://stackoverflow.com/a/25885326/2900840
      var el = $(e.target);

      // Help to close other popovers: http://stackoverflow.com/a/34320956/2900840
      el.data("bs.popover")._activeTrigger.click = false;
      $('.share').not(el).popover('hide');

      //build share_url and popover html content
      var dataset_share_id = $(el).attr('value');
      var current_url = window.location.href;
      var current_page = current_url.lastIndexOf("dataset_manager.html");
      var dataset_permalink = current_url.substring(0, current_page) + 'p?s=' + dataset_share_id;
      console.log(dataset_permalink);

      var html = "<p>Include this URL in publications so others can view this dataset.</p>" +
          		"<div class='btn-toolbar'>" +
              "<input type='text' name='dataset_permalink' id='dataset_permalink' class='form-control' value='" + dataset_permalink + "'>" +
          		"<button class='btn btn-default cancel_share pull-right' data-dismiss='popover' value='cancel_share'>Done</button>" +
          		"</div>";

      //insert content into share popover
      el.attr('data-content', html);
  });

	//Hides popover with cancel button: http://stackoverflow.com/a/28106695/2900840
	$('[data-toggle="popover"]').each(function () {
        var button = $(this);
		button.popover().on('shown.bs.popover', function() {
		    $(button.data('bs.popover').tip).find('[data-dismiss="popover"]').on('click', function () {
		        button.popover('hide');
		    });
		});
	});

	$( ".actions" ).on("click", "button.removefromprofile, button.change_access, button.curate, button.view, button.delete", function(e) {
      // get the dataset id from siblings since button is inserted after ajax
      var siblings_dataset_id= $(this).siblings().attr('value');
      var dataset_div = '#' + siblings_dataset_id;

      if ( $(e.target).hasClass('curate') ) {
          window.location.href = `./dataset_curator.html#/dataset/${$(this).attr('value')}/displays`;
      }

      if ( $(e.target).hasClass('removefromprofile') ) {

    			// Just in case this button were to be enabled with 'site default' profile selected
    			if ($('#selected_layout').val() == '0' || $('#selected_layout').val() == 'Site default') {

      				// Display error message. Cannot edit site default
      				$('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
      					'<button type="button" class="close close-alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
      					'<p class="alert-message"><strong>Oops! </strong>You cannot remove this dataset from the Site Default Profile.</p></div>').show();
    			} else {
    				// get layout id
    				var layout_id = $('#selected_layout').val();

    				$(dataset_div).css('opacity', '0.5'); //fade until ajax success/error
    				session_id = Cookies.get('gear_session_id'); //get user id for sql qry
    				$.ajax({
    					cache: false,
    					url: './cgi/remove_dataset_from_profile.cgi',
    					data: { 'session_id': session_id, 'dataset_id': siblings_dataset_id, 'layout_id': layout_id },
    					dataType: 'json',
    					success: function(data, textStatus, jqXHR) {
      						// if dataset was removed
      						if ( data['dataset'] ) {
        						$(dataset_div).hide('fade', '700'); //hide if success (if deleted)

                                //remove dataset from listed_datasets so it can be added back to the profile without reloading
                                if ( $.inArray(siblings_dataset_id, listed_datasets) >= -1 ) {
                                    var i = listed_datasets.indexOf(siblings_dataset_id);
                                    listed_datasets.splice(i,1);
                                }

                                //remove owner from listed_owners
                                var owner = $('#owner_name_' + siblings_dataset_id).text().replace('Owner: ', '');
                                var access = $('#access_' + siblings_dataset_id).text().replace('Access: ', '');
                                for ( var i=0, len=listed_owners.length; i<len; i++ ) {
                                    if ( listed_owners[i]['owner'] == owner && listed_owners[i]['access'] == access ) {
                                        listed_owners.splice(i,1);
                                        break; //only remove one
                                    }
                                }

                                update_profile_share();
      						}
      						// if dataset was not removed
      						if ( data['error'] ) {
  				                $(dataset_div).css('opacity', '1');
  		  			            $('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
      	  							                       '<button type="button" class="close close-alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
      	  							                       '<p class="alert-message"><strong>Oops! </strong> ' + data["error"] + '</p></div>').show();
      						}
    					},
    					error: function (jqXHR, textStatus, errorThrown) {
    			    		$(dataset_div).css('opacity', '1');
	                        display_error_bar(jqXHR.status + ' ' + errorThrown.name);
    					}
    				}); //end ajax
    			} // end else
	  } // end if removefromprofile

  	  // Change access level of dataset
      $('.private, .public').click(function() {
      		var id = dataset_div.replace("#", "");
      		var access = $(this).attr('value');
      		$(this).popover('hide');

      		$(dataset_div).css('opacity', '0.5'); //fade until ajax success/error
  		    session_id = Cookies.get('gear_session_id'); // gets the user id for the sql qry
    			$.ajax({
      				cache: false,
      				url : './cgi/change_dataset_access.cgi',
      				type: "POST",
      				data : { 'session_id': session_id, 'dataset_id': id, 'access': access },
      				dataType:"json",
      				success: function(data, textStatus, jqXHR) {

      					// if user owns dataset
      					if ( data['dataset'] ) {
  						      $('#access_' + id).text("Access: " + data['dataset'][0]['access']);
						        $(dataset_div).css('opacity', '1');
  						      $('#access_' + id).effect('highlight', 'slow');
      					}
      					// if user does NOT own dataset
      					if ( data['error'] ) {

  						      // unfade dataset and alert user access cannot be changed, they dont own the dataset
						        $(dataset_div).css('opacity', '1');
  					        $('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
        							'<button type="button" class="close close-alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
        							'<p class="alert-message"><strong>Oops! </strong> ' + data["error"] + '</p></div>').show();
      					}
      					id = null;
      				},
      				error: function (jqXHR, textStatus, errorThrown) {
              		$(dataset_div).css('opacity', '1');
                  display_error_bar(jqXHR.status + ' ' + errorThrown.name);
      				}
    			}); // end ajax
  	  }); // end access level

      $('.confirm_delete').click(function() {
      		$(this).popover('hide');
      		var id = dataset_div.replace("#", "");
      		$(dataset_div).css('opacity', '0.5'); // fade while ajax is made

      		session_id = Cookies.get('gear_session_id');
          $.ajax({
              url: './cgi/mark_dataset_for_removal.cgi',
      				type: "POST",
      				data : { 'session_id': session_id, 'dataset_id': id },
      				dataType:"json",
      				success: function(data, textStatus, jqXHR) {

        					// user owns dataset
        					if ( data['success'] == 1 ) {
        	    				$(dataset_div).hide('fade', '700'); //hide if success (if deleted)
        					}
        					// user does NOT own dataset
        					else {
    						      $(dataset_div).css('opacity', '1');
            					$('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
          							'<button type="button" class="close close-alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
          							'<p class="alert-message"><strong>Oops! </strong> ' + data["error"] + '</p></div>').show();
        					}
        					id = null;
      				},
      				error: function (jqXHR, textStatus, errorThrown) {
    	        		$(dataset_div).css('opacity', '1');
                  display_error_bar(jqXHR.status + ' ' + errorThrown.name);
              }
          }); //end ajax for .confirm_delete
          dataset_id = null;
	   }); // end .confirm_delete
	}); // end .actions btns

	// Close any open popovers and open current. URL: http://www.bootply.com/61233
	$('[data-toggle=popover]').on('click', function (e) {
   		$('[data-toggle=popover]').not(this).popover('hide');
	});
}

function update_profile_share() {
    // check if profile contains any datasets owned by someone else
    // if found, disable #btn_share_layout
    var multiple_owners = false;
    var is_private = false;
    for (var i=0, len=listed_owners.length; i<len; i++) {
        if ( listed_owners[i]['owner'] != CURRENT_USER.user_name && listed_owners[i]['access'] == 'Private' ) {
            multiple_owners = true;
            // if ( listed_owners[i]['access'] == 'Private' ) {
                // is_private = true;
                // break;
            // }
            is_private = true;
            break;
        }
    }
    if ( multiple_owners == true && is_private == true ) {
        $('#btn_share_layout').prop('disabled', true);
    } else if ( $('#selected_layout').val() == '0' || $('#selected_layout').val() == 'Site default' ) {
        $('#btn_share_layout').prop('disabled', true);
    } else {
        $('#btn_share_layout').prop('disabled', false);
    }
} //end update_profile_share

//remove warning labels in modal on focus
$('#dataset_title').focus(function(){
	$('#label-warning-title').remove();
});
$('#dataset_desc').focus(function(){
	$('#label-warning-desc').remove();
});

function check_required_fields() {
	var pass = true;

	// test whether required fields were filled in
	if ( !$('#dataset_title').val() ) {
		$('#label_title_input').append(' <span id="label-warning-title" class="label label-warning">Oops. This is required.</span>');
		$('#dataset_title').effect('highlight', 'slow');
		pass=false;
	}

    if ( !$('#dataset_desc').val() ) {
        $('#label_desc_input').append(' <span id="label-warning-desc" class="label label-warning">Oops. This is required.</span>');
        $('#dataset_desc').effect('highlight', 'slow');
        pass=false;
    }

	return pass;

}

// X-EDITABLE - TITLE
$(document).on('click', '.editable-title', function() {
  	var infoType = 'title';
  	var dataset_id = $(this).data('dataset-id');
  	var oldValue = $(this).text();
  	var session_id = Cookies.get('gear_session_id');
  	var data = {};

  	$('#title_' + dataset_id).editable({
    		mode: 'inline',
    		type: 'text',
            tpl: '<input type="text" style="width: 800px;">',
    		pk: dataset_id,
    		params: function(params) {
      			var data = {};
      			data['id'] = params.pk; //dataset_id
      			data['field'] = infoType; //for CGI to determine the field to change
      			data['new_value'] = params.value; // db_column_name = 'new value'
      			data['session_id'] = session_id;
      			return data;
    		},
    		url: 'cgi/save_datasetinfo_changes.cgi',
    		title: 'Enter Title',
    		success: function(data, newValue, oldValue) {
      			if (data['dataset']) {
        				$(this).editable('destroy');
        				$(this).text(newValue);
        				$(this).effect('highlight', 'slow');
      			}

      			// user doesn't own dataset or unexpected field was edited
      			if (data['error']) {
        				$(this).editable('destroy');
        				$(this).text(oldValue);
        				$('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
        					'<button type="button" class="close close-alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
        					'<p class="alert-message"><strong>Oops! </strong> ' + data["error"] + '</p></div>').show();
      			}
    		},
    		error: function(data, newValue, oldValue) {
      			$(this).text(oldValue);
      			if (data['error']) {
        				$(this).editable('destroy');
        				$(this).text(oldValue);
        				$('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
        					'<button type="button" class="close close-alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
        					'<p class="alert-message"><strong>Oops! </strong> ' + data["error"] + '</p></div>').show();
      			}
    		}
  	});
});

// X-EDITABLE - DESCRIPTION
$(document).on('click', '.editable-ldesc', function() {
  	var infoType = 'ldesc';
  	var dataset_id = $(this).data('dataset-id');
  	var oldValue = $(this).text();
  	var session_id = Cookies.get('gear_session_id');
  	var data = {};
  	$('#ldesc_' + dataset_id).editable({
    		mode: 'inline',
    		type: 'textarea',
            tpl: '<textarea type="text" style="width: 800px;">',
            inputclass: 'editable-textarea-block',
    		pk: dataset_id,
    		params: function(params) {
      			var data = {};
      			data['id'] = params.pk; //dataset_id
      			data['field'] = infoType; //for CGI to determine the field to change
      			data['new_value'] = params.value; // db_column_name = 'new value'
      			data['session_id'] = session_id;
      			return data;
    		},
    		url: 'cgi/save_datasetinfo_changes.cgi',
    		title: 'Enter Description',
    		success: function(data, newValue, oldValue) {
      			if (data['dataset']) {
        				$(this).editable('destroy');
        				$(this).text(newValue);
        				$(this).effect('highlight', 'slow');
      			}

      			// user doesn't own dataset or unexpected field was edited
      			if (data['error']) {
        				$(this).editable('destroy');
        				$(this).text(oldValue);
        				$('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
        					'<button type="button" class="close close-alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
        					'<p class="alert-message"><strong>Oops! </strong> ' + data["error"] + '</p></div>').show();
      			}
    		},
    		error: function(data, newValue, oldValue) {
      			$(this).text(oldValue);
      			if (data['error']) {
        				$(this).editable('destroy');
        				$(this).text(oldValue);
        				$('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
        					'<button type="button" class="close close-alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
        					'<p class="alert-message"><strong>Oops! </strong> ' + data["error"] + '</p></div>').show();
      			}
    		}
  	});
});

// X-EDITABLE - TAGS
$(document).on('click', '.editable-tag', function() {
  	var infoType = 'tag';
  	var dataset_id = $(this).data('dataset-id');
  	var oldValue = $(this).text();
  	var session_id = Cookies.get('gear_session_id');
  	var data = {};

    // get_tag_list(dataset_id);
  	$('#tag_list_' + dataset_id).editable({
    		mode: 'inline',
    		type: 'text',
            tpl: '<input type="text" id="tag_input_'+dataset_id+'" value="'+oldValue+'" style="width:600px;">',
    		pk: dataset_id,
        clear: false,
    		params: function(params) {
      			var data = {};
      			data['id'] = params.pk; //dataset_id
      			data['field'] = infoType; //for CGI to determine the field to change
      			data['new_value'] = params.value.replace(/ /g, ''); //makes comma-separated string for CGI script
            data['session_id'] = session_id;
      			return data;
    		},
    		url: 'cgi/save_datasetinfo_changes.cgi',
    		title: 'Enter Description',
    		success: function(data, newValue, oldValue) {
      			if (data['dataset']) {
                // Remove any duplicate tags
                var tag_list = [];
                var raw_tag_list = newValue.split(',');
                for ( var t=0, len=raw_tag_list.length; t<len; t++ ) {
                  if (tag_list.indexOf(raw_tag_list[t]) <= 0){
                    tag_list.push(raw_tag_list[t]);
                  }
                }
                $(this).editable('destroy');
                $(this).text(tag_list.join(', '));
        				$(this).effect('highlight', 'slow');
      			}

      			// user doesn't own dataset or unexpected field was edited
      			if (data['error']) {
        				$(this).editable('destroy');
        				$(this).text(oldValue);
        				$('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
        					'<button type="button" class="close close-alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
        					'<p class="alert-message"><strong>Oops! </strong> ' + data["error"] + '</p></div>').show();
      			}
    		},
    		error: function(data, newValue, oldValue) {
      			$(this).text(oldValue);
      			if (data['error']) {
        				$(this).editable('destroy');
        				$(this).text(oldValue);
        				$('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
        					'<button type="button" class="close close-alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
        					'<p class="alert-message"><strong>Oops! </strong> ' + data["error"] + '</p></div>').show();
      			}
    		}
  	});
    get_tag_list(dataset_id);
});

// X-EDITABLE - PUBMED ID
$(document).on('click', '.editable-pubmed', function() {
  	var infoType = 'pubmed';
  	var dataset_id = $(this).data('dataset-id');
  	var oldValue = $('#pubmed_'+dataset_id).text();
  	var session_id = Cookies.get('gear_session_id');
  	var data = {};
  	$('#pubmed_' + dataset_id).editable({
    		mode: 'inline',
    		type: 'text',
    		pk: dataset_id,
    		params: function(params) {
      			var data = {};
      			data['id'] = params.pk; //dataset_id
      			data['field'] = infoType; //for CGI to determine the field to change
      			data['new_value'] = params.value; // db_column_name = 'new value'
      			data['session_id'] = session_id;
      			//console.log("pudmed: " + JSON.stringify(data));
      			return data;
    		},
    		url: 'cgi/save_datasetinfo_changes.cgi',
    		title: 'Enter PubMed ID',
    		success: function(data, newValue, oldValue) {
      			if (data['dataset']) {
        				$(this).editable('destroy');
                $(this).text(newValue);
                $('#pubmed_url_'+dataset_id).attr('href', 'http://www.ncbi.nlm.nih.gov/pubmed/' + newValue);
                if (newValue == '') {
                    $('#pubmed_url_'+dataset_id).text('');
                } else {
                    $('#pubmed_url_'+dataset_id).text('link');
                }
        				$(this).effect('highlight', 'slow');
      			}

      			// user doesn't own dataset or unexpected field was edited
      			if (data['error']) {
        				$(this).editable('destroy');
                $('#pubmed_url_'+dataset_id).attr('href', 'http://www.ncbi.nlm.nih.gov/pubmed/' + oldValue);
                if (oldValue == '') {
                    $('#pubmed_url_'+dataset_id).text('');
                } else {
                    $('#pubmed_url_'+dataset_id).text('link');
                }
        				$('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
        					'<button type="button" class="close close-alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
        					'<p class="alert-message"><strong>Oops! </strong> ' + data["error"] + '</p></div>').show();
      			}
    		},
    		error: function(data, newValue, oldValue) {
      			$(this).text(oldValue);
        			if (data['error']) {
        				$(this).editable('destroy');
        				$(this).text(oldValue);
        				$('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
        					'<button type="button" class="close close-alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
        					'<p class="alert-message"><strong>Oops! </strong> ' + data["error"] + '</p></div>').show();
      			}
    		}
  	});
});

// X-EDITABLE - GEO ID
$(document).on('click', '.editable-geo', function() {
  	var infoType = 'geo';
  	var dataset_id = $(this).data('dataset-id');
  	var oldValue = $('#geo_'+dataset_id).text();
  	var session_id = Cookies.get('gear_session_id');
  	var data = {};
  	$('#geo_' + dataset_id).editable({
    		mode: 'inline',
    		type: 'text',
    		pk: dataset_id,
    		params: function(params) {
      			var data = {};
      			data['id'] = params.pk; //dataset_id
      			data['field'] = infoType; //for CGI to determine the field to change
      			data['new_value'] = params.value; // db_column_name = 'new value'
      			data['session_id'] = session_id;
      			//console.log("pudmed: " + JSON.stringify(data));
      			return data;
    		},
    		url: 'cgi/save_datasetinfo_changes.cgi',
    		title: 'Enter GEO ID',
    		success: function(data, newValue, oldValue) {
      			if (data['dataset']) {
        				$(this).editable('destroy');
                $(this).text(newValue);
                $('#geo_url_'+dataset_id).attr('href','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=' + newValue);
                if (newValue == '') {
                    $('#geo_url_'+dataset_id).text('');
                } else {
                    $('#geo_url_'+dataset_id).text('link');
                }
        				$(this).effect('highlight', 'slow');
      			}

      			// user doesn't own dataset or unexpected field was edited
      			if (data['error']) {
        				$(this).editable('destroy');
        				$('#geo_url_'+dataset_id).attr('href','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=' + oldValue);
                if (oldValue == '') {
                    $('#geo_url_'+dataset_id).text('');
                } else {
                    $('#geo_url_'+dataset_id).text('link');
                }
        				$('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
        					'<button type="button" class="close close-alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
        					'<p class="alert-message"><strong>Oops! </strong> ' + data["error"] + '</p></div>').show();
      			}
    		},
    		error: function(data, newValue, oldValue) {
      			$(this).text(oldValue);
      			if (data['error']) {
        				$(this).editable('destroy');
        				$(this).text(oldValue);
        				$('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
        					'<button type="button" class="close close-alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
        					'<p class="alert-message"><strong>Oops! </strong> ' + data["error"] + '</p></div>').show();
      			}
    		}
  	});
});

// X-EDITABLE - SCHEMATIC IMAGE
$(document).on('click', '.editable-schematic', function() {
  	var infoType = 'schematic';
  	var dataset_id = $(this).data('dataset-id');
  	var oldValue = $(this).attr('src');
  	var session_id = Cookies.get('gear_session_id');
  	var data = {};

  	//move file to 'uploads'
  	$('#schematic_upload').fileupload({
    		url: 'uploads/',
    		dataType: 'json',
    		done: function (e, data) {
            var filename = null
            $.each(data.result.files, function (index, file) {
                filename = file.name
            });

            //remove preexisting thumbnail and filename
            $('#thumbnail_c').remove();
            $('#schematic_upload').attr('value', '');

            //apply new file to input (for x-editable)
            $('#schematic_upload').attr('value', filename);

            //show thumbnail of new schematic
            var schematic_filepath = 'uploads/files/' + filename;
            $('#schematic_upload').after('<div id="thumbnail_c" class="row">' +
            		'<img id="schematic_thumbnail" src="' + schematic_filepath +
            		'" style="height:120px; margin-top: 7px;" alt="uploaded image"/></div>');
    		}
  	});

	$('#schematic_' + dataset_id).editable({
  		mode: 'inline',
  		type: 'text',
  		tpl: '<input id="schematic_upload" type="file" name="files[]">',
  		showbuttons: 'bottom',
  		clear: false,
  		pk: dataset_id,
  		params: function(params) {
    			var data = {};
    			data['id'] = params.pk; //dataset_id
    			data['field'] = infoType; //for CGI to determine the field to change
    			data['session_id'] = session_id;

    			// set params.value to latest filename selected
    			params.value = $('#schematic_upload').attr('value');

    			// For Chrome (was inserting '\fakepath\' before actual filename)
    			if ( params.value.indexOf('\\fakepath\\') >= 0) {
      				var new_filename = params.value;
      				new_filename = new_filename.split("\\fakepath\\")[1]; //remove Chrome inserted substring
      				data['new_value'] = new_filename; // db_column_name = 'new value'
      				return data;
    			} else {
      				data['new_value'] = params.value;
      				return data;
    			}
  		},
  		title: 'Select schematic file',
  		url: 'cgi/save_datasetinfo_changes.cgi',
  		display: function(value, response) { return false;},
  		success: function(data, newValue, oldValue) {
    			$('#schematic_' + dataset_id).empty();
    			if (data['dataset']) {

      				// For Chrome (was inserting '\fakepath\' before actual filename)
      				$('#schematic_' + dataset_id).editable('destroy');
      				if ( newValue.indexOf('\\fakepath\\') >= 0) {
      					     newValue = newValue.split("\\fakepath\\")[1]; //remove Chrome
      				}
      				$(this).attr('src', 'uploads/files/' + newValue);
    			}

    			// user doesn't own dataset or unexpected field was edited
    			if (data['error']) {
      				$(this).editable('destroy');
      				$(this).attr('src', oldValue);
      				$('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
      					'<button type="button" class="close close-alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
      					'<p class="alert-message"><strong>Oops! </strong> ' + data["error"] + '</p></div>').show();
    			}
  		},
  		error: function(data, newValue, oldValue) {
    			//alert("ooops");
    			$(this).text(oldValue);
    			if (data['error']) {
      				console.log(data['error']);
      				$(this).editable('destroy');
      				$(this).attr('src', oldValue);
      				$('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
      					'<button type="button" class="close close-alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
      					'<p class="alert-message"><strong>Oops! </strong>Sorry, an error has occurred. Please try again.</p></div>').show();
    			}
  		}
  	});
});

//Sort datasets by date uploaded
var sort_date_click_count = null;
$(document).on('click', 'button#btn_sort_date', function() {
    sort_date_click_count +=1;
    // console.log(sort_date_click_count, 'click!');

    // odd number sorts newest --> oldest
    if (sort_date_click_count % 2) {

        //Nice, Quick Source: http://stackoverflow.com/a/7211762/2900840
        $("div.dataset").sort(function(a,b){
            return new Date($(a).attr("data-date-added")) > new Date($(b).attr("data-date-added"));
          }).each(function(){
              $("div#dataset_list_c").prepend(this);
          });

        //change button glyphicon
        $('button#btn_sort_date').replaceWith('<button type="button" class="btn btn-default btn-sm" ' +
              'id="btn_sort_date" data-toggle="tooltip" data-placement="bottom" ' +
              'title="Sort datasets by upload date">' +
              '<span class="glyphicon glyphicon-sort-by-attributes-alt" aria-hidden="true"></span> Date</button>');

    // even number sorts oldest --> newest
    } else {
        //Nice, Quick Source: http://stackoverflow.com/a/7211762/2900840
        $("div.dataset").sort(function(a,b){
            return new Date($(a).attr("data-date-added")) < new Date($(b).attr("data-date-added"));
          }).each(function(){
              $("div#dataset_list_c").prepend(this);
          });

        //change button glyphicon
        $('button#btn_sort_date').replaceWith('<button type="button" class="btn btn-default btn-sm" ' +
              'id="btn_sort_date" data-toggle="tooltip" data-placement="bottom" ' +
              'title="Sort datasets by upload date">' +
              '<span class="glyphicon glyphicon-sort-by-attributes" aria-hidden="true"></span> Date</button>');
    }
});

//Sort datasets by title
var sort_title_click_count = null;
$(document).on('click', 'button#btn_sort_title', function() {
    sort_title_click_count +=1;
    // console.log(sort_title_click_count, 'click!');

    // odd number sorts A --> Z
    if (sort_title_click_count % 2) {

        // complex string sorting: http://www.w3schools.com/js/js_array_sort.asp
        $('div.dataset').sort(function(a,b) {
            var a_title = a.getAttribute('data-title').toLowerCase();
            var b_title = b.getAttribute('data-title').toLowerCase();
            if (a_title < b_title) {return -1;}
            if (a_title > b_title) {return 1;}
            return 0;
        }).appendTo( $('div#dataset_list_c') );

        //change button to Z->A glyphicon
        $('button#btn_sort_title').replaceWith('<button type="button" class="btn btn-default btn-sm" ' +
              'id="btn_sort_title" data-toggle="tooltip" data-placement="bottom" ' +
              'title="Sort datasets by title">' +
              '<span class="glyphicon glyphicon-sort-by-alphabet-alt" aria-hidden="true"></span> Title</button>');

        $('.dataset');

    // even number sorts Z --> A
    } else {

        // complex string sorting: http://www.w3schools.com/js/js_array_sort.asp
        $('div.dataset').sort(function(a,b) {
            var a_title = a.getAttribute('data-title').toLowerCase();
            var b_title = b.getAttribute('data-title').toLowerCase();
            if (a_title < b_title) {return 1;}
            if (a_title > b_title) {return -1;}
            return 0;
        }).appendTo( $('div#dataset_list_c') );

        //change button to A->Z glyphicon
        $('button#btn_sort_title').replaceWith('<button type="button" class="btn btn-default btn-sm" ' +
              'id="btn_sort_title" data-toggle="tooltip" data-placement="bottom" ' +
              'title="Sort datasets by upload title">' +
              '<span class="glyphicon glyphicon-sort-by-alphabet" aria-hidden="true"></span> Title</button>');
    }
});

// Filter datasets
$(document).on('change', '.filterby', function() {

    if ( $(".filterby:checked").length == 0 ) {
        // show all if no filters are checked
        $("div.dataset").show();
    } else {
        // one or filters checked
        $("div.dataset").hide();

        // collect the filters that are checked
        var filter_organisms = [];
        var filter_dtypes = [];
        $(".filterby:checked").each(function(i, item) {
            if ( $(item).data('category') == 'dtype' ) {
                filter_dtypes.push(item.value);

                // Combine 'image-static' dtype with 'linegraph-standard'
                // TODO (only til chick BP & Utricle are reloaded)
                if (item.value == "linegraph-standard") {
                    filter_dtypes.push("image-static");
                }

                // Combine 'image-static-standard' dtype with 'violin-standard'
                if (item.value == "violin-standard") {
                    filter_dtypes.push("image-static-standard");
                }
            }
            if ( $(item).data('category') == 'organism' ) {
                filter_organisms.push(item.value);
            }
        });

        var dtypes_cnt = filter_dtypes.length;
        var organisms_cnt = filter_organisms.length;
        $("div.dataset").each(function(i, el) {
            var div = $(el);

            switch(true) {
                case (dtypes_cnt > 0 && organisms_cnt == 0):
                    if ( $.inArray(div.data('dtype'), filter_dtypes) > -1 ) {
                        div.show();
                        break;
                    }
                case (dtypes_cnt == 0 && organisms_cnt >= 0):
                    if ( $.inArray(div.data('organism'), filter_organisms) > -1 ) {
                        div.show();
                        break;
                    }
                case (dtypes_cnt > 0 && organisms_cnt > 0):
                    if ( $.inArray(div.data('dtype'), filter_dtypes) > -1 && $.inArray(div.data('organism'), filter_organisms) > -1) {
                        div.show();
                        break;
                    }
                default:
                    div.hide();
            }
        }); // end each div.dataset
    }

    // Add count of each filter
    $(".filterby").each(function(i, item) {
        var count_id = item.value.toLowerCase() + '-count';

        if ( $(item).data('category') == 'dtype' ) {

            if (item.value == 'linegraph-standard') {
                // Combine 'image-static' dtype with 'linegraph-standard'
                // TODO (only til chick BP & Utricle are reloaded)
                $('#' + count_id).text( "(" + $('div.dataset:visible[data-dtype="' + item.value + '"], div.dataset:visible[data-dtype="image-static"]').length.toString() + ")" );
            } else if (item.value == 'violin-standard') {
                // Combine 'image-static-standard' dtype with 'violin-standard'
                $('#' + count_id).text( "(" + $('div.dataset:visible[data-dtype="' + item.value + '"], div.dataset:visible[data-dtype="image-static-standard"]').length.toString() + ")" );
            } else {
              $('#' + count_id).text( "(" + $('div.dataset:visible[data-dtype="' + item.value + '"]').length.toString() + ")" );
            }
        }
        if ( $(item).data('category') == 'organism' ) {
            $('#' + count_id).text( "(" + $('div.dataset:visible[data-organism="' + item.value + '"]').length.toString() + ")" );
        }

        if ( $("#" + count_id).text().indexOf('(0)') !== -1 ) {
            $("#" + count_id).parent().addClass('no_results');
        } else {
            $("#" + count_id).parent().removeClass('no_results');
        }
    });

    //switch back to list view (if needed)
    if ( $('#dataset_arrangement_c:visible') ) {
        $('#btn_list_view').trigger('click');
    }
});

// Generate a list of existing tags from database
// Autocompletes in TAG input as twitter-like tokens
function get_tag_list(dataset_id) {
    console.log('Getting tags...');
  	$.get('./cgi/get_tag_list.cgi', function(data) {
    		if (data['success'] == 1) {
    		    // put all tags into a list
            var tokenized_tags = [];
            for (i=0; i < data['tags'].length; i++) {
                tokenized_tags.push(data['tags'][i]['label']);
            }

            //initialize tokenfield for tags
      			$('#tag_input_' + dataset_id).tokenfield({
      			    autocomplete: {
      			        source: tokenized_tags,
                    ShowAutocompleteOnFocus: true
      			    },
                delimiter: [',', ' ']
			       });
             console.log("Tags received.");
    		} else {
    		    console.log("Handle a failed report from the CGI");
    		}
  	})
  	.fail(function() {
  		  console.log("Something broke. Unable to generate tag list");
  	});
} //end get_tag_list()

// Set schematics max-height to match height of metadata
function set_schematics_maxheight() {
    $('div.schematic-c').each(function() {
        var h = $(this).siblings('div.dataset-info-c').height().toString() + 'px';
        $(this).children().css("max-height", h);
    });
}
