var sort_note_date_click_count = null;
var sort_note_title_click_count = null;
var popout_dataset_title = null;
var popout_window = false;
var parent_window = null;
var notewindow_id = null;

var note_location = 'right'; //default. can also be 'left' or 'popout'
var slide_direction = ''; //'right' for right side panel. 'left' for left side panel.
if ( $('#note_side_panel_c').hasClass('dock-right') ) {
    slide_direction = 'right';
} else {
    slide_direction = 'left';
}

$(document).ready(function() {
    // Initially load side panel for notes
    $('#note_side_panel_c').load('./include/sidepanel_notes.html');

    //Load if full-screen notes URL
    if (document.URL.indexOf("?notewindow_id=") >= 0 ) {
        slide_direction = 'left';
        popout_window = window;
        popout_window.window.opener.popout_window = window;

        $('#note_fullwindow_c').load('./include/fullwindow_notes.html');
        notewindow_id = getUrlParameter('notewindow_id');

        //Get notes for the popout
        get_dataset_notes(notewindow_id, 'popout');

        // Format sidepanel for popout window
        // Apply sidepanel css
        $('#note_side_panel_c').removeClass('dock-right dock-left').addClass('dock-popout');
        // Change sidepanel options - popin only
        $('#sidepanel_options_docked').hide();
        $('#sidepanel_options_popped').show();

        // Hide non-side panel stuff
        $('#note_panel_closer').hide();
        $('#navigation_bar, #page_name, #controls, #main_panel, #site_intro_c, #sidebar, footer').hide();

        // show fullwindow text editing area
        $('#note_fullwindow_c').show();

        // show everything at once. Gives a cleaner show.
        $('#body_c').show();
    } else {
        // show everything at once. Gives a cleaner show.
        $('#body_c').show();

        //save for later. Helps with interwindow communication
        parent_window = window;
    }
});

function get_dataset_notes(dataset_id, scope) {
    //scope = 'sidepanel' || 'popout'
    $.ajax({
        url: './cgi/get_dataset_notes.cgi',
        type: 'POST',
        data: { 'session_id': session_id, 'dataset_id': dataset_id, 'scope': scope },
        dataType: 'json',
        success: function(data, textStatus, jqXHR) {
            if (data['success'] == 1) {

                if (data['notes'].length > 0) {
                    for (var i=0,len=data['notes'].length; i<len; i++) {

                        //make date_last_changed nicer looking
                        data['notes'][i]['date_formatted'] = formatDate(data['notes'][i]['date_last_changed']);

                        // Make a preview of the note
                        var ldesc_preview = "";
                        if (data['notes'][i]['ldesc'] == null) {
                            ldesc_preview = "...";
                        } else if (data['notes'][i]['ldesc'].length > 50 ) {
                            ldesc_preview = data['notes'][i]['ldesc'].substring(0, 49) + "...";
                        } else {
                            ldesc_preview = data['notes'][i]['ldesc'] + "...";
                        }
                        data['notes'][i]['ldesc_preview'] = ldesc_preview;

                        // Insert line-breaks into dataset descriptions (for public notes that user does not own)
                        if ( data['notes'][i]['is_owner'] == 0 ) {
                            var text = data['notes'][i]['ldesc'].replace(/(\r\n|\n\r|\r|\n)/g, '<br />');
                            data['notes'][i]['ldesc'] = text;
                        }

                    }

                    // Populate sidepanel heading with dataset's title
                    if ( scope == 'popout' ) {
                        // Add the dataset title to popout header
                        $('#popout_dataset_title').text( data['notes'][0]['dataset_title'] );
                    } else {
                        $('#dataset_title').text( data['notes'][0]['dataset_title'] );
                    }

                    // Populate the note side panel
                    var notepanelbodyViewTmpl = $.templates('#tmpl_notetable_body');
                    var notepanelbodyViewHtml = notepanelbodyViewTmpl.render(data['notes']);
                    $('#note_table_body').html(notepanelbodyViewHtml);

                    // Populate the note side panel
                    var editnoteViewTmpl = $.templates('#tmpl_edit_note');
                    var editnoteViewHtml = editnoteViewTmpl.render(data['notes']);
                    $('#edit_notes_c').html(editnoteViewHtml);

                } else {
                    var no_notes_msg =
                        '<br />' +
                        '<h3 style="text-align:center;">No notes were found.</h3>' +
                        '<br /><br /><br /><br /><br /><br />' +
                        '<h2 style="text-align:center;">Create a new note below.</h2>'
                    $('#note_table_body').html(no_notes_msg);

                    // Populate popout heading with dataset's title
                    if (scope == 'popout') {
                        $('#popout_dataset_title').text( data['dataset_title'][0]['dataset_title'] );
                    } else {
                        $('#dataset_title').text( data['dataset_title'][0]['dataset_title'] );
                    }
                }

                // Finally, display the note_side_panel
                if (scope == 'popout') {
                    $('#note_dataset_title_c').hide();
                    $('#sidepanel_options_docked').hide();

                    $('#sidepanel_options_popped').show();
                    $('#note_side_panel_c').show();
                } else {
                    $('#note_side_panel_c').show('slide', {direction: slide_direction}, 400);
                }

            } else {
                console.log(data['error']);

                $('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
                  '<button type="button" class="close close-alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
                  '<p class="alert-message">' + data['error'] + '</p>' +
                  '</div>').show();
            }
        },
        error: function(jqXHR, textStatus, errorThrown) {
            display_error_bar(jqXHR.status + ' ' + errorThrown.name);
        }
    }); //end ajax
} //end get_dataset_notes()

function formatDate(date) {
    // Formats a date to "Today", "Yesterday", or "DD/MM/YYYY" and returns it
    var date_str = '';

    var today = new Date();
    var date_last_changed = new Date(date);

    var time_today = (today.getHours() * 3600) + (today.getMinutes() * 60) + today.getSeconds(); //number of seconds we are into the day
    var sum_time_yesterday = time_today + 86400; // sum of a day (86400 seconds) and the part of day that has passed
    var time_delta = Math.floor( (today - date_last_changed) / 1000 ); // length of time between note and now (in seconds)

    // When was note posted? Format it as "Today", "Yesterday", or "DD/MM/YYYY"
    if ( time_delta < time_today ) {
        date_str = "Today";
    } else if ( time_delta > time_today && time_delta < sum_time_yesterday ) {
        date_str = "Yesterday";
    } else if ( time_delta > sum_time_yesterday ) {
        var day = date_last_changed.getDate();
        var month = date_last_changed.getMonth() + 1;
        var year = date_last_changed.getFullYear();

        date_str = day + '/' + month + '/' + year;
    } else {
        //Format any missed date (This should not happen though)
        var day = date_last_changed.getDate();
        var month = date_last_changed.getMonth() + 1;
        var year = date_last_changed.getFullYear();

        date_str = day + '/' + month + '/' + year;
    }
    return date_str;

} // end formatDate()

// Launches the side panel
$(document).on('click', 'li.note_launcher', function(e) {
    var dataset_id = $(this).attr('data-dataset-id');

    // Prepopulate the form dataset_id
    $('input#dataset_id').val(dataset_id);

    // Set popout window's url
    // remove '#' from URL if present
    if (window.location.href.indexOf('#') ) {
        var popout_url = window.location.href.split('#')[0] + '?notewindow_id=' + dataset_id;
    } else {
        var popout_url = window.location + '?notewindow_id=' + dataset_id;
    }
    $('#movenote2popout').attr('data-popout-url', popout_url);

    //// Load notes. ////
    // If popout has been closed or hasn't been opened yet
    // I'm not sure why this works better than if(popout_window)
    if ( popout_window.closed  || !popout_window) {

        //Ensures parent popout_window becomes 'null', not sure why
        doUnloadPopout();
        // Get the notes for that dataset
        get_dataset_notes(dataset_id, 'sidepanel');

        // Hide edit_note and add_note panels (if open)
        if ( $('#edit_notes_c:visible') ) {
            $('#edit_notes_c').hide('slide', {direction: slide_direction}, 400);
            $('.edit_form').hide();
        }
        if ( $('.unable2edit:visible') ) {
            $('.unable2edit').hide('slide', {direction: slide_direction}, 400);
        }
        if ( $('div#add_note_c:visible') ) {
            $('div#add_note_c').hide('slide', {direction: slide_direction}, 400);
            $('input#add_note_title, textarea').val('');
        }
    } else {
        // Prevent the URL change from initializing reset_parent_window()
        dontUnloadPopout();

        // Refresh Popout with the dataset's notes
        popout_window.location.href = popout_url;
        popout_window.focus();
    }
});

// Toggle Docking locations (left, right, popout)
$(document).on('click', '#movenote2left', function(e) {
    $('#note_side_panel_c').removeClass('dock-right dock-popout').addClass('dock-left');
    popout_window = false;
    $('a#main_logo').css({'height': '50px'});

    slide_direction = 'left';
});
$(document).on('click', '#movenote2right', function(e) {
    $('#note_side_panel_c').removeClass('dock-left dock-popout').addClass('dock-right');
    popout_window = false;
    $('a#main_logo').css({'height': '75px'});

    slide_direction = 'right';
});
$(document).on('click', '#movenote2popout', function(e) {
    // Open note popout window and create variable to later use as a way of refreshing the popout
    // Assign reset_parent_window() to run if browser's 'close' is clicked
    //Source: http://stackoverflow.com/a/15425536/2900840
    popout_window = window.open( $('#movenote2popout').attr('data-popout-url'), '', 'fullscreen=yes, scrollbars=yes' );
    popout_window.onbeforeunload = unloadPage;

    // Hide note_panel and change class
    $('#note_side_panel_c').hide('slide', {direction: slide_direction}, 400);
});
// Only unload page if closed. Source: http://stackoverflow.com/a/18310601/2900840
var resetParentOnClose = true;
function unloadPage() {
    if (resetParentOnClose) {
        reset_parent_window('close');
    }
    if (window.closed) {
        reset_parent_window('close');
    }
}
function dontUnloadPopout() {
    //Don't reset parent_window
    resetParentOnClose = false;
}
function doUnloadPopout() {
    //Do reset parent_window
    resetParentOnClose = true;
}

$(document).on('click', '#movenote2popin', function(e) {
    // former code moved here...
    reset_parent_window('open');
});

function reset_parent_window(sidepanel_status) {
    // Resets the parent window and sidepanel to onload/default
    // sidepanel_status = 'open' || 'close'. Have sidepanel open or closed after reset.

    // Get parent window
    var p_opener = popout_window.opener;
    // Reset parent's popout
    p_opener.popout_window = false;

    //reset sidepanel position and css
    var parent_sidepanel = p_opener.document.getElementById('note_side_panel_c');
    $( parent_sidepanel ).removeClass('dock-left dock-popout')
      .addClass('dock-right')
      .show();
    p_opener.slide_direction = 'right';

    //reset gEAR logo (if user popped out from left-side)
    var parent_logo = p_opener.document.getElementById('main_logo');
    $( parent_logo ).css('height', '75px');

    // Open parent's sidepanel?
    if ( sidepanel_status == 'open' ) {
        p_opener.get_dataset_notes(notewindow_id, 'sidepanel');
    } else {
        var p_sidepanel = p_opener.document.getElementById('note_side_panel_c');
        $(p_sidepanel).hide('slide', {direction: slide_direction}, 400);

        var p_sidepanel_table = p_opener.document.getElementById('note_table_body');
        $(p_sidepanel_table).empty();
    }

    if (popout_window) {
        popout_window.close();
    }
}

// Hide and reset the note panel
$(document).on('click', '#note_panel_closer', function(e) {
    $('#note_side_panel_c').hide('slide', {direction: slide_direction}, 400);
    $('#note_table_body').empty();
});

$(document).on('click', '.note_preview', function(e) {
    var note_id = $(this).attr('data-note-id');

    // Hide notes that might already be showing
    $('#edit_notes_c').hide();
    $('.able2edit, .unable2edit').hide();

    if ( popout_window ) {
        // Move and adjust edit_notes_c for full window
        $('#note_fullwin_body_c').append( $('#edit_notes_c') );
        $('#edit_notes_c').css('padding-right', '30px');
    }

    // Display contents first, then slide them in view
    $('#view_note_c_' + note_id).show();
    $('#edit_notes_c').show('slide', {direction: slide_direction}, 400);
    $('#no_selected_note_c').hide();
});

$(document).on('click', '.close-public-note', function(e) {
    $('#no_selected_note_c').show();
    $('#edit_notes_c').hide('slide', {direction: slide_direction}, 400);
    $('.unable2edit').hide();
});

$(document).on('click', '.btn_cancel_edit_note', function(e) {
    var note_id = $(this).attr('data-note-id');
    $('#no_selected_note_c').show();

    // Slide contents out, then hide them
    $('#edit_notes_c').hide('slide', {direction: slide_direction}, 400);
    $('#view_note_c_' + note_id).hide();
});

$(document).on('click', 'button#add_note', function(e){
    if ( popout_window ) {
      // Move and adjust add_note_c for full window
      $('#note_fullwin_body_c').append( $('#add_note_c') );
      $('#add_note_c').css('padding-right', '30px');
    }

    // Use timestamp as a backup for title
    var d = new Date().toString();
    d = d.replace(/\:\d+\s.+/g, ''); //remove seconds and everything after
    var placeholder = 'Note from ' + d;
    $('#add_note_title').attr( 'placeholder', placeholder );

    // Slide the contents into view
    $('div#add_note_c').show('slide', {direction: slide_direction}, 400);
});

//Save a new note or the changes made to existing note
$(document).on('click', 'button#save_new_note, button.btn_save_edit_note', function(e){
    // What change is user performing on the note?
    var el = e.target;
    var scope = '';
    if( $(el).hasClass('btn_save_edit_note') ) {
        scope = 'edit';
    } else {
        scope = 'new';
    }

    //Where is the change occurring? 'sidepanel' vs 'popout'
    var note_location = '';
    if ( popout_window ) {
        note_location = 'popout';
    } else {
        note_location = 'sidepanel';
    }

    // Collect note information
    var note_id = '';
    var text = '';
    var ldesc = '';
    var access_level = '';

    // Get dataset_id
    var dataset_id = '';
    if ( popout_window ) {
        dataset_id = notewindow_id;
    } else {
        dataset_id = $('input#dataset_id').val();
    }
    //If saving changes to existing note...
    if( scope == 'edit' ) {
        note_id = $(el).attr('data-note-id');
        text = $('input#edit_note_title_' + note_id).val();
        if ( text.length > 0 ) {
            title = text;
        } else {
            title =  $('input#edit_note_title_' + note_id).attr('placeholder');
        }
        ldesc = $('#edit_note_text_' + note_id).prop('value');
        access_level = $('#edit_note_access_level_' + note_id).val();

    // Creating a new note...
    } else {
        text = $('input#add_note_title').val();
        if ( text.length > 0 ) {
            title = text;
        } else {
            title =  $('input#add_note_title').attr('placeholder');
        }
        ldesc = $('#add_note_text').prop('value');
        access_level = $('#add_note_access_level').val();
    }

    $.ajax({
        url: './cgi/apply_note_changes.cgi',
        type: 'POST',
        data: {
            'session_id': session_id, 'dataset_id': dataset_id, 'title': title,
            'ldesc': ldesc, 'access_level': access_level,
            'scope': scope, 'note_id': note_id
        },
        dataType: 'json',
        success: function(data, textStatus, jqXHR) {
            if (data['success'] == 1) {

                // Get the notes for that dataset
                get_dataset_notes(dataset_id, note_location);

                if ( scope == 'edit') {
                    // Hide edit panel for that note
                    $('div#edit_notes_c').hide('slide', {direction: slide_direction}, 400);
                } else {
                    // Hide add note panel and clear it
                    $('div#add_note_c').hide('slide', {direction: slide_direction}, 400);
                    $('input#add_note_title, textarea').val('');
                }

            } else {
                console.log(data['error']);
                display_error_bar(data['error']);
            }
        },
        error: function(jqXHR, textStatus, errorThrown) {
            display_error_bar(jqXHR.status + ' ' + errorThrown.name);
        }
    }); //end ajax
});

// Delete user-owned note
$(document).on('click', '.btn_remove_note', function(e) {
    var note_id = $(e.target).attr('data-note-id');
    var scope = 'remove';

    $.ajax({
        url: './cgi/apply_note_changes.cgi',
        type: 'POST',
        data: {'session_id': session_id, 'note_id': note_id, 'scope': scope},
        dataType: 'json',
        success: function(data, textStatus, jqXHR) {
            if (data['success'] == 1) {
                // Hide edit_note_c. Then hide the deleted note's row
                $('div#edit_notes_c').hide('slide', {direction: slide_direction}, 400);
                $('#view_note_c_' + note_id).hide();
                $('#row_' + note_id).hide('fade',{}, 400);
            } else {
                display_error_bar(data['error']);
            }
        },
        error: function(jqXHR, textStatus, errorThrown) {
            display_error_bar(jqXHR.status + ' ' + errorThrown.name);
        }
    }); //end .ajax
});

$(document).on('click', 'button#cancel_new_note', function(e){
    // Hide the New Note form
    $('div#add_note_c').hide('slide', {direction: slide_direction}, 400);

    //clear the form
    $('input#add_note_title, textarea').val('');
});

// Toggle Hide/Show Public notes
$(document).on('click', '#btn_hide_public', function(e) {
    //Hide public notes
    $('.note_preview[data-access="Public"]').hide();

    //change button icon
    $('#btn_hide_public').replaceWith('<button type="button" class="btn btn-xs" ' +
          'id="btn_show_public" data-toggle="tooltip" data-placement="bottom" ' +
          'title="Show public notes">' +
          '<span class="fa fa-eye" aria-hidden="true"></span> Public</button>');
});
$(document).on('click', '#btn_show_public', function(e) {
    //Show public notes
    $('.note_preview[data-access="Public"]').show();

    //change button icon
    $('#btn_show_public').replaceWith('<button type="button" class="btn btn-xs" ' +
          'id="btn_hide_public" data-toggle="tooltip" data-placement="bottom" ' +
          'title="Hide public notes">' +
          '<span class="fa fa-eye-slash" aria-hidden="true"></span> Public</button>');
});

//Sort datasets by date uploaded
$(document).on('click', 'button#btn_sort_note_date', function() {
    sort_note_date_click_count +=1;

    // odd number sorts newest --> oldest
    if (sort_note_date_click_count % 2) {

        //Nice, Quick Source: http://stackoverflow.com/a/7211762/2900840
        $(".note_preview").sort(function(a,b){
            return new Date($(a).attr("data-date-changed")) > new Date($(b).attr("data-date-changed"));
          }).each(function(){
              $("#note_table_body").prepend(this);
          });

        //change button icon
        $('button#btn_sort_note_date').replaceWith('<button type="button" class="btn btn-xs" ' +
              'id="btn_sort_note_date" data-toggle="tooltip" data-placement="bottom" ' +
              'title="Sort notes by date">' +
              '<span class="fa fa-sort-amount-desc" aria-hidden="true"></span> Date</button>');

    // even number sorts oldest --> newest
    } else {
        //Nice, Quick Source: http://stackoverflow.com/a/7211762/2900840
        $(".note_preview").sort(function(a,b){
            return new Date($(a).attr("data-date-changed")) < new Date($(b).attr("data-date-changed"));
          }).each(function(){
              $("#note_table_body").prepend(this);
          });

        //change button icon
        $('button#btn_sort_note_date').replaceWith('<button type="button" class="btn btn-xs" ' +
              'id="btn_sort_note_date" data-toggle="tooltip" data-placement="bottom" ' +
              'title="Sort notes by date">' +
              '<span class="fa fa-sort-amount-asc" aria-hidden="true"></span> Date</button>');
    }
});

//Sort notes by title
$(document).on('click', 'button#btn_sort_note_title', function() {
    sort_note_title_click_count +=1;

    // odd number sorts A --> Z
    if (sort_note_title_click_count % 2) {

        // complex string sorting: http://www.w3schools.com/js/js_array_sort.asp
        $('.note_preview').sort(function(a,b) {
            var a_title = a.getAttribute('data-title').toLowerCase();
            var b_title = b.getAttribute('data-title').toLowerCase();
            if (a_title < b_title) {return -1;}
            if (a_title > b_title) {return 1;}
            return 0;
        }).appendTo( $('#note_table_body') );

        //change button to Z->A icon
        $('button#btn_sort_note_title').replaceWith('<button type="button" class="btn btn-xs pull-right" ' +
              'id="btn_sort_note_title" data-toggle="tooltip" data-placement="bottom" ' +
              'title="Sort notes by title">' +
              '<span class="fa fa-sort-alpha-desc" aria-hidden="true"></span> Title</button>');

        $('.dataset');

    // even number sorts Z --> A
    } else {

        // complex string sorting: http://www.w3schools.com/js/js_array_sort.asp
        $('.note_preview').sort(function(a,b) {
            var a_title = a.getAttribute('data-title').toLowerCase();
            var b_title = b.getAttribute('data-title').toLowerCase();
            if (a_title < b_title) {return 1;}
            if (a_title > b_title) {return -1;}
            return 0;
        }).appendTo( $('#note_table_body') );

        //change button to A->Z icon
        $('button#btn_sort_note_title').replaceWith('<button type="button" class="btn btn-xs pull-right" ' +
              'id="btn_sort_note_title" data-toggle="tooltip" data-placement="bottom" ' +
              'title="Sort notes by title">' +
              '<span class="fa fa-sort-alpha-asc" aria-hidden="true"></span> Title</button>');
    }
});

$('.note_access').popover({
	animation: true,
	trigger: 'click',
	// title: "Delete profile",
	content:
      '<div id="public_note_warning_add" class="alert alert-warning public_note_warning" role="alert">' +
          '<b>Caution!</b> Public notes are visible to everyone. Your name will also be displayed.' +
      '</div>',
	html: true,
	placement: 'auto'
});

// #add_note_access_level
// Source show popover for specific <select><option>: http://stackoverflow.com/a/13129651/2900840
$(document).on('click', '.note_access', function(e) {
    $('.note_access').popover('destroy');

    if ( $(e.target).attr('id') == 'add_note_access_level' ) {
        if ( $(e.target).val() == "1" || $(e.target).val() == "Public" ) {
            $('#add_note_access_level').popover({
                trigger: "manual",
                placement: "auto",
              	html: true
            }).popover('show');
        } else {
            $('#add_note_access_level').popover('destroy');
        }
    } else {
        var note_id = $(e.target).attr('data-note-id');
        if ( $(e.target).val() == "1" || $(e.target).val() == "Public" ) {
            $('#edit_note_access_level_' + note_id).popover({
                trigger: "manual",
                placement: "auto",
              	html: true
            }).popover('show');
        } else {
            $('#edit_note_access_level_' + note_id).popover('destroy');
        }
    }
});
