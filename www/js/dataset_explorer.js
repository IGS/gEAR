var first_search = true;
var animation_time = 200;
var current_layout_dataset_ids = [];
var dataset_id_to_delete = null;

// Values are recent, search, or a profile share ID
var current_dataset_list_label = 'recent';

// toggles whether we are in profile management mode
var mgmt_mode = false;

window.onload=function() {
    // check if the user is already logged in
    check_for_login();

    $('#selected_layout').on('change', function() {
        update_layout_dataset_id_list();
        update_view_buttons();
    });

    $("#search_clear").click(function(){
        $("#search_terms").val('');
        submit_search();
    });

    $('#search_terms').keyup(function() {
        if ($("#search_terms").val().length > 0) {
            $("#search_clear").show();
        } else {
            $("#search_clear").hide();
        }
    });

    $('#submit_search').submit(function(event) {
        event.preventDefault(); 
        submit_search();
    });

    $('#sort_by').on('change', function() {
        submit_search();
    });

    $(document).on('click', 'button.curate', function(e){
        window.location = "./dataset_curator.html#/dataset/" + $(this).val()  + "/displays";
    });

    $(document).on('click', 'button.view_dataset', function(e){
        window.location = "./p?s=" + $(this).val();        
    });

    // Generic function to handle all collapsable menus
    // h.expandable_control is clicked and looks for plus/minus icons as siblings
    // and an .expandable_target as a direct child
    $(document).on('click', "h4.expandable_control", function() {
        var exblock = $(this).siblings(".expandable_target")[0];
        if ($(exblock).is(":visible")) {
            $(this).children(".fa-plus").show();
            $(this).children(".fa-minus").hide();
            $(exblock).hide(animation_time);

            if ($(this).siblings(".profile_control").length) {
                $(".profile_control").hide();
                $("#btn_arrangement_view").hide();
                mgmt_mode = false;
            }
        } else {
            $(this).children(".fa-plus").hide();
            $(this).children(".fa-minus").show();
            $(exblock).show(animation_time);

            if ($(this).siblings(".profile_control").length) {
                $(".profile_control").show();

                if ($('#selected_layout').find(':selected').data('is_domain') == "0") {
                    $("#btn_arrangement_view").show();
                }
                mgmt_mode = true;
            }
        }
    });

    $("#initial_instructions_bar").on('click', function() {
        if ($("#initial_instructions_body").is(":visible")) {
            $("#initial_instructions_body").hide(animation_time);
        } else {
            $("#initial_instructions_body").show(animation_time);
        }
    });

    $("#initial_instructions_closer i").on('click', function() {
        $("#initial_instructions_c").hide();
    });

    // Generic function to handle the facet selector choices
    //  For any ul.controls_filter_options the list elements can have a class="selected"
    //  The groups of <li> also have one/top li with class="all_selector" which
    //  toggles the rest of them off since no filter is applied.
    $(document).on('click', "ul.controls_filter_options li", function() {
        // if the one clicked is the all_selector then highlight it and unclick the rest
        if ($(this).hasClass('all_selector')) {
            if (! $(this).hasClass('selected')) {
                $(this).addClass('selected');
            }

            $(this).siblings().removeClass('selected');
        } else {
            if (! $(this).hasClass('selected')) {
                // If turning on, make sure all_selector is off                
                $(this).parent().find("li.all_selector").removeClass('selected');

                // If this selection group has the 'only_one' option deselect the rest
                if ($(this).parent().hasClass('only_one')) {
                    $(this).siblings().removeClass('selected');
                }
                
                $(this).addClass('selected');
            } else {
                // If turning off, make sure at least one other option is selected, else set
                //  set all_selector on
                $(this).removeClass('selected');

                if ($(this).parent().children("li.selected").length == 0) {
                    $(this).parent().find("li.all_selector").addClass('selected');
                }
            }
        }

        submit_search();
    });

     // keep the add layout button disabled unless a full new name has been entered
    $(document).on('keyup', '#new_layout_name', function() {
        if ($('#new_layout_name').val().length == 0) {
            $('#confirm_layout_add').prop('disabled', true);
        } else {
            $('#confirm_layout_add').prop('disabled', false);
        }
    });

    $(document).on('click', 'span.dataset-expander', function() {
        var dataset_id = $(this).data('dataset-id');
        var selector_base = "#result_dataset_id_" + dataset_id;

        if ($(selector_base + " .expandable-view").hasClass('expanded-view-hidden')) {
            $(selector_base + " .expandable-view").removeClass('expanded-view-hidden');
            $(selector_base + " .dataset-expander i").removeClass('fa-expand');
            $(selector_base + " .dataset-expander i").addClass('fa-compress');
        } else {
            $(selector_base + " .expandable-view").addClass('expanded-view-hidden');
            $(selector_base + " .dataset-expander i").removeClass('fa-compress');
            $(selector_base + " .dataset-expander i").addClass('fa-expand');
        }
    });
    
    $(document).on('click', 'button.add2profile', function() {
        dataset_id = $(this).attr('value');
        add_to_profile(dataset_id);
    });

    $(document).on('click', 'button.edit_dataset', function() {
        var dataset_id = $(this).data('dataset-id');
        var selector_base = "#result_dataset_id_" + dataset_id;

        // Show editable versions where there are some and hide the display versions
        $(selector_base + " .is-editable").hide();
        $(selector_base + " .editable-version").show();

        // Make sure the view is expanded
        if ($(selector_base + " .expandable-view").hasClass('expanded-view-hidden')) {
            $(selector_base + " span.dataset-expander").click();
        }
    });

    $(document).on('click', 'button.removefromprofile', function() {
        dataset_id = $(this).attr('value');
        remove_from_profile(dataset_id);
    });

    $(document).on('click', 'button.share_dataset', function() {
        share_id = $(this).attr('value');
        var current_url = window.location.href;
        var current_page = current_url.lastIndexOf("dataset_explorer.html");
        var share_url = current_url.substring(0, current_page) + 'p?s=' + share_id;
        var dataset_id = $(this).data('dataset-id');

        if (copyToClipboard(share_url)) {
            show_dataset_action_note(dataset_id, "URL copied to clipboard");
        } else {
            show_dataset_action_note(dataset_id, "Failed to copy to clipboard. URL: " + share_url);
        }
    });
};  // end window onloads

// setup alert-container to give 'hide' function to alert button
// alert doesn't exist at window load therefore it needs to inherit from the parent div
$( ".alert-container" ).on("click", "button.close-alert", function() {
	 $( ".alert-container" ).hide();
});

$('#btn_add_layout').popover({
	animation: true,
	trigger: 'click',
	title: "Add profile",
	content: "<p>Enter a name for a profile to add</p>" +
		"<div class='btn-toolbar' style='width:250px'>" +
        "<input type='text' name='new_layout_name' id='new_layout_name' class='form-control' placeholder='Profile name' maxlength='30'>" +
		"<button id='confirm_layout_add' class='btn btn-default btn-primary confirm_layout_add' data-dismiss='popover' disabled>Add</button>" +
		"<button id='cancel_layout_add' class='btn btn-default cancel_add' value='cancel_add'>Cancel</button>" +
		"</div>",
		html: true,
		placement: 'auto',
    container: 'body'
}).removeClass('actions');

$(document).on('click', '#cancel_layout_add', function() {
    $('#btn_add_layout').popover('hide');
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
                $('#selected_layout').val(data['layout_id']);
                show_layout_action_note("Profile created");
                update_add_remove_buttons();
                update_view_buttons();
            }
        },
        error: function (jqXHR, textStatus, errorThrown) {
            console.log('textStatus= ', textStatus);
            console.log('errorThrown= ', errorThrown);
            display_error_bar(jqXHR.status + ' ' + errorThrown.name);
        }
    });
});

// Popover for these is created within process_search_results()

$(document).on('click', '#cancel_dataset_delete', function() {
    dataset_id_to_delete = null;
    $('.delete_dataset').popover('hide');
});

$(document).on('click', '.confirm_dataset_delete', function() {
    $('.delete_dataset').popover('hide');

    session_id = Cookies.get('gear_session_id');

    $.ajax({
        url : './cgi/remove_dataset.cgi',
        type: "POST",
        data : { 'session_id': session_id, 'dataset_id': dataset_id_to_delete },
        dataType:"json",
        success: function(data, textStatus, jqXHR) {
            if (data['success'] == 1) {
                $("#result_dataset_id_" + dataset_id_to_delete).fadeOut("slow", function() {
                    $("#result_count").html( $("#result_count").html() - 1  );
                    $(this).remove();
                });
                dataset_id_to_delete = null;
            } else {
                display_error_bar(data['error']);
            }
        },
        error: function (jqXHR, textStatus, errorThrown) {
            console.log('textStatus= ', textStatus);
            console.log('errorThrown= ', errorThrown);
            display_error_bar(jqXHR.status + ' ' + errorThrown.name);
        }
    }); //end ajax for .confirm_delete
});

$(document).on('click', '.edit_dataset_cancel', function() {
    var dataset_id = $(this).data('dataset-id');
    var selector_base = "#result_dataset_id_" + dataset_id;

    // Show editable versions where there are some and hide the display versions
    $(selector_base + " .editable-version").hide();
    $(selector_base + " .is-editable").show();

    // Reset any unsaved/edited values
    var visibility = $(selector_base + "_visibility").data("original-val");
    $(selector_base + "_visibility").val(visibility);

    var title = $(selector_base + "_editable_title").data("original-val");
    $(selector_base + "_editable_title").val(title);

    var pubmed_id = $(selector_base + "_editable_pubmed_id").data("original-val");
    $(selector_base + "_editable_pubmed_id").val(pubmed_id);

    var geo_id = $(selector_base + "_editable_geo_id").data("original-val");
    $(selector_base + "_editable_geo_id").val(geo_id);

    var ldesc = $(selector_base + "_editable_ldesc").data("original-val");
    $(selector_base + "_editable_ldesc").val(ldesc);
});

$(document).on('click', '.edit_dataset_save', function() {
    session_id = Cookies.get('gear_session_id');
    var dataset_id = $(this).data('dataset-id');
    var selector_base = "#result_dataset_id_" + dataset_id;
    var new_visibility = $(selector_base + "_visibility").val();
    var new_title = $(selector_base + "_editable_title").val();
    var new_pubmed_id = $(selector_base + "_editable_pubmed_id").val();
    var new_geo_id = $(selector_base + "_editable_geo_id").val();
    var new_ldesc = $(selector_base + "_editable_ldesc").val();

    $.ajax({
        url : './cgi/save_datasetinfo_changes.cgi',
        type: "POST",
        data : { 'session_id': session_id,
                 'dataset_id': dataset_id,
                 'visibility': new_visibility,
                 'title': new_title,
                 'pubmed_id': new_pubmed_id,
                 'geo_id': new_geo_id,
                 'ldesc': new_ldesc
               },
        dataType:"json",
        success: function(data, textStatus, jqXHR) {
            if ( data['success'] == 1 ) {
                // Update the UI for the new values
                $(selector_base + "_visibility").data("original-val", new_visibility);

                var visibility_html = ''
                if (new_visibility == 'public') {
                    visibility_html = '<h3><span class="badge badge-light">Public dataset</span></h3>'
                } else {
                    visibility_html = '<h3><span class="badge badge-danger">Private dataset</span></h3>';
                }
                $(selector_base + "_display_visibility").html(visibility_html);
                
                $(selector_base + "_editable_title").data("original-val", new_title);
                $(selector_base + "_display_title").html(new_title);
                
                $(selector_base + "_editable_pubmed_id").data("original-val", new_pubmed_id);
                var pubmed_html = "<span class='att_label'>Pubmed</span>";

                if (new_pubmed_id) {
                    pubmed_html += "<a href='https://pubmed.ncbi.nlm.nih.gov/" + new_pubmed_id + "' target='_blank'> " + new_pubmed_id;
                    pubmed_html += " <i class='fa fa-external-link'></i></a>";
                } else {
                    pubmed_html = "Not given";
                }

                $(selector_base + "_display_pubmed_id").html(pubmed_html);

                $(selector_base + "_editable_geo_id").data("original-val", new_geo_id);
                var geo_html = "<span class='att_label'>GEO ID</span>";

                if (new_geo_id) {
                    geo_html += "<a href='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=" + new_geo_id + "' target='_blank'> " + new_geo_id;
                    geo_html += " <i class='fa fa-external-link'></i></a>";
                } else {
                    geo_html = "Not given";
                }

                $(selector_base + "_display_geo_id").html(geo_html);

                $(selector_base + "_editable_ldesc").data("original-val", new_ldesc);
                $(selector_base + "_display_ldesc").html(new_ldesc);
                
                // Put interface back to view mode.
                $(selector_base + " .editable-version").hide();
                $(selector_base + " .is-editable").show();
            }
        },
        error: function (jqXHR, textStatus, errorThrown) {
            console.log('textStatus= ', textStatus);
            console.log('errorThrown= ', errorThrown);
            display_error_bar(jqXHR.status + ' ' + errorThrown.name);
        }
    }); //end ajax
});   

$('#btn_delete_layout').popover({
		animation: true,
		trigger: 'click',
		title: "Delete profile",
		content: "<p>Are you sure you want to delete this profile?</p>" +
		"<div class='btn-toolbar' style='width:250px'>" +
		"<button class='btn btn-default btn-danger confirm_layout_delete' data-dismiss='popover'>Delete</button>" +
		"<button id='cancel_layout_delete' class='btn btn-default cancel_delete' value='cancel_delete'>Cancel</button>" +
		"</div>",
		html: true,
		placement: 'auto',
    container: 'body'
}).removeClass('actions');

$(document).on('click', '#cancel_layout_delete', function() {
    $('#btn_delete_layout').popover('hide');
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
                
                // Reset the select box to be the first option
                $('#selected_layout').val($("#selected_layout option:first").val());

                // Delete the other value
                $('option[value="' + old_value + '"]', $('#selected_layout')).remove();

                update_add_remove_buttons();
                update_view_buttons();
            }
        },
        error: function (jqXHR, textStatus, errorThrown) {
            console.log('textStatus= ', textStatus);
            console.log('errorThrown= ', errorThrown);
            display_error_bar(jqXHR.status + ' ' + errorThrown.name);
        }
    }); //end ajax for .confirm_delete
});

$(document).on('click', 'ul.layout_links li', function() {
    window.location = "./p?l=" + $(this).data('profile-share-id');  
});

$("#btn_arrangement_view").click(function(e) {
    // make sure the current list of datasets matches the currently-selected profile
    if (current_dataset_list_label != $('#selected_layout').find(':selected').data('share-id')) {
        $("#btn_view_layout_datasets").trigger('click');
    }
    
    $("#btn_arrangement_view").addClass('active');
    $("#btn_list_view_compact").removeClass('active');
    $("#btn_list_view_expanded").removeClass('active');

    $("#dataset_list_c").hide();
    $("#dataset_arrangement_c").show();
});

$("#btn_list_view_compact").click(function(e) {
    $("#btn_arrangement_view").removeClass('active');
    $("#btn_list_view_compact").addClass('active');
    $("#btn_list_view_expanded").removeClass('active');

    $("#dataset_arrangement_c").hide();
    $("#dataset_list_c").show();
    
    // find all elements with class 'expandable-view' and make sure they also have 'expanded-view-hidden'
    $(".expandable-view").each(function() {
        $(this).addClass("expanded-view-hidden");
    });
});

$("#btn_list_view_expanded").click(function(e) {
    $("#btn_arrangement_view").removeClass('active');
    $("#btn_list_view_compact").removeClass('active');
    $("#btn_list_view_expanded").addClass('active');

    $("#dataset_arrangement_c").hide();
    $("#dataset_list_c").show();
    
    // find all elements with class 'expandable-view' and make sure they also have 'expanded-view-hidden'
    $(".expandable-view").each(function() {
        $(this).removeClass("expanded-view-hidden");
    });
});

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
                $("#arrangement_save_notification").show();
                setTimeout(function() {
                    $("#arrangement_save_notification").fadeOut();
                    $("#btn_list_view_compact").trigger('click');
                }, 800);
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

$("#btn_set_primary_layout").click(function(e) {
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
            show_layout_action_note("Profile set as primary");
        },
        error: function (jqXHR, textStatus, errorThrown) {
            //alert("Failure!  Status: (" + textStatus + ") Error: (" + errorThrown + ")");
		    console.log('textStatus= ', textStatus);
		    console.log('errorThrown= ', errorThrown);
	        display_error_bar(jqXHR.status + ' ' + errorThrown.name);
        }
    }); //end ajax for search
});

$("#btn_share_layout").click(function(e) {
    var current_url = window.location.href;
    var current_page = current_url.lastIndexOf("dataset_explorer.html");
    var share_url = current_url.substring(0, current_page) + 'p?l=' + $('#selected_layout').find(':selected').data('share-id');

    if (copyToClipboard(share_url)) {
        show_layout_action_note("URL copied to clipboard");
    } else {
        show_layout_action_note("Failed to copy to clipboard. URL: " + share_url);
    }
});

$("#btn_view_layout_datasets").click(function(e) {
    $.ajax({
        url : './cgi/search_datasets.cgi',
        type: "POST",
        data : { 'session_id': session_id,
                 'layout_share_id': $('#selected_layout').find(':selected').data('share-id'),
                 'sort_by': $("#sort_by").val()
               },
        dataType:"json",
        success: function(data, textStatus, jqXHR) {
            process_search_results(data, ' in profile');
            current_dataset_list_label = $('#selected_layout').find(':selected').data('share-id');
        },
        error: function (jqXHR, textStatus, errorThrown) {
	        console.log('textStatus= ', textStatus);
	        console.log('errorThrown= ', errorThrown);
            display_error_bar(jqXHR.status + ' ' + errorThrown.name);
        }
    }); //end ajax for search 
});

function add_to_profile(dataset_id) {
    session_id = Cookies.get('gear_session_id');

    $.ajax({
        url : './cgi/add_dataset_to_layout.cgi',
        type: "POST",
        data : { 'session_id': session_id, 'layout_id': $('#selected_layout').val(), 'dataset_id': dataset_id },
        dataType:"json",
        success: function(data, textStatus, jqXHR) {
            if ( data['success'] ) {
                // switch the buttons
                $("#result_dataset_id_" + dataset_id + " button.add2profile").prop("disabled", true);
                $("#result_dataset_id_" + dataset_id + " button.removefromprofile").prop("disabled", false);
                show_dataset_action_note(dataset_id, "Dataset added to profile");
                current_layout_dataset_ids.push(dataset_id);
            } else {
                display_error_bar(data['error']);
            }
        },
        error: function (jqXHR, textStatus, errorThrown) {
            console.log('textStatus= ', textStatus);
            console.log('errorThrown= ', errorThrown);
            display_error_bar(jqXHR.status + ' ' + errorThrown.name);
        }
    });
}

function build_filter_string(group_name, att_name, crit) {
    // Builds a comma-separated search string based on the selected options
    //  in one of the filter option blocks
    if ($("#" + group_name + " ul li.selected").not(".all_selector").length) {
        var dbvals = [];

        $("#" + group_name + " ul li.selected").not(".all_selector").each(function() {
            dbvals.push($(this).data('dbval'));
        });

        crit[att_name] = dbvals.join(",");
    }
}

function copyToClipboard(text) {
    // https://stackoverflow.com/a/59594066
    if (window.clipboardData && window.clipboardData.setData) {
        // IE specific code path to prevent textarea being shown while dialog is visible.
        return clipboardData.setData("Text", text); 

    } else if (document.queryCommandSupported && document.queryCommandSupported("copy")) {
        var textarea = document.createElement("textarea");
        textarea.textContent = text;
        textarea.style.position = "fixed";  // Prevent scrolling to bottom of page in MS Edge.
        document.body.appendChild(textarea);
        textarea.select();
        try {
            return document.execCommand("copy");  // Security exception may be thrown by some browsers.
        } catch (ex) {
            console.warn("Copy to clipboard failed.", ex);
            return false;
        } finally {
            document.body.removeChild(textarea);
        }
    }
}

function display_error_bar(msg) {
    $('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
      '<button type="button" class="close close-alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
      '<p class="alert-message">' +
      '<strong>Fail. </strong> Sorry, something went wrong.  Please contact us with this message if you need help.' +
      '</p>' +
      '<p style="text-align: center;">(<em>Error: ' + msg + '</em>)</p>' +
      '</div>').show();
}

function load_initial_results() {
    params = { 'session_id': session_id,
               'sort_by': $("#sort_by").val()
             };

    var initial_search_terms = getUrlParameter('search_terms');
    if (initial_search_terms) {
        $("#search_clear").show();
        $("#search_terms").val(initial_search_terms);
        params['search_terms'] = initial_search_terms;
    } else {
        params['custom_list'] = 'most_recent';
    }
    
    $.ajax({
        url : './cgi/search_datasets.cgi',
        type: "POST",
        data : params,
        dataType:"json",
        success: function(data, textStatus, jqXHR) {
            process_search_results(data, ' most recent');
        },
        error: function (jqXHR, textStatus, errorThrown) {
	        console.log('textStatus= ', textStatus);
	        console.log('errorThrown= ', errorThrown);
            display_error_bar(jqXHR.status + ' ' + errorThrown.name);
        }
    }); //end ajax for search    
}

function load_organism_list() {
    session_id = Cookies.get('gear_session_id');
    
    $.ajax({
        url : './cgi/get_organism_list.cgi',
        type: "GET",
        data : {},
        dataType:"json",
        success: function(data, textStatus, jqXHR) {
            var ListTmpl = $.templates("#organism_list_tmpl");
            var ListHtml = ListTmpl.render(data['organisms']);
            $("#organism_choices").append(ListHtml);
        },
        error: function (jqXHR, textStatus, errorThrown) {
            display_error_bar(jqXHR.status + ' ' + errorThrown.name);
        }
    });
}

function load_preliminary_data() {
    /*
      Loads all the parts of the page which need initial calls from the server, such as 
      database-driven select boxes.
    */
    load_organism_list();
    load_user_layouts();
    load_initial_results();
}

function load_user_layouts() {
    session_id = Cookies.get('gear_session_id');

    if (session_id) {
        $.ajax({
            url : './cgi/get_user_layouts.cgi',
            type: "POST",
            data : { 'session_id': session_id, 'no_public': 0 },
            dataType:"json",
            success: function(data, textStatus, jqXHR) {
                // clear layout selector
                $('#site_profile_list').empty();
                $('#user_profile_list').empty();

                var site_profiles = [];
                var user_profiles = [];

                data['layouts'].forEach(function(p) {
                    if (p["is_domain"] == 1) {
                        site_profiles.push(p);
                    } else {
                        user_profiles.push(p);
                    }
                });

                var userListTmpl = $.templates("#layout_list_tmpl");
                var userListHtml = userListTmpl.render(user_profiles);
                $("#user_profile_list").html(userListHtml);

                var siteListTmpl = $.templates("#layout_list_tmpl");
                var siteListHtml = siteListTmpl.render(site_profiles);
                $("#site_profile_list").html(siteListHtml);

                $('#selected_layout').val(data['selected']);
                current_primary_layout = data['selected'];
                update_layout_dataset_id_list();
                update_view_buttons();
            },
            error: function (jqXHR, textStatus, errorThrown) {
                display_error_bar(jqXHR.status + ' ' + errorThrown.name);
            }
        });
    }
}

function process_search_results(data, result_label) {
    var domain_profile_selected = false;
    if ($("#selected_layout").find(':selected').data("is-domain") == "1") {
        domain_profile_selected = true;
    }
    
    for (const dataset of data['datasets']) {
        // make date_added nicer looking
        dataset['date_formatted'] = new Date(dataset['date_added']);
        dataset['date_formatted'] = dataset['date_formatted'].toDateString();

        var i;
        for (i = 0; i < dataset['layouts'].length; i++) {
            dataset['layouts'][i] = JSON.parse(dataset['layouts'][i]);
        }
        
        // some fields should never have whitespace
        if (dataset['pubmed_id'] !== 'None' && dataset['pubmed_id']) {
            dataset['pubmed_id'] = dataset['pubmed_id'].trim();
        } else {
            dataset['pubmed_id'] = 0;
        }

        if (dataset['geo_id'] !== 'None' && dataset['geo_id']) {
            dataset['geo_id'] = dataset['geo_id'].trim();
        } else {
            dataset['geo_id'] = 0;
        }

        if (dataset['is_public'] === null) {
            dataset['is_public'] = 0;
        }

        // set grid width styling for arrangement view
        if (dataset['grid_width'] == 12) {
            dataset['style_width'] = '588px';
        } else if (dataset['grid_width'] == 8) {
            dataset['style_width'] = '392px';
        } else { // grid_width == 4
            dataset['style_width'] = '196px';
        }
    } //end for

    // For the list view
    var resultsViewTmpl = $.templates("#dataset_results_view_tmpl");
    var resultsViewHtml = resultsViewTmpl.render(data['datasets']);
    $("#dataset_list_results_view_c").html(resultsViewHtml);

    // For the arrangement view
    var arrangementViewTmpl = $.templates("#dataset_arrangement_view_tmpl");
    var arrangementViewHtml = arrangementViewTmpl.render(data['datasets']);
    $("#dataset_arrangement").html(arrangementViewHtml);

    update_add_remove_buttons();

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

    // show the counts
    $("#result_count").html(data['datasets'].length);
    $("#result_label").html(result_label);

    if (mgmt_mode) {
        $(".profile_control").show();
    } else {
        $(".profile_control").hide();        
    }

    // Respect whatever the current view is
    if ($("#btn_list_view_expanded").hasClass('active')) {
        // find all elements with class 'expandable-view' and make sure they also have 'expanded-view-hidden'
        $(".expandable-view").each(function() {
            $(this).removeClass("expanded-view-hidden");
        });
    }

    $('.delete_dataset').popover({
		animation: true,
		trigger: 'click',
		title: "Delete dataset",
		content: "<p>Are you sure you want to delete this dataset?</p>" +
		    "<div class='btn-toolbar' style='width:250px'>" +
            "<span id='dataset_id_to_delete'></span>" +
		    "<button class='btn btn-default btn-danger confirm_dataset_delete' data-dismiss='popover'>Delete</button>" +
		    "<button id='cancel_dataset_delete' class='btn btn-default cancel_delete' value='cancel_delete'>Cancel</button>" +
		    "</div>",
		html: true,
		placement: 'auto',
        container: 'body'
    }).on('show.bs.popover', function(e) {
        // e.target is the popover trigger..
        dataset_id_to_delete = $(e.target).val();
    });
}

function remove_from_profile(dataset_id) {
    session_id = Cookies.get('gear_session_id');

    $.ajax({
        url : './cgi/remove_dataset_from_layout.cgi',
        type: "POST",
        data : { 'session_id': session_id, 'layout_id': $('#selected_layout').val(), 'dataset_id': dataset_id },
        dataType:"json",
        success: function(data, textStatus, jqXHR) {
            if ( data['success'] ) {
                // switch the buttons
                $("#result_dataset_id_" + dataset_id + " button.add2profile").prop("disabled", false);
                $("#result_dataset_id_" + dataset_id + " button.removefromprofile").prop("disabled", true);
                show_dataset_action_note(dataset_id, "Dataset removed from profile");
                
                // Javascript needs a better remove by value for arrays
                current_layout_dataset_ids = current_layout_dataset_ids.filter(val => val !== dataset_id);

                $("#result_dataset_id_" + dataset_id).fadeOut("slow", function() {
                    $("#result_count").html( $("#result_count").html() - 1  );
                    $(this).remove();
                });
            } else {
                display_error_bar(data['error']);
            }
        },
        error: function (jqXHR, textStatus, errorThrown) {
            console.log('textStatus= ', textStatus);
            console.log('errorThrown= ', errorThrown);
            display_error_bar(jqXHR.status + ' ' + errorThrown.name);
        }
    });
}

function show_dataset_action_note(dataset_id, msg) {
    var note_selector = "#result_dataset_id_" + dataset_id + " span.dataset_action_note";
    $(note_selector).html(msg).show();
    setTimeout(function() {
        $(note_selector).fadeOut().empty();
    }, 5000);
}

function show_layout_action_note(msg) {
    $("#layout_action_note").html(msg);
}

function submit_search() {
    // clear out any previous results:
    $("#dataset_list_results_view_c").empty();

    // If this is the first time searching with terms, set the sort by to relevance
    if ($("#search_terms").val() && first_search) {
        $("#sort_by").val('relevance');
        first_search = false;
    }
    
    var search_criteria = {
        'session_id': session_id,
        'search_terms': $("#search_terms").val(),
        'sort_by': $("#sort_by").val()
    };

    // collect the filter options the user defined
    build_filter_string('controls_organism', 'organism_ids', search_criteria);
    build_filter_string('controls_dtype', 'dtypes', search_criteria);
    build_filter_string('controls_date_added', 'date_added', search_criteria);
    build_filter_string('controls_ownership', 'ownership', search_criteria);

    $.ajax({
        url : './cgi/search_datasets.cgi',
        type: "POST",
        data : search_criteria,
        dataType:"json",
        success: function(data, textStatus, jqXHR) {
            process_search_results(data, ' results');
            current_dataset_list_label = 'search';
        },
        error: function (jqXHR, textStatus, errorThrown) {
	        console.log('textStatus= ', textStatus);
	        console.log('errorThrown= ', errorThrown);
            display_error_bar(jqXHR.status + ' ' + errorThrown.name);
        }
    }); //end ajax for search    
};

function update_add_remove_buttons() {
    // Iterates through each of the datasets in #dataset_list_results_view_c and updates
    //  the add to / remove from layout button and delete button
    $(".dataset_list_element_c").each(function() {
        var dataset_id = $(this).data("dataset-id");

        var domain_profile_selected = false;
        if ($("#selected_layout").find(':selected').data("is-domain") == "1") {
            domain_profile_selected = true;
        }

        // The ability to edit and delete and dataset are currently paired
        if (CURRENT_USER.id == $(this).find("button.delete_dataset").data('owner-id')) {
            $(this).find("button.delete_dataset").show();
            $(this).find("button.edit_dataset").show();
        } else {
            $(this).find("button.delete_dataset").hide();
            $(this).find("button.edit_dataset").hide();
        }
        
        if (domain_profile_selected) {
            $(this).find("button.add2profile").prop("disabled", true);
            $(this).find("button.removefromprofile").prop("disabled", true);
        } else {
            if (current_layout_dataset_ids.includes(dataset_id)) {
                $(this).find("button.add2profile").prop("disabled", true);
                $(this).find("button.removefromprofile").prop("disabled", false);
            } else {
                $(this).find("button.add2profile").prop("disabled", false);
                $(this).find("button.removefromprofile").prop("disabled", true);
            }
        }
    });
}

function update_layout_dataset_id_list() {
    $.ajax({
        url : './cgi/get_users_layout_members.cgi',
        type: "POST",
        data : { 'session_id': session_id, 'layout_id': $('#selected_layout').val() },
        dataType:"json",
        success: function(data, textStatus, jqXHR) {
            current_layout_dataset_ids = [];
            
            for (var lm of data['layout_members']) {
                lm = JSON.parse(lm);
                current_layout_dataset_ids.push(lm['dataset_id']);
            }

            update_add_remove_buttons();
        },
        error: function (jqXHR, textStatus, errorThrown) {
	        console.log('textStatus= ', textStatus);
	        console.log('errorThrown= ', errorThrown);
            display_error_bar(jqXHR.status + ' ' + errorThrown.name);
        }
    });
}

function update_view_buttons() {
    // Controls the on/off and visibility status of the compact/expanded/arrangement view buttons
    if ($('#selected_layout').find(':selected').data('is-domain') == "1") {
        $("#btn_arrangement_view").hide();
        
        // is arrangement the current view?  If so, set to compact
        if ($("#btn_arrangement_view").hasClass('active')) {
            $("#btn_list_view_compact").trigger('click');
        }
    } else {
        $("#btn_arrangement_view").show();
    }
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

