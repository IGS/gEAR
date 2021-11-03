var first_search = true;
var animation_time = 200;
var gc_id_to_delete = null;

window.onload=function() {
    // check if the user is already logged in
    check_for_login();

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

    $(document).on('blur', '#new_cart_label', function(e){
        if (! $(this).val()) {
            $(this).addClass("input-validation-error");
        } else {
            $(this).removeClass("input-validation-error");
        }
    });

    $(document).on('click', 'button.edit_gc', function() {
        var gc_id = $(this).data('gc-id');
        var selector_base = "#result_gc_id_" + gc_id;

        // copy the organism selection list for this row
        $("#result_gc_id_" + gc_id + "_editable_organism_id").html(
            $("#new_cart_organism_id").html()
        );

        // set the current value as selected
        $("#result_gc_id_" + gc_id + "_editable_organism_id").val(
            $("#result_gc_id_" + gc_id + "_editable_organism_id").data('original-val')
        );

        if ($("#result_gc_id_" + gc_id + "_editable_visibility").data('is-public')) {
            $("#result_gc_id_" + gc_id + "_editable_visibility").bootstrapToggle('on');
        } else {
            $("#result_gc_id_" + gc_id + "_editable_visibility").bootstrapToggle('off');
        }

        // Show editable versions where there are some and hide the display versions
        $(selector_base + " .is-editable").hide();
        $(selector_base + " .editable-version").show();

        // Make sure the view is expanded
        if ($(selector_base + " .expandable-view").hasClass('expanded-view-hidden')) {
            $(selector_base + " span.gc-expander").click();
        }
    });

    $(document).on('click', 'button.view_gc', function(e){
        window.location = "./p?c=" + $(this).val();
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
};  // end window onloads

// Popover for these is created within process_search_results()
$(document).on('click', '#cancel_gc_delete', function() {
    gc_id_to_delete = null;
    $('.delete_gc').popover('hide');
});

$(document).on('click', '.edit_gc_cancel', function() {
    var gc_id = $(this).data('gc-id');
    var selector_base = "#result_gc_id_" + gc_id;

    // Show editable versions where there are some and hide the display versions
    $(selector_base + " .editable-version").hide();
    $(selector_base + " .is-editable").show();

    // Reset any unsaved/edited values
    var visibility = $(selector_base + "_editable_visibility").data("original-val");
    $(selector_base + "_editable_visibility").val(visibility);

    var title = $(selector_base + "_editable_title").data("original-val");
    $(selector_base + "_editable_title").val(title);

    var org_id = $(selector_base + "_editable_organism_id").data("original-val");
    $(selector_base + "_editable_organism_id").val(org_id);
});

$(document).on('click', '.confirm_gc_delete', function() {
    $('.delete_gc').popover('hide');

    session_id = Cookies.get('gear_session_id');

    $.ajax({
        url : './cgi/remove_gene_cart.cgi',
        type: "POST",
        data : { 'session_id': session_id, 'gene_cart_id': gc_id_to_delete },
        dataType:"json",
        success: function(data, textStatus, jqXHR) {
            if (data['success'] == 1) {
                $("#result_gc_id_" + gc_id_to_delete).fadeOut("slow", function() {
                    $("#result_count").html( $("#result_count").html() - 1  );
                    $(this).remove();
                });
                gc_id_to_delete = null;
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

$(document).on('click', '.download_gc', function() {
    /*
      Reformats the <ul> containing the gene symbols into a text file with one gene
      per row.
     */
    var gc_id = $(this).data('gc-id');
    var file_contents = '';

    $("#" + gc_id + "_gene_list li").each(function(idx, li) {
        var gene_sym = $(li).html();
        file_contents += gene_sym + "\n";
    });

    var element = document.createElement("a");
    element.setAttribute(
        "href",
        "data:text/tab-separated-values;charset=utf-8," + encodeURIComponent(file_contents)
    );
    element.setAttribute("download", "gene_cart." + $(this).data("gc-share-id") + ".tsv");
    element.style.display = "none";
    document.body.appendChild(element);
    element.click();
    document.body.removeChild(element);
});

$(document).on('click', '.gc_unweighted_gene_list_toggle', function() {
    var gc_id = $(this).data('gc-id');
    var gene_list = $("#" + gc_id + "_gene_list");

    // see if the .gene_list is visible and toggle
    if (gene_list.is(":visible")) {
        gene_list.hide();
        $(this).addClass('btn-outline-secondary');
        $(this).removeClass('btn-secondary');
    } else {
        gene_list.show();
        $(this).removeClass('btn-outline-secondary');
        $(this).addClass('btn-secondary');
    }
});

$(document).on('click', '.gc_weighted_gene_list_toggle', function() {
    var gc_id = $(this).data('gc-id');
    var gene_list = $("#" + gc_id + "_gene_list");
    var share_id = $(this).data('gc-share-id');

    $("#btn_gc_" + gc_id + "_preview").hide();
    $("#btn_gc_" + gc_id + "_loading").show();

    $.ajax({
        url : './cgi/get_weighted_gene_cart_preview.cgi',
        type: "POST",
        data : {'share_id': share_id},
        dataType:"json",
        success: function(data, textStatus, jqXHR) {
            process_weighted_gc_list(gc_id, data['preview_json']);
        },
        error: function (jqXHR, textStatus, errorThrown) {
	        console.log('textStatus= ', textStatus);
	        console.log('errorThrown= ', errorThrown);
            display_error_bar(jqXHR.status + ' ' + errorThrown.name);
        }
    }); //end ajax
});

$(document).on('click', '.gc_weighted_gene_list_hider', function() {
    var gc_id = $(this).data('gc-id');
    $("#" + gc_id + '_gene_list').hide(animation_time);
    $("#btn_gc_" + gc_id + "_loading").hide();
    $("#btn_gc_" + gc_id + "_preview").show();
    $("#" + gc_id + "_gene_table").html('');
    $(this).hide();
});

$(document).on('click', 'button.share_gc', function() {
    share_id = $(this).attr('value');
    var current_url = window.location.href;
    var current_page = current_url.lastIndexOf("gene_cart_manager.html");
    var share_url = current_url.substring(0, current_page) + 'p?c=' + share_id;
    var gc_id = $(this).data('gc-id');

    if (copyToClipboard(share_url)) {
        show_gc_action_note(gc_id, "URL copied to clipboard");
    } else {
        show_gc_action_note(gc_id, "Failed to copy to clipboard. URL: " + share_url);
    }
});

$("#btn_create_cart_toggle").click(function(e) {
    if ($("#add_cart_panel").is(":visible")) {
        $("#add_cart_panel").hide();
        $("#gc_viewport").show(animation_time);
        $("#view_controls").show(animation_time);

        $("#btn_create_cart_toggle").html('Create new cart');
        reset_add_form();
    } else {
        $("#view_controls").hide();
        $("#gc_viewport").hide();
        $("#add_cart_panel").show(animation_time);
        $('#new_cart_is_public').bootstrapToggle('off');

        $("#btn_create_cart_toggle").html('Cancel cart creation');
    }
});

$("#btn_gc_paste_unweighted_list").click(function(e) {
    $("#new_cart_unweighted_header").addClass('bg-primary');
    $("#new_cart_unweighted_header").css('color', 'white');
    $("#btn_gc_upload_unweighted_list").addClass('disabled');
    $("#btn_gc_upload_weighted_list").addClass('disabled');
    $("#new_cart_pasted_genes_c").show();
    $("#new_cart_form_c").show(animation_time);
    $("#new_cart_upload_type").val('pasted_genes');
    $("#file_upload_c").hide();
});

$("#btn_gc_upload_unweighted_list").click(function(e) {
    $("#new_cart_unweighted_header").addClass('bg-primary');
    $("#new_cart_unweighted_header").css('color', 'white');
    $("#btn_gc_upload_weighted_list").addClass('disabled');
    $("#btn_gc_paste_unweighted_list").addClass('disabled');
    $("#new_cart_form_c").show(animation_time);
    $("#new_cart_upload_type").val('uploaded-unweighted');
    $("#file_upload_c").show();
});

$("#btn_gc_upload_weighted_list").click(function(e) {
    $("#new_cart_weighted_header").addClass('bg-primary');
    $("#new_cart_weighted_header").css('color', 'white');
    $("#btn_gc_upload_unweighted_list").addClass('disabled');
    $("#btn_gc_paste_unweighted_list").addClass('disabled');
    $("#new_cart_form_c").show(animation_time);
    $("#new_cart_upload_type").val('uploaded-weighted');
    $("#file_upload_c").show();
});

$("#btn_new_cart_cancel").click(function(e) {
    $("#btn_create_cart_toggle").trigger('click');
});

$('#new_cart_data').on('submit', function(e) {
    // disable button and show indicator that it's loading
    $("#btn_new_cart_save").hide();
    $("#btn_new_cart_saving").show();

    // check required fields
    if (! $("#new_cart_label").val()) {
        $("#new_cart_label").addClass("input-validation-error");
        $("#btn_new_cart_save").show();
        $("#btn_new_cart_saving").hide();
        return false;
    } else {
        $("#new_cart_label").removeClass("input-validation-error");
    }
    
    session_id = Cookies.get('gear_session_id');
    var is_public = ($("#new_cart_is_public").prop('checked') ? 1 : 0);

    var formData = new FormData($(this)[0]);
    formData.append('is_public', is_public);
    formData.append('session_id', session_id);

    // https://stackoverflow.com/questions/45594504/upload-file-and-json-data-in-the-same-post-request-using-jquery-ajax
    var gc = new GeneCart();
    gc.add_cart_to_db_from_form(gene_cart_saved, formData);

    return false;
});

$("#btn_new_cart_save").click(function(e) {
    $('#new_cart_data').trigger("submit");
});


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

function display_error_bar(msg) {
    $('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
      '<button type="button" class="close close-alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
      '<p class="alert-message">' +
      '<strong>Fail. </strong> Sorry, something went wrong.  Please contact us with this message if you need help.' +
      '</p>' +
      '<p style="text-align: center;">(<em>Error: ' + msg + '</em>)</p>' +
      '</div>').show();
}

function gene_cart_saved() {
    $("#btn_create_cart_toggle").trigger('click');
    submit_search();
    reset_add_form();
}

function load_preliminary_data() {
    /*
      Loads all the parts of the page which need initial calls from the server, such as 
      database-driven select boxes.
    */
    load_organism_list()
    $("#your_gene_cart_filter").trigger('click');
}

function load_organism_list() {
    $.ajax({
        url : './cgi/get_organism_list.cgi',
        type: "GET",
        data : {},
        dataType:"json",
        success: function(data, textStatus, jqXHR) {
            var ListTmpl = $.templates("#organism_list_tmpl");
            var ListHtml = ListTmpl.render(data['organisms']);
            $("#organism_choices").append(ListHtml);

            var selectTmpl = $.templates("#organism_select_tmpl");
            var selectHtml = selectTmpl.render(data['organisms']);
            $("#new_cart_organism_id").append(selectHtml);
        },
        error: function (jqXHR, textStatus, errorThrown) {
            display_error_bar(jqXHR.status + ' ' + errorThrown.name);
        }
    });
}

function process_search_results(data, result_label) {
    // each element keeps getting sent as a string
    for (i = 0; i < data['gene_carts'].length; i++) {
        data['gene_carts'][i] = JSON.parse(data['gene_carts'][i]);
        data['gene_carts'][i]['date_added'] = new Date(data['gene_carts'][i]['date_added']);
        data['gene_carts'][i]['date_added'] = data['gene_carts'][i]['date_added'].toDateString();
    }

    var resultsViewTmpl = $.templates("#gc_results_view_tmpl");
    var resultsViewHtml = resultsViewTmpl.render(data['gene_carts']);
    $("#gc_list_results_view_c").html(resultsViewHtml);

    update_gene_cart_list_buttons();

    $("#result_count").html(data['gene_carts'].length);
    $("#result_label").html(result_label);

    $('.delete_gc').popover({
		animation: true,
		trigger: 'click',
		title: "Delete gene cart",
		content: "<p>Are you sure you want to delete this gene cart?</p>" +
		    "<div class='btn-toolbar' style='width:250px'>" +
            "<span id='gc_id_to_delete'></span>" +
		    "<button class='btn btn-default btn-danger confirm_gc_delete' data-dismiss='popover'>Delete</button>" +
		    "<button id='cancel_gc_delete' class='btn btn-default cancel_delete' value='cancel_delete'>Cancel</button>" +
		    "</div>",
		html: true,
		placement: 'auto',
        container: 'body',
        sanitize: false
    }).on('show.bs.popover', function(e) {
        // e.target is the popover trigger..
        gc_id_to_delete = $(e.target).val();
    });
}

function process_weighted_gc_list(gc_id, jdata) {
    $("#btn_gc_" + gc_id + "_loading").hide();
    $("#btn_gc_" + gc_id + "_hider").show();

    // This creates a table with classes dataframe and weighted-list
    $("#" + gc_id + '_gene_table').html(jdata);
    $("#" + gc_id + '_gene_table').show(animation_time);
}

function reset_add_form() {
    $("#btn_new_cart_saving").hide();
    $("#btn_new_cart_save").show();
    
    $("#new_cart_label").val('');
    $("#new_cart_ldesc").val('');
    $("#new_cart_pasted_genes").val('');
    $("#new_cart_file").val('');
    $('#new_cart_is_public').bootstrapToggle('off');

    $("#new_cart_unweighted_header").removeClass('bg-primary');
    $("#new_cart_unweighted_header").css('color', 'black');

    $("#new_cart_weighted_header").removeClass('bg-primary');
    $("#new_cart_weighted_header").css('color', 'black');
    
    $("#btn_gc_paste_unweighted_list").removeClass('disabled');
    $("#btn_gc_upload_unweighted_list").removeClass('disabled');
    $("#btn_gc_upload_weighted_list").removeClass('disabled');

    $("#new_cart_form_c").hide();
    $("#new_cart_pasted_genes_c").hide();
}

function show_gc_action_note(gc_id, msg) {
    var note_selector = "#result_gc_id_" + gc_id + " span.gc_action_note";
    $(note_selector).html(msg).show();
    setTimeout(function() {
        $(note_selector).fadeOut().empty();
    }, 5000);
}

function submit_search() {
    // clear out any previous results:
    $("#gc_list_results_view_c").empty();

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
    build_filter_string('controls_date_added', 'date_added', search_criteria);
    build_filter_string('controls_ownership', 'ownership', search_criteria);

    $.ajax({
        url : './cgi/search_gene_carts.cgi',
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

function update_gene_cart_list_buttons() {
    // Iterates through each of the carts in the result list and makes sure the appropriate
    //  buttons are visible/hidden
    $(".gc_list_element_c").each(function() {
        var gc_id = $(this).data("gc-id");

        // The ability to edit and delete and dataset are currently paired
        if (CURRENT_USER.id == $(this).find("button.delete_gc").data('owner-id')) {
            $(this).find("button.delete_gc").show();
            $(this).find("button.edit_gc").show();
        } else {
            $(this).find("button.delete_gc").hide();
            $(this).find("button.edit_gc").hide();
        }
    });
}
