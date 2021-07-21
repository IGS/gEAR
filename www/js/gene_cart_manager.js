var first_search = true;
var animation_time = 200;

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

function load_preliminary_data() {
    /*
      Loads all the parts of the page which need initial calls from the server, such as 
      database-driven select boxes.
    */
    //load_user_layouts();
    //load_initial_results();
}

function process_search_results(data, result_label) {
    // each element keeps getting sent as a string
    for (i = 0; i < data['gene_carts'].length; i++) {
        data['gene_carts'][i] = JSON.parse(data['gene_carts'][i]);
    }
    
    // For the list view
    var resultsViewTmpl = $.templates("#gc_results_view_tmpl");
    var resultsViewHtml = resultsViewTmpl.render(data['gene_carts']);
    $("#gc_list_results_view_c").html(resultsViewHtml);
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
    build_filter_string('controls_date_added', 'date_added', search_criteria);
    build_filter_string('controls_ownership', 'ownership', search_criteria);

    $.ajax({
        url : './cgi/search_gene_carts.cgi',
        type: "POST",
        data : search_criteria,
        dataType:"json",
        success: function(data, textStatus, jqXHR) {
            console.log(data);
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
