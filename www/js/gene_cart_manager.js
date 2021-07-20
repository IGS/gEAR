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

function display_error_bar(msg) {
    $('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
      '<button type="button" class="close close-alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
      '<p class="alert-message">' +
      '<strong>Fail. </strong> Sorry, something went wrong.  Please contact us with this message if you need help.' +
      '</p>' +
      '<p style="text-align: center;">(<em>Error: ' + msg + '</em>)</p>' +
      '</div>').show();
}


