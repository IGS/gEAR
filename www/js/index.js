var search_results = [];
var layouts = []; //populated on load by loads_layouts();
var gene_carts = [];

// key dataset_id, value = Snap paths
var svgs = {};

var SCROLLBAR_DRAWN = false;
var GO_TERM_SCROLLBAR_DRAWN = false;
var AT_FIRST_MATCH_RECORD = false;
var AT_LAST_MATCH_RECORD  = false;
var PREVIOUS_SELECTED_RECORD_NUM = null;
var SCORING_METHOD = 'gene';
var SELECTED_GENE = null;

var share_id = null; //from permalink
var permalinked_dataset_id = null; //holds dataset_id obtained from load_dataset_frames()

var annotation_panel = new FunctionalAnnotationPanel();
var dataset_collection_panel = new DatasetCollectionPanel();

var search_result_postselection_functions = [];

window.onload=function() {
    // check if the user is already logged in
    check_for_login();

    // Was a permalink found?
    share_id = getUrlParameter('share_id');
    scope = "permalink";
    if (!share_id) {
        // layout_id is a share_id for the profile layout
        share_id = getUrlParameter('layout_id');
        scope = "profile";
    }
    if (share_id) {
        //hide site_into and display the permalink message
        $('#intro_content').hide();
        $('#viewport_intro').children().hide();
        $('#searching_indicator_c').hide();

        $('#leftbar_main').show();
        $('#permalink_intro_c').show();
        // validate the share_id. runs load_dataset_frames() on success
        validate_permalink(share_id, scope);
    } else {
        get_index_info();
        load_gene_carts();
    }

    // Was help_id found?
    var help_id = getUrlParameter('help_id');
    if (help_id) {
        validate_help_id(help_id);
    }

    var permalinked_gene_symbol = getUrlParameter('gene_symbol');
    var permalinked_gsem = getUrlParameter('gene_symbol_exact_match');
    var permalinked_multigene_plots = getUrlParameter('multigene_plots');

    if (permalinked_gene_symbol) {
        $("#search_gene_symbol_intro").val(permalinked_gene_symbol);

        if (permalinked_gsem && permalinked_gsem === "1") {
            set_exact_match('on');
        }

        if (permalinked_multigene_plots && permalinked_multigene_plots === "1") {
            set_multigene_plots('on');
        }

        sleep(1000).then(() => {
            $('#intro_search_icon').trigger('click');
            // clear any open tooltips
            $('[data-toggle="tooltip"], .tooltip').tooltip("hide");
        })
    }

    // The search button starts out disabled, make sure it gets re-enabled.
    $("button#submit_search").prop( "disabled", false );

    $('#exact_match_icon').click(function() {
        // handle if it was already in the on position
        if ( $('#exact_match_input').prop("checked") ) {
            set_exact_match('off');
        } else {
            // else it was off, so turn it on
            set_exact_match('on');
        }
    });


    $('#multigene_plots_input').change(function() {
        if ( $('#multigene_plots_input').prop("checked") ) {
            $('#search_results_c').hide();
        } else {
            // else it was off, so turn it on
            $('#search_results_c').show();
        }
    });

    $('#intro_search_form').on('submit', function(e) {
        // TODO: It makes sense to remove/destroy those elements we aren't showing after a search
        e.preventDefault();
        $("#search_gene_symbol").val( $("#search_gene_symbol_intro").val());
        $('#intro_content').hide();

        $("#leftbar_main").show();
        $("#viewport_main").show();

        // fire the true search button
        $("#submit_search").trigger( "click" );
    });

    $('#intro_search_icon').click(function() {
        $('#intro_search_form').submit();
    });

    $('#dataset_search_form').on('submit', function(e) {
        e.preventDefault();
        window.location.replace("./dataset_explorer.html?search_terms=" + encodeURI($('#search_dataset_intro').val()));
    });

    $('#launcher_manual').click(function() {
        window.location.replace('./manual.html');
    });

    $('#launcher_expression_uploader').click(function() {
        window.location.replace('./upload_dataset.html');
    });

    $('#launcher_epigenetic_uploader').click(function() {
        window.location.replace('./upload_epigenetic_data.html');
    });

    // add post-page load listeners
    $( "#dataset_zoomed_zoom_out_control" ).click(function() {
        zoom_out_dataset();
    });

    $(document).on('click', '.domain_choice_c', function() {
        dataset_collection_panel.set_layout($(this).data('profile-id'), $(this).data('profile-label'), true);
        share_id = $(this).data('profile-share-id');
    });

    $( document ).on("click", ".scope_choice", function() {
        SCORING_METHOD = $(this).data('choice');
        if (SELECTED_GENE !== null) {
            select_search_result($(SELECTED_GENE));
        }
    });

    // track the mouse movement so we can display scoring tooltips
    $( document ).on( "mousemove", function( event ) {
        // Positioning for dataset_grid tips
        // Why is this pixel adjustment necessary?
        xpos = event.pageX - 240;
        ypos = event.pageY - 130 - 30;
        $("#tip").css("left", xpos + "px" );
        $("#tip").css("top" , ypos + "px" );
    });
};

function get_index_info() {
    $.ajax({
        url: './cgi/get_index_info.cgi',
        type: 'GET',
        dataType: 'json',
        success: function(data, textStatus, jqXHR) {

            $('#stats_dataset_count').text(data['dataset_count'])
            $('#stats_user_count').text(data['user_count'])
        },
        error: function(jqXHR, textStatus, errorThrown) {
            display_error_bar(jqXHR.status + ' ' + errorThrown.name, 'Error getting index info.');
        }
    });
};

//Check help_id is valid. For Forgotten Password
function validate_help_id(help_id) {
    $.ajax({
        url: './cgi/validate_help_id.cgi',
        type: 'POST',
        data: {'help_id': help_id},
        dataType: 'json',
        success: function(data, textStatus, jqXHR) {
            if ( data['success'] == 1 ) {
                // Add help_id to form
                $('#user_help_id').val(help_id);

                // Greet user by name (a subtle confirmation they know it's their account)
                if (data['user_name'].length > 0) {
                    var user_first_name = ' ' + data['user_name'].split(' ')[0] + '!';
                    $('#forgot_password_user_name').text(user_first_name);
                }

                $('#forgot_password_modal').modal('show');
            } else {
                // Invalid help_id, display invalid message
                $("#valid_forgot_pass_modal_body_c").hide();
                $("#save_user_new_pass").hide();
                $("#invalid_forgot_pass_modal_body_c").show();

                $('#forgot_password_modal').modal('show');
            }
        },
        error: function(jqXHR, textStatus, errorThrown) {
            display_error_bar(jqXHR.status + ' ' + errorThrown.name, 'Error validating help ID');
        }
    });
};

// Disable 2nd password input until 1st is populated
$(document).on('keydown', 'input#user_new_pass_1', function(){
    if ( $(this).val().length > 1 ) {
        $("input#user_new_pass_2").prop("disabled", false);
    } else {
        $("input#user_new_pass_2").prop("disabled", true);
    }
});

// Disable Save password button until 1st and 2nd inputs match
$(document).on('keyup', 'input#user_new_pass_2', function(){
    var pass_1 = $('input#user_new_pass_1').val();

    if ( $('input#user_new_pass_2').val() == pass_1 ) {
        $('button#save_user_new_pass').prop('disabled', false);
    } else {
        $('button#save_user_new_pass').prop('disabled', true);
    }
});

// Submit new password
$(document).on('click', 'button#save_user_new_pass', function(){
    // Hide password form and show waiting
    $('#valid_forgot_pass_modal_body_c').hide();
    $('#forgot_pass_modal_footer').hide();
    $('#saving_forgot_pass_modal_body_c').show();

    var help_id = $('input#user_help_id').val();
    var new_password = $('input#user_new_pass_2').val();
    $.ajax({
        url: './cgi/save_user_account_changes.cgi',
        type: 'POST',
        data: { 'help_id': help_id, 'new_password': new_password, 'scope': 'password'},
        dataType: 'json',
        success: function(data, textStatus, jqXHR) {
            if (data['success'] == 1) {
                // Hide waiting and show success
                $('#saving_forgot_pass_modal_body_c').hide();
                $('#success_forgot_pass_modal_body_c').show();

                // Redirect to home page
                setInterval(function(){
                    window.location.replace('./index.html');
                }, 2000);
            } else {
                $('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
                    '<button type="button" class="close close-alert" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
                    '<p class="alert-message"><strong>Oops! </strong> ' + data["error"] + '</p></div>').show();
            }
        },
        error: function (jqXHR, textStatus, errorThrown) {
            display_error_bar(jqXHR.status + ' ' + errorThrown.name, 'Error saving new user password');
        }
    });//end ajax
});

function validate_permalink(share_id, scope) {
    $.ajax({
        url : './cgi/validate_share_id.cgi',
        type: "POST",
        data : { 'share_id': share_id, 'scope': scope },
        dataType:"json",
        success: function(data, textStatus, jqXHR) {
            if ( data['success'] == 1 ) {
                // query the db and load the images, including permalink dataset
                dataset_collection_panel.load_frames({ share_id });

            } else {
                // query the db and load the images
                dataset_collection_panel.load_frames();
                $('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
                    '<button type="button" class="close close-alert" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
                    '<p class="alert-message"><strong>Oops! </strong> ' + data["error"] + '</p></div>').show();
            }
        },
        error: function (jqXHR, textStatus, errorThrown) {
            display_error_bar(jqXHR.status + ' ' + errorThrown.name, 'Error validating share ID');
        }
    });
}

function load_layouts() {
    var d = new $.Deferred();
    var session_id = Cookies.get('gear_session_id');
    var layout_share_id = getUrlParameter('layout_id');

    // Temporary hack for Heller lab
    if (layout_share_id == '8d38b600' || layout_share_id == 'afd2eb77') {
        $("#intro_selected_profile_warning").show();
    }

    //organize user and domain profiles in x-editable format
    $.ajax({
        url: './cgi/get_user_layouts.cgi',
        type: 'post',
        data: { 'session_id': session_id, 'layout_share_id': layout_share_id },
        dataType: 'json',
        success: function(data, textStatus, jqXHR) {
            /*
              Priority of displayed profile:
                0.  Passed layout ID via layout_id URL parameter
                1.  Cookie value
                2.  User's DB-saved value (when they go to a machine, and there's no cookie)
                3.  Admin's active domain
             */

            var domain_profiles = [];
            var user_profiles = [];

            var active_layout_id = null;
            var active_layout_label = null;

            // Pass through once to sort domains from user profiles AND see if it matches a shared layout
            $.each(data['layouts'], function(i, item){
                if ( item['is_domain'] == 1 ) {
                    domain_profiles.push({value: item['id'], text: item['label'], share_id: item['share_id'] });
                } else {
                    user_profiles.push({value: item['id'], text: item['label'], share_id: item['share_id']  });
                }

                if (item['share_id'] == layout_share_id) {
                    active_layout_id = item['id'];
                    active_layout_label = item['label'];
                    share_id = item['share_id'];
                }
            });

            // pass through again and look for one set by a cookie
            if (active_layout_id == null) {
                $.each(data['layouts'], function(i, item) {
                    if (item['label'] == CURRENT_USER.profile) {
                        active_layout_id = item['id'];
                        active_layout_label = item['label'];
                        share_id = item['share_id'];
                        return false;
                    }
                });
            }

            // pass through again and look for one set as current by the user
            if (active_layout_id == null) {
                $.each(data['layouts'], function(i, item) {
                    if ( item['is_domain'] == 0 && item['is_current'] == 1 ) {
                        active_layout_id = item['id'];
                        active_layout_label = item['label'];
                        share_id = item['share_id'];
                        return false;
                    }
                });
            }

            // pass through again if no active layout was found for user and choose the admin's
            if (active_layout_id == null) {
                $.each(data['layouts'], function(i, item) {
                    if ( item['is_domain'] == 1 && item['is_current'] == 1 ) {
                        active_layout_id = item['id'];
                        active_layout_label = item['label'];
                        share_id = item['share_id'];
                        return false;
                    }
                });
            }

            var formattedData = [{text: 'Public profiles', children: domain_profiles}];

            if (user_profiles.length > 0) {
                formattedData.push({text: 'Your profiles', children: user_profiles});
            }

            var layout_items_tmpl = $.templates("#tmpl_intro_search_layouts");

            var domain_items_html = layout_items_tmpl.render(domain_profiles);
            $(domain_items_html).insertAfter('#intro_domains_header');

            var profile_items_html = layout_items_tmpl.render(user_profiles);
            $(profile_items_html).insertAfter('#intro_profiles_header');

            //Serves as source for #selected_profile editable
            layouts = formattedData;

            dataset_collection_panel.set_layout(active_layout_id, active_layout_label, true);

            d.resolve();
        },
        error: function (jqXHR, textStatus, errorThrown) {
            display_error_bar(jqXHR.status + ' ' + errorThrown.name, 'Error loading layouts.');
            d.fail();
        }
    });

    d.promise();
}

function load_gene_carts() {
  var d = new $.Deferred();
  var session_id = Cookies.get('gear_session_id');

  if (!session_id) {
      //User is not logged in. Hide gene carts container
      $("#selected_gene_cart_c").hide();
      d.resolve();
  } else {
      $("#selected_gene_cart_c").show(); //Show if hidden
      $.ajax({
        url: './cgi/get_user_gene_carts.cgi',
        type: 'post',
        data: { 'session_id': session_id },
        dataType: 'json',
        success: function(data, textStatus, jqXHR){ //source https://stackoverflow.com/a/20915207/2900840
            var user_gene_carts = [];
            var formattedData = [];

            if (data['gene_carts'].length > 0) {
                //User has some profiles
                $.each(data['gene_carts'], function(i, item){
                    user_gene_carts.push({value: item['id'], text: item['label'] });

                });

                formattedData = [{text: 'My gene carts', children: user_gene_carts }];
            } else {
                $("#selected_gene_cart_c").hide();
            }

            //Serves as source for #selected_gene_cart editable
            gene_carts = formattedData;
            d.resolve();
        },
        error: function (jqXHR, textStatus, errorThrown) {
            display_error_bar(jqXHR.status + ' ' + errorThrown.name);
            d.fail();
        }
      });
  }
  d.promise();
}

// Hide option menu when scope is changed.
$(document).on('click', '.scope_choice', function(){
    $('#toggle_options').popover('hide');
});

// Handle direct documentation links
$(document).on('click', '#doc-link-choices li', function(){
    window.location.replace("./manual.html?doc=" + $(this).data('doc-link'));
});

function populate_search_result_list(data) {
    // so we can display in sorted order.  javascript sucks like that.
    sorted_gene_syms = [];

    for (var key in data) {
        if (data.hasOwnProperty(key)) {
            sorted_gene_syms.push(key);
        }
    }

    sorted_gene_syms.sort();
    sorted_gene_syms_len = sorted_gene_syms.length

    var items = [];

    for (i = 0; i < sorted_gene_syms_len; i++) {
        gene_symbol = sorted_gene_syms[i];

        // Build search result html
        var gene_result_html = '<a class="list-group-item" data-gene_symbol="' + gene_symbol + '" href="#">' + gene_symbol;

        gene_result_html += '</a>';
        items.push(gene_result_html);
    }

    if (items.length == 0) {
        $('#search_results').text('No results found');
        $('#search_result_count').text(items.length);
    } else {
        $('#search_results').append( items.join('') );

        // the value here needs to match the max in gene_search.cgi
        if (items.length == 100) {
            $('#search_result_count').text('max:' + items.length);
        } else {
            $('#search_result_count').text(items.length);
        }
    }
}

var lastCall = 0;
function select_search_result(elm) {
    //TODO Prevents this function from being double-called by #gene_search_form.submit()
    var callTime = new Date().getTime();
    if (callTime - lastCall <= 500) {
      return false;
    }
    lastCall = callTime;

    SELECTED_GENE = $(elm);
    gene_sym = $(elm).data("gene_symbol");

    // remove coloring from other result links
    $('.list-group-item-active').removeClass('list-group-item-active');
    $(elm).addClass('list-group-item-active');

    var current_record_number = 0
    var match_record_number = null

    annotation_panel.annotation = search_results[gene_sym];
    annotation_panel.display_first_organism();

    // hide the intro, show the search result box
    if( $('#site_intro_c').is(':visible') ) {
        $('#site_intro_c').hide({easing: 'fade', duration: 400});
        $('#recent_updates_c').hide({easing: 'fade', duration: 400, complete: show_search_result_info_box});
    }

    dataset_collection_panel.update_by_search_result(search_results[gene_sym]);

    // call any plugin functions
    search_result_postselection_functions.forEach(function(f) {f()})
}

function isNumeric(n) {
  return !isNaN(parseFloat(n)) && isFinite(n);
}

function show_search_result_info_box() {
    if( permalinked_dataset_id ) {
        // show links_out and gene_annot with zoom_on
        $('#links_out_c, #gene_details_c').addClass('search_result_c').removeClass('search_result_c_DISABLED').show('fade', {}, 400);
        $('#dataset_zoomed_c').show('fade', {}, 400);
    } else {
        // $('div.search_result_c').show('fade', {}, 400);
        $('.search_result_c_DISABLED').addClass('search_result_c').removeClass('search_result_c_DISABLED');
        // $('div.search_result_c').toggleClass('search_result_c_DISABLED');
    }
}

function find_share_id_by_pk(pk) {
    // Given the layout primary key, retrieve the share_id from the array of previously generated layouts

    // layouts -> public/user profiles -> children -> indiv. layouts -> value/text/share_id

    for (i = 0; i < layouts.length; i++) {
        let found_idx = layouts[i].children.findIndex(function(v) {
            return v.value === pk;
        })
        if (found_idx > 0)
            return layouts[i].children[found_idx].share_id;
    }

    // Failsafe if primary key was not found (though this should not happen)
    console.log("Cannot find share_id for layout with id " + pk + " in list of layouts");
    return "undefined";
}

$('#search_results').on("click", "a", function(e) {
    e.preventDefault(); //prevent page scrolling to top
    $(this).blur(); //removes focus so active's purple coloring can show
    select_search_result(this);
});

// Warn user if no datasets in profile
$( "#search_gene_symbol").focus(function(){
    if (dataset_collection_panel.datasets.length < 1) {
        $("#search_gene_symbol").popover('show');
    } else {
        $("#search_gene_symbol").popover('hide');
    }
});
$("#search_gene_symbol").blur(function() {
    $("#search_gene_symbol").popover('hide');
});

// Popover for warning user that profile lacks datasets.
$('#search_gene_symbol').popover({
  	animation: true,
  	trigger: 'manual',
    container: 'body',
  	content: "<div class='text-center' style='width:250px;'>" +
        "<div class='alert alert-warning text-center'>" +
        "<p><span class='fa fa-exclamation'></span> <b>No datasets in current profile</b></p>" +
      	"</div>" +
        "<p>To search a gene, add a dataset to your layout profile in the <a href='./dataset_manager.html'>Dataset Manager</a>.</p>" +
        "</div>",
  	html: true,
  	placement: 'right'
});

$("#gene_search_form").submit(function( event ) {
    $("#viewport_intro").hide();
    $("#viewport_main").show();

    // determine if searching for exact matches (convert from bool to 1 or 0)
    $("#exact_match").val( Number($('#exact_match_input').prop("checked")) );
    $("#multigene_plots").val( Number($('#multigene_plots_input').prop("checked")) );


    $('#recent_updates_c').hide();
    $('#searching_indicator_c').show();

    var formData = $("#gene_search_form").serializeArray();
    //console.log(formData);

    // split on combination of space and comma (individually or both together.)
    var gene_symbol_array = $("#search_gene_symbol").val().split(/[\s,]+/);
    // Remove duplicates in gene search if they exist
    var uniq_gene_symbols = gene_symbol_array.filter((value, index, self) => self.indexOf(value) === index);
    var curated_searched_gene_symbols = uniq_gene_symbols.join(',');

    history.pushState(
        // State Info
        {
            'layout_id': share_id,
            'gene_symbol': curated_searched_gene_symbols,
            'gene_symbol_exact_match': $("#exact_match").val(),
            'multigene_plots': $("#multigene_plots").val()
        },
        // State title
        "Gene search",
        // URL
        "/index.html?layout_id=" + share_id
            + "&gene_symbol=" + encodeURIComponent(curated_searched_gene_symbols)
            + "&gene_symbol_exact_match=" + $("#exact_match").val()
            + "&multigene_plots=" + $("#multigene_plots").val()
    )

    $('#search_results').empty();
    // show search results
    $('#search_results_c').removeClass('search_result_c_DISABLED');

    $.ajax({
        url : './cgi/search_genes.py',
        type: "POST",
        data : formData,
        dataType:"json",
        success: function(data, textStatus, jqXHR) {
        	// reset search_results
        	search_results = data;

            populate_search_result_list(data);
            $('#searching_indicator_c').hide();
            $('#intro_content').hide('fade', {}, 400, function() {
                // auto-select the first match.  first <a class="list-group-item"
                first_thing = $('#search_results a.list-group-item').first();
                select_search_result(first_thing);
            });

            // http://manos.malihu.gr/jquery-custom-content-scroller/
            // The author of this utility was wonderfully responsive and helpful
            if (SCROLLBAR_DRAWN == false) {
                $("#search_results_scrollbox").mCustomScrollbar({
                	theme: '3d-thick-dark',
                    scrollButtons:{ enable:true },
                    // we need to disable the keyboard scrolling so our custom indicators can work
                    keyboard:{ enable:false }
                });

                // Change height after initializing scroller.
                $("#search_results_scrollbox").css({'height': 'calc(90vh - 260px)'});

                SCROLLBAR_DRAWN = true;
            } else {
                $("#search_results_scrollbox").mCustomScrollbar("update");
            }

            return false;
        },
        error: function (jqXHR, textStatus, errorThrown) {
            $('#searching_indicator_c').hide();

            // Error occurred
      			if ( $('#search_gene_symbol').val().length < 1 ) {
                // No gene symbol entered
        				$('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
        					'<button type="button" class="close close-alert" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
        					'<p class="alert-message"><strong>Oops! </strong>No gene symbol was entered. Enter a gene symbol and try again.</p></div>').show();

      			} else if ( dataset_ids_loaded.length == 0) {
                // No datasets in current layout profile
        				$('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
        					'<button type="button" class="close close-alert" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
        					'<p class="alert-message"><strong>Oops! </strong>No datasets were found in the current layout profile.</p><p>To add datasets to a profile or choose a different profile, go to the <a href="./dataset_manager.html" class="alert-link">Dataset Manager</a>.</p></div>').show();

            } else {
                // Some other error occurred
              	display_error_bar(jqXHR.status + ' ' + errorThrown.name);
      		}
        }
    });

    return false;  // keeps the page from not refreshing

});

// controls to enable user scrolling of results with mouse arrow
scrolling_results = false

$('body').click(function(event) {
    if (!$(event.target).closest('#search_results_c').length) {
        scrolling_results = false
    } else {
        scrolling_results = true
    };
});

$(document).keydown(function(event) {
    // don't do anything unless we're scrolling results
    if (scrolling_results == true) {
        // this makes sure the browser doesn't scroll the window
        event.preventDefault();

        switch (event.keyCode) {
            // up key
            case 38:
            if (AT_FIRST_MATCH_RECORD == false) {
                select_search_result($(SELECTED_GENE).prev())
            }
            break;

            // down key
            case 40:
            if (AT_LAST_MATCH_RECORD == false) {
                select_search_result($(SELECTED_GENE).next())
            }
            break;
        }
    }
});

// Gene Details collapse css changes
$(document).on('click', '#gene_details_header, #gene_collapse_btn', function() {
    if ($('#gene_collapse_btn').attr('aria-expanded') == 'false') {
        //Details is collapsed
        $('#gene_details_header').css({
            'border-bottom-left-radius': '4px',
            'border-bottom-right-radius': '4px'
        });

        //Change '-' back to '+' button
        $('#gene_collapse_btn').replaceWith('<span id="gene_collapse_btn" class="fa ' +
            'fa-minus pull-right" title="Show gene information" ' +
            'data-toggle="collapse"' +
            'data-target="#gene_details_info" aria-expanded="true" ' +
            'aria-controls="gene_details_info"></span>');
    } else {
      //Details is expanded
      $('#gene_details_header').css({
          'border-bottom-left-radius': '0px',
          'border-bottom-right-radius': '0px'
      });

      //Change '+' to '-'
      $('#gene_collapse_btn').replaceWith('<span id="gene_collapse_btn" class="fa ' +
          'fa-plus pull-right" title="Hide gene information" ' +
          'data-toggle="collapse" ' +
          'data-target="#gene_details_info" aria-expanded="false" ' +
          'aria-controls="gene_details_info"></span>');
    }
});

if (window.location.href.indexOf("manual.html") === -1) { //Without this the manual links break.
    // Allows user to change layout profile without having to go to dataset Manager
    // User not logged in can only change domain profiles
    $('#selected_profile').editable({
        mode: 'inline',
        type: 'select',
        pk: function() {
            //return the layout_id
            return $('.editable-input select option:selected').attr('value');
        },
        params: function(params) {
            //data to be sent to 'url'
            var data = {};
            data['layout_id'] = parseInt(params.pk);
            data['session_id'] =  Cookies.get('gear_session_id');
            return data;
        },
        url: function(params) {
            /*
              When a user changes the profile the following needs to happen

              - Store the new layout ID
              - set share_id so new state can be pushed to history (in update_datasetframes_generesults)
              - Fetch and draw the new DatasetCollectionPanel
              - If a user is logged in, set a cookie with the new layout ID
             */
            dataset_collection_panel.set_layout(params.layout_id, $('.editable-input select option:selected').text(), true);
            share_id = find_share_id_by_pk(params.layout_id)
            update_datasetframes_generesults();
        },
        title: 'Click to change',
        tpl: '<select style="width:160px"></select>',
        showbuttons: false,
        source: function() {
            return layouts; //populated by load_layouts()
        },
        sourceCache: true
    });

    // Load user's gene carts
    $('#selected_gene_cart').editable({
        mode: 'inline',
        type: 'select',
        pk: function() {
            //return the gene_cart_id
            return $('.editable-input select option:selected').attr('value');
        },
        params: function(params) {
            //data to be sent to 'url'
            var data = {};
            data['gene_cart_id'] = parseInt(params.pk);
            data['session_id'] =  Cookies.get('gear_session_id');
            return data;
        },
        url: function(params) {
            var d = new $.Deferred(); //Causes editable to wait until results are returned

            //User is not logged in
            if (!params.session_id) {
                d.resolve();
            } else {
                //User is logged in
                $("#search_gene_symbol").prop("disabled", true);
                $("#selected_gene_cart_loading_c").show();

                //Get the gene cart members and populate the gene symbol search bar
                $.ajax({
                    url: './cgi/get_gene_cart_members.cgi',
                    type: 'post',
                    data: params,
                    success: function(data, newValue, oldValue) {
                        if (data['success'] == 1) {
                            // Append gene symbols to search bar
                            gene_symbols = ''

                            //format gene symbols into search string
                            $.each(data['gene_symbols'], function(i, item){
                                gene_symbols += item['label'] + ' ';
                            });

                            $('#search_gene_symbol').val(gene_symbols);

                            // determine if searching for exact matches
                            if ( $('#exact_match_input').prop("checked") == false ) {
                                $('#exact_match_input').bootstrapToggle('on');
                            }

                        } else {
                            $('#selected_gene_cart').text(oldValue);
                            $('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
                              '<button type="button" class="close close-alert" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
                              '<p class="alert-message"><strong>Oops! </strong> ' + data["error"] + '</p></div>').show();
                        }
                        $("#search_gene_symbol").prop("disabled", false);
                        $("#selected_gene_cart_loading_c").hide();
                        d.resolve();
                    }
                });
            }
            return d.promise();
        },
        title: 'Click to change',
        tpl: '<select style="width:160px"></select>',
        showbuttons: false,
        source: function() {
            return gene_carts; //populated by load_gene_carts()
        },
        sourceCache: true
    });
}

// automatically reloads dataset grid and resubmits gene search
function update_datasetframes_generesults() {
    function resubmit_gene_search() {
        $('#gene_search_form').trigger('submit');
    }

    $.when( resubmit_gene_search() ).done(function(){
        // auto-select the first match.  first <a class="list-group-item"
        first_thing = $('#search_results a.list-group-item').first();
        select_search_result(first_thing);
    });
}

function set_exact_match(mode) {
    if (mode == 'on') {
        $('#exact_match_input').bootstrapToggle('on');
        $("#exact_match_icon img").attr("src", "img/arrow_target_selected.png");
        $("#exact_match_icon img").attr('data-original-title', "Exact match (currently on)").tooltip('show');
    } else if (mode == 'off') {
        $('#exact_match_input').bootstrapToggle('off');
        $("#exact_match_icon img").attr("src", "img/arrow_target_unselected.png");
        $("#exact_match_icon img").attr('data-original-title', "Exact match (currently off)").tooltip('show');
    }
}

function set_multigene_plots(mode) {
    if (mode == 'on') {
        $('#search_results_c').hide();
        $('#multigene_plots_input').bootstrapToggle('on');
    } else if (mode == 'off') {
        $('#search_results_c').show();
        $('#multigene_plots_input').bootstrapToggle('off');
    }
}
