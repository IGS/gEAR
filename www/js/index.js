let search_results = [];

// key dataset_id, value = Snap paths
const svgs = {};

let SCROLLBAR_DRAWN = false;
const GO_TERM_SCROLLBAR_DRAWN = false;
const AT_FIRST_MATCH_RECORD = false;
const AT_LAST_MATCH_RECORD  = false;
const PREVIOUS_SELECTED_RECORD_NUM = null;
let SCORING_METHOD = 'gene';
let SELECTED_GENE = null;

let share_id = null; //from permalink - dataset share ID
let layout_id = null; //from permalink - profile grid layout ID
const gene_cart_id = null; //from permalink - gene cart share ID
const permalinked_dataset_id = null; //holds dataset_id obtained from load_dataset_frames()
let multigene = false;  // Is this a multigene search?

const annotation_panel = new FunctionalAnnotationPanel();
const dataset_collection_panel = new DatasetCollectionPanel();

const profile_tree = new ProfileTree({treeDiv: '#profile_tree'});
const selected_profile_tree = new ProfileTree({treeDiv: '#selected_profile_tree'});

const gene_cart_tree = new GeneCartTree({treeDiv: '#selected_gene_cart_tree'});

const search_result_postselection_functions = [];

window.onload=() => {
    // check if the user is already logged in
    check_for_login();

    gene_cart_share_id = getUrlParameter('gene_cart_share_id');
    load_gene_carts(gene_cart_share_id);

    // Was a permalink found?
    share_id = getUrlParameter('share_id');
    scope = "permalink";

    if (share_id) {
        //hide site_into and display the permalink message
        $('#intro_content').hide();
        $('#viewport_intro').children().hide();
        $('#searching_indicator_c').hide();

        $('#leftbar_main').show();
        $('#permalink_intro_c').show();

        // validate the share_id. runs load_dataset_frames() on success
        validate_permalink(scope);
    } else {
        // layout_id is a share_id for the profile layout
        layout_id = getUrlParameter('layout_id');
        scope = "profile";
        get_index_info();

        if (document.URL.includes("index.html") ||
        window.location.pathname == '/' ) {
            load_layouts();
        }
    }

    // Was help_id found?
    const help_id = getUrlParameter('help_id');
    if (help_id) {
        validate_help_id(help_id);
    }

    const permalinked_gene_symbol = getUrlParameter('gene_symbol');
    const permalinked_gsem = getUrlParameter('gene_symbol_exact_match');
    const permalinked_multigene_plots = getUrlParameter('multigene_plots');
    multigene = (permalinked_multigene_plots && permalinked_multigene_plots === "1")

    if (permalinked_gene_symbol) {
        $("#search_gene_symbol_intro").val(permalinked_gene_symbol);

        if (permalinked_gsem && permalinked_gsem === "1") {
            set_exact_match('on');
        }

        if (multigene) {
            set_multigene_plots('on');
            set_exact_match('on');  // Multigene searches are exact.  This should speed up search_genes.py
        }

        sleep(1000).then(() => {
            $('#intro_search_icon').trigger('click');
        })
    } else if (share_id) {
        $('#permalink_intro_c').show();
    }

    // The search button starts out disabled, make sure it gets re-enabled.
    $("button#submit_search").prop( "disabled", false );

    $('#exact_match_icon').click(() => {
        // handle if it was already in the on position
        if ( $('#exact_match_input').prop("checked") ) {
            set_exact_match('off');
        } else {
            // else it was off, so turn it on
            set_exact_match('on');
        }
    });

    // If MG search icon on front page is clicked
    $('#multigene_search_icon').click(() => {
        // handle if it was already in the on position
        if ( $('#multigene_plots_input').prop("checked") ) {
            set_multigene_plots('off');
        } else {
            // else it was off, so turn it on
            set_multigene_plots('on');
        }
    });

    $('#intro_search_form').on('submit', (e) => {
        // TODO: It makes sense to remove/destroy those elements we aren't showing after a search
        e.preventDefault();
        $("#search_gene_symbol").val( $("#search_gene_symbol_intro").val());
        $('#intro_content').hide();

        $("#leftbar_main").show();
        $("#viewport_main").show();

        // fire the true search button
        $("#submit_search").trigger( "click" );
    });

    $('#intro_search_icon').click(() => {
        $('#intro_search_form').submit();
    });

    $('#dataset_search_form').on('submit', (e) => {
        e.preventDefault();
        window.location.replace(`./dataset_explorer.html?search_terms=${encodeURI($('#search_dataset_intro').val())}`);
    });

    $('#launcher_manual').click(() => {
        window.location.replace('./manual.html');
    });

    $('#launcher_expression_uploader').click(() => {
        window.location.replace('./upload_dataset.html');
    });

    $('#launcher_epigenetic_uploader').click(() => {
        window.location.replace('./upload_epigenetic_data.html');
    });

    $('.tool-launcher').click(function() {
        if ($(this).data('tool-name') == 'comparison') {
            window.location.replace('./compare_datasets.html');

        } else if ($(this).data('tool-name') == 'workbench') {
            window.location.replace('./analyze_dataset.html');
        } else if ($(this).data('tool-name') == 'mg_curator') {
            window.location.replace('./multigene_curator.html');
        }
    });

    // add post-page load listeners
    $( "#dataset_zoomed_zoom_out_control" ).click(() => {
        zoom_out_dataset();
    });

    // If a ProfileTree element is selected from the results page,
    // adjust state history and other results page things
    // These are adjustments that are normally made when the "search" button is hit
    $(document).on('change', '#selected_profile', () => {
        //TODO: get rid of redundancy with #search_gene_form.submit()

        // split on combination of space and comma (individually or both together.)
        const gene_symbol_array = $("#search_gene_symbol").val().split(/[\s,]+/);
        // Remove duplicates in gene search if they exist
        const uniq_gene_symbols = gene_symbol_array.filter((value, index, self) => self.indexOf(value) === index);
        const curated_searched_gene_symbols = uniq_gene_symbols.join(',');

        // determine if searching for exact matches (convert from bool to 1 or 0)
        $("#exact_match").val( Number($('#exact_match_input').prop("checked")) );
        $("#multigene_plots").val( Number($('#multigene_plots_input').prop("checked")) );

        // Update multigene toggle so correct grid widths are loaded.
        multigene = ($("#multigene_plots").val() && $("#multigene_plots").val() === "1")

        // Update search history
        add_state_history(curated_searched_gene_symbols);

        $("#too_many_genes_warning").hide();
        $('#search_result_count').text('');
        if (multigene) {
            // MG enabled
            $('#search_results_scrollbox').hide();
            $('#multigene_search_indicator').show();
            // Show warning if too many genes are entered
            if (uniq_gene_symbols.length > 10) {
                $("#too_many_genes_warning").text(`There are currently ${uniq_gene_symbols.length} genes to be searched and plotted. This can be potentially slow. Also be aware that with some plots, a high number of genes can make the plot congested or unreadable.`);
                $("#too_many_genes_warning").show();
            }
        } else {
            // MG disabled
            $('#search_results_scrollbox').show();
            $('#multigene_search_indicator').hide();
        }

        // Adjust num_genes badge (this does not use the result from search_genes.py so that may mismatch if "exact" is not chosen)
        // TODO: Actually run search_genes.py to get annotation information
        $('#search_result_count').text(uniq_gene_symbols.length);
    });

    // If a ProfileTree element is selected, this is changed and the new layout is set
    // NOTE: I don't think #search_param_profile needs to be a trigger
    $(document).on('change', '#search_param_profile, #selected_profile', function() {
        dataset_collection_panel.set_layout($(this).data('profile-id'), $(this).data('profile-label'), true, multigene);
        layout_id = $(this).data('profile-share-id');
    });

    $( document ).on("click", ".scope_choice", function() {
        SCORING_METHOD = $(this).data('choice');
        if (SELECTED_GENE !== null) {
            select_search_result($(SELECTED_GENE));
        }
    });

    // track the mouse movement so we can display scoring tooltips
    $( document ).on( "mousemove", (event) => {
        // Positioning for dataset_grid tips
        // Why is this pixel adjustment necessary?
        xpos = event.pageX - 240;
        ypos = event.pageY - 130 - 30;
        $("#tip").css("left", `${xpos}px` );
        $("#tip").css("top" , `${ypos}px` );
    });

    // Create observer to watch if user changes (ie. successful login does not refresh page)
    // See: https://developer.mozilla.org/en-US/docs/Web/API/MutationObserver

    // Select the node that will be observed for mutations
    const target_node = document.getElementById('loggedin_controls');
    // Create an observer instance linked to the callback function
    const observer = new MutationObserver(reload_trees);
    // For the "config" settings, do not monitor the subtree of nodes as that will trigger the callback multiple times.
    // Just seeing #loggedin_controls go from hidden (not logged in) to shown (logged in) is enough to trigger.
    observer.observe(target_node, { attributes: true });
};

function get_index_info() {
    $.ajax({
        url: './cgi/get_index_info.cgi',
        type: 'GET',
        dataType: 'json',
        success(data) {

            $('#stats_dataset_count').text(data.dataset_count);
            $('#stats_user_count').text(data.user_count);
        },
        error(jqXHR, _textStatus, errorThrown) {
            display_error_bar(`${jqXHR.status} ${errorThrown.name}`, 'Error getting index info.');
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
        success(data) {
            if ( data.success == 1 ) {
                // Add help_id to form
                $('#user_help_id').val(help_id);

                // Greet user by name (a subtle confirmation they know it's their account)
                if (data.user_name.length > 0) {
                    const user_first_name = ` ${data.user_name.split(' ')[0]}!`;
                    $('#forgot_password_user_name').text(user_first_name);
                }
            } else {
                // Invalid help_id, display invalid message
                $("#valid_forgot_pass_modal_body_c").hide();
                $("#save_user_new_pass").hide();
                $("#invalid_forgot_pass_modal_body_c").show();
            }
            $('#forgot_password_modal').modal('show');
        },
        error(jqXHR, _textStatus, errorThrown) {
            display_error_bar(`${jqXHR.status} ${errorThrown.name}`, 'Error validating help ID');
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

function validate_permalink(scope) {

    // Works for dataset or layout-based share IDs, which is differentiated by scope
    $.ajax({
        url : './cgi/validate_share_id.cgi',
        type: "POST",
        data : { 'share_id': share_id, 'scope': scope },
        dataType:"json",
        success: function(data, textStatus, jqXHR) {
            if ( data['success'] != 1 ) {
                $('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
                    '<button type="button" class="close close-alert" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
                    '<p class="alert-message"><strong>Oops! </strong> ' + data["error"] + '</p></div>').show();
            }

            if (scope == 'permalink') {
                dataset_collection_panel.load_frames({share_id});
            }
        },
        error: function (jqXHR, textStatus, errorThrown) {
            display_error_bar(jqXHR.status + ' ' + errorThrown.name, 'Error validating share ID');
        }
    });
}

function load_layouts() {
    const d = new $.Deferred();
    const session_id = Cookies.get('gear_session_id');
    const layout_share_id = getUrlParameter('layout_id');

    // Temporary hack for Heller lab
    if (layout_share_id == '8d38b600' || layout_share_id == 'afd2eb77') {
        $("#intro_selected_profile_warning").show();
    }

    //organize user and domain profiles in a tree format
    $.ajax({
        url: './cgi/get_user_layouts.cgi',
        type: 'post',
        async: false,
        data: { 'session_id': session_id, 'layout_share_id': layout_share_id },
        dataType: 'json',
        success(data, textStatus, jqXHR) {
            /*
              Priority of displayed profile:
                0.  Passed layout ID via layout_id URL parameter
                1.  Cookie value
                2.  User's DB-saved value (when they go to a machine, and there's no cookie)
                3.  Admin's active domain
             */

            const domain_profiles = [];
            const user_profiles = [];

            let active_layout_id = null;
            let active_layout_label = null;

            // Pass through once to sort domains from user profiles AND see if it matches a shared layout
            $.each(data['layouts'], (i, item) => {
                if ( item['is_domain'] == 1 ) {
                    domain_profiles.push({value: item['id'], text: item['label'], share_id: item['share_id'] });
                } else {
                    user_profiles.push({value: item['id'], text: item['label'], share_id: item['share_id']  });
                }

                if (item['share_id'] == layout_share_id) {
                    active_layout_id = item['id'];
                    active_layout_label = item['label'];
                    layout_id = item['share_id'];
                }
            });

            // Generate the tree structure for the layouts
            profile_tree.domainProfiles = domain_profiles;
            profile_tree.userProfiles = user_profiles;
            profile_tree.generateTree();
            selected_profile_tree.domainProfiles = domain_profiles;
            selected_profile_tree.userProfiles = user_profiles;
            selected_profile_tree.generateTree();

            // pass through again and look for one set by a cookie
            if (active_layout_id == null) {
                $.each(data.layouts, (_i, item) => {
                    if (item.label == CURRENT_USER.profile) {
                        active_layout_id = item.id;
                        active_layout_label = item.label;
                        layout_id = item.share_id;
                        return false;
                    }
                });
            }

            // pass through again and look for one set as current by the user
            if (active_layout_id == null) {
                $.each(data.layouts, (_i, item) => {
                    if ( item['is_domain'] == 0 && item.is_current == 1 ) {
                        active_layout_id = item.id;
                        active_layout_label = item.label;
                        layout_id = item.share_id;
                        return false;
                    }
                });
            }

            // pass through again if no active layout was found for user and choose the admin's
            if (active_layout_id == null) {
                $.each(data.layouts, (_i, item) => {
                    if ( item['is_domain'] == 1 && item.is_current == 1 ) {
                        active_layout_id = item.id;
                        active_layout_label = item.label;
                        layout_id = item.share_id;
                        return false;
                    }
                });
            }

            dataset_collection_panel.set_layout(active_layout_id, active_layout_label, false, multigene);

            d.resolve();
        },
        error(jqXHR, _textStatus, errorThrown) {
            profile_tree.generateTree();
            selected_profile_tree.generateTree();
            display_error_bar(`${jqXHR.status} ${errorThrown.name}`, 'Error loading layouts.');
            d.fail();
        }
    });

    d.promise();
}

function load_gene_carts(cart_share_id) {
    const d = new $.Deferred();
    const session_id = Cookies.get('gear_session_id');

    if (!session_id) {
        //User is not logged in. Hide gene carts container
        $("#selected_gene_cart_c").hide();
        gene_cart_tree.generateTree();
        d.resolve();
    } else {
        $("#selected_gene_cart_c").show(); //Show if hidden
        $.ajax({
        url: './cgi/get_user_gene_carts.cgi',
        type: 'post',
            data: { 'session_id': session_id, 'share_id': cart_share_id },
        dataType: 'json',
        success(data, textStatus, jqXHR) { //source https://stackoverflow.com/a/20915207/2900840
            const carts = {};
            let permalink_cart_id = null
            let permalink_cart_label = null
            const cart_types = ['domain', 'user', 'group', 'shared', 'public'];
            let carts_found = false;

            for (const ctype of cart_types) {
                carts[ctype] = [];

                if (data[`${ctype}_carts`].length > 0) {
                    carts_found = true;

                    //User has some profiles
                    $.each(data[`${ctype}_carts`], (_i, item) => {
                        // If cart permalink was passed in, retrieve gene_cart_id for future use.
                        if (cart_share_id && item.share_id == cart_share_id) {
                            permalink_cart_id = item.id;
                            permalink_cart_label = item.label;
                        }

                        carts[ctype].push({value: item.id, text: item.label });
                    });
                }
            }

            gene_cart_tree.domainGeneCarts = carts.domain;
            gene_cart_tree.userGeneCarts = carts.user;
            gene_cart_tree.groupGeneCarts = carts.group;
            gene_cart_tree.sharedGeneCarts = carts.shared;
            gene_cart_tree.publicGeneCarts = carts.public;
            gene_cart_tree.generateTree();

            if (! carts_found ) {
                $("#selected_gene_cart_c").hide();
            }

            // If gene_cart_permalink was provided:
            // 1) Set the value in the gene cart tree and gene search bar
            // 2) Trigger change event to populate the gene search bar with the genes
            // 3) Click the "search" button on the front page (happens via another process that actually calls load_frames)
            // 4) show sidebar stuff in the display panel
            if (permalink_cart_id) {
                $("#selected_gene_cart").text(permalink_cart_label);
                $("#selected_gene_cart").val(permalink_cart_id);
                $("#selected_gene_cart").trigger('change');

                //hide site_into and display the permalink message
                $('#intro_content').hide();
                $('#viewport_intro').children().hide();
                $('#searching_indicator_c').hide();

                $('#leftbar_main').show();
                $('#permalink_intro_c').show();
            }

            d.resolve();
        },
        error(jqXHR, textStatus, errorThrown) {
            gene_cart_tree.generateTree();
            display_error_bar(`${jqXHR.status} ${errorThrown.name}`, "Gene carts not sucessfully loaded.");
            d.fail();
        }
        });
    }
    d.promise();
}

// If user changes, update genecart/profile trees
async function reload_trees(){
    // Update dataset and genecart trees in parallel
    // Works if they were not populated or previously populated
    await Promise.all([load_layouts(), load_gene_carts(gene_cart_share_id)]);
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

    for (const key in data) {
        if (data.hasOwnProperty(key)) {
            sorted_gene_syms.push(key);
        }
    }

    sorted_gene_syms.sort();
    sorted_gene_syms_len = sorted_gene_syms.length

    const items = [];

    for (i = 0; i < sorted_gene_syms_len; i++) {
        gene_symbol = sorted_gene_syms[i];

        // Build search result html
        let gene_result_html = `<a class="list-group-item" data-gene_symbol="${gene_symbol}" href="#">${gene_symbol}`;

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
            $('#search_result_count').text(`max:${items.length}`);
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

$("#gene_search_form").submit((event) => {
    $("#viewport_intro").hide();
    $("#viewport_main").show();

    // clear any open tooltips
    $('[data-toggle="tooltip"], .tooltip').tooltip("hide");

    // determine if searching for exact matches (convert from bool to 1 or 0)
    $("#exact_match").val( Number($('#exact_match_input').prop("checked")) );
    $("#multigene_plots").val( Number($('#multigene_plots_input').prop("checked")) );

    $('#recent_updates_c').hide();
    $('#searching_indicator_c').show();

    const formData = $("#gene_search_form").serializeArray();

    // split on combination of space and comma (individually or both together.)
    const gene_symbol_array = $("#search_gene_symbol").val().split(/[\s,]+/);
    // Remove duplicates in gene search if they exist
    const uniq_gene_symbols = gene_symbol_array.filter((value, index, self) => self.indexOf(value) === index);
    const curated_searched_gene_symbols = uniq_gene_symbols.join(',');

    // Update multigene toggle so correct grid widths are loaded.
    multigene = ($("#multigene_plots").val() && $("#multigene_plots").val() === "1")

    add_state_history(curated_searched_gene_symbols);

    $("#too_many_genes_warning").hide();
    $('#search_result_count').text('');
    if (multigene) {
        // MG enabled
        $('#search_results_scrollbox').hide();
        $('#multigene_search_indicator').show();
        // Show warning if too many genes are entered
        if (uniq_gene_symbols.length > 10) {
            $("#too_many_genes_warning").text(`There are currently ${uniq_gene_symbols.length} genes to be searched and plotted. This can be potentially slow. Also be aware that with some plots, a high number of genes can make the plot congested or unreadable.`);
            $("#too_many_genes_warning").show();
        }
    } else {
        // MG disabled
        $('#search_results_scrollbox').show();
        $('#multigene_search_indicator').hide();
    }

    $('#search_results').empty();
    // show search results
    $('#search_results_c').removeClass('search_result_c_DISABLED');

    dataset_collection_panel.load_frames({share_id, multigene});

    $.ajax({
        url : './cgi/search_genes.py',
        type: "POST",
        data : formData,
        dataType:"json",
        success(data, textStatus, jqXHR) {
        	// reset search_results
        	search_results = data;
            populate_search_result_list(data);
            $('#searching_indicator_c').hide();
            $('#intro_content').hide('fade', {}, 400, () => {
                if ($('#multigene_plots').val() == 1){
                    dataset_collection_panel.update_by_all_results(uniq_gene_symbols);
                } else {
                    // auto-select the first match.  first <a class="list-group-item"
                    const first_thing = $('#search_results a.list-group-item').first();
                    select_search_result(first_thing);
                }
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
        error(jqXHR, textStatus, errorThrown) {
            $('#searching_indicator_c').hide();

            // Error occurred
            if ( $('#search_gene_symbol').val().length < 1 ) {
            // No gene symbol entered
                    $('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
                        '<button type="button" class="close close-alert" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
                        '<p class="alert-message"><strong>Oops! </strong>No gene symbol was entered. Enter a gene symbol and try again.</p></div>').show();


            } else if ( dataset_collection_panel.datasets.length == 0) {
            // No datasets in current layout profile
                    $('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
                        '<button type="button" class="close close-alert" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
                        '<p class="alert-message"><strong>Oops! </strong>No datasets were found in the current layout profile.</p><p>To add datasets to a profile or choose a different profile, go to the <a href="./dataset_manager.html" class="alert-link">Dataset Manager</a>.</p></div>').show();

            } else {
                // Some other error occurred
              	display_error_bar(`${jqXHR.status} ${errorThrown.name}`, "Could not successfully search for genes in database.");
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

// When a gene cart is selected, populate the gene search bar with its members
$('#selected_gene_cart').change( function() {
    let geneCartId = $(this).val();
    const params = { session_id: session_id, gene_cart_id: geneCartId };
    const d = new $.Deferred(); // Causes editable to wait until results are returned

    if (typeof session_id !== 'undefined') {
        // Get the gene cart members and populate the gene symbol search bar
        $.ajax({
        url: './cgi/get_gene_cart_members.cgi',
        type: 'post',
        data: params,
        success: function (data, newValue, oldValue) {
            if (data.success === 1) {
                const gene_symbols_array = []
                // format gene symbols into search string
                $.each(data.gene_symbols, function (i, item) {
                    gene_symbols_array.push(item.label);
                });
                //deduplicate gene cart
                const dedup_gene_symbols_array = [...new Set(gene_symbols_array)]

                gene_symbols = dedup_gene_symbols_array.join(' ')
                $('#search_gene_symbol_intro').val(gene_symbols);
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
    } else {
        d.resolve();
    }
    return d.promise();
});

// Set the state history based on current conditions
function add_state_history(gene_symbols) {
    const state_info = {
        'gene_symbol_exact_match': $("#exact_match").val(),
        'multigene_plots': $("#multigene_plots").val()
    };

    let state_url = "/index.html?"
                + `&gene_symbol_exact_match=${$("#exact_match").val()}`
                + `&multigene_plots=${$("#multigene_plots").val()}`;

    if (layout_id) {
        state_info.layout_id = layout_id;
        state_url += `&layout_id=${layout_id}`;
    }

    if (gene_symbols) {
        state_info.gene_symbol = gene_symbols;
        state_url += `&gene_symbol=${gene_symbols}`;
    }

    if (getUrlParameter('gene_cart_share_id')) {
        state_info.gene_cart_share_id = getUrlParameter('gene_cart_share_id');
        state_url += `&gene_cart_share_id=${getUrlParameter('gene_cart_share_id')}`;
    }

    // SAdkins - Should we have a separate history state for dataset share IDs?
    history.pushState(
        // State Info
        state_info,
        // State title
        "Gene search",
        // URL
        state_url
    )
}

// automatically reloads dataset grid and resubmits gene search
function update_datasetframes_generesults() {
    function resubmit_gene_search() {
        return;
        // SAdkins - 1/6/22 - commenting out since it loads double the frames
        //$('#gene_search_form').trigger('submit');
    }

    $.when( resubmit_gene_search() ).done(() => {
        if ($('#multigene_plots').val() == 1){
            // split on combination of space and comma (individually or both together.)
            const gene_symbol_array = $("#search_gene_symbol").val().split(/[\s,]+/);
            // Remove duplicates in gene search if they exist
            const uniq_gene_symbols = gene_symbol_array.filter((value, index, self) => self.indexOf(value) === index);
            dataset_collection_panel.update_by_all_results(uniq_gene_symbols);
        } else {
            // auto-select the first match.  first <a class="list-group-item"
            const first_thing = $('#search_results a.list-group-item').first();
            select_search_result(first_thing);
        }
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
        $('#multigene_plots_input').bootstrapToggle('on');  // Toggles gene results display things upon change
        $("#multigene_search_icon i").attr('data-original-title', "Multigene displays enabled. Click to search for single-gene displays ").tooltip('show');
        $("#multigene_search_icon i").addClass("fa-gears");
        $("#multigene_search_icon i").removeClass("fa-gear");
    } else if (mode == 'off') {
        $('#multigene_plots_input').bootstrapToggle('off');
        $("#multigene_search_icon i").attr('data-original-title', "Single-gene displays enabled. Click to search for multigene displays").tooltip('show');
        $("#multigene_search_icon i").addClass("fa-gear");
        $("#multigene_search_icon i").removeClass("fa-gears");
    }
}

