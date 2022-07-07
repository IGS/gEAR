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

let dataset_id = null; //from permalink - dataset share ID
let layout_id = null; //from permalink - profile grid layout ID
let gene_cart_id = null; //from permalink - gene cart share ID
let multigene = false;  // Is this a multigene search?
let exact_match = true; // Set on by default
let projection = false;

const annotation_panel = new FunctionalAnnotationPanel();
const dataset_collection_panel = new DatasetCollectionPanel();
const { controller } = dataset_collection_panel;

/*
Tree properties for constructor:
treeDiv - Element to generate the tree structure on
storedValElt - Element to store text, vals, and data properties on (if not in a treeDiv descendant "dropdown-toggle" element)
*/

const profile_tree = new ProfileTree({treeDiv: '#profile_tree'});
const selected_profile_tree = new ProfileTree({treeDiv: '#selected_profile_tree'});

const gene_cart_tree = new GeneCartTree({treeDiv: '#gene_cart_tree', storedValElt: '#search_param_gene_cart'});
const selected_gene_cart_tree = new GeneCartTree({treeDiv: '#selected_gene_cart_tree', storedValElt: '#selected_gene_cart'});

//const projection_source_tree = new ProjectionSourceTree({treeDiv: '#projection_source_tree'});

const search_result_postselection_functions = [];

$(document).on("handle_page_loading", () => {

    // Ensure "exact match" and "multigene" tooltips work upon page load
    $('#intro_search_div [data-toggle="tooltip"]').tooltip();

    // Was a permalink found?
    dataset_id = getUrlParameter('share_id');
    scope = "permalink";

    if (dataset_id) {
        //hide site_into and display the permalink message
        $('#intro_content').hide();
        $('#viewport_intro').children().hide();
        $('#searching_indicator_c').hide();
        $('#selected_profile_container').hide();

        $('#leftbar_main').show();
        $('#permalink_intro_c').show();

        show_gene();

        // validate the dataset_id. runs load_frames() on success
        validate_permalink(scope);
    } else {
        // layout_id is a share_id for the profile layout
        layout_id = getUrlParameter('layout_id');
        scope = "profile";
        get_index_info();
    }

    // Was help_id found?
    const help_id = getUrlParameter('help_id');
    if (help_id) {
        validate_help_id(help_id);
    }

    gene_cart_id = getUrlParameter('gene_cart_share_id');
    if (gene_cart_id) {
        console.info(`Gene cart share ID found: ${gene_cart_id}`);
        $('#intro_search_icon').trigger('click');
    }

    const permalinked_gsem = getUrlParameter('gene_symbol_exact_match');
    if (permalinked_gsem !== (null || undefined)
        && permalinked_gsem === "0") {
        exact_match = false;
    }
    set_exact_match(exact_match, false);

    const permalinked_multigene_plots = getUrlParameter('multigene_plots');
    multigene = (permalinked_multigene_plots !== (null || undefined)
        && permalinked_multigene_plots === "1");
    set_multigene_plots(multigene, false);

    if (multigene) {
        exact_match = multigene;
        set_exact_match(exact_match, false);  // Multigene searches are exact.  This should speed up search_genes.py
    }

    // If gene symbols were provided (via either URL param method), click search button.
    const permalinked_gene_symbol = getUrlParameter('gene_symbol');
    if (permalinked_gene_symbol) {
        $("#search_gene_symbol_intro").val(permalinked_gene_symbol);

        console.info(`Permalinked gene symbols found: ${permalinked_gene_symbol}`);
        $('#intro_search_icon').trigger('click');
    } else if (dataset_id) {
        $('#permalink_intro_c').show();
    }

    // Repopulate projection information... projection_source URL loaded earlier
    const permalinked_is_pca = getUrlParameter('is_pca')
    const permalinked_projection_patterns = getUrlParameter('projection_patterns');

    // Only apply if both are present
    const permalinked_projection_id = getUrlParameter('projection_source');
    if (permalinked_projection_id) {
        let selected_projections_string;
        if (permalinked_projection_patterns) {
            selected_projections_string = permalinked_projection_patterns;
        } else {
            // jQuery .map() returns a jQuery object, not an array, so use .get() to convert to array
            const selected_projections =  $('.js-projection-pattern-elts-check').map(function() {
                return $(this).data('label');
            }).get();
            selected_projections_string = selected_projections.join(',');
        }

        // Patterns applied after HTML renders
        $("#search_gene_symbol_intro").val(selected_projections_string);
        // Check boxes for the elements that were found in the URL
        selected_projections_string.split(',').forEach((pattern) => {
            $(`.js-projection-pattern-elts-check[data-label="${pattern}"]`).prop('checked', true);
        });

        $("#is_pca").prop('checked', permalinked_is_pca === "1");

        // Correct tab is active
        $("#projection_tab").click();
        projection = true;
        console.info(`Projection ID found: ${permalinked_projection_id}`);
        $('#intro_search_icon').trigger('click');
    }

    // TODO: I think everything beyond this point could actually be in a window.onload event or set global.

    // The search button starts out disabled, make sure it gets re-enabled.
    $("button#submit_search").prop( "disabled", false );

    // If exact match icon is clicked, toggle parameters
    $('.js-exact-match').click(() => {
        exact_match = !exact_match;
        set_exact_match(exact_match);
    });

    // If MG search icon is clicked, toggle parameters
    $('.js-multigene').click(() => {
        multigene = !multigene;
        set_multigene_plots(multigene);
    });

    // If multi-pattern set, toggle multigene
    $('input[name="projection_display_mode"]').change(() => {
        multigene = $('#multi_pattern').is(':checked');
        set_multigene_plots(multigene, false);
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

    // If a ProfileTree element is selected, this is changed and the new layout is set
    // NOTE: I don't think #search_param_profile needs to be a trigger
    $(document).on('change', '#search_param_profile, #selected_profile', function() {
        dataset_collection_panel.set_layout($(this).data('profile-id'), $(this).data('profile-label'), true, multigene);
        layout_id = $(this).data('profile-share-id');

        // If a ProfileTree element is selected from the results page,
        // adjust state history and other results page things
        // These are adjustments that are normally made when the "search" button is hit
        if (this.id === "selected_profile") {

            // split on combination of space and comma (individually or both together.)
            const gene_symbol_array = $("#search_gene_symbol").val().split(/[\s,]+/);
            // Remove duplicates in gene search if they exist
            const uniq_gene_symbols = gene_symbol_array.filter((value, index, self) => self.indexOf(value) === index);
            const curated_searched_gene_symbols = uniq_gene_symbols.join(',');

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
        }
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

    // But we need to wait for navigation_bar to load first (in common.js) so do some polling
    // See: https://stackoverflow.com/q/38881301

    // Select the node that will be observed for mutations
    const target_node = document.getElementById('loggedin_controls');
    const safer_node = document.getElementById("navigation_bar");   // Empty div until loaded
    // Create an observer instance linked to the callback function
    const observer = new MutationObserver(function(_mutationList, _observer) {
        if (target_node) {
            load_all_trees();
            this.disconnect();  // Don't need to reload once the trees are updated
        }
    });
    // For the "config" settings, do not monitor the subtree of nodes as that will trigger the callback multiple times.
    // Just seeing #loggedin_controls go from hidden (not logged in) to shown (logged in) is enough to trigger.
    observer.observe(target_node || safer_node , { attributes: true });
});

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

// TODO: If search is clicked before trees are loaded, gene is searched twice as submission happens twice.
// Disable search button until trees are loaded.
$(document).on("build_jstrees", async () => await load_all_trees());

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
    const pass_1 = $('input#user_new_pass_1').val();

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

    const help_id = $('input#user_help_id').val();
    const new_password = $('input#user_new_pass_2').val();
    $.ajax({
        url: './cgi/save_user_account_changes.cgi',
        type: 'POST',
        data: { 'help_id': help_id, 'new_password': new_password, 'scope': 'password'},
        dataType: 'json',
        success: function(data, textStatus, jqXHR) {
            if (data.success == 1) {
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
                    '<p class="alert-message"><strong>Oops! </strong> ' + data.error + '</p></div>').show();
            }
        },
        error: function (jqXHR, textStatus, errorThrown) {
            display_error_bar(`${jqXHR.status} ${errorThrown.name}`, 'Error saving new user password');
        }
    });//end ajax
});

function validate_permalink(scope) {
    // Works for dataset or layout-based share IDs, which is differentiated by scope
    $.ajax({
        url : './cgi/validate_share_id.cgi',
        type: "POST",
        data : { 'share_id': dataset_id, 'scope': scope },
        dataType:"json",
        success: function(data, textStatus, jqXHR) {
            if ( data.success != 1 ) {
                $('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
                    '<button type="button" class="close close-alert" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
                    '<p class="alert-message"><strong>Oops! </strong> ' + data.error + '</p></div>').show();
            }

            if (scope == 'permalink') {
                dataset_collection_panel.load_frames({dataset_id});
            }
        },
        error: function (jqXHR, textStatus, errorThrown) {
            display_error_bar(`${jqXHR.status} ${errorThrown.name}`, 'Error validating share ID');
        }
    });
}

async function load_layouts() {
    const layout_share_id = getUrlParameter('layout_id');
    let active_layout_id = null;
    let active_layout_label = null;

    // Temporary hack for Heller lab
    if (layout_share_id == '8d38b600' || layout_share_id == 'afd2eb77') {
        $("#intro_selected_profile_warning").show();
    }

    //organize user and domain profiles in a tree format
    await $.ajax({
        url: './cgi/get_user_layouts.cgi',
        type: 'post',
        data: { 'session_id': session_id, 'layout_share_id': layout_share_id },
        dataType: 'json'
    }).done((data) => {
        /*
            Priority of displayed profile:
            0.  Passed layout ID via layout_id URL parameter
            1.  Cookie value
            2.  User's DB-saved value (when they go to a machine, and there's no cookie)
            3.  Admin's active domain
            */
        const layouts = {};
        const layout_types = ['domain', 'user', 'group', 'shared']

        for (const ltype of layout_types) {
            layouts[ltype] = [];

            $.each(data[`${ltype}_layouts`], (_i, item) => {
                layouts[ltype].push({value: item['id'],
                                        text: item['label'],
                                        share_id: item['share_id'],
                                        folder_id: item['folder_id'],
                                        folder_parent_id: item['folder_parent_id'],
                                        folder_label: item['folder_label']
                                    });

                if (item['share_id'] == layout_share_id) {
                    active_layout_id = item.id;
                    active_layout_label = item.label;
                    layout_id = item.share_id;
                }
            });
        }

        // Generate the tree structure for the layouts
        profile_tree.domainProfiles = layouts.domain;
        profile_tree.userProfiles = layouts.user;
        profile_tree.groupProfiles = layouts.group;
        profile_tree.sharedProfiles = layouts.shared;
        selected_profile_tree.domainProfiles = layouts.domain;
        selected_profile_tree.userProfiles = layouts.user;
        selected_profile_tree.groupProfiles = layouts.group;
        selected_profile_tree.sharedProfiles = layouts.shared;

        // pass through again and look for one set by a cookie
        if (active_layout_id == null) {
            for (const ltype of layout_types) {
                $.each(data[`${ltype}_layouts`], (_i, item) => {
                    if (item.label == CURRENT_USER.profile) {
                        active_layout_id = item.id;
                        active_layout_label = item.label;
                        layout_id = item.share_id;
                        return false;
                    }
                });
            }
        }

        // pass through again and look for one set as current by the user
        if (active_layout_id == null) {
            for (const ltype of ['user', 'group']) {
                $.each(data[`${ltype}_layouts`], (_i, item) => {
                    if ( item['is_domain'] == 0 && item.is_current == 1 ) {
                        active_layout_id = item.id;
                        active_layout_label = item.label;
                        layout_id = item.share_id;
                        return false;
                    }
                });
            }
        }

        // pass through again if no active layout was found for user and choose the admin's
        if (active_layout_id == null) {
            $.each(data['domain_layouts'], (_i, item) => {
                if ( item['is_domain'] == 1 && item.is_current == 1 ) {
                    active_layout_id = item.id;
                    active_layout_label = item.label;
                    layout_id = item.share_id;
                    return false;
                }
            });
        }

        dataset_collection_panel.set_layout(active_layout_id, active_layout_label, false, multigene);

    }).fail((jqXHR, textStatus, errorThrown) => {
        display_error_bar(`${jqXHR.status} ${errorThrown.name}`, 'Error loading layouts.');
    });

    profile_tree.generateTree();
    selected_profile_tree.generateTree();

}

async function load_gene_carts() {
    let carts_found = false;
    let permalink_cart_id = null
    let permalink_cart_label = null
    $("#selected_gene_cart_c").prop("disabled", false);
    const cart_share_id = getUrlParameter('gene_cart_share_id');

    if (!session_id) {
        //User is not logged in. Hide gene carts container
        $("#selected_gene_cart_c").prop("disabled", true);
        gene_cart_tree.generateTree();
        selected_gene_cart_tree.generateTree();
        return;
    }
    await $.ajax({
        url: './cgi/get_user_gene_carts.cgi',
        type: 'post',
        data: { 'session_id': session_id, 'share_id': cart_share_id },
        dataType: 'json'
    }).done((data, textStatus, jqXHR) => {
        const carts = {};
        const cart_types = ['domain', 'user', 'group', 'shared', 'public'];
        for (const ctype of cart_types) {
            carts[ctype] = [];

            if (data[`${ctype}_carts`].length > 0) {
                carts_found = true;

                $.each(data[`${ctype}_carts`], (_i, item) => {
                    // If cart permalink was passed in, retrieve gene_cart_id for future use.
                    if (cart_share_id && item.share_id == cart_share_id) {
                        permalink_cart_id = item.id;
                        permalink_cart_label = item.label;
                    }

                    carts[ctype].push({value: item.id,
                                        text: item.label,
                                        folder_id: item.folder_id,
                                        folder_label: item.folder_label,
                                        folder_parent_id: item.folder_parent_id
                                        });
                });
            }
        }

        gene_cart_tree.domainGeneCarts = carts.domain;
        gene_cart_tree.userGeneCarts = carts.user;
        gene_cart_tree.groupGeneCarts = carts.group;
        gene_cart_tree.sharedGeneCarts = carts.shared;
        gene_cart_tree.publicGeneCarts = carts.public;
        selected_gene_cart_tree.domainGeneCarts = carts.domain;
        selected_gene_cart_tree.userGeneCarts = carts.user;
        selected_gene_cart_tree.groupGeneCarts = carts.group;
        selected_gene_cart_tree.sharedGeneCarts = carts.shared;
        selected_gene_cart_tree.publicGeneCarts = carts.public;

    }).fail((jqXHR, textStatus, errorThrown) => {
        display_error_bar(`${jqXHR.status} ${errorThrown.name}`, "Gene carts not sucessfully loaded.");
    });
    if (! carts_found ) {
        $("#selected_gene_cart_c").prop("disabled", true);
    }
    gene_cart_tree.generateTree();
    selected_gene_cart_tree.generateTree();

    // If gene_cart_permalink was provided:
    // 1) Set the value in the gene cart tree and gene search bar
    // 2) Trigger change event to populate the gene search bar with the genes
    // 3) (Outside of function) Search button is clicked
    // 4) (Outside of function) Show sidebar stuff in the display panel
    if (permalink_cart_id) {
        // This will also change selected_gene_cart via the "change" event trigger
        $("#search_param_gene_cart").text(permalink_cart_label);
        $("#search_param_gene_cart").val(permalink_cart_id);
        $("#search_param_gene_cart").trigger('change');
    }
}

async function load_weighted_gene_carts(cart_share_id) {
    let permalink_cart_id = null
    let permalink_cart_label = null

    if (!session_id) {
        console.info("User is not logged in. Weighted gene carts not loaded.");
        return [permalink_cart_id, permalink_cart_label];
    }

    await $.ajax({
        url: './cgi/get_user_gene_carts.cgi',
        type: 'post',
        data: { 'session_id': session_id, "cart_type": "weighted-list" },
        dataType: 'json'
    }).done((data, textStatus, jqXHR) => {
        const carts = {};

        const cart_types = ['domain', 'user', 'group', 'shared', 'public'];

        for (const ctype of cart_types) {
            carts[ctype] = [];
            if (data[`${ctype}_carts`].length > 0) {
                $.each(data[`${ctype}_carts`], (_i, item) => {
                    const share_id = `cart.${item.share_id}`; // normalizing name for easy filepath retrieval

                    // If cart permalink was passed in, retrieve gene_cart_id for future use.
                    if (cart_share_id && share_id == cart_share_id) {
                        permalink_cart_id = share_id;
                        permalink_cart_label = item.label;
                    }


                    carts[ctype].push({value: share_id,    // Use share ID as it is used in the cart file's basename
                                        text: item.label,
                                        folder_id: item.folder_id,
                                        folder_label: item.folder_label,
                                        folder_parent_id: item.folder_parent_id
                                        });
                });
            }
        }

        // Tree is generated in `load_pattern_tree`
        projection_source_tree.domainGeneCarts = carts.domain;
        projection_source_tree.userGeneCarts = carts.user;
        projection_source_tree.groupGeneCarts = carts.group;
        projection_source_tree.sharedGeneCarts = carts.shared;
        projection_source_tree.publicGeneCarts = carts.public;

    })
    .fail((jqXHR, textStatus, errorThrown) => {
        display_error_bar(`${jqXHR.status} ${errorThrown.name}`, "Weighted gene carts not sucessfully retrieved.");
    });
    return [permalink_cart_id, permalink_cart_label];
}

async function populate_pattern_selection(projection_source) {
    let permalink_projection_id = null
    let permalink_projection_label = null

    await $.ajax({
        type: "POST",
        url: "./cgi/get_projection_pattern_list.cgi",
        //data: {session_id},
        dataType: "json",
    }).done((data) => {
        const patterns_list = [];

        $.each(data, (_i, item) => {
            if (projection_source && item.id == projection_source) {
                permalink_projection_id = item.id;
                permalink_projection_label = item.title;
            }
            patterns_list.push({value: item.id,
                 text: item.title
            });
        });

        // Tree is generated in `load_pattern_tree`
        projection_source_tree.projectionPatterns = patterns_list;
    }).fail((jqXHR, textStatus, errorThrown) => {
        display_error_bar(`${jqXHR.status} ${errorThrown.name}`,`Failed to populate patterns list`);
    });
    return [permalink_projection_id, permalink_projection_label];
}

async function load_pattern_tree() {
    const projection_id = getUrlParameter('projection_source');
    const values = await Promise.allSettled([load_weighted_gene_carts(projection_id), populate_pattern_selection(projection_id)])
        .catch((err) => {
            console.error(err);
        });

    // If projection info was in URL, one of the above should have the JSTree element returned
    projection_source_tree.generateTree();

    // If a permalink was provided, set the value in the tree and search bar
    for (const val of values) {
        if (val.value[0]) {
            $("#projection_source").text(val.value[1]);
            $("#projection_source").val(val.value[0]);
            // At this point, the tree is generated but loading data attributes to the storedValElt does not occur until a node is selected.
            // So we need to manually set the data attribute for the first-pass.
            const tree_leaf = projection_source_tree.treeData.find(e => e.id === $("#projection_source").val());
            $("#projection_source").data("scope", tree_leaf.scope);
            $("#projection_source").trigger('change');
            return;
        }
    }
}

// If user changes, update genecart/profile trees
async function load_all_trees(){
    // Update dataset and genecart trees in parallel
    // Works if they were not populated or previously populated
    await Promise.allSettled([load_layouts(), load_gene_carts()])//, load_pattern_tree()])
        .catch((err) => {
            console.error(err)
        });
    console.info("Trees loaded");

    // NOTE: This will trigger again if the MutationObserver catches a login, but that may be acceptable.
    $(document).trigger("handle_page_loading");

}

// Hide option menu when scope is changed.
$(document).on('click', '.scope_choice', function(){
    $('#toggle_options').popover('hide');
});

// Handle direct documentation links
$(document).on('click', '#doc-link-choices li', function(){
    window.location.replace(`./manual.html?doc=${$(this).data('doc-link')}`);
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
    const callTime = new Date().getTime();
    if (callTime - lastCall <= 500) {
      return false;
    }
    lastCall = callTime;

    SELECTED_GENE = $(elm);
    gene_sym = $(elm).data("gene_symbol");

    // remove coloring from other result links
    $('.list-group-item-active').removeClass('list-group-item-active');
    $(elm).addClass('list-group-item-active');

    if (!projection) {
        annotation_panel.annotation = search_results[gene_sym];
        annotation_panel.autoselect_organism();
    }

    // hide the intro, show the search result box
    if( $('#site_intro_c').is(':visible') ) {
        $('#site_intro_c').hide({easing: 'fade', duration: 400});
        $('#recent_updates_c').hide({easing: 'fade', duration: 400, complete: show_search_result_info_box});
    }

    dataset_collection_panel.update_by_search_result(search_results[gene_sym]);

    // call any plugin functions
    search_result_postselection_functions.forEach((f) => {f()})
}

function isNumeric(n) {
    return !isNaN(parseFloat(n)) && isFinite(n);
}

function show_search_result_info_box() {
    if( dataset_id ) {
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

    dataset_collection_panel.reset_abort_controller();

    let draw=true;
    if (projection) {
        draw=false;

        dataset_collection_panel.datasets.forEach((dataset) => {
            dataset.draw({gene_symbol: $(this).data("gene_symbol")});
        });
    }
    select_search_result(this, draw_display=draw);
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
    show_gene();

    // re-initialize any open tooltips
    $('[data-toggle="tooltip"], .tooltip').tooltip();
    // Ensure correct tooltip text is displayed
    set_exact_match(exact_match, false);
    set_multigene_plots(multigene, false);

    $('#recent_updates_c').hide();
    $('#searching_indicator_c').show();

    const formData = $("#gene_search_form").serializeArray();

    // split on combination of space and comma (individually or both together.)
    const gene_symbol_array = $("#search_gene_symbol").val().split(/[\s,]+/);
    // Remove duplicates in gene search if they exist
    const uniq_gene_symbols = gene_symbol_array.filter((value, index, self) => self.indexOf(value) === index);
    const curated_searched_gene_symbols = uniq_gene_symbols.join(',');

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

    dataset_collection_panel.load_frames({dataset_id, multigene});

    // Add Exact Match param
    formData.push({"name": "exact_match", "value" :Number(exact_match)});

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
                if (multigene){
                    dataset_collection_panel.update_by_all_results(uniq_gene_symbols);
                } else {
                    // auto-select the first match.  first <a class="list-group-item"
                    const first_thing = $('#search_results a.list-group-item').first();
                    select_search_result(first_thing);
                }
            });

            set_scrollbar_props();
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

$("#projection_search_form").submit((event) => {
    $("#viewport_intro").hide();
    $("#viewport_main").show();
    show_projection();

    // re-initialize any open tooltips
    $('[data-toggle="tooltip"], .tooltip').tooltip();

    // If front page MG was enabled, ensure correct radio button is selected
    if (multigene) {
        $('#multi_pattern').click();
    }

    $('#recent_updates_c').hide();
    $('#searching_indicator_c').show();

    // Get selected projections and add as state
    const selected_projections = [];
    $('.js-projection-pattern-elts-check:checked').each(function() {
        // Needs to be Object so it can be the same structure as "search_genes.py" so it fits nicesly in populate_search_result_list()
        const label = $(this).data('label');
        selected_projections[label] = label;
    });
    const selected_projections_string = Object.keys(selected_projections).join(',');

    if (multigene) {
        // MG enabled
        $('#search_results_scrollbox').hide();
        $('#multigene_search_indicator').show();
    } else {
        // MG disabled
        $('#search_results_scrollbox').show();
        $('#multigene_search_indicator').hide();
    }

    $('#search_results').empty();
    // show search results
    $('#search_results_c').removeClass('search_result_c_DISABLED');

    projection_source = $("#set_of_patterns").val() ? $("#set_of_patterns").val() : null;
    projection = true;

    // Add the patterns source to the history.
    add_state_history(selected_projections_string, projection_source);

    dataset_collection_panel.load_frames({dataset_id, multigene, projection});

    // Run ProjectR for the chosen pattern
    if (projection_source) {
        const scope = $("#projection_source").data('scope');

        search_results = selected_projections;

        // Implementing search_genes.py results without the CGI execution
        // ? Can we use the DIMRED_meta file to get annotation info?
        populate_search_result_list(selected_projections);
        $('#searching_indicator_c').hide();
        $('#intro_content').hide('fade', {}, 400);
        // auto-select the first match.  first <a class="list-group-item"
        const first_thing = $('#search_results a.list-group-item').first();
        select_search_result(first_thing, draw_display=false);
        set_scrollbar_props();

        dataset_collection_panel.reset_abort_controller();

        dataset_collection_panel.datasets.forEach(dataset => {
            dataset.run_projectR(projection_source, is_pca, scope)
                .then(() => {

                    if (dataset.projection_id) {
                        if (multigene) {
                            // 'entries' is array of gene_symbols
                            dataset.draw_mg({ gene_symbols: Object.keys(selected_projections) });
                        } else {
                            dataset.draw({ gene_symbol: first_thing.data('gene_symbol')});
                        }
                    } else {
                        if (dataset.display) dataset.display.clear_display();
                        dataset.show_no_match();
                    }

                })
                .catch(error => console.error(error));
        })
    }
    return false;   // keeps the page from not refreshing
})

$("#set_of_patterns").on('change', () => {
    if ($('#set_of_patterns').val()) {
        $.ajax({
            type: "POST",
            url: "./cgi/get_pattern_element_list.cgi",
            async: false,
            data: {
                'file_name': $('#set_of_patterns').val(),
            },
            dataType: "json",
            success: (data) => {
                const pattern_elements_tmpl = $.templates("#pattern_elements_tmpl");
                const pattern_elements_html = pattern_elements_tmpl.render(data);
                $("#projection_pattern_elements").html(pattern_elements_html);
                $("#projection_pattern_elements_c").show();

                // Check boxes for the elements that were found in the URL
                if (projection_patterns) {
                    projection_patterns.split(',').forEach((pattern) => {
                        $(`.js-projection-pattern-elts-check[data-label="${pattern}"]`).prop('checked', true);
                    });
                }

            },
            error(xhr, status, msg) {
                console.error(`Failed to load dataset list because msg: ${msg}`);
            }
        });
    }
});

// controls to enable user scrolling of results with mouse arrow
scrolling_results = false

$('body').click((event) => {
    scrolling_results = !$(event.target).closest('#search_results_c').length ? false : true;;
});

$(document).keydown((event) => {
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
$(document).on('click', '#gene_details_header, #gene_collapse_btn', () => {
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

// If #search_param_gene_cart (front page) was the setter, set #selected_gene_cart
$('#search_param_gene_cart').change(() => {
    $('.js-gene-cart').text($('#search_param_gene_cart').text());
    $('.js-gene-cart').val($('#search_param_gene_cart').val());
    $('.js-gene-cart-div').css("visibility", "visible");
});

// When a gene cart is selected, populate the gene search bar with its members
$('.js-gene-cart').change( function() {
    const gene_cart_id = $(this).val();
    const params = { session_id, gene_cart_id };
    const d = new $.Deferred(); // Causes editable to wait until results are returned

    if (typeof session_id !== 'undefined') {
        // Get the gene cart members and populate the gene symbol search bar
        $.ajax({
            url: './cgi/get_gene_cart_members.cgi',
            async: false,
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
                    $('.js-gene-symbols').val(gene_symbols);

                    // determine if searching for exact matches
                    exact_match = true;
                    set_exact_match(exact_match);
                    $(".js-exact-match img").tooltip('hide'); // Do not want it to autoshow since it most likely is not hovered over right now
                } else {
                    $('.js-gene-cart').text(oldValue);
                    $('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
                        '<button type="button" class="close close-alert" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
                        '<p class="alert-message"><strong>Oops! </strong> ' + data.error + '</p></div>').show();
                }
                d.resolve();
            }
        });
    } else {
        d.resolve();
    }
    return d.promise();
});

// Set the state history based on current conditions
function add_state_history(searched_entities, projection_source=null) {
    const state_info = {
        'multigene_plots': Number(multigene)
    };

    let state_url = "/index.html?"
        + `multigene_plots=${state_info.multigene_plots}`;

    // Currently dataset share id and layout id URL params are mutually exclusive
    if (dataset_id) {
        state_info.share_id = dataset_id;
        state_url += `&share_id=${dataset_id}`;
    }

    if (layout_id) {
        state_info.layout_id = layout_id;
        state_url += `&layout_id=${layout_id}`;
    }

    if (projection_source) {
        state_info.projection_patterns = searched_entities;
        state_url += `&projection_patterns=${searched_entities}`;
        state_info.projection_source = projection_source;
        state_url += `&projection_source=${projection_source}`;

    } else {
        // gene_cart_id automatically enables the exact_match and populates gene symbols,
        // so let's not crowd up the history with that.
        if (gene_cart_id) {
            state_info.gene_cart_share_id = gene_cart_id;
            state_url += `&gene_cart_share_id=${gene_cart_id}`;
        } else {
            state_info.gene_symbol_exact_match = Number(exact_match);
            state_url += `&gene_symbol_exact_match=${state_info.gene_symbol_exact_match}`;
            if (searched_entities) {
                state_info.gene_symbol = searched_entities;
                state_url += `&gene_symbol=${searched_entities}`;
            }
        }
    }

    // SAdkins - Should we have a separate history state for dataset share IDs?
    history.pushState(
        // State Info
        state_info,
        // State title
        projection_source ? "Projection Pattern Search" : "Gene search",
        // URL
        state_url
    )
}

function set_scrollbar_props() {
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
}

// automatically reloads dataset grid and resubmits gene search
function update_datasetframes_generesults() {
    function resubmit_gene_search() {
        return;
        // SAdkins - 1/6/22 - commenting out since it loads double the frames
        //$('#gene_search_form').trigger('submit');
    }

    $.when( resubmit_gene_search() ).done(() => {
        if (multigene){
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

function set_exact_match(is_enabled, show_tooltip=true) {
    let action = "hide";
    if (show_tooltip) {
        action = "show";
    }

    if (is_enabled) {
        $(".js-exact-match img").attr("src", "img/arrow_target_selected.png");
        $(".js-exact-match:visible img").attr('data-original-title', "Exact match (currently on).").tooltip(action); // Show is added so the old version of the tooltip does not persist
    } else {
        $(".js-exact-match img").attr("src", "img/arrow_target_unselected.png");
        $(".js-exact-match:visible img").attr('data-original-title', "Exact match (currently off).").tooltip(action);
    }
}

function set_multigene_plots(is_enabled, show_tooltip=true) {
    let action = "hide";
    if (show_tooltip) {
        action = "show";
    }

    if (is_enabled) {
        $(".js-multigene:visible img").attr('data-original-title', "Multigene displays enabled. Click to search for single-gene displays.").tooltip(action);
        $(".js-multigene img").attr('src', 'img/icons/multi-dna.svg');
    } else {
        $(".js-multigene:visible img").attr('data-original-title', "Single-gene displays enabled. Click to search for multigene displays.").tooltip(action);
        $(".js-multigene img").attr('src', 'img/icons/single-dna.svg');
    }
}

// Show gene-related container and options
function show_gene() {
    $('#gene_search_div').show();
    $('#projection_search_div').hide();
    $('#submit_search_projection').hide();
    $("#gene_tab").addClass("active");
    $("#projection_tab").removeClass("active");

}

// Show projection-related container and options
function show_projection() {
    $('#projection_search_div').show();
    $('#submit_search_projection').show();
    $('#gene_search_div').hide();
    $("#projection_tab").addClass("active");
    $("#gene_tab").removeClass("active");
}


// Events to select and deselect all projection pattern checkboxes
$(document).on("click", "#projection_pattern_select_all", () => {
    $('.js-projection-pattern-elts-check').prop('checked', true);
});

$(document).on("click", "#projection_pattern_deselect_all", () => {
    $('.js-projection-pattern-elts-check').prop('checked', false);
});

// SAdkins - 5/6/22 - Moved these functions to top level so they are loaded quicker, so triggers work without needing a timeout period before.
$('#intro_search_form').on('submit', (e) => {
    // TODO: It makes sense to remove/destroy those elements we aren't showing after a search
    $('#intro_content').hide();

    $("#leftbar_main").show();
    $("#viewport_main").show();


    // fire the true search button, to submit the true form
    if (projection) {
        $("#submit_search_projection").trigger( "click" );
    } else {
        $("#search_gene_symbol").val( $("#search_gene_symbol_intro").val());
        $("#submit_search").trigger( "click" );
    }
    return false;   // prevent the default action
});

// Search from front page is clicked
$('#intro_search_icon').click((e) => {
    $('#intro_search_form').submit();
});

// Search from results page is clicked
$('#submit_search').click((e) => {
    // Reset some stuff before submission, so it does not show while AJAX stuff is happening
    $('#search_results').empty();
    $('#search_result_count').empty();
    $('#searching_indicator_c').show();

    // Scope selection
    $('#toggle_options').show();

    $('#gene_search_form').submit();
})

// Display curations using projections instead of genes
$('#submit_search_projection').click((e) => {
    // Reset some stuff before submission, so it does not show while AJAX stuff is happening
    $('#search_results').empty();
    $('#search_result_count').empty();
    $('#searching_indicator_c').show();

    // Scope selection
    $('#toggle_options').hide();  // Not sure if this is relevant for projections
    $('#projection_search_form').submit();
})