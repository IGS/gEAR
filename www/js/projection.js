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
let multipattern = false;  // Is this a multipattern search?

const dataset_collection_panel = new DatasetCollectionPanel();

/*
Tree properties for constructor:
treeDiv - Element to generate the tree structure on
storedValElt - Element to store text, vals, and data properties on (if not in a treeDiv descendant "dropdown-toggle" element)
*/

const selected_profile_tree = new ProfileTree({treeDiv: '#selected_profile_tree'});
const projection_source_tree = new ProjectionSourceTree({treeDiv: '#projection_source_tree'});

const search_result_postselection_functions = [];

$(document).on("handle_page_loading", () => {
    // Was a permalink found?
    dataset_id = getUrlParameter('share_id');
    let scope = "permalink";

    if (dataset_id) {
        hide_functional_panel();

        // validate the dataset_id. runs load_frames() on success
        validate_permalink(scope);
    } else {
        // layout_id is a share_id for the profile layout
        layout_id = getUrlParameter('layout_id');
        scope = "profile";
    }

    // Was help_id found?
    const help_id = getUrlParameter('help_id');
    if (help_id) {
        validate_help_id(help_id);
    }

    const permalinked_multipattern_plots = getUrlParameter('multipattern_plots');
    multipattern = (permalinked_multipattern_plots !== (null || undefined)
        && permalinked_multipattern_plots === "1");
    set_multipattern_plots(multipattern, false);

    // Repopulate projection information... projection_source URL loaded earlier
    const permalinked_projection_algo = getUrlParameter('projection_algo')
    const permalinked_projection_patterns = getUrlParameter('projection_patterns');

    if (permalinked_projection_algo) {
        $(`input[name='projection_algo'][value='${permalinked_projection_algo}']`).prop("checked",true);
    }

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

        // Check boxes for the elements that were found in the URL
        selected_projections_string.split(',').forEach((pattern) => {
            const escaped_pattern = $.escapeSelector(pattern);
            $(`.js-projection-pattern-elts-check[data-label="${escaped_pattern}"]`).prop('checked', true);
        });
    }

    $("#submit_search_projection").trigger( "click" );

    // If MG search icon is clicked, toggle parameters
    $('.js-multipattern').click(() => {
        multipattern = !multipattern;
        set_multipattern_plots(multipattern);
    });

    // If multi-pattern set, toggle multipattern
    $('input[name="projection_display_mode"]').change(() => {
        multipattern = $('#multi_pattern').is(':checked');
        set_multipattern_plots(multipattern, false);
    });

    // add post-page load listeners
    $( "#dataset_zoomed_zoom_out_control" ).click(() => {
        zoom_out_dataset();
    });

    // If a ProfileTree element is selected, this is changed and the new layout is set
    $(document).on('change', '#selected_profile', function() {
        dataset_collection_panel.set_layout($(this).data('profile-id'), $(this).data('profile-label'), true, multipattern, true);
        layout_id = $(this).data('profile-share-id');
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

//Check help_id is valid. For Forgotten Password
function validate_help_id(help_id) {
    $.ajax({
        url: './cgi/validate_help_id.cgi',
        type: 'POST',
        data: {'help_id': help_id},
        dataType: 'json'
    }).done((data) => {
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
    }).fail((jqXHR, textStatus, errorThrown) => {
        display_error_bar(`${jqXHR.status} ${errorThrown.name}`, 'Error validating help ID');
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
$(document).on('keyup', 'input#user_new_pass_2', () => {
    const pass_1 = $('input#user_new_pass_1').val();

    if ( $('input#user_new_pass_2').val() == pass_1 ) {
        $('button#save_user_new_pass').prop('disabled', false);
    } else {
        $('button#save_user_new_pass').prop('disabled', true);
    }
});

// Submit new password
$(document).on('click', 'button#save_user_new_pass', () => {
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
        dataType: 'json'
    }).done((data) => {
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
    }).fail((jqXHR, textStatus, errorThrown) => {
        display_error_bar(`${jqXHR.status} ${errorThrown.name}`, 'Error saving new user password');
    });//end ajax
});

function validate_permalink(scope) {
    // Works for dataset or layout-based share IDs, which is differentiated by scope
    $.ajax({
        url : './cgi/validate_share_id.cgi',
        type: "POST",
        data : { 'share_id': dataset_id, 'scope': scope },
        dataType:"json"
    }).done((data) => {
        if ( data.success != 1 ) {
            $('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
                '<button type="button" class="close close-alert" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
                '<p class="alert-message"><strong>Oops! </strong> ' + data.error + '</p></div>').show();
        }

        if (scope == 'permalink') {
            dataset_collection_panel.load_frames({dataset_id});
        }
    }).fail((jqXHR, textStatus, errorThrown) => {
        display_error_bar(`${jqXHR.status} ${errorThrown.name}`, 'Error validating share ID');
    });
}

async function load_layouts() {
    const layout_share_id = getUrlParameter('layout_id');

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
        let active_layout_id = null;
        let active_layout_label = null;
        const layout_types = ['domain', 'user', 'group', 'shared', 'public']

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

        // NOTE: Need to deep copy the carts so the same node is not referenced in multiple trees.
        selected_profile_tree.domainProfiles = deepCopy(layouts.domain);
        selected_profile_tree.userProfiles = deepCopy(layouts.user);
        selected_profile_tree.groupProfiles = deepCopy(layouts.group);
        selected_profile_tree.sharedProfiles = deepCopy(layouts.shared);
        selected_profile_tree.publicProfiles = deepCopy(layouts.public);
        selected_profile_tree.folders = deepCopy(data.folders);

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

        dataset_collection_panel.set_layout(active_layout_id, active_layout_label, false, multipattern);
    }).fail((jqXHR, textStatus, errorThrown) => {
        display_error_bar(`${jqXHR.status} ${errorThrown.name}`, 'Error loading layouts.');
    });

    selected_profile_tree.generateTree();
}

async function load_pattern_carts() {
    let permalink_cart_id = null
    let permalink_cart_label = null

    const cart_share_id = getUrlParameter('projection_source');

    if (!session_id) {
        console.info("User is not logged in. Patterns not loaded.");
        projection_source_tree.generateTree();
        return;
    }

    await $.ajax({
        url: './cgi/get_user_gene_carts.cgi',
        type: 'post',
        data: { 'session_id': session_id, "group_by_type":true, 'share_id': cart_share_id},
        dataType: 'json'
    }).done((data, textStatus, jqXHR) => {

        const gctypes = ["unweighted-list", "weighted-list"];
        const cart_types = ['domain', 'user', 'group', 'shared', 'public'];

        for (const gctype of gctypes) {
            const carts = {};
            for (const ctype of cart_types) {
                carts[ctype] = [];
                if (data[gctype][`${ctype}_carts`].length > 0) {
                    $.each(data[gctype][`${ctype}_carts`], (_i, item) => {

                        // If cart permalink was passed in, retrieve gene_cart_id for future use.
                        if (cart_share_id && item.share_id == cart_share_id) {
                            permalink_cart_id = item.share_id;
                            permalink_cart_label = item.label;
                        }


                        carts[ctype].push({value: item.share_id,    // Use share ID as it is used in the cart file's basename
                                            text: item.label,
                                            folder_id: item.folder_id,
                                            folder_label: item.folder_label,
                                            folder_parent_id: item.folder_parent_id,
                                            gctype
                                            });
                    });
                }
            }

            projection_source_tree[gctype].domainGeneCarts = carts.domain;
            projection_source_tree[gctype].userGeneCarts = carts.user;
            projection_source_tree[gctype].groupGeneCarts = carts.group;
            projection_source_tree[gctype].sharedGeneCarts = carts.shared;
            projection_source_tree[gctype].publicGeneCarts = carts.public;
            projection_source_tree[gctype].folders = data.folders || [];
        }

    })
    .fail((jqXHR, textStatus, errorThrown) => {
        display_error_bar(`${jqXHR.status} ${errorThrown.name}`, "Weighted gene carts not sucessfully retrieved.");
    });

    // If projection info was in URL, one of the above should have the JSTree element returned
    projection_source_tree.generateTree();

    // If a permalink was provided, set the value in the tree and search bar
    if (permalink_cart_id) {
        $("#projection_source").text(permalink_cart_label);
        $("#projection_source").val(permalink_cart_id);
        // At this point, the tree is generated but loading data attributes to the storedValElt does not occur until a node is selected.
        // So we need to manually set the data attribute for the first-pass.
        const tree_leaf = projection_source_tree.treeData.find(e => e.id === `genecart__${$("#projection_source").val()}`);
        $("#projection_source").data("gctype", tree_leaf.gctype);
        $("#projection_source").trigger('change');
    }
}

// If user changes, update genecart/profile trees
async function load_all_trees(){
    // Update dataset and genecart trees in parallel
    // Works if they were not populated or previously populated
    try {
        await Promise.allSettled([load_layouts(), load_pattern_carts()]);
    } catch (err) {
        console.error(err)
    }

    // NOTE: This will trigger again if the MutationObserver catches a login, but that may be acceptable.
    $(document).trigger("handle_page_loading");

}

// Handle direct documentation links
$(document).on('click', '#doc-link-choices li', function(){
    window.location.replace(`./manual.html?doc=${$(this).data('doc-link')}`);
});

// Sort numerically in situations where a alphabetical string is followed by a number
// i.e. PC1, PC2, etc.
// NOTE: Ignores the leading string altogether, so this still applies even if that is not consistent.
const customNumericSort = function (a, b) {
    return (Number(a.match(/(\d+)$/g)[0]) - Number((b.match(/(\d+)$/g)[0])));
}

function populate_search_result_list(data) {
    clear_search_result_info();
    // so we can display in sorted order.  javascript sucks like that.
    sorted_gene_syms = [];

    for (const key in data) {
        if (data.hasOwnProperty(key)) {
            sorted_gene_syms.push(key);
        }
    }

    // Source - https://stackoverflow.com/a/29180576
    sorted_gene_syms.sort(customNumericSort);
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
function select_search_result(elm, draw_display=true) {
    //TODO Prevents this function from being double-called by #gene_search_form.submit()
    // This might be replaced by my addition of DatasetCollectionPanel.search_performed attribute
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

    if (! multipattern) {
        // Get to top up- and down-regulated genes for each pattern if they exist.
        const escaped_gene_sym = $.escapeSelector(gene_sym);
        const top_up = $(`.js-projection-pattern-elts-check[data-label=${escaped_gene_sym}]`).data('top-up') || undefined;
        const top_down = $(`.js-projection-pattern-elts-check[data-label=${escaped_gene_sym}]`).data('top-down') || undefined;

        if (! (top_up === undefined)) {
            $("#highly_expressed_genes_card .card-header").text(`Pattern ${gene_sym}`);
            $("#highly_expressed_genes_card #top_up_genes .card-text").text(top_up);
            $("#highly_expressed_genes_card #top_down_genes .card-text").text(top_down);
            $("#highly_expressed_genes_card").show();
            $('#functional_not_supported_alert').hide();    // Hide functional support panel to clean up some screen real estate
        }
    }

    if (draw_display) {
        dataset_collection_panel.update_by_search_result(search_results[gene_sym]);
    }

    // call any plugin functions
    search_result_postselection_functions.forEach((f) => {f()})
}

function isNumeric(n) {
    return !isNaN(parseFloat(n)) && isFinite(n);
}

$('#search_results').on("click", "a", function(e) {
    e.preventDefault(); //prevent page scrolling to top
    $(this).blur(); //removes focus so active's purple coloring can show

    dataset_collection_panel.reset_abort_controller();

    dataset_collection_panel.datasets.forEach((dataset) => {
        if (dataset.projection_id)
        dataset.draw({gene_symbol: $(this).data("gene_symbol")});
    });
    select_search_result(this, draw_display=false);
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

$("#projection_search_form").submit((event) => {
    // ! Needs to be in 'submit' event to ensure "Enter" triggera are handled correctly as well
    // Reset some stuff before submission, so it does not show while AJAX stuff is happening
    clear_search_result_info();

    hide_functional_panel();

    // re-initialize any open tooltips
    $('[data-toggle="tooltip"], .tooltip').tooltip();

    // If front page MG was enabled, ensure correct radio button is selected
    if (multipattern) {
        $('#multi_pattern').click();
    } else {
        // If single gene, select all genes
        $(`.js-projection-pattern-elts-check`).prop('checked', true);
    }

    const projection_algorithm = $('[name="projection_algo"]:checked').val();

    // Get selected projections and add as state
    const selected_projections = {};
    $('.js-projection-pattern-elts-check:checked').each(function() {
        // Needs to be Object so it can be the same structure as "search_genes.cgi" so it fits nicesly in populate_search_result_list()
        const label = $(this).data('label');
        selected_projections[label] = label;
    });
    const selected_projections_string = Object.keys(selected_projections).join(',');

    // MG disabled
    $('#search_results_scrollbox').show();
    $('#multipattern_search_indicator').hide();
    if (multipattern) {
        // MG enabled
        $('#search_results_scrollbox').hide();
        $('#multipattern_search_indicator').show();
    }

    // show search results
    $('#search_results_c').removeClass('search_result_c_DISABLED');

    const projection_source = $("#projection_source").val() ? $("#projection_source").val() : null;
    projection = true;

    // Run ProjectR for the chosen pattern
    if (projection_source) {
        // Add the patterns source to the history.
        add_state_history(selected_projections_string, projection_source, projection_algorithm);

        dataset_collection_panel.load_frames({dataset_id, multipattern, projection});
        const gctype = $("#projection_source").data('gctype');

        search_results = selected_projections;

        // Implementing search_genes.cgi results without the CGI execution
        populate_search_result_list(selected_projections);
        // auto-select the first match.  first <a class="list-group-item"
        const first_thing = $('#search_results a.list-group-item').first();
        select_search_result(first_thing, draw_display=false);
        set_scrollbar_props();

        dataset_collection_panel.reset_abort_controller();
        dataset_collection_panel.datasets.map((dataset) => {
            run_projection(dataset, projection_source, projection_algorithm, gctype, selected_projections, first_thing);
        });
        return false;  // keeps the page from not refreshing
    }

    $("#functional_not_supported_alert").hide();    // otherwise it shows way at the bottom of the screen
    return false;   // keeps the page from not refreshing
})

$("#projection_source").on('change', (_event) => {
    // Hide previous genecart pattern list of results
    $('#search_results_scrollbox').hide();
    // Empty these since the old pattern vals may not be relevent to the current pattern.
    clear_search_result_info();
    $("#search_gene_symbol").empty();

    const gctype = $("#projection_source").data("gctype");

    $("#binary_algo_form_check").hide();
    $('#multi_pattern_group').show();
    if (gctype === "unweighted-list") {
        $("#binary_algo_form_check").show();
        $('#multi_pattern_group').hide();
    }

    $.ajax({
        type: "POST",
        url: "./cgi/get_pattern_element_list.cgi",
        async: false,   // No clue why this works but async/wait does not.. maybe it's the onchange event?
        data: {
            'source_id': $('#projection_source').val(),
            'scope': gctype
        },
        dataType: "json"
    }).done((data) => {
        const pattern_elements_tmpl = $.templates("#pattern_elements_tmpl");
        const pattern_elements_html = pattern_elements_tmpl.render(data);
        $("#projection_pattern_elements").html(pattern_elements_html);

        // If only one pattern, disable multi-pattern
        $("#multi_pattern").prop( "disabled", false )
        if (!data.length) {
            $("#single_pattern").click();
            $("#multi_pattern").prop( "disabled", true )
        }

        // Only show if multipattern is enabled. All projections are used in single gene search.
        if (multipattern) {
            $("#projection_pattern_elements_c").show();
        }
    }).fail((jqXHR, textStatus, errorThrown) => {
        display_error_bar(`${jqXHR.status} ${errorThrown.name}`, 'Error getting list of patterns from source.');
    });
});

$("#highly_expressed_genes_card .card-footer button").on("click", (_event) => {
    const pattern_id = $(SELECTED_GENE).data("gene_symbol");

    $.ajax({
        type: "POST",
        url: "./cgi/get_pattern_weighted_genes.cgi",
        async: false,   // No clue why this works but async/wait does not.. maybe it's the onchange event?
        data: {
            'source_id': $('#projection_source').val(),
            pattern_id
        },
        dataType: "json"
    }).done((data) => {
        let html_stream = "<table>";
        // Gather genes and weights and show in new page
        for (const row of data) {
            html_stream += `<tr><td>${row["gene"]}</td><td>${row["weight"]}</td></tr>`;
        }
        html_stream += "</table>"
        const tab = window.open('about:blank', '_blank');
        tab.document.write(html_stream);
        tab.document.close();

    }).fail((jqXHR, textStatus, errorThrown) => {
        display_error_bar(`${jqXHR.status} ${errorThrown.name}`, `Error getting gene weights for pattern ${pattern_id} from source.`);
    });

})

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

        const draw = projection ? false : true;

        switch (event.keyCode) {
            // up key
            case 38:
            if (AT_FIRST_MATCH_RECORD == false) {
                const new_gene = $(SELECTED_GENE).prev();
                // For projections, we do not want to go through checks where the organism_id is checked
                // (since there isn't one).  So we draw the display and select the search result independently
                if (!draw) {
                    dataset_collection_panel.datasets.forEach((dataset) => {
                        dataset.draw({gene_symbol: $(new_gene).data("gene_symbol")});
                    });
                }

                select_search_result($(SELECTED_GENE).prev(), draw_display=draw);
            }
            break;

            // down key
            case 40:
            if (AT_LAST_MATCH_RECORD == false) {
                const new_gene = $(SELECTED_GENE).next();
                if (!draw) {
                    dataset_collection_panel.datasets.forEach((dataset) => {
                        dataset.draw({gene_symbol: $(new_gene).data("gene_symbol")});
                    });
                }

                select_search_result($(SELECTED_GENE).next(), draw_display=draw);
            }
            break;
        }
    }
});

// If #search_param_gene_cart (front page) was the setter, set #selected_gene_cart
$('#search_param_gene_cart').change(() => {
    $('.js-gene-cart').text($('#search_param_gene_cart').text());
    $('.js-gene-cart').val($('#search_param_gene_cart').val());
    $('.js-gene-cart-div').css("visibility", "visible");
});

async function run_projection(dataset, projection_source, projection_algorithm, gctype, selected_projections, first_thing) {
    try {
        await dataset.run_projectR(projection_source, projection_algorithm, gctype);
        if (dataset.projection_id) {
            if (multipattern) {
                // 'entries' is array of gene_symbols
                dataset.draw_mg({ gene_symbols: Object.keys(selected_projections) });
            } else {
                dataset.draw({ gene_symbol: first_thing.data('gene_symbol') });
            }
        } else {
            if (dataset.active_display)
                dataset.active_display.clear_display();
            dataset.show_no_match();
        }

    } catch(error) { console.error(error)};
}


// Set the state history based on current conditions
function add_state_history(searched_entities, projection_source=null, projection_algorithm=null) {
    const state_info = {
        'multipattern_plots': Number(multipattern)
    };

    let state_url = "/projection.html?"
        + `multipattern_plots=${state_info.multipattern_plots}`;

    // Currently dataset share id and layout id URL params are mutually exclusive (or neither is used)
    if (dataset_id) {
        state_info.share_id = dataset_id;
        state_url += `&share_id=${dataset_id}`;
    }

    if (layout_id) {
        state_info.layout_id = layout_id;
        state_url += `&layout_id=${layout_id}`;
    }

    // If "transfer learning button on front page was clicked, this will initially be undefined
    if (projection_source) {
        state_info.projection_source = projection_source;
        state_url += `&projection_source=${projection_source}`;
    }

    if (projection_algorithm) {
        state_info.projection_algo = projection_algorithm;
        state_url += `&projection_algo=${state_info.projection_algo}`;
    }
    // If a single-projection display, do not bother adding the list of patterns
    if (state_info.multipattern_plots == 1) {
        state_info.projection_patterns = searched_entities;
        state_url += `&projection_patterns=${searched_entities}`;
    }

    // SAdkins - Should we have a separate history state for dataset share IDs?
    history.pushState(
        // State Info
        state_info,
        // State title
        "Projection Pattern Search",
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

function update_datasetframes_projections() {
    // TODO: Refactor this to not repeat so much of the projection search submit event

    // Assumes that projection source is still present.
    const projection_source = $("#projection_source").val();
    if (!projection_source) {
        return;
    }

    if (!multipattern) {
        // If single gene, select all genes
        $(`.js-projection-pattern-elts-check`).prop('checked', true);
    }

    // Get selected projections and add as state
    const selected_projections = {};
    $('.js-projection-pattern-elts-check:checked').each(function() {
        // Needs to be Object so it can be the same structure as "search_genes.cgi" so it fits nicesly in populate_search_result_list()
        const label = $(this).data('label');
        selected_projections[label] = label;
    });

    // Implementing search_genes.cgi results without the CGI execution
    populate_search_result_list(selected_projections);

    const gctype = $("#projection_source").data('gctype');
    const first_thing = $('#search_results a.list-group-item').first();
    select_search_result(first_thing, draw_display=false);

    const projection_algorithm = $('[name="projection_algo"]:checked').val();

    // Run ProjectR for the chosen pattern
    dataset_collection_panel.datasets.map((dataset) => {
        run_projection(dataset, projection_source, projection_algorithm, gctype, selected_projections, first_thing);
    });
}

function set_multipattern_plots(is_enabled, show_tooltip=true) {
    const action = show_tooltip ? "show" : "hide";

    if (is_enabled) {
        $(".js-multipattern:visible").attr('data-original-title', "multipattern displays enabled. Click to search for single-gene displays.").tooltip(action);
        $(".js-multipattern img").attr('src', 'img/icons/multi-dna.svg');
        // projection tab stuff
        $("#projection_pattern_deselect_all").click();
        $("#projection_pattern_elements_c").show();
        return;
    }
    $(".js-multipattern:visible").attr('data-original-title', "Single-gene displays enabled. Click to search for multipattern displays.").tooltip(action);
    $(".js-multipattern img").attr('src', 'img/icons/single-dna.svg');
    // projection tab stuff
    $("#projection_pattern_select_all").click();
    $("#projection_pattern_elements_c").hide();
}

function hide_functional_panel() {
    $('#functional_not_supported_alert').show();
    $('#highly_expressed_genes_card').hide();
}

function show_functional_panel() {
    $('#functional_not_supported_alert').hide();
    $('#highly_expressed_genes_card').hide();
}

function clear_search_result_info() {
    $('#search_results').empty();
    $('#search_result_count').empty();
}

// Events to select and deselect all projection pattern checkboxes
$(document).on("click", "#projection_pattern_select_all", () => {
    $('.js-projection-pattern-elts-check').prop('checked', true);
});

$(document).on("click", "#projection_pattern_deselect_all", () => {
    $('.js-projection-pattern-elts-check').prop('checked', false);
});

// Display curations using projections instead of genes
$('#submit_search_projection').click((e) => {
    $('#projection_search_form').submit();
})