var impc_data;
var omim_data;

function deafness_plugin_gene_change() {
    var gene_symbol = SELECTED_GENE.text().toLowerCase();

    // reset all
    $('#deafness_gene_mouse img').attr('src', './img/icons/org-1-outline-64.svg');
    $('#deafness_gene_human img').attr('src', './img/icons/org-2-outline-64.svg');
    $('#deafness_gene_human_putative img').attr('src', './img/icons/org-2-unknown-outline-64.svg');
    $('#deafness_gene_mouse').attr('data-title', 'No data available');
    $('#deafness_gene_mouse').attr('data-title', 'No data available');
    $('#deafness_gene_mouse').attr('data-title', 'No data available');

    if (impc_data.hasOwnProperty(gene_symbol)) {
        $('#deafness_gene_mouse img').attr('src', './img/icons/org-1-dark-64.svg');
        $('#deafness_gene_mouse').attr('data-title', impc_data[gene_symbol]['on_hover']);

        var links_tmpl = $.templates("#tmpl_deafness_resource_links");
        var links_html = links_tmpl.render(impc_data[gene_symbol]['links']);
        $("#deafness_popover_links").html(links_html);
        $('#deafness_gene_mouse').attr("data-popover", $("#deafness_popover_links").html())
    }

    if (omim_data.hasOwnProperty(gene_symbol)) {
        $('#deafness_gene_human img').attr('src', './img/icons/org-2-dark-64.svg');
        $('#deafness_gene_human').attr('data-title', 
                                       omim_data[gene_symbol]['phenotypes'].join(' - '));

        var links_tmpl = $.templates("#tmpl_deafness_resource_links");
        var links_html = links_tmpl.render(omim_data[gene_symbol]['links']);
        $("#deafness_popover_links").html(links_html);
        $('#deafness_gene_human').attr("data-popover", $("#deafness_popover_c").html())
    }

    $('button.icon-deafness-gene').each(function() {
        $(this).popover("dispose").popover({    
            content : $(this).attr("data-popover"),
            placement : 'auto',
            title : 'Deafness gene info',
            trigger: 'focus',
            html: true
        }).tooltip("dispose").tooltip({    
            placement : 'top',  
            title : $(this).attr("data-title")         
        }).on('show.bs.popover', function() {
            $(this).tooltip('hide')
        })
    });
}

$(document).ready(function() {
    // place the plugin HTML content
    $("#annotation_plugin_c").html(
        $("#deafness_gene_annotation_html_c").html()
    );

    // then destroy the initial placement
    document.getElementById("deafness_gene_annotation_html_c").remove();

    // get the data files
    fetch("./plugins/deafness_gene_annotation/impc_data.json")
        .then(
            function(response) {
                if (response.status !== 200) {
                    console.log('Looks like there was a problem. Status Code: ' +
                                response.status);
                    return;
                }

                // Examine the text in the response
                response.json().then(function(data) {
                    impc_data = data;
                });
            }
        )
        .catch(function(err) {
            console.log('Fetch Error :-S', err);
        });

    fetch("./plugins/deafness_gene_annotation/omim_data.json")
        .then(
            function(response) {
                if (response.status !== 200) {
                    console.log('Looks like there was a problem. Status Code: ' +
                                response.status);
                    return;
                }

                // Examine the text in the response
                response.json().then(function(data) {
                    omim_data = data;
                });
            }
        )
        .catch(function(err) {
            console.log('Fetch Error :-S', err);
        });

    search_result_postselection_functions.push(deafness_plugin_gene_change);
});
