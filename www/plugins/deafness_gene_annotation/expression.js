let impc_data;
let omim_data;
let hl_data;
let mgi_data;

function deafness_plugin_gene_change() {
    const gene_symbol = SELECTED_GENE.text().toLowerCase();

    // disable button#deafness_gene_mouse, button#deafness_gene_human, button#deafness_gene_human_putative
    document.getElementById("deafness_gene_mouse").disabled = true;
    document.getElementById("deafness_gene_human").disabled = true;
    document.getElementById("deafness_gene_human_putative").disabled = true;

    // reset all
    document.querySelector('#deafness_gene_mouse img').setAttribute('src', './img/icons/org-1-outline-64.svg');
    document.querySelector('#deafness_gene_human img').setAttribute('src', './img/icons/org-2-outline-64.svg');
    document.querySelector('#deafness_gene_human_putative img').setAttribute('src', './img/icons/org-2-unknown-outline-64.svg');
    document.querySelector('#deafness_gene_mouse').setAttribute('data-title', 'No data available');
    document.querySelector('#deafness_gene_human').setAttribute('data-title', 'No data available');
    document.querySelector('#deafness_gene_human_putative').setAttribute('data-title', 'No data available');
    document.querySelector('#deafness_popover_phenotypes').innerHTML = 'No data available';

    let mouse_deafness_links = [];

    if (mgi_data.hasOwnProperty(gene_symbol)) {
        mouse_deafness_links = mouse_deafness_links.concat(mgi_data[gene_symbol]['links']);
        document.querySelector('#deafness_gene_mouse').setAttribute('data-title', 'Click for links and phenotypes');

        document.querySelector('#deafness_popover_phenotypes').innerHTML = '';

        for (link of mouse_deafness_links) {

        }

        const phenotypes_tmpl = $.templates("#tmpl_deafness_phenotypes");
        const phenotypes_html = phenotypes_tmpl.render(mgi_data[gene_symbol]['phenotypes']);
        $("#deafness_popover_phenotypes").html(phenotypes_html);
    }

    if (impc_data.hasOwnProperty(gene_symbol)) {
        $('#deafness_gene_mouse').attr('data-title', impc_data[gene_symbol]['on_hover']);
        mouse_deafness_links = mouse_deafness_links.concat(impc_data[gene_symbol]['links']);
    }

    if (mouse_deafness_links.length > 0) {
        $("button#deafness_gene_mouse").attr("disabled", false);
        $('#deafness_gene_mouse img').attr('src', './img/icons/org-1-dark-64.svg');
        var links_tmpl = $.templates("#tmpl_deafness_resource_links");
        var links_html = links_tmpl.render(mouse_deafness_links);
        $("#deafness_popover_links").html(links_html);

        $('#deafness_gene_mouse').attr("data-popover", $("#deafness_popover_c").html())
    }

    // this is here after mouse so that the phenotypes don't carry over to other organisms
    $('#deafness_popover_phenotypes').html('No data available');

    if (omim_data.hasOwnProperty(gene_symbol)) {
        $("button#deafness_gene_human").attr("disabled", false);
        $('#deafness_gene_human img').attr('src', './img/icons/org-2-dark-64.svg');
        $('#deafness_gene_human').attr('data-title',
                                       omim_data[gene_symbol]['phenotypes'].join(' - '));

        var links_tmpl = $.templates("#tmpl_deafness_resource_links");
        var links_html = links_tmpl.render(omim_data[gene_symbol]['links']);
        $("#deafness_popover_links").html(links_html);
        $('#deafness_gene_human').attr("data-popover", $("#deafness_popover_c").html())
    }

    if (hl_data.hasOwnProperty(gene_symbol)) {
        $("button#deafness_gene_human_putative").attr("disabled", false);
        $('#deafness_gene_human_putative img').attr('src', './img/icons/org-2-unknown-dark-64.svg');
        $('#deafness_gene_human_putative').attr('data-title',
                                                hl_data[gene_symbol]['locus']);

        var links_tmpl = $.templates("#tmpl_deafness_resource_links");
        var links_html = links_tmpl.render(hl_data[gene_symbol]['links']);
        $("#deafness_popover_links").html(links_html);
        $('#deafness_gene_human_putative').attr("data-popover", $("#deafness_popover_c").html())
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


// First make the external resource links column smaller
document.getElementById("annotation-panel-external-links-c").classList.remove("is-9");
document.getElementById("annotation-panel-external-links-c").classList.add("is-5");

// Add the contents of #deafness_gene_annotation_html_c after the first column under #annotation-panel
const oselector_c = document.getElementById('annotation-panel-organism-selector-c');
const deafness_gene_annotation_html = document.getElementById('deafness_gene_annotation_html_c').innerHTML;
oselector_c.insertAdjacentHTML('afterend', deafness_gene_annotation_html);

// then destroy the initial placement
document.getElementById("deafness_gene_annotation_html_c").remove();

// get the data files
fetch("./plugins/deafness_gene_annotation/impc_data.json")
    .then(
        (response) => {
            if (response.status !== 200) {
                console.error('Looks like there was a problem. Status Code: ' +
                            response.status);
                return;
            }

            // Examine the text in the response
            response.json().then((data) => {
                impc_data = data;
            });
        }
    )
    .catch((err) => {
        console.error('Fetch Error :-S', err);
    });

fetch("./plugins/deafness_gene_annotation/omim_data.json")
    .then(
        (response) => {
            if (response.status !== 200) {
                console.error('Looks like there was a problem. Status Code: ' +
                            response.status);
                return;
            }

            // Examine the text in the response
            response.json().then((data) => {
                omim_data = data;
            });
        }
    )
    .catch((err) => {
        console.error('Fetch Error :-S', err);
    });

fetch("./plugins/deafness_gene_annotation/hl_loci.json")
    .then(
        (response) => {
            if (response.status !== 200) {
                console.error('Looks like there was a problem. Status Code: ' +
                            response.status);
                return;
            }

            // Examine the text in the response
            response.json().then((data) => {
                hl_data = data;
            });
        }
    )
    .catch((err) => {
        console.error('Fetch Error :-S', err);
    });

fetch("./plugins/deafness_gene_annotation/mgi_data.json")
    .then(
        (response) => {
            if (response.status !== 200) {
                console.error('Looks like there was a problem. Status Code: ' +
                            response.status);
                return;
            }

            // Examine the text in the response
            response.json().then((data) => {
                mgi_data = data;
            });
        }
    )
    .catch(function(err) {
        console.error('Fetch Error :-S', err);
    });

//search_result_postselection_functions.push(deafness_plugin_gene_change);
