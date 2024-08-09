let impc_data;
let omim_data;
let hl_data;
let mgi_data;

function add_dga_button_listeners() {
    document.querySelectorAll('.btn-dm').forEach(item => {
        item.addEventListener('click', event => {
            const dm = item.closest('.dropdown');

            if (dm.classList.contains('is-active')) {
                dm.classList.remove('is-active');
            } else {
                dm.classList.add('is-active');
            }
        });
    });
}


function deafness_plugin_gene_change() {
    const gene_symbol = currently_selected_gene_symbol.toLowerCase();

    // remove the active class from all dropdowns
    document.querySelectorAll('.dropdown').forEach(item => {
        item.classList.remove('is-active');
    });

    // Disable the buttons (do we really need to?)
    document.getElementById("btn-deafness-gene-mouse").disabled = true;
    document.getElementById("btn-deafness-gene-human").disabled = true;
    document.getElementById("btn-deafness-gene-human-putative").disabled = true;

    // Reset the images to the outline states
    document.getElementById("img-deafness-gene-mouse").src = "./img/icons/org-1-outline-64.svg";
    document.getElementById("img-deafness-gene-human").src = "./img/icons/org-2-outline-64.svg";
    document.getElementById("img-deafness-gene-human-putative").src = "./img/icons/org-2-unknown-outline-64.svg";

    // Clear any existing lists
    document.querySelector(".dm-deafness-phenotypes").innerHTML = "";
    document.querySelector(".dm-deafness-links").innerHTML = "";

    const phenotype_template = document.getElementById('tmpl-deafness-phenotype');
    const link_template = document.getElementById('tmpl-deafness-resource-link');

    if (mgi_data.hasOwnProperty(gene_symbol)) {
        for (const phenotype of mgi_data[gene_symbol]['phenotypes']) {
            const span = phenotype_template.content.cloneNode(true);
            span.querySelector('span').innerHTML = phenotype;
            document.getElementById('dm-deafness-gene-mouse-phenotypes').appendChild(span);
        }

        for (const link of mgi_data[gene_symbol]['links']) {
            const a = link_template.content.cloneNode(true);
            a.querySelector('a').innerHTML = link['label'];
            a.querySelector('a').href = link['url'];
            document.getElementById('dm-deafness-gene-mouse-links').appendChild(a);
        }

        document.getElementById("img-deafness-gene-mouse").src = "./img/icons/org-1-dark-64.svg";
        document.getElementById("btn-deafness-gene-mouse").disabled = false;
    }
    
    return;

    

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
document.getElementById("annotation-panel-external-links-c").classList.remove("is-12");
document.getElementById("annotation-panel-external-links-c").classList.add("is-7");

// Add the contents of #deafness_gene_annotation_html_c after the first column under #annotation-panel
const external_links_c = document.getElementById('annotation-panel-external-links-c');
const deafness_gene_annotation_html = document.getElementById('deafness_gene_annotation_html_c').innerHTML;
external_links_c.insertAdjacentHTML('afterend', deafness_gene_annotation_html);

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
    
geneChangeCallbacks.push(deafness_plugin_gene_change);
geneChangeCallbacks.push(add_dga_button_listeners);
