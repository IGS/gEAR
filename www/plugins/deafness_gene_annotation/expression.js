let impc_data;
let omim_data;
let hl_data;
let mgi_data;

function add_dga_button_listeners() {
    document.querySelectorAll('.btn-dm').forEach(item => {
        item.addEventListener('click', event => {
            let dm = item.closest('.dropdown');

            if (dm.classList.contains('is-active')) {
                dm.classList.remove('is-active');
            } else {
                dm.classList.add('is-active');
            }

            document.querySelectorAll('.btn-dm').forEach(btn => {
                if (btn !== item) {
                    dm = btn.closest('.dropdown');
                    dm.classList.remove('is-active');
                }
            });
        });
    });
}


function deafness_plugin_gene_change(selected_gene) {
    const gene_symbol = selected_gene.toLowerCase();

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
    for (const list of document.querySelectorAll(".dm-deafness-phenotypes")) {
        list.innerHTML = "";
    }

    for (const list of document.querySelectorAll(".dm-deafness-links")) {
        list.innerHTML = "";
    }

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

    if (omim_data.hasOwnProperty(gene_symbol)) {
        for (const phenotype of omim_data[gene_symbol]['phenotypes']) {
            const span = phenotype_template.content.cloneNode(true);
            span.querySelector('span').innerHTML = phenotype;
            document.getElementById('dm-deafness-gene-human-phenotypes').appendChild(span);
        }

        for (const link of omim_data[gene_symbol]['links']) {
            const a = link_template.content.cloneNode(true);
            a.querySelector('a').innerHTML = link['label'];
            a.querySelector('a').href = link['url'];
            document.getElementById('dm-deafness-gene-human-links').appendChild(a);
        }

        document.getElementById("img-deafness-gene-human").src = "./img/icons/org-2-dark-64.svg";
        document.getElementById("btn-deafness-gene-human").disabled = false;
    }

    if (hl_data.hasOwnProperty(gene_symbol)) {
        const locus = hl_data[gene_symbol]['locus'];
        const span = phenotype_template.content.cloneNode(true);
        span.querySelector('span').innerHTML = locus;
        document.getElementById('dm-deafness-gene-human-putative-loci').appendChild(span);

        for (const link of hl_data[gene_symbol]['links']) {
            const a = link_template.content.cloneNode(true);
            a.querySelector('a').innerHTML = link['label'];
            a.querySelector('a').href = link['url'];
            document.getElementById('dm-deafness-gene-human-putative-links').appendChild(a);
        }

        document.getElementById("img-deafness-gene-human-putative").src = "./img/icons/org-2-dark-64.svg";
        document.getElementById("btn-deafness-gene-human-putative").disabled = false;
    }
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

window.geneChangeCallbacks = [deafness_plugin_gene_change]

add_dga_button_listeners();
