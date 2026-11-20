/* Various functions and constants that will be used across multiple pages */

import { chromium, webkit, firefox, devices } from 'playwright';

const iPhone = devices['iPhone 12'];
const pixel = devices['Pixel 5'];
export let browsers = [
    {name: "chromium", launch: chromium, browserContext: null},
    {name: "webkit", launch: webkit, browserContext: null},
    {name: "firefox", launch: firefox, browserContext: null},
    {name: "iPhone", launch: webkit, browserContext: {...iPhone}},
    {name: "pixel", launch: chromium, browserContext: {...pixel}}
];

// If BROWSER is set then only run that browser
if (process.env.BROWSER) {
    browsers = browsers.filter(b => b.name === process.env.BROWSER);
}

// Check if devel is local or on server
const isLocal = process.env.LOCAL?.toLowerCase() === "true";

// Control debug output to make it easier to see what's going on
const isDebug = process.env.DEBUG?.toLowerCase() === "true";

export const gearBase = isLocal ? "http://localhost:8080" : "https://devel.umgear.org";

export let browser, context;


const datasetId = "00000000-0000-0000-0000-000000000000";   // dummy uuid for P2, mouse, scRNA-seq, cochlea (Hertzano/Ament)


/**
 * Sets up the browser context for testing.
 * @param {number} browserIndex - The index of the browser to use.
 * @param {Object} browser - The browser object.
 * @param {Object} context - The context object.
 * @returns {Object} - The page object.
 */
export const setupBrowserContext = async (browserIndex) => {
    // Runs before each test
    browser = await browsers[browserIndex].launch.launch({  // uncomment for debugging
        headless: !isDebug,
        slowMo: isDebug ? 1000 : undefined,
    });
    context = await browser.newContext(browsers[browserIndex].browserContext);
    return await context.newPage()
}

/**
 * Closes the browser context and browser after each test.
 * @param {puppeteer.Browser} browser - The Puppeteer browser instance.
 * @param {puppeteer.BrowserContext} context - The Puppeteer browser context instance.
 * @returns {Promise<void>} - A promise that resolves when the browser context and browser are closed.
 */
export const teardownBrowserContext = async () => {
    // Runs after each test
    await context.close();
    await browser.close();
}

/**
 * Performs a mock login on the given page using the provided gear URL.
 * @param {Page} page - The page object to perform the login on.
 * @param {string} gearUrl - The gear URL to navigate to after login.
 * @returns {Promise<void>} - A promise that resolves once the login is completed.
 */
export const login = async (page, gearUrl) => {
    // Mock login
    await page.route(`${gearBase}/cgi/login.v2.cgi`, async route => {
        const json = {
            session_id : "1234567890ABCDEF",  // validation checks for >10 length
            user_name : "Test Armstrong",
            is_admin : 0,
            gear_default_domain : "Hearing (site default)",
            default_profile_share_id : null
        };
        await route.fulfill({ json });
        // Login refreshes and get_session_info.v2.cgi is called which overrides the user_name
    });

    await page.route(`${gearBase}/cgi/get_session_info.v2.cgi`, async route => {
        const json = {'email': "test@testing.com"
        , 'user_name': "Test Armstrong"
        , 'is_admin': 0
        , 'id':1
        , 'institution': "Test Institution"
        , 'colorblind_mode': 0
        , 'updates_wanted': 0
        , 'success':1
        }
        await route.fulfill({ json });
    });


    // Hover over the login button
    await page.locator("css=.navbar-item.not-logged-in").hover();

    await page.locator("css=[name='user-email']").fill("test");
    await page.locator("css=[name='user-password']").fill("secret_test");
    await page.locator("css=#submit-login").click();    // page refreshes

    // gear_session_id is set in the cookie
    //await context.addCookies({name: "gear_session_id", value: "1234", domain: gearBase, path: "/"});

    // Wait for the page to load
    await page.waitForURL(gearUrl);
}

/**
 * Simulates a login failure by mocking the login process and overriding session information.
 * @param {Page} page - The page object representing the browser page.
 * @param {string} gearUrl - The URL of the gear.
 * @returns {Promise<void>} - A promise that resolves when the login failure process is complete.
 */
export const loginFailure = async (page, gearUrl) => {
    // Mock login
    await page.route(`${gearBase}/cgi/login.v2.cgi`, async route => {
        const json = {
            session_id : "-1"
        };
        await route.fulfill({ json });
        // Login refreshes and get_session_info.v2.cgi is called which overrides the user_name
    });

    await page.route(`${gearBase}/cgi/get_session_info.v2.cgi`, async route => {
        const json = {
            "email": null,
            "user_name": null,
            'success':0
        }
        await route.fulfill({ json });
    });

    // Hover over the login button
    await page.locator("css=.navbar-item.not-logged-in").hover();

    await page.locator("css=[name='user-email']").fill("test");
    await page.locator("css=[name='user-password']").fill("secret_test");
    await page.locator("css=#submit-login").click();    // page refreshes
    await page.waitForURL(gearUrl);
}

/**
 * Mocks the getOrganismList function by intercepting the network request and returning a predefined JSON response.
 * @param {Object} page - The page object.
 * @param {string} organism - The organism name.
 * @returns {Promise<void>} - A promise that resolves when the network request is intercepted and fulfilled.
 */
export const mockGetOrganismList = async (page) => {
    await page.route(`${gearBase}/cgi/get_organism_list.cgi`, async route => {
        const json = {
            "organisms": [
                {
                    "id": 5,
                    "label": "Chicken",
                    "genus": "Gallus",
                    "species": "gallus",
                    "strain": null,
                    "taxon_id": 9031,
                    "has_icon": 1
                },
                {
                    "id": 2,
                    "label": "Human",
                    "genus": "Homo",
                    "species": "sapiens",
                    "strain": "sapiens",
                    "taxon_id": 9606,
                    "has_icon": 1
                },
                {
                    "id": 7,
                    "label": "Marmoset",
                    "genus": "Callithrix",
                    "species": "jacchus",
                    "strain": "jacchus",
                    "taxon_id": 9483,
                    "has_icon": 1
                },
                {
                    "id": 1,
                    "label": "Mouse",
                    "genus": "Mus",
                    "species": "musculus",
                    "strain": null,
                    "taxon_id": 10090,
                    "has_icon": 1
                },
                {
                    "id": 6,
                    "label": "Rat",
                    "genus": "Rattus",
                    "species": "norvegicus",
                    "strain": null,
                    "taxon_id": 10116,
                    "has_icon": 1
                },
                {
                    "id": 3,
                    "label": "Zebrafish",
                    "genus": "Danio",
                    "species": "rerio",
                    "strain": null,
                    "taxon_id": 7955,
                    "has_icon": 1
                }
            ]
        };
        await route.fulfill({ json });
    });
}
/**
 * Mocks the behavior of the getDatasetList function.
 * @param {Object} page - The page object.
 * @returns {Promise<void>}
 */
export const mockGetDatasetList = async (page) => {
    await page.route(`${gearBase}/cgi/get_h5ad_dataset_list.cgi`, async route => {
        const json = {
            "shared_with_user": {
                "datasets": []
            },
            "user": {
                "datasets": [
                    {
                        "id": "7812a487-932b-32f7-2de7-33dd3155c849",
                        "owner_id": 3,
                        "title": "P0, mouse, RNA-seq, cochlea, hair cells compared with rest of cochlear duct (Groves)",
                        "organism_id": 1,
                        "pubmed_id": "25855195",
                        "geo_id": null,
                        "is_public": 1,
                        "ldesc": "(P0) Atoh1&lt;super&gt;A1GFP/A1GFP&lt;/super&gt; mice which express GFP in their hair cells were used. Cochlear ducts were dissociated and hair cells were separated from the rest of the cells in the cochlear duct. Gene expression was measured using RNA-seq. (Cai T, Jen HI, Kang H, Klisch TJ, Zoghbi HY and Groves AK (2015) Characterization of the transcriptome of nascent hair cells and identification of direct targets of the Atoh1 transcription factor. J Neurosci. 2015 Apr 8;35(14):5870-83)",
                        "date_added": "2016-02-09 09:27:49",
                        "dtype": "svg-expression",
                        "schematic_image": "datasets_uploaded/7812a487-932b-32f7-2de7-33dd3155c849.jpg",
                        "share_id": "2c88b23c-b0f6-4f92-b88d-791c4de63527",
                        "math_default": "raw",
                        "marked_for_removal": 0,
                        "load_status": "completed",
                        "has_h5ad": 1,
                        "platform_id": null,
                        "instrument_model": null,
                        "library_selection": null,
                        "library_source": null,
                        "library_strategy": null,
                        "contact_email": null,
                        "contact_institute": null,
                        "contact_name": null,
                        "annotation_source": null,
                        "plot_default": null,
                        "annotation_release": null,
                        "gene_count": null,
                        "obs_count": null,
                        "has_tarball": 0,
                        "displays": [],
                        "tags": [],
                        "layouts": [],
                        "links": []
                    }
                ]
            },
            "public": {
                "datasets": [
                    {
                        "id": "00000000-0000-0000-0000-000000000000",
                        "owner_id": 3,
                        "title": "P2, mouse, scRNA-seq, cochlea (Hertzano/Ament)",
                        "organism_id": 1,
                        "pubmed_id": "None",
                        "geo_id": "None",
                        "is_public": 1,
                        "ldesc": "Cochleae of P2 C57BL/6N wild type mice were dissociated. Gene expression was measured using the 10x scRNA-seq platform. This is an unpublished dataest from the laboratories of Drs. Hertzano and Ament from the University of Maryland School of Medicine (UMSOM).\n\nThis dataset has a limited number of cells. Hair cells were identified as a single cluster. Inner and outer hair cells cluster together. ",
                        "date_added": "2019-02-01 15:29:49",
                        "dtype": "single-cell-rnaseq",
                        "schematic_image": null,
                        "share_id": "2050da31",
                        "math_default": "raw",
                        "marked_for_removal": 0,
                        "load_status": "completed",
                        "has_h5ad": 1,
                        "platform_id": null,
                        "instrument_model": null,
                        "library_selection": null,
                        "library_source": null,
                        "library_strategy": null,
                        "contact_email": null,
                        "contact_institute": null,
                        "contact_name": null,
                        "annotation_source": null,
                        "plot_default": null,
                        "annotation_release": null,
                        "gene_count": null,
                        "obs_count": null,
                        "has_tarball": 0,
                        "displays": [],
                        "tags": [],
                        "layouts": [],
                        "links": []
                    }
                ]
            }
        }
        await route.fulfill({ json });
    });
};

/**
 * Mocks the behavior of the get_dataset_displays.cgi API endpoint.
 * @param {Object} page - The page object.
 * @returns {Promise<void>} - A promise that resolves when the mocking is complete.
 */
export const mockGetDatasetDisplays = async (page) => {
    await page.route(`${gearBase}/cgi/get_dataset_displays.cgi`, async route => {
        const json = {
                "user": [
                    {
                        "id": 2416,
                        "dataset_id": "00000000-0000-0000-0000-000000000000",
                        "label": "ugly test plot",
                        "plot_type": "scatter",
                        "plotly_config": {
                            "x_axis": "cell_type",
                            "y_axis": "raw_value",
                            "hide_x_labels": false,
                            "hide_y_labels": false,
                            "color_name": "cell_type",
                            "hide_legend": false,
                            "jitter": false,
                            "marker_size": "8",
                            "color_palette": "purp",
                            "reverse_palette": false,
                            "order": {
                                "cell_type": [
                                    "Crabp1+_Cells",
                                    "Dividing_Glial_Cells",
                                    "Dividing_Mes_Cells",
                                    "Fst+_Cells",
                                    "Glial_Cells",
                                    "Hair_Cells",
                                    "Medial_Interdental_Cells",
                                    "Mes_Cells_2/3",
                                    "Mes_Cells_1",
                                    "Oc90+_Cells",
                                    "Pf4+_Cells",
                                    "Pillar_Cells",
                                    "Rgs5+_Cells",
                                    "Supporting_Cells",
                                    "Vascular_Cells"
                                ]
                            },
                            "colors": {
                                "Crabp1+_Cells": "#1f77b4",
                                "Dividing_Glial_Cells": "#ff7f0e",
                                "Dividing_Mes_Cells": "#2ca02c",
                                "Fst+_Cells": "#d62728",
                                "Glial_Cells": "#9467bd",
                                "Hair_Cells": "#8c564b",
                                "Medial_Interdental_Cells": "#e377c2",
                                "Mes_Cells_2/3": "#7f7f7f",
                                "Mes_Cells_1": "#bcbd22",
                                "Oc90+_Cells": "#17becf",
                                "Pf4+_Cells": "#16537e",
                                "Pillar_Cells": "#b3590a",
                                "Rgs5+_Cells": "#1f701f",
                                "Supporting_Cells": "#961b1c",
                                "Vascular_Cells": "#684884"
                            },
                            "vlines": [],
                            "gene_symbol": "Six1"
                        }
                    },
                    {
                        "id": 2417,
                        "dataset_id": "00000000-0000-0000-0000-000000000000",
                        "label": "ugly multigene plot",
                        "plot_type": "heatmap",
                        "plotly_config": {
                            "primary_col": "cell_type",
                            "colorscale": "bluered",
                            "reverse_colorscale": false,
                            "distance_metric": "euclidean",
                            "matrixplot": true,
                            "center_around_zero": false,
                            "cluster_obs": false,
                            "cluster_genes": true,
                            "flip_axes": false,
                            "hide_obs_labels": false,
                            "hide_gene_labels": false,
                            "clusterbar_fields": [
                                "cell_type"
                            ],
                            "obs_filters": {
                                "cell_type": [
                                    "Crabp1+_Cells",
                                    "Dividing_Glial_Cells",
                                    "Dividing_Mes_Cells",
                                    "Fst+_Cells",
                                    "Glial_Cells",
                                    "Hair_Cells",
                                    "Medial_Interdental_Cells",
                                    "Mes_Cells_2/3",
                                    "Mes_Cells_1",
                                    "Oc90+_Cells",
                                    "Pf4+_Cells",
                                    "Pillar_Cells",
                                    "Rgs5+_Cells",
                                    "Supporting_Cells",
                                    "Vascular_Cells"
                                ]
                            },
                            "sort_order": {
                                "cell_type": [
                                    "Hair_Cells",
                                    "Supporting_Cells",
                                    "Pillar_Cells",
                                    "Mes_Cells_1",
                                    "Mes_Cells_2/3",
                                    "Dividing_Mes_Cells",
                                    "Fst+_Cells",
                                    "Glial_Cells",
                                    "Dividing_Glial_Cells",
                                    "Medial_Interdental_Cells",
                                    "Vascular_Cells",
                                    "Crabp1+_Cells",
                                    "Oc90+_Cells",
                                    "Pf4+_Cells",
                                    "Rgs5+_Cells"
                                ]
                            },
                            "gene_symbols": [
                                "Acbd7",
                                "Cabp2",
                                "Cd164l2",
                                "Cib2",
                                "Gfi1",
                                "Gng8",
                                "Myo6",
                                "Myo7a",
                                "Pcp4",
                                "Pou3f4",
                                "Pou4f3",
                                "Ptprq",
                                "Pvalb",
                                "Shtn1",
                                "Smpx",
                                "Sox2",
                                "Tmc1"
                            ]
                        }
                    }
                ],
                "owner": [
                    {
                        "id": 2423,
                        "dataset_id": "00000000-0000-0000-0000-000000000000",
                        "label": "volcano test case",
                        "plot_type": "volcano",
                        "plotly_config": {
                            "de_test_algo": "t-test",
                            "annot_nonsignificant": true,
                            "pvalue_threshold": "0.05",
                            "adj_pvals": true,
                            "lower_logfc_threshold": "-1",
                            "upper_logfc_threshold": "1",
                            "ref_condition": "cell_type;-;Supporting_Cells",
                            "query_condition": "cell_type;-;Hair_Cells",
                            "obs_filters": {},
                            "gene_symbols": [
                                "Acbd7",
                                "Cabp2",
                                "Cd164l2",
                                "Cib2",
                                "Gng8",
                                "Myo6",
                                "Myo7a",
                                "Pcp4",
                                "Pou4f3",
                                "Ptprq",
                                "Pvalb",
                                "Shtn1",
                                "Smpx",
                                "Tmc1"
                            ]
                        }
                    },
                    {
                        "id": 2425,
                        "dataset_id": "00000000-0000-0000-0000-000000000000",
                        "label": "static_tsne_test",
                        "plot_type": "tsne_static",
                        "plotly_config": {
                            "x_axis": "X_tsne_1",
                            "y_axis": "X_tsne_2",
                            "flip_x": false,
                            "flip_y": false,
                            "colorize_legend_by": "cell_type",
                            "max_columns": "3",
                            "skip_gene_plot": false,
                            "horizontal_legend": false,
                            "marker_size": null,
                            "order": {},
                            "obs_filters": {},
                            "colors": {
                                "Crabp1+_Cells": "#1f77b4",
                                "Dividing_Glial_Cells": "#ff7f0e",
                                "Dividing_Mes_Cells": "#2ca02c",
                                "Fst+_Cells": "#d62728",
                                "Glial_Cells": "#9467bd",
                                "Hair_Cells": "#8c564b",
                                "Medial_Interdental_Cells": "#e377c2",
                                "Mes_Cells_2/3": "#7f7f7f",
                                "Mes_Cells_1": "#bcbd22",
                                "Oc90+_Cells": "#17becf",
                                "Pf4+_Cells": "#16537e",
                                "Pillar_Cells": "#b3590a",
                                "Rgs5+_Cells": "#1f701f",
                                "Supporting_Cells": "#961b1c",
                                "Vascular_Cells": "#684884"
                            },
                            "gene_symbol": "Pou4f3"
                        }
                    }
                ]
        }
        await route.fulfill({ json });
    });
};

/**
 * Mocks the getDatasetGenes function.
 * @param {Page} page - The page object.
 * @param {string} datasetId - The dataset ID.
 * @returns {Promise<void>}
 */
export const mockGetDatasetGenes = async (page) => {
    await page.route(`${gearBase}/api/h5ad/${datasetId}/genes`, async route => {
        const json = {
            "success": 1,
            "gene_symbols": [
                "Atoh1",
                "Pou4f3",
                "Sox2",
            ]
        }
        await route.fulfill({ json });
    });
}

/**
 * Mocks the getDatasetAnalyses function.
 * @param {Page} page - The page object.
 * @param {string} datasetId - The dataset ID.
 * @returns {Promise<void>}
 */
export const mockGetDatasetAnalyses = async (page) => {
    await page.route(`${gearBase}/api/h5ad/${datasetId}/analyses`, async route => {
        const json = {
            "success": 1,
            "public":[],
            "private":[]
        }
        await route.fulfill({ json });
    });
}

/**
 * Mocks the API endpoint for retrieving available display types for a dataset.
 * @param {Page} page - The Puppeteer page object.
 * @param {string} datasetId - The ID of the dataset.
 * @returns {Promise<void>} - A promise that resolves when the API endpoint is mocked.
 */
export const mockGetDatasetAvailableDisplayTypes = async (page) => {
    await page.route(`${gearBase}/api/h5ad/${datasetId}/availableDisplayTypes`, async route => {
        const json = {
            "scatter": true,
            "tsne_static": true,
            "umap_static": true,
            "pca_static": true,
            "tsne/umap_dynamic": true,
            "bar": true,
            "violin": true,
            "line": true,
            "svg": true
        }
        await route.fulfill({ json });
    });
}

/**
 * Mocks the API endpoint for retrieving available multigene display types for a dataset.
 * @param {Page} page - The page object.
 * @returns {Promise<void>}
 */
export const mockGetDatasetAvailableMultigeneDisplayTypes = async (page) => {
    await page.route(`${gearBase}/api/h5ad/${datasetId}/mg_availableDisplayTypes`, async route => {
        const json = {
            "dotplot": true,
            "heatmap": true,
            "mg_violin": true,
            "volcano": true,
            "quadrant": true

        }
        await route.fulfill({ json });
    });
}

/**
 * Mocks the response for the h5ad info API call.
 * @param {Object} page - The page object.
 * @returns {Promise<void>}
 */
export const mockGetDatasetH5adInfo = async (page) => {
    await page.route(`${gearBase}/api/h5ad/${datasetId}`, async route => {
        const json = {
            "success": 1,
            "num_obs": 4194,
            "obs_columns": [
                "cell_type",
                "barcode",
                "tSNE_1",
                "tSNE_2",
                "louvain",
                "cell_type_colors",
                "X_tsne_1",
                "X_tsne_2"
            ],
            "obs_levels": {
                "cell_type": [
                    "Crabp1+_Cells",
                    "Dividing_Glial_Cells",
                    "Dividing_Mes_Cells",
                    "Fst+_Cells",
                    "Glial_Cells",
                    "Hair_Cells",
                    "Medial_Interdental_Cells",
                    "Mes_Cells_2/3",
                    "Mes_Cells_1",
                    "Oc90+_Cells",
                    "Pf4+_Cells",
                    "Pillar_Cells",
                    "Rgs5+_Cells",
                    "Supporting_Cells",
                    "Vascular_Cells"
                ],
                "louvain": [
                    "Crabp1+_Cells",
                    "Dividing_Glial_Cells",
                    "Dividing_Mes_Cells",
                    "Fst+_Cells",
                    "Glial_Cells",
                    "Hair_Cells",
                    "Medial_Interdental_Cells",
                    "Mes_Cells_2/3",
                    "Mes_Cells_1",
                    "Oc90+_Cells",
                    "Pf4+_Cells",
                    "Pillar_Cells",
                    "Rgs5+_Cells",
                    "Supporting_Cells",
                    "Vascular_Cells"
                ],
                "cell_type_colors": [
                    "#00B0F6",
                    "#00BA38",
                    "#00BCD8",
                    "#00BF7D",
                    "#00C0AF",
                    "#6BB100",
                    "#619CFF",
                    "#A3A500",
                    "#B983FF",
                    "#C99800",
                    "#E76BF3",
                    "#E58700",
                    "#F8766D",
                    "#FD61D1",
                    "#FF67A4"
                ]
            },
            "has_replicates": 1
        }
        await route.fulfill({ json });
    });
}

/**
 * Mocks the response for getting dataset aggregations.
 * @param {Page} page - The page object.
 * @returns {Promise<void>}
 */
export const mockGetDatasetAggregations = async (page) => {
    await page.route(`${gearBase}/api/h5ad/${datasetId}/aggregations`, async route => {
        const json = {
            "success": 1,
            "aggregations": [
                {
                    "name": "cell_type",
                    "count": "4194",
                    "items": [
                        {
                            "name": "Mes_Cells_2/3",
                            "count": 1032
                        },
                        {
                            "name": "Mes_Cells_1",
                            "count": 1009
                        },
                        {
                            "name": "Glial_Cells",
                            "count": 398
                        },
                        {
                            "name": "Medial_Interdental_Cells",
                            "count": 382
                        },
                        {
                            "name": "Supporting_Cells",
                            "count": 339
                        },
                        {
                            "name": "Crabp1+_Cells",
                            "count": 252
                        },
                        {
                            "name": "Dividing_Glial_Cells",
                            "count": 115
                        },
                        {
                            "name": "Dividing_Mes_Cells",
                            "count": 113
                        },
                        {
                            "name": "Oc90+_Cells",
                            "count": 113
                        },
                        {
                            "name": "Hair_Cells",
                            "count": 101
                        },
                        {
                            "name": "Fst+_Cells",
                            "count": 85
                        },
                        {
                            "name": "Pillar_Cells",
                            "count": 70
                        },
                        {
                            "name": "Pf4+_Cells",
                            "count": 63
                        },
                        {
                            "name": "Rgs5+_Cells",
                            "count": 62
                        },
                        {
                            "name": "Vascular_Cells",
                            "count": 60
                        }
                    ]
                },
                {
                    "name": "louvain",
                    "count": "4194",
                    "items": [
                        {
                            "name": "Mes_Cells_2/3",
                            "count": 1032
                        },
                        {
                            "name": "Mes_Cells_1",
                            "count": 1009
                        },
                        {
                            "name": "Glial_Cells",
                            "count": 398
                        },
                        {
                            "name": "Medial_Interdental_Cells",
                            "count": 382
                        },
                        {
                            "name": "Supporting_Cells",
                            "count": 339
                        },
                        {
                            "name": "Crabp1+_Cells",
                            "count": 252
                        },
                        {
                            "name": "Dividing_Glial_Cells",
                            "count": 115
                        },
                        {
                            "name": "Dividing_Mes_Cells",
                            "count": 113
                        },
                        {
                            "name": "Oc90+_Cells",
                            "count": 113
                        },
                        {
                            "name": "Hair_Cells",
                            "count": 101
                        },
                        {
                            "name": "Fst+_Cells",
                            "count": 85
                        },
                        {
                            "name": "Pillar_Cells",
                            "count": 70
                        },
                        {
                            "name": "Pf4+_Cells",
                            "count": 63
                        },
                        {
                            "name": "Rgs5+_Cells",
                            "count": 62
                        },
                        {
                            "name": "Vascular_Cells",
                            "count": 60
                        }
                    ]
                }
            ],
            "total_count": 4194
        }
        await route.fulfill({ json });
    })
}