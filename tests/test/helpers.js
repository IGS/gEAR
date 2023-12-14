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
const testEnv = process.env.TEST_ENV?.toLowerCase();

// Control debug output to make it easier to see what's going on
const isDebug = process.env.DEBUG?.toLowerCase() === "true";

// Set the gEAR URL based on the testing environment
const isLocal = testEnv === "local";
const isNode = testEnv === "node";
const isDevel = testEnv === "devel";
let gearBase;
switch (true) {
    case isLocal:
        console.log("Testing locally");
        gearBase = "http://localhost:8080";
        break;

    case isNode:
        console.log("Testing on node");
        gearBase = "http://localhost:3000";
        break;

    case isDevel:
        console.log("Testing on devel");
        gearBase = "https://devel.umgear.org";
        break;

    default:
        console.error("Unknown testing environment");
        throw new Error("Unknown testing environment");
        break;
}

export let browser, context;

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

export const mockGetOrganismList = async (page, organism) => {
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
export const mockGetDatasetList = async (page) => {
    await page.route(`${gearBase}/cgi/get_h5ad_dataset_list.cgi`, async route => {
        const json = {
            "shared_with_user": {
                "datasets": []
            },
            "user": {
                "datasets": [
                    {
                        "id": "dcfb4818-5302-83b5-4cc1-e586a5f81f74",
                        "owner_id": 662,
                        "title": "RNAseq, Deiters' cells,pillar cells, Inner Hair Cells and Outer Hair Cells,workshop (He)",
                        "organism_id": 1,
                        "pubmed_id": "30327589",
                        "geo_id": "GSE111347",
                        "is_public": null,
                        "ldesc": "RNA extracted from Deiters' cells,pillar cells, Inner Hair Cells and Outer Hair Cells. three  biological replicates, each with two technical repeats for Deiter and pillar cells, two biological replicates of IHCs and three biological replicates of OHCs, each with two technical repeats. They was seuqenced by Illumina HiSeq 2500.",
                        "date_added": "2022-05-25 15:21:30",
                        "dtype": "bulk-rnaseq",
                        "schematic_image": null,
                        "share_id": "711a8330",
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
                        "id": "1b12dde9-1762-7564-8fbd-1b07b750505f",
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