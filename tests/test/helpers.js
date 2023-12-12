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