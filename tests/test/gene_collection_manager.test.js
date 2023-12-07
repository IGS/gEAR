#!/usr/bin/env mocha

/*
    Unit tests for gene_collection_manager.js
*/

// Some setup help by https://publishing-project.rivendellweb.net/testing-front-end-with-mocha-and-playwright-2/

import "mocha" // to set global variables like describe, it, etc.
import { chromium, webkit, firefox, devices } from 'playwright';
import { expect } from "playwright/test";

const iPhone = devices['iPhone 12'];
const pixel = devices['Pixel 5'];
let browsers = [
    {name: "chromium", launch: chromium, browserContext: null},
    {name: "webkit", launch: webkit, browserContext: null},
    {name: "firefox", launch: firefox, browserContext: null},
    {name: "iPhone", launch: webkit, browserContext: {...iPhone}},
    {name: "pixel", launch: chromium, browserContext: {...pixel}}
];

let browser, page, context;

if (process.env.BROWSER) {
    browsers = browsers.filter(b => b.name === process.env.BROWSER);
}

// Check if devel is local or on server
const isLocal = process.env.LOCAL;
const userEmail = process.env.EMAIL;
const userPass = process.env.PASSWORD;

const skipLoginTests = !userEmail || !userPass;


const gearBase = isLocal ? "http://localhost:8080" : "https://devel.umgear.org";
const gearUrl = `${gearBase}/gene_collection_manager.html`;

const login = async (page) => {
    // Hover over the login button
    await page.locator("css=.navbar-item.not-logged-in").hover();

    await page.locator("css=[name='user-email']").fill(userEmail);
    await page.locator("css=[name='user-password']").fill(userPass);
    await page.locator("css=#submit-login").click();    // page refreshes
    await page.waitForURL(gearUrl);
}

// NOTE https://mochajs.org/#arrow-functions - If using "this.timeout" then don't use arrow functions

// NOTE async/await is not needed for finding the locator but rather for the interactions (i.e. click, type, etc.)
//  and some assertions that retry until they pass or timeout.

// NOTE Playwright prefers getting elements by role, name, or text since that is what users see.
//  However sometimes it's just easier to use the locator method to grab the css selector

describe('Gene Collection Manager', function () {

    //this.timeout(10000);    // default is 2000

    let browserIndex = 0;

    beforeEach("Setting up", async () => {

        // Runs before each test
        browser = await browsers[browserIndex].launch.launch({  // uncomment for debugging
            //headless: false,
            //slowMo: 1000
        });
        context = await browser.newContext(browsers[browserIndex].browserContext);
        page = await context.newPage();
        await page.goto(gearUrl);
    });

    afterEach("Closing", async () => {
        // Runs after each test
        await context.close();
        await browser.close();
    });

    for (const b of browsers) {

        describe(`Browser: ${b.name}`, () => {

            after("Switching browsers", () => {
                // Runs once after all tests in this block
                browserIndex++;
            });

            it('Checks that page is returned', () => {
                // Ensure response has table element
                const controlFacet = page.locator("css=#controls_ownership");
                expect(controlFacet).toBeDefined();
            });

            describe("New Gene Collection", () => {
                it("should not appear if not logged in", async () => {
                    // Check that the form is hidden
                    await expect(page.getByRole("button", {name: "Create new gene collection"})).not.toBeVisible();
                });

                it('should validate #new_collection_label input on blur', async function () {
                    if (skipLoginTests) this.skip();

                    await login(page);
                    await page.getByRole("button", {name: "Create new gene collection"}).click();
                    await page.getByRole("button", {name: "Paste list of genes"}).click();

                    // Trigger the blur event
                    const label = await page.locator("css=#new_collection_label");
                    await label.click();
                    await page.locator("css=#new_collection_ldesc").click();

                    // Check that the correct class was added/removed from #new_collection_label
                    // Check that the helper text was added/removed
                    await expect(label).toHaveClass(/is-danger/);
                    await expect(page.getByText("Please enter a value")).toBeVisible();
                });

                it('should validate #new_collection_label input on save', async function () {
                    if (skipLoginTests) this.skip();

                    await login(page);
                    await page.getByRole("button", {name: "Create new gene collection"}).click();
                    await page.getByRole("button", {name: "Paste list of genes"}).click();

                    // Trigger the click event
                    await page.locator("css=#btn_new_collection_save").click();

                    await expect(page.locator("css=#new_collection_label")).toHaveClass(/is-danger/);

                    const parent = page.locator("css=.control:has(#new_collection_label)")
                        .filter({ has: page.getByText("Please enter a value") });

                    await expect(parent).toBeVisible();
                });

                it("should validate #new_collection_pasted_genes input on save", async function () {
                    if (skipLoginTests) this.skip();

                    await login(page);
                    await page.getByRole("button", {name: "Create new gene collection"}).click();
                    await page.getByRole("button", {name: "Paste list of genes"}).click();

                    // Trigger the click event
                    await page.locator("css=#btn_new_collection_save").click();

                    await expect(page.locator("css=#new_collection_label")).toHaveClass(/is-danger/);

                    const parent = page.locator("css=.control:has(#new_collection_pasted_genes)")
                        .filter({ has: page.getByText("Please enter a value") });

                    await expect(parent).toBeVisible();
                });

                it("should validate #new_collection_file input on save", async function () {
                    if (skipLoginTests) this.skip();

                    await login(page);
                    await page.getByRole("button", {name: "Create new gene collection"}).click();
                    await page.locator("css=#btn_gc_upload_unweighted_list").click();

                    // Trigger the click event
                    await page.locator("css=#btn_new_collection_save").click();

                    await expect(page.locator("css=.file:has(#new_collection_file)")).toHaveClass(/is-danger/);
                    await expect(page.getByText("Please select a file")).toBeVisible();
                });


                it('should validate organism select on save', async function () {
                    if (skipLoginTests) this.skip();

                    await login(page);
                    await page.getByRole("button", {name: "Create new gene collection"}).click();
                    await page.getByRole("button", {name: "Paste list of genes"}).click();

                    // Trigger the click event
                    await page.locator("css=#btn_new_collection_save").click();

                    await expect(page.locator("css=.select:has(#new_collection_organism_id)")).toHaveClass(/is-danger/);
                    await expect(page.getByText("Please select an organism")).toBeVisible();
                });

                it("should hide #new_collection_form_c on cancel", async function () {
                    if (skipLoginTests) this.skip();
                    await login(page);
                    await page.getByRole("button", {name: "Create new gene collection"}).click();
                    await page.getByRole("button", {name: "Paste list of genes"}).click();

                    // Trigger the click event
                    await page.locator("css=#btn_new_collection_cancel").click();

                    // Check that the form is hidden
                    await expect(page.locator("css=#new_collection_form_c")).not.toBeVisible();
                });
            });

            describe("Existing Gene Collection", () => {

            });

            describe("Gene Collection Search", () => {
                it('should clear search terms on #search_clear click', async () => {
                    // Set the #search_terms value
                    await page.locator("css=#search_terms").fill("test");
                    await page.keyboard.up("t")
                    // Trigger the click event
                    await page.locator("css=#search_clear").click();
                    // Check that the #search_terms value was cleared
                    await expect(page.locator("css=#search_terms")).toHaveValue("");
                });

                it('should search after enter clicked for #search_terms', async () => {
                    // Trigger the keyup event
                    await page.locator("css=#search_terms").type("test");
                    await page.keyboard.up("Enter");
                    // Check that submitSearch() was called if the key was "Enter"
                    await expect(page.locator("css=#gc_count_label")).toContainText("results");
                });

                it('should update search results on #sort_by change', async () => {
                    // Perform a search
                    await page.locator("css=#search_terms").fill("test");
                    await page.keyboard.up("Enter");
                    // Trigger the change event
                    await page.locator("css=#sort_by").selectOption("Gene cart title");
                    // Check that the search results were updated
                    await expect(page.locator("css=#gc_count_label")).toContainText("results");
                });
            });

        });
    };
});