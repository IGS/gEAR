#!/usr/bin/env mocha

/*
    Unit tests for gene_collection_manager.js
*/

// Some setup help by https://publishing-project.rivendellweb.net/testing-front-end-with-mocha-and-playwright-2/

import { describe } from 'node:test';
import { chromium, webkit, firefox, devices } from 'playwright';
import { expect } from "playwright/test";

const iPhone = devices['iPhone 12'];
const pixel = devices['Pixel 5'];
const browsers = {chromium, webkit, firefox, iPhone, pixel};

let browser, page, context;

const browserName = process.env.BROWSER || 'chromium';

// Check if devel is local or on server
const isLocal = process.env.LOCAL;
const userEmail = "sadkins@som.umaryland.edu";
const userPass = process.env.PASSWORD;

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

// NOTE https://mochajs.org/#arrow-functions

// NOTE async/await is not needed for finding the locator but rather for the interactions (i.e. click, type, etc.)
//  and some assertions that retry until they pass or timeout.

describe('Gene Collection Manager', function () {

    beforeEach(async function () {
        // Runs before each test
        browser = await browsers[browserName].launch({  // uncomment for debugging
            //headless: false,
            //slowMo: 1000
        });
        context = await browser.newContext();
        page = await context.newPage();
        await page.goto(gearUrl);
    });

    afterEach(async function () {
        // Runs after each test
        await browser.close();
    });

    it('Checks that page is returned', function () {
        // Ensure response has table element
        const controlFacet = page.locator("css=#controls_ownership");
        expect(controlFacet).toBeDefined();
    });

    describe("New Gene Collection", function () {
        it('should validate #new_collection_label input on blur', async function () {
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
    });

    describe("Existing Gene Collection", function () {

    });

    describe("Gene Collection Search", function () {
        it('should clear search terms on #search_clear click', async function () {
            // Set the #search_terms value
            await page.locator("css=#search_terms").fill("test");
            await page.keyboard.up("t")
            // Trigger the click event
            await page.locator("css=#search_clear").click();
            // Check that the #search_terms value was cleared
            await expect(page.locator("css=#search_terms")).toHaveValue("");
        });

        it('should search after enter clicked for #search_terms', async function () {
            // Trigger the keyup event
            await page.locator("css=#search_terms").type("test");
            await page.keyboard.up("Enter");
            // Check that submitSearch() was called if the key was "Enter"
            await expect(page.locator("css=#gc_count_label")).toContainText("results");
        });

        it('should update search results on #sort_by change', async function () {
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