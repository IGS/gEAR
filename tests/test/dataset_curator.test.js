#!/usr/bin/env mocha

/*
    Unit tests for dataset_curator.js
*/

// Some setup help by https://publishing-project.rivendellweb.net/testing-front-end-with-mocha-and-playwright-2/

import "mocha" // to set global variables like describe, it, etc.
import { expect } from "playwright/test";

import { browsers, gearBase, login, loginFailure, setupBrowserContext, teardownBrowserContext } from "./helpers.js";

let page;

const gearUrl = `${gearBase}/dataset_curator.html`;


describe('Dataset Curator', function () {
    this.retries(3);

    this.timeout(10000);    // default is 2000

    let browserIndex = 0;

    beforeEach("Setting up", async () => {

        page = await setupBrowserContext(browserIndex);
        await page.goto(gearUrl);
    });

    afterEach("Closing", async () => {
        await teardownBrowserContext();
    });

    for (const b of browsers) {

        describe(`Browser: ${b.name}`, () => {

        });
    }
});