#!/usr/bin/env mocha

/*
    Unit tests for dataset_curator.js
*/

// Some setup help by https://publishing-project.rivendellweb.net/testing-front-end-with-mocha-and-playwright-2/

import "mocha" // to set global variables like describe, it, etc.
import { expect } from "playwright/test";

import { browsers, gearBase, login, loginFailure, setupBrowserContext, teardownBrowserContext } from "./helpers.js";

import { mockGetDatasetList } from "./helpers.js";


let page;

const gearUrl = `${gearBase}/multigene_curator.html`;


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

            it("should populate the dataset list", async () => {
                it("should populate the dataset list", async () => {
                    await mockGetDatasetList(page);
                    const datasetInput = page.locator("css=#dataset_query");
                    const datasetTree = page.locator("css=#dataset_tree");
                    expect(await datasetTree).toBeVisible();
                    expect(await datasetTree.locator("css=.wb-row")).toHaveCount(3);    // 3 categories
                    await datasetInput.fill("Hertzano");
                    expect(await page.getByText("(Hertzano/Ament)")).toBeVisible();

                });

            });

            describe("gene selection", () => {

                it("should append to the genes list when a gene is manually selected", async () => {});

                it("should append to the genes list when a genecart is selected", async () => {});

                it("should remove a gene from the genes list when the remove button is clicked", async () => {});

                it("should remove all genes when Clear is clicked", async () => {});
            });

            describe("plotting", () => {

                it("should plot a heatmap", async () => {});

                it("should plot a matrixplot", async () => {});

                it("should plot a dotplot", async () => {});

                it("should plot a violin plot", async () => {});

                it("should plot a stacked violin plot", async () => {});

                it("should plot a volcano plot", async () => {});

                it("should plot a quadrant plot", async () => {});

            });

            describe("options settings", () => {

                describe("common options", () => {

                });

                describe("plot-specific options", () => {

                });

            });

            describe("saving", () => {

            });


            describe("logged in", () => {

                beforeEach("Logging in", async () => {
                    await login(page, gearUrl);
                });

                it('should login', async () => {
                    expect(await page.title()).toEqual("gEAR multigene viewer page");
                });

                it("should prevent user from deleting a display owned by another user", async () => {});

                it("should be able to delete a display owned by the user", async () => {});

                it("Should be able to switch default displays", async () => {});

                describe("saving", () => {

                    it("should save a new display", async () => {});

                    it("should save a new display as a new default", async () => {});
                });

            })

            describe("not logged in", () => {

                it('should fail to login with incorrect password', async () => {
                    await loginFailure(page, gearUrl);
                    // Check that the error message is displayed
                    const incorrectPw = page.getByText("Incorrect password");
                    await expect(incorrectPw).toBeVisible();
                    expect(await page.title()).toEqual("gEAR multigene viewer page");
                });

                it("should prevent user from selecting a default display", async () => {})

                it("should prevent user from saving a display", async () => {});

                it("should prevent user from deteting a display", async () => {});

            })

        });
    }
});