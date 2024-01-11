#!/usr/bin/env mocha

/*
    Unit tests for multigene_curator.js
*/

// Some setup help by https://publishing-project.rivendellweb.net/testing-front-end-with-mocha-and-playwright-2/

import "mocha" // to set global variables like describe, it, etc.
import { expect } from "playwright/test";

import { browsers, gearBase, login, loginFailure, setupBrowserContext, teardownBrowserContext } from "./helpers.js";
import { mockGetDatasetAnalyses, mockGetDatasetAvailableMultigeneDisplayTypes, mockGetDatasetGenes, mockGetDatasetAggregations } from "./helpers.js";

import { mockGetDatasetList } from "./helpers.js";


let page;

const gearUrl = `${gearBase}/multigene_curator.html`;


describe('Multigene Curator', function () {
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
                    await expect(datasetTree).toBeVisible();
                    await expect(datasetTree.locator("css=.wb-row")).toHaveCount(4);    // 3 categories + root
                    await datasetInput.fill("Hertzano");
                    await expect(page.getByText("P2, mouse, scRNA-seq, cochlea")).toBeVisible();

                });

            });

            describe("gene selection", () => {

                it("should append to the genes list when a gene is manually selected", async () => {});

                it("should append to the genes list when a genecart is selected", async () => {});

                it("should remove a gene from the genes list when the remove button is clicked", async () => {});

                it("should remove all genes when Clear is clicked", async () => {});
            });

            describe("plotting", () => {

                beforeEach("Selecting a dataset", async () => {
                    await mockGetDatasetList(page);
                    const datasetInput = page.locator("css=#dataset_query");
                    await datasetInput.fill("Hertzano");
                    const datasetTree = page.locator("css=#dataset_tree");
                    await datasetTree.getByText("P2, mouse, scRNA-seq, cochlea (Hertzano/Ament)").click();

                    await mockGetDatasetGenes(page);
                    await mockGetDatasetAnalyses(page);
                    await mockGetDatasetAvailableMultigeneDisplayTypes(page);

                    await page.getByText("Curate new display").click();

                    // Should now be on Plot Type and Analysis select part
                    await page.getByText("Choose how to plot").click();
                });

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