#!/usr/bin/env mocha

/*
    Unit tests for dataset_curator.js
*/

// Some setup help by https://publishing-project.rivendellweb.net/testing-front-end-with-mocha-and-playwright-2/

import "mocha" // to set global variables like describe, it, etc.
import { expect } from "playwright/test";

import { browsers, gearBase, login, loginFailure, setupBrowserContext, teardownBrowserContext } from "./helpers.js";
import { mockGetDatasetAnalyses, mockGetDatasetAvailableDisplayTypes, mockGetDatasetGenes, mockGetDatasetAggregations } from "./helpers.js";

import { mockGetDatasetList } from "./helpers.js";

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

            it("should populate the dataset list", async () => {
                await mockGetDatasetList(page);
                const datasetInput = page.locator("css=#dataset_query");
                const datasetTree = page.locator("css=#dataset_tree");
                await expect(datasetTree).toBeVisible();
                await expect(datasetTree.locator("css=.wb-row")).toHaveCount(4);    // 3 categories + root
                await datasetInput.fill("Hertzano");
                await expect(page.getByText("P2, mouse, scRNA-seq, cochlea")).toBeVisible();

            });

            it("Plot button should be disabled if a dataset is not selected", async () => {

            });

            it("Plot button should be disabled if a plot type is not selected", async () => {

            });

            it("Plot button should be disabled if a gene is not selected", async () => {
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
                    await mockGetDatasetAvailableDisplayTypes(page);

                    await page.getByText("Curate new display").click();

                    // Should now be on Plot Type and Analysis select part
                    await page.getByText("Choose how to plot").click();
                });

                it.only("should plot a scatterplot", async () => {
                    const plotSelect = page.locator("css=#plot_type_select_c .nice-select-dropdown .list");
                    await plotSelect.getByText("Scatter").click();

                    await mockGetDatasetAggregations(page);

                    await page.getByText("Please select a gene").click();
                    const geneSelect = page.locator("css=#gene_select_c .nice-select-dropdown .list");
                    await geneSelect.getByText("Pou4f3").click();
                });

                it("should plot a bar plot", async () => {
                    const plotSelect = page.locator("css=#plot_type_select_c .nice-select-dropdown .list");
                    await plotSelect.getByText("Bar").click();

                    await mockGetDatasetAggregations(page);

                    await page.getByText("Please select a gene").click();
                    const geneSelect = page.locator("css=#gene_select_c .nice-select-dropdown .list");
                    await geneSelect.getByText("Pou4f3").click();
                });

                it("should plot a line plot", async () => {
                    const plotSelect = page.locator("css=#plot_type_select_c .nice-select-dropdown .list");
                    await plotSelect.getByText("Line").click();

                    await mockGetDatasetAggregations(page);

                    await page.getByText("Please select a gene").click();
                    const geneSelect = page.locator("css=#gene_select_c .nice-select-dropdown .list");
                    await geneSelect.getByText("Pou4f3").click();
                });

                it("should plot a violin plot", async () => {
                    const plotSelect = page.locator("css=#plot_type_select_c .nice-select-dropdown .list");
                    await plotSelect.getByText("Violin").click();

                    await mockGetDatasetAggregations(page);

                    await page.getByText("Please select a gene").click();
                    const geneSelect = page.locator("css=#gene_select_c .nice-select-dropdown .list");
                    await geneSelect.getByText("Pou4f3").click();
                });

                it("After plotting a violin plot with no color series, ensure color selection dropdown is disabled", async () => {

                });

                it("should plot a tSNE static plot", async () => {
                    const plotSelect = page.locator("css=#plot_type_select_c .nice-select-dropdown .list");
                    await plotSelect.getByText("tSNE static").click();

                    await mockGetDatasetAggregations(page);

                    await page.getByText("Please select a gene").click();
                    const geneSelect = page.locator("css=#gene_select_c .nice-select-dropdown .list");
                    await geneSelect.getByText("Pou4f3").click();
                });

                it("should plot a UMAP static plot", async () => {
                    const plotSelect = page.locator("css=#plot_type_select_c .nice-select-dropdown .list");
                    await plotSelect.getByText("UMAP static").click();

                    await mockGetDatasetAggregations(page);

                    await page.getByText("Please select a gene").click();
                    const geneSelect = page.locator("css=#gene_select_c .nice-select-dropdown .list");
                    await geneSelect.getByText("Pou4f3").click();
                });

                it("tSNE dynamic plot should correctly replot colors in plot when editing parameters after plotting", async () => {

                });

            });

            // TODO: Save SVG template as artifact for "local" testing
            it.skip("should plot a SVG plot", async () => {
                // Get dataset with an SVG image
                await mockGetDatasetList(page);
                const datasetInput = page.locator("css=#dataset_query");
                await datasetInput.fill("Hertzano");
                const datasetTree = page.locator("css=#dataset_tree");
                await datasetTree.getByText("P0, mouse, RNA-seq, hair cells vs epithelial non-hair cells (Hertzano)").click();
                await mockGetDatasetGenes(page);
                await mockGetDatasetAnalyses(page);
                await mockGetDatasetAvailableDisplayTypes(page);

                await page.getByText("Curate new display").click();

                // Should now be on Plot Type and Analysis select part
                await page.getByText("Choose how to plot").click();
                const plotSelect = page.locator("css=#plot_type_select_c .nice-select-dropdown .list");
                await plotSelect.getByText("SVG image").click();

                await mockGetDatasetAggregations(page);

                await page.getByText("Please select a gene").click();
                const geneSelect = page.locator("css=#gene_select_c .nice-select-dropdown .list");
                await geneSelect.getByText("Pou4f3").click();

                await page.locator("css=#plot_btn").click();

                // Check that the SVG image is displayed
                const svg = page.locator("css=#plot_container svg");
                await expect(svg).toBeVisible();
            });

            describe("options settings", () => {

                describe("common options", () => {

                });

                describe("plot-specific options", () => {

                });

            });


            describe("logged in", () => {

                beforeEach("Logging in", async () => {
                    await login(page, gearUrl);
                });

                it('should login', async () => {
                    expect(await page.title()).toEqual("gEAR dataset curator page");
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
                    expect(await page.title()).toEqual("gEAR dataset curator page");
                });

                it("should prevent user from selecting a default display", async () => {})

                it("should prevent user from saving a display", async () => {});

                it("should prevent user from deteting a display", async () => {});


            })

        });
    }
});