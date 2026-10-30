#!/usr/bin/env mocha

/*
    Unit tests for cell_cell.js (Cell-Cell Communication viewer for DominoSignal outputs)
*/

// Some setup help by https://publishing-project.rivendellweb.net/testing-front-end-with-mocha-and-playwright-2/

import "mocha" // to set global variables like describe, it, etc.
import { expect } from "playwright/test";

import { browsers, gearBase, setupBrowserContext, teardownBrowserContext } from "./helpers.js";

let page;

const gearUrl = `${gearBase}/cell_cell.html`;

const DOMINO_FILE = "domino_test.parquet";
const DOMINO_FILE_NO_INFO = "domino_no_info.parquet";           // no matching .json file
const DOMINO_FILE_PARTIAL_INFO = "domino_partial_info.parquet"; // .json file without an author
const CLUSTERS = ["B cells", "Macrophages_Pos", "Pericytes", "T cells"];
const DOMINO_OPTIONS = {
    "tfs": ["Atf3", "Foxp1"],
    "ligands": ["Gdf11", "Tgfb1"],
    "receptors": ["Antxr1", "Egfr"],
    "sending_clusters": CLUSTERS,
    "receiving_clusters": CLUSTERS,
};

// Dataset details from each output's matching .json file, as returned by get_domino_options.cgi
const DATASET_INFO = {
    [DOMINO_FILE]: {
        "title": "A <i>test</i> study of cell-cell communication",   // markup should be shown as text, not HTML
        "author": "A Researcher",
        "year": "2024",
        "doi": "10.1000/test.2024",
        "pmid": null,
        "gear_url": "umgear.org",
        "description": "Test dataset",
    },
    [DOMINO_FILE_NO_INFO]: null,
    [DOMINO_FILE_PARTIAL_INFO]: {
        "title": "An incomplete record", "author": null, "year": "2023",
        "doi": null, "pmid": "12345", "gear_url": null, "description": null,
    },
};

// Form fields sent to get_domino_plots.cgi, one object per request (reset before each test)
let plotRequests = [];

/**
 * Gets a field from a request's multipart form data (as sent via convertToFormData).
 * @param {Request} request - The Playwright request.
 * @param {string} name - The form field name.
 * @returns {string|null} - The field value, or null if not sent.
 */
const getFormField = (request, name) => {
    const match = (request.postData() || "").match(new RegExp(`name="${name}"\\r\\n\\r\\n([^\\r]*)`));
    return match ? match[1] : null;
};

/**
 * Mocks get_domino_options.cgi: the list of DominoSignal outputs and, when a file is passed, its parameter options.
 * @param {Page} page - The page object.
 * @returns {Promise<void>}
 */
const mockGetDominoOptions = async (page) => {
    await page.route(`${gearBase}/cgi/get_domino_options.cgi`, async route => {
        const json = { "success": 1, "files": Object.keys(DATASET_INFO) };
        const filename = getFormField(route.request(), "filename");
        if (filename) {
            json.options = DOMINO_OPTIONS;
            json.dataset_info = DATASET_INFO[filename];
        }
        await route.fulfill({ json });
    });
};

// Small figures standing in for the real dot and river plots
const DOT_FIGURE = {
    data: [{ type: "scatter", mode: "markers", x: ["Pericytes->Macrophages_Pos"], y: ["Gdf11->Antxr1"], marker: { color: [1.2], size: 10 } }],
    layout: { title: { text: "Geometric Mean of Ligand and Receptor Expression" }, height: 500 },
};
const RIVER_FIGURE = {
    data: [{ type: "sankey", node: { label: ["Pericytes", "Gdf11", "Antxr1", "Macrophages_Pos"] }, link: { source: [0, 1, 2], target: [1, 2, 3], value: [1, 1, 1] } }],
    layout: { title: { text: "Cell-cell communication" }, height: 500 },
};

// The selection each plot needs, mirroring get_domino_plots.cgi's "can't generate" messages
const PLOT_REQUIREMENTS = {
    dot_tf: { fields: ["selected_tf"], message: "Select a TF to generate this plot" },
    dot_clusters: { fields: ["sending_clusters", "receiving_clusters"], message: "Select sending and/or receiving clusters to generate this plot" },
    dot_clusters_filtered: { fields: ["sending_clusters", "receiving_clusters"], message: "Select sending and/or receiving clusters to generate this plot" },
    river_tf: { fields: ["selected_tf"], message: "Select a TF to generate this plot" },
    river_ligand: { fields: ["selected_ligand"], message: "Select a ligand to generate this plot" },
    river_receptor: { fields: ["selected_receptor"], message: "Select a receptor to generate this plot" },
    river_selected: { fields: [], message: null },
    river_selected_notf: { fields: [], message: null },
};

/**
 * Mocks get_domino_plots.cgi.  Returns plot JSON (or a PDF when image_format=pdf), recording each request's form fields.
 * @param {Page} page - The page object.
 * @returns {Promise<void>}
 */
const mockGetDominoPlots = async (page) => {
    await page.route(`${gearBase}/cgi/get_domino_plots.cgi`, async route => {
        const request = route.request();
        const fields = Object.fromEntries(
            [...(request.postData() || "").matchAll(/name="([^"]+)"\r\n\r\n([^\r]*)/g)].map((m) => [m[1], m[2]])
        );
        plotRequests.push(fields);

        const plotTypes = JSON.parse(fields.plots || "[]");

        if (fields.image_format === "pdf") {
            await route.fulfill({
                status: 200,
                contentType: "application/pdf",
                headers: { "Content-Disposition": `attachment; filename="domino_test_${plotTypes[0]}.pdf"` },
                body: Buffer.from("%PDF-1.4 test"),
            });
            return;
        }

        const json = { "success": 1, "plots": {}, "messages": {} };
        for (const plotType of plotTypes) {
            const { fields: needed, message } = PLOT_REQUIREMENTS[plotType];
            if (needed.length && !needed.some((name) => fields[name])) {
                json.plots[plotType] = null;
                json.messages[plotType] = message;
            } else {
                json.plots[plotType] = plotType.startsWith("dot") ? DOT_FIGURE : RIVER_FIGURE;
            }
        }
        await route.fulfill({ json });
    });
};

/**
 * Chooses the mocked DominoSignal output and waits for its parameters to load.
 * @param {Page} page - The page object.
 * @returns {Promise<void>}
 */
const chooseDominoFile = async (page, filename = DOMINO_FILE) => {
    await page.getByTestId("domino-file-select").selectOption(filename);
    await expect(page.getByTestId("parameters-panel")).toBeVisible();
    await expect(page.getByTestId("show-visualizations-button")).toBeEnabled();
};

/**
 * The searchable dropdown widget (nice-select2) shown in place of a hidden native select.
 * @param {Page} page - The page object.
 * @param {string} testId - The data-testid of the native select.
 * @returns {Locator}
 */
const dropdown = (page, testId) => page.getByTestId(testId).locator("xpath=following-sibling::div[contains(@class, 'nice-select')][1]");

/**
 * Selects values in a multiple dropdown by clicking them in its list, as a user would, then closes it.
 * @param {Page} page - The page object.
 * @param {string} testId - The data-testid of the native select.
 * @param {string[]} values - The values to select.
 * @returns {Promise<void>}
 */
const selectInDropdown = async (page, testId, values) => {
    const widget = dropdown(page, testId);
    await widget.click();
    for (const value of values) {
        await widget.locator(`li.option[data-value="${value}"]`).click();
    }
    // Clicking outside closes a multiple dropdown
    await page.locator("#page-header-label").click();
};

/**
 * The values currently selected in a multiple dropdown (tracked by NiceSelect via "selected" attributes).
 * @param {Page} page - The page object.
 * @param {string} testId - The data-testid of the native select.
 * @returns {Promise<string[]>}
 */
const getSelectedValues = async (page, testId) => {
    return await page.getByTestId(testId).evaluate((select) => Array.from(select.querySelectorAll("option[selected]")).map((option) => option.value));
};

/**
 * Sets which plots are requested, by ticking only the given plot checkboxes.
 * @param {Page} page - The page object.
 * @param {string[]} plotTypes - The plot types to request.
 * @returns {Promise<void>}
 */
const choosePlots = async (page, plotTypes) => {
    for (const checkbox of await page.locator("[data-testid^='plot-toggle-']").all()) {
        await checkbox.setChecked(plotTypes.includes(await checkbox.getAttribute("value")));
    }
};

describe('Cell-Cell Communication', function () {
    this.retries(3);
    this.timeout(20000);    // default is 2000

    let browserIndex = 0;

    beforeEach("Setting up", async () => {
        plotRequests = [];
        page = await setupBrowserContext(browserIndex);

        // The output list is fetched as the page loads, so mock before navigating
        await mockGetDominoOptions(page);
        await mockGetDominoPlots(page);
        await page.goto(gearUrl);
    });

    afterEach("Closing", async () => {
        await teardownBrowserContext();
    });

    for (const b of browsers) {

        describe(`Browser: ${b.name}`, () => {

            after("Switching browsers", () => {
                // Runs once after all tests in this block
                browserIndex++;
            });

            it("should list the DominoSignal outputs, with parameters hidden until one is chosen", async () => {
                const fileSelect = page.getByTestId("domino-file-select");
                await expect(fileSelect.locator("option")).toHaveText(["Choose a DominoSignal output", ...Object.keys(DATASET_INFO)]);
                await expect(page.getByTestId("parameters-panel")).toBeHidden();
                await expect(page.getByTestId("dataset-info")).toBeHidden();
                await expect(page.getByTestId("results-initial-notification")).toBeVisible();
                await expect(page.getByTestId("show-visualizations-button")).toBeDisabled();
            });

            it("should load the parameter options for the chosen output", async () => {
                await chooseDominoFile(page);

                // Single selects have a "None" option first
                await expect(page.getByTestId("tf-select").locator("option")).toHaveText(["None", ...DOMINO_OPTIONS.tfs]);
                await expect(page.getByTestId("river-ligands-select").locator("option")).toHaveText(DOMINO_OPTIONS.ligands);
                await expect(page.getByTestId("dot-sending-clusters-select").locator("option")).toHaveText(CLUSTERS);
                await expect(page.getByTestId("copy-share-link-button")).toBeEnabled();
            });

            it("should show the chosen output's dataset details below the output list", async () => {
                await chooseDominoFile(page);
                await expect(page.getByTestId("dataset-info-missing")).toBeHidden();

                await expect(page.getByTestId("dataset-info-citation")).toHaveText('"A <i>test</i> study of cell-cell communication", A Researcher, 2024');
                await expect(page.getByTestId("dataset-info-ids")).toHaveText("DOI: 10.1000/test.2024 | PMID: N/A");
                await expect(page.getByTestId("dataset-info-url")).toHaveText("gEAR Dataset: umgear.org");
                await expect(page.getByTestId("dataset-info-description")).toHaveText("Description: Test dataset");

                await expect(page.getByTestId("dataset-info-ids").getByRole("link", { name: "10.1000/test.2024" })).toHaveAttribute("href", "https://doi.org/10.1000/test.2024");
                await expect(page.getByTestId("dataset-info-url").getByRole("link")).toHaveAttribute("href", "https://umgear.org/");
                await expect(page.getByTestId("dataset-info-citation").locator("i")).toHaveCount(0);  // shown as text, not HTML
            });

            it("should show a disclaimer when an output's dataset details are missing or incomplete", async () => {
                // No .json file
                await chooseDominoFile(page, DOMINO_FILE_NO_INFO);
                await expect(page.getByTestId("dataset-info-missing")).toBeVisible();
                await expect(page.getByTestId("dataset-info-missing")).toHaveText("Dataset information has not been provided for this DominoSignal output.");
                await expect(page.getByTestId("dataset-info-details")).toBeHidden();

                // .json file without an author
                await chooseDominoFile(page, DOMINO_FILE_PARTIAL_INFO);
                await expect(page.getByTestId("dataset-info-missing")).toBeVisible();
                await expect(page.getByTestId("dataset-info-details")).toBeHidden();

                // Back to an output with complete details
                await chooseDominoFile(page, DOMINO_FILE);
                await expect(page.getByTestId("dataset-info-missing")).toBeHidden();
                await expect(page.getByTestId("dataset-info-details")).toBeVisible();
            });

            it("should collapse and expand a parameter section with its arrow", async () => {
                await chooseDominoFile(page);
                const toggle = page.getByTestId("dot-options-toggle");

                await toggle.click();
                await expect(page.getByTestId("dot-options-section")).toBeHidden();
                await expect(toggle).toHaveAttribute("aria-expanded", "false");

                await toggle.click();
                await expect(page.getByTestId("dot-options-section")).toBeVisible();
                await expect(toggle).toHaveAttribute("aria-expanded", "true");
            });

            it("should start with the advanced plotting options collapsed, and show the dot size as a percentage", async () => {
                await chooseDominoFile(page);
                await expect(page.getByTestId("advanced-options-section")).toBeHidden();

                await page.getByTestId("advanced-options-toggle").click();
                await expect(page.getByTestId("advanced-options-section")).toBeVisible();
                await expect(page.getByTestId("dot-size-value")).toHaveText("100%");

                await page.getByTestId("dot-size-slider").fill("2");
                await expect(page.getByTestId("dot-size-value")).toHaveText("200%");
            });

            it("should swap the dot plot sending and receiving clusters", async () => {
                await chooseDominoFile(page);
                await selectInDropdown(page, "dot-sending-clusters-select", ["Pericytes", "T cells"]);
                await selectInDropdown(page, "dot-receiving-clusters-select", ["Macrophages_Pos"]);

                await page.getByTestId("swap-dot-clusters-button").click();
                expect(await getSelectedValues(page, "dot-sending-clusters-select")).toEqual(["Macrophages_Pos"]);
                expect(await getSelectedValues(page, "dot-receiving-clusters-select")).toEqual(["Pericytes", "T cells"]);
                await expect(dropdown(page, "dot-receiving-clusters-select").locator(".current")).toHaveText(["Pericytes", "T cells"]);
            });

            it("should copy the dot plot clusters to the river plot options, replacing them", async () => {
                await chooseDominoFile(page);
                await selectInDropdown(page, "dot-sending-clusters-select", ["Pericytes"]);
                await selectInDropdown(page, "dot-receiving-clusters-select", ["Macrophages_Pos"]);
                await selectInDropdown(page, "river-sending-clusters-select", ["B cells"]);

                await page.getByTestId("copy-clusters-to-river-button").click();
                expect(await getSelectedValues(page, "river-sending-clusters-select")).toEqual(["Pericytes"]);
                expect(await getSelectedValues(page, "river-receiving-clusters-select")).toEqual(["Macrophages_Pos"]);
            });

            it("should remove a selected value with its chip's remove button", async () => {
                await chooseDominoFile(page);
                await selectInDropdown(page, "river-ligands-select", ["Gdf11", "Tgfb1"]);

                await dropdown(page, "river-ligands-select").getByRole("button", { name: "Remove Gdf11" }).click();
                expect(await getSelectedValues(page, "river-ligands-select")).toEqual(["Tgfb1"]);
                await expect(dropdown(page, "river-ligands-select")).not.toHaveClass(/open/);  // the dropdown stays closed
            });

            it("should show each requested plot in its own tab, in order", async () => {
                await chooseDominoFile(page);
                await selectInDropdown(page, "dot-sending-clusters-select", ["Pericytes"]);
                await choosePlots(page, ["dot_tf", "dot_clusters", "river_selected"]);
                await page.getByTestId("show-visualizations-button").click();

                const tabs = page.getByTestId("results-tabs").getByRole("tab");
                await expect(tabs).toHaveText(["Dot: Selected TF", "Dot: Selected clusters", "River: Custom selection"]);

                // No TF was selected, so that tab shows why, and the first tab with a plot opens instead
                await expect(page.getByTestId("plot-tab-dot_tf").locator(".plot-tab-warning")).toBeVisible();
                await expect(page.getByTestId("plot-tab-dot_clusters")).toHaveClass(/is-active/);
                await expect(page.getByTestId("plot-pane-dot_clusters").locator(".js-plotly-plot")).toBeVisible();

                await page.getByTestId("plot-tab-dot_tf").click();
                await expect(page.getByTestId("plot-pane-dot_tf").getByTestId("plot-card-message")).toHaveText("Select a TF to generate this plot");

                // The request carried the selections
                const request = plotRequests.at(-1);
                expect(JSON.parse(request.plots)).toEqual(["dot_tf", "dot_clusters", "river_selected"]);
                expect(JSON.parse(request.sending_clusters)).toEqual(["Pericytes"]);
            });

            it("should draw a tab's plot when it is first opened, and move between tabs with the arrow keys", async () => {
                await chooseDominoFile(page);
                await selectInDropdown(page, "dot-sending-clusters-select", ["Pericytes"]);
                await choosePlots(page, ["dot_clusters", "river_selected"]);
                await page.getByTestId("show-visualizations-button").click();

                await expect(page.getByTestId("plot-pane-dot_clusters").locator(".js-plotly-plot")).toBeVisible();
                await expect(page.getByTestId("plot-pane-river_selected").locator(".js-plotly-plot")).toHaveCount(0);  // not drawn yet

                await page.getByTestId("plot-tab-dot_clusters").getByRole("tab").focus();
                await page.keyboard.press("ArrowRight");
                await expect(page.getByTestId("plot-tab-river_selected")).toHaveClass(/is-active/);
                await expect(page.getByTestId("plot-pane-river_selected").locator(".js-plotly-plot")).toBeVisible();
                await expect(page.getByTestId("plot-pane-dot_clusters")).toBeHidden();
            });

            it("should send the thresholds and advanced plotting options with the request", async () => {
                await chooseDominoFile(page);
                await selectInDropdown(page, "dot-sending-clusters-select", ["Pericytes"]);
                await choosePlots(page, ["dot_clusters_filtered", "river_selected"]);
                await page.getByTestId("exp-range-threshold-input").fill("0.7");
                await page.getByTestId("river-apply-thresholds-checkbox").check();
                await page.getByTestId("advanced-options-toggle").click();
                await page.getByTestId("dot-size-slider").fill("2");
                await page.getByTestId("dot-white-background-checkbox").check();
                await page.getByTestId("show-visualizations-button").click();

                await expect(page.getByTestId("plot-tab-dot_clusters_filtered")).toBeVisible();
                const request = plotRequests.at(-1);
                expect(request.exp_range_threshold).toBe("0.7");
                expect(request.river_apply_thresholds).toBe("1");
                expect(request.dot_size_scale).toBe("2");
                expect(request.dot_white_background).toBe("1");
            });

            it("should restore selections and generate plots from a share link", async () => {
                const shareUrl = new URL(gearUrl);
                shareUrl.searchParams.append("domino_file", DOMINO_FILE);
                shareUrl.searchParams.append("selected_tf", "Foxp1");
                shareUrl.searchParams.append("sending_clusters", "Pericytes,T cells");
                shareUrl.searchParams.append("plots", "dot_tf,dot_clusters");
                shareUrl.searchParams.append("dot_size_scale", "1.5");
                await page.goto(shareUrl.toString());

                await expect(page.getByTestId("plot-tab-dot_tf")).toBeVisible();
                await expect(page.getByTestId("tf-select")).toHaveValue("Foxp1");
                expect(await getSelectedValues(page, "dot-sending-clusters-select")).toEqual(["Pericytes", "T cells"]);
                await expect(page.getByTestId("plot-toggle-dot_clusters")).toBeChecked();
                await expect(page.getByTestId("plot-toggle-river_selected")).not.toBeChecked();
                await expect(page.getByTestId("dot-size-value")).toHaveText("150%");

                // The URL is left as it was
                expect(page.url()).toBe(shareUrl.toString());
            });

            it("should download a plot as a PDF when chosen, and as a PNG otherwise", async () => {
                await chooseDominoFile(page);
                await selectInDropdown(page, "dot-sending-clusters-select", ["Pericytes"]);
                await choosePlots(page, ["dot_clusters"]);
                await page.getByTestId("show-visualizations-button").click();

                const cameraButton = page.getByTestId("plot-pane-dot_clusters").locator(".modebar-btn[data-title='Download plot']");
                await expect(cameraButton).toHaveCount(1);

                // PNG by default (drawn in the browser)
                let downloadPromise = page.waitForEvent("download");
                await cameraButton.click({ force: true });  // the toolbar only shows on hover
                expect((await downloadPromise).suggestedFilename()).toBe("domino_test_dot_clusters.png");

                // PDF from the server when the box is ticked
                await page.getByTestId("advanced-options-toggle").click();
                await page.getByTestId("download-pdf-checkbox").check();
                downloadPromise = page.waitForEvent("download");
                await cameraButton.click({ force: true });
                expect((await downloadPromise).suggestedFilename()).toBe("domino_test_dot_clusters.pdf");

                const request = plotRequests.at(-1);
                expect(request.image_format).toBe("pdf");
                expect(JSON.parse(request.plots)).toEqual(["dot_clusters"]);
            });
        });
    }
});
