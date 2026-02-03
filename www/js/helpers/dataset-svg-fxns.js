"use strict";

import * as d3 from "https://cdn.jsdelivr.net/npm/d3@7/+esm";
import { getCurrentUser } from "../common.v2.js";

// for those fields which have no reading, a specific value is sometimes put in instead
// These are colored a neutral color
const NA_FIELD_PLACEHOLDER = -0.012345679104328156;
const NA_FIELD_COLOR = '#808080';

/**
 * Color the SVG based on the chart data and plot configuration.
 * @param {Object} chartData - The data for the chart.
 * @param {string} datasetId - The ID of the dataset.
* @param {Object} plotConfig - The configuration for the plot.
 * @param {string} containerInfo - Various container information
 * @param {string} [svgScoringMethod="user_defined"] - The scoring method for coloring the SVG.
 * @returns {Promise<void>} - A promise that resolves when the SVG is colored.
 */
export const colorSVG = async (chartData, datasetId, plotConfig, containerInfo, svgScoringMethod="user_defined") => {
    // I found adding the mid color for the colorblind mode  skews the whole scheme towards the high color
    const colorblindMode = getCurrentUser()?.colorblind_mode || false;
    const colors = plotConfig?.colors || {};
    if (!colors) {
        throw new Error("No colors found in plotConfig.");
    }
    const { lowColor, midColor, highColor } = getConfigColors(colors, colorblindMode);

    // Fill in tissue classes with the expression colors
    const {data: expression} = chartData;

    const [scoringMethodUsed, score] = getScoreObject(chartData, svgScoringMethod);

    const imageContainer = containerInfo.imageContainer
    if (!imageContainer) {
        throw new Error("No imageContainer found in containerInfo to render SVG into.");
    }

    const containerId = containerInfo.containerId;
    if (!containerId) {
        throw new Error("No containerId found in containerInfo to render SVG into.");
    }

    // Load SVG file and set up the window
    const cardImage = document.querySelector(`${imageContainer}`);

    // create a legend div
    const legendDiv = document.createElement('div');
    legendDiv.classList.add('legend');
    legendDiv.style.zIndex = 1;
    legendDiv.style.height = "40px";    // match the viewbox height of child
    cardImage.append(legendDiv);

    // create a svg div (CSS for margin-top will now work nicely)
    const svgDiv = document.createElement('div');
    svgDiv.classList.add('svg');
    // higher z-index so we can mouseover the svg
    svgDiv.style.zIndex = 2;
    svgDiv.style.height = "calc(100% - 40px)";
    cardImage.append(svgDiv);

    const snap = Snap(svgDiv);
    const svg_path = `datasets_uploaded/${datasetId}.svg`;

    await Snap.load(svg_path, async (path) => {
        await snap.append(path);
        const svg = snap.select(`${imageContainer} svg`);
        svg.attr({
            width: "100%"
        });

        // TODO: Set viewbar just like the legend.
        // TODO: Set at bottom of card-image

        // Get all paths, circles, rects, and ellipses
        const paths = svg.selectAll(`path,circle,rect,ellipse`);

        // Rename path IDs to include the containerId
        paths.forEach(path => {
            path.attr('id', `${containerId}-${path.attr('id')}`);
        });

        const singleColorScaleMethods = ['gene', 'dataset', 'user_defined'];
        if (singleColorScaleMethods.includes(scoringMethodUsed)) {
            createSingleScaleColorMapping(paths, svgDiv, expression, score, lowColor, midColor, highColor);
            let title;
            switch (scoringMethodUsed) {
                case 'gene':
                    title = "Gene-level expression";
                    break;
                case 'dataset':
                    title = "Dataset-level expression";
                    break;
                case 'user_defined':
                    title = "Owner-defined expression";
                    break;
                default:
                    // Should never happen
                    title = "Expression";
            }
            // Draw the legend
            drawSVGLegend([lowColor, midColor, highColor], containerInfo, title, score);
            return;
        }
        if (scoringMethodUsed === 'tissue') {
            const geneSymbol = chartData.scores.gene.gene
            createTissueColorMapping(paths, svgDiv, expression, score, lowColor, midColor, highColor, containerInfo, geneSymbol);
        } else {
            throw new Error(`Invalid svg scoring method ${scoringMethodUsed}.`);
        }

    });
}

/**
 * Draws a legend for a SVG image
 *
 * @param {Array} colors - Array of colors for the legend (low, mid, high).
 * @param {Object} containerInfo - The container information object.
 * @param {string} title - The title for the legend.
 * @param {Object} score - The score object containing the minimum and maximum values.
 */
export const drawSVGLegend = (colors, containerInfo, title, score) => {
    const [ lowColor, midColor, highColor ] = colors;

    const containerId = containerInfo.containerId;

    const card = document.querySelector(`${containerInfo.outerContainer}`);
    const node = document.querySelector(`#${containerId} .legend`);

    const width = node.getBoundingClientRect().width;

    // Create our legend svg
    const legend = d3.select(node)  // returns document.documentElement
        .append('svg')
        .style('position', 'relative')
        .style('width', '100%')
        .style("height", "40px")    // Without a fixed height, the box is too tall and prevents mouseover of the svg image
        .attr('viewBox', `0 0 ${width} 40`)
        .attr('class', 'svg-gradient-container');
    const defs = legend.append('defs');
    // Define our gradient shape
    const linearGradient = defs
        .append('linearGradient')
        .attr('id', `${containerId}-linear-gradient`)
        .attr('x1', '0%')
        .attr('y1', '0%')
        .attr('x2', '100%')
        .attr('y2', '0%');

    const { min, max } = score;
    const range33 = ((max - min) / 3) + min;
    const range66 = (2 * (max - min) / 3) + min;

    // Create the gradient points for either three- or two-color gradients
    if (midColor) {
        // This means we've got a good distribution of min under 0 and max above
        //  it, so we can do a proper three-color range
        // midpoint offset calculation, so the mid color is at 0
        let midOffset = (Math.abs(min) / (max + Math.abs(min))) * 100;
        if (min >= 0 || max <= 0) {
            midOffset = 50  // not diverging, just put in the middle
        }

        linearGradient
            .append('stop')
            .attr('offset', '0%')
            .attr('stop-color', lowColor);
        linearGradient
            .append('stop')
            .attr('offset', `${midOffset}%`)
            .attr('stop-color', midColor);
        linearGradient
            .append('stop')
            .attr('offset', '100%')
            .attr('stop-color', highColor);
    } else {
        linearGradient
            .append('stop')
            .attr('offset', '0%')
            .attr('stop-color', lowColor);
        linearGradient
            .append('stop')
            .attr('offset', '100%')
            .attr('stop-color', highColor);
    }

    // Draw the rectangle using the linear gradient
    legend
        .append('rect')
        .attr('width', "50%")
        .attr('y', 15)
        .attr('x', "25%")
        .attr('height', 10) // quarter of viewport height
        .style(
            'fill',
            `url(#${containerId}-linear-gradient)`
        );

    // Define the x-axis range
    const xScale = d3
        .scaleLinear()
        .domain([min, max])
        .range([0, width / 2]).clamp(true);

    const xAxis = d3
        .axisBottom(xScale)
        .tickValues([min, range33, range66, max])

    // Add the x-axis to the legend
    legend
        .append('g')
        .attr('class', 'axis')
        .attr('transform', `translate(${width / 4}, 22)`)   // start quarter from left, and 10 px below rectangle
        .attr("stroke", "black")
        .call(xAxis);

    // Make tick marks black
    legend.selectAll('.tick line')
        .attr('stroke', 'black');

    // Hide upper tick bar
    legend.selectAll(".domain ")
        .attr("opacity", 0);

    // Add title
    legend
        .append('text')
        .attr('x', "50%")
        .attr('y', 12)
        .attr('text-anchor', 'middle')
        .attr('font-size', '10px')
        .attr('font-weight', 'bold')
        .text(title);

    // Ensure axis is responsive
    window.addEventListener('resize', () => {
        legend.attr('viewbox', `0 0 ${card.getBoundingClientRect().width} 40`);

        // purge old axis
        legend.select('.axis').remove();

        // redraw axis
        // TODO: this is hacky, but it works
        const xScale = d3
            .scaleLinear()
            .domain([min, max])
            .range([0, card.getBoundingClientRect().width / 2]).clamp(true);

        const xAxis = d3
            .axisBottom(xScale)
            .tickValues([min, range33, range66, max])

        // Add the x-axis to the legend
        legend
            .append('g')
            .attr('class', 'axis')
            .attr('transform', `translate(${card.getBoundingClientRect().width / 4}, 22)`)   // start quarter from left, and 22 px below rectangle
            .attr("stroke", "black")
            .call(xAxis);

        // Make tick marks black
        legend.selectAll('.tick line')
            .attr('stroke', 'black');

        // Hide upper tick bar
        legend.selectAll(".domain ")
            .attr("opacity", 0);
    });
}


/**
 * Retrieves the score object for a given scoring method from the chart data.
 *
 * If the requested scoring method is "user_defined" but it does not exist in `chartData.scores`,
 * the function falls back to the "gene" scoring method and logs a warning.
 *
 * @param {Object} chartData - The chart data containing scores.
 * @param {Object} chartData.scores - An object mapping scoring methods to their respective scores.
 * @param {string} svgScoringMethod - The scoring method to retrieve ("user_defined", "gene", etc.).
 * @returns {[string, any]} A tuple where the first element is the scoring method used,
 *   and the second element is the corresponding score object.
 */
const getScoreObject = (chartData, svgScoringMethod) => {
    if (svgScoringMethod === "user_defined" && !chartData.scores.hasOwnProperty("user_defined")) {
        console.warn("user_defined scoring method requested but not found in chartData.scores. Falling back to gene scoring method.");
        return [ "gene", chartData.scores.gene ];
    }
    return [svgScoringMethod, chartData.scores[svgScoringMethod]];
}

const getConfigColors = (colors, colorblindMode) => {
    const lowColor = colorblindMode ? 'rgb(254, 232, 56)' : colors.low_color;
    const midColor = colorblindMode ? null : colors.mid_color
    const highColor = colorblindMode ? 'rgb(0, 34, 78)' : colors.high_color;
    return {lowColor, midColor, highColor};
}

/**
 * Renders SVG paths with a single color scale based on expression values.
 * Supports two-color and three-color gradients, and handles NA values.
 * Attaches tooltips to each path showing the tissue name and formatted expression score.
 *
 * @param {Array} paths - Array of SnapSVG path objects to color.
 * @param {Object} svgDiv - The SVG container div for appending tooltips.
 * @param {Object} expression - Object mapping tissue names to expression values.
 * @param {Object} score - Object containing min and max values for the color scale.
 * @param {string} lowColor - Color for the low end of the scale.
 * @param {string} midColor - Color for the middle of the scale (optional for three-color gradient).
 * @param {string} highColor - Color for the high end of the scale.
 */
const createSingleScaleColorMapping = (paths, svgDiv, expression, score, lowColor, midColor, highColor) => {
    const { min, max } = score;
    let color = null;

    // are we doing a three- or two-color gradient?
    if (midColor) {
        let mid = 0;    // Diverging color scale
        if (min >= 0 || max <= 0) {
            // If not diverging, just set the mid color halfway
            mid = (min + max) / 2;
        }
        // We have a good value range, do the three-color
        color = d3
            .scaleLinear()
            .domain([min, mid, max])
            .range([lowColor, midColor, highColor])
            .clamp(true);   // prevent colors outside the range
    } else {
        color = d3
            .scaleLinear()
            .domain([min, max])
            .range([lowColor, highColor])
            .clamp(true);
    }

    // NOTE: This must use the SnapSVG API Set.forEach function to iterate
    paths.forEach(path => {
        const tissue_classes = path.node.className.baseVal.split(' ');
        tissue_classes.forEach(tissue => {
            if (!expression.hasOwnProperty(tissue)) {
                return;
            }

            if (expression[tissue] == NA_FIELD_PLACEHOLDER) {
                path.attr('fill', NA_FIELD_COLOR);
            } else {
                path.attr('fill', color(expression[tissue]));
            }

            // log-transfom the expression score
            const math = "raw";
            let score;
            // Apply math transformation to expression score
            if (math == 'log2') {
                score = d3.format('.2f')(Math.log2(expression[tissue]));
            } else if (math == 'log10') {
                score = d3.format('.2f')(Math.log10(expression[tissue]));
            } else {
                //math == 'raw'
                score = d3.format('.2f')(expression[tissue]);
            }

            const tooltip = document.createElement('div');
            tooltip.classList.add('tooltip');
            tooltip.style.position = 'relative';
            tooltip.style.fontSize = "12px";
            tooltip.style.bottom = `${0}px`;
            tooltip.style.left = `${0}px`;
            tooltip.style.backgroundColor = 'white';
            tooltip.style.opacity = 0.8;
            tooltip.style.color = 'black';
            tooltip.style.padding = '5px';
            tooltip.style.border = '1px solid gray';
            tooltip.style.zIndex = 3;

            // Place tissue in score in a nice compact tooltip
            const tooltipText = `<strong>${tissue}</strong>: ${score}`;
            tooltip.innerHTML = tooltipText;

            // Add mouseover and mouseout events to create and destroy the tooltip
            path.mouseover(() => {
                // Add tooltip to the bottom-left of the SVG
                svgDiv.appendChild(tooltip);

            });
            path.mouseout(() => {
                svgDiv.querySelector('.tooltip').remove();
            });

        });
    });
}

/**
 * Renders a color scale for tissue expression values on SVG paths, adds tooltips and legends for each tissue.
 *
 * @param {Array} paths - Array of D3 selection SVG path elements representing tissues.
 * @param {Object} svgDiv - The SVG container div for appending tooltips.
 * @param {Object} expression - Object mapping tissue names to expression values for the selected gene.
 * @param {Object} score - Object mapping tissue names to min/max expression values for color scaling.
 * @param {string} lowColor - Color representing the lowest expression value.
 * @param {string} midColor - (Optional) Color representing the mid expression value (used for diverging scales).
 * @param {string} highColor - Color representing the highest expression value.
 * @param {Object} containerInfo - Information about the SVG container.
 * @param {string} geneSymbol - Gene symbol being visualized.
 */
const createTissueColorMapping = (paths, svgDiv, expression, score, lowColor, midColor, highColor, containerInfo, geneSymbol) => {
    // tissue scope
    const tissues = Object.keys(score);
    const color = {};

    if (midColor) {

        tissues.forEach(tissue => {
            let {
                min,
                max
            } = score[tissue];

            let mid = 0;    // Diverging color scale
            if (min >= 0 || max <= 0) {
                // If not diverging, just set the mid color halfway
                mid = (min + max) / 2;
            }

            // We have a good value range, do the three-color
            color[tissue] = d3
                .scaleLinear()
                .domain([min, mid, max])
                .range([lowColor, midColor, highColor])
                .clamp(true);   // prevent colors outside the range
        });
    } else {
        tissues.forEach(tissue => {
            let {
                min,
                max
            } = score[tissue];

            color[tissue] = d3
                .scaleLinear()
                .domain([min, max])
                .range([lowColor, highColor])
                .clamp(true);
        });
    }
    // Add instructions to hover over a tissue class to see the legend
    const instructions = document.createElement('div');
    instructions.textContent = `Hover over a tissue to see ${geneSymbol} expression compared to all genes only for this tissue`;
    instructions.style.position = 'relative';
    instructions.style.top = `${0}px`;
    instructions.style.left = `${0}px`;
    instructions.style.backgroundColor = 'white';
    instructions.style.color = 'black';
    instructions.style.fontWeight = 'bold';
    instructions.style.padding = '5px';
    instructions.style.zIndex = 3;
    instructions.style.alignContent = 'center';

    const legendDiv = document.querySelector(`#${containerInfo.containerId} .legend`);
    legendDiv.appendChild(instructions);

    paths.forEach(path => {
        const tissue_classes = path.node.className.baseVal.split(' ');
        tissue_classes.forEach(tissue => {
            if (!(tissue && color[tissue])) {
                return;
            }
            const color_scale = color[tissue];
            if (expression[tissue] == NA_FIELD_PLACEHOLDER) {
                path.attr('fill', NA_FIELD_COLOR);
            } else {
                path.attr('fill', color_scale(expression[tissue]));
            }

            // log-transfom the expression score
            const math = "raw";
            let expressionScore;
            // Apply math transformation to expression score
            if (math == 'log2') {
                expressionScore = d3.format('.2f')(Math.log2(expression[tissue]));
            } else if (math == 'log10') {
                expressionScore = d3.format('.2f')(Math.log10(expression[tissue]));
            } else {
                //math == 'raw'
                expressionScore = d3.format('.2f')(expression[tissue]);
            }

            const tooltip = document.createElement('div');
            tooltip.classList.add('tooltip');
            tooltip.style.position = 'relative';
            tooltip.style.fontSize = "12px";
            tooltip.style.bottom = `${0}px`;
            tooltip.style.left = `${0}px`;
            tooltip.style.backgroundColor = 'white';
            tooltip.style.opacity = 0.8;
            tooltip.style.color = 'black';
            tooltip.style.padding = '5px';
            tooltip.style.border = '1px solid gray';
            tooltip.style.zIndex = 3;

            // Place tissue in score in a nice compact tooltip
            const tooltipText = `<strong>${tissue}</strong>: ${expressionScore}`;
            tooltip.innerHTML = tooltipText;

            // Add mouseover and mouseout events to create and destroy the tooltip
            path.mouseover(() => {
                // clear legend
                legendDiv.replaceChildren();

                // Add tooltip to the bottom-left of the SVG
                svgDiv.appendChild(tooltip);

                const tissueScore = {min: score[tissue].min, max: score[tissue].max};

                const title = "Tissue-level expression";
                // draw legend for this tissue class
                drawSVGLegend([lowColor, midColor, highColor], containerInfo, title, tissueScore);
            });
            path.mouseout(() => {
                svgDiv.querySelector('.tooltip').remove();
                // clear legend
                legendDiv.replaceChildren();

                legendDiv.appendChild(instructions);
            });

        });
    });
}