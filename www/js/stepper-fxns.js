'use strict';

// Stepper Functions
// Functions ro reduce redundancy in common actions for stepper components
// This is for https://github.com/octoshrimpy/bulma-o-steps/tree/master

/**
 * Fails a step with the specified href by adding the "is-active" and "is-danger" classes to the corresponding anchor element.
 * @param {string} selectorHref - The href value of the anchor element to be failed.
 */
const failStepWithHref = (selectorHref) => {
    document.querySelector(`a[href='${selectorHref}']`).parentElement.classList.add("is-active")
    document.querySelector(`a[href='${selectorHref}']`).classList.add("is-danger")
}

/**
 * Adds the "is-active" and "is-danger" classes to the element with the specified selector.
 * @param {string} selector - The CSS selector of the element to add the classes to.
 */
const failStep = (selector) => {
    document.querySelector(selector).parentElement.classList.add("is-active")
    document.querySelector(selector).classList.add("is-danger")
}

/**
 * Passes a step by removing the "is-active" and "is-dashed" classes from the anchor element with the specified href,
 * and adds the "is-hollow" class to it.
 *
 * @param {string} selectorHref - The href value of the anchor element to modify.
 */
const passStepWithHref = (selectorHref) => {
    document.querySelector(`a[href='${selectorHref}']`).parentElement.classList.remove("is-active", "is-dashed")
    document.querySelector(`a[href='${selectorHref}']`).classList.remove("is-danger")
    document.querySelector(`a[href='${selectorHref}']`).classList.add("is-hollow")
}

/**
 * Removes the "is-active" and "is-dashed" classes from the element with the specified selector,
 * and adds the "is-hollow" class.
 *
 * @param {string} selector - The CSS selector of the element to modify.
 */
const passStep = (selector) => {
    document.querySelector(selector).parentElement.classList.remove("is-active", "is-dashed")
    document.querySelector(selector).classList.remove("is-danger")
    document.querySelector(selector).classList.add("is-hollow")

}

/**
 * Adds CSS classes to the specified selector hrefs to control their appearance.
 * @param {string[]} selectorHrefs - An array of selector hrefs.
 * @param {string|null} activeSelectorHref - The active selector href.
 */
const openNextStepHrefs = (selectorHrefs, activeSelectorHref=null, clickActive=false) => {

    // if length of selectorHrefs is 1, then the step is the active step
    if (selectorHrefs.length === 1) {
        activeSelectorHref = selectorHrefs[0]
    }

    for (const href of selectorHrefs) {
        document.querySelector(`a[href='${href}']`).parentElement.classList.add("is-dashed")
        if (href !== activeSelectorHref) {
            document.querySelector(`a[href='${href}']`).classList.add("is-hollow")
        }
    }

    if (activeSelectorHref) {
        document.querySelector(`a[href='${activeSelectorHref}']`).parentElement.classList.add("is-active")
        // click the active step
        if (clickActive) {
            document.querySelector(`a[href='${activeSelectorHref}']`).click()
        }
    }

}

/**
 * Opens the next steps in a stepper.
 *
 * @param {string[]} selectors - An array of CSS selectors for the steps to open.
 * @param {string|null} activeSelector - The CSS selector for the currently active step.
 */
const openNextSteps = (selectors, activeSelector=null) => {

    // if length of selectors is 1, then the step is the active step
    if (selectors.length === 1) {
        activeSelector = selectors[0]
    }

    for (const selector of selectors) {
        document.querySelector(selector).parentElement.classList.add("is-dashed");

        if (selector !== activeSelector) {
            document.querySelector(selector).classList.add("is-hollow")
        }
    }

    if (activeSelector) {
        document.querySelector(activeSelector).parentElement.classList.add("is-active")
    }
}