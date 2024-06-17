'use strict';

// Stepper Functions
// Functions ro reduce redundancy in common actions for stepper components
// This is for https://github.com/octoshrimpy/bulma-o-steps/tree/master

// States
// Completed steps
// parent.is-dashed: The past steps that have been completed
// i.mdi-check: The icon for the completed steps

// Active step
// parent.is-active: Line marker up to the current step
// is-light: The current step
// i.mdi-pencil: The icon for the active step

// Failed step
// is-danger: Failed step
// parent.is-active: Line marker up to the failed step


/**
 * Marks a step as failed by adding appropriate CSS classes to the step element.
 * @param {string} selectorHref - The href attribute value of the step element to mark as failed.
 */
const failStepWithHref = (selectorHref) => {
    document.querySelector(`a[href='${selectorHref}']`).parentElement.classList.add("is-active")
    document.querySelector(`a[href='${selectorHref}']`).classList.remove("is-light")
    document.querySelector(`a[href='${selectorHref}']`).classList.add("is-danger")
    document.querySelector(`a[href='${selectorHref}'] i`).classList.remove("mdi-check")
}


/**
 * Marks a step as failed by adding appropriate CSS classes.
 * @param {string} selector - The CSS selector for the step element.
 */
const failStep = (selector) => {
    document.querySelector(selector).parentElement.classList.add("is-active")
    document.querySelector(selector).classList.remove("is-light")
    document.querySelector(selector).classList.add("is-danger")
    document.querySelector(`a[href='${selector}'] i`).classList.remove("mdi-check")
}

/**
 * Marks a step as passed by adding appropriate CSS classes to the step element.
 * @param {string} selectorHref - The href attribute value of the step element to mark as passed.
 */
const passStepWithHref = (selectorHref) => {
    document.querySelector(`a[href='${selectorHref}']`).parentElement.classList.remove("is-active")
    document.querySelector(`a[href='${selectorHref}']`).parentElement.classList.add("is-dashed")
    document.querySelector(`a[href='${selectorHref}']`).classList.remove("is-danger", "is-light")
    document.querySelector(`a[href='${selectorHref}'] i`).classList.remove("mdi-pencil")
    document.querySelector(`a[href='${selectorHref}'] i`).classList.add("mdi-check")
}

/**
 * Marks a step as passed by adding appropriate CSS classes.
 * @param {string} selector - The CSS selector for the step element.
 */
const passStep = (selector) => {
    document.querySelector(selector).parentElement.classList.remove("is-active")
    document.querySelector(selector).parentElement.classList.add("is-dashed")
    document.querySelector(selector).classList.remove("is-danger", "is-light")
    document.querySelector(`${selector} i`).classList.remove("mdi-pencil")
    document.querySelector(`${selector} i`).classList.add("mdi-check")
}

/**
 * Opens the next step hrefs and applies appropriate classes to the elements.
 * @param {string[]} selectorHrefs - An array of hrefs for the selectors.
 * @param {string|null} activeSelectorHref - The href of the active selector. Defaults to null.
 * @param {boolean} clickActive - Specifies whether to click the active step. Defaults to false.
 */
const openNextStepHrefs = (selectorHrefs, activeSelectorHref=null, clickActive=false) => {

    // if length of selectorHrefs is 1, then the step is the active step
    if (selectorHrefs.length === 1) {
        activeSelectorHref = selectorHrefs[0]
    }

    for (const href of selectorHrefs) {
        document.querySelector(`a[href='${href}']`).classList.add("is-light")
        document.querySelector(`a[href='${href}'] i`).classList.add("mdi-pencil")
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
 * @param {string[]} selectors - An array of CSS selectors for the steps.
 * @param {string|null} activeSelector - The CSS selector for the active step. If null, the first step will be considered active.
 */
const openNextSteps = (selectors, activeSelector=null) => {

    // if length of selectors is 1, then the step is the active step
    if (selectors.length === 1) {
        activeSelector = selectors[0]
    }

    for (const selector of selectors) {
        document.querySelector(selector).classList.add("is-light");
        document.querySelector(`${selector} i`).classList.add("mdi-pencil")
    }

    if (activeSelector) {
        document.querySelector(activeSelector).parentElement.classList.add("is-active")
    }
}