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

// Irreversible (blocked) steps
// is-dark: The blocked step
// i.mdi-cancel: The icon for the irreversible steps
// parent.is-dashed turned off

// I add `.steps:not(.is-hidden)` to the selectors to avoid selecting hidden steps
// in the single-cell workbench, there are two steppers that have the same selectors, but one is hidden

// NOTE: These functions make the assumption that the active selector or href is in the same element as the .steps-marker


/**
 * Blocks a step with the specified href.
 * @param {string} selectorHref - The href of the step to be blocked.
 */
const blockStepWithHref = (selectorHref) => {
    document.querySelector(`.steps:not(.is-hidden) a[href='${selectorHref}']`).parentElement.classList.remove("is-dashed", "is-active");
    document.querySelector(`.steps:not(.is-hidden) a[href='${selectorHref}']`).classList.add("is-dark");
    document.querySelector(`.steps:not(.is-hidden) a[href='${selectorHref}'] i`).classList.remove("mdi-check");
    document.querySelector(`.steps:not(.is-hidden) a[href='${selectorHref}'] i`).classList.add("mdi-cancel");
}

/**
 * Blocks a step by modifying the CSS classes of the specified selector.
 * @param {string} selector - The CSS selector of the step element.
 */
const blockStep = (selector) => {
    document.querySelector(`.steps:not(.is-hidden) ${selector}`).parentElement.classList.remove("is-dashed", "is-active");
    document.querySelector(`.steps:not(.is-hidden) ${selector}`).classList.add("is-dark");
    document.querySelector(`.steps:not(.is-hidden) ${selector} i`).classList.remove("mdi-check");
    document.querySelector(`.steps:not(.is-hidden) ${selector} i`).classList.add("mdi-cancel");
}

/**
 * Marks a step as failed by adding appropriate CSS classes to the step element.
 * @param {string} selectorHref - The href attribute value of the step element to mark as failed.
 */
const failStepWithHref = (selectorHref) => {
    document.querySelector(`.steps:not(.is-hidden) a[href='${selectorHref}']`).parentElement.classList.add("is-active")
    document.querySelector(`.steps:not(.is-hidden) a[href='${selectorHref}']`).classList.remove("is-light")
    document.querySelector(`.steps:not(.is-hidden) a[href='${selectorHref}']`).classList.add("is-danger")
    document.querySelector(`.steps:not(.is-hidden) a[href='${selectorHref}'] i`).classList.remove("mdi-check")
}


/**
 * Marks a step as failed by adding appropriate CSS classes.
 * @param {string} selector - The CSS selector for the step element.
 */
const failStep = (selector) => {
    document.querySelector(`.steps:not(.is-hidden) ${selector}`).parentElement.classList.add("is-active")
    document.querySelector(`.steps:not(.is-hidden) ${selector}`).classList.remove("is-light")
    document.querySelector(`.steps:not(.is-hidden) ${selector}`).classList.add("is-danger")
    document.querySelector(`.steps:not(.is-hidden) ${selector} i`).classList.remove("mdi-check")
}

/**
 * Marks a step as passed by adding appropriate CSS classes to the step element.
 * @param {string} selectorHref - The href attribute value of the step element to mark as passed.
 */
const passStepWithHref = (selectorHref) => {
    document.querySelector(`.steps:not(.is-hidden) a[href='${selectorHref}']`).parentElement.classList.remove("is-active")
    document.querySelector(`.steps:not(.is-hidden) a[href='${selectorHref}']`).parentElement.classList.add("is-dashed")
    document.querySelector(`.steps:not(.is-hidden) a[href='${selectorHref}']`).classList.remove("is-danger", "is-light")
    document.querySelector(`.steps:not(.is-hidden) a[href='${selectorHref}'] i`).classList.remove("mdi-pencil")
    document.querySelector(`.steps:not(.is-hidden) a[href='${selectorHref}'] i`).classList.add("mdi-check")
}

/**
 * Marks a step as passed by adding appropriate CSS classes.
 * @param {string} selector - The CSS selector for the step element.
 */
const passStep = (selector) => {
    document.querySelector(`.steps:not(.is-hidden) ${selector}`).parentElement.classList.remove("is-active")
    document.querySelector(`.steps:not(.is-hidden) ${selector}`).parentElement.classList.add("is-dashed")
    document.querySelector(`.steps:not(.is-hidden) ${selector}`).classList.remove("is-danger", "is-light")
    document.querySelector(`.steps:not(.is-hidden) ${selector} i`).classList.remove("mdi-pencil")
    document.querySelector(`.steps:not(.is-hidden) ${selector} i`).classList.add("mdi-check")
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
        document.querySelector(`.steps:not(.is-hidden) a[href='${href}']`).classList.add("is-light")
        document.querySelector(`.steps:not(.is-hidden) a[href='${href}'] i`).classList.add("mdi-pencil")
    }

    if (activeSelectorHref) {
        document.querySelector(`.steps:not(.is-hidden) a[href='${activeSelectorHref}']`).parentElement.classList.add("is-active")
        // click the active step
        if (clickActive) {
            document.querySelector(`.steps:not(.is-hidden) a[href='${activeSelectorHref}']`).click()
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
        document.querySelector(`.steps:not(.is-hidden) ${selector}`).classList.add("is-light");
        document.querySelector(`.steps:not(.is-hidden) ${selector} i`).classList.add("mdi-pencil")
    }

    if (activeSelector) {
        document.querySelector(`.steps:not(.is-hidden) ${activeSelector}`).parentElement.classList.add("is-active")
    }
}

/**
 * Resets the stepper by removing classes and icons from each step element.
 * Optionally, sets a specific step as active based on the provided href.
 *
 * @param {string|null} activeSelectorHref - The href of the step to set as active. Default is null.
 */
const resetStepperWithHrefs = (activeSelectorHref=null) => {
    const steps = document.querySelectorAll(".steps:not(.is-hidden) .steps-marker")
    for (const step of steps) {
        step.classList.remove("is-light", "is-danger");
        step.parentElement.classList.remove("is-active", "is-dashed");
        step.querySelector("i").classList.remove("mdi-pencil", "mdi-check");
    }


    if (activeSelectorHref) {
        // pass all steps before the active step
        const selectorHrefs = Array.from(steps).map(step => step.getAttribute("href"));
        const activeIndex = selectorHrefs.indexOf(activeSelectorHref)

        if (activeIndex > 0) {
            // if found, pass all steps before the active step
            for (let i = 0; i < activeIndex; i++) {
                passStepWithHref(selectorHrefs[i])
            }
            openNextStepHrefs([activeSelectorHref])
        }
    } else {
        // make the first step active
        document.querySelector(".steps:not(.is-hidden) .steps-marker").parentElement.classList.add("is-active")
    }
}

/**
 * Resets the stepper by removing classes and icons from each step element.
 * Optionally, sets a specific step as active.
 *
 * @param {string|null} activeSelector - The selector for the step to set as active. If null, no step will be set as active.
 */
const resetStepper = (activeSelector=null) => {
    const steps = document.querySelectorAll(".steps:not(.is-hidden) .steps-marker")
    for (const step of steps) {
        step.classList.remove("is-light", "is-danger", "is-dark");
        step.parentElement.classList.remove("is-active", "is-dashed");
        step.querySelector("i").classList.remove("mdi-pencil", "mdi-check", "mdi-cancel");
    }

    if (activeSelector) {
        // TODO: need to test this
        // Select all steps before the active step, and pass them
        const beforeSteps = document.querySelectorAll(`.steps:not(.is-hidden) .steps-marker ~ ${activeSelector}`);
        for (const step of beforeSteps) {
            passStep(step);
        }
        openNextSteps([activeSelector]);

    } else {
        // make the first step active
        document.querySelector(".steps:not(.is-hidden) .steps-marker").parentElement.classList.add("is-active")
    }

}
