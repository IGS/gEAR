document.addEventListener('DOMContentLoaded', () => {
    // Add listeners for any elements which have generic close controls
    (document.querySelectorAll('.notification .delete') || []).forEach(($delete) => {
        const $notification = $delete.parentNode;

        $delete.addEventListener('click', () => {
            $notification.parentNode.removeChild($notification);
        });
    });
});

// If "step" sections are clicked, expand that section and collaps the others
const jsSteps = document.getElementsByClassName("js-steps");

const resetSteps = (event) => {
    // Make every step uniform
    for (const jsStep of jsSteps) {
        jsStep.classList.remove("has-background-light")
        jsStep.classList.add("has-background-white-bis");

        const rightIcons = jsStep.querySelectorAll("span.is-pulled-right i");
        rightIcons.forEach((elt) => {
            elt.classList.remove("mdi-chevron-up");
            elt.classList.add("mdi-chevron-down");
        });

        //TODO: Collapse all
    }

    // Reset this section
    const currentStep = event.target.closest(".js-steps");
    currentStep.classList.remove("has-background-white-bis");
    currentStep.classList.add("has-background-light");
    currentStep.querySelectorAll("span.is-pulled-right i")[0].classList.remove("mdi-chevron-down");
    currentStep.querySelectorAll("span.is-pulled-right i")[0].classList.add("mdi-chevron-up");
    //TODO: Expand all
}

for (const jsStep of jsSteps) {
    jsStep.addEventListener("click", resetSteps);
}