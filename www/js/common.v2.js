document.addEventListener('DOMContentLoaded', () => {
    // Add listeners for any elements which have generic close controls
    // NOTE: ".delete" has Bulma CSS associated, which overwrites mdi icon classes
    (document.querySelectorAll('.notification .delete') || []).forEach(($delete) => {
        // const $notification = $delete.parentNode;
        // SAdkins - modifying to find ".notification" in case .delete is not immediate child. Otherwise partial removal could happen.
        const $notification = $delete.closest(".notification");

        $delete.addEventListener('click', () => {
            $notification.parentNode.removeChild($notification);
        });
    });
});

// If "step" sections are clicked, expand that section and collaps the others
const jsSteps = document.getElementsByClassName("js-step");
const resetSteps = (event) => {
    // Make every step uniform
    for (const jsStep of jsSteps) {
        // Reset the bg color
        jsStep.classList.remove("has-background-light")
        jsStep.classList.add("has-background-white-bis");

        // Reset the rightmost arrow
        const rightIcon = jsStep.querySelector("span.is-pulled-right i");
        rightIcon.classList.remove("mdi-chevron-up");
        rightIcon.classList.add("mdi-chevron-down");

        // Collapse all "non-header" parts of step
        const collapsableElt = jsStep.querySelector(".js-step-collapsable");
        collapsableElt.style.display = "none";
    }

    // Modify this step to look active
    const currentStep = event.target.closest(".js-step");
    currentStep.classList.remove("has-background-white-bis");
    currentStep.classList.add("has-background-light");
    currentStep.querySelector("span.is-pulled-right i").classList.remove("mdi-chevron-down");
    currentStep.querySelector("span.is-pulled-right i").classList.add("mdi-chevron-up");
    //TODO: Expand all
    const collapsableElt = currentStep.querySelector(".js-step-collapsable");
    collapsableElt.style.display = "block";

}

for (const jsStep of jsSteps) {
    jsStep.addEventListener("click", resetSteps);
}