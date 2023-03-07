/*
Pagination
code inspired by - https://webdesign.tutsplus.com/tutorials/pagination-with-vanilla-javascript--cms-41896
TODO: Consider doing this server-side to reduce load if datasets are huge
*/

const paginationNumbers = $("#pagination-numbers");
const nextButton = $("#results-next-button");
const prevButton = $("#results-prev-button");
const resultsItems = $("#results tr");

const paginationLimit = 20;
let pageCount = Math.ceil(resultsItems.length / paginationLimit);    // Number of pages to render
let currentPage = 1;
let datasetsSelected = 0;

const disableButton = (button) => button.css("visibility", "hidden");   // Use "visibility" to preserve space in layout
const enableButton = (button) => button.css("visibility", "visible");

const closeModal = (modal) => modal.removeClass("is-active");
const openModal = (modal) => modal.addClass("is-active");

// Situationally enable/disable Prev/Next buttons
const handlePageButtonsStatus = () => {
    if (currentPage === 1) {
        disableButton(prevButton);
    } else {
        enableButton(prevButton);
    }

    if (pageCount === currentPage) {
        disableButton(nextButton);
    } else {
        enableButton(nextButton);
    }
};

// Update which page is current
const handleActivePageNumber = () => {
    $(".pagination-link").each( (i, button) => {
        const buttonItem = $(button);
        const buttonText = buttonItem.text();
        buttonItem.removeClass("is-current");
        buttonItem.removeAttr("aria-current");
        buttonItem.attr("aria-label", `Goto page ${buttonText}`);
        const pageIndex = Number(buttonItem.attr("page-index"));
        if (pageIndex == currentPage) {
            buttonItem.addClass("is-current");
            buttonItem.attr("aria-current", "page");
            buttonItem.attr("aria-label", `Page ${buttonText}`);
        }
    });
};

// Add pagination page numbers to the container
const appendPageNumber = (index) => {
    const pageNumber = $('<li></li>');
    pageNumber.html(`<a class="pagination-link" page-index="${index}">${index}</a>`);
    pageNumber.attr("aria-label", "Goto page " + index);
    paginationNumbers.append(pageNumber);
};

  // Update the pagination container
const getPaginationNumbers = () => {
    for (let i = 1; i <= pageCount; i++) {
        appendPageNumber(i);
    }
};

// Set the active page
const setCurrentPage = (pageNum) => {
    currentPage = pageNum;

    handleActivePageNumber();
    handlePageButtonsStatus();

    const prevRange = (pageNum - 1) * paginationLimit;
    const currRange = pageNum * paginationLimit;

    // Restrict results to only N items from the total results
    resultsItems.each( (i, el) => {
        $(el).hide();
        if (i >= prevRange && i < currRange) {
            $(el).show();
        }
    });
};

// Reset pagination page rendering
const resetPagination = () => {
    pageCount = Math.ceil(resultsItems.length / paginationLimit);    // Number of pages to render
    getPaginationNumbers();
    setCurrentPage(1);

    prevButton.click( () => {
        setCurrentPage(currentPage - 1);
    });

    nextButton.click( () => {
        setCurrentPage(currentPage + 1);
    });

    // Make "page number" buttons interactive
    $(".pagination-link").each((i, button) => {
        const pageIndex = Number($(button).attr("page-index"));
        if (pageIndex) {
            $(button).click( () => {
                setCurrentPage(pageIndex);
            });
        }
    });

    $("#pagination-limit").text(paginationLimit > resultsItems.length ? resultsItems.length : paginationLimit);
    $("#results-total").text(resultsItems.length);
}

/* Dropdown panel */
$(".dropdown-trigger div").click(function() {
    const dropdownElt = $(this).closest(".dropdown");
    // If this element is already active, close it.
    if ($(dropdownElt).hasClass("is-active")) {
        $(dropdownElt).removeClass("is-active");
        return
    }
    // Set is-active to the newly clicked dropdown
    $(".dropdown").each((i, el) => $(el).removeClass("is-active"));
    $(dropdownElt).addClass("is-active");
});

/* Axios fetch commands */

const fetchNemoArchiveFiles = () => {

    const data = {
        filters: {"op":"and","content":[{"op":"in","content":{"field":"file.subtype","value":["counts"]}}]}
        , from: 1
        , size: 20
        , sort: "file_id:asc"
    }

    axios({
        method: 'get',
        url: 'https://portal.nemoarchive.org/api/gql/_mapping',
        //url: 'https://portal.nemoarchive.org/api/files',
        //data: data
        }).then(function (response) {
            console.log(response.data)
        });
}

$('#results .selection').click(function() {
    this.checked ? datasetsSelected++ : datasetsSelected--;  // not jQuery $(this)
    $("#selected-datasets").text(datasetsSelected);
    $("#dataset-selection-submit").prop("disabled", false);

    if (!datasetsSelected) {
        // 0 datasets... cannot continue
        $("#dataset-selection-submit").prop("disabled", true);
    }

    if (datasetsSelected > 1) {
        $("#name-profile").show();
    } else {
        $("#name-profile").hide();
        $("#name-profile").empty();
    }

});


/* Modal operations */
$(".js-modal-trigger").click(function() {
    const target = $(this).data('target');  // modal element ID
    openModal($(`#${target}`));
});

$(".modal-background, .modal-close").click(function() {
    closeModal($(this).closest(".modal"));
});

$(document).ready( () => {
    $("#name-profile").hide();
    fetchNemoArchiveFiles();
    resetPagination();
    $("#selected-datasets").text(datasetsSelected);
});