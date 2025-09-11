'use strict';

let currentPaperIndex = 0;
let papers = [];
let pubmed_url_base = 'https://pubmed.ncbi.nlm.nih.gov/';


window.onload=function() {
    // Set the page title
    document.getElementById('page-header-label').textContent = 'HRP Landing page';

    // Fetch papers from JSON file
    fetch('/landing/hrp/papers.json')
    .then(response => response.json())
    .then(data => {
      papers = data;
      if (papers.length > 0) {
      displayPaper(currentPaperIndex); // Set the initial paper
      }
    })
    .catch(error => console.error('Error fetching papers:', error));

    // Event listener for the next button
    document.getElementById('next-paper-btn').addEventListener('click', () => {
        currentPaperIndex = (currentPaperIndex + 1) % papers.length;
        displayPaper(currentPaperIndex);
    });
};

const displayPaper = (index) => {
    const paper = papers[index];
    if (!paper) return;

    document.getElementById('paper-title').textContent = paper.title;
    document.getElementById('paper-authors').textContent = paper.authors;
    document.getElementById('paper-journal').innerHTML = `${paper.journal} <a href="${paper.doi}" target="_blank">doi:${paper.doi}</a>`;

    setLink('landing-link', paper.landingLink);
    setLink('read-paper-link', pubmed_url_base + paper.pubmedID);
    
    setLink('view-data-link', paper.viewDataLink);
    setLink('download-data-link', paper.downloadLink);
    
    document.getElementById('paper-image').src = paper.image;
    document.getElementById('paper-image').alt = paper.imageAlt || 'Figure';
};

const handlePageSpecificLoginUIUpdates = async (event) => {
    
}

const setLink = (elementId, link) => {
    const element = document.getElementById(elementId);
    if (link) {
        element.href = link;
        element.parentElement.style.display = 'block'; // Show the link
    } else {
        element.parentElement.style.display = 'none'; // Hide the link
    }
};