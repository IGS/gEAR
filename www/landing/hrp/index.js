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

      // Filter out papers without titles (usually just the template entry at the end)
      papers = papers.filter(paper => paper.title && paper.title.trim() !== '');

      // Update the counter display
      document.getElementById('total-paper-count').textContent = papers.length;
      document.getElementById('current-paper-count').textContent = currentPaperIndex + 1;

      if (papers.length > 0) {
        displayPaper(currentPaperIndex); // Set the initial paper
      }

      // shuffle all papers except the first one
      for (let i = papers.length - 1; i > 1; i--) {
        const j = Math.floor(Math.random() * (i - 1)) + 1; // random index from 1 to i
        [papers[i], papers[j]] = [papers[j], papers[i]]; // swap elements
      }
    })
    .catch(error => console.error('Error fetching papers:', error));

    // Event listener for the next button
    document.getElementById('next-paper-btn').addEventListener('click', () => {
        currentPaperIndex = (currentPaperIndex + 1) % papers.length;
        document.getElementById('current-paper-count').textContent = currentPaperIndex + 1;
        displayPaper(currentPaperIndex);
    });

    // Event listener for the previous button
    document.getElementById('prev-paper-btn').addEventListener('click', () => {
        currentPaperIndex = (currentPaperIndex - 1 + papers.length) % papers.length;
        document.getElementById('current-paper-count').textContent = currentPaperIndex + 1;
        displayPaper(currentPaperIndex);
    });
};

const displayPaper = (index) => {
    const paper = papers[index];
    if (!paper) return;

    document.getElementById('paper-title').textContent = paper.title;
    document.getElementById('paper-authors').textContent = paper.authors;
    document.getElementById('paper-journal').textContent = paper.journal;

    if (paper.doi) {
        document.getElementById('paper-doi').textContent = `DOI: ${paper.doi}`;
    }

    setLink('landing-link', paper.landingLink);
    setLink('read-paper-link', pubmed_url_base + paper.pubmedId);
    
    setLink('view-data-link', paper.dataLink);
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