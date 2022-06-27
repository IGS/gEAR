window.onload=function() {
    // This matches the choices in #list-tab
    //  'basic' corresponds to '#list-basic-choice', 'curation' to '#list-curation-choice', etc.
    const  doc_link = getUrlParameter('doc');
    if (doc_link) {
        $(`#list-${doc_link}-choice`).trigger('click');
    }
    deferIframe();
}

// Defer iFrame loading until after the page has loaded.  Speeds up loading of navbar and other things
// Source - https://www.annacantrell.com/how-to/how-to-defer-loading-iframe-content/
function deferIframe() {
    const iframeElem = document.getElementsByTagName('iframe');
    for ( let i = 0; i < iframeElem.length; i++ ) {
      if(iframeElem[i].getAttribute('data-src')) {
        iframeElem[i].setAttribute('src',iframeElem[i].getAttribute('data-src'));
      }
    }
  }