const puppeteer = require('puppeteer');

// (async () => {
//   console.log("test-bp-10K start");
//   const browser = await puppeteer.launch();
//   const page = await browser.newPage();
//   page.on('console', msg => console.log('PAGE LOG:', msg.text()));
  
//   const start = Date.now();
//   await page.goto('http://localhost:8081/components/epiviz-charts/performance/tests/test-bp-10K.html', 
//         {waitUntil: 'networkidle0'});
//   const end = Date.now();
//   console.log("before loading Page : " + start);
//   console.log("After Page Load : " + end);
//   console.log("Total Elapsed Time : " + (end - start));
//   const metrics = await page.metrics();
//   console.log(JSON.stringify(metrics));
//   await page.screenshot({path: 'test.png'});
//   await browser.close();
//   console.log("test-bp-10K end");
// })();

// (async () => {
//   console.log("test-bp-100K start");
//   const browser = await puppeteer.launch();
//   const page = await browser.newPage();
//   page.on('console', msg => console.log('PAGE LOG:', msg.text()));
  
//   const start = Date.now();
//   await page.goto('http://localhost:8081/components/epiviz-charts/performance/tests/test-20K.html', 
//         {waitUntil: 'networkidle0'});
//   const end = Date.now();
//   console.log("before loading Page : " + start);
//   console.log("After Page Load : " + end);
//   console.log("Total Elapsed Time : " + (end - start));
//   const metrics = await page.metrics();
//   console.log(JSON.stringify(metrics));
//   await page.screenshot({path: 'test.png'});
//   await browser.close();
//   console.log("test-bp-100K end");
// })();

(async () => {
  console.log("test-bp-1M start");
  const browser = await puppeteer.launch();
  const page = await browser.newPage();
  page.on('console', msg => console.log('PAGE LOG:', msg.text()));
  
  const start = Date.now();
  await page.goto('http://localhost:8081/components/epiviz-charts/performance/tests/test-30K.html', 
        {waitUntil: 'networkidle0'});
  const end = Date.now();
  console.log("before loading Page : " + start);
  console.log("After Page Load : " + end);
  console.log("Total Elapsed Time : " + (end - start));
  const metrics = await page.metrics();
  console.log(JSON.stringify(metrics));
  await page.screenshot({path: 'test.png'});
  await browser.close();
  console.log("test-bp-1M end");
})();

(async () => {
  console.log("test-bp-100M start");
  const browser = await puppeteer.launch();
  const page = await browser.newPage();
  page.on('console', msg => console.log('PAGE LOG:', msg.text()));
  
  const start = Date.now();
  await page.goto('http://localhost:8081/components/epiviz-charts/performance/tests/test-40K.html', 
        {waitUntil: 'networkidle0'});
  const end = Date.now();
  console.log("before loading Page : " + start);
  console.log("After Page Load : " + end);
  console.log("Total Elapsed Time : " + (end - start));
  const metrics = await page.metrics();
  console.log(JSON.stringify(metrics));
  await page.screenshot({path: 'test.png'});
  await browser.close();
  console.log("test-bp-100M end");
})();