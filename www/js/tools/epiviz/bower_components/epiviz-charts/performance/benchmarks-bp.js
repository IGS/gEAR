const puppeteer = require("puppeteer");

async function multiple_run_benchmark(html, nruns, file_path) {
  test_runs = [];
  datafields = ["domInteractive", "domComplete"];

  for (var i = 0; i < nruns; i++) {
    const browser = await puppeteer.launch({
      args: ["--no-sandbox", "--disable-setuid-sandbox"]
    });
    const page = await browser.newPage();

    var chart_draw_time = null;
    await page.exposeFunction("onCustomEvent", e => {
      chart_draw_time = e.detail;
    });

    function listenFor(type) {
      return page.evaluateOnNewDocument(type => {
        document.addEventListener(type, e => {
          window.onCustomEvent({ type, detail: e.detail });
        });
      }, type);
    }
    await listenFor("total_draw_time");

    var responseMap = {};
    var responseTimingMap = {};
    const client = await page.target().createCDPSession();
    await client.send("Network.enable");

    //  Track requests
    client.on("Network.requestWillBeSent", event => {
      var requestId = event.requestId;
      var url = event.request.url;

      if (
        url.indexOf("?requestId=") != -1 &&
        url.indexOf("&action=getValues") != -1
      ) {
        responseMap[requestId] = {
          requestSent: [],
          responseReceived: [],
          dataReceived: [],
          requestFinished: [],
          requestSize: []
        };
        responseMap[requestId].requestSent.push(Date.now());
      }
    });

    //  data chunks received
    client.on("Network.dataReceived", event => {
      var requestId = event.requestId;

      if (responseMap[requestId]) {
        responseMap[requestId].dataReceived.push(Date.now());
      }
    });

    client.on("Network.responseReceived", event => {
      var response = event.response;
      var requestId = event.requestId;

      if (responseMap[requestId]) {
        responseMap[requestId].responseReceived.push(Date.now());
        responseTimingMap[requestId] = response;
      }
    });

    client.on("Network.loadingFinished", event => {
      var requestId = event.requestId;

      if (responseMap[requestId]) {
        responseMap[requestId].requestFinished.push(Date.now());
        responseMap[requestId].requestSize.push(event.encodedDataLength);
      }
    });

    await page.goto(html, { waitUntil: "networkidle0" });
    const performanceTiming = JSON.parse(
      await page.evaluate(() => JSON.stringify(window.performance.timing))
    );

    await client.detach();

    // await page.screenshot({ path: file_path + "-" + i + ".png" });
    await browser.close();

    const navigationStart = performanceTiming.navigationStart;

    var total_http_time = 0;
    var total_latency_time = 0;
    var total_request_size = 0;

    Object.keys(responseMap).forEach(function(reqId) {
      total_http_time +=
        responseMap[reqId].requestFinished - responseMap[reqId].requestSent;
      total_latency_time +=
        responseMap[reqId].responseReceived - responseMap[reqId].requestSent;
      total_request_size += responseMap[reqId].requestSize[0];
    });

    const extractedData = {
      total_draw_time: chart_draw_time.total_draw_time,
      total_http_time: total_http_time / Object.keys(responseMap).length,
      total_latency_time: total_latency_time / Object.keys(responseMap).length,
      total_request_size: total_request_size / Object.keys(responseMap).length
    };

    datafields.forEach(name => {
      extractedData[name] = performanceTiming[name] - navigationStart;
    });

    test_runs.push(extractedData);
  }

  return test_runs;
}

(async () => {
  console.log(
    '{ "10K": ' +
      JSON.stringify(
        await multiple_run_benchmark(
          "http://localhost:8081/components/epiviz-charts/performance/tests/test-bp-10K.html",
          10,
          "./screenshots/test10k"
        )
      ) +
      ","
  );
  console.log(
    ' "100K": ' +
      JSON.stringify(
        await multiple_run_benchmark(
          "http://localhost:8081/components/epiviz-charts/performance/tests/test-bp-100K.html",
          10,
          "./screenshots/test100k"
        )
      ) +
      ","
  );
  console.log(
    ' "1M": ' +
      JSON.stringify(
        await multiple_run_benchmark(
          "http://localhost:8081/components/epiviz-charts/performance/tests/test-bp-1M.html",
          10,
          "./screenshots/test1M"
        )
      ) +
      ","
  );
  console.log(
    ' "10M": ' +
      JSON.stringify(
        await multiple_run_benchmark(
          "http://localhost:8081/components/epiviz-charts/performance/tests/test-bp-10M.html",
          10,
          "./screenshots/test10M"
        )
      ) +
      ","
  );
  console.log(
    ' "100M": ' +
      JSON.stringify(
        await multiple_run_benchmark(
          "http://localhost:8081/components/epiviz-charts/performance/tests/test-bp-100M.html",
          10,
          "./screenshots/test100M"
        )
      ) +
      ","
  );
  console.log(
    ' "chr": ' +
      JSON.stringify(
        await multiple_run_benchmark(
          "http://localhost:8081/components/epiviz-charts/performance/tests/test-bp-chr.html",
          10,
          "./screenshots/testchr"
        )
      ) +
      "}"
  );
})();
