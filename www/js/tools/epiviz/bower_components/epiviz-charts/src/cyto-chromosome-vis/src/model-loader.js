
(function (cyto_chr, d3) {

  var defaultDataURLs = {
    "400" : "ideogram_9606_GCF_000001305.14_400_V1",
    "550" : "ideogram_9606_GCF_000001305.14_550_V1",
    "850" : "ideogram_9606_GCF_000001305.14_850_V1",
    "1200" : "ideogram_9606_GCF_000001305.13_1200_v1"
  };

  var baseDir = 'data/';

  function CacheInstance() {
    this.status = "notloaded";
    this.cache = [];
  }

  var dataCache = {
    "400" : new CacheInstance,
    "550" : new CacheInstance,
    "850" : new CacheInstance,
    "1200" : new CacheInstance
  };

  var callQueue = [];

  function loadData(file, res, cb) {

    var c = dataCache[res];
    if (c.cache.length === 0) {
      if (c.status === "loading") {
        callQueue.push({
          res: res,
          cb: cb
        });

        return;
      } else if (c.status === "notloaded") {
        c.status = "loading";
        d3.tsv(file, function(d) {
          c.cache = d;
          c.status = "loaded";
          cb(d);

          while(callQueue.length > 0) {
            var cbq = callQueue.shift();
            cbq.cb(d);
          }

        });
      }
    } else {
      cb(c.cache);
    }
  }

  function getChromosomeData(chr, resolution, cb) {

    chr = chr || '1';
    resolution = resolution || "550";

    var fileName = defaultDataURLs[resolution];

    loadData(baseDir + fileName, resolution, function (d) {

      if(d) {
        var filteredResults = filterByChromosome(d, chr);
        cb(filteredResults);
      }

    });
  }

  function filterByChromosome(data, chr) {
    var newAry = [];
    for(var i = 0; i < data.length; i++) {
      if (data[i]['#chromosome'] === chr) {
        newAry.push(data[i]);
      }
    }
    return newAry;
  }

  cyto_chr.modelLoader = {
    load: getChromosomeData,
    setDataDir: function(d) {baseDir = d;},
    getDataDir: function() {return baseDir;}
  };

})(cyto_chr || {}, d3);