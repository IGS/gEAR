/**
 * Created with JetBrains PhpStorm.
 * User: florin
 * Date: 4/12/13
 * Time: 11:12 PM
 * To change this template use File | Settings | File Templates.
 */


(
  function(i,s,o,g,r,a,m) {
    i['GoogleAnalyticsObject'] = r;
    i[r] = i[r] ||
      function() {
       (i[r].q = i[r].q || []).push(arguments)
      },
      i[r].l = 1 * new Date();
    a = s.createElement(o),
      m = s.getElementsByTagName(o)[0];
    a.async = 1;
    a.src = g;
    m.parentNode.insertBefore(a,m)
  }
)(window, document, 'script','//www.google-analytics.com/analytics.js','ga');

ga('create', 'UA-40113910-1', 'umd.edu');
ga('send', 'pageview');
