// Load Mermaid from CDN and render every <pre class="mermaid"> block.
// pymdownx.superfences emits the wrapper; mermaid handles the SVG.
//
// We use the ESM build and call run() once on DOMContentLoaded so the
// MkDocs Material instant-loading nav (which swaps content without
// reload) doesn't leave un-rendered diagrams behind.
(function () {
  function inject() {
    if (window.__mermaidLoaded) return;
    window.__mermaidLoaded = true;
    var script = document.createElement('script');
    script.type = 'module';
    script.textContent =
      "import mermaid from 'https://cdn.jsdelivr.net/npm/mermaid@10/dist/mermaid.esm.min.mjs';" +
      "mermaid.initialize({startOnLoad: false, theme: 'base', " +
      "themeVariables: {fontFamily: 'IBM Plex Sans, sans-serif'," +
      "fontSize: '13px', primaryColor: '#fafaf7', primaryTextColor: '#0e0e10'," +
      "primaryBorderColor: '#0e0e10', lineColor: '#0e0e10', secondaryColor: '#f3efe7'," +
      "tertiaryColor: '#fff'}});" +
      "window.__mermaid = mermaid;" +
      "window.__renderMermaid = function () {" +
      "  document.querySelectorAll('pre.mermaid').forEach(function (el) {" +
      "    if (el.getAttribute('data-processed') === 'true') return;" +
      "    var src = el.textContent;" +
      "    var div = document.createElement('div');" +
      "    div.className = 'mermaid';" +
      "    div.textContent = src;" +
      "    el.replaceWith(div);" +
      "  });" +
      "  mermaid.run({querySelector: 'div.mermaid'});" +
      "};" +
      "window.__renderMermaid();";
    document.head.appendChild(script);
  }
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', inject);
  } else {
    inject();
  }
  // Material's instant nav swaps content without reload — re-render
  // diagrams on every navigation by listening to its document$ stream.
  if (typeof document$ !== 'undefined') {
    document$.subscribe(function () {
      if (window.__renderMermaid) window.__renderMermaid();
    });
  }
})();
