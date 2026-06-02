// Load Mermaid from CDN and render every <pre class="mermaid"> block.
// pymdownx.superfences emits the wrapper; mermaid handles the SVG.
//
// We use the ESM build and call run() once on DOMContentLoaded so the
// MkDocs Material instant-loading nav (which swaps content without
// reload) doesn't leave un-rendered diagrams behind.
//
// Theme choice: 'neutral' (clean light palette, respects per-node
// inline color attributes and classDef overrides). We do NOT set
// primaryTextColor in themeVariables, because doing so makes Mermaid
// use it as the default text color for ALL nodes, ignoring per-node
// `color:` attributes on dark-fill boxes and leaving dark text on
// dark backgrounds (unreadable). Diagrams that need dark fills use
// `classDef accent fill:#1f8a4c,color:#fff` / `class N accent` to
// pin readable text per-node.
(function () {
  function inject() {
    if (window.__mermaidLoaded) return;
    window.__mermaidLoaded = true;
    var script = document.createElement('script');
    script.type = 'module';
    script.textContent =
      "import mermaid from 'https://cdn.jsdelivr.net/npm/mermaid@10/dist/mermaid.esm.min.mjs';" +
      "mermaid.initialize({startOnLoad: false, theme: 'neutral', " +
      "themeVariables: {fontFamily: 'IBM Plex Sans, sans-serif'," +
      "fontSize: '13px', lineColor: '#0e0e10'}});" +
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
