/* GenAIRR mockup — shared JS for visualizers */

/* ----------------------------------------
   Sequence rendering helper
   tokens: array of { c: 'A', s: 'v'|'d'|'j'|'np'|'mut'|'n'|'pcr'|'del'|'ins'|'umi' }
   ---------------------------------------- */
function renderSeq(tokens) {
  return tokens.map(t => `<span class="b b-${t.s}">${t.c}</span>`).join('');
}

/* Make a uniform random N-mer */
function rndSeq(len, seed = 1) {
  const bases = ['A','C','G','T'];
  let s = '';
  let x = seed;
  for (let i = 0; i < len; i++) {
    x = (x * 9301 + 49297) % 233280;
    s += bases[Math.floor((x / 233280) * 4)];
  }
  return s;
}

/* ----------------------------------------
   Canonical example sequence used across the site
   so the user sees the SAME molecule transformed by
   each phase. This builds intuition.
   ---------------------------------------- */
const CANONICAL = (() => {
  const v   = 'GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCGTGGTCCAGCCTGGGAGGTCC';
  const np1 = 'CCGTA';
  const d   = 'GGGTATAGCAGCTCG';
  const np2 = 'TACCG';
  const j   = 'ACTACTGGTACTTCGATCTC';
  return { v, np1, d, np2, j };
})();

function tokensFor(stage, opts = {}) {
  const { v, np1, d, np2, j } = CANONICAL;
  const out = [];
  const push = (str, s) => { for (const c of str) out.push({ c, s }); };

  // Stage 1: V only
  if (stage >= 1) push(v, 'v');
  // Stage 2: + NP1 + D
  if (stage >= 2) { push(np1, 'np'); push(d, 'd'); }
  // Stage 3: + NP2 + J  → complete recombination
  if (stage >= 3) { push(np2, 'np'); push(j, 'j'); }

  // Stage 4: somatic hypermutation
  if (stage >= 4 && opts.mutate !== false) {
    const mutPos = [9, 18, 26, 38, 71, 78];
    const mutBases = ['T','C','A','G','C','T'];
    let i = 0;
    for (const p of mutPos) {
      if (out[p] && out[p].s !== 'np') {
        out[p] = { c: mutBases[i++ % mutBases.length], s: 'mut' };
      }
    }
  }

  // Stage 5: corruption — 5'/3' loss, N-bases, a single insertion
  if (stage >= 5) {
    // 5' loss: trim first 8 bases
    out.splice(0, 8);
    // 3' loss: trim last 4 bases
    out.splice(out.length - 4, 4);
    // scatter a few low-quality N-bases
    [10, 25, 50, 70].forEach(p => { if (out[p]) out[p] = { c: 'N', s: 'n' }; });
    // a single inserted base
    if (out[35]) out.splice(35, 0, { c: 'A', s: 'ins' });
  }

  return out;
}

/* ----------------------------------------
   Pipeline scrubber
   ---------------------------------------- */
function setupScrubber(rootSelector) {
  const root = document.querySelector(rootSelector);
  if (!root) return;
  const stops = root.querySelectorAll('.scrubber-stop');
  const labels = root.querySelectorAll('.scrubber-labels span');
  const seqEl = root.querySelector('.scrubber-seq');
  const meta = root.querySelector('.scrubber-meta');

  const stages = [
    { name: 'Empty',     desc: 'Pre-rearrangement',                  stage: 0, code: '' },
    { name: 'Recombine', desc: 'V·NP1·D·NP2·J assembled',            stage: 3, code: '.recombine()' },
    { name: 'Mutate',    desc: 'S5F somatic hypermutation applied',  stage: 4, code: '.mutate(model="s5f", count=(5, 25))' },
    { name: 'Corrupt',   desc: "5'/3' loss · N-bases · indels",      stage: 5,
      code:
        '.corrupt_5prime_loss(length=(5, 30))\n' +
        '.corrupt_3prime_loss(length=(0, 15))\n' +
        '.corrupt_ns(count=(0, 5))\n' +
        '.corrupt_indels(count=(0, 2), insertion_prob=0.5)' },
  ];

  let active = 1;

  function render() {
    stops.forEach((s, i) => {
      s.classList.toggle('active', i === active);
      s.classList.toggle('passed', i < active);
    });
    labels.forEach((l, i) => l.classList.toggle('active', i === active));
    const stg = stages[active];
    const tokens = tokensFor(stg.stage, { mutate: stg.stage >= 4 });
    if (stg.stage === 0) {
      seqEl.innerHTML = '<span style="color:var(--ink-muted)">— no sequence yet —</span>';
    } else {
      seqEl.innerHTML = renderSeq(tokens);
    }

    // metadata
    const len = tokens.length;
    const muts = tokens.filter(t => t.s === 'mut').length;
    const ns = tokens.filter(t => t.s === 'n').length;
    const ins = tokens.filter(t => t.s === 'ins').length;
    meta.innerHTML = `
      <div class="m"><div class="k">phase</div><div class="v">${stg.name}</div></div>
      <div class="m"><div class="k">length</div><div class="v">${len || 0}</div></div>
      <div class="m"><div class="k">mutations</div><div class="v">${muts}</div></div>
      <div class="m"><div class="k">N-bases</div><div class="v">${ns}</div></div>
      <div class="m"><div class="k">indels</div><div class="v">${ins}</div></div>
    `;

    const codeBlock = root.parentElement.querySelector('.scrubber-code code');
    if (codeBlock) {
      // chain methods are indented to 7 spaces to match the lesson code template
      codeBlock.textContent = stages.slice(1, active + 1)
        .map(s => s.code.split('\n').map(l => '       ' + l).join('\n'))
        .join('\n');
    }
  }

  stops.forEach((s, i) => s.addEventListener('click', () => { active = i; render(); }));
  render();
}

/* ----------------------------------------
   V(D)J recombination workbench (hero)
   ---------------------------------------- */
function setupBench(rootSelector) {
  const root = document.querySelector(rootSelector);
  if (!root) return;

  // ---- realistic allele names sampled from human IGH ----
  const V_LIB = [
    'IGHV1-69*06','IGHV3-23*04','IGHV4-39*02','IGHV3-30*18','IGHV1-2*02',
    'IGHV3-7*01','IGHV4-59*01','IGHV3-48*03','IGHV1-18*01','IGHV2-5*10',
    'IGHV3-33*05','IGHV5-51*01','IGHV3-15*07','IGHV1-46*02','IGHV4-34*01',
  ];
  const D_LIB = [
    'IGHD3-22*01','IGHD2-21*02','IGHD3-10*01','IGHD6-13*01','IGHD2-2*02',
    'IGHD1-26*01','IGHD4-17*01','IGHD3-3*01','IGHD5-12*01','IGHD2-15*01',
    'IGHD6-19*01','IGHD3-9*01','IGHD1-1*01','IGHD7-27*01',
  ];
  const J_LIB = [
    'IGHJ4*02','IGHJ6*03','IGHJ5*02','IGHJ3*02','IGHJ1*01',
    'IGHJ2*01','IGHJ4*01','IGHJ6*02','IGHJ5*01',
  ];

  // ---- DOM refs ----
  const seedEl  = root.querySelector('#bench-seed');
  const stepNum = root.querySelector('#bench-step-num');
  const stepName = root.querySelector('#bench-step-name');
  const stepFn  = root.querySelector('#bench-step-fn');
  const rackV = root.querySelector('#rack-v');
  const rackD = root.querySelector('#rack-d');
  const rackJ = root.querySelector('#rack-j');
  const rackColV = root.querySelector('[data-rack="v"]');
  const rackColD = root.querySelector('[data-rack="d"]');
  const rackColJ = root.querySelector('[data-rack="j"]');
  const ruler = root.querySelector('#bench-ruler');
  const track = root.querySelector('#bench-track');
  const prov  = root.querySelector('#bench-prov');
  const fields = root.querySelector('#bench-fields');
  const fieldByKey = {};
  fields.querySelectorAll('.bf').forEach(el => {
    const k = el.querySelector('.bf-k').textContent.trim();
    fieldByKey[k] = el.querySelector('.bf-v');
  });

  // ---- helpers ----
  const $ = (sel, parent = root) => parent.querySelector(sel);
  const sleep = ms => new Promise(r => setTimeout(r, ms));
  let cancelled = false;

  function setStep(num, name, fn) {
    stepNum.textContent = String(num).padStart(2, '0');
    stepName.textContent = name;
    stepFn.textContent = fn;
  }

  function setField(key, value) {
    const el = fieldByKey[key];
    if (!el) return;
    el.textContent = value;
    el.classList.remove('is-fresh');
    void el.offsetWidth;
    el.classList.add('is-fresh');
  }

  function fillRack(rackEl, lib, picked) {
    rackEl.innerHTML = '';
    // show 5 names: picked in middle if possible
    const idx = lib.indexOf(picked);
    const start = Math.max(0, Math.min(lib.length - 5, idx - 2));
    for (let i = start; i < start + 5 && i < lib.length; i++) {
      const item = document.createElement('div');
      item.className = 'rack-item';
      item.textContent = lib[i];
      if (lib[i] === picked) item.dataset.picked = '1';
      rackEl.appendChild(item);
    }
  }

  async function flickerPick(rackEl, lib, picked) {
    fillRack(rackEl, lib, picked);
    const items = [...rackEl.querySelectorAll('.rack-item')];
    // brief flicker — short enough to survive throttling
    for (let i = 0; i < 4 && !cancelled; i++) {
      items.forEach(it => it.classList.remove('is-flicker'));
      const it = items[Math.floor(Math.random() * items.length)];
      it.classList.add('is-flicker');
      await sleep(120);
    }
    items.forEach(it => it.classList.remove('is-flicker'));
    // settle on the picked one
    const target = items.find(it => it.dataset.picked === '1');
    if (target) target.classList.add('is-picked');
  }

  function drawRuler(totalLen) {
    ruler.innerHTML = '';
    const w = ruler.clientWidth;
    const stops = [0, Math.round(totalLen * 0.25), Math.round(totalLen * 0.5),
                   Math.round(totalLen * 0.75), totalLen];
    stops.forEach((bp, i) => {
      const t = document.createElement('div');
      t.className = 'tick';
      t.style.left = `${(bp / totalLen) * 100}%`;
      if (i === 0) t.style.transform = 'translateX(0)';
      else if (i === stops.length - 1) t.style.transform = 'translateX(-100%)';
      else t.style.transform = 'translateX(-50%)';
      t.textContent = `${bp}`;
      ruler.appendChild(t);
    });
  }

  function clearTrack() {
    track.innerHTML = '';
    prov.innerHTML = '';
    prov.classList.remove('is-shown');
  }

  // segments are placed by % offsets relative to a 'planned total length'
  // we keep a running model of segments [{key, start, len, label, finalLen}]
  function renderSegments(segs, total, opts = {}) {
    track.querySelector('.molecule-empty')?.remove();
    // track existing nodes by key
    const existing = {};
    track.querySelectorAll('.seg').forEach(el => existing[el.dataset.key] = el);
    segs.forEach(s => {
      let el = existing[s.key];
      if (!el) {
        el = document.createElement('div');
        el.className = `seg seg-${s.key}`;
        el.dataset.key = s.key;
        el.innerHTML = `<span class="seg-label">${s.label}</span>`;
        track.appendChild(el);
      }
      el.style.left  = `${(s.start / total) * 100}%`;
      el.style.width = `${(s.len   / total) * 100}%`;
      el.style.opacity = s.len > 0 ? 1 : 0;
    });
    // remove gone segments
    Object.keys(existing).forEach(k => {
      if (!segs.find(s => s.key === k)) existing[k].remove();
    });
  }

  function renderProvenance(segs, total) {
    prov.innerHTML = '';
    segs.forEach(s => {
      if (s.len <= 0) return;
      const pv = document.createElement('div');
      const cls = s.key.startsWith('np') ? 'pv-np' : `pv-${s.key}`;
      pv.className = `pv ${cls}`;
      pv.style.left  = `${(s.start / total) * 100}%`;
      pv.style.width = `${(s.len   / total) * 100}%`;
      pv.textContent = s.provLabel || s.label;
      prov.appendChild(pv);
    });
    prov.classList.add('is-shown');
  }

  async function showTrimChip(text, leftPct) {
    const chip = document.createElement('div');
    chip.className = 'trim-chip';
    chip.textContent = text;
    chip.style.left = `calc(${leftPct}% - 30px)`;
    chip.style.top = '-20px';
    track.appendChild(chip);
    await sleep(20);
    chip.classList.add('is-shown');
    await sleep(900);
    chip.classList.remove('is-shown');
    await sleep(250);
    chip.remove();
  }

  // ---- the main loop ----
  async function runOnce(seed) {
    cancelled = false;
    if (seedEl) seedEl.textContent = seed;

    // pick alleles deterministically from seed
    const pickV = V_LIB[seed % V_LIB.length];
    const pickD = D_LIB[(seed * 7) % D_LIB.length];
    const pickJ = J_LIB[(seed * 13) % J_LIB.length];

    // realistic untrimmed lengths
    const vLenRaw = 296;
    const dLenRaw = 17 + (seed % 8);   // 17–24
    const jLenRaw = 50 + (seed % 6);   // 50–55
    const np1Len = 3 + (seed % 7);     // 3–9
    const np2Len = 2 + ((seed * 5) % 8); // 2–9
    const vTrim3 = 1 + (seed % 4);     // 1–4
    const dTrim5 = 1 + ((seed * 3) % 5);
    const dTrim3 = 1 + ((seed * 2) % 4);
    const jTrim5 = 2 + ((seed * 11) % 6);

    const vLen = vLenRaw - vTrim3;
    const dLen = dLenRaw - dTrim5 - dTrim3;
    const jLen = jLenRaw - jTrim5;
    const total = vLen + np1Len + dLen + np2Len + jLen;
    const productive = (total % 3 === 0) || ((seed * 7) % 3 !== 0);

    // reset
    [rackV, rackD, rackJ].forEach(r => r.innerHTML = '');
    [rackColV, rackColD, rackColJ].forEach(c => c.classList.remove('is-active', 'is-done'));
    clearTrack();
    track.innerHTML = '<div class="molecule-empty">— pre-rearrangement —</div>';
    Object.values(fieldByKey).forEach(el => { el.textContent = '—'; el.classList.remove('is-fresh'); });
    drawRuler(total);

    // step 1 — sample V
    setStep(1, 'sample V allele from germline library', 'SimulateSequence.sample_v()');
    rackColV.classList.add('is-active');
    await flickerPick(rackV, V_LIB, pickV);
    if (cancelled) return;
    setField('v_call', pickV);
    let segs = [{ key: 'v', start: 0, len: vLenRaw, label: 'V', provLabel: `v_call · ${vLenRaw} bp` }];
    renderSegments(segs, total);
    await sleep(550);

    // step 2 — trim V 3'
    setStep(2, `trim V 3' end (${vTrim3} bp removed)`, 'sample_v_trim_3()');
    showTrimChip(`v_trim_3 = ${vTrim3}`, (vLen / total) * 100);
    segs = [{ key: 'v', start: 0, len: vLen, label: 'V', provLabel: `v_call · ${vLen} bp` }];
    renderSegments(segs, total);
    setField('trim_3 v · trim_5/3 d · trim_5 j', `${vTrim3} · ${dTrim5}/${dTrim3} · ${jTrim5}`);
    await sleep(900);
    rackColV.classList.add('is-done');
    rackColV.classList.remove('is-active');

    // step 3 — generate NP1
    setStep(3, `generate NP1 region (${np1Len} bp, Markov chain)`, 'generate_NP1()');
    segs.push({ key: 'np1', start: vLen, len: np1Len, label: 'NP1', provLabel: `np1 · ${np1Len}` });
    renderSegments(segs, total);
    await sleep(800);

    // step 4 — sample D
    setStep(4, 'sample D allele', 'SimulateSequence.sample_d()');
    rackColD.classList.add('is-active');
    await flickerPick(rackD, D_LIB, pickD);
    if (cancelled) return;
    setField('d_call', pickD);
    segs.push({ key: 'd', start: vLen + np1Len, len: dLenRaw, label: 'D', provLabel: `d_call · ${dLenRaw}` });
    renderSegments(segs, total);
    await sleep(550);

    // step 5 — trim D both ends
    setStep(5, `trim D ends (5'=${dTrim5}, 3'=${dTrim3})`, 'sample_d_trim_5() / sample_d_trim_3()');
    const dStart = vLen + np1Len + dTrim5;
    showTrimChip(`d_trim_5 = ${dTrim5}`, (dStart / total) * 100);
    segs = segs.filter(s => s.key !== 'd');
    segs.push({ key: 'd', start: dStart, len: dLen, label: 'D', provLabel: `d_call · ${dLen}` });
    renderSegments(segs, total);
    await sleep(900);
    rackColD.classList.add('is-done');
    rackColD.classList.remove('is-active');

    // step 6 — generate NP2
    setStep(6, `generate NP2 region (${np2Len} bp)`, 'generate_NP2()');
    segs.push({ key: 'np2', start: dStart + dLen, len: np2Len, label: 'NP2', provLabel: `np2 · ${np2Len}` });
    renderSegments(segs, total);
    setField('np1 / np2', `${np1Len} bp / ${np2Len} bp`);
    await sleep(700);

    // step 7 — sample + trim J
    setStep(7, 'sample J allele', 'SimulateSequence.sample_j()');
    rackColJ.classList.add('is-active');
    await flickerPick(rackJ, J_LIB, pickJ);
    if (cancelled) return;
    setField('j_call', pickJ);
    const jStart = dStart + dLen + np2Len;
    segs.push({ key: 'j', start: jStart, len: jLen, label: 'J', provLabel: `j_call · ${jLen}` });
    renderSegments(segs, total);
    await sleep(550);

    // step 8 — reveal ground truth
    setStep(8, 'ground-truth annotations populated · 47 fields', 'SimulationContainer → AIRR record');
    setField('productive', productive ? 'true' : 'false');
    renderProvenance(segs, total);
    await sleep(2400);
  }

  async function loop() {
    let seed = 42;
    while (!cancelled && root.isConnected) {
      await runOnce(seed);
      if (cancelled) break;
      seed = (seed * 1103515245 + 12345) % 2147483647;
      seed = seed % 9000 + 100; // keep 3-4 digit
      await sleep(800);
    }
  }

  // Pause when off-screen so it doesn't burn CPU on long pages
  let started = false;
  function kickoff() {
    if (started) return;
    started = true;
    loop();
  }

  // Render a complete static "final state" immediately so first paint is rich
  // (in case the tab is throttled / not yet visible)
  function renderStatic() {
    fillRack(rackV, V_LIB, V_LIB[0]);
    fillRack(rackD, D_LIB, D_LIB[0]);
    fillRack(rackJ, J_LIB, J_LIB[0]);
    rackV.querySelector('[data-picked]')?.classList.add('is-picked');
    rackD.querySelector('[data-picked]')?.classList.add('is-picked');
    rackJ.querySelector('[data-picked]')?.classList.add('is-picked');
    const vL = 290, np1L = 5, dL = 13, np2L = 4, jL = 47;
    const total = vL + np1L + dL + np2L + jL;
    drawRuler(total);
    const segs = [
      { key: 'v',   start: 0,                   len: vL,   label: 'V',   provLabel: `v_call · ${vL}` },
      { key: 'np1', start: vL,                  len: np1L, label: 'NP1', provLabel: `np1 · ${np1L}` },
      { key: 'd',   start: vL + np1L,           len: dL,   label: 'D',   provLabel: `d_call · ${dL}` },
      { key: 'np2', start: vL + np1L + dL,      len: np2L, label: 'NP2', provLabel: `np2 · ${np2L}` },
      { key: 'j',   start: vL + np1L + dL + np2L, len: jL, label: 'J',   provLabel: `j_call · ${jL}` },
    ];
    track.innerHTML = '';
    renderSegments(segs, total);
    renderProvenance(segs, total);
    setStep(8, 'ground-truth annotations populated · 47 fields', 'SimulationContainer → AIRR record');
    setField('v_call', V_LIB[0]);
    setField('d_call', D_LIB[0]);
    setField('j_call', J_LIB[0]);
    setField('np1 / np2', `${np1L} bp / ${np2L} bp`);
    setField('trim_3 v · trim_5/3 d · trim_5 j', `2 · 1/2 · 3`);
    setField('productive', 'true');
  }
  renderStatic();

  // Start the animated loop only when visible (saves CPU + avoids throttle weirdness)
  const rect = root.getBoundingClientRect();
  if (rect.top < window.innerHeight && rect.bottom > 0 && !document.hidden) {
    setTimeout(kickoff, 1200);  // brief pause on the static state, then start animating
  } else if ('IntersectionObserver' in window) {
    const io = new IntersectionObserver((entries) => {
      entries.forEach(e => { if (e.isIntersecting && !document.hidden) kickoff(); });
    }, { threshold: 0.1 });
    io.observe(root);
    document.addEventListener('visibilitychange', () => {
      if (!document.hidden) kickoff();
    });
  } else {
    setTimeout(kickoff, 1200);
  }
}

// Legacy assembly animation — used by lesson-1 (kept for backward compat)
function setupAssembly(rootSelector, opts = {}) {
  const root = document.querySelector(rootSelector);
  if (!root) return;
  const rows = {
    v:   root.querySelector('[data-stage="v"]   .seq'),
    np1: root.querySelector('[data-stage="np1"] .seq'),
    d:   root.querySelector('[data-stage="d"]   .seq'),
    np2: root.querySelector('[data-stage="np2"] .seq'),
    j:   root.querySelector('[data-stage="j"]   .seq'),
    final: root.querySelector('[data-stage="final"] .seq'),
  };
  const allRows = root.querySelectorAll('.assembly-row');
  const stepName = root.querySelector('.step-name');
  const { v, np1, d, np2, j } = CANONICAL;
  const toks = (s, type) => [...s].map(c => ({ c, s: type }));
  const phases = [
    { key: 'v',   label: 'sample_v + trim_v',  html: renderSeq(toks(v, 'v')) },
    { key: 'np1', label: 'generate_NP1',       html: renderSeq(toks(np1, 'np')) },
    { key: 'd',   label: 'sample_d + trim_d',  html: renderSeq(toks(d, 'd')) },
    { key: 'np2', label: 'generate_NP2',       html: renderSeq(toks(np2, 'np')) },
    { key: 'j',   label: 'sample_j + trim_j',  html: renderSeq(toks(j, 'j')) },
    { key: 'final', label: 'assemble (concatenate ASeq)', html: renderSeq([
        ...toks(v, 'v'), ...toks(np1, 'np'), ...toks(d, 'd'), ...toks(np2, 'np'), ...toks(j, 'j')
      ]) },
  ];
  let i = 0;
  Object.values(rows).forEach(el => { if (el) el.innerHTML = ''; });
  allRows.forEach(r => r.classList.add('is-pending'));
  function step() {
    if (i >= phases.length) {
      setTimeout(() => {
        i = 0;
        Object.values(rows).forEach(el => { if (el) el.innerHTML = ''; });
        allRows.forEach(r => { r.classList.add('is-pending'); r.classList.remove('is-active'); });
        step();
      }, 2200);
      return;
    }
    const p = phases[i];
    const row = root.querySelector(`[data-stage="${p.key}"]`);
    if (row) {
      row.classList.remove('is-pending');
      allRows.forEach(r => r.classList.remove('is-active'));
      row.classList.add('is-active');
    }
    if (rows[p.key]) rows[p.key].innerHTML = p.html;
    if (stepName) stepName.innerHTML = `step <b>${p.label}</b>`;
    i++;
    setTimeout(step, 1100);
  }
  setTimeout(step, 400);
}

/* ----------------------------------------
   S5F mutation heatmap
   ---------------------------------------- */
function setupHeatmap(rootSelector) {
  const root = document.querySelector(rootSelector);
  if (!root) return;
  // 5 regions × 30 cells. Higher in CDR1/CDR2/CDR3 (S5F-realistic).
  const regions = [
    { name: 'FWR1', n: 25, base: 1 },
    { name: 'CDR1', n: 8,  base: 3 },
    { name: 'FWR2', n: 17, base: 1 },
    { name: 'CDR2', n: 8,  base: 4 },
    { name: 'FWR3', n: 38, base: 1 },
    { name: 'CDR3', n: 14, base: 5 },
    { name: 'FWR4', n: 10, base: 1 },
  ];
  let html = '';
  for (const r of regions) {
    const cells = Array.from({ length: r.n }).map((_, i) => {
      // jitter intensity around base
      const j = Math.floor(Math.random() * 3);
      const intensity = Math.max(0, Math.min(5, r.base + j - 1));
      return `<div class="cell" data-h="${intensity}"></div>`;
    }).join('');
    html += `
      <div class="heatmap-row">
        <div class="lab">${r.name}</div>
        <div class="heatmap-cells">${cells}</div>
      </div>`;
  }
  html += `<div class="heatmap-axis"><span>5'</span><span>nucleotide position</span><span>3'</span></div>`;
  root.innerHTML = html;
}

/* ----------------------------------------
   Corruption visualizer
   ---------------------------------------- */
function setupCorruption(rootSelector) {
  const root = document.querySelector(rootSelector);
  if (!root) return;
  const out = root.querySelector('.corruption-output');
  const baselineEl = root.querySelector('.corruption-baseline');
  const inputs = root.querySelectorAll('input[type=range]');
  const valEls = root.querySelectorAll('.t-val');

  const baseline = (() => {
    const { v, np1, d, np2, j } = CANONICAL;
    const out = [];
    [...v].forEach(c => out.push({ c, s: 'v' }));
    [...np1].forEach(c => out.push({ c, s: 'np' }));
    [...d].forEach(c => out.push({ c, s: 'd' }));
    [...np2].forEach(c => out.push({ c, s: 'np' }));
    [...j].forEach(c => out.push({ c, s: 'j' }));
    return out;
  })();
  baselineEl.innerHTML = renderSeq(baseline);

  function recompute() {
    const trim5 = parseInt(inputs[0].value, 10);
    const trim3 = parseInt(inputs[1].value, 10);
    const indelP = parseInt(inputs[2].value, 10);
    const nP = parseInt(inputs[3].value, 10);
    valEls[0].textContent = trim5;
    valEls[1].textContent = trim3;
    valEls[2].textContent = indelP + '%';
    valEls[3].textContent = nP + '%';

    let toks = baseline.map(t => ({ ...t }));
    if (trim5) toks.splice(0, trim5);
    if (trim3) toks.splice(toks.length - trim3, trim3);

    if (indelP > 0) {
      const n = Math.round(toks.length * indelP / 100);
      for (let k = 0; k < n; k++) {
        const p = Math.floor((k * 17 + 5) % toks.length);
        if (k % 2) toks.splice(p, 1);  // del
        else toks.splice(p, 0, { c: 'A', s: 'ins' });
      }
    }
    if (nP > 0) {
      const n = Math.round(toks.length * nP / 100);
      for (let k = 0; k < n; k++) {
        const p = Math.floor((k * 23 + 11) % toks.length);
        if (toks[p]) toks[p] = { c: 'N', s: 'n' };
      }
    }
    out.innerHTML = renderSeq(toks);

    const stats = root.querySelector('.corruption-stats');
    if (stats) {
      stats.innerHTML = `
        <span><b>${toks.length}</b> bp</span>
        <span><b>${toks.filter(t => t.s === 'n').length}</b> N's</span>
        <span><b>${toks.filter(t => t.s === 'ins').length}</b> insertions</span>
        <span><b>${baseline.length - toks.length + toks.filter(t => t.s === 'ins').length}</b> deletions</span>
      `;
    }
  }
  inputs.forEach(i => i.addEventListener('input', recompute));
  recompute();
}


/* ----------------------------------------
   Specimen opener — annotated museum card
   ---------------------------------------- */
function setupSpecimen(rootSelector) {
  const root = document.querySelector(rootSelector);
  if (!root) return;

  // Each specimen = one record showing one stage of the pipeline.
  // Coordinates are in BASE PAIRS — % positions are derived at render time.
  // Region [a, b] is half-open: a inclusive, b exclusive. Lengths add up.
  //
  // Realistic IGH proportions (post-recombination, pre-corruption):
  //   V ≈ 290 bp  ·  NP1 ≈ 5  ·  D ≈ 14  ·  NP2 ≈ 4  ·  J ≈ 47
  // Total ≈ 360 bp. Corruption shortens the read.

  const SPECIMENS = [
    {
      id: '00042', seed: 42, locus: 'human_igh', level: '★ baseline',
      headline: 'Recombination',
      narrate: 'A clean V–D–J rearrangement. Every segment, every junction nucleotide, every coordinate — recorded by construction.',
      length: 360,
      regions: [
        { kind: 'v',   label: 'V',   start: 0,   end: 290 },
        { kind: 'np1', label: 'N1',  start: 290, end: 295 },
        { kind: 'd',   label: 'D',   start: 295, end: 309 },
        { kind: 'np2', label: 'N2',  start: 309, end: 313 },
        { kind: 'j',   label: 'J',   start: 313, end: 360 },
      ],
      muts: [], ns: [], drift: false,
      callouts: {
        top: [
          { anchor: 145, k: 'v_call',  v: 'IGHV3-23*01', klass: 'is-truth' },
          { anchor: 302, k: 'd_call',  v: 'IGHD3-10*01', klass: 'is-cobalt' },
          { anchor: 336, k: 'j_call',  v: 'IGHJ4*02',    klass: 'is-plum' },
        ],
        bot: [
          { anchor: 145, k: 'v_sequence_end', v: '290' },
          { anchor: 302, k: 'junction_aa',    v: 'CARDVPYAFDIW', klass: 'is-truth' },
          { anchor: 336, k: 'productive',     v: 'True',         klass: 'is-truth' },
        ],
      },
    },
    {
      id: '00043', seed: 43, locus: 'human_igh', level: '★★ mutated',
      headline: 'Maturation',
      narrate: 'S5F somatic hypermutation hits the V region. Eleven substitutions — every position logged, every base-change recoverable.',
      length: 360,
      regions: [
        { kind: 'v',   label: 'V',   start: 0,   end: 290 },
        { kind: 'np1', label: 'N1',  start: 290, end: 295 },
        { kind: 'd',   label: 'D',   start: 295, end: 309 },
        { kind: 'np2', label: 'N2',  start: 309, end: 313 },
        { kind: 'j',   label: 'J',   start: 313, end: 360 },
      ],
      // bp positions inside V (mostly), one in J — biologically appropriate
      muts: [34, 67, 92, 118, 143, 171, 198, 224, 251, 273, 339],
      ns: [], drift: false,
      callouts: {
        top: [
          { anchor: 67,  k: 'mutation',     v: 'A→G @ 67',  klass: 'is-mut' },
          { anchor: 198, k: 'n_mutations',  v: '11',         klass: 'is-mut' },
          { anchor: 339, k: 'mutation',     v: 'C→T @ 339',  klass: 'is-mut' },
        ],
        bot: [
          { anchor: 145, k: 'v_identity',   v: '0.962',      klass: 'is-truth' },
          { anchor: 273, k: 'shm_model',    v: 's5f',        klass: 'is-cobalt' },
          { anchor: 339, k: 'mutation_rate',v: '0.031',      klass: 'is-mut' },
        ],
      },
    },
    {
      id: '00044', seed: 44, locus: 'human_igh', level: '★★★ corrupted',
      headline: 'Lab effects',
      narrate: "5' loss truncates the V. Three N-bases scatter into the read. The truth survives intact — only the observation is degraded.",
      // 28 bp lost from 5', 12 bp lost from 3' → 320 bp remaining
      length: 320,
      regions: [
        { kind: 'lost', label: "5'×",  start: 0,   end: 28  },
        { kind: 'v',    label: 'V',    start: 28,  end: 290 },
        { kind: 'np1',  label: 'N1',   start: 290, end: 295 },
        { kind: 'd',    label: 'D',    start: 295, end: 309 },
        { kind: 'np2',  label: 'N2',   start: 309, end: 313 },
        { kind: 'j',    label: 'J',    start: 313, end: 348 },
        { kind: 'lost', label: "3'×",  start: 348, end: 360 },
      ],
      muts: [67, 143, 224],
      ns:   [89, 207, 275],
      lostBp: { p5: 28, p3: 12 },
      drift: false,
      callouts: {
        top: [
          { anchor: 14,  k: "corrupt_5prime_loss", v: '28 bp',  klass: 'is-warn' },
          { anchor: 207, k: 'corrupt_ns',          v: '3 bases', klass: 'is-warn' },
          { anchor: 354, k: "corrupt_3prime_loss", v: '12 bp',  klass: 'is-warn' },
        ],
        bot: [
          { anchor: 89,  k: 'sequence_length',     v: '320 bp', klass: 'is-cobalt' },
          { anchor: 224, k: 'truth_v_call',        v: 'IGHV3-23*01', klass: 'is-truth' },
          { anchor: 336, k: 'productive',          v: 'True',  klass: 'is-truth' },
        ],
      },
    },
    {
      id: '00045', seed: 45, locus: 'human_igh', level: '★★★★ provenance',
      headline: 'Aligner drift',
      narrate: 'Heavy mutation pushes the best-match aligner call away from truth. GenAIRR records both — the sampled allele AND the evidence-based call. You can score the drift directly.',
      length: 348,
      regions: [
        { kind: 'lost', label: "5'×",  start: 0,   end: 18  },
        { kind: 'v',    label: 'V',    start: 18,  end: 290 },
        { kind: 'np1',  label: 'N1',   start: 290, end: 295 },
        { kind: 'd',    label: 'D',    start: 295, end: 309 },
        { kind: 'np2',  label: 'N2',   start: 309, end: 313 },
        { kind: 'j',    label: 'J',    start: 313, end: 360 },
      ],
      muts: [42, 58, 79, 96, 117, 138, 159, 184, 211, 247, 268, 339],
      ns: [],
      lostBp: { p5: 18, p3: 0 },
      drift: { from: 145, to: 200 }, // bp anchors of truth_v_call vs v_call
      callouts: {
        top: [
          { anchor: 145, k: 'truth_v_call',     v: 'IGHV3-23*01', klass: 'is-truth' },
          { anchor: 200, k: 'v_call (best match)', v: 'IGHV3-30*04', klass: 'is-warn' },
          { anchor: 336, k: 'truth_j_call',     v: 'IGHJ4*02',    klass: 'is-truth' },
        ],
        bot: [
          { anchor: 145, k: 'v_identity',       v: '0.892', klass: 'is-cobalt' },
          { anchor: 247, k: 'n_mutations',      v: '12',     klass: 'is-mut' },
          { anchor: 339, k: 'aligner_disagrees',v: 'true',   klass: 'is-warn' },
        ],
      },
    },
  ];

  // bp → % helper, used by every renderer
  const pct = (bp, total) => (bp / total) * 100;

  // Bases for the dim base-letter strip. Each cell maps back to a bp
  // position (cell i ↔ bp = i / cellCount * length), so mutation/N marks
  // on the strip line up with marks on the ribbon above.
  function buildBases(spec, cellCount) {
    const cycle = ['G','A','T','C','C','G','T','A','C','G','A','T','G','C','A','T'];
    const out = [];
    const mutSet = new Set(spec.muts);
    const nSet = new Set(spec.ns);
    const lostP5 = spec.lostBp ? spec.lostBp.p5 : 0;
    const lostP3End = spec.length - (spec.lostBp ? spec.lostBp.p3 : 0);
    for (let i = 0; i < cellCount; i++) {
      const bp = Math.round(i / cellCount * spec.length);
      // round each cell to nearest mut/n bp within ±2 bp
      const isMut = [...mutSet].some(m => Math.abs(m - bp) <= 2);
      const isN   = [...nSet].some(n => Math.abs(n - bp) <= 2);
      const isLost = bp < lostP5 || bp >= lostP3End;
      let cls = '';
      let ch = cycle[(i + spec.seed) % cycle.length];
      if (isLost)   { cls = 'b-lost'; ch = '·'; }
      else if (isN) { cls = 'b-n';    ch = 'N'; }
      else if (isMut){ cls = 'b-mut'; }
      out.push(`<span class="b ${cls}">${ch}</span>`);
    }
    return out.join('');
  }

  // refs
  const idEl = document.getElementById('spec-id');
  const locEl = document.getElementById('spec-locus');
  const lvlEl = document.getElementById('spec-level');
  const headlineEl = document.getElementById('spec-headline');
  const narrateEl = document.getElementById('spec-narrate');
  const dots = root.parentElement.querySelectorAll('.sp-dot');
  const ribbon = document.getElementById('ribbon');
  const ruler = document.getElementById('ribbon-ruler');
  const basesEl = document.getElementById('ribbon-bases');
  const topC = document.getElementById('callouts-top');
  const botC = document.getElementById('callouts-bot');
  const prevBtn = document.getElementById('spec-prev');
  const nextBtn = document.getElementById('spec-next');
  const regenBtn = document.getElementById('spec-regen');

  let active = 0;
  let timer = null;

  function drawRuler(totalLen) {
    ruler.innerHTML = '';
    const stops = [0, Math.round(totalLen * 0.25), Math.round(totalLen * 0.5),
                   Math.round(totalLen * 0.75), totalLen];
    stops.forEach((bp, i) => {
      const t = document.createElement('div');
      t.className = 'tick';
      t.style.left = `${(bp / totalLen) * 100}%`;
      if (i === 0) t.style.transform = 'translateX(0)';
      else if (i === stops.length - 1) t.style.transform = 'translateX(-100%)';
      else t.style.transform = 'translateX(-50%)';
      t.textContent = `${bp}`;
      ruler.appendChild(t);
    });
  }

  function renderSpecimen(spec) {
    // header
    idEl.textContent = spec.id;
    locEl.textContent = `${spec.locus} · seed=${spec.seed}`;
    lvlEl.textContent = spec.level;

    // dots
    dots.forEach((d, i) => {
      d.classList.toggle('is-active', i === active);
      d.classList.toggle('is-passed', i < active);
    });

    // ribbon regions — diff in place to keep transitions smooth
    const existingRegions = {};
    ribbon.querySelectorAll('.ribbon-region').forEach(el => {
      existingRegions[el.dataset.idx] = el;
    });
    spec.regions.forEach((r, idx) => {
      let el = existingRegions[idx];
      if (!el) {
        el = document.createElement('div');
        el.className = 'ribbon-region';
        el.dataset.idx = idx;
        ribbon.appendChild(el);
      }
      el.className = `ribbon-region is-${r.kind}`;
      el.dataset.idx = idx;
      const leftPct = pct(r.start, spec.length);
      const widthPct = pct(r.end - r.start, spec.length);
      el.style.left = `${leftPct}%`;
      el.style.width = `${widthPct}%`;
      el.textContent = widthPct >= 4 ? r.label : '';
    });
    // remove extras
    Object.entries(existingRegions).forEach(([idx, el]) => {
      if (Number(idx) >= spec.regions.length) el.remove();
    });

    // mutation + N-base marks (positions are in bp)
    ribbon.querySelectorAll('.ribbon-mut, .ribbon-n').forEach(el => el.remove());
    spec.muts.forEach(bp => {
      const m = document.createElement('div');
      m.className = 'ribbon-mut';
      m.style.left = `${pct(bp, spec.length)}%`;
      ribbon.appendChild(m);
    });
    spec.ns.forEach(bp => {
      const n = document.createElement('div');
      n.className = 'ribbon-n';
      n.style.left = `${pct(bp, spec.length)}%`;
      ribbon.appendChild(n);
    });

    // make sure scan line exists
    if (!ribbon.querySelector('.ribbon-scan')) {
      const scan = document.createElement('div');
      scan.className = 'ribbon-scan';
      ribbon.appendChild(scan);
    }

    // ruler
    drawRuler(spec.length);

    // bases strip — sized to fill width
    const cellCount = 96;
    basesEl.innerHTML = buildBases(spec, cellCount);

    // callouts — clean column layout: each callout occupies a slot
    // (left/center/right when N=3), with a hairline that bends from the
    // slot center to the actual bp anchor on the ribbon.
    function renderCallouts(side, list, container) {
      const existing = [...container.querySelectorAll('.callout')];
      existing.forEach(el => el.classList.remove('is-shown'));

      setTimeout(() => {
        container.innerHTML = '';
        const N = list.length;
        list.forEach((c, i) => {
          const el = document.createElement('div');
          el.className = `callout callout-slot-${i}`;
          // Slot column positions: 1/(2N), 3/(2N), 5/(2N) ... — even spacing
          const slotPct = (i + 0.5) / N * 100;
          // Anchor (where the hairline meets the ribbon) in % of container
          const anchorPct = pct(c.anchor, spec.length);
          el.style.setProperty('--slot-x',   `${slotPct}%`);
          el.style.setProperty('--anchor-x', `${anchorPct}%`);
          if (c.klass) el.classList.add(`callout-${c.klass}`);
          el.innerHTML = `
            <span class="callout-jog"></span>
            <span class="callout-stub"></span>
            <span class="callout-k">${c.k}</span>
            <span class="callout-v ${c.klass || ''}">${c.v}</span>`;
          // Compute the jog (horizontal stub between slot and anchor).
          // Width = |anchor - slot|, positioned at min(anchor, slot).
          const jog = el.querySelector('.callout-jog');
          const diff = anchorPct - slotPct;
          const minPct = Math.min(slotPct, anchorPct);
          // jog is positioned absolutely inside the .callout — offset from
          // its own slot-x by (minPct - slotPct), width = |diff|.
          // Convert from container-% to a left value in px-equivalent of
          // container width by re-anchoring relative to the callout center.
          jog.style.left  = `calc(${minPct - slotPct}% + 50%)`;
          jog.style.width = `${Math.abs(diff)}%`;
          // Reframe: callout's own coord system has 0 = left edge of card,
          // but we want absolute % of the .callouts container. The card's
          // left edge sits at slot-x - 50% of card-width. Easier: re-set
          // the .callout-jog's left/width relative to the .callouts container
          // by anchoring with custom props the CSS already understands.
          jog.style.setProperty('--jog-from', `${minPct}%`);
          jog.style.setProperty('--jog-w', `${Math.abs(diff)}%`);
          container.appendChild(el);
        });
        requestAnimationFrame(() => {
          container.querySelectorAll('.callout').forEach((el, i) => {
            setTimeout(() => el.classList.add('is-shown'), i * 80);
          });
        });
      }, existing.length ? 220 : 0);
    }
    renderCallouts('top', spec.callouts.top, topC);
    renderCallouts('bot', spec.callouts.bot, botC);

    // narration update
    if (headlineEl) {
      headlineEl.textContent = spec.headline || '';
      headlineEl.classList.remove('is-shown');
      void headlineEl.offsetWidth;
      headlineEl.classList.add('is-shown');
    }
    if (narrateEl) {
      narrateEl.textContent = spec.narrate || '';
      narrateEl.classList.remove('is-shown');
      void narrateEl.offsetWidth;
      narrateEl.classList.add('is-shown');
    }

    // drift indicator: dashed link between two top callouts
    const driftLine = document.getElementById('drift-line');
    if (driftLine) {
      if (spec.drift) {
        const fromPct = pct(spec.drift.from, spec.length);
        const toPct   = pct(spec.drift.to,   spec.length);
        driftLine.style.setProperty('--from-x', `${fromPct}%`);
        driftLine.style.setProperty('--to-x',   `${toPct}%`);
        driftLine.classList.add('is-shown');
      } else {
        driftLine.classList.remove('is-shown');
      }
    }
  }

  function show(i) {
    active = ((i % SPECIMENS.length) + SPECIMENS.length) % SPECIMENS.length;
    renderSpecimen(SPECIMENS[active]);
  }

  function startCycle() {
    stopCycle();
    timer = setInterval(() => show(active + 1), 7000);
  }
  function stopCycle() {
    if (timer) { clearInterval(timer); timer = null; }
  }

  prevBtn.addEventListener('click', () => { show(active - 1); startCycle(); });
  nextBtn.addEventListener('click', () => { show(active + 1); startCycle(); });
  regenBtn.addEventListener('click', () => { show(active); startCycle(); });

  // pause cycle on hover for inspection
  root.addEventListener('mouseenter', stopCycle);
  root.addEventListener('mouseleave', startCycle);

  show(0);
  startCycle();
}
