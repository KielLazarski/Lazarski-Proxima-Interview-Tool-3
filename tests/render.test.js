/* Render smoke test.
 *
 * There is no browser in this environment, so this drives the page builders
 * against a DOM shim small enough to be obviously correct. It will not catch
 * a layout problem, but it does catch the things that actually break a build:
 * a page builder throwing, a field renamed in the cache but not in the page,
 * a route producing nothing.
 *
 * Run: node tests/render.test.js
 */

'use strict';

const fs = require('fs');
const path = require('path');
const vm = require('vm');

const ROOT = path.resolve(__dirname, '..');
const CACHE = path.join(ROOT, 'cache', 'web');

/* ---- a very small DOM ------------------------------------------------- */

function makeNode(tag) {
  return {
    tagName: String(tag).toUpperCase(),
    nodeType: 1,
    children: [],
    childNodes: [],
    attributes: {},
    style: {},
    classList: {
      _set: new Set(),
      add(c) { this._set.add(c); },
      remove(c) { this._set.delete(c); },
      contains(c) { return this._set.has(c); }
    },
    _text: '',
    _html: '',
    dataset: {},
    get firstChild() { return this.childNodes[0] || null; },
    get lastChild() { return this.childNodes[this.childNodes.length - 1] || null; },
    getBoundingClientRect() { return { left: 0, top: 0, width: 600, height: 500 }; },
    set className(v) { this.attributes.class = v; },
    get className() { return this.attributes.class || ''; },
    set textContent(v) { this._text = String(v); },
    get textContent() {
      if (this._text) return this._text;
      return this.childNodes.map((c) => c.textContent || '').join('');
    },
    set innerHTML(v) { this._html = v; this.children = []; this.childNodes = []; },
    get innerHTML() { return this._html; },
    setAttribute(k, v) { this.attributes[k] = String(v); },
    getAttribute(k) { return this.attributes[k]; },
    appendChild(c) {
      this.childNodes.push(c);
      if (c.nodeType === 1) this.children.push(c);
      c.parentNode = this;
      return c;
    },
    // Handlers are kept rather than discarded so a test can fire one. A
    // button is two things -- a label and what it does -- and a fake DOM
    // that drops the second can only ever check the first. That is how
    // "Assess this pose" shipped as a button that merely navigated.
    addEventListener(ev, fn) { (this._on || (this._on = {}))[ev] = fn; },
    removeEventListener(ev) { if (this._on) delete this._on[ev]; },
    click() { if (this._on && this._on.click) this._on.click({
      preventDefault() {}, stopPropagation() {} }); },
    getContext() { return makeCtx(); },
    get clientWidth() { return 600; },
    get clientHeight() { return 500; },
    setPointerCapture() {},
    releasePointerCapture() {}
  };
}

function makeCtx() {
  const noop = () => {};
  return {
    setTransform: noop, clearRect: noop, beginPath: noop, moveTo: noop,
    lineTo: noop, stroke: noop, arc: noop, fill: noop, fillText: noop,
    save: noop, restore: noop,
    set strokeStyle(v) {}, set fillStyle(v) {}, set lineWidth(v) {},
    set font(v) {}
  };
}

function makeDocument() {
  const doc = {
    body: makeNode('body'),
    readyState: 'complete',
    createElement: (t) => makeNode(t),
    createElementNS: (ns, t) => makeNode(t),
    createTextNode: (t) => ({ nodeType: 3, textContent: String(t) }),
    addEventListener() {},
    _byId: {},
    getElementById(id) {
      if (!this._byId[id]) this._byId[id] = makeNode('div');
      return this._byId[id];
    }
  };
  return doc;
}

/* ---- harness ---------------------------------------------------------- */

function loadApp(data) {
  const doc = makeDocument();
  const sandbox = {
    Math, JSON, console, Date, Map, Set, Number, String, Array, Object, Error,
    Int32Array, Float64Array, Float32Array, Uint8Array, Int8Array,
    Uint16Array, Uint32Array, isNaN, parseFloat, parseInt, Infinity,
    document: doc,
    performance: { now: () => Date.now() },
    requestAnimationFrame: (fn) => fn(),
    setTimeout: (fn) => fn(),
    location: { hash: '#/benchmark' },
    devicePixelRatio: 1
  };
  sandbox.window = sandbox;
  sandbox.globalThis = sandbox;
  sandbox.window.addEventListener = () => {};
  sandbox.window.scrollTo = () => {};

  vm.createContext(sandbox);

  const manifest = JSON.parse(
    fs.readFileSync(path.join(ROOT, 'src', 'manifest.json'), 'utf8')
  );
  // Everything except startup, which boots on load.
  manifest.modules
    .filter((m) => !m.includes('50-startup'))
    .forEach((rel) => {
      const src = fs.readFileSync(path.join(ROOT, 'src', rel), 'utf8');
      vm.runInContext(src, sandbox, { filename: rel });
    });

  sandbox.GF.DATA = data;
  return sandbox;
}

function loadCache() {
  const index = JSON.parse(fs.readFileSync(path.join(CACHE, 'index.json'), 'utf8'));
  const systems = {};
  index.systems.forEach((s) => {
    systems[s.entry] = JSON.parse(
      fs.readFileSync(path.join(CACHE, s.entry + '.json'), 'utf8')
    );
  });
  const potential = JSON.parse(
    fs.readFileSync(path.join(CACHE, 'crossval.json'), 'utf8')
  ).potential;
  return { index, systems, potential };
}

let failures = 0;
function findByText(node, re) {
  if (!node || node.nodeType !== 1) return null;
  for (const c of node.childNodes) {
    const hit = findByText(c, re);
    if (hit) return hit;
  }
  return re.test(node.textContent || '') ? node : null;
}

function check(label, ok, detail) {
  if (!ok) failures++;
  console.log(`  [${ok ? 'PASS' : 'FAIL'}] ${label}${detail ? '  ' + detail : ''}`);
}

function textOf(node) {
  return node.textContent || '';
}

function main() {
  if (!fs.existsSync(path.join(CACHE, 'index.json'))) {
    console.error('no cache; run python3 sidecar/export_web.py first');
    process.exit(2);
  }

  const data = loadCache();
  const app = loadApp(data);
  const GF = app.GF;

  console.log(`render check: ${Object.keys(data.systems).length} systems`);

  // --- benchmark -------------------------------------------------------
  let bench = null;
  try { bench = GF.pages.benchmark(); } catch (e) { console.log('    ' + e.stack); }
  check('benchmark page builds', bench !== null);

  if (bench) {
    const txt = textOf(bench);
    check('headline question present', txt.includes('frustration mark'));
    check(
      'every entry appears',
      data.index.systems.every((s) => txt.includes(s.entry)),
      data.index.systems.length + ' entries'
    );
    check('warns about the stand-in potential',
      data.index.potential_published || txt.includes('stand-in'));
    check('reports a sign test', txt.includes('sign test'));
  }

  // --- each system ------------------------------------------------------
  let built = 0, firstErr = null;
  Object.keys(data.systems).forEach((entry) => {
    try {
      const page = GF.pages.system(entry);
      const txt = textOf(page);
      if (txt.includes(entry) && txt.length > 200) built++;
    } catch (e) {
      if (!firstErr) firstErr = entry + ': ' + e.message;
    }
  });
  check(
    'every system page builds',
    built === Object.keys(data.systems).length,
    built + '/' + Object.keys(data.systems).length +
      (firstErr ? '  first error: ' + firstErr : '')
  );

  // --- method -----------------------------------------------------------
  let method = null;
  try { method = GF.pages.method(); } catch (e) { console.log('    ' + e.stack); }
  check('method page builds', method !== null);
  if (method) {
    const txt = textOf(method);
    check('provenance includes the potential', txt.includes(data.index.potential));
    check('states the oracle caveat', txt.includes('already contains the answer'));
  }

  // --- the maths the pages rely on --------------------------------------
  const diffs = data.index.systems.map((s) => s.result.difference);
  const sign = GF.stat.signTest(diffs);
  check(
    'sign test is sane',
    sign.n === diffs.length && sign.pGreater >= 0 && sign.pGreater <= 1,
    `${sign.negative}/${sign.n}, p=${sign.pGreater.toExponential(2)}`
  );

  const c = GF.frustrationColour(-2, 1);
  check('colour scale returns rgba', /^rgba\(/.test(c), c);

  // --- glossary ----------------------------------------------------------
  const terms = Object.keys(GF.help.TERMS);
  check('glossary is populated', terms.length >= 12, terms.length + ' terms');
  const complete = terms.filter((k) => {
    const t = GF.help.TERMS[k];
    return t.title && t.eli5 && t.concise && t.detail;
  });
  check(
    'every term has all three depths',
    complete.length === terms.length,
    complete.length + '/' + terms.length
  );
  const distinct = terms.filter((k) => {
    const t = GF.help.TERMS[k];
    return t.eli5 !== t.concise && t.concise !== t.detail;
  });
  check('the three depths differ', distinct.length === terms.length);

  // Terms referenced by the pages must exist, or a "?" silently does nothing.
  const pagesSrc = fs.readFileSync(
    path.join(ROOT, 'src', 'js', '05-pages.js'), 'utf8'
  );
  const referenced = new Set();
  const re = /help\.attach\([^,]+,\s*'([a-z0-9_]+)'/g;
  let m;
  while ((m = re.exec(pagesSrc)) !== null) referenced.add(m[1]);
  const unknown = Array.from(referenced).filter((k) => !GF.help.TERMS[k]);
  check(
    'every attached term exists in the glossary',
    unknown.length === 0,
    referenced.size + ' referenced' + (unknown.length ? '; missing: ' + unknown.join(', ') : '')
  );

  // --- viewer ------------------------------------------------------------
  const firstEntry = data.index.systems[0].entry;
  const vsys = data.systems[firstEntry];
  const canvas = makeNode('canvas');
  const viewer = new GF.Viewer(canvas, vsys, function () {});
  viewer.draw();
  const legend = viewer.legend();
  check('viewer legend names both proteins', legend.length >= 5 &&
    legend[0].label === vsys.assignment.anchor_name &&
    legend[1].label === vsys.assignment.partner_name,
    legend.length + ' rows');
  check('viewer legend gives the chains different colours',
    legend[0].swatch !== legend[1].swatch);
  const described = viewer.describe(0);
  check('contact description is complete',
    !!(described.a && described.b && described.verdict &&
       Number.isFinite(described.frustration)),
    described.a + ' — ' + described.b + ' · ' + described.verdict);
  let ifaceNodes = 0;
  for (let i = 0; i < viewer._interfaceResidue.length; i++) {
    if (viewer._interfaceResidue[i]) ifaceNodes++;
  }
  check('interface residues are marked for drawing', ifaceNodes > 0,
    ifaceNodes + ' residues');

  // --- the viewer on a run blob, which has no frustration behind it ------
  // Docking results reach this viewer too, and they are coordinates and
  // nothing else. Every style that reads contact energies or exit vectors
  // has to be absent rather than broken.
  {
    const full = data.systems[firstEntry];
    const lean = JSON.parse(JSON.stringify(full));
    delete lean.contacts;
    delete lean.exit_vectors;
    delete lean.result;
    delete lean.residues.rho;
    lean.drawing_only = true;
    lean.missing = ['contacts', 'exit_vectors', 'result', 'residues.rho'];

    let v = null, threw = null;
    try {
      v = new GF.Viewer(makeNode('canvas'), lean, function () {});
      v.draw();
    } catch (e) { threw = e; }
    check('viewer builds a system with no frustration data',
      !threw, threw ? threw.message : '');

    if (v) {
      const ids = v.styles().map((s) => s.id);
      check('it offers only the styles it can draw',
        ids.indexOf('contacts') < 0 && ids.indexOf('vectors') < 0 &&
        ids.indexOf('ribbon') >= 0 && ids.indexOf('sidechains') >= 0,
        ids.join(', '));
      check('it does not open on a style it cannot draw',
        v.style !== 'contacts' && v.style !== 'vectors', v.style);

      // The interface has to come from the blob, since there are no
      // contacts to derive it from, and it must not come back empty.
      let marked = 0;
      for (let i = 0; i < v._interfaceResidue.length; i++) {
        if (v._interfaceResidue[i]) marked++;
      }
      check('it finds the interface from the blob instead', marked > 0,
        marked + ' residues');

      // And every drawable style must survive being selected.
      let bad = null;
      v.styles().forEach((st) => {
        try { v.style = st.id; v.draw(); } catch (e) { bad = st.id + ': ' + e.message; }
      });
      check('every offered style draws on a run blob', !bad, bad || '');

      // Selecting a style it cannot draw must be inert, not fatal: the
      // style is set from a saved preference in some paths.
      let crashed = null;
      try { v.style = 'contacts'; v.draw(); v.style = 'vectors'; v.draw(); }
      catch (e) { crashed = e.message; }
      check('an unavailable style is inert rather than fatal',
        !crashed, crashed || '');
    }
  }

  // --- the landing page leads with the work ------------------------------
  // It used to open on a wall of install buttons, which is the wrong first
  // thing when, most of the time, nothing needs installing.
  if (GF.pages.landing && GF.pages.setup) {
    const env = (missingPkgs) => ({
      python: { version: '3.13.5', platform: 'Windows', machine: 'AMD64' },
      packages: [
        { key: 'numpy', present: true, required: true, version: '2.1.0',
          stage: 'everything', why: 'geometry' },
        { key: 'scipy', present: true, required: true, version: '1.14.0',
          stage: 'everything', why: 'stats' },
        { key: 'rdkit', present: missingPkgs.indexOf('rdkit') < 0,
          required: false, stage: 'chemistry', why: 'molecules' }
      ],
      stages: [
        { key: 'assess', label: 'Assess a complex', ready: true, missing: [],
          better_with_missing: [], note: 'n' },
        { key: 'dock', label: 'Dock a new pair', ready: true, missing: [],
          better_with_missing: ['lightdock'], note: 'n' },
        { key: 'chemistry', label: 'Grow chemistry',
          ready: missingPkgs.indexOf('rdkit') < 0,
          missing: missingPkgs.indexOf('rdkit') < 0 ? [] : ['rdkit'],
          better_with_missing: [], note: 'n' }
      ],
      runs: []
    });

    const show = (missingPkgs) => {
      GF.api.state.live = true;
      GF.api.state.env = env(missingPkgs);
      let node = null, err = null;
      try { node = GF.pages.landing(); } catch (e) { err = e; }
      return { node: node, err: err, text: node ? textOf(node) : '' };
    };

    let r = show([]);
    check('landing builds', !!r.node && !r.err, r.err ? r.err.message : '');
    check('it leads with what you can do',
      r.text.indexOf('What you can do') >= 0 &&
      r.text.indexOf('What you can do') < r.text.indexOf('Setup'),
      'the stage list must come before the setup pointer');
    check('the install buttons are no longer on the front page',
      !/What is installed/.test(r.text) && !/\bInstall\b/.test(r.text));
    check('with nothing missing the pointer is quiet',
      /Everything these need is installed/.test(r.text));

    r = show(['rdkit']);
    check('a blocked stage is called out by name',
      /rdkit not installed/.test(r.text) &&
      /grow chemistry unavailable/i.test(r.text), r.text.slice(0, 120));
    check('and it does not shout about a stage that merely runs better',
      !/lightdock not installed/i.test(r.text));

    let sp = null, e2 = null;
    try { sp = GF.pages.setup(); } catch (e) { e2 = e; }
    check('setup page builds', !!sp && !e2, e2 ? e2.message : '');
    if (sp) {
      const t = textOf(sp);
      check('setup holds the packages', /Packages/.test(t) && /rdkit/.test(t));
      check('setup names the Python it installs into',
        /3\.13\.5/.test(t) && /installs into whichever Python/.test(t));
    }
    GF.api.state.live = false;
    GF.api.state.env = null;
  }

  // --- a failed install must say what would fix it ----------------------
  // The message people see when pip fails is the whole of what they have to
  // work from, so it is worth a test. This is the verbatim output from a
  // real failure on Windows with Python 3.13.
  const PIP_FAIL = [
    '  note: This error originates from a subprocess, and is likely not a problem with pip.',
    '  ERROR: Failed building wheel for prody',
    'Failed to build lightdock freesasa prody',
    '',
    'ERROR: Failed to build installable wheels for some pyproject.toml based projects (lightdock, freesasa, prody)'
  ].join('\n');

  if (GF.landing && GF.landing.diagnose) {
    const onWindows = (fn) => {
      GF.api.state.env = {
        python: { version: '3.13.5', platform: 'Windows', machine: 'AMD64' }
      };
      try { return fn(); } finally { GF.api.state.env = null; }
    };

    const node = GF.landing.diagnose(PIP_FAIL);
    const said = node ? textOf(node) : '';
    check('a failed install is diagnosed', !!node);
    check('it names the Python as the cause, version or not',
      /Python 3\.\d+|this Python/.test(said) && !/Python \./.test(said),
      said.slice(0, 54));
    check('it names the packages that failed',
      said.includes('prody') && said.includes('freesasa'));
    check('it does not suggest re-running what just failed',
      !/run the same pip command/i.test(said));

    // The case the advice was rewritten for. LightDock has no Windows build
    // at any Python version and its own page does not claim Windows support,
    // so an older Python is not a fix and saying so would waste an hour.
    const win = onWindows(() => textOf(GF.landing.diagnose(PIP_FAIL)));
    check('on Windows it reads the running Python', /Python 3\.13/.test(win));
    check('it sends an unbuildable lightdock to WSL, not to an older Python',
      /wsl --install/.test(win) && /build-essential/.test(win) &&
      !/py -3\.\d+ -m venv/.test(win));
    check('the WSL path says how the page is reached from Windows',
      /localhost/.test(win));
    check('it says which stages still work without lightdock',
      /Docking is the only stage/.test(win) && /RDKit/.test(win));

    // A failure with no lightdock in it is the ordinary version mismatch,
    // and there the download probe is the right next step -- one line per
    // package, because pip stops at the first failure.
    const SOME = 'ERROR: Failed to build installable wheels for some ' +
      'pyproject.toml based projects (freesasa, prody)';
    const other = onWindows(() => textOf(GF.landing.diagnose(SOME)));
    check('a plain version mismatch gets the download probe',
      /--only-binary/.test(other));
    check('the probe is one line per package',
      (other.match(/--only-binary/g) || []).length === 2);
    check('and an older Python that has a build',
      /py -3\.11 -m venv/.test(other));

    check('an unrecognised failure gets no invented explanation',
      GF.landing.diagnose('something else entirely') === null);
    const perm = textOf(GF.landing.diagnose('ERROR: Permission denied') || makeNode('p'));
    check('a permissions failure is told apart from a missing wheel',
      /venv/.test(perm) && !/Nothing prebuilt/.test(perm));
  }

  // --- the docking screen, in every engine state ------------------------
  // Added after this screen threw in front of a user. The whole page is
  // decided by what /api/presets said, so building it without that payload
  // tests the word "Loading" and nothing else.
  if (GF.pages.dock && GF.dock) {
    const ldPreset = (key) => ({
      key: key, label: key, swarms: 400, glowworms: 100, steps: 60,
      minutes_low: 20, minutes_high: 80, note: 'a note', measured: null
    });
    const fftPreset = (key) => ({
      key: key, label: key, engine: 'builtin', rotations: 3000,
      spacing: 1.6, keep: 200, focus: [80, 30], minutes_estimate: [4, 11]
    });
    const keys = ['quick', 'standard', 'thorough'];
    const payload = (lightdockHere, measured) => ({
      presets: keys.map(ldPreset),
      default: 'standard',
      calibrated: false,
      engine_default: lightdockHere ? 'lightdock' : 'builtin',
      engines: {
        lightdock: { available: lightdockHere, presets: keys.map(ldPreset) },
        builtin: { available: true, presets: keys.map(fftPreset),
                   measured: measured || null }
      }
    });

    let built = null;
    try { built = GF.pages.dock('run-under-test'); }
    catch (e) { console.log('    ' + e.stack); }
    check('dock screen builds before the presets arrive', !!built);

    const draw = (state) => {
      Object.keys(state).forEach((k) => { GF.dock.state[k] = state[k]; });
      try { GF.dock.render(); return null; } catch (e) { return e; }
    };

    // The case that broke: LightDock absent, so there is one engine and the
    // picker is null.
    let err = draw({ presets: payload(false, null), engine: 'builtin',
                     chosen: 'quick' });
    check('dock screen builds with only the built-in engine',
      !err, err ? err.message : '');
    let said = built ? textOf(built) : '';
    check('it does not offer a choice there is no choice about',
      !/LightDock/.test(said) || !/Built-in/.test(said));
    check('it admits the recovery rate is unmeasured',
      /not measured/i.test(said));

    // With a measurement, it must carry the caveat that makes the number
    // readable rather than just the number.
    err = draw({ presets: payload(false, {
      recovered: 1, designable: 0, n: 6, cutoff: 5, preset: 'quick',
      bound_bound: true, largest_partner: 2423 }) });
    said = built ? textOf(built) : '';
    check('a measured engine shows its rate', !err && /1 of 6/.test(said),
      err ? err.message : '');
    check('and says the benchmark was bound-bound', /bound shape/.test(said));

    // Both engines present: now the picker should appear.
    err = draw({ presets: payload(true, null), engine: 'lightdock',
                 chosen: 'standard' });
    said = built ? textOf(built) : '';
    check('with both engines it offers the choice',
      !err && /LightDock/.test(said) && /Built-in/.test(said),
      err ? err.message : '');
    // The two engines have never been run against each other here, so the
    // screen must describe them rather than rank them. "Stronger" and
    // "weaker" shipped on these buttons for a while, standing in for a
    // comparison nobody had made -- the one kind of claim this tool exists
    // to not make. If a comparison is ever measured, say the number here
    // instead of reaching for the adjective again.
    check('but does not call either engine stronger without a comparison',
      !/stronger|weaker|better of the two/i.test(said),
      said.slice(said.search(/stronger|weaker|better of the two/i) - 60, 120));

    // And the poses table, which is where the ranking change landed.
    err = draw({
      presets: payload(false, null), engine: 'builtin', chosen: 'quick',
      state: 'done',
      result: {
        engine: 'builtin',
        poses: [
          { label: 'fft-0000', face_burial: 0.2, contacts: 90, clashes: 1,
            ligand_clashes: 0, restraint_satisfaction: 0.5, viable: true },
          { label: 'fft-0001', face_burial: 0.9, contacts: 70, clashes: 2,
            ligand_clashes: 0, restraint_satisfaction: 0.4, viable: true }
        ]
      }
    });
    said = built ? textOf(built) : '';
    check('the poses table builds', !err, err ? err.message : '');
    check('the built-in engine keeps its own order, not face burial',
      said.indexOf('fft-0000') < said.indexOf('fft-0001'),
      'the pose with less face buried must stay first');

    // The button has to run the gate, not just change the hash. It shipped
    // doing only the second, so the next screen said nothing had been
    // assessed -- true, and a button that names what it does not do.
    const btn = findByText(built, /^Assess this pose$/);
    check('the dock screen offers to assess a pose', !!btn);
    if (btn) {
      const calls = [];
      const realJson = GF.api.json;
      GF.api.json = (url, body) => {
        calls.push([url, body]);
        return { then: () => ({ catch: () => {} }) };
      };
      try { btn.click(); } finally { GF.api.json = realJson; }
      check('clicking it asks the server to assess, not just navigates',
        calls.length === 1, 'calls: ' + JSON.stringify(calls));
      check('and it says which pose',
        calls.length === 1 && /\/assess$/.test(calls[0][0])
        && calls[0][1] && calls[0][1].pose === 'fft-0000',
        JSON.stringify(calls[0] || null));
    }
  }

  // --- the assessment panel ---------------------------------------------
  // The thresholds behind every verdict were set on crystal structures,
  // and the engine lands 6-8 A out while a glue cavity tolerates ~2, so a
  // verdict on a docked pose is a different claim from a verdict on a
  // solved one. The panel has to say which it is.
  if (GF.assess && GF.assess.panel) {
    const assessment = (pose) => {
      const a = JSON.parse(JSON.stringify(
        data.systems[firstEntry].result || {}));
      a.verdict = a.verdict || { verdict: 'borderline', headline: 'Mixed',
                                 n_checks: 4, n_pass: 2, n_fail: 1 };
      if (pose) a.pose = pose;
      return a;
    };

    let built = null, perr = null;
    try { built = GF.assess.panel(assessment(null), null, {}); }
    catch (e) { perr = e; }
    check('the assessment panel builds for a solved complex',
      !perr && !!built, perr ? perr.stack : '');
    let said = built ? textOf(built) : '';
    check('and claims nothing about docking when none was done',
      !/Judged on a docked pose/.test(said));

    try { built = GF.assess.panel(assessment('fft-0000'), null, {}); perr = null; }
    catch (e) { perr = e; built = null; }
    said = built ? textOf(built) : '';
    check('a verdict on a docked pose says so',
      !perr && /Judged on a docked pose/.test(said), perr ? perr.stack : '');
    check('and names the pose it judged', /fft-0000/.test(said));
    check('and warns the verdict may be about the docking',
      /where the search put protein B/.test(said));
    // The warning has to carry the measurement, not an adjective. A
    // caveat with no number in it gets read as boilerplate and skipped.
    check('and gives the measured miss distance, not a vague caveat',
      /2\.6/.test(said) && /20/.test(said), said.slice(0, 200));
  }

  // --- the run viewer panels --------------------------------------------
  // These are what a docking result is actually read through, so they get
  // built here against a run-shaped blob rather than trusted.
  if (GF.view) {
    const runBlob = JSON.parse(JSON.stringify(data.systems[firstEntry]));
    delete runBlob.contacts;
    delete runBlob.exit_vectors;
    delete runBlob.result;
    delete runBlob.residues.rho;
    runBlob.drawing_only = true;
    runBlob.missing = ['contacts', 'exit_vectors', 'result', 'residues.rho'];
    runBlob.label = 'fft-0000';

    let sp = null, e1 = null;
    try { sp = GF.view.panel(runBlob, { title: 'The placement' }); }
    catch (e) { e1 = e; }
    check('structure panel builds for a run', !sp === false && !e1,
      e1 ? e1.message : '');
    if (sp) {
      const said = textOf(sp.node);
      check('it names the pose and both proteins',
        /fft-0000/.test(said) && said.length > 40, said.slice(0, 60));
      check('it says why two styles are missing',
        /not computed for a docking result/.test(said));
      // Swapping pose must not rebuild the viewer from nothing.
      let e2 = null;
      try {
        const other = JSON.parse(JSON.stringify(runBlob));
        other.label = 'fft-0007';
        sp.setSystem(other);
      } catch (e) { e2 = e; }
      check('switching pose redraws without throwing', !e2,
        e2 ? e2.message : '');
      check('and the new pose is named', /fft-0007/.test(textOf(sp.node)));
    }

    let cp = null, e3 = null;
    try { cp = GF.view.cavityPanel(runBlob); } catch (e) { e3 = e; }
    check('cavity panel builds for a run', !!cp && !e3, e3 ? e3.message : '');
    if (cp) {
      check('it reports the pocket it found',
        /Å³|no cavity/i.test(textOf(cp.node)), textOf(cp.node).slice(0, 70));
    }

    // A pose with no pocket beside the ligand is a real outcome and must
    // draw as one rather than throwing on an absent cavity.
    const noCav = JSON.parse(JSON.stringify(runBlob));
    noCav.cavity = { points: 0, note: 'no cavity detected at this site' };
    let e4 = null, cp2 = null;
    try { cp2 = GF.view.cavityPanel(noCav); } catch (e) { e4 = e; }
    check('an empty cavity is a result, not a crash', !!cp2 && !e4,
      e4 ? e4.message : '');
    if (cp2) {
      check('and it says so in those terms',
        /nowhere for added chemistry/.test(textOf(cp2.node)));
    }
  }

  // --- a page newer than the server answering it -------------------------
  // The situation after unzipping a new version over a running one: new
  // pages, old routes, symptoms everywhere and none of them naming the
  // cause.
  if (GF.api && GF.api.serverIsStale) {
    const built = app.GF.BUILT;
    app.GF.BUILT = 2000;

    GF.api.state.env = { started: 1000 };
    check('a server older than the page is detected',
      GF.api.serverIsStale() === true);

    GF.api.state.env = { started: 3000 };
    check('a server newer than the page is not', !GF.api.serverIsStale());

    // The normal case: built moments before the server started, or the
    // other way round. Neither is a mismatch worth a banner.
    GF.api.state.env = { started: 1990 };
    check('a few seconds either way is not a mismatch',
      !GF.api.serverIsStale(), 'built 2000, started 1990');

    GF.api.state.env = {};
    check('with nothing to compare it says no', !GF.api.serverIsStale());
    app.GF.BUILT = built;
    GF.api.state.env = null;
  }

  // --- the built artefacts -----------------------------------------------
  // Two builds with opposite jobs. The served page is handed out by
  // glueforge.py and fetches the validation set on demand, so inlining it
  // would make every start pay four megabytes nobody asked for. The
  // standalone page is the export: it has to open on a machine with no
  // Python and no network, so the data must be inside it. Checking only one
  // of them lets the other regress into the wrong shape unnoticed.
  const served = path.join(ROOT, 'dist', 'glueforge.html');
  const alone = path.join(ROOT, 'dist', 'glueforge-standalone.html');

  const noMarkers = (h) =>
    !h.includes('__STYLES__') && !h.includes('__SCRIPTS__') &&
    !h.includes('__DATA__');
  const noNetwork = (h) =>
    !/<(script|link)[^>]+(src|href)=["']https?:/i.test(h);

  if (fs.existsSync(served)) {
    const html = fs.readFileSync(served, 'utf8');
    const kb = Math.round(html.length / 1024);
    check('served build has no unreplaced markers', noMarkers(html));
    check('served build leaves the cache on the server',
      html.includes('GF.DATA = null'), kb + ' kB');
    check('served build stays small', kb < 600, kb + ' kB');
    check('served build needs no network', noNetwork(html));
  }

  if (fs.existsSync(alone)) {
    const html = fs.readFileSync(alone, 'utf8');
    check('standalone build has no unreplaced markers', noMarkers(html));
    check('standalone build inlines the cache', html.includes('GF.DATA = {'));
    check('standalone build needs no network', noNetwork(html));
    const entries = Object.keys(
      JSON.parse(fs.readFileSync(path.join(CACHE, 'index.json'), 'utf8'))
    ).length ? JSON.parse(
      fs.readFileSync(path.join(CACHE, 'index.json'), 'utf8')
    ).systems.map((s) => s.entry) : [];
    const missing = entries.filter((e) => !html.includes('"' + e + '"'));
    check('standalone build carries every system', missing.length === 0,
      entries.length + ' systems' +
      (missing.length ? ', missing ' + missing.join(', ') : ''));
  }

  // --- cache integrity ---------------------------------------------------
  let badJson = [];
  fs.readdirSync(CACHE).filter((f) => f.endsWith('.json')).forEach((f) => {
    const raw = fs.readFileSync(path.join(CACHE, f), 'utf8');
    if (/[:,[]\s*(NaN|Infinity|-Infinity)\b/.test(raw)) badJson.push(f);
    try { JSON.parse(raw); } catch (e) { badJson.push(f + ' (' + e.message.slice(0, 40) + ')'); }
  });
  check('cache files are valid JSON', badJson.length === 0, badJson.join(', '));

  // --- exit vectors ------------------------------------------------------
  const withVectors = Object.keys(data.systems).filter((k) => {
    const ev = data.systems[k].exit_vectors;
    return ev && ev.vectors && ev.vectors.length;
  });
  check(
    'most systems yield growth points',
    withVectors.length >= Object.keys(data.systems).length - 1,
    withVectors.length + '/' + Object.keys(data.systems).length
  );

  const ranked = withVectors.every((k) => {
    const v = data.systems[k].exit_vectors.vectors;
    for (let i = 1; i < v.length; i++) {
      if (v[i].addressable > v[i - 1].addressable + 1e-9) return false;
    }
    return true;
  });
  check('growth points are ranked by addressable frustration', ranked);

  const exposureDiverges = withVectors.some((k) => {
    const v = data.systems[k].exit_vectors.vectors;
    const byExposure = v.slice().sort((a, b) => b.relative_sasa - a.relative_sasa);
    return byExposure[0].atom_name !== v[0].atom_name;
  });
  check(
    'exposure and reach disagree somewhere',
    exposureDiverges,
    'the most exposed atom is not always the best vector'
  );

  const ligOk = withVectors.every((k) => {
    const l = data.systems[k].ligand;
    return l.sasa_total >= l.sasa_total_ternary - 1e-6 &&
           l.buried_by_partner >= -1e-6 && l.buried_by_partner <= 1 + 1e-6;
  });
  check('presented face shrinks in the ternary', ligOk);

  // --- drawing styles ----------------------------------------------------
  const styles = GF.Viewer.STYLES.map((s) => s.id);
  check('four drawing styles offered', styles.length === 4, styles.join(', '));
  let drew = 0;
  styles.forEach((id) => {
    try { viewer.style = id; viewer.draw(); drew++; } catch (e) {
      console.log('    ' + id + ': ' + e.message);
    }
  });
  check('every style draws without throwing', drew === styles.length, drew + '/' + styles.length);
  viewer.style = 'contacts';

  const hasCA = Object.keys(data.systems).every((k) => {
    const r = data.systems[k].residues;
    return r.cax && r.cax.length === r.x.length && r.ss && r.ss.length === r.x.length;
  });
  check('CA trace and secondary structure exported for every system', hasCA);

  // --- cavity ------------------------------------------------------------
  const withCavity = Object.keys(data.systems).filter((k) => {
    const c = data.systems[k].cavity;
    return c && c.points > 0;
  });
  check(
    'every system has a cavity',
    withCavity.length === Object.keys(data.systems).length,
    withCavity.length + '/' + Object.keys(data.systems).length
  );

  let cavBuilt = 0, cavErr = null;
  withCavity.forEach((k) => {
    try {
      const page = GF.pages.cavity(k);
      const txt = textOf(page);
      if (txt.includes(k) && txt.length > 200) cavBuilt++;
    } catch (e) {
      if (!cavErr) cavErr = k + ': ' + e.message;
    }
  });
  check('every cavity page builds', cavBuilt === withCavity.length,
    cavBuilt + '/' + withCavity.length + (cavErr ? '  first error: ' + cavErr : ''));

  const fieldsOk = withCavity.every((k) => {
    const c = data.systems[k].cavity;
    return ['electrostatic', 'donor', 'acceptor', 'hydrophobic', 'from_anchor']
      .every((f) => Array.isArray(c[f]) && c[f].length === c.points);
  });
  check('all five fields present and correctly sized', fieldsOk);

  const splitOk = withCavity.every((k) => {
    const c = data.systems[k].cavity;
    return c.anchor_lined_share > 0.3 && c.anchor_lined_share < 0.95;
  });
  check(
    'the cavity is mostly anchor-lined everywhere',
    splitOk,
    withCavity.map((k) => Math.round(data.systems[k].cavity.anchor_lined_share * 100)).join('/') + '%'
  );

  // The ligand sits inside its own cavity, which is the detector's sanity check.
  const inside = withCavity.filter((k) => {
    const c = data.systems[k].cavity, l = data.systems[k].ligand;
    let hit = 0;
    for (let i = 0; i < l.x.length; i++) {
      let best = 1e9;
      for (let j = 0; j < c.x.length; j++) {
        const d = Math.hypot(l.x[i] - c.x[j], l.y[i] - c.y[j], l.z[i] - c.z[j]);
        if (d < best) best = d;
      }
      if (best < 2.0) hit++;
    }
    return hit / l.x.length > 0.5;
  });
  check(
    'the ligand lies inside its own cavity',
    inside.length >= withCavity.length - 3,
    inside.length + '/' + withCavity.length + ' (shallow grooves may fail)'
  );

  const withSens = withCavity.filter((k) => data.systems[k].cavity.sensitivity);
  if (withSens.length) {
    const degrades = withSens.every((k) => {
      const rows = data.systems[k].cavity.sensitivity.rows;
      const near = rows.find((r) => r.shift_angstrom === 1.0);
      const far = rows.find((r) => r.shift_angstrom === 5.0);
      return !near || !far || near.jaccard >= far.jaccard;
    });
    check('cavity degrades monotonically with partner error', degrades,
      withSens.length + ' systems with a budget');
  } else {
    console.log('  . no sensitivity data yet (run run_cavity.py --sensitivity)');
  }

  const cavityStyles = ['electrostatic', 'donor', 'acceptor', 'hydrophobic', 'from_anchor'];
  const vsys2 = data.systems[withCavity[0]];
  const cv = new GF.CavityView(makeNode('canvas'), vsys2, function () {});
  let cdrew = 0;
  cavityStyles.forEach((f) => {
    try { cv.field = f; cv.draw(); cdrew++; } catch (e) {
      console.log('    ' + f + ': ' + e.message);
    }
  });
  check('cavity view draws every field', cdrew === cavityStyles.length,
    cdrew + '/' + cavityStyles.length);

  console.log(failures === 0 ? '\nall checks passed' : `\n${failures} check(s) failed`);
  process.exit(failures === 0 ? 0 : 1);
}

main();
