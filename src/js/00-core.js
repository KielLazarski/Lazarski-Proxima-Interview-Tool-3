/* GlueForge — namespace, data access and small shared helpers.
 *
 * Every module hangs off one global, GF, and they run in the order given by
 * src/manifest.json. The built file inlines the cache as GF.DATA, so the
 * dashboard opens with no server and no network.
 */

(function (GF) {
  'use strict';

  /* ---- numbers -------------------------------------------------------- */

  const fmt = {
    /** Signed fixed-point, for values where direction is the point. */
    signed(v, dp) {
      if (v === null || v === undefined || Number.isNaN(v)) return '—';
      const s = v.toFixed(dp === undefined ? 3 : dp);
      return v >= 0 ? '+' + s : s;
    },
    fixed(v, dp) {
      if (v === null || v === undefined || Number.isNaN(v)) return '—';
      return v.toFixed(dp === undefined ? 3 : dp);
    },
    /** p-values: small ones as an exponent rather than a row of zeros. */
    p(v) {
      if (v === null || v === undefined || Number.isNaN(v)) return '—';
      if (v < 1e-4) return v.toExponential(1);
      return v.toFixed(4);
    },
    int(v) {
      return v === null || v === undefined ? '—' : String(v);
    },
    bytes(n) {
      return n < 1024 ? n + ' B'
        : n < 1048576 ? (n / 1024).toFixed(0) + ' kB'
        : (n / 1048576).toFixed(1) + ' MB';
    }
  };

  /* ---- statistics ------------------------------------------------------ */

  const stat = {
    mean(a) {
      if (!a.length) return NaN;
      let s = 0;
      for (let i = 0; i < a.length; i++) s += a[i];
      return s / a.length;
    },
    median(a) {
      if (!a.length) return NaN;
      const b = Array.from(a).sort((x, y) => x - y);
      const m = b.length >> 1;
      return b.length % 2 ? b[m] : 0.5 * (b[m - 1] + b[m]);
    },
    sd(a) {
      const n = a.length;
      if (n < 2) return NaN;
      const mu = stat.mean(a);
      let acc = 0;
      for (let i = 0; i < n; i++) { const d = a[i] - mu; acc += d * d; }
      return Math.sqrt(acc / (n - 1));
    },
    /** Histogram counts over a fixed range, for the distribution plot. */
    histogram(values, lo, hi, bins) {
      const counts = new Int32Array(bins);
      const w = (hi - lo) / bins;
      for (let i = 0; i < values.length; i++) {
        let k = Math.floor((values[i] - lo) / w);
        if (k < 0) k = 0;
        if (k >= bins) k = bins - 1;
        counts[k]++;
      }
      return { counts, lo, hi, width: w };
    },
    /** Two-sided sign test, used for the across-system meta-analysis. */
    signTest(values) {
      const n = values.filter((v) => v !== 0 && Number.isFinite(v)).length;
      const k = values.filter((v) => v < 0).length;
      let p = 0;
      for (let x = k; x <= n; x++) p += binomPmf(x, n, 0.5);
      return { n, negative: k, pGreater: Math.min(1, p) };
    }
  };

  function binomPmf(k, n, p) {
    return Math.exp(
      lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1) +
      k * Math.log(p) + (n - k) * Math.log(1 - p)
    );
  }

  /** Lanczos approximation; plenty for binomial terms at these sizes. */
  function lgamma(z) {
    const g = [
      676.5203681218851, -1259.1392167224028, 771.32342877765313,
      -176.61502916214059, 12.507343278686905, -0.13857109526572012,
      9.9843695780195716e-6, 1.5056327351493116e-7
    ];
    if (z < 0.5) {
      return Math.log(Math.PI / Math.sin(Math.PI * z)) - lgamma(1 - z);
    }
    z -= 1;
    let x = 0.99999999999980993;
    for (let i = 0; i < g.length; i++) x += g[i] / (z + i + 1);
    const t = z + g.length - 0.5;
    return 0.5 * Math.log(2 * Math.PI) + (z + 0.5) * Math.log(t) - t + Math.log(x);
  }

  /* ---- colour --------------------------------------------------------- */

  /** Diverging scale for frustration: warm where frustrated, cool where
   *  minimally frustrated, near-neutral in between. Mid-point is 0, and the
   *  conventional cut-offs at -1 and +0.78 sit well inside the ramp. */
  function frustrationColour(f, alpha) {
    const t = Math.max(-1, Math.min(1, f / 2.5));
    const a = alpha === undefined ? 1 : alpha;
    if (t < 0) {
      const u = -t;
      const r = Math.round(118 + u * 82);
      const g = Math.round(116 - u * 52);
      const b = Math.round(112 - u * 58);
      return `rgba(${r},${g},${b},${a})`;
    }
    const u = t;
    const r = Math.round(118 - u * 68);
    const g = Math.round(116 + u * 12);
    const b = Math.round(112 + u * 76);
    return `rgba(${r},${g},${b},${a})`;
  }

  /* ---- data access ---------------------------------------------------- */

  const data = {
    index() { return GF.DATA && GF.DATA.index; },
    potential() { return GF.DATA && GF.DATA.potential; },
    system(entry) { return GF.DATA && GF.DATA.systems && GF.DATA.systems[entry]; },
    entries() {
      const ix = data.index();
      return ix ? ix.systems.map((s) => s.entry) : [];
    },
    /** Systems grouped by biological pair, so replicates are never shown as
     *  independent evidence without saying so. */
    groups() {
      const ix = data.index();
      if (!ix) return [];
      const by = new Map();
      ix.systems.forEach((s) => {
        const key = s.group + (s.form ? ' · ' + s.form : '');
        if (!by.has(key)) by.set(key, { label: key, group: s.group, form: s.form, rows: [] });
        by.get(key).rows.push(s);
      });
      return Array.from(by.values());
    }
  };

  /* ---- dom ------------------------------------------------------------ */

  function el(tag, attrs, children) {
    const node = document.createElement(tag);
    if (attrs) {
      Object.keys(attrs).forEach((k) => {
        if (k === 'class') node.className = attrs[k];
        else if (k === 'text') node.textContent = attrs[k];
        else if (k === 'html') node.innerHTML = attrs[k];
        else if (k.startsWith('on') && typeof attrs[k] === 'function') {
          node.addEventListener(k.slice(2).toLowerCase(), attrs[k]);
        } else if (attrs[k] !== null && attrs[k] !== undefined) {
          node.setAttribute(k, attrs[k]);
        }
      });
    }
    (children || []).forEach((c) => {
      if (c === null || c === undefined) return;
      node.appendChild(typeof c === 'string' ? document.createTextNode(c) : c);
    });
    return node;
  }

  /** Give every wide table its own horizontal scroll.
   *
   *  A seven-column table of measurements has no honest narrow layout: you
   *  either scroll it or you drop columns, and dropping columns from a
   *  results table is how numbers get quoted without their controls. So it
   *  scrolls. Done here, once, rather than at each of the dozen places a
   *  table is built.
   */
  function wrapTables(root) {
    if (!root || !root.querySelectorAll) return root;
    root.querySelectorAll('table.grid').forEach((t) => {
      const parent = t.parentNode;
      if (!parent || (parent.classList && parent.classList.contains('tscroll'))) return;
      const box = el('div', { class: 'tscroll' });
      parent.insertBefore(box, t);
      box.appendChild(t);
    });
    return root;
  }

  GF.fmt = fmt;
  GF.stat = stat;
  GF.data = data;
  GF.el = el;
  GF.wrapTables = wrapTables;
  GF.frustrationColour = frustrationColour;
})(globalThis.GF || (globalThis.GF = {}));
