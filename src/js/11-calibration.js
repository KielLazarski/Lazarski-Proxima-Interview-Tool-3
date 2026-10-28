/* GlueForge — where the gate's thresholds come from.
 *
 * Every number the assessment screen judges against is reachable from here,
 * and what it shows is the raw thing: all 14 crystal poses and all 210
 * plausible wrong poses, plotted on the same axis with the cut marked.
 *
 * The reason to build this rather than quote an AUC is that a threshold you
 * cannot see the working for is a threshold nobody can argue with, and the
 * arguable ones are the only kind worth having. On this plot you can see
 * immediately that "lined by both proteins" has two almost disjoint clouds
 * while "mean enclosure" has one cloud with a line through it — which is
 * exactly why one of them carries four times the weight of the other in a
 * verdict, and why three real glue systems come out borderline.
 */

(function (GF) {
  'use strict';

  const el = GF.el;
  const api = GF.api;

  const NS = 'http://www.w3.org/2000/svg';
  const svg = (tag, attrs) => {
    const n = document.createElementNS(NS, tag);
    Object.keys(attrs || {}).forEach((k) => {
      if (attrs[k] !== null && attrs[k] !== undefined) {
        n.setAttribute(k, attrs[k]);
      }
    });
    return n;
  };

  let cache = null;

  /* ---- one metric, with every point that calibrated it ------------------ */

  function strip(key, spec, systems) {
    const W = 680, H = 108, L = 8, R = 14, T = 26, B = 24;

    const nat = systems.map((s) => s.native[key])
                       .filter(Number.isFinite);
    const dec = [];
    systems.forEach((s) => s.decoys.forEach((d) => {
      if (Number.isFinite(d[key])) dec.push(d[key]);
    }));
    if (!nat.length || !dec.length) return null;

    const all = nat.concat(dec);
    let lo = Math.min.apply(null, all);
    let hi = Math.max.apply(null, all);
    const pad = (hi - lo) * 0.05 || 1;
    lo -= pad; hi += pad;
    const x = (v) => L + (W - L - R) * (v - lo) / (hi - lo);

    const g = svg('svg', {
      viewBox: `0 0 ${W} ${H}`, class: 'strip',
      role: 'img',
      'aria-label': spec.label + ': ' + nat.length + ' crystal poses and ' +
                    dec.length + ' wrong poses'
    });

    // The shaded side is the one the gate treats as glue-like.
    const cut = x(spec.threshold);
    g.appendChild(svg('rect', {
      x: spec.direction > 0 ? cut : L,
      y: T - 12,
      width: Math.max(0, spec.direction > 0 ? (W - R - cut) : (cut - L)),
      height: H - T - B + 24,
      class: 'strip-pass'
    }));

    // Decoys first so natives draw over them.
    dec.forEach((v, i) => g.appendChild(svg('circle', {
      cx: x(v), cy: H - B + 2 + ((i * 7) % 11) * 0.45,
      r: 1.9, class: 'pt decoy'
    })));
    nat.forEach((v, i) => g.appendChild(svg('circle', {
      cx: x(v), cy: T + ((i * 5) % 9) * 0.8,
      r: 3.1, class: 'pt native'
    })));

    g.appendChild(svg('line', {
      x1: cut, x2: cut, y1: T - 14, y2: H - B + 10, class: 'cut'
    }));

    // Place the label on whichever side has room. Anchoring it by the
    // metric's direction alone clips it whenever the cut sits near an edge,
    // which is exactly where the interesting thresholds end up.
    const roomRight = (W - R) - cut;
    const labelLeft = roomRight < 70;
    const lab = svg('text', {
      x: cut + (labelLeft ? -5 : 5), y: 13,
      class: 'cutlab',
      'text-anchor': labelLeft ? 'end' : 'start'
    });
    lab.textContent = 'cut ' + fmtNum(spec.threshold);
    g.appendChild(lab);

    [[L, lo], [W - R, hi]].forEach(([px, v], i) => {
      const t = svg('text', {
        x: px, y: H - 4, class: 'axlab',
        'text-anchor': i ? 'end' : 'start'
      });
      t.textContent = fmtNum(v);
      g.appendChild(t);
    });

    return g;
  }

  function fmtNum(v) {
    if (!Number.isFinite(v)) return '—';
    if (Math.abs(v) >= 100) return v.toFixed(0);
    if (Math.abs(v) >= 1) return v.toFixed(2);
    return v.toFixed(3);
  }

  function metricBlock(key, spec, systems, weightPct) {
    const plot = strip(key, spec, systems);
    const block = el('section', { class: 'calib', id: 'm-' + key }, [
      el('div', { class: 'calib-head' }, [
        el('span', { class: 'calib-name', text: spec.label }),
        el('span', { class: 'calib-auc mono', text: 'AUC ' + spec.auc.toFixed(3) }),
        el('span', { class: 'calib-weight', text: weightPct + '% of the verdict' })
      ]),
      plot,
      el('div', { class: 'calib-meta mono', text:
        'cut at ' + fmtNum(spec.threshold) +
        '  ·  catches ' + (100 * spec.sensitivity).toFixed(0) + '% of real glues' +
        '  ·  rejects ' + (100 * spec.specificity).toFixed(0) + '% of wrong poses' +
        '  ·  ' + (spec.direction > 0 ? 'higher is better' : 'lower is better') })
    ]);

    // The interpretation, which is the point of showing the plot at all.
    let read;
    if (spec.auc >= 0.95) {
      read = 'Two almost separate clouds. This measurement decides most of ' +
             'what the gate says, and deserves to.';
    } else if (spec.auc >= 0.85) {
      read = 'Mostly separate, with a tail of wrong poses reaching into the ' +
             'real range.';
    } else {
      read = 'One cloud with a line through it. Real glue sites and wrong ' +
             'poses overlap here, so this can support a verdict but should ' +
             'never carry one.';
    }
    block.appendChild(el('p', { class: 'calib-read', text: read }));
    return block;
  }

  /* ---- the page --------------------------------------------------------- */

  function render(host, data) {
    host.innerHTML = '';
    const metrics = data.metrics || {};
    const systems = data.systems || [];

    host.appendChild(el('div', { class: 'page' }, [
      el('h1', { text: 'Where the thresholds come from' }),
      el('p', { class: 'lede' }, [
        'Every cut the assessment screen judges against, with the evidence ' +
        'that set it. Filled circles are the ' + systems.length + ' solved ' +
        'glue systems. Small open circles are the ' +
        systems.reduce((n, s) => n + s.decoys.length, 0) + ' plausible wrong ' +
        'poses of those same partners — the ones that bury the presented ' +
        'face best, which is to say the ones a docking run would rank ' +
        'highest.'
      ]),
      el('p', { class: 'sub' }, [
        'Calibrating on the solved systems alone is the easy mistake and it ' +
        'produces a gate that passes whatever it is shown: every example is ' +
        'a positive, so no threshold drawn from them can fail anything. The ' +
        'wrong poses are what make this a test.'
      ])
    ]));

    const usable = Object.entries(metrics).filter(([, m]) => m.auc >= 0.70);
    const total = usable.reduce((t, [, m]) => t + (m.auc - 0.5), 0) || 1;

    usable
      .sort((a, b) => b[1].auc - a[1].auc)
      .forEach(([key, spec]) => {
        host.appendChild(metricBlock(
          key, spec, systems,
          Math.round(100 * (spec.auc - 0.5) / total)
        ));
      });

    const weak = Object.entries(metrics).filter(([, m]) => m.auc < 0.70);
    if (weak.length) {
      host.appendChild(el('section', { class: 'panel' }, [
        el('h2', { text: 'Excluded' }),
        el('p', { class: 'sub', text:
          'Below AUC 0.70 a measurement is close enough to a coin flip that ' +
          'including it would add confidence without adding evidence.' }),
        el('ul', { class: 'caveats' }, weak.map(([k, m]) =>
          el('li', { text: m.label + ' — AUC ' + m.auc.toFixed(3) })))
      ]));
    }

    host.appendChild(el('section', { class: 'panel' }, [
      el('h2', { text: 'What this cannot tell you' }),
      el('p', { class: 'sub' }, [
        'Fourteen systems is a small set and they are not a random sample of ' +
        'anything — they are the glue complexes that have been solved, which ' +
        'selects for the ones that crystallised. The wrong poses are wrong ' +
        'placements of the right partner, so they test whether a pose is ' +
        'right and not whether a pair is glueable at all. A verdict here is ' +
        'evidence about this interface, not a prediction about the project.'
      ]),
      el('p', { class: 'fnote mono', text:
        'seed ' + data.seed + ' · ' + data.decoys_per_system +
        ' decoys kept per system' })
    ]));

    GF.wrapTables(host);
  }

  function page(anchor) {
    const host = el('div', {});
    const draw = (d) => {
      cache = d;
      render(host, d);
      if (anchor) {
        // The router appends this node after page() returns, so a scroll
        // scheduled now runs against a node that is not in the document yet.
        // Wait for it to be attached, then scroll.
        let tries = 0;
        const go = () => {
          const n = document.getElementById('m-' + anchor);
          if (n && n.isConnected) {
            n.scrollIntoView({ block: 'start', behavior: 'auto' });
            n.classList.add('flash');
          } else if (tries++ < 20) {
            requestAnimationFrame(go);
          }
        };
        requestAnimationFrame(go);
      }
    };

    if (cache) {
      draw(cache);
    } else if (GF.EXPORT && GF.EXPORT.calibration) {
      // An exported copy carries its own calibration, so every threshold
      // stays one click from its evidence with no server anywhere.
      draw(GF.EXPORT.calibration);
    } else {
      host.appendChild(el('div', { class: 'page' }, [
        el('h1', { text: 'One moment' }),
        el('p', { class: 'lede', text: 'Loading the calibration.' })
      ]));
      api.call('/api/calibration').then(draw).catch((err) => {
        host.innerHTML = '';
        host.appendChild(el('div', { class: 'page' }, [
          el('h1', { text: 'No calibration' }),
          el('p', { class: 'lede', text: err.message }),
          el('p', { class: 'note' }, [
            'Build it with ',
            el('code', { text: 'python3 sidecar/run_gate.py' }),
            '.'
          ])
        ]));
      });
    }
    return host;
  }

  GF.pages = GF.pages || {};
  GF.pages.calibration = page;
})(globalThis.GF || (globalThis.GF = {}));
