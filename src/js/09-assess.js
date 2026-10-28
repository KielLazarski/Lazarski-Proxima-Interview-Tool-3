/* GlueForge — the assessment screen.
 *
 * The gate. It answers "is this interface worth doing chemistry into" with no
 * chemistry done, and it is allowed to say no.
 *
 * The design problem here is that a verdict is a word and the evidence is a
 * handful of numbers, and the word is what people remember. So no number
 * appears alone: each one is drawn on a track showing where the validated
 * glue systems sit and where plausible wrong poses reach, with this value
 * marked on it. A borderline call then looks borderline instead of arriving
 * as a shrug, and a pass that scrapes in looks like a pass that scrapes in.
 *
 * Metrics that failed to separate right poses from wrong ones are shown too,
 * under their own heading, with their measured value and the reason they are
 * not in the verdict. Hiding them would conceal that the most intuitive
 * metric — how big the pocket is — is the one that carries no signal.
 */

(function (GF) {
  'use strict';

  const el = GF.el;
  const fmt = GF.fmt;

  // `.ok`, `.warn` and `.bad` are block-callout classes in this stylesheet:
  // they carry padding, a border and a background. Reusing them as state
  // modifiers turns every status chip into a bordered box. Prefixed.
  const WORDS = {
    pass: { label: 'Pass', cls: 'v-ok' },
    borderline: { label: 'Borderline', cls: 'v-warn' },
    fail: { label: 'Fail', cls: 'v-bad' }
  };

  /* ---- one metric, drawn ------------------------------------------------ */

  /** A track showing the decoy range, the validated range, and this value.
   *
   *  Both ranges are drawn on one axis because the interesting cases are the
   *  ones where they nearly touch — and a number quoted without them is a
   *  number nobody can argue with.
   */
  function track(c) {
    const pts = [c.value, c.native_median, c.native_min,
                 c.decoy_median, c.decoy_p90].filter(Number.isFinite);
    let lo = Math.min.apply(null, pts);
    let hi = Math.max.apply(null, pts);
    const pad = (hi - lo) * 0.12 || Math.abs(hi) * 0.12 || 1;
    lo -= pad; hi += pad;
    const at = (v) => (100 * (v - lo) / (hi - lo));

    // For a higher-is-better metric the decoy band runs from the bottom up to
    // where decoys reach, and the validated band from the weakest validated
    // system upward.
    const up = c.direction > 0;
    const decoyFrom = up ? 0 : at(c.decoy_p90);
    const decoyTo = up ? at(c.decoy_p90) : 100;
    const natFrom = up ? at(c.native_min) : 0;
    const natTo = up ? 100 : at(c.native_min);

    const row = (cls, from, to, label) => el('div', { class: 'lane' }, [
      el('span', { class: 'lane-label', text: label }),
      el('span', { class: 'lane-track' }, [
        el('span', {
          class: 'band ' + cls,
          style: `left:${from}%;width:${Math.max(0.8, to - from)}%`
        })
      ])
    ]);

    return el('div', { class: 'track' }, [
      row('native', natFrom, natTo, 'validated'),
      row('decoy', decoyFrom, decoyTo, 'wrong poses'),
      el('div', { class: 'lane marker-lane' }, [
        el('span', { class: 'lane-label', text: 'this one' }),
        el('span', { class: 'lane-track' }, [
          el('span', {
            class: 'tick',
            style: `left:${at(c.native_median)}%`,
            title: 'median of the validated systems'
          }),
          el('span', {
            class: 'marker ' + WORDS[c.verdict].cls,
            style: `left:${at(c.value)}%`,
            title: 'this interface'
          })
        ])
      ])
    ]);
  }

  function checkRow(c) {
    const w = WORDS[c.verdict];
    const v = Math.abs(c.value) >= 100 ? c.value.toFixed(0) : c.value.toFixed(3);
    return el('div', { class: 'check ' + w.cls }, [
      el('div', { class: 'check-head' }, [
        el('span', { class: 'check-label', text: c.label }),
        el('span', { class: 'check-value mono', text: v + (c.unit ? ' ' + c.unit : '') }),
        el('span', { class: 'check-verdict ' + w.cls, text: w.label })
      ]),
      track(c),
      el('div', { class: 'check-note', text: c.note }),
      el('div', { class: 'check-meta mono' }, [
        el('span', { text:
          'validated ' + c.native_min.toFixed(3) + '–' +
          (c.direction > 0 ? '↑' : '↓') +
          '   decoys reach ' + c.decoy_p90.toFixed(3) +
          '   AUC ' + c.auc.toFixed(2) + '   ' }),
        // Why this cut-off is always one click from the number it governs.
        el('a', { href: '#/calibration/' + c.key, text: 'why this cut-off' })
      ])
    ]);
  }

  /* ---- the screen ------------------------------------------------------- */

  function panel(assessment, onContinue, opts) {
    const v = assessment.verdict || {};
    const w = WORDS[v.verdict] || WORDS.borderline;
    const host = el('div', {});

    host.appendChild(el('div', { class: 'page' }, [
      el('h1', { text: 'Assessment' }),
      el('p', { class: 'lede' }, [
        'Whether this interface has a pocket worth designing chemistry into. ' +
        'No chemistry is done to answer it, and the answer is often no — ' +
        'glue discovery is mostly pair selection.'
      ])
    ]));

    // What was judged. A verdict on a docked pose and a verdict on a solved
    // complex are different claims, and the thresholds behind them were
    // calibrated on crystal structures. Saying which one this is matters
    // more than the verdict itself. On the six bound-bound cases in
    // cache/recovery, the pose this engine puts first was 2.6 A from the
    // crystal placement once and 20-35 A out the other five times, while
    // the cavity at a glue interface tolerates about 2 A. So a verdict on
    // a docked pose is usually describing the docking, not the biology.
    if (v.pose || assessment.pose) {
      host.appendChild(el('div', { class: 'callout note judged-on' }, [
        el('strong', { text: 'Judged on a docked pose' }),
        ' — ',
        el('code', { text: String(v.pose || assessment.pose) }),
        ', not a solved structure. The thresholds below come from crystal '
        + 'complexes. On the six benchmark cases the pose this engine ranks '
        + 'first landed within 2.6 Å of the crystal placement once and 20 Å '
        + 'or more out the other five times, and a cavity here tolerates '
        + 'about 2 Å — so read any result on a docked pose as being first '
        + 'of all about where the search put protein B.'
      ]));
    }

    // The verdict, with its count, so "pass" is visibly 4-of-4 rather than a
    // word that arrived from somewhere.
    host.appendChild(el('section', { class: 'verdict ' + w.cls }, [
      el('div', { class: 'verdict-word', text: w.label }),
      el('div', { class: 'verdict-head', text: v.headline || '' }),
      v.n_checks
        ? el('div', { class: 'verdict-count', text:
            v.n_pass + ' of ' + v.n_checks + ' checks in the validated range' +
            (v.n_fail ? ', ' + v.n_fail + ' inside the decoy range' : '') })
        : null
    ]));

    if (!v.calibrated) {
      host.appendChild(el('div', { class: 'callout warn' }, [
        el('strong', { text: 'Not calibrated. ' }),
        'There is no gate calibration on disk, so the measurements below ' +
        'have nothing to be compared against. Run ',
        el('code', { text: 'sidecar/run_gate.py' }),
        ' to build one.'
      ]));
    }

    const checks = v.checks || [];
    if (checks.length) {
      const sec = el('section', { class: 'panel' }, [
        el('h2', { text: 'What was measured' }),
        el('p', { class: 'sub' }, [
          'Each value is shown against the range the 14 validated glue ' +
          'systems occupy and the range that plausible wrong poses of the ' +
          'same partners reach. Where those two overlap, the measurement ' +
          'cannot separate them and the check says so.'
        ])
      ]);
      checks.forEach((c) => sec.appendChild(checkRow(c)));
      host.appendChild(sec);
    }

    const ex = v.excluded || [];
    if (ex.length) {
      const sec = el('section', { class: 'panel' }, [
        el('h2', { text: 'Measured, and deliberately not in the verdict' }),
        el('p', { class: 'sub' }, [
          'These were computed and then excluded because they do not ' +
          'separate right poses from wrong ones. They are here rather than ' +
          'hidden, because the most intuitive of them is the one that ' +
          'carries the least signal.'
        ])
      ]);
      const t = el('table', { class: 'grid' });
      t.appendChild(el('thead', {}, [el('tr', {}, [
        el('th', { text: 'Metric' }),
        el('th', { text: 'Value' }),
        el('th', { text: 'AUC' }),
        el('th', { text: 'Why not' })
      ])]));
      const b = el('tbody');
      ex.forEach((e) => b.appendChild(el('tr', {}, [
        el('td', { text: e.label }),
        el('td', { class: 'val mono', text: String(e.value) }),
        el('td', { class: 'val mono', text: e.auc.toFixed(2) }),
        el('td', { class: 'dim', text: e.why })
      ])));
      t.appendChild(b);
      sec.appendChild(t);
      host.appendChild(sec);
    }

    (v.caveats || []).forEach((c) => {
      host.appendChild(el('div', { class: 'callout', text: c }));
    });

    // The gate is bypassable by design. What it is not is quiet about it:
    // taking the override is recorded on the run and travels with anything
    // exported from it.
    const next = el('section', { class: 'panel' }, [
      el('h2', { text: 'Chemistry' })
    ]);
    const ready = GF.api.stageReady('chemistry');
    if (v.verdict === 'pass') {
      next.appendChild(el('p', { class: 'sub', text:
        'This interface looks like the validated ones. Growing chemistry ' +
        'into it is a reasonable next step.' }));
    } else {
      next.appendChild(el('div', { class: 'callout warn' }, [
        el('strong', { text: 'The gate did not pass. ' }),
        v.verdict === 'fail'
          ? 'Nothing stops you continuing, but the pocket here looks like ' +
            'the ones wrong poses produce. Compounds grown into it would be ' +
            'optimised against a site that may not exist.'
          : 'The evidence does not separate this interface from a wrong ' +
            'pose. Continuing is reasonable if you have a reason to trust ' +
            'the pose that this screen cannot see.'
      ]));
      next.appendChild(el('p', { class: 'fnote', text:
        'Continuing past a failed gate is recorded on the run and travels ' +
        'with anything exported from it, so a result built on it says so ' +
        'wherever it ends up.' }));
    }
    next.appendChild(el('button', {
      class: 'btn' + (v.verdict === 'pass' ? ' primary' : ' ghost'),
      text: v.verdict === 'pass' ? 'Grow chemistry' : 'Continue anyway',
      disabled: ready ? null : 'disabled',
      onclick: onContinue || null
    }));
    if (!ready) {
      next.appendChild(el('p', { class: 'fnote', text:
        'Chemistry needs rdkit. Install it from the landing page.' }));
    }
    host.appendChild(next);

    if (opts && opts.runId && GF.api.state.live) {
      host.appendChild(el('section', { class: 'panel' }, [
        el('h2', { text: 'Share this' }),
        el('p', { class: 'sub' }, [
          'Saves as one HTML file that opens on a machine with no Python and ' +
          'nothing installed. The calibration travels inside it, so whoever ' +
          'you send it to can click through to the evidence behind every ' +
          'threshold rather than taking the verdict on trust.'
        ]),
        el('a', {
          class: 'btn ghost',
          href: '/api/run/' + encodeURIComponent(opts.runId) + '/export',
          download: '',
          text: 'Export as a single file'
        }),
        el('p', { class: 'fnote', text:
          'About 4 MB — most of it is the 14-system validation set that makes ' +
          '"why this cut-off" work offline.' })
      ]));
    }

    return GF.wrapTables(host);
  }

  /* ---- an exported copy -------------------------------------------------- */

  /** What an exported file opens on.
   *
   *  It has to stand on its own: whoever receives it was not here when the
   *  run happened and has no way to ask. So it leads with what was measured
   *  and against what, rather than with a verdict and a date.
   */
  function exported() {
    const e = GF.EXPORT || {};
    const run = e.run || {};
    const host = el('div', {});

    if (!run.assessment) {
      host.appendChild(el('div', { class: 'page' }, [
        el('h1', { text: 'Nothing in this file' }),
        el('p', { class: 'lede', text:
          'This looks like an exported GlueForge page but it carries no ' +
          'assessment.' })
      ]));
      return host;
    }

    const slots = run.slots || {};
    const a = slots.A || {};
    const b = slots.B || {};
    const face = run.face || {};

    host.appendChild(el('div', { class: 'page' }, [
      el('h1', { text: 'Exported assessment' }),
      el('p', { class: 'lede' }, [
        'A finished GlueForge run, saved as one file. Nothing here needs a ' +
        'server; the thresholds and the evidence behind them travel inside ' +
        'it.'
      ])
    ]));

    const facts = el('div', { class: 'facts' });
    const fact = (label, value, note, cls) => facts.appendChild(
      el('div', { class: 'fact' }, [
        el('div', { class: 'flabel', text: label }),
        el('div', { class: 'fval' + (cls ? ' ' + cls : ''), text: value }),
        note ? el('div', { class: 'fnote', text: note }) : null
      ]));
    fact('Protein A', a.entry || a.name || '—', 'holds the ligand', 'name');
    fact('Protein B', b.entry || b.name || '—', 'recruited', 'name');
    fact('Ligand', (run.ligand && run.ligand.comp_id) || '—',
         face.ligand_heavy_atoms ? face.ligand_heavy_atoms + ' heavy atoms' : null);
    fact('Exposed face', face.presented_area !== undefined
         ? face.presented_area.toFixed(1) + ' Å²' : '—',
         face.buried_fraction !== undefined
         ? (100 * face.buried_fraction).toFixed(0) + '% buried' : null);
    host.appendChild(el('section', { class: 'panel' }, [facts]));

    host.appendChild(panel(run.assessment, null, null));

    host.appendChild(el('p', { class: 'fnote' }, [
      'run ', el('code', { text: run.id || '?' }),
      run.exported ? ' · exported ' + new Date(e.exported * 1000).toISOString().slice(0, 10) : '',
      ' · GlueForge ', el('code', { text: e.version || '?' })
    ]));
    return GF.wrapTables(host);
  }

  GF.pages = GF.pages || {};
  GF.pages.exported = exported;

  GF.assess = { panel: panel, track: track, WORDS: WORDS, exported: exported };
})(globalThis.GF || (globalThis.GF = {}));
