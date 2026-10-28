/* GlueForge — the docking screen.
 *
 * Choose how long to spend, watch it happen, then look at what came back.
 *
 * Two things here are deliberate and both are about not overselling a result.
 *
 * Each preset shows its measured pose-recovery rate, or says it has none.
 * The run times are estimates from LightDock's own published figure; the
 * recovery rates have to be measured on a machine that has LightDock
 * installed, and until that has happened the preset says so. A preset
 * showing an invented recovery rate next to a real run time is worse than
 * one admitting ignorance, because the invented number is the one that gets
 * quoted.
 *
 * And poses are ranked by how much of the presented face they bury, not by
 * docking score. A glue interface is often not especially favourable; what
 * makes it a glue interface is that the ligand is in it. Both numbers are
 * shown, because where they disagree is the interesting part.
 */

(function (GF) {
  'use strict';

  const el = GF.el;
  const api = GF.api;
  const fmt = GF.fmt;

  const S = {
    runId: null,
    run: null,
    presets: null,
    engine: 'builtin',
    chosen: 'standard',
    job: null,
    lines: [],
    state: null,
    error: null,
    result: null,
    selected: 0,
    // The 3D panels, and which pose they are currently showing. Kept here
    // rather than rebuilt on every render: re-creating the viewer would
    // throw away the rotation the user had set, which on a structure is
    // most of what they were doing.
    panels: null,
    shownPose: null,
    viewError: null,
    assessing: false,
    assessError: null
  };

  /* ---- choosing a preset ------------------------------------------------ */

  /** Which engine's presets are on screen. The two search differently and
   *  are measured separately, so they never share a card or a number. */
  function engineInfo() {
    const e = S.presets && S.presets.engines;
    if (!e) return null;
    return e[S.engine] || null;
  }

  function presetCard(p, measured) {
    const on = S.chosen === p.key;
    const lo = p.minutes_estimate ? p.minutes_estimate[0] : p.minutes_low;
    const hi = p.minutes_estimate ? p.minutes_estimate[1] : p.minutes_high;
    const card = el('button', {
      class: 'preset' + (on ? ' on' : ''),
      onclick: () => { S.chosen = p.key; render(); }
    }, [
      el('div', { class: 'preset-name', text: p.label }),
      el('div', { class: 'preset-time mono', text: lo + '–' + hi + ' min' }),
      p.note ? el('div', { class: 'preset-note', text: p.note }) : null
    ]);

    // The built-in engine's measurement covers the engine, not each preset,
    // so it is shown once above the grid rather than copied onto three cards
    // as if three things had been measured.
    const m = p.measured !== undefined ? p.measured : null;
    if (measured !== false) {
      card.appendChild(el('div', { class: 'preset-rec' + (m ? '' : ' unknown') },
        m
          ? [
              el('strong', { text: m.recovered + '/' + m.n }),
              el('span', { text: ' recovered the crystal pose within ' +
                (m.cutoff || 5) + ' Å' }),
              m.designable !== undefined
                ? el('span', { class: 'dim', text:
                    '  ·  ' + m.designable + '/' + m.n + ' close enough to design into' })
                : null
            ]
          : [
              el('span', { text: 'recovery rate not measured yet' })
            ]
      ));
    }

    card.appendChild(el('div', { class: 'preset-spec mono', text:
      p.engine === 'builtin'
        ? p.rotations + ' orientations · ' + p.spacing + ' Å grid · ' +
          (p.focus ? p.focus[0] + ' revisited' : 'one sweep')
        : p.swarms + ' swarms · ' + p.glowworms + ' glowworms · ' +
          p.steps + ' steps' }));
    return card;
  }

  /** Pick an engine, when there is a choice to make.
   *
   *  There usually is not: without LightDock the built-in engine is all
   *  there is, which is every Windows machine. The row is drawn when both
   *  exist, and describes each by how it searches rather than ranking them.
   *  They have never been run against each other on this benchmark -- only
   *  the built-in one is measured -- so "stronger" and "weaker" were words
   *  standing in for a comparison nobody had made, which is the one kind of
   *  claim this tool is supposed to not make.
   */
  function engineChooser() {
    const engines = S.presets && S.presets.engines;
    if (!engines || !engines.lightdock || !engines.lightdock.available) {
      return null;
    }
    const row = el('div', { class: 'engines' });
    [['lightdock', 'LightDock',
      'Swarm docking. Not measured here, so its rate on this benchmark ' +
      'is unknown.'],
     ['builtin', 'Built-in',
      'Correlation docking. Measured: see the rate below. Slows down on ' +
      'large partners.']
    ].forEach(([key, label, note]) => {
      row.appendChild(el('button', {
        class: 'engine' + (S.engine === key ? ' on' : ''),
        onclick: () => {
          S.engine = key;
          const list = engines[key].presets;
          if (!list.some((p) => p.key === S.chosen)) S.chosen = list[0].key;
          render();
        }
      }, [
        el('div', { class: 'engine-name', text: label }),
        el('div', { class: 'engine-note', text: note })
      ]));
    });
    return row;
  }

  function chooser() {
    const sec = el('section', { class: 'panel' }, [
      el('h2', { text: 'How long to spend' })
    ]);

    if (!S.presets) {
      sec.appendChild(el('p', { class: 'sub', text: 'Loading…' }));
      return sec;
    }

    // Guarded, because `engineChooser` returns null when there is only one
    // engine to choose from -- which is every Windows machine, LightDock
    // having no build for it. `el()` skips nulls inside a children array, so
    // the pattern is safe everywhere it is used that way and not here, and
    // an unguarded appendChild(null) throws. This screen had no render test,
    // so nothing caught it until it was in someone's hands.
    const engines = engineChooser();
    if (engines) sec.appendChild(engines);

    const info = engineInfo();
    const builtin = S.engine === 'builtin';

    sec.appendChild(el('p', { class: 'sub' }, [
      builtin
        ? 'Times are for the largest partner in the benchmark, 6500 atoms. ' +
          'The grid is sized by the partner’s reach and the cost goes as its ' +
          'cube, so a small protein is much quicker than the figure shown ' +
          'and a very large one is not.'
        : 'Times are estimates from LightDock’s own figure of one to two ' +
          'minutes per ten steps per swarm, on eight cores. Restraints cut ' +
          'the swarm count a long way, so a real run is usually faster than ' +
          'the upper end.'
    ]));

    if (builtin) {
      const m = info && info.measured;
      sec.appendChild(el('div', { class: 'callout ' + (m ? 'note' : 'warn') },
        m
          ? [
              el('strong', { text: m.recovered + ' of ' + m.n +
                ' benchmark systems recovered. ' }),
              'Within ' + (m.cutoff || 5) + ' Å of the crystal pose, ' +
              m.designable + ' of them close enough to design into. ' +
              'Measured on the ' + (m.preset || 'quick') + ' setting' +
              (m.bound_bound
                ? ', and on pairs taken from solved complexes, so both ' +
                  'halves are already in their bound shape — real use is ' +
                  'harder than this number.'
                : '.')
            ]
          : [
              el('strong', { text: 'Recovery rate not measured here. ' }),
              'Nothing on this machine has checked how often this engine ' +
              'finds the crystal pose. Running ',
              el('code', { text: 'python sidecar/run_builtin_recovery.py' }),
              ' fills it in from the benchmark systems.'
            ]
      ));
    } else if (!S.presets.calibrated) {
      sec.appendChild(el('div', { class: 'callout warn' }, [
        el('strong', { text: 'No recovery rates yet. ' }),
        'Nothing has measured how often each of these finds the crystal ' +
        'pose. Running ',
        el('code', { text: 'python sidecar/run_recovery.py --all-presets' }),
        ' fills these in from the 14 benchmark systems. Until then the ' +
        'difference between the presets is time spent, not accuracy known.'
      ]));
    }

    const grid = el('div', { class: 'presets' });
    const list = info ? info.presets : S.presets.presets;
    list.forEach((p) => grid.appendChild(presetCard(p, !builtin)));
    sec.appendChild(grid);

    const face = S.run && S.run.face;
    if (face && !face.directional) {
      sec.appendChild(el('div', { class: 'callout warn' }, [
        el('strong', { text: 'Weakly restrained. ' }),
        'This ligand sits in a groove rather than a pocket, so no residue ' +
        'could honestly be made a requirement and the whole rim went in as ' +
        'permitted. The search is narrowed less than usual; a longer preset ' +
        'buys more here than it would otherwise.'
      ]));
    }

    const r = (S.run && S.run.restraintSummary) ||
              (face && face.restraints) || null;
    if (r) {
      sec.appendChild(el('p', { class: 'fnote', text:
        'Restraint: ' + r.n_active + ' active, ' + r.n_passive + ' passive, ' +
        'from the ' + (face ? face.n_rim : '?') + ' residues lining the face.' }));
    }

    sec.appendChild(el('button', {
      class: 'btn primary',
      text: 'Start docking',
      disabled: S.state === 'running' ? 'disabled' : null,
      onclick: start
    }));
    if (S.error) {
      sec.appendChild(el('p', { class: 'drop-err', text: S.error }));
    }
    return sec;
  }

  /* ---- while it runs ---------------------------------------------------- */

  function progress() {
    if (!S.job && S.state !== 'running') return null;
    const sec = el('section', { class: 'panel' }, [
      el('h2', { text: S.state === 'running' ? 'Docking' : 'Docking output' }),
      el('p', { class: 'sub', text:
        S.state === 'running'
          ? 'This takes a while. You can close the page — the run keeps going ' +
            'and the results land in the run folder either way.'
          : 'What LightDock said.' })
    ]);
    const log = el('pre', { class: 'dep-log open' });
    log.textContent = S.lines.join('\n');
    sec.appendChild(log);
    return sec;
  }

  /* ---- what came back --------------------------------------------------- */

  /** The poses on offer, in the order the table shows them. The viewer and
   *  the table have to agree about what row 3 means, so they ask here. */
  function selectable() {
    const res = S.result;
    if (!res || !res.poses || !res.poses.length) return [];
    const rows = res.poses.slice();
    if (res.engine !== 'builtin') {
      rows.sort((a, b) => b.face_burial - a.face_burial);
    }
    const viable = rows.filter((p) => p.viable);
    return (viable.length ? viable : rows).slice(0, 25);
  }

  function poses() {
    const res = S.result;
    if (!res || !res.poses || !res.poses.length) return null;

    const prefiltered = res.engine === 'builtin';
    const rows = res.poses;
    const viable = rows.filter((p) => p.viable);
    const inPocket = rows.filter((p) => p.ligand_clashes > 2).length;
    // One list, shared with the viewer, so the highlighted row and the
    // structure on screen are the same pose. They were computed separately
    // at first and drifted the moment the ordering rule changed.
    const shown = selectable();

    const sec = el('section', { class: 'panel' }, [
      el('h2', { text: 'Poses' }),
      el('p', { class: 'sub' }, [
        prefiltered
          ? 'In the engine’s own order. Every pose here already had to ' +
            'reach the presented face to be returned at all, so ranking ' +
            'them by how much of it they bury sorts by how many atoms each ' +
            'piles on rather than by which is right. Face burial is still ' +
            'shown, and where it disagrees with the order is worth a look.'
          : 'Ranked by how much of the presented face each one buries, not ' +
            'by docking score. A glue interface is often not especially ' +
            'favourable; what makes it one is that the ligand is in it. ' +
            'Both numbers are here, and where they disagree is worth a look.'
      ])
    ]);

    sec.appendChild(el('p', { class: 'fnote', text:
      rows.length + ' poses, ' + viable.length + ' viable' +
      (inPocket
        ? ' · ' + inPocket + ' had settled into the ligand’s own pocket ' +
          'and were rejected'
        : '') }));

    const t = el('table', { class: 'grid' });
    t.appendChild(el('thead', {}, [el('tr', {}, [
      el('th', { text: '' }),
      el('th', { text: 'Pose' }),
      el('th', { text: 'Face buried' }),
      el('th', { text: 'Restraints met' }),
      el('th', { text: 'Contacts' }),
      el('th', { text: 'Clashes' }),
      el('th', { text: 'On ligand' })
    ])]));
    const b = el('tbody');
    shown.forEach((p, i) => {
      const tr = el('tr', { class: i === S.selected ? 'on' : '' }, [
        el('td', {}, [el('input', {
          type: 'radio', name: 'pose', checked: i === S.selected ? '' : null,
          onchange: () => { S.selected = i; render(); }
        })]),
        el('td', { class: 'mono', text: (p.label || '').slice(0, 26) }),
        el('td', { class: 'val', text: (100 * p.face_burial).toFixed(1) + '%' }),
        el('td', { class: 'val', text:
          (100 * (p.restraint_satisfaction || 0)).toFixed(0) + '%' }),
        el('td', { class: 'val', text: String(p.contacts) }),
        el('td', { class: 'val', text: String(p.clashes) }),
        el('td', { class: 'val', text: String(p.ligand_clashes) })
      ]);
      b.appendChild(tr);
    });
    t.appendChild(b);
    sec.appendChild(t);

    sec.appendChild(el('div', { class: 'callout' }, [
      el('strong', { text: 'Pose quality limits everything after this. ' }),
      'The cavity at the interface tolerates about 2 Å of error in the ' +
      'partner’s placement and falls apart by 3. A pose you cannot trust ' +
      'to that will not give a cavity worth designing into, however good its ' +
      'score looks.'
    ]));

    // It runs the gate on the selected pose and then goes to the result.
    // It used to only navigate, which meant the button named the thing it
    // did not do and the next screen said nothing had been assessed --
    // because nothing had.
    const go = el('button', {
      class: 'btn primary',
      text: S.assessing ? 'Assessing…' : 'Assess this pose'
    });
    go.disabled = S.assessing ? true : null;
    go.addEventListener('click', () => {
      const list = selectable();
      const chosen = list[Math.min(S.selected, list.length - 1)];
      if (!chosen) return;
      S.assessing = true;
      S.assessError = null;
      render();
      api.json('/api/run/' + S.runId + '/assess', { pose: chosen.label })
        .then(() => { location.hash = '#/assess/' + S.runId; })
        .catch((err) => {
          S.assessing = false;
          S.assessError = err && err.message ? err.message : String(err);
          render();
        });
    });
    sec.appendChild(go);
    if (S.assessError) {
      sec.appendChild(el('div', { class: 'callout warn' }, [
        el('strong', { text: 'Could not assess that pose. ' }), S.assessError
      ]));
    }
    return sec;
  }

  /* ---- server ----------------------------------------------------------- */

  async function start() {
    S.error = null;
    S.lines = [];
    S.state = 'running';
    render();
    try {
      const out = await api.json('/api/run/' + S.runId + '/dock',
                                 { preset: S.chosen, engine: S.engine });
      S.job = out.job;
      api.follow(out.job, (line) => {
        S.lines.push(line);
        if (S.lines.length > 600) S.lines.splice(0, 200);
        const log = document.querySelector('.dep-log.open');
        if (log) { log.textContent = S.lines.join('\n'); log.scrollTop = log.scrollHeight; }
      }, async (done) => {
        S.state = done && done.code === 0 ? 'done' : 'failed';
        if (S.state === 'failed') {
          S.error = 'Docking exited with code ' +
            (done ? done.code : '?') + '. The output above says why.';
        }
        await refresh();
        render();
      });
    } catch (err) {
      S.state = 'failed';
      S.error = err.message;
      // The server usually says what to do about it already; only add the
      // hint when it has not, or the user reads the same sentence twice.
      if (err.missing && !/landing page/i.test(S.error)) {
        S.error += ' Install it from the landing page, then reopen this run.';
      }
      render();
    }
  }

  async function refresh() {
    try {
      const d = await api.call('/api/run/' + S.runId + '/dock');
      S.state = d.state === 'running' ? 'running' : (S.state || d.state);
      if (d.result) S.result = d.result;
      S.run = await api.runs.get(S.runId);
    } catch (err) { /* leave what we have */ }
  }

  /* ---- assembly ---------------------------------------------------------- */

  let host = null;

  /** Draw the selected pose, fetching it if it is not the one on screen. */
  function showPose(label) {
    if (!label || S.shownPose === label) return;
    S.shownPose = label;
    S.viewError = null;
    GF.view.forRun(S.runId, label).then((blob) => {
      if (S.shownPose !== label) return;   // selection moved on while waiting
      if (!S.panels) {
        S.panels = {
          structure: GF.view.panel(blob, { title: 'The placement' }),
          cavity: GF.view.cavityPanel(blob)
        };
      } else {
        S.panels.structure.setSystem(blob);
        S.panels.cavity.setSystem(blob);
      }
      render();
    }).catch((err) => {
      if (S.shownPose !== label) return;
      S.viewError = err && err.message ? err.message : String(err);
      render();
    });
  }

  function render() {
    if (!host) return;
    host.innerHTML = '';
    host.appendChild(el('div', { class: 'page' }, [
      el('h1', { text: 'Dock' }),
      el('p', { class: 'lede' }, [
        'Place protein B against protein A, restrained to the face the ' +
        'ligand presents. That restraint is the whole of what this tool ' +
        'contributes to the search — it comes from the ligand, which is ' +
        'information ordinary docking does not have.'
      ]),
      S.runId ? el('p', { class: 'fnote mono', text: 'run ' + S.runId }) : null
    ]));

    host.appendChild(chooser());
    const p = progress();
    if (p) host.appendChild(p);
    const r = poses();
    if (r) host.appendChild(r);

    // The structure, under the table it is a row of. Asked for whenever the
    // selection and what is on screen have diverged -- which includes the
    // first time results arrive, so the picture appears on its own rather
    // than waiting for someone to click the row that is already selected.
    const list = selectable();
    if (list.length) {
      const want = (list[Math.min(S.selected, list.length - 1)] || {}).label;
      showPose(want);
      if (S.viewError) {
        host.appendChild(el('div', { class: 'callout warn' }, [
          el('strong', { text: 'Could not draw this pose. ' }),
          S.viewError
        ]));
      } else if (S.panels) {
        host.appendChild(S.panels.structure.node);
        host.appendChild(S.panels.cavity.node);
      } else {
        host.appendChild(el('p', { class: 'sub', text: 'Drawing…' }));
      }
    }
    GF.wrapTables(host);
  }

  function page(runId) {
    host = el('div', {});
    if (S.runId !== runId) {
      S.runId = runId;
      S.job = null; S.lines = []; S.state = null;
      S.error = null; S.result = null; S.selected = 0;
    }
    render();

    Promise.all([
      S.presets ? Promise.resolve(S.presets) : api.call('/api/presets'),
      api.runs.get(runId).catch(() => null)
    ]).then(([presets, run]) => {
      S.presets = presets;
      if (presets && presets.engine_default) S.engine = presets.engine_default;
      if (presets && presets.default && !S.chosen) S.chosen = presets.default;
      S.run = run;
      return refresh();
    }).then(render).catch(() => render());

    return host;
  }

  GF.pages = GF.pages || {};
  GF.pages.dock = page;

  /** A seam for the render test, and the reason it exists is worth stating.
   *
   *  This screen's shape is decided entirely by what `/api/presets` returned
   *  -- which engines exist, which of them are installed -- and none of that
   *  is reachable by calling `pages.dock` alone, because the first render
   *  happens before the request resolves. So a test that only builds the
   *  page exercises the loading state and nothing else, and the screen can
   *  be broken in every state that matters while the suite stays green.
   *
   *  That is not hypothetical. The picker is null when there is only one
   *  engine to choose between, which is every Windows machine, and appending
   *  it unguarded threw. Nobody found it until it was in front of someone.
   */
  GF.dock = { state: S, render: render };
})(globalThis.GF || (globalThis.GF = {}));
