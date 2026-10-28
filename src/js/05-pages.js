/* GlueForge — the pages.
 *
 * Three of them. The benchmark page is the deliverable; the system page
 * exists so you can see how a number was arrived at; the method page records
 * what was computed with what.
 *
 * Every metric carries a glossary affordance (see 01-help.js) because the
 * numbers here are only meaningful if you know which of them is an effect
 * size, which is a raw difference, and which is controlled for what.
 */

(function (GF) {
  'use strict';

  const el = GF.el, fmt = GF.fmt, help = GF.help;

  /* ===================================================================== */
  /* Benchmark                                                             */
  /* ===================================================================== */

  function benchmarkPage() {
    const ix = GF.data.index();
    const root = el('div', { class: 'page' });

    const diffs = ix.systems.map((s) => s.result.difference);
    const strat = ix.systems
      .map((s) => s.stratified)
      .filter((v) => v !== null && v !== undefined);
    const sign = GF.stat.signTest(diffs);
    const signStrat = GF.stat.signTest(strat);

    root.appendChild(el('h1', { text: 'Does frustration mark where the glue sits?' }));
    root.appendChild(el('p', {
      class: 'lede',
      text:
        'Each row is one crystal ternary complex with the glue deleted. The ' +
        'number is the difference in mean frustration between interface ' +
        'contacts next to where the glue sat and the rest of that interface. ' +
        'Negative is the direction the hypothesis predicts.'
    }));

    if (!ix.potential_published) {
      root.appendChild(el('div', {
        class: 'warn',
        text:
          'Running on the stand-in potential (' + ix.potential + '), not the ' +
          'published AWSEM tables. Every number below is a statement about ' +
          'this potential, not about AWSEM. Drop the --gamma flag to use the ' +
          'bundled AWSEM_2015.json instead.'
      }));
    } else {
      root.appendChild(el('div', {
        class: 'note',
        text:
          'Published AWSEM-MD (2015) tables, contact and burial terms. An ' +
          'earlier run on a stand-in potential gave the same direction with ' +
          'roughly twice the effect size, so the sign is robust to the ' +
          'potential and the magnitude is not \u2014 see Method.'
      }));
    }

    const cards = el('div', { class: 'cards' });
    cards.appendChild(statCard(
      sign.negative + ' / ' + sign.n,
      'in the predicted direction',
      'sign test p = ' + fmt.p(sign.pGreater),
      'difference'
    ));
    cards.appendChild(statCard(
      fmt.signed(GF.stat.mean(diffs), 3),
      'mean difference',
      'raw, across all structures',
      'difference'
    ));
    cards.appendChild(statCard(
      fmt.signed(GF.stat.mean(strat), 3),
      'burial-controlled',
      signStrat.negative + ' / ' + signStrat.n + ' negative',
      'burial_controlled'
    ));
    root.appendChild(cards);

    const table = el('table', { class: 'grid' });
    table.appendChild(el('thead', {}, [
      el('tr', {}, [
        el('th', { text: 'Entry' }),
        el('th', { text: 'Ligand' }),
        headCell('Interface', 'interface_contacts'),
        headCell('Prox / distal', 'proximal_distal'),
        headCell('Δ raw', 'difference'),
        headCell('Δ burial-ctrl', 'burial_controlled'),
        headCell('rank-biserial', 'rank_biserial'),
        headCell('gradient ρ', 'gradient'),
        headCell('p', 'p_permutation')
      ])
    ]));

    const body = el('tbody');
    GF.data.groups().forEach((g) => {
      const gdiffs = g.rows.map((r) => r.result.difference);
      const allNeg = gdiffs.every((d) => d < 0);
      body.appendChild(el('tr', { class: 'group' }, [
        el('td', { colspan: '10' }, [
          el('span', { class: 'gname', text: g.label }),
          el('span', {
            class: 'gmeta',
            text:
              g.rows.length + ' structure' + (g.rows.length === 1 ? '' : 's') +
              ' · mean ' + fmt.signed(GF.stat.mean(gdiffs), 3) +
              (allNeg ? ' · all in predicted direction' : '')
          })
        ])
      ]));

      g.rows.forEach((s) => {
        const r = s.result;
        const tr = el('tr', { class: 'row' });
        tr.appendChild(el('td', {}, [
          el('a', { href: '#/system/' + s.entry, class: 'entry', text: s.entry })
        ]));
        tr.appendChild(el('td', { class: 'mono', text: s.ligand }));
        tr.appendChild(el('td', { class: 'num', text: fmt.int(s.n_interface) }));
        tr.appendChild(el('td', {
          class: 'num dim', text: s.n_proximal + ' / ' + s.n_distal
        }));
        tr.appendChild(numCell(r.difference));
        tr.appendChild(numCell(s.stratified));
        tr.appendChild(el('td', {
          class: 'num dim', text: fmt.signed(r.rank_biserial, 3)
        }));
        tr.appendChild(el('td', {
          class: 'num ' + (s.gradient_p !== null && s.gradient_p < 0.05 ? 'sig' : 'dim'),
          text: s.gradient_rho === null || s.gradient_rho === undefined
            ? '—' : fmt.signed(s.gradient_rho, 3)
        }));
        tr.appendChild(el('td', {
          class: 'num ' + (r.p_value < 0.05 ? 'sig' : 'dim'),
          text: fmt.p(r.p_value)
        }));
        body.appendChild(tr);
      });
    });
    table.appendChild(body);
    root.appendChild(table);

    root.appendChild(el('div', { class: 'reading' }, [
      el('p', {
        text:
          'Read the groups, not the rows. The eight CDK12 structures are one ' +
          'interface with eight different ligands, so they are eight ' +
          'measurements of a stable effect rather than eight independent ' +
          'tests of the hypothesis. Across biological systems the count is ' +
          'much smaller than fourteen.'
      }),
      el('p', {
        text:
          'The two DCAF15 forms are the most informative comparison here: ' +
          'the same pair, different crystal forms, and they disagree in sign. ' +
          'That bounds reproducibility more honestly than any p-value in the ' +
          'table.'
      }),
      el('p', {
        text:
          'The gradient column is the threshold-free version of the same ' +
          'question: the rank correlation between a contact\'s frustration ' +
          'and its distance to the ligand, with no cut-off involved. It ' +
          'agrees with the split comparison, which means the result is not ' +
          'an artefact of where the 6 Å line was drawn.'
      })
    ]));

    return root;
  }

  function headCell(text, slug) {
    const th = el('th', { class: 'num' }, [el('span', { text: text })]);
    return help.attach(th, slug);
  }

  function statCard(value, label, note, slug) {
    const card = el('div', { class: 'card' }, [
      el('div', { class: 'cval', text: value })
    ]);
    const lab = el('div', { class: 'clabel' }, [el('span', { text: label })]);
    help.attach(lab, slug);
    card.appendChild(lab);
    card.appendChild(el('div', { class: 'cnote', text: note }));
    return card;
  }

  function numCell(v) {
    if (v === null || v === undefined) {
      return el('td', { class: 'num dim', text: '—' });
    }
    const td = el('td', { class: 'num val', text: fmt.signed(v, 3) });
    td.style.color = GF.frustrationColour(v < 0 ? -1.6 : 1.6, 1);
    return td;
  }

  /* ===================================================================== */
  /* One system                                                            */
  /* ===================================================================== */

  function systemPage(entry) {
    const sys = GF.data.system(entry);
    const root = el('div', { class: 'page' });
    if (!sys) {
      root.appendChild(el('p', { text: 'No cached data for ' + entry }));
      return root;
    }

    const a = sys.assignment;
    root.appendChild(el('h1', { text: entry + ' · ' + a.ligand }));

    const lede = el('p', { class: 'lede' }, [
      el('span', { text: 'Anchor ' }),
      el('strong', { text: a.anchor_name }),
      el('span', { text: ' (chain ' + a.anchor_chain + '), partner ' }),
      el('strong', { text: a.partner_name }),
      el('span', { text: ' (chain ' + a.partner_chain + '). The ligand puts ' }),
      el('strong', { text: Math.round(a.anchor_share * 100) + '%' }),
      el('span', {
        text:
          ' of its protein contacts into the anchor, which is how the roles ' +
          'were assigned — nothing here was configured by hand.'
      })
    ]);
    help.attach(lede, 'anchor_share', { inline: true });
    root.appendChild(lede);

    if (sys.cavity && sys.cavity.points) {
      root.appendChild(el('p', { class: 'crosslink' }, [
        el('a', {
          href: '#/cavity/' + entry,
          class: 'entry',
          text: 'Open the cavity at this site \u2192'
        }),
        el('span', {
          class: 'dim',
          text:
            '  ' + Math.round(sys.cavity.summary.volume) + ' \u00c5\u00b3, ' +
            Math.round(sys.cavity.anchor_lined_share * 100) + '% anchor-lined'
        })
      ]));
    }

    const split = el('div', { class: 'split' });

    /* --- left: viewer -------------------------------------------------- */
    const left = el('div', { class: 'panel' });
    const vhead = el('div', { class: 'phead' }, [
      el('h3', { text: 'Structure' })
    ]);
    help.attach(vhead, 'viewer');
    left.appendChild(vhead);

    const canvas = el('canvas', { class: 'viewer' });
    left.appendChild(canvas);
    const readout = el('div', { class: 'readout', text: 'Hover a contact to inspect it.' });
    left.appendChild(readout);
    const legend = el('div', { class: 'vlegend' });
    left.appendChild(legend);
    const controls = el('div', { class: 'controls' });
    left.appendChild(controls);
    split.appendChild(left);

    /* --- right: numbers ------------------------------------------------ */
    const right = el('div', { class: 'panel' });
    right.appendChild(el('h3', { text: 'Result' }));
    right.appendChild(resultTable(sys));

    const hhead = el('div', { class: 'phead' }, [
      el('h3', { text: 'Frustration across the interface' })
    ]);
    help.attach(hhead, 'histogram');
    right.appendChild(hhead);

    const histBox = el('div', { class: 'hist' });
    right.appendChild(histBox);
    const liveBox = el('div', { class: 'live' });
    right.appendChild(liveBox);
    split.appendChild(right);

    root.appendChild(split);

    const vectorHost = el('div', { class: 'vectors-section' });
    root.appendChild(vectorHost);

    requestAnimationFrame(function () {
      const v = new GF.Viewer(canvas, sys, function (hit) {
        renderReadout(readout, hit);
      });
      v.draw();
      renderLegend(legend, v);
      const redrawPlot = function (live) {
        drawPlot(histBox, sys, live, function (k) {
          v.highlight = k;
          v.draw();
        });
      };
      buildControls(controls, v, sys, redrawPlot, liveBox, legend);
      redrawPlot(null);
      buildVectorPanel(vectorHost, sys, v, legend);
      window.addEventListener('resize', function () { v.draw(); });
    });

    return root;
  }

  function renderLegend(host, viewer) {
    host.innerHTML = '';
    viewer.legend().forEach(function (row) {
      const sw = el('span', { class: 'lsw' });
      sw.style.background = row.swatch;
      host.appendChild(el('div', { class: 'lrow' }, [
        sw,
        el('span', { class: 'llabel', text: row.label }),
        el('span', { class: 'lnote', text: row.note })
      ]));
    });
  }

  function renderReadout(host, hit) {
    host.innerHTML = '';
    if (!hit) {
      host.textContent = 'Hover a contact to inspect it.';
      host.classList.remove('on');
      return;
    }
    host.classList.add('on');
    const chip = el('span', { class: 'verdict', text: hit.verdict });
    chip.style.color = GF.frustrationColour(hit.frustration, 1);
    host.appendChild(el('span', { class: 'mono rpair', text: hit.a + ' — ' + hit.b }));
    host.appendChild(chip);
    host.appendChild(el('span', {
      class: 'mono dim',
      text:
        'F ' + fmt.signed(hit.frustration, 2) +
        ' · ' + fmt.fixed(hit.distance, 1) + ' Å apart' +
        ' · ' + fmt.fixed(hit.ligandDistance, 1) + ' Å from ligand' +
        (hit.iface ? (hit.proximal ? ' · glue-proximal' : ' · interface, distal')
                   : ' · within one chain')
    }));
  }

  function resultTable(sys) {
    const r = sys.result, c = r.comparison;
    const t = el('table', { class: 'kv' });
    const rows = [
      ['Interface contacts', fmt.int(r.contacts_interface), 'interface_contacts'],
      ['Glue-proximal / distal', r.contacts_proximal + ' / ' + r.contacts_distal, 'proximal_distal'],
      ['Mean F, proximal', fmt.signed(c.mean_proximal, 3), 'frustration'],
      ['Mean F, distal', fmt.signed(c.mean_distal, 3), 'frustration'],
      ['Difference', fmt.signed(c.difference, 3), 'difference'],
      ['Burial-controlled difference', fmt.signed(r.stratified.pooled_difference, 3), 'burial_controlled'],
      ['Rank-biserial', fmt.signed(c.rank_biserial, 3), 'rank_biserial'],
      ['p (permutation)', fmt.p(c.p_value), 'p_permutation'],
      ['Mean density, proximal', fmt.fixed(r.mean_density_proximal, 2), 'density'],
      ['Mean density, distal', fmt.fixed(r.mean_density_distal, 2), 'density'],
      ['Gradient, rho(F, distance)',
       r.gradient ? fmt.signed(r.gradient.spearman_rho, 3) : '—', 'gradient'],
      ['Gradient p',
       r.gradient ? fmt.p(r.gradient.p_value) : '—', 'gradient'],
      ['Sign holds 4 to 10 Å',
       r.gradient ? (r.gradient.sweep_all_same_sign ? 'yes' : 'no') : '—',
       'gradient']
    ];
    rows.forEach(function (kv) {
      const k = el('td', { class: 'k' }, [el('span', { text: kv[0] })]);
      if (kv[2]) help.attach(k, kv[2]);
      t.appendChild(el('tr', {}, [
        k, el('td', { class: 'v mono', text: kv[1] })
      ]));
    });
    return t;
  }

  /* --- controls --------------------------------------------------------- */

  function buildControls(host, viewer, sys, redrawPlot, liveBox, legendHost) {
    const styleLab = el('div', { class: 'clab' }, [el('span', { text: 'Style' })]);
    help.attach(styleLab, 'style');
    const styleSeg = el('div', { class: 'seg' });
    // Only the styles this system can be drawn in. A docking result has no
    // frustration analysis behind it, so there are no per-contact energies
    // to colour by and no exit vectors to trace; offering those buttons
    // would promise a view that cannot be produced.
    viewer.styles().forEach(function (st) {
      const b = el('button', {
        class: 'seg-btn' + (st.id === viewer.style ? ' on' : ''),
        text: st.label
      });
      b.addEventListener('click', function () {
        viewer.style = st.id;
        viewer.activeVector = -1;
        viewer.draw();
        renderLegend(legendHost, viewer);
        Array.from(styleSeg.children).forEach(function (c) { c.classList.remove('on'); });
        b.classList.add('on');
      });
      styleSeg.appendChild(b);
    });
    host.appendChild(el('div', { class: 'cgroup' }, [styleLab, styleSeg]));

    host.appendChild(el('div', { class: 'cgroup' }, [
      el('div', { class: 'clab', text: 'Show' }),
      el('div', { class: 'seg' }, [
        segButton('Interface only', 'interface', viewer, true),
        segButton('Glue-proximal', 'proximal', viewer, false),
        segButton('Everything', 'all', viewer, false)
      ])
    ]));

    const decoys = el('select', { class: 'sel' });
    [200, 400, 800, 1600].forEach(function (n) {
      decoys.appendChild(el('option', { value: n, text: n }));
    });
    decoys.value = '800';

    const flavour = el('select', { class: 'sel' });
    flavour.appendChild(el('option', { value: 'configurational', text: 'configurational' }));
    flavour.appendChild(el('option', { value: 'mutational', text: 'mutational' }));

    const run = el('button', { class: 'btn', text: 'Recompute in page' });
    const reset = el('button', { class: 'btn ghost', text: 'Back to cached' });

    const dLab = el('div', { class: 'clab' }, [el('span', { text: 'Decoys per contact' })]);
    help.attach(dLab, 'decoys');
    const fLab = el('div', { class: 'clab' }, [el('span', { text: 'Decoy mode' })]);
    help.attach(fLab, 'flavour');
    const rLab = el('div', { class: 'clab' }, [el('span', { text: 'Run it' })]);
    help.attach(rLab, 'recompute');

    host.appendChild(el('div', { class: 'cgrid' }, [
      el('div', {}, [dLab, decoys]),
      el('div', {}, [fLab, flavour]),
      el('div', {}, [rLab, el('div', { class: 'seg' }, [run, reset])])
    ]));

    host.appendChild(el('p', {
      class: 'hint',
      text:
        'The settings above change how the frustration score is estimated, ' +
        'not what is drawn. More decoys steadies the number; configurational ' +
        'is the stricter of the two modes. Press Recompute to run the engine ' +
        'live rather than reading the cached answer — the "?" next to each ' +
        'explains more.'
    }));
    host.appendChild(el('p', {
      class: 'hint',
      text:
        'Drag to rotate, scroll to zoom, hover a contact to inspect it. The ' +
        'magenta dots are the ligand; it is drawn for reference and never ' +
        'enters the calculation, because the potential has no small-molecule ' +
        'term.'
    }));

    run.addEventListener('click', function () {
      run.disabled = true;
      run.textContent = 'computing…';
      setTimeout(function () {
        const out = recompute(sys, Number(decoys.value), flavour.value);
        viewer.live = out.frustration;
        viewer.draw();
        redrawPlot(out.frustration);
        liveBox.innerHTML = '';
        liveBox.appendChild(el('div', { class: 'livemsg' }, [
          el('strong', { text: 'Recomputed in the page. ' }),
          el('span', {
            text:
              out.ms + ' ms for ' + out.n + ' contacts at ' + decoys.value +
              ' decoys, ' + flavour.value + '. '
          }),
          el('span', {
            class: 'dim',
            text:
              'Largest change against the cached run: ' +
              out.maxDelta.toFixed(4) + '.'
          })
        ]));
        run.disabled = false;
        run.textContent = 'Recompute in page';
      }, 16);
    });

    reset.addEventListener('click', function () {
      viewer.live = null;
      viewer.draw();
      redrawPlot(null);
      liveBox.innerHTML = '';
    });
  }

  function segButton(label, mode, viewer, active) {
    const b = el('button', { class: 'seg-btn' + (active ? ' on' : ''), text: label });
    b.addEventListener('click', function () {
      viewer.mode = mode;
      viewer.draw();
      Array.from(b.parentNode.children).forEach(function (c) {
        c.classList.remove('on');
      });
      b.classList.add('on');
    });
    return b;
  }

  function recompute(sys, nDecoys, flavour) {
    const res = sys.residues;
    const alphabet = GF.data.potential().alphabet;
    const idx = new Map(alphabet.map(function (a, k) { return [a, k]; }));
    const chains = new Map();
    const chainIds = Int32Array.from(res.chain, function (c) {
      if (!chains.has(c)) chains.set(c, chains.size);
      return chains.get(c);
    });

    const site = {
      x: Float64Array.from(res.x),
      y: Float64Array.from(res.y),
      z: Float64Array.from(res.z),
      aa: Int32Array.from(res.one, function (c) {
        const v = idx.get(c);
        return v === undefined ? 0 : v;
      }),
      chainIds: chainIds
    };

    const t0 = performance.now();
    const out = GF.frustration.compute(site, GF.data.potential(), {
      nDecoys: nDecoys, flavour: flavour
    });
    const ms = Math.round(performance.now() - t0);

    let maxDelta = 0;
    const cached = sys.contacts.f;
    if (out.frustration.length === cached.length) {
      for (let k = 0; k < cached.length; k++) {
        const d = Math.abs(out.frustration[k] - cached[k]);
        if (d > maxDelta) maxDelta = d;
      }
    } else {
      maxDelta = NaN;
    }
    return { frustration: out.frustration, ms: ms, n: out.frustration.length, maxDelta: maxDelta };
  }


  /* --- step 6: exit vectors and what they can reach ---------------------
   *
   * This is the bridge from the benchmark to chemistry. Exposure says where
   * you *can* grow; reach says where it is worth growing. They disagree
   * often -- the most solvent-exposed atom on a ligand frequently points
   * away from the interface entirely.
   */

  function buildVectorPanel(host, sys, viewer, legendHost) {
    host.innerHTML = '';
    const ev = sys.exit_vectors;
    const lig = sys.ligand;

    const head = el('div', { class: 'phead' }, [
      el('h3', { text: 'Growth points on ' + sys.assignment.ligand })
    ]);
    help.attach(head, 'exit_vectors');
    host.appendChild(head);

    if (!ev || !ev.vectors.length) {
      host.appendChild(el('p', {
        class: 'lede',
        text:
          'No growth points passed the exposure threshold for this ligand. ' +
          'That happens when the ligand is almost fully enclosed by the ' +
          'anchor, which leaves nothing presented for chemistry to extend ' +
          'into the interface.'
      }));
      return;
    }

    if (lig && lig.sasa_total !== undefined) {
      const facts = el('div', { class: 'facts' });
      facts.appendChild(factCell(
        fmt.fixed(lig.sasa_total, 0) + ' Å²', 'presented face',
        'exposed on the anchor alone', 'presented_face'
      ));
      facts.appendChild(factCell(
        Math.round(lig.buried_by_partner * 100) + '%', 'buried by the partner',
        fmt.fixed(lig.sasa_total_ternary, 0) + ' Å² left in the ternary',
        'presented_face'
      ));
      facts.appendChild(factCell(
        String(ev.vectors.length), 'growth points',
        'exposed heavy atoms with an outward normal', 'exit_vectors'
      ));
      host.appendChild(facts);
    }

    host.appendChild(el('p', {
      class: 'lede',
      text:
        'Each row is an exposed ligand atom and the direction a substituent ' +
        'would project from it. "Reaches" counts frustrated interface ' +
        'contacts within ' + ev.max_reach + ' Å and ' + ev.max_angle_deg +
        '° of that direction. Select a row to isolate it in the structure — ' +
        'switch the view to Exit vectors to see them drawn.'
    }));

    const table = el('table', { class: 'grid vectors' });
    table.appendChild(el('thead', {}, [
      el('tr', {}, [
        el('th', { text: 'Atom' }),
        headCell('Exposure', 'presented_face'),
        headCell('Reaches', 'exit_vectors'),
        headCell('Mean distance', 'exit_vectors'),
        headCell('Addressable', 'addressable')
      ])
    ]));

    const best = Math.max.apply(null, ev.vectors.map(function (v) {
      return v.addressable;
    })) || 1;

    const body = el('tbody');
    ev.vectors.forEach(function (v, idx) {
      const tr = el('tr', { class: 'row vrow' });
      tr.appendChild(el('td', { class: 'mono', text: v.atom_name + ' (' + v.element + ')' }));
      tr.appendChild(el('td', { class: 'num dim', text: fmt.fixed(v.relative_sasa, 2) }));
      tr.appendChild(el('td', { class: 'num', text: String(v.n_reachable) }));
      tr.appendChild(el('td', {
        class: 'num dim',
        text: v.n_reachable ? fmt.fixed(v.mean_reach_distance, 1) + ' Å' : '—'
      }));

      const fill = el('span', { class: 'barfill' });
      fill.style.width = Math.round(100 * v.addressable / best) + '%';
      tr.appendChild(el('td', { class: 'num' }, [
        el('span', { class: 'barval mono', text: fmt.fixed(v.addressable, 2) }),
        el('span', { class: 'bartrack' }, [fill])
      ]));

      tr.addEventListener('click', function () {
        const already = viewer.activeVector === idx;
        viewer.activeVector = already ? -1 : idx;
        if (!already && viewer.style !== 'vectors') {
          viewer.style = 'vectors';
          renderLegend(legendHost, viewer);
        }
        viewer.draw();
        Array.from(body.children).forEach(function (r) { r.classList.remove('sel'); });
        if (!already) tr.classList.add('sel');
      });
      body.appendChild(tr);
    });
    table.appendChild(body);
    host.appendChild(table);

    host.appendChild(el('p', {
      class: 'hint',
      text:
        'Reach is a distance-and-angle feasibility test, not a docking ' +
        'result and not a synthesis plan. It says which frustrated contacts ' +
        'a substituent off this atom could plausibly touch; whether any ' +
        'makeable molecule actually does is the next question, and it needs ' +
        'a chemistry toolkit this build does not carry.'
    }));
  }

  function factCell(value, label, note, slug) {
    const cell = el('div', { class: 'fact' }, [
      el('div', { class: 'fval', text: value })
    ]);
    const lab = el('div', { class: 'flabel' }, [el('span', { text: label })]);
    if (slug) help.attach(lab, slug);
    cell.appendChild(lab);
    cell.appendChild(el('div', { class: 'fnote', text: note }));
    return cell;
  }

  /* --- the comparison plot ---------------------------------------------
   *
   * A strip plot rather than a histogram. The groups here run from about 25
   * to 190 contacts, and 28 bins over 25 points is mostly binning noise --
   * which is what the first version looked like. One dot per contact shows
   * the sample size honestly, needs no vertical scale to misread, and lets
   * you hover an individual contact and find it in the structure.
   */

  function drawPlot(host, sys, live, onPick) {
    host.innerHTML = '';
    const ct = sys.contacts, res = sys.residues;
    const source = live || ct.f;
    const cutoff = sys.proximal_cutoff || 6.0;

    const pts = [];
    for (let k = 0; k < ct.i.length; k++) {
      if (!ct.interface[k]) continue;
      pts.push({
        k: k,
        f: source[k],
        d: Math.min(res.dlig[ct.i[k]], res.dlig[ct.j[k]]),
        prox: !!ct.proximal[k]
      });
    }
    if (!pts.length) return;

    const W = 440, H = 320, padL = 46, padR = 16, padT = 20, padB = 46;
    const fs = pts.map(function (p) { return p.f; });
    const ds = pts.map(function (p) { return p.d; });
    const lo = Math.min(-2, Math.floor(Math.min.apply(null, fs)));
    const hi = Math.max(2, Math.ceil(Math.max.apply(null, fs)));
    const dMax = Math.max(cutoff * 2, Math.ceil(Math.max.apply(null, ds)));

    const X = function (v) {
      return padL + ((v - lo) / (hi - lo)) * (W - padL - padR);
    };
    const Y = function (v) {
      return H - padB - (v / dMax) * (H - padT - padB);
    };

    const ns = 'http://www.w3.org/2000/svg';
    const svg = document.createElementNS(ns, 'svg');
    svg.setAttribute('viewBox', '0 0 ' + W + ' ' + H);
    svg.setAttribute('class', 'histsvg');

    function line(x1, y1, x2, y2, opacity, dash, colour) {
      const l = document.createElementNS(ns, 'line');
      l.setAttribute('x1', x1); l.setAttribute('y1', y1);
      l.setAttribute('x2', x2); l.setAttribute('y2', y2);
      l.setAttribute('stroke', colour || 'currentColor');
      l.setAttribute('stroke-opacity', opacity);
      if (dash) l.setAttribute('stroke-dasharray', dash);
      return l;
    }
    function text(x, y, str, cls, anchor) {
      const t = document.createElementNS(ns, 'text');
      t.setAttribute('x', x); t.setAttribute('y', y);
      t.setAttribute('class', cls);
      if (anchor) t.setAttribute('text-anchor', anchor);
      t.textContent = str;
      return t;
    }

    // The band below the cut-off is the glue-proximal group, so shade it
    // rather than leaving the split to a legend.
    const band = document.createElementNS(ns, 'rect');
    band.setAttribute('x', padL);
    band.setAttribute('y', Y(cutoff));
    band.setAttribute('width', W - padL - padR);
    band.setAttribute('height', Math.max(Y(0) - Y(cutoff), 0));
    band.setAttribute('fill', GF.frustrationColour(-1.4, 0.07));
    svg.appendChild(band);

    // Vertical reference lines: the two conventional frustration cut-offs.
    [[-1, 'frustrated'], [0, null], [0.78, 'minimally']].forEach(function (ref) {
      if (ref[0] < lo || ref[0] > hi) return;
      svg.appendChild(line(X(ref[0]), padT, X(ref[0]), H - padB,
        ref[0] === 0 ? 0.3 : 0.15, '3 3'));
      if (ref[1]) {
        svg.appendChild(text(X(ref[0]) + (ref[0] < 0 ? -4 : 4), padT - 6,
          ref[1], 'axis', ref[0] < 0 ? 'end' : 'start'));
      }
    });

    // The proximal / distal boundary.
    svg.appendChild(line(padL, Y(cutoff), W - padR, Y(cutoff), 0.5, '4 3'));
    svg.appendChild(text(W - padR, Y(cutoff) - 5,
      cutoff + ' Å cut-off', 'axis', 'end'));

    const tip = el('div', { class: 'tip' });

    // Group means, drawn as vertical rules inside each band.
    [
      { rows: pts.filter(function (p) { return p.prox; }),
        y0: Y(cutoff), y1: Y(0), colour: GF.frustrationColour(-1.6, 1) },
      { rows: pts.filter(function (p) { return !p.prox; }),
        y0: Y(dMax), y1: Y(cutoff), colour: 'rgba(150,146,140,1)' }
    ].forEach(function (g) {
      if (!g.rows.length) return;
      const mu = GF.stat.mean(g.rows.map(function (p) { return p.f; }));
      svg.appendChild(line(X(mu), g.y0, X(mu), g.y1, 0.85, null, g.colour));
      svg.appendChild(text(X(mu), (g.y0 + g.y1) / 2 - 4,
        'mean ' + fmt.signed(mu, 2), 'meanlab', 'middle'));
    });

    pts.forEach(function (p) {
      const dot = document.createElementNS(ns, 'circle');
      dot.setAttribute('cx', X(p.f));
      dot.setAttribute('cy', Y(p.d));
      dot.setAttribute('r', 3.2);
      dot.setAttribute('fill', p.prox
        ? GF.frustrationColour(-1.4, 0.85)
        : 'rgba(140,136,130,0.7)');
      dot.setAttribute('class', 'dot');
      dot.addEventListener('mousemove', function (e) {
        const info = sysContactInfo(sys, p.k, p.f);
        tip.innerHTML = '';
        tip.appendChild(el('div', { class: 'tiphead', text: info.pair }));
        tip.appendChild(el('div', { class: 'mono',
          text: 'F ' + fmt.signed(info.f, 2) + ' · ' + info.verdict }));
        tip.appendChild(el('div', { class: 'dim', text: info.detail }));
        tip.style.display = 'block';
        const box = host.getBoundingClientRect();
        tip.style.left = Math.min(e.clientX - box.left + 10, box.width - 196) + 'px';
        tip.style.top = (e.clientY - box.top + 12) + 'px';
        if (onPick) onPick(p.k);
      });
      dot.addEventListener('mouseleave', function () {
        tip.style.display = 'none';
        if (onPick) onPick(-1);
      });
      svg.appendChild(dot);
    });

    // --- axes -------------------------------------------------------------
    svg.appendChild(line(padL, H - padB, W - padR, H - padB, 0.35));
    svg.appendChild(line(padL, padT, padL, H - padB, 0.35));

    for (let v = Math.ceil(lo); v <= hi; v += 1) {
      svg.appendChild(line(X(v), H - padB, X(v), H - padB + 4, 0.35));
      svg.appendChild(text(X(v), H - padB + 15,
        (v > 0 ? '+' : '') + v, 'axis', 'middle'));
    }
    const step = dMax <= 12 ? 2 : dMax <= 30 ? 5 : 10;
    for (let v = 0; v <= dMax; v += step) {
      svg.appendChild(line(padL - 4, Y(v), padL, Y(v), 0.35));
      svg.appendChild(text(padL - 7, Y(v) + 3.5, String(v), 'axis', 'end'));
    }

    svg.appendChild(text(padL + (W - padL - padR) / 2, H - 6,
      'frustration index  —  more frustrated to the left', 'axis', 'middle'));
    const yl = text(0, 0, 'distance to ligand (Å)', 'axis', 'middle');
    yl.setAttribute('transform',
      'translate(13,' + (padT + (H - padT - padB) / 2) + ') rotate(-90)');
    svg.appendChild(yl);

    host.appendChild(svg);
    host.appendChild(tip);
    host.appendChild(el('p', {
      class: 'hint',
      text:
        'Each dot is one interface contact. Horizontal position is how ' +
        'frustrated it is; vertical position is how close it sits to the ' +
        'ligand. Everything in the shaded band is the glue-proximal group. ' +
        'If the hypothesis holds, the shaded band leans left of the rest.'
    }));
  }

  function sysContactInfo(sys, k, f) {
    const ct = sys.contacts, r = sys.residues;
    const a = ct.i[k], b = ct.j[k];
    return {
      pair: r.one[a] + r.seq[a] + ' ' + r.chain[a] + ' — ' + r.one[b] + r.seq[b] + ' ' + r.chain[b],
      f: f,
      verdict: f > 0.78 ? 'minimally frustrated' : f < -1 ? 'highly frustrated' : 'neutral',
      detail: fmt.fixed(ct.d[k], 1) + ' Å apart · ' +
        fmt.fixed(Math.min(r.dlig[a], r.dlig[b]), 1) + ' Å from the ligand'
    };
  }

  /* ===================================================================== */
  /* Method and provenance                                                 */
  /* ===================================================================== */

  function methodPage() {
    const ix = GF.data.index();
    const pot = GF.data.potential();
    const root = el('div', { class: 'page' });
    root.appendChild(el('h1', { text: 'Method and provenance' }));

    root.appendChild(el('p', {
      class: 'lede',
      text:
        'This page exists because the pipeline mixes a live browser engine, ' +
        'a cached sidecar run and a choice of potential. A surprising number ' +
        'is unreadable unless you can tell which part of the stack produced it.'
    }));

    const t = el('table', { class: 'kv wide' });
    [
      ['Arm', ix.arm + ' — interface geometry from the crystal ternary with the glue deleted'],
      ['Potential', pot.name + (pot.published ? ' (published)' : ' (stand-in, not AWSEM)')],
      ['Provenance', pot.provenance],
      ['Decoys per contact', String(ix.n_decoys)],
      ['Proximity cut-off', ix.proximal_cutoff + ' Å from any ligand heavy atom'],
      ['Cache generated', ix.generated],
      ['Systems', String(ix.systems.length)],
      ['r_min / r_max', pot.params.r_min + ' / ' + pot.params.r_max + ' Å'],
      ['r_minII / r_maxII', pot.params.r_minII + ' / ' + pot.params.r_maxII + ' Å'],
      ['eta, eta_sigma', pot.params.eta + ', ' + pot.params.eta_sigma],
      ['rho_0', String(pot.params.rho_0)],
      ['k_contact', pot.params.k_contact + ' kJ/mol']
    ].forEach(function (kv) {
      t.appendChild(el('tr', {}, [
        el('td', { class: 'k', text: kv[0] }),
        el('td', { class: 'v', text: kv[1] })
      ]));
    });
    root.appendChild(t);

    root.appendChild(el('h3', { text: 'How much the potential matters' }));
    root.appendChild(el('p', {
      class: 'lede',
      text:
        'This pipeline ran for a while on a transparent stand-in potential ' +
        'built from hydropathy and charge, because the published tables were ' +
        'not to hand. Swapping them in \u2014 and adding the burial term, ' +
        'which the stand-in lacked \u2014 is the clearest sensitivity test ' +
        'this project has. The direction survived; the magnitude did not.'
    }));
    const pt = el('table', { class: 'grid' });
    pt.appendChild(el('thead', {}, [
      el('tr', {}, [
        el('th', { text: '' }),
        el('th', { class: 'num', text: 'stand-in' }),
        el('th', { class: 'num', text: 'AWSEM 2015' })
      ])
    ]));
    const pb = el('tbody');
    [
      ['Split comparison, direction', '12/14, p = 0.0065', '12/14, p = 0.0065'],
      ['Split comparison, mean effect', '\u22120.332', '\u22120.141'],
      ['Burial-stratified', '12/14, p = 0.0065', '10/14, p = 0.090'],
      ['Gradient, direction', '13/14, p = 0.0009', '13/14, p = 0.0009'],
      ['Gradient, mean rho', '+0.227', '+0.164'],
      ['Individually significant gradients', '11/14', '8/14']
    ].forEach(function (r) {
      pb.appendChild(el('tr', { class: 'row' }, [
        el('td', { text: r[0] }),
        el('td', { class: 'num dim mono', text: r[1] }),
        el('td', { class: 'num mono', text: r[2] })
      ]));
    });
    pt.appendChild(pb);
    root.appendChild(pt);
    root.appendChild(el('p', {
      class: 'hint',
      text:
        'The two potentials agree at r = 0.868 on which systems show the ' +
        'effect, so they measure the same thing with different sensitivity. ' +
        'Reproduce the left column with --gamma builtin on any script.'
    }));

    root.appendChild(el('h3', { text: 'Glossary' }));
    root.appendChild(el('p', {
      class: 'lede',
      text:
        'Every term below also appears behind the "?" next to the number it ' +
        'describes, at three depths.'
    }));
    const gl = el('div', { class: 'glossary' });
    Object.keys(help.TERMS).forEach(function (slug) {
      const term = help.TERMS[slug];
      const title = el('div', { class: 'gtitle' }, [el('span', { text: term.title })]);
      help.attach(title, slug);
      gl.appendChild(el('div', { class: 'gterm' }, [
        title,
        el('div', { class: 'gtext', text: term.concise })
      ]));
    });
    root.appendChild(gl);

    root.appendChild(el('h3', { text: 'What this cannot tell you' }));
    const ul = el('ul', { class: 'caveats' });
    [
      'The oracle arm uses interface geometry from the structure that ' +
      'already contains the answer. It tests whether the hypothesis is true ' +
      'at all, not whether a pipeline could find it blind.',
      'Frustration is defined over amino-acid pairs. The ligand never enters ' +
      'the calculation; it only labels which contacts count as proximal.',
      'Replicate structures of one interface are not independent evidence. ' +
      'The benchmark page groups them for that reason.',
      'The proximity cut-off has now been swept from 4 to 10 Å and the ' +
      'threshold-free rank correlation agrees with it, so the result is not ' +
      'an artefact of that choice. What remains untested is the contact ' +
      'cut-off and the decoy count.',
      'Every system here is a DDB1-based CRL4 complex, so a signal common to ' +
      'all of them may be about that machinery rather than about glue ' +
      'interfaces in general.',
      'The cavity pages are computed from crystal ternary complexes, so the ' +
      'partner is placed perfectly. The error budget says what happens when ' +
      'it is not; the docking that would produce a real placement is not ' +
      'built here.',
      'Cavity electrostatics use formal charges with a distance-dependent ' +
      'dielectric. They are directionally right and quantitatively crude.'
    ].forEach(function (s) { ul.appendChild(el('li', { text: s })); });
    root.appendChild(ul);

    return root;
  }

  GF.pages = { benchmark: benchmarkPage, system: systemPage, method: methodPage };
})(globalThis.GF || (globalThis.GF = {}));
