/* GlueForge — the upload screen.
 *
 * Two structures in, a presented face out. Everything is parsed the moment it
 * lands rather than when docking starts, because the failure modes here are
 * mundane and common — the wrong chain, a ligand that is really a buffer
 * component, a file that is actually a nucleic acid — and every one of them
 * is cheap to catch now and expensive to catch twenty minutes into a run.
 *
 * The presented face is shown before anything is committed to. It is the
 * restraint the docking run will use, so if it looks wrong the time to find
 * out is while it is still free to change.
 */

(function (GF) {
  'use strict';

  const el = GF.el;
  const api = GF.api;
  const fmt = GF.fmt;

  /** Live screen state. One run at a time; the server holds the truth. */
  const S = {
    run: null,        // the run.json payload
    busy: {},         // slot -> true while uploading
    error: {},        // slot -> message
    faceError: null,
    faceBusy: false,
    chainBusy: false,
    chainError: null,
    examples: null,
    examplesAsked: false,
    examplesNote: null,
    examplesError: null,
    exampleBusy: false
  };

  /** Keep only these chains of a slot, and redraw.
   *
   *  The server rewrites protein A to the selection and discards anything
   *  measured against the old one -- the face, the verdict, a docking run.
   *  That is said out loud when it happens rather than left for someone to
   *  notice a number went missing.
   */
  function applyChains(slot, keep) {
    if (S.chainBusy) return;
    S.chainBusy = true;
    S.chainError = null;
    render();
    api.json('/api/run/' + S.run.id + '/chains',
             { slot: slot, keep: keep }).then((res) => {
      S.chainBusy = false;
      if (S.run && S.run.slots) S.run.slots[slot] = res.structure;
      (res.invalidated || []).forEach((k) => {
        if (S.run) delete S.run[k];
      });
      S.lastInvalidated = res.invalidated || [];
      return api.runs.get(S.run.id).then((r) => { S.run = r; });
    }).then(render).catch((err) => {
      S.chainBusy = false;
      S.chainError = err && err.message ? err.message : String(err);
      render();
    });
  }

  function reset() {
    S.run = null;
    S.busy = {};
    S.error = {};
    S.faceError = null;
    S.faceBusy = false;
    S.assessBusy = false;
    S.assessError = null;
    S.roles = null;
  }

  /* ---- drop zones ------------------------------------------------------ */

  const SLOT_COPY = {
    A: {
      title: 'Protein A',
      sub: 'the one holding the small molecule',
      why: 'The ligand has to already be bound here. Its exposed face is ' +
           'what aims the docking run and what chemistry grows from.'
    },
    B: {
      title: 'Protein B',
      sub: 'the one to be recruited',
      why: 'No ligand needed. This is the protein the glue has to pull in.'
    }
  };

  function dropZone(slot) {
    const copy = SLOT_COPY[slot];
    const have = S.run && S.run.slots && S.run.slots[slot];
    const zone = el('div', {
      class: 'drop' + (have ? ' filled' : '') + (S.busy[slot] ? ' busy' : '')
    });

    const input = el('input', {
      type: 'file',
      accept: '.pdb,.cif,.mmcif,.ent',
      style: 'display:none'
    });
    input.addEventListener('change', () => {
      if (input.files && input.files[0]) send(slot, input.files[0]);
    });

    const head = el('div', { class: 'drop-head' }, [
      el('span', { class: 'drop-slot', text: slot }),
      el('div', {}, [
        el('div', { class: 'drop-title', text: copy.title }),
        el('div', { class: 'drop-sub', text: copy.sub })
      ])
    ]);
    zone.appendChild(head);

    if (S.busy[slot]) {
      zone.appendChild(el('p', { class: 'drop-state', text: 'Reading…' }));
    } else if (have) {
      zone.appendChild(summary(slot, have));
    } else {
      zone.appendChild(el('p', { class: 'drop-why', text: copy.why }));
      zone.appendChild(el('div', { class: 'drop-cta' }, [
        el('button', {
          class: 'btn',
          text: 'Choose a file',
          onclick: () => input.click()
        }),
        el('span', { class: 'drop-or', text: 'or drop a .pdb or .cif here' })
      ]));
    }

    if (S.error[slot]) {
      zone.appendChild(el('p', { class: 'drop-err', text: S.error[slot] }));
    }

    zone.appendChild(input);

    ['dragenter', 'dragover'].forEach((ev) =>
      zone.addEventListener(ev, (e) => {
        e.preventDefault();
        zone.classList.add('over');
      }));
    ['dragleave', 'drop'].forEach((ev) =>
      zone.addEventListener(ev, (e) => {
        e.preventDefault();
        zone.classList.remove('over');
      }));
    zone.addEventListener('drop', (e) => {
      const f = e.dataTransfer && e.dataTransfer.files && e.dataTransfer.files[0];
      if (f) send(slot, f);
    });

    return zone;
  }

  /* ---- what the parse found -------------------------------------------- */

  function summary(slot, st) {
    const box = el('div', { class: 'parsed' });

    box.appendChild(el('div', { class: 'parsed-file' }, [
      el('span', { class: 'mono', text: st.name }),
      el('span', { class: 'dim', text: fmt.bytes(st.bytes) }),
      el('button', {
        class: 'btn ghost small',
        text: 'Replace',
        onclick: () => {
          delete S.error[slot];
          if (S.run && S.run.slots) delete S.run.slots[slot];
          render();
        }
      })
    ]));

    const facts = el('div', { class: 'facts' });
    const fact = (label, value, note, cls) => facts.appendChild(
      el('div', { class: 'fact' }, [
        el('div', { class: 'flabel', text: label }),
        el('div', { class: 'fval' + (cls ? ' ' + cls : ''), text: value }),
        note ? el('div', { class: 'fnote', text: note }) : null
      ]));
    fact('Entry', st.entry || '—', null, 'name');
    fact('Chains', String(st.n_chains));
    fact('Residues', String(st.n_residues));
    fact('Hetero groups', String((st.ligands || []).length),
         st.waters_removed ? st.waters_removed + ' waters dropped' : null);
    box.appendChild(facts);

    // Chains, biggest first, with how close each comes to the molecule and
    // a tick to keep it. A deposited structure is whatever crystallised:
    // CRBN arrives bound to DDB1, which sits 34 Å from the glue site and
    // belongs to a different question. Without a way to drop it, the file
    // cannot be used at all.
    const kept = (st.kept_chains && st.kept_chains.slice()) ||
      (st.chains || []).map((c) => c.id);
    const ch = el('table', { class: 'grid chains' });
    ch.appendChild(el('thead', {}, [el('tr', {}, [
      el('th', { text: slot === 'A' ? 'Use' : '' }),
      el('th', { text: 'Chain' }),
      el('th', { text: 'Residues' }),
      el('th', { text: 'From the ligand' }),
      el('th', { text: 'Description' })
    ])]));
    const cb = el('tbody');
    (st.chains || []).forEach((c) => {
      const on = kept.indexOf(c.id) >= 0;
      const tick = el('input', { type: 'checkbox', checked: on ? '' : null });
      tick.addEventListener('change', () => {
        const next = (st.chains || []).filter((x) => {
          const box2 = cb.querySelector('tr[data-chain="' + x.id + '"] input');
          return box2 ? box2.checked : false;
        }).map((x) => x.id);
        applyChains(slot, next);
      });
      const near = c.ligand_distance === null || c.ligand_distance === undefined
        ? '—'
        : c.ligand_distance.toFixed(1) + ' Å' +
          (c.ligand_contacts ? ' · touches it' : '');
      cb.appendChild(el('tr', { 'data-chain': c.id,
                                class: on ? '' : 'chain-off' }, [
        el('td', {}, slot === 'A' ? [tick] : []),
        el('td', { class: 'mono', text: c.id }),
        el('td', { class: 'val', text: String(c.residues) }),
        el('td', { class: 'val', text: near }),
        el('td', { text: c.name || '—' })
      ]));
    });
    ch.appendChild(cb);
    box.appendChild(ch);

    if (slot === 'A') {
      box.appendChild(el('p', { class: 'fnote', text:
        'Untick a chain to leave it out of protein A. A chain that sits far '
        + 'from the ligand is usually a crystallisation partner rather than '
        + 'part of the question — though one that packs against the site is '
        + 'worth keeping, because a docked protein must not land where it '
        + 'already is.' }));
      if (S.chainBusy) {
        box.appendChild(el('p', { class: 'sub', text: 'Rewriting protein A…' }));
      }

      // The chains taken out of protein A are usually the obvious protein
      // B -- dropping them is what makes room to dock something, and the
      // thing most worth docking first is the protein that was there.
      const dropped = st.dropped_chains || [];
      const haveB = !!((S.run.slots || {}).B);
      if (dropped.length && !haveB && !S.chainBusy) {
        const use = el('button', { class: 'btn',
          text: 'Use chain ' + dropped.join(', ') + ' as protein B' });
        use.addEventListener('click', () => {
          applyChains('A', (st.kept_chains || []).slice(), true);
        });
        box.appendChild(el('div', { class: 'callout note' }, [
          el('p', { class: 'sub', text:
            'Chain ' + dropped.join(', ') + ' is out of protein A. Docking '
            + 'it back is the one case where you already know the answer, '
            + 'so it is worth doing once to see what the search finds.' }),
          use
        ]));
      }
      if (S.chainError) {
        box.appendChild(el('div', { class: 'callout warn' }, [
          el('strong', { text: 'Could not change the chains. ' }), S.chainError
        ]));
      }
      box.appendChild(ligandPicker(st));
    }
    return box;
  }

  /* ---- choosing the ligand --------------------------------------------- */

  function ligandPicker(st) {
    const ligs = st.ligands || [];
    if (!ligs.length) {
      return el('div', { class: 'callout bad' }, [
        el('strong', { text: 'No ligand in this file. ' }),
        'The method starts from a small molecule already bound to protein A, ' +
        'so without one there is nothing to grow from and nothing to aim a ' +
        'docking run at. If the ligand is in a separate file you will need to ' +
        'place it into the pocket first.'
      ]);
    }

    const chosen = (S.run && S.run.ligand) || null;
    const box = el('div', { class: 'ligands' }, [
      el('div', { class: 'llabel', text: 'The bound molecule' })
    ]);

    if (ligs.length > 1) {
      box.appendChild(el('p', { class: 'sub' }, [
        ligs.length + ' hetero groups here, largest first. Obvious buffer and ' +
        'cryoprotectant components are already filtered out, but a ' +
        'crystallisation additive can still look like a ligand — pick the one ' +
        'you mean.'
      ]));
    }

    const list = el('div', { class: 'lrows' });
    ligs.forEach((l) => {
      const on = chosen && chosen.comp_id === l.comp_id &&
                 (chosen.chain == null || chosen.chain === l.chain) &&
                 (chosen.seq_id == null || chosen.seq_id === l.seq_id);
      list.appendChild(el('button', {
        class: 'lig' + (on ? ' on' : ''),
        onclick: () => chooseLigand(l)
      }, [
        el('span', { class: 'lig-id mono', text: l.comp_id }),
        el('span', { class: 'lig-meta', text:
          l.heavy_atoms + ' heavy atoms · chain ' + l.chain + ' · ' + l.seq_id }),
        on ? el('span', { class: 'lig-on', text: 'measuring this one' }) : null
      ]));
    });
    box.appendChild(list);

    if (S.faceBusy) {
      box.appendChild(el('p', { class: 'drop-state', text:
        'Measuring the presented face…' }));
    }
    if (S.faceError) {
      box.appendChild(el('p', { class: 'drop-err', text: S.faceError }));
    }
    return box;
  }

  /* ---- the restraint --------------------------------------------------- */

  /** The rim, split the way LightDock treats it.
   *
   *  Active restraints must be satisfied; passive ones are permitted. Showing
   *  them as one undifferentiated list would hide the only decision being
   *  made here. A residue within van der Waals reach of an exposed ligand
   *  atom is certainly lining the opening; one seven angstroms away might
   *  merely be nearby, and making it mandatory over-constrains the run.
   */
  function rimPanel(face) {
    const byLabel = {};
    (face.rim || []).forEach((r) => { byLabel[r.label] = r; });
    const rst = face.restraints || null;

    const chip = (r, active) => el('span', {
      class: 'chip mono' + (active ? ' active' : ''),
      title: (r
        ? r.min_distance + ' Å from an exposed ligand atom · ' +
          r.sasa + ' Å² exposed'
        : '') + (active ? ' · must be at the interface'
                        : ' · allowed at the interface'),
      text: r ? r.comp_id + r.seq_id : '?'
    });

    const box = el('div', { class: 'rim' });

    if (!rst) {
      box.appendChild(el('div', { class: 'llabel', text:
        'Docking will be restrained to these residues, nearest first' }));
      box.appendChild(el('div', { class: 'chips' },
        (face.rim || []).slice(0, 40).map((r) => chip(r, false))));
      return box;
    }

    box.appendChild(el('div', { class: 'llabel', text: 'The docking restraint' }));
    box.appendChild(el('p', { class: 'sub' }, [
      'These go to the docking run as residue restraints. It discards the ' +
      'starting regions that are nowhere near them, which is the whole of ' +
      'what this tool contributes to the search — the restraint comes from ' +
      'the ligand, which is information ordinary docking does not have.'
    ]));

    if (rst.n_active) {
      box.appendChild(el('div', { class: 'rimgroup' }, [
        el('div', { class: 'rimhead' }, [
          el('strong', { text: rst.n_active + ' active' }),
          el('span', { class: 'dim', text:
            'furthest along the face’s outward normal — the pose has to ' +
            'touch these' })
        ]),
        el('div', { class: 'chips' },
          rst.active.map((l) => chip(byLabel[l], true)))
      ]));
      box.appendChild(el('p', { class: 'fnote' }, [
        'Chosen by how far each residue sits along the direction the partner ' +
        'must approach from, not by how close it is to the ligand. Across ' +
        'the 14 validated systems, closeness to the ligand does not predict ' +
        'which rim residues the real partner touches (r = −0.07) and this ' +
        'does (r = +0.57).'
      ]));
    } else {
      box.appendChild(el('div', { class: 'callout warn' }, [
        el('strong', { text: 'Nothing is required, only permitted. ' }),
        'This face has no single approach direction, so there is no honest ' +
        'way to pick residues the partner must touch. The whole rim goes in ' +
        'as permitted, which still narrows where docking looks but narrows ' +
        'it less.'
      ]));
    }

    if (rst.n_passive) {
      box.appendChild(el('div', { class: 'rimgroup' }, [
        el('div', { class: 'rimhead' }, [
          el('strong', { text: rst.n_passive + ' passive' }),
          el('span', { class: 'dim', text:
            'further out — permitted at the interface, not required' })
        ]),
        el('div', { class: 'chips' },
          rst.passive.slice(0, 40).map((l) => chip(byLabel[l], false)))
      ]));
      if (rst.passive.length > 40) {
        box.appendChild(el('p', { class: 'fnote', text:
          'Showing 40 of ' + rst.passive.length + '.' }));
      }
    }
    return box;
  }

  /* ---- the presented face ---------------------------------------------- */

  function facePanel() {
    const face = S.run && S.run.face;
    if (!face) return null;

    const heading = el('h2', { text: 'The presented face' });
    if (GF.help) GF.help.attach(heading, 'presented_face');

    const panel = el('section', { class: 'panel' }, [
      heading,
      el('p', { class: 'note' }, [
        'Measured against protein A alone, which is the only way the number ' +
        'means anything: in an assembled complex the partner is already ' +
        'sitting on this very surface, so every atom would read as buried.'
      ])
    ]);

    const facts = el('div', { class: 'facts' });
    const fact = (label, value, note) => facts.appendChild(
      el('div', { class: 'fact' }, [
        el('div', { class: 'flabel', text: label }),
        el('div', { class: 'fval', text: value }),
        note ? el('div', { class: 'fnote', text: note }) : null
      ]));

    fact('Exposed area', fmt.fixed(face.presented_area, 1) + ' Å²',
         'across ' + face.n_exposed + ' heavy atoms');
    fact('Buried', (100 * face.buried_fraction).toFixed(1) + '%',
         'of the molecule’s free surface');
    fact('Rim residues', String(face.n_rim),
         'protein A residues lining the opening');
    if (face.directionality !== undefined) {
      fact('Face shape', face.directionality.toFixed(2),
           face.directional ? 'points one way — a pocket'
                            : 'wraps the ligand — a groove');
    }
    fact('Anchor chain', (face.anchor_chains || []).join(', ') || '—');
    panel.appendChild(facts);

    panel.appendChild(el('p', { class: 'sub hyd' }, [
      'Hydrogens are discarded whether or not the file has them, so this ' +
      'area is comparable between structures. Buried fraction is the more ' +
      'robust of the two numbers — it is stable to about ±0.05 across the ' +
      'two conventions, where absolute area is not.'
    ]));

    (face.warnings || []).forEach((w) => {
      panel.appendChild(el('div', { class: 'callout warn', text: w }));
    });

    if (face.n_exposed === 0) {
      panel.appendChild(el('div', { class: 'callout bad' }, [
        el('strong', { text: 'Nothing to work with. ' }),
        'A fully enclosed ligand cannot be touched by a second protein, so ' +
        'this pair has no glue to make regardless of what docking returns.'
      ]));
      return panel;
    }

    // The rim is the restraint. Showing it as the actual residue list, split
    // the way a docking run will treat it, is the only way to tell whether it
    // landed on the right surface — and the split matters, because active and
    // passive restraints do different things.
    panel.appendChild(rimPanel(face));

    const exposed = el('table', { class: 'grid' });
    exposed.appendChild(el('thead', {}, [el('tr', {}, [
      el('th', { text: 'Atom' }),
      el('th', { text: 'Element' }),
      el('th', { text: 'Exposed Å²' }),
      el('th', { text: 'Share of free' })
    ])]));
    const eb = el('tbody');
    (face.exposed || []).forEach((a) => {
      eb.appendChild(el('tr', {}, [
        el('td', { class: 'mono', text: a.name }),
        el('td', { class: 'mono', text: a.element }),
        el('td', { class: 'val', text: fmt.fixed(a.sasa, 1) }),
        el('td', { class: 'val', text: (100 * a.relative_sasa).toFixed(0) + '%' })
      ]));
    });
    exposed.appendChild(eb);
    panel.appendChild(el('details', { class: 'fold' }, [
      el('summary', { text: 'The exposed atoms, one per row' }),
      exposed
    ]));

    return panel;
  }

  /* ---- worked examples --------------------------------------------------- */

  /** Solved complexes that ship with the tool, with their own measurements.
   *
   *  Every number shown here was produced by running the same code the tool
   *  runs, not written alongside it. An example that described itself in
   *  terms its own measurement contradicted would be worse than no example,
   *  because it is the first thing someone reads and the thing they
   *  calibrate every later number against.
   */
  function examplePanel() {
    if (S.examplesError) {
      return el('div', { class: 'callout warn' }, [
        el('strong', { text: 'Could not load the examples. ' }),
        S.examplesError
      ]);
    }
    if (!S.examples) {
      if (!S.examplesAsked) {
        S.examplesAsked = true;
        api.call('/api/examples').then((d) => {
          S.examples = d.examples || [];
          S.examplesNote = d.note || null;
          render();
        }).catch((err) => {
          // A dead server is not an examples problem, and saying so inside
          // this panel while the rest of the page carries on looking live
          // is the wrong thing to show. `api.call` has already marked the
          // session offline; re-routing puts the whole page into the state
          // that explains it.
          if (err && err.status === 0) { GF.route(); return; }
          S.examplesError = err && err.message ? err.message : String(err);
          render();
        });
      }
      return null;
    }
    if (!S.examples.length) {
      return S.examplesNote
        ? el('p', { class: 'fnote', text: S.examplesNote })
        : null;
    }

    const sec = el('section', { class: 'panel' }, [
      el('h2', { text: 'Or start from a solved complex' }),
      el('p', { class: 'sub' }, [
        'Five ternary complexes that ship with GlueForge, already split '
        + 'into protein A with its ligand and protein B on its own. Every '
        + 'number below was measured by the tool itself. They are different '
        + 'from each other on purpose: the benchmark holds eight CDK12 '
        + 'entries and four DCAF15 ones, and shipping them all would be a '
        + 'longer list showing less.'
      ])
    ]);

    const grid = el('div', { class: 'examples' });
    S.examples.forEach((x) => {
      const bad = x.verdict !== 'pass';
      const card = el('div', { class: 'example' + (bad ? ' example-warn' : '') });
      card.appendChild(el('div', { class: 'ex-head' }, [
        el('code', { class: 'ex-entry', text: x.entry }),
        el('span', { class: 'ex-lig mono', text: x.ligand }),
        el('span', { class: 'ex-verdict v-' + (bad ? 'warn' : 'ok'),
                     text: x.verdict || '—' })
      ]));
      card.appendChild(el('div', { class: 'ex-title', text: x.title }));
      card.appendChild(el('p', { class: 'ex-shows', text: x.shows }));
      card.appendChild(el('div', { class: 'ex-nums mono', text:
        x.presented_area + ' Å² face · ' +
        Math.round(100 * x.buried_fraction) + '% buried · ' +
        Math.round(x.cavity_volume) + ' Å³ cavity' }));
      const btn = el('button', { class: 'btn', text: 'Load this' });
      btn.addEventListener('click', () => loadExample(x.id));
      card.appendChild(btn);
      grid.appendChild(card);
    });
    sec.appendChild(grid);
    if (S.exampleBusy) {
      sec.appendChild(el('p', { class: 'sub', text: 'Loading…' }));
    }
    return sec;
  }

  function loadExample(id) {
    if (S.exampleBusy) return;
    S.exampleBusy = true;
    S.examplesError = null;
    render();
    api.json('/api/run/' + S.run.id + '/example', { example: id })
      .then(() => api.runs.get(S.run.id))
      .then((r) => { S.run = r; S.exampleBusy = false; render(); })
      .catch((err) => {
        S.exampleBusy = false;
        S.examplesError = err && err.message ? err.message : String(err);
        render();
      });
  }

  /* ---- already a complex ----------------------------------------------- */

  /** When protein A's file already holds both proteins.
   *
   *  Asking for a second upload that is already on disk is the kind of thing
   *  that makes a tool feel like paperwork. If the ligand touches two chains,
   *  this file is a ternary complex and can go straight to the gate — the
   *  only thing missing is which chain is the anchor, and that is a question
   *  with a sensible default.
   */
  function complexPanel() {
    const face = S.run && S.run.face;
    if (!face) return null;
    const chains = face.contacting_chains || [];
    if (chains.length < 2) return null;

    const roles = S.roles || { anchor: chains[0], partner: chains[1] };
    const panel = el('section', { class: 'panel' }, [
      el('h2', { text: 'This file already holds two proteins at the ligand' }),
      el('p', { class: 'sub' }, [
        'The ligand touches ' + chains.length + ' chains, so this is an ' +
        'assembled complex. Chain ' + chains[0] + ' is taken as the anchor ' +
        'because it buries most of the ligand. That choice decides which ' +
        'surface the face above is measured against, so it is not cosmetic: ' +
        'change it and the face is re-measured.'
      ]),
      // Being an assembled complex used to end the matter, and that was
      // wrong. A solved ternary complex is also the best possible starting
      // point for docking a *different* partner at the same site -- the
      // ligand is already in its pocket and the face is already measured.
      // The only thing in the way is the protein that is currently there.
      el('p', { class: 'sub' }, [
        'You can assess it as it stands, or drop chain ' +
        chains[1] + ' above and dock something else against the same face. ' +
        'Both are reasonable; the first asks whether this complex is a glue, ' +
        'the second asks what else that face could recruit.'
      ])
    ]);

    const pick = (which, other) => {
      const sel = el('select', { class: 'sel' });
      chains.forEach((c) => {
        const o = el('option', { value: c, text: 'chain ' + c });
        if (roles[which] === c) o.selected = true;
        sel.appendChild(o);
      });
      sel.addEventListener('change', () => {
        const next = Object.assign({}, roles);
        next[which] = sel.value;
        if (next[other] === sel.value) {
          next[other] = chains.find((c) => c !== sel.value) || next[other];
        }
        S.roles = next;
        // Changing the anchor changes which surface the face is measured
        // against, so the face shown above has to be recomputed. Leaving a
        // stale number on screen next to a changed choice is how people end
        // up quoting the wrong one.
        remeasure(next.anchor);
      });
      return sel;
    };

    panel.appendChild(el('div', { class: 'roles' }, [
      el('label', {}, [
        el('span', { class: 'llabel', text: 'Anchor' }),
        el('span', { class: 'dim', text: 'holds the ligand' }),
        pick('anchor', 'partner')
      ]),
      el('label', {}, [
        el('span', { class: 'llabel', text: 'Partner' }),
        el('span', { class: 'dim', text: 'the recruited one' }),
        pick('partner', 'anchor')
      ])
    ]));

    if (S.assessError) {
      panel.appendChild(el('p', { class: 'drop-err', text: S.assessError }));
    }
    panel.appendChild(el('button', {
      class: 'btn primary',
      text: S.assessBusy ? 'Measuring…' : 'Assess this complex',
      disabled: S.assessBusy ? 'disabled' : null,
      onclick: () => assess(roles)
    }));
    return panel;
  }

  /** Re-measure the presented face against a different anchor chain. */
  async function remeasure(anchor) {
    const lig = (S.run && S.run.ligand) || null;
    if (!lig) { render(); return; }
    S.faceBusy = true;
    S.faceError = null;
    render();
    try {
      await api.runs.face(S.run.id, {
        ligand: lig.comp_id,
        ligand_chain: lig.chain,
        ligand_seq: lig.seq_id,
        anchor_chains: [anchor]
      });
      S.run = await api.runs.get(S.run.id);
    } catch (err) {
      S.faceError = err.message;
    }
    S.faceBusy = false;
    render();
  }

  async function assess(roles) {
    S.assessBusy = true;
    S.assessError = null;
    render();
    try {
      await api.json('/api/run/' + S.run.id + '/assess', {
        slot: 'A',
        anchor_chains: [roles.anchor],
        partner_chains: [roles.partner]
      });
      S.assessBusy = false;
      location.hash = '#/assess/' + S.run.id;
      return;
    } catch (err) {
      S.assessError = err.message;
    }
    S.assessBusy = false;
    render();
  }

  /* ---- next step ------------------------------------------------------- */

  function nextStep() {
    const run = S.run;
    const haveA = run && run.slots && run.slots.A;
    const haveB = run && run.slots && run.slots.B;
    const face = run && run.face;
    const dock = api.stage('dock');

    const reasons = [];
    if (!haveA) reasons.push('protein A');
    if (!haveB) reasons.push('protein B');
    if (haveA && !face) reasons.push('a ligand to measure');

    const panel = el('section', { class: 'panel' }, [
      el('h2', { text: 'Dock them' })
    ]);

    if (reasons.length) {
      panel.appendChild(el('p', { class: 'sub', text:
        'Still needed: ' + reasons.join(', ') + '.' }));
      panel.appendChild(el('button', {
        class: 'btn primary', text: 'Dock', disabled: 'disabled'
      }));
      return panel;
    }

    if (dock && !dock.ready) {
      panel.appendChild(el('div', { class: 'callout warn' }, [
        'Docking needs ',
        el('code', { text: dock.missing.join(', ') }),
        '. Install it from the ',
        el('a', { href: '#/', text: 'landing page' }),
        ' and come back — the two structures and the face stay where they are.'
      ]));
      // The screen is still worth reaching: it shows what the run would cost
      // and what restraint it would use, which is how you decide whether
      // installing is worth it.
      panel.appendChild(el('a', {
        class: 'btn ghost', href: '#/dock/' + run.id,
        text: 'See what docking would do'
      }));
      return panel;
    }

    panel.appendChild(el('p', { class: 'sub' }, [
      'The restraint is the rim above. A docking run is scored on how much of ' +
      'the presented face it buries, not on how good the interface looks ' +
      'overall, because an interface that does not cover the ligand is not a ' +
      'glue interface however favourable it is.'
    ]));
    panel.appendChild(el('a', {
      class: 'btn primary', href: '#/dock/' + run.id, text: 'Dock'
    }));
    return panel;
  }

  /* ---- server calls ---------------------------------------------------- */

  async function ensureRun() {
    if (S.run) return S.run;
    S.run = await api.runs.create();
    return S.run;
  }

  async function send(slot, file) {
    delete S.error[slot];
    S.busy[slot] = true;
    render();
    try {
      const run = await ensureRun();
      const res = await api.runs.upload(run.id, slot, file);
      S.run = await api.runs.get(run.id);
      // One obvious ligand and no ambiguity: measure it without being asked.
      if (slot === 'A') {
        const ligs = (res.structure && res.structure.ligands) || [];
        if (ligs.length === 1) {
          S.busy[slot] = false;
          render();
          await chooseLigand(ligs[0]);
          return;
        }
      }
    } catch (err) {
      S.error[slot] = err.message;
      if (err.missing) {
        S.error[slot] += ' — install it from the landing page.';
      }
    }
    S.busy[slot] = false;
    render();
  }

  async function chooseLigand(l) {
    S.faceError = null;
    S.faceBusy = true;
    render();
    try {
      await api.runs.face(S.run.id, {
        ligand: l.comp_id,
        ligand_chain: l.chain,
        ligand_seq: l.seq_id
      });
      S.run = await api.runs.get(S.run.id);
    } catch (err) {
      S.faceError = err.message;
    }
    S.faceBusy = false;
    render();
  }

  /* ---- assembly -------------------------------------------------------- */

  let host = null;

  function render() {
    if (!host) return;
    host.innerHTML = '';

    host.appendChild(el('div', { class: 'page' }, [
      el('h1', { text: 'Screen a pair' }),
      el('p', { class: 'lede' }, [
        'Protein A with its ligand, protein B on its own. Both are read as ' +
        'soon as they land so anything wrong shows up here rather than later.'
      ]),
      S.run
        ? el('p', { class: 'fnote mono', text: 'run ' + S.run.id })
        : null
    ]));

    host.appendChild(el('div', { class: 'drops' }, [
      dropZone('A'), dropZone('B')
    ]));

    // Offered only before anything is uploaded. Once there are files on the
    // run, a list of examples is a way to lose them by accident.
    const slots = (S.run && S.run.slots) || {};
    if (!slots.A && !slots.B) {
      const ex = examplePanel();
      if (ex) host.appendChild(ex);
    }

    const face = facePanel();
    if (face) host.appendChild(face);
    const cx = complexPanel();
    if (cx) host.appendChild(cx);
    // A file that is already a complex does not need a second upload, so the
    // docking step is not the next thing to offer.
    if (!cx) host.appendChild(nextStep());
    GF.wrapTables(host);
  }

  function upload(runId) {
    host = el('div', {});
    if (runId && (!S.run || S.run.id !== runId)) {
      reset();
      api.runs.get(runId).then((r) => { S.run = r; render(); })
        .catch((err) => {
          host.innerHTML = '';
          host.appendChild(el('div', { class: 'page' }, [
            el('h1', { text: 'That run is gone' }),
            el('p', { class: 'lede', text: err.message }),
            el('a', { class: 'btn', href: '#/new', text: 'Start a new one' })
          ]));
        });
    } else {
      if (!runId) reset();
      render();
    }
    return host;
  }

  /** The assessment screen for a run that has one. */
  function assessPage(runId) {
    const host = el('div', {});
    const draw = (run) => {
      host.innerHTML = '';
      if (!run || !run.assessment) {
        host.appendChild(el('div', { class: 'page' }, [
          el('h1', { text: 'Nothing assessed yet' }),
          el('p', { class: 'lede', text:
            'This run has no assessment on it. Go back and measure the ' +
            'complex first.' }),
          el('a', { class: 'btn', href: '#/run/' + runId, text: 'Back to the run' })
        ]));
        return;
      }
      host.appendChild(GF.assess.panel(run.assessment, () => {
        location.hash = '#/chemistry/' + runId;
      }, { runId: runId }));

      // The complex the verdict was reached about, and the cavity the
      // verdict is largely about. Fetched after the panel is on screen
      // rather than before it: the assessment is the thing the user came
      // for, and making them wait on a structure export to read it would be
      // the wrong way round.
      const slot = el('div', {});
      host.appendChild(slot);
      GF.view.forRun(runId).then((blob) => {
        slot.appendChild(GF.view.panel(blob, {
          title: 'The complex as assessed' }).node);
        slot.appendChild(GF.view.cavityPanel(blob).node);
      }).catch((err) => {
        slot.appendChild(el('p', { class: 'fnote', text:
          'The structure could not be drawn: ' +
          (err && err.message ? err.message : String(err)) }));
      });

      host.appendChild(el('p', { class: 'fnote' }, [
        'run ', el('code', { text: runId }), ' · ',
        el('a', { href: '#/run/' + runId, text: 'back to the structures' })
      ]));
    };

    if (S.run && S.run.id === runId && S.run.assessment) {
      draw(S.run);
    } else {
      host.appendChild(el('div', { class: 'page' }, [
        el('h1', { text: 'One moment' }),
        el('p', { class: 'lede', text: 'Loading the assessment.' })
      ]));
      api.runs.get(runId).then((r) => { S.run = r; draw(r); })
        .catch((err) => {
          host.innerHTML = '';
          host.appendChild(el('div', { class: 'page' }, [
            el('h1', { text: 'That run is gone' }),
            el('p', { class: 'lede', text: err.message })
          ]));
        });
    }
    return host;
  }

  GF.pages = GF.pages || {};
  GF.pages.upload = upload;
  GF.pages.assess = assessPage;
})(globalThis.GF || (globalThis.GF = {}));
