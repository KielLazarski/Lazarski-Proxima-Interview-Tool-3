/* GlueForge — where to grow, and into what.
 *
 * The last screen. It takes the growth points computed during assessment and
 * says, for each one, how much room is in front of it, whether that room
 * belongs to the partner or back to the anchor, and what kind of group the
 * environment there would suit.
 *
 * Two things this screen tries hard not to do.
 *
 * It does not present a character description as a prescription. The
 * electrostatics behind "this would suit a basic group" are a Coulomb sum
 * over formal charges with a distance-dependent dielectric and no explicit
 * solvent. That is enough to say a region leans positive or negative and
 * nowhere near enough to rank two real compounds, so the wording stays at
 * the level the measurement supports.
 *
 * And it shows what each field value is made of rather than only the value.
 * The best growth point on the CDK12 system sits beside an arginine and
 * wants an acceptor, which looks like a contradiction until you see the
 * contributions: a close arginine pushing one way, a handful of acidic
 * residues pushing the other, and a net that lands between them. The
 * breakdown is the difference between a result and a puzzle.
 */

(function (GF) {
  'use strict';

  const el = GF.el;
  const api = GF.api;

  /* ---- one growth point -------------------------------------------------- */

  function bar(value, max, cls) {
    const pct = Math.max(0, Math.min(100, 100 * value / max));
    return el('span', { class: 'minibar' }, [
      el('span', { class: 'minibar-fill ' + (cls || ''),
                   style: 'width:' + pct + '%' })
    ]);
  }

  function drivers(g) {
    if (!g.field_drivers || !g.field_drivers.length) return null;
    const box = el('div', { class: 'drivers' }, [
      el('div', { class: 'drivers-head', text: 'what makes the field' })
    ]);
    const rows = el('div', { class: 'drivers-rows' });
    g.field_drivers.forEach((d) => {
      rows.appendChild(el('div', {
        class: 'driver' + (d.remainder ? ' rest' : '')
      }, [
        el('span', { class: 'driver-name mono', text: d.residue }),
        el('span', {
          class: 'driver-val mono ' + (d.contribution > 0 ? 'pos' : 'neg'),
          text: (d.contribution > 0 ? '+' : '') + d.contribution.toFixed(2)
        })
      ]));
    });
    box.appendChild(rows);
    box.appendChild(el('div', { class: 'driver net' }, [
      el('span', { class: 'driver-name mono', text: 'net' }),
      el('span', {
        class: 'driver-val mono ' + (g.electrostatic > 0 ? 'pos' : 'neg'),
        text: (g.electrostatic > 0 ? '+' : '') + g.electrostatic.toFixed(2)
      })
    ]));
    return box;
  }

  function pointRow(g, i) {
    const usable = g.reachable_volume >= 8;
    const node = el('div', { class: 'gpoint' + (usable ? '' : ' thin') }, [
      el('div', { class: 'gpoint-head' }, [
        el('span', { class: 'gpoint-rank mono', text: String(i + 1) }),
        el('span', { class: 'gpoint-atom mono', text: g.name }),
        el('span', { class: 'gpoint-el', text: g.element }),
        el('span', { class: 'gpoint-suggest', text: g.suggests })
      ])
    ]);

    const facts = el('div', { class: 'gfacts' }, [
      el('div', {}, [
        el('span', { class: 'glabel', text: 'room in front' }),
        el('span', { class: 'gval mono', text: g.reachable_volume.toFixed(0) + ' Å³' }),
        bar(g.reachable_volume, 120, 'room')
      ]),
      el('div', {}, [
        el('span', { class: 'glabel', text: 'facing the partner' }),
        el('span', { class: 'gval mono', text: (100 * g.toward_partner).toFixed(0) + '%' }),
        bar(g.toward_partner, 1, 'partner')
      ]),
      el('div', {}, [
        el('span', { class: 'glabel', text: 'exposed' }),
        el('span', { class: 'gval mono', text: (100 * g.relative_sasa).toFixed(0) + '%' }),
        bar(g.relative_sasa, 1, 'exposed')
      ])
    ]);
    node.appendChild(facts);

    const d = drivers(g);
    if (d) {
      const fold = el('details', { class: 'fold' }, [
        el('summary', { text:
          'field ' + (g.electrostatic > 0 ? '+' : '') +
          g.electrostatic.toFixed(2) + ' — what makes it' })
      ]);
      fold.appendChild(d);
      if (g.partner_contacts && g.partner_contacts.length) {
        fold.appendChild(el('p', { class: 'fnote', text:
          'Partner residues within reach: ' + g.partner_contacts.join(', ') +
          '. These are the nearest, which is not the same as the ones ' +
          'driving the field above.' }));
      }
      node.appendChild(fold);
    }
    return node;
  }

  /* ---- the screen -------------------------------------------------------- */

  function render(host, run) {
    host.innerHTML = '';
    const a = run.assessment || {};
    const growth = a.growth || {};
    const points = growth.points || [];
    const summary = growth.summary || {};
    const verdict = (a.verdict || {}).verdict;

    host.appendChild(el('div', { class: 'page' }, [
      el('h1', { text: 'Chemistry' }),
      el('p', { class: 'lede' }, [
        'Where the molecule could be grown, and what the space in front of ' +
        'each place would suit. Ranked by room and by whether that room ' +
        'belongs to the partner — growing back into the protein the ligand ' +
        'is already bound to buys no new interface.'
      ])
    ]));

    if (verdict && verdict !== 'pass') {
      host.appendChild(el('div', { class: 'callout warn' }, [
        el('strong', { text: 'The gate did not pass. ' }),
        'These growth points are computed from a cavity the assessment could ' +
        'not distinguish from one a wrong pose would produce. They are ' +
        'arithmetic about this structure either way; whether the structure ' +
        'is right is the open question. ',
        el('a', { href: '#/assess/' + run.id, text: 'Back to the verdict' })
      ]));
    }

    if (!points.length) {
      host.appendChild(el('section', { class: 'panel' }, [
        el('h2', { text: 'Nothing to grow from' }),
        el('p', { class: 'sub', text: summary.note ||
          'No part of the ligand is exposed enough to grow from.' })
      ]));
      return;
    }

    // Toward the partner and back into the anchor are different kinds of
    // opportunity, and a single ranked list invites the wrong read: a point
    // with plenty of room but facing the anchor rises on room alone and
    // looks like a glue opportunity when it is an affinity one.
    const TOWARD = 0.35;
    const toward = points.filter((g) => g.toward_partner >= TOWARD);
    const inward = points.filter((g) => g.toward_partner < TOWARD);

    const sec = el('section', { class: 'panel' }, [
      el('h2', { text: 'Toward the partner' }),
      el('p', { class: 'sub' }, [
        'These face the protein being recruited, so chemistry here adds new ' +
        'interface. That is what makes a glue. ' + (summary.note || '')
      ])
    ]);

    // Two claims live here and only one of them is supported. Letting the
    // supported one lend its credibility to the other is the easiest way to
    // mislead someone with true statements.
    sec.appendChild(el('div', { class: 'callout' }, [
      el('strong', { text: 'Which atoms are listed is validated. The order ' +
                           'they are listed in is not. ' }),
      'Across the eight CDK12 analogues, how exposed an atom is predicts ' +
      'whether that campaign changed it — rank correlation −0.50, and −0.43 ' +
      'with distance from the molecule’s centre partialled out. That is the ' +
      'test for getting onto this list.'
    ]));
    sec.appendChild(el('div', { class: 'callout warn' }, [
      el('strong', { text: 'The ranking points away from where that ' +
                           'campaign actually grew. ' }),
      'Ordering these by room, direction and enclosure correlates +0.37 ' +
      'with being conserved — the tool puts highest the atoms the chemists ' +
      'left alone. The reason looks structural rather than fixable: the ' +
      'atoms they did change sit closer to the partner, in space a rigid ' +
      'cavity does not show as free. Real elaboration reshapes an interface; ' +
      'this measures the room already there. Treat the order as a ' +
      'description of available space, not as advice about what to make ' +
      'next.'
    ]));

    if (toward.length) {
      toward.forEach((g, i) => sec.appendChild(pointRow(g, i)));
    } else {
      sec.appendChild(el('div', { class: 'callout warn' }, [
        el('strong', { text: 'None. ' }),
        'Every exposed atom points back into the anchor. Chemistry here ' +
        'could improve how tightly the molecule holds protein A, but none ' +
        'of it would reach protein B, so none of it would make a glue.'
      ]));
    }
    host.appendChild(sec);

    if (inward.length) {
      const sec2 = el('section', { class: 'panel' }, [
        el('h2', { text: 'Back into the anchor' }),
        el('p', { class: 'sub' }, [
          'These have room too, but it belongs to the protein the molecule ' +
          'is already bound to. Growing here can buy affinity for protein A ' +
          'and will not add interface with protein B. Listed separately ' +
          'because ranking them together makes a roomy inward vector look ' +
          'like an opportunity it is not.'
        ])
      ]);
      inward.forEach((g, i) => sec2.appendChild(pointRow(g, i)));
      host.appendChild(sec2);
    }

    // What this screen cannot do yet, said plainly rather than by omission.
    const ready = api.stageReady('chemistry');
    host.appendChild(el('section', { class: 'panel' }, [
      el('h2', { text: 'Building molecules' }),
      ready
        ? el('p', { class: 'sub', text:
            'RDKit is installed. Attaching fragments at these vectors and ' +
            'scoring the products is the next step and is not built yet.' })
        : el('div', { class: 'callout warn' }, [
            el('strong', { text: 'Needs rdkit. ' }),
            'Everything above is geometry and needs no chemistry toolkit — ' +
            'deciding which vector is worth the attempt does not require ' +
            'building anything. Turning a vector into an actual molecule ' +
            'does. Install it from the ',
            el('a', { href: '#/', text: 'landing page' }), '.'
          ])
    ]));

    GF.wrapTables(host);
  }

  function page(runId) {
    const host = el('div', {});
    host.appendChild(el('div', { class: 'page' }, [
      el('h1', { text: 'One moment' }),
      el('p', { class: 'lede', text: 'Loading the growth points.' })
    ]));
    api.runs.get(runId).then((run) => {
      if (!run.assessment) {
        host.innerHTML = '';
        host.appendChild(el('div', { class: 'page' }, [
          el('h1', { text: 'Not assessed yet' }),
          el('p', { class: 'lede', text:
            'Growth points come from the cavity the assessment measures, so ' +
            'the complex has to be assessed first.' }),
          el('a', { class: 'btn', href: '#/run/' + runId, text: 'Back to the run' })
        ]));
        return;
      }
      render(host, run);
    }).catch((err) => {
      host.innerHTML = '';
      host.appendChild(el('div', { class: 'page' }, [
        el('h1', { text: 'That run is gone' }),
        el('p', { class: 'lede', text: err.message })
      ]));
    });
    return host;
  }

  GF.pages = GF.pages || {};
  GF.pages.chemistry = page;
})(globalThis.GF || (globalThis.GF = {}));
