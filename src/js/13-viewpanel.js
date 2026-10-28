/* GlueForge — a 3D panel for live runs.
 *
 * The validation browser has had a viewer since the frustration work, but it
 * is built into that page and wired to the cached benchmark data. A docking
 * result wants the same picture and almost none of the surrounding
 * apparatus: no frustration plot, no decoy histogram, no live potential
 * recompute — those read measurements a run does not have.
 *
 * So this is the viewer with the structural half kept and the analytical
 * half dropped: rotate, zoom, pick a style, see which residues are at the
 * interface and where the ligand sits. The same `GF.Viewer`, which already
 * knows to offer only the styles a given system can be drawn in.
 *
 * It also carries the cavity, because for this tool the cavity *is* the
 * result — a glue has to fill a pocket at the new interface, and a pose with
 * no pocket at it is not a lead however well it docked.
 */

(function (GF) {
  'use strict';

  const el = GF.el;
  const api = GF.api;

  /** Colour ramps, matching the cavity page so one legend means one thing
   *  across the tool. */
  const FIELDS = [
    { id: 'none', label: 'Shape' },
    { id: 'electrostatic', label: 'Charge' },
    { id: 'hydrophobic', label: 'Greasiness' },
    { id: 'enclosure', label: 'Enclosure' },
    { id: 'from_anchor', label: 'Which protein lines it' }
  ];

  function segmented(options, current, onPick) {
    const seg = el('div', { class: 'seg' });
    options.forEach(function (o) {
      const b = el('button', {
        class: 'seg-btn' + (o.id === current ? ' on' : ''), text: o.label
      });
      b.addEventListener('click', function () {
        Array.from(seg.children).forEach(function (c) {
          c.classList.remove('on');
        });
        b.classList.add('on');
        onPick(o.id);
      });
      seg.appendChild(b);
    });
    return seg;
  }

  /** Draw the cavity points into their own canvas, sharing the viewer's
   *  rotation so the two pictures agree. */
  function drawCavity(canvas, sys, field) {
    const cav = sys.cavity;
    const ctx = canvas.getContext('2d');
    const w = canvas.clientWidth || 600;
    const h = canvas.clientHeight || 340;
    const dpr = window.devicePixelRatio || 1;
    canvas.width = w * dpr; canvas.height = h * dpr;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, w, h);
    if (!cav || !cav.points) return;

    const vals = field === 'none' ? null : cav[field];
    let lo = 0, hi = 1;
    if (vals) {
      lo = Math.min.apply(null, vals);
      hi = Math.max.apply(null, vals);
      if (hi - lo < 1e-6) hi = lo + 1;
    }

    // One frame of reference for the cavity and the ligand together: the
    // whole point of the picture is where the molecule sits in the pocket,
    // and scaling them separately would destroy exactly that.
    const lig = sys.ligand;
    let max = 1;
    for (let i = 0; i < cav.x.length; i++) {
      max = Math.max(max, Math.abs(cav.x[i]), Math.abs(cav.y[i]),
        Math.abs(cav.z[i]));
    }
    for (let i = 0; i < lig.x.length; i++) {
      max = Math.max(max, Math.abs(lig.x[i]), Math.abs(lig.y[i]),
        Math.abs(lig.z[i]));
    }
    const cx = cav.x.reduce(function (a, b) { return a + b; }, 0) / cav.x.length;
    const cy = cav.y.reduce(function (a, b) { return a + b; }, 0) / cav.y.length;
    const cz = cav.z.reduce(function (a, b) { return a + b; }, 0) / cav.z.length;
    const span = Math.max(12, max - Math.min(Math.abs(cx), Math.abs(cy)));
    const s = 0.42 * Math.min(w, h) / span;

    const ax = canvas._ax || -0.35, ay = canvas._ay || 0.6;
    function project(x, y, z) {
      const X0 = x - cx, Y0 = y - cy, Z0 = z - cz;
      const ca = Math.cos(ay), sa = Math.sin(ay);
      const X = X0 * ca + Z0 * sa;
      let Z = -X0 * sa + Z0 * ca;
      const cb = Math.cos(ax), sb = Math.sin(ax);
      const Y = Y0 * cb - Z * sb;
      Z = Y0 * sb + Z * cb;
      return [w / 2 + X * s, h / 2 - Y * s, Z];
    }

    const order = [];
    for (let i = 0; i < cav.x.length; i++) order.push(i);
    const proj = order.map(function (i) {
      return project(cav.x[i], cav.y[i], cav.z[i]);
    });
    order.sort(function (a, b) { return proj[a][2] - proj[b][2]; });

    for (let q = 0; q < order.length; q++) {
      const i = order[q];
      const p = proj[i];
      let col = 'rgba(150,150,160,0.5)';
      if (vals) {
        const t = (vals[i] - lo) / (hi - lo);
        if (field === 'electrostatic') {
          col = t < 0.5
            ? 'rgba(' + Math.round(70 + 300 * t) + ',110,210,0.72)'
            : 'rgba(220,' + Math.round(200 - 200 * (t - 0.5)) + ',90,0.72)';
        } else if (field === 'from_anchor') {
          col = 'rgba(' + Math.round(90 + 150 * t) + ',' +
            Math.round(190 - 60 * t) + ',' + Math.round(230 - 120 * t) + ',0.72)';
        } else {
          col = 'rgba(' + Math.round(90 + 150 * t) + ',' +
            Math.round(90 + 120 * t) + ',' + Math.round(220 - 120 * t) + ',0.7)';
        }
      }
      ctx.fillStyle = col;
      ctx.beginPath();
      ctx.arc(p[0], p[1], 2.3, 0, 6.2832);
      ctx.fill();
    }

    // The molecule, drawn over the pocket it has to fill.
    ctx.strokeStyle = 'rgba(245,240,230,0.95)';
    ctx.fillStyle = 'rgba(245,240,230,0.95)';
    ctx.lineWidth = 1.6;
    const lp = [];
    for (let i = 0; i < lig.x.length; i++) {
      lp.push(project(lig.x[i], lig.y[i], lig.z[i]));
    }
    for (let i = 0; i < lp.length; i++) {
      for (let j = i + 1; j < lp.length; j++) {
        const dx = lig.x[i] - lig.x[j], dy = lig.y[i] - lig.y[j],
          dz = lig.z[i] - lig.z[j];
        if (dx * dx + dy * dy + dz * dz > 3.0 * 3.0) continue;
        ctx.beginPath();
        ctx.moveTo(lp[i][0], lp[i][1]);
        ctx.lineTo(lp[j][0], lp[j][1]);
        ctx.stroke();
      }
    }
    for (let i = 0; i < lp.length; i++) {
      ctx.beginPath();
      ctx.arc(lp[i][0], lp[i][1], 2.2, 0, 6.2832);
      ctx.fill();
    }
  }

  function bindDrag(canvas, onMove) {
    let drag = null;
    canvas.addEventListener('pointerdown', function (e) {
      drag = { x: e.clientX, y: e.clientY,
               ax: canvas._ax || -0.35, ay: canvas._ay || 0.6 };
      canvas.setPointerCapture(e.pointerId);
    });
    canvas.addEventListener('pointermove', function (e) {
      if (!drag) return;
      canvas._ay = drag.ay + (e.clientX - drag.x) * 0.01;
      canvas._ax = drag.ax + (e.clientY - drag.y) * 0.01;
      onMove();
    });
    canvas.addEventListener('pointerup', function (e) {
      drag = null;
      canvas.releasePointerCapture(e.pointerId);
    });
  }

  /** Build a panel for one already-fetched blob.
   *
   *  Returns { node, viewer, setSystem } so a caller holding several poses
   *  can swap the structure without rebuilding the controls underneath it.
   */
  function panel(sys, opts) {
    opts = opts || {};
    const node = el('section', { class: 'panel viewpanel' });
    const head = el('div', { class: 'phead' }, [
      el('h3', { text: opts.title || 'Structure' })
    ]);
    node.appendChild(head);

    const sub = el('p', { class: 'sub' });
    node.appendChild(sub);

    const canvas = el('canvas', { class: 'viewer' });
    node.appendChild(canvas);
    const legend = el('div', { class: 'vlegend' });
    node.appendChild(legend);
    const controls = el('div', { class: 'controls' });
    node.appendChild(controls);

    let viewer = null;

    function drawLegend() {
      legend.innerHTML = '';
      viewer.legend().forEach(function (row) {
        const sw = el('span', { class: 'lsw' });
        sw.style.background = row.swatch;
        legend.appendChild(el('div', { class: 'lrow' }, [
          sw,
          el('span', { class: 'llabel', text: row.label }),
          el('span', { class: 'lnote', text: row.note })
        ]));
      });
    }

    function buildControls() {
      controls.innerHTML = '';
      controls.appendChild(el('div', { class: 'cgroup' }, [
        el('div', { class: 'clab' }, [el('span', { text: 'Style' })]),
        segmented(viewer.styles(), viewer.style, function (id) {
          viewer.style = id;
          viewer.draw();
          drawLegend();
        })
      ]));
      // Said once, plainly, rather than left for someone to notice that two
      // buttons they saw on the Validation page are not here.
      if (sys.drawing_only) {
        controls.appendChild(el('p', { class: 'fnote', text:
          'Contact energies and exit vectors are not computed for a docking '
          + 'result, so the styles that draw them are not offered. The '
          + 'Validation page has them for the benchmark systems.' }));
      }
    }

    function setSystem(next) {
      sys = next;
      viewer = new GF.Viewer(canvas, sys, function () {});
      const a = sys.assignment || {};
      sub.textContent = (sys.label ? sys.label + ' · ' : '') +
        (a.anchor_name || 'anchor') + ' and ' + (a.partner_name || 'partner') +
        ', with ' + (a.ligand || 'the ligand') + ' in place.' +
        ' Drag to turn, scroll to zoom.';
      buildControls();
      drawLegend();
      viewer.draw();
    }

    setSystem(sys);
    return { node: node, setSystem: setSystem,
             viewer: function () { return viewer; } };
  }

  /** The cavity at the interface, with the molecule sitting in it. */
  function cavityPanel(sys) {
    const node = el('section', { class: 'panel viewpanel' });
    node.appendChild(el('div', { class: 'phead' }, [
      el('h3', { text: 'The cavity, with the molecule in it' })
    ]));

    const cav = sys.cavity || {};
    if (!cav.points) {
      node.appendChild(el('p', { class: 'sub', text:
        'No cavity was found at this interface. That is a result, not a '
        + 'failure to draw one: a pose with no pocket beside the ligand has '
        + 'nowhere for added chemistry to go.' }));
      return { node: node, setSystem: function () {} };
    }

    const sub = el('p', { class: 'sub' });
    node.appendChild(sub);
    const canvas = el('canvas', { class: 'viewer' });
    node.appendChild(canvas);
    const controls = el('div', { class: 'controls' });
    node.appendChild(controls);

    let field = 'none';
    function redraw() { drawCavity(canvas, sys, field); }

    controls.appendChild(el('div', { class: 'cgroup' }, [
      el('div', { class: 'clab' }, [el('span', { text: 'Colour by' })]),
      segmented(FIELDS, field, function (id) { field = id; redraw(); })
    ]));

    function setSystem(next) {
      sys = next;
      const c = sys.cavity || {};
      const s = c.summary || {};
      sub.textContent = c.points
        ? Math.round(s.volume || 0) + ' Å³ over ' + c.points + ' points'
          + (c.anchor_lined_share !== undefined
            ? ', ' + Math.round(c.anchor_lined_share * 100) + '% of it lined '
              + 'by the anchor' : '')
          + '. Drag to turn.'
        : 'No cavity at this interface.';
      redraw();
    }

    bindDrag(canvas, redraw);
    setSystem(sys);
    return { node: node, setSystem: setSystem };
  }

  /** Fetch a run's view and hand back both panels.
   *
   *  `pose` names a docked placement; without it the run is drawn as it was
   *  uploaded, which is what an assessed complex wants.
   */
  function forRun(runId, pose) {
    return api.call('/api/run/' + encodeURIComponent(runId) + '/view' +
      (pose ? '?pose=' + encodeURIComponent(pose) : ''));
  }

  GF.view = {
    panel: panel,
    cavityPanel: cavityPanel,
    forRun: forRun,
    drawCavity: drawCavity
  };
})(window.GF);
