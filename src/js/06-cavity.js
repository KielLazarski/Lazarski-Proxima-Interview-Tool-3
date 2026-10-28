/* GlueForge — the cavity page.
 *
 * The composite pocket at a glued interface, its pharmacophore fields, and
 * how much partner-placement error it survives. This is the design target:
 * unlike frustration, a molecule placed here is scored directly rather than
 * through the geometry it happens to induce.
 *
 * The page leads with the error budget because that is the number which
 * decides whether a design campaign is viable, and with the anchor/partner
 * split because that is why the problem is better conditioned than it looks.
 */

(function (GF) {
  'use strict';

  const el = GF.el, fmt = GF.fmt, help = GF.help;

  const FIELDS = [
    { id: 'electrostatic', label: 'Electrostatic', unit: 'kcal/mol/e',
      note: 'negative attracts a cation' },
    { id: 'donor', label: 'H-bond donors', unit: '', note: 'protein donors nearby' },
    { id: 'acceptor', label: 'H-bond acceptors', unit: '', note: 'protein acceptors nearby' },
    { id: 'hydrophobic', label: 'Hydrophobic', unit: '', note: 'apolar contact' },
    { id: 'from_anchor', label: 'Lining', unit: '', note: '1 = anchor, 0 = partner' }
  ];

  function cavityPage(entry) {
    const sys = GF.data.system(entry);
    const root = el('div', { class: 'page' });
    if (!sys || !sys.cavity || !sys.cavity.points) {
      root.appendChild(el('h1', { text: 'Cavity · ' + entry }));
      root.appendChild(el('p', {
        class: 'lede',
        text: sys
          ? 'No cavity was detected at this site. That happens where the ' +
            'glue sits in a surface groove rather than an enclosed pocket.'
          : 'No cached data for ' + entry
      }));
      return root;
    }

    const cav = sys.cavity;
    const s = cav.summary || {};
    const a = sys.assignment;

    root.appendChild(el('h1', { text: 'Cavity · ' + entry + ' · ' + a.ligand }));
    root.appendChild(el('p', {
      class: 'lede',
      text:
        'The pocket a molecule has to fill: part of it is ' + a.anchor_name +
        '’s, part is completed by ' + a.partner_name + '. Unlike ' +
        'frustration, this objective has the molecule in the calculation.'
    }));

    /* --- headline ------------------------------------------------------ */
    const cards = el('div', { class: 'cards' });
    cards.appendChild(card(fmt.fixed(s.volume, 0) + ' Å³', 'volume',
      cav.total_points + ' grid points', 'cavity_volume'));
    cards.appendChild(card(fmt.fixed(s.enclosure_quality, 2), 'enclosure quality',
      s.enclosure_quality >= 0.4 ? 'an enclosed pocket' : 'closer to a surface groove',
      'enclosure_quality'));
    cards.appendChild(card(Math.round(cav.anchor_lined_share * 100) + '%',
      'anchor-lined', 'known exactly from the co-crystal', 'cavity_split'));
    cards.appendChild(card(fmt.signed(s.mean_electrostatic, 2), 'mean field',
      s.mean_electrostatic < 0 ? 'wants a cation or donor' : 'wants an anion',
      'cavity_field'));
    root.appendChild(cards);

    /* --- the error budget ---------------------------------------------- */
    if (cav.sensitivity && cav.sensitivity.rows) {
      const h = el('div', { class: 'phead' }, [
        el('h3', { text: 'How much docking error this cavity survives' })
      ]);
      help.attach(h, 'error_budget');
      root.appendChild(h);
      root.appendChild(el('p', {
        class: 'lede',
        text:
          'The partner is displaced rigidly; the anchor and the ligand stay ' +
          'put. Shape degrades quickly, the fields far less so, and the ' +
          'anchor-lined part of the cavity barely moves at all.'
      }));
      root.appendChild(budgetTable(cav.sensitivity.rows));
      root.appendChild(el('p', {
        class: 'hint',
        text:
          'Rigid-body, so no induced fit. A real partner would relax into a ' +
          'displaced position, which should make the true degradation gentler ' +
          'than this.'
      }));
    }

    /* --- fields --------------------------------------------------------- */
    const fh = el('div', { class: 'phead' }, [
      el('h3', { text: 'Pharmacophore fields' })
    ]);
    help.attach(fh, 'cavity_field');
    root.appendChild(fh);

    const split = el('div', { class: 'split' });
    const left = el('div', { class: 'panel' });
    const canvas = el('canvas', { class: 'viewer' });
    left.appendChild(canvas);
    const readout = el('div', { class: 'readout', text: 'Hover a cavity point.' });
    left.appendChild(readout);
    const controls = el('div', { class: 'controls' });
    left.appendChild(controls);
    split.appendChild(left);

    const right = el('div', { class: 'panel' });
    right.appendChild(fieldTable(cav, s));
    split.appendChild(right);
    root.appendChild(split);

    requestAnimationFrame(function () {
      const v = new CavityView(canvas, sys, function (hit) {
        readout.innerHTML = '';
        if (!hit) {
          readout.textContent = 'Hover a cavity point.';
          readout.classList.remove('on');
          return;
        }
        readout.classList.add('on');
        readout.appendChild(el('span', { class: 'mono rpair', text: hit.label }));
        readout.appendChild(el('span', { class: 'mono dim', text: hit.detail }));
      });
      v.draw();
      buildCavityControls(controls, v);
      window.addEventListener('resize', function () { v.draw(); });
    });

    return root;
  }

  function card(value, label, note, slug) {
    const c = el('div', { class: 'card' }, [el('div', { class: 'cval', text: value })]);
    const lab = el('div', { class: 'clabel' }, [el('span', { text: label })]);
    if (slug) help.attach(lab, slug);
    c.appendChild(lab);
    c.appendChild(el('div', { class: 'cnote', text: note }));
    return c;
  }

  function budgetTable(rows) {
    const t = el('table', { class: 'grid' });
    t.appendChild(el('thead', {}, [
      el('tr', {}, [
        el('th', { text: 'Partner error' }),
        el('th', { class: 'num', text: 'Shape (Jaccard)' }),
        el('th', { class: 'num', text: 'r(electrostatic)' }),
        el('th', { class: 'num', text: 'r(donor)' }),
        el('th', { class: 'num', text: 'Anchor side kept' }),
        el('th', { class: 'num', text: 'Partner side kept' })
      ])
    ]));
    const body = el('tbody');
    rows.forEach(function (r) {
      const tr = el('tr', { class: 'row' });
      tr.appendChild(el('td', { class: 'mono', text: r.shift_angstrom + ' Å' }));
      [r.jaccard, r.elec, r.donor, r.anchor_keep, r.partner_keep].forEach(function (v, i) {
        const td = el('td', { class: 'num' + (v === null ? ' dim' : '') });
        td.textContent = v === null ? '—' : fmt.fixed(v, 2);
        if (v !== null && i === 0) {
          td.style.color = v > 0.7 ? 'var(--ok)' : v > 0.5 ? 'var(--ink)' : 'var(--frustrated)';
          td.style.fontWeight = '600';
        }
        tr.appendChild(td);
      });
      body.appendChild(tr);
    });
    t.appendChild(body);
    return t;
  }

  function fieldTable(cav, s) {
    const t = el('table', { class: 'kv' });
    const rows = [
      ['Grid points', cav.total_points + ' (' + cav.points + ' drawn)'],
      ['Spacing', cav.spacing + ' Å'],
      ['Composite fraction', fmt.fixed(s.composite_fraction, 2)],
      ['Anchor-lined', Math.round(cav.anchor_lined_share * 100) + '%'],
      ['Partner-lined', Math.round(cav.partner_lined_share * 100) + '%'],
      ['Electrostatic range', fmt.signed(s.electrostatic_range ? s.electrostatic_range[0] : 0, 1) +
        ' to ' + fmt.signed(s.electrostatic_range ? s.electrostatic_range[1] : 0, 1)],
      ['Donor-rich points', Math.round((s.donor_rich_fraction || 0) * 100) + '%'],
      ['Acceptor-rich points', Math.round((s.acceptor_rich_fraction || 0) * 100) + '%'],
      ['Mean hydrophobic', fmt.fixed(s.mean_hydrophobic, 2)]
    ];
    rows.forEach(function (kv) {
      t.appendChild(el('tr', {}, [
        el('td', { class: 'k', text: kv[0] }),
        el('td', { class: 'v mono', text: kv[1] })
      ]));
    });
    return t;
  }

  function buildCavityControls(host, view) {
    const lab = el('div', { class: 'clab' }, [el('span', { text: 'Colour by' })]);
    help.attach(lab, 'cavity_field');
    const seg = el('div', { class: 'seg' });
    FIELDS.forEach(function (f, i) {
      const b = el('button', {
        class: 'seg-btn' + (i === 0 ? ' on' : ''), text: f.label
      });
      b.addEventListener('click', function () {
        view.field = f.id;
        view.draw();
        Array.from(seg.children).forEach(function (c) { c.classList.remove('on'); });
        b.classList.add('on');
      });
      seg.appendChild(b);
    });
    host.appendChild(el('div', { class: 'cgroup' }, [lab, seg]));
    host.appendChild(el('p', {
      class: 'hint',
      text:
        'Drag to rotate, scroll to zoom. Each dot is a cavity grid point ' +
        'coloured by the selected field; the magenta dots are the ligand, ' +
        'shown for reference. A generative model would be conditioned on ' +
        'these fields rather than on the shape alone — they are the part ' +
        'that survives docking error.'
    }));
  }

  /* --- the cavity viewer ------------------------------------------------ */

  function CavityView(canvas, sys, onHover) {
    this.canvas = canvas;
    this.ctx = canvas.getContext('2d');
    this.sys = sys;
    this.cav = sys.cavity;
    this.onHover = onHover || function () {};
    this.field = 'electrostatic';
    this.ax = -0.35;
    this.ay = 0.6;
    this.scale = 1;
    this._drag = null;
    this._proj = null;
    this._hover = -1;
    let r = 1;
    for (let i = 0; i < this.cav.x.length; i++) {
      r = Math.max(r, Math.hypot(this.cav.x[i], this.cav.y[i], this.cav.z[i]));
    }
    this.radius = r;
    this._bind();
  }

  CavityView.prototype._bind = function () {
    const self = this;
    this.canvas.addEventListener('pointerdown', function (e) {
      self._drag = { x: e.clientX, y: e.clientY, ax: self.ax, ay: self.ay };
      self.canvas.setPointerCapture(e.pointerId);
    });
    this.canvas.addEventListener('pointermove', function (e) {
      if (self._drag) {
        self.ay = self._drag.ay + (e.clientX - self._drag.x) * 0.01;
        self.ax = self._drag.ax + (e.clientY - self._drag.y) * 0.01;
        self.draw();
        return;
      }
      const r = self.canvas.getBoundingClientRect();
      self._pick(e.clientX - r.left, e.clientY - r.top);
    });
    this.canvas.addEventListener('pointerleave', function () {
      if (self._hover !== -1) { self._hover = -1; self.draw(); self.onHover(null); }
    });
    this.canvas.addEventListener('pointerup', function (e) {
      self._drag = null;
      self.canvas.releasePointerCapture(e.pointerId);
    });
    this.canvas.addEventListener('wheel', function (e) {
      e.preventDefault();
      self.scale *= e.deltaY > 0 ? 0.92 : 1.08;
      self.scale = Math.max(0.3, Math.min(6, self.scale));
      self.draw();
    }, { passive: false });
  };

  CavityView.prototype._pick = function (mx, my) {
    if (!this._proj) return;
    let best = -1, bestD = 7;
    for (let i = 0; i < this._proj.x.length; i++) {
      const d = Math.hypot(mx - this._proj.x[i], my - this._proj.y[i]);
      if (d < bestD) { bestD = d; best = i; }
    }
    if (best !== this._hover) {
      this._hover = best;
      this.draw();
      if (best === -1) { this.onHover(null); return; }
      const c = this.cav;
      this.onHover({
        label: 'cavity point ' + best,
        detail:
          'elec ' + fmt.signed(c.electrostatic[best], 2) +
          ' · donor ' + fmt.fixed(c.donor[best], 2) +
          ' · acceptor ' + fmt.fixed(c.acceptor[best], 2) +
          ' · phobic ' + fmt.fixed(c.hydrophobic[best], 2) +
          ' · ' + Math.round(c.from_anchor[best] * 100) + '% anchor-lined' +
          ' · enclosed on ' + c.enclosure[best] + '/3 axes'
      });
    }
  };

  /** Diverging for signed fields, sequential for the rest. */
  CavityView.prototype.colour = function (i, alpha) {
    const c = this.cav;
    const f = this.field;
    if (f === 'electrostatic') {
      const lo = c.summary.electrostatic_range[0];
      const hi = c.summary.electrostatic_range[1];
      const span = Math.max(Math.abs(lo), Math.abs(hi), 1e-6);
      return GF.frustrationColour(-(c.electrostatic[i] / span) * 2.5, alpha);
    }
    if (f === 'from_anchor') {
      const t = c.from_anchor[i];
      const A = GF.Viewer.CHAIN_A, B = GF.Viewer.CHAIN_B;
      return 'rgba(' +
        Math.round(B[0] + (A[0] - B[0]) * t) + ',' +
        Math.round(B[1] + (A[1] - B[1]) * t) + ',' +
        Math.round(B[2] + (A[2] - B[2]) * t) + ',' + alpha + ')';
    }
    const v = c[f];
    let max = 1e-6;
    for (let k = 0; k < v.length; k++) max = Math.max(max, v[k]);
    const t = Math.min(1, v[i] / max);
    return 'rgba(' + Math.round(70 + 150 * t) + ',' +
      Math.round(90 + 90 * t) + ',' + Math.round(160 - 40 * t) + ',' + alpha + ')';
  };

  CavityView.prototype.draw = function () {
    const c = this.ctx, cv = this.canvas;
    const dpr = (typeof window !== 'undefined' && window.devicePixelRatio) || 1;
    const w = cv.clientWidth, h = cv.clientHeight;
    if (cv.width !== w * dpr || cv.height !== h * dpr) {
      cv.width = w * dpr; cv.height = h * dpr;
    }
    c.setTransform(dpr, 0, 0, dpr, 0, 0);
    c.clearRect(0, 0, w, h);

    const cav = this.cav;
    const n = cav.x.length;
    const s = (Math.min(w, h) * 0.40 / this.radius) * this.scale;
    const cx = w / 2, cy = h / 2;
    const ca = Math.cos(this.ay), sa = Math.sin(this.ay);
    const cb = Math.cos(this.ax), sb = Math.sin(this.ax);
    const self = this;

    function proj(x, y, z) {
      const X = x * ca + z * sa;
      let Z = -x * sa + z * ca;
      const Y = y * cb - Z * sb;
      Z = y * sb + Z * cb;
      return [cx + X * s, cy - Y * s, Z];
    }

    const px = new Float32Array(n), py = new Float32Array(n), pz = new Float32Array(n);
    const order = [];
    for (let i = 0; i < n; i++) {
      const p = proj(cav.x[i], cav.y[i], cav.z[i]);
      px[i] = p[0]; py[i] = p[1]; pz[i] = p[2];
      order.push(i);
    }
    this._proj = { x: px, y: py };
    order.sort(function (a, b) { return pz[a] - pz[b]; });

    for (let q = 0; q < order.length; q++) {
      const i = order[q];
      c.fillStyle = i === this._hover
        ? 'rgba(255,255,255,0.95)'
        : self.colour(i, 0.55);
      c.beginPath();
      c.arc(px[i], py[i], i === this._hover ? 5 : 3.2, 0, Math.PI * 2);
      c.fill();
    }

    const lig = this.sys.ligand;
    if (lig && lig.x.length) {
      for (let i = 0; i < lig.x.length; i++) {
        const p = proj(lig.x[i], lig.y[i], lig.z[i]);
        c.fillStyle = 'rgba(204,68,136,0.95)';
        c.beginPath();
        c.arc(p[0], p[1], 2.6, 0, Math.PI * 2);
        c.fill();
      }
    }

    const bar = 10 * s;
    c.strokeStyle = 'rgba(140,136,130,0.7)';
    c.lineWidth = 1;
    c.beginPath();
    c.moveTo(14, h - 16); c.lineTo(14 + bar, h - 16); c.stroke();
    c.fillStyle = 'rgba(140,136,130,0.9)';
    c.font = '10px ui-monospace, SFMono-Regular, Menlo, monospace';
    c.fillText('10 Å', 14, h - 22);
  };

  GF.pages = GF.pages || {};
  GF.pages.cavity = cavityPage;
  GF.CavityView = CavityView;
})(globalThis.GF || (globalThis.GF = {}));
