/* GlueForge — structure viewer.
 *
 * Four drawing styles, because they answer different questions. Contacts is
 * the one the science needs: it shows which pairs are frustrated and where
 * they sit relative to the ligand. The ribbon styles are for orienting
 * yourself in the fold before reading that.
 *
 *   contacts   backbone trace + contact chords, interface residues as nodes
 *   ribbon     smoothed CA spline, width by secondary structure
 *   sidechains ribbon + sidechain sticks on interface residues
 *   vectors    contacts + the ligand's exit vectors and what they reach
 *
 * Four things carry meaning, kept on separate visual channels so they cannot
 * be confused: chain identity is hue of the protein, interface is whether a
 * residue is drawn as a node, frustration is hue of the contact chords, and
 * the ligand is magenta, a colour used for nothing else.
 *
 * Secondary structure is estimated from CA geometry, not DSSP. It is for
 * drawing only and the page says so.
 */

(function (GF) {
  'use strict';

  const CHAIN_A = [47, 125, 140];    // teal, the anchor
  const CHAIN_B = [168, 118, 47];    // ochre, the partner
  const LIGAND = 'rgba(204,68,136,0.95)';
  const VECTOR = 'rgba(120,200,160,0.95)';

  const STYLES = [
    { id: 'contacts', label: 'Contacts' },
    { id: 'ribbon', label: 'Ribbon' },
    { id: 'sidechains', label: 'Ribbon + side chains' },
    { id: 'vectors', label: 'Exit vectors' }
  ];

  function rgba(c, a) {
    return 'rgba(' + c[0] + ',' + c[1] + ',' + c[2] + ',' + a + ')';
  }

  function rotate(x, y, z, ax, ay) {
    const ca = Math.cos(ay), sa = Math.sin(ay);
    const X = x * ca + z * sa;
    let Z = -x * sa + z * ca;
    const cb = Math.cos(ax), sb = Math.sin(ax);
    const Y = y * cb - Z * sb;
    Z = y * sb + Z * cb;
    return [X, Y, Z];
  }

  function Viewer(canvas, system, onHover) {
    this.canvas = canvas;
    this.ctx = canvas.getContext('2d');
    this.sys = system;
    this.onHover = onHover || function () {};
    this.ax = -0.35;
    this.ay = 0.6;
    this.scale = 1;
    // Contacts is the best view when there is contact data. A blob built
    // from coordinates alone has none, and defaulting to a style that
    // cannot draw would open the viewer on an empty canvas.
    this.style = this.has('contacts') ? 'contacts' : 'ribbon';
    this.filter = 'interface';     // interface | proximal | all
    this.live = null;
    this.highlight = -1;           // a contact to pick out, driven from the plot
    this.activeVector = -1;
    this._drag = null;
    this._proj = null;
    this._hover = -1;
    this._interfaceResidue = this._markInterfaceResidues();
    this._bind();
    this.fit();
  }

  Viewer.STYLES = STYLES;

  /** Does this system carry the data a style needs?
   *
   *  Two kinds of system reach this viewer. A benchmark entry has the whole
   *  frustration analysis behind it -- per-contact energies, exit vectors --
   *  and every style draws. A docking result is coordinates and nothing
   *  else, because the analysis takes eight hundred decoys and a published
   *  potential and is not run while someone waits.
   *
   *  So the styles that read those are not offered for it. The alternative
   *  would be to ship zeros and draw them, which is worse than an absent
   *  button: a flat grey contact map looks like a measurement that came out
   *  uniform rather than one that was never made.
   */
  Viewer.prototype.has = function (what) {
    const s = this.sys;
    if (what === 'contacts') {
      return !!(s.contacts && s.contacts.i && s.contacts.i.length);
    }
    if (what === 'vectors') {
      return !!(s.exit_vectors && s.exit_vectors.vectors &&
                s.exit_vectors.vectors.length);
    }
    return true;
  };

  /** The styles this system can actually be drawn in. */
  Viewer.prototype.styles = function () {
    const self = this;
    return STYLES.filter(function (st) {
      if (st.id === 'contacts') return self.has('contacts');
      if (st.id === 'vectors') return self.has('vectors');
      return true;
    });
  };

  Viewer.prototype._markInterfaceResidues = function () {
    const ct = this.sys.contacts;
    const n = this.sys.residues.x.length;
    const flag = new Uint8Array(n);
    // Without contacts the interface cannot be derived, so it is carried in
    // the blob instead -- computed from the coordinates by whoever built it.
    if (!this.has('contacts')) {
      const given = this.sys.residues.iface;
      if (given) for (let i = 0; i < n; i++) flag[i] = given[i];
      return flag;
    }
    for (let k = 0; k < ct.i.length; k++) {
      if (!ct.interface[k]) continue;
      const v = ct.proximal[k] ? 2 : 1;
      if (v > flag[ct.i[k]]) flag[ct.i[k]] = v;
      if (v > flag[ct.j[k]]) flag[ct.j[k]] = v;
    }
    return flag;
  };

  Viewer.prototype._bind = function () {
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

  Viewer.prototype.fit = function () {
    const r = this.sys.residues;
    let max = 1;
    for (let i = 0; i < r.x.length; i++) {
      const d = Math.sqrt(r.x[i] * r.x[i] + r.y[i] * r.y[i] + r.z[i] * r.z[i]);
      if (d > max) max = d;
    }
    this.radius = max;
  };

  Viewer.prototype.frustrationAt = function (k) {
    return this.live ? this.live[k] : this.sys.contacts.f[k];
  };

  Viewer.prototype.showsContacts = function () {
    return (this.style === 'contacts' || this.style === 'vectors') &&
      this.has('contacts');
  };

  Viewer.prototype.visible = function (k) {
    const ct = this.sys.contacts;
    if (this.filter === 'interface') return !!ct.interface[k];
    if (this.filter === 'proximal') return !!(ct.interface[k] && ct.proximal[k]);
    return true;
  };

  Viewer.prototype._pick = function (mx, my) {
    if (!this._proj || !this.showsContacts()) return;
    const ct = this.sys.contacts, p = this._proj;
    let best = -1, bestD = 9;
    for (let k = 0; k < ct.i.length; k++) {
      if (!this.visible(k)) continue;
      const a = ct.i[k], b = ct.j[k];
      const d = pointToSegment(mx, my, p.x[a], p.y[a], p.x[b], p.y[b]);
      if (d < bestD) { bestD = d; best = k; }
    }
    if (best !== this._hover) {
      this._hover = best;
      this.draw();
      this.onHover(best === -1 ? null : this.describe(best));
    }
  };

  Viewer.prototype.describe = function (k) {
    const ct = this.sys.contacts, r = this.sys.residues;
    const a = ct.i[k], b = ct.j[k];
    const f = this.frustrationAt(k);
    return {
      index: k,
      a: r.one[a] + r.seq[a] + ' ' + r.chain[a],
      b: r.one[b] + r.seq[b] + ' ' + r.chain[b],
      distance: ct.d[k],
      frustration: f,
      energy: ct.e[k],
      iface: !!ct.interface[k],
      proximal: !!ct.proximal[k],
      ligandDistance: Math.min(r.dlig[a], r.dlig[b]),
      verdict: f > 0.78 ? 'minimally frustrated'
             : f < -1 ? 'highly frustrated' : 'neutral'
    };
  };

  function pointToSegment(px, py, x1, y1, x2, y2) {
    const dx = x2 - x1, dy = y2 - y1;
    const len2 = dx * dx + dy * dy;
    let t = len2 ? ((px - x1) * dx + (py - y1) * dy) / len2 : 0;
    t = Math.max(0, Math.min(1, t));
    return Math.hypot(px - (x1 + t * dx), py - (y1 + t * dy));
  }

  /* ---- spline ---------------------------------------------------------- */

  /** Catmull-Rom through the CA positions, so the trace reads as a fold
   *  rather than a polyline of straight segments. */
  function spline(points, segments) {
    if (points.length < 2) return points.slice();
    const out = [];
    for (let i = 0; i < points.length - 1; i++) {
      const p0 = points[i === 0 ? 0 : i - 1];
      const p1 = points[i];
      const p2 = points[i + 1];
      const p3 = points[i + 2 >= points.length ? points.length - 1 : i + 2];
      for (let s = 0; s < segments; s++) {
        const t = s / segments, t2 = t * t, t3 = t2 * t;
        out.push([
          0.5 * ((2 * p1[0]) + (-p0[0] + p2[0]) * t +
            (2 * p0[0] - 5 * p1[0] + 4 * p2[0] - p3[0]) * t2 +
            (-p0[0] + 3 * p1[0] - 3 * p2[0] + p3[0]) * t3),
          0.5 * ((2 * p1[1]) + (-p0[1] + p2[1]) * t +
            (2 * p0[1] - 5 * p1[1] + 4 * p2[1] - p3[1]) * t2 +
            (-p0[1] + 3 * p1[1] - 3 * p2[1] + p3[1]) * t3),
          p1[2]
        ]);
      }
    }
    out.push(points[points.length - 1]);
    return out;
  }

  /* ---- drawing --------------------------------------------------------- */

  Viewer.prototype.draw = function () {
    const c = this.ctx, cv = this.canvas;
    const dpr = (typeof window !== 'undefined' && window.devicePixelRatio) || 1;
    const w = cv.clientWidth, h = cv.clientHeight;
    if (cv.width !== w * dpr || cv.height !== h * dpr) {
      cv.width = w * dpr; cv.height = h * dpr;
    }
    c.setTransform(dpr, 0, 0, dpr, 0, 0);
    c.clearRect(0, 0, w, h);

    const sys = this.sys, res = sys.residues;
    const n = res.x.length;
    const s = (Math.min(w, h) * 0.42 / this.radius) * this.scale;
    const cx = w / 2, cy = h / 2;
    const self = this;

    function project(x, y, z) {
      const p = rotate(x, y, z, self.ax, self.ay);
      return [cx + p[0] * s, cy - p[1] * s, p[2]];
    }

    const px = new Float32Array(n), py = new Float32Array(n), pz = new Float32Array(n);
    for (let i = 0; i < n; i++) {
      const p = project(res.x[i], res.y[i], res.z[i]);
      px[i] = p[0]; py[i] = p[1]; pz[i] = p[2];
    }
    this._proj = { x: px, y: py, z: pz };

    const anchor = sys.assignment.anchor_chain;
    const ribbon = this.style === 'ribbon' || this.style === 'sidechains';

    if (ribbon) this._drawRibbon(c, project, anchor);
    else this._drawTrace(c, px, py, anchor);

    if (this.showsContacts()) this._drawContacts(c, px, py, pz);
    if (this.style === 'sidechains') this._drawSidechains(c, project, px, py, anchor);
    if (!ribbon) this._drawNodes(c, px, py, anchor);

    this._drawLigand(c, project);
    if (this.style === 'vectors') this._drawVectors(c, project, px, py);
    this._drawScaleBar(c, s, h);
  };

  Viewer.prototype._drawTrace = function (c, px, py, anchor) {
    const res = this.sys.residues, n = res.x.length;
    const focus = this.filter !== 'all';
    c.lineWidth = 1.1;
    for (let i = 1; i < n; i++) {
      if (res.chain[i] !== res.chain[i - 1]) continue;
      if (res.seq[i] - res.seq[i - 1] > 4) continue;
      const near = this._interfaceResidue[i] || this._interfaceResidue[i - 1];
      c.strokeStyle = rgba(
        res.chain[i] === anchor ? CHAIN_A : CHAIN_B,
        focus ? (near ? 0.70 : 0.20) : 0.45
      );
      c.beginPath();
      c.moveTo(px[i - 1], py[i - 1]);
      c.lineTo(px[i], py[i]);
      c.stroke();
    }
  };

  /** Smoothed CA spline, stroked wider through helices and strands. Drawn in
   *  depth order so the fold occludes itself the way a solid object would. */
  Viewer.prototype._drawRibbon = function (c, project, anchor) {
    const res = this.sys.residues;
    const n = res.x.length;
    if (!res.cax) { this._drawTrace(c, this._proj.x, this._proj.y, anchor); return; }

    const runs = [];
    let start = 0;
    for (let i = 1; i <= n; i++) {
      const broken = i === n ||
        res.chain[i] !== res.chain[i - 1] ||
        res.seq[i] - res.seq[i - 1] > 4;
      if (broken) {
        if (i - start >= 2) runs.push([start, i]);
        start = i;
      }
    }

    const segs = [];
    runs.forEach(function (run) {
      const pts = [];
      for (let i = run[0]; i < run[1]; i++) {
        pts.push(project(res.cax[i], res.cay[i], res.caz[i]));
      }
      const sm = spline(pts, 6);
      const perStep = 6;
      for (let k = 1; k < sm.length; k++) {
        const srcIdx = run[0] + Math.min(
          Math.floor((k - 1) / perStep), run[1] - run[0] - 1
        );
        segs.push({
          a: sm[k - 1], b: sm[k],
          chain: res.chain[srcIdx],
          ss: res.ss ? res.ss[srcIdx] : 'C',
          iface: !!this._interfaceResidue[srcIdx],
          z: 0.5 * (sm[k - 1][2] + sm[k][2])
        });
      }
    }, this);

    segs.sort(function (a, b) { return a.z - b.z; });
    c.lineCap = 'round';
    c.lineJoin = 'round';
    for (let k = 0; k < segs.length; k++) {
      const sg = segs[k];
      const width = sg.ss === 'H' ? 5.0 : sg.ss === 'E' ? 3.8 : 1.9;
      const alpha = this.filter !== 'all' ? (sg.iface ? 0.95 : 0.45) : 0.8;
      c.strokeStyle = rgba(sg.chain === anchor ? CHAIN_A : CHAIN_B, alpha);
      c.lineWidth = width;
      c.beginPath();
      c.moveTo(sg.a[0], sg.a[1]);
      c.lineTo(sg.b[0], sg.b[1]);
      c.stroke();
    }
    c.lineCap = 'butt';
  };

  Viewer.prototype._drawSidechains = function (c, project, px, py, anchor) {
    const sc = this.sys.sidechains;
    if (!sc || !sc.res || !sc.res.length) return;
    const res = this.sys.residues;
    c.lineWidth = 1.2;
    for (let k = 0; k < sc.res.length; k++) {
      const ri = sc.res[k];
      const p = project(sc.x[k], sc.y[k], sc.z[k]);
      const prox = this._interfaceResidue[ri] === 2;
      c.strokeStyle = rgba(
        res.chain[ri] === anchor ? CHAIN_A : CHAIN_B, prox ? 0.9 : 0.45
      );
      c.beginPath();
      c.moveTo(px[ri], py[ri]);
      c.lineTo(p[0], p[1]);
      c.stroke();
      c.fillStyle = rgba(
        res.chain[ri] === anchor ? CHAIN_A : CHAIN_B, prox ? 0.9 : 0.45
      );
      c.beginPath();
      c.arc(p[0], p[1], 1.5, 0, Math.PI * 2);
      c.fill();
    }
  };

  Viewer.prototype._drawContacts = function (c, px, py, pz) {
    if (!this.has('contacts')) return;
    const ct = this.sys.contacts;
    const order = [];
    for (let k = 0; k < ct.i.length; k++) if (this.visible(k)) order.push(k);
    order.sort(function (a, b) {
      return (pz[ct.i[a]] + pz[ct.j[a]]) - (pz[ct.i[b]] + pz[ct.j[b]]);
    });

    const active = (this.activeVector >= 0 && this.has('vectors'))
      ? new Set(this.sys.exit_vectors.vectors[this.activeVector].reachable)
      : null;

    for (let q = 0; q < order.length; q++) {
      const k = order[q];
      const a = ct.i[k], b = ct.j[k];
      const emph = ct.interface[k] && ct.proximal[k];
      const isHover = k === this._hover || k === this.highlight;
      const dimmed = active && !active.has(k);
      c.strokeStyle = isHover
        ? 'rgba(255,255,255,0.95)'
        : GF.frustrationColour(
            this.frustrationAt(k),
            dimmed ? 0.08 : (emph ? 0.95 : 0.35)
          );
      c.lineWidth = isHover ? 3 : (emph ? 2 : 1);
      c.beginPath();
      c.moveTo(px[a], py[a]);
      c.lineTo(px[b], py[b]);
      c.stroke();
    }
  };

  Viewer.prototype._drawNodes = function (c, px, py, anchor) {
    const res = this.sys.residues, n = res.x.length;
    for (let i = 0; i < n; i++) {
      const flag = this._interfaceResidue[i];
      if (!flag) continue;
      const col = res.chain[i] === anchor ? CHAIN_A : CHAIN_B;
      c.fillStyle = rgba(col, flag === 2 ? 0.95 : 0.5);
      c.beginPath();
      c.arc(px[i], py[i], flag === 2 ? 3.4 : 2.2, 0, Math.PI * 2);
      c.fill();
    }
  };

  Viewer.prototype._drawLigand = function (c, project) {
    const lig = this.sys.ligand;
    if (!lig || !lig.x.length) return;
    for (let i = 0; i < lig.x.length; i++) {
      const p = project(lig.x[i], lig.y[i], lig.z[i]);
      c.fillStyle = LIGAND;
      c.beginPath();
      c.arc(p[0], p[1], 3.1, 0, Math.PI * 2);
      c.fill();
    }
  };

  /** Exit vectors as arrows off the ligand, scaled by how much frustrated
   *  contact each one could reach. */
  Viewer.prototype._drawVectors = function (c, project) {
    const ev = this.sys.exit_vectors;
    if (!ev || !ev.vectors.length) return;
    const best = Math.max.apply(null, ev.vectors.map(function (v) {
      return v.addressable;
    })) || 1;

    ev.vectors.forEach(function (v, idx) {
      const len = 3.0 + 7.0 * (v.addressable / best);
      const from = project(v.position[0], v.position[1], v.position[2]);
      const to = project(
        v.position[0] + v.direction[0] * len,
        v.position[1] + v.direction[1] * len,
        v.position[2] + v.direction[2] * len
      );
      const active = this.activeVector === idx;
      const strength = v.addressable / best;
      c.strokeStyle = active
        ? 'rgba(255,255,255,0.95)'
        : 'rgba(120,200,160,' + (0.25 + 0.7 * strength).toFixed(2) + ')';
      c.lineWidth = active ? 3 : 1.6;
      c.beginPath();
      c.moveTo(from[0], from[1]);
      c.lineTo(to[0], to[1]);
      c.stroke();
      c.fillStyle = c.strokeStyle;
      c.beginPath();
      c.arc(to[0], to[1], active ? 3.6 : 2.4, 0, Math.PI * 2);
      c.fill();
    }, this);
  };

  Viewer.prototype._drawScaleBar = function (c, s, h) {
    const bar = 20 * s;
    c.strokeStyle = 'rgba(140,136,130,0.7)';
    c.lineWidth = 1;
    c.beginPath();
    c.moveTo(14, h - 16);
    c.lineTo(14 + bar, h - 16);
    c.stroke();
    c.fillStyle = 'rgba(140,136,130,0.9)';
    c.font = '10px ui-monospace, SFMono-Regular, Menlo, monospace';
    c.fillText('20 Å', 14, h - 22);
  };

  Viewer.prototype.legend = function () {
    const a = this.sys.assignment;
    const rows = [
      { swatch: rgba(CHAIN_A, 0.95), label: a.anchor_name, note: 'anchor · chain ' + a.anchor_chain },
      { swatch: rgba(CHAIN_B, 0.95), label: a.partner_name, note: 'partner · chain ' + a.partner_chain },
      { swatch: LIGAND, label: a.ligand, note: 'ligand' }
    ];
    if (this.showsContacts()) {
      rows.push({
        swatch: GF.frustrationColour(-2, 0.95),
        label: 'frustrated contact', note: 'F < −1'
      });
      rows.push({
        swatch: GF.frustrationColour(2, 0.95),
        label: 'minimally frustrated', note: 'F > +0.78'
      });
    }
    if (this.style === 'vectors') {
      rows.push({ swatch: VECTOR, label: 'exit vector', note: 'length = reach' });
    }
    if (this.style === 'ribbon' || this.style === 'sidechains') {
      rows.push({
        swatch: 'transparent', label: 'wide = helix or strand',
        note: 'estimated from CA'
      });
    }
    return rows;
  };

  GF.Viewer = Viewer;
  GF.Viewer.CHAIN_A = CHAIN_A;
  GF.Viewer.CHAIN_B = CHAIN_B;
})(globalThis.GF || (globalThis.GF = {}));
