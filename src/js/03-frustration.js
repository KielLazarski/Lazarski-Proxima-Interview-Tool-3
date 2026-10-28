/* GlueForge — the frustration engine, browser side.
 *
 * This is the one expensive step that recomputes live rather than being read
 * from the cache, which is what lets the page strip the glue, move a residue
 * or swap a model and show the frustration map change.
 *
 * It mirrors sidecar/glueforge/frustration.py line for line. The two are
 * checked against each other by tests/crossval.test.js, which only works
 * because both draw decoys from the same mulberry32 stream in the same order:
 * for each contact, for each decoy, the identity pair first and then the
 * resampled environment.
 *
 * Sign convention, same as the reference:
 *   F = (mean(E_decoy) - E_native) / std(E_decoy)
 * Large positive F is minimally frustrated; negative F is frustrated.
 */

(function (GF) {
  'use strict';

  const MINIMALLY_FRUSTRATED = 0.78;
  const HIGHLY_FRUSTRATED = -1.0;

  /* ---- deterministic PRNG -------------------------------------------- */

  /** mulberry32. Its state advances by a fixed increment, so the k-th draw
   *  is a pure function of the seed and k — which is what lets the Python
   *  side produce the same stream vectorised. */
  function mulberry32(seed) {
    let a = seed | 0;
    return function () {
      a = (a + 0x6d2b79f5) | 0;
      let t = a;
      t = Math.imul(t ^ (t >>> 15), 1 | t);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  /* ---- switching functions ------------------------------------------- */

  function theta(r, rMin, rMax, eta) {
    return (
      0.25 *
      (1 + Math.tanh(eta * (r - rMin))) *
      (1 + Math.tanh(eta * (rMax - r)))
    );
  }

  /** AWSEM burial term: how well a residue type suits its local density.
   *  Three density shells, each a tanh window. Returns 0 when the potential
   *  carries no burial table, which is the case for the stand-in. */
  function burialEnergy(pot, aa, rho) {
    const gb = pot.burial;
    if (!gb) return 0;
    const bp = pot.burial_params;
    const kappa = bp.kappa, lo = bp.ro_min, hi = bp.ro_max;
    let acc = 0;
    for (let k = 0; k < 3; k++) {
      acc += gb[aa][k] * (
        Math.tanh(kappa * (rho - lo[k])) + Math.tanh(kappa * (hi[k] - rho))
      );
    }
    return -0.5 * pot.params.k_contact * acc;
  }

  function sigmaWater(rhoI, rhoJ, etaSigma, rho0) {
    return (
      0.25 *
      (1 - Math.tanh(etaSigma * (rhoI - rho0))) *
      (1 - Math.tanh(etaSigma * (rhoJ - rho0)))
    );
  }

  /* ---- geometry ------------------------------------------------------- */

  /** Pairwise CB-CB distances as a flat Float64Array of n*n. */
  function distanceMatrix(x, y, z) {
    const n = x.length;
    const d = new Float64Array(n * n);
    for (let i = 0; i < n; i++) {
      const xi = x[i], yi = y[i], zi = z[i];
      for (let j = i + 1; j < n; j++) {
        const dx = xi - x[j], dy = yi - y[j], dz = zi - z[j];
        const r = Math.sqrt(dx * dx + dy * dy + dz * dz);
        d[i * n + j] = r;
        d[j * n + i] = r;
      }
    }
    return d;
  }

  /** rho_i, the AWSEM local density over the direct-well window. Pairs are
   *  skipped when they are close in sequence *within the same chain*; an
   *  inter-chain pair always counts, since sequence separation is undefined
   *  across chains. */
  function localDensity(dist, chainIds, params) {
    const n = chainIds.length;
    const rho = new Float64Array(n);
    const minSep = params.min_sequence_separation_rho;
    for (let i = 0; i < n; i++) {
      let acc = 0;
      for (let j = 0; j < n; j++) {
        if (i === j) continue;
        if (chainIds[i] === chainIds[j] && Math.abs(i - j) < minSep) continue;
        acc += theta(dist[i * n + j], params.r_min, params.r_max, params.eta);
      }
      rho[i] = acc;
    }
    return rho;
  }

  /* ---- the engine ----------------------------------------------------- */

  const DEFAULTS = {
    nDecoys: 800,
    seed: 20261026,
    minSepContact: 2,
    contactCutoff: 9.5,
    flavour: 'configurational'
  };

  /**
   * @param {{x:Float64Array|number[], y:*, z:*, aa:Int32Array|number[],
   *          chainIds:Int32Array|number[]}} site
   *        aa holds indices into the potential's 20-letter alphabet.
   * @param {object} pot  potential, as written by the sidecar.
   * @param {object} [opts]
   */
  function compute(site, pot, opts) {
    const cfg = Object.assign({}, DEFAULTS, opts || {});
    const p = pot.params;
    const n = site.aa.length;
    if (n < 2) throw new Error('need at least two residues');

    const dist = distanceMatrix(site.x, site.y, site.z);
    const rho = localDensity(dist, site.chainIds, p);

    // Contact list, upper triangle, in the same order numpy's nonzero gives:
    // row-major, so i ascending and j ascending within i.
    const ii = [], jj = [], dd = [];
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        const sameChain = site.chainIds[i] === site.chainIds[j];
        if (sameChain && j - i < cfg.minSepContact) continue;
        const r = dist[i * n + j];
        if (r >= cfg.contactCutoff) continue;
        ii.push(i); jj.push(j); dd.push(r);
      }
    }
    const m = ii.length;
    if (m === 0) throw new Error('no contacts found; check the selection');

    const gd = pot.direct, gw = pot.water, gp = pot.protein;
    const k = p.k_contact;

    const energy = new Float64Array(m);
    const th = new Float64Array(m);
    const th2 = new Float64Array(m);
    const sw = new Float64Array(m);

    for (let c = 0; c < m; c++) {
      const r = dd[c];
      th[c] = theta(r, p.r_min, p.r_max, p.eta);
      th2[c] = theta(r, p.r_minII, p.r_maxII, p.eta);
      sw[c] = sigmaWater(rho[ii[c]], rho[jj[c]], p.eta_sigma, p.rho_0);
      const a = site.aa[ii[c]], b = site.aa[jj[c]];
      energy[c] = -k * (
        th[c] * gd[a][b] +
        sw[c] * th2[c] * gw[a][b] +
        (1 - sw[c]) * th2[c] * gp[a][b]
      ) + burialEnergy(pot, a, rho[ii[c]]) + burialEnergy(pot, b, rho[jj[c]]);
    }

    const configurational = cfg.flavour === 'configurational';
    const drawsPer = configurational ? 5 : 2;
    const rand = mulberry32(cfg.seed);
    const nd = cfg.nDecoys;

    const frustration = new Float64Array(m);
    const decoyMean = new Float64Array(m);
    const decoyStd = new Float64Array(m);
    const buf = new Float64Array(nd);

    for (let c = 0; c < m; c++) {
      let sum = 0;
      for (let dIdx = 0; dIdx < nd; dIdx++) {
        const ua = rand(), ub = rand();
        let a = Math.floor(ua * n); if (a > n - 1) a = n - 1;
        let b = Math.floor(ub * n); if (b > n - 1) b = n - 1;
        const ia = site.aa[a], ib = site.aa[b];

        let t1, t2, s, rhoA, rhoB;
        if (configurational) {
          const ur = rand(), ui = rand(), uj = rand();
          let ri = Math.floor(ur * m); if (ri > m - 1) ri = m - 1;
          let qi = Math.floor(ui * n); if (qi > n - 1) qi = n - 1;
          let qj = Math.floor(uj * n); if (qj > n - 1) qj = n - 1;
          const rr = dd[ri];
          t1 = theta(rr, p.r_min, p.r_max, p.eta);
          t2 = theta(rr, p.r_minII, p.r_maxII, p.eta);
          s = sigmaWater(rho[qi], rho[qj], p.eta_sigma, p.rho_0);
          // Burial follows the resampled environment too.
          rhoA = rho[qi]; rhoB = rho[qj];
        } else {
          t1 = th[c]; t2 = th2[c]; s = sw[c];
          rhoA = rho[ii[c]]; rhoB = rho[jj[c]];
        }

        const e = -k * (
          t1 * gd[ia][ib] +
          s * t2 * gw[ia][ib] +
          (1 - s) * t2 * gp[ia][ib]
        ) + burialEnergy(pot, ia, rhoA) + burialEnergy(pot, ib, rhoB);
        buf[dIdx] = e;
        sum += e;
      }
      const mu = sum / nd;
      let varAcc = 0;
      for (let dIdx = 0; dIdx < nd; dIdx++) {
        const dv = buf[dIdx] - mu;
        varAcc += dv * dv;
      }
      const sd = Math.sqrt(varAcc / nd); // population sd, matching numpy's default
      decoyMean[c] = mu;
      decoyStd[c] = sd;
      frustration[c] = (mu - energy[c]) / (sd < 1e-12 ? 1e-12 : sd);
    }

    return {
      i: Int32Array.from(ii),
      j: Int32Array.from(jj),
      distance: Float64Array.from(dd),
      energy,
      decoyMean,
      decoyStd,
      frustration,
      density: rho,
      flavour: cfg.flavour,
      potential: pot.name,
      config: cfg
    };
  }

  function classify(f) {
    if (f > MINIMALLY_FRUSTRATED) return 1;
    if (f < HIGHLY_FRUSTRATED) return -1;
    return 0;
  }

  GF.frustration = {
    compute,
    classify,
    mulberry32,
    theta,
    sigmaWater,
    burialEnergy,
    distanceMatrix,
    localDensity,
    MINIMALLY_FRUSTRATED,
    HIGHLY_FRUSTRATED,
    DEFAULTS
  };
})(globalThis.GF || (globalThis.GF = {}));
