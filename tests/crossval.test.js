/* Cross-engine check.
 *
 * The architecture claims the browser's frustration engine is validated
 * against the Python reference. This is that check, and it is exact rather
 * than statistical: both engines draw decoys from the same mulberry32 stream
 * in the same order, so disagreement means a real difference in the maths,
 * not sampling noise.
 *
 * Run: node tests/crossval.test.js
 */

'use strict';

const fs = require('fs');
const path = require('path');
const vm = require('vm');

const ROOT = path.resolve(__dirname, '..');
const FIXTURE = path.join(ROOT, 'cache', 'web', 'crossval.json');

// Tolerance: these are IEEE754 doubles through tanh and sqrt on both sides,
// so agreement should be at the last few bits, not merely "close".
const TOL_ENERGY = 1e-9;
const TOL_FRUSTRATION = 1e-9;
const TOL_DENSITY = 1e-9;

function loadEngine() {
  const sandbox = { globalThis: {}, Math, console };
  sandbox.globalThis = sandbox;
  const src = fs.readFileSync(
    path.join(ROOT, 'src', 'js', '03-frustration.js'), 'utf8'
  );
  vm.createContext(sandbox);
  vm.runInContext(src, sandbox, { filename: '03-frustration.js' });
  if (!sandbox.GF || !sandbox.GF.frustration) {
    throw new Error('engine did not register itself on GF');
  }
  return sandbox.GF.frustration;
}

function maxAbsDiff(a, b) {
  if (a.length !== b.length) {
    return { diff: Infinity, at: -1, note: `length ${a.length} vs ${b.length}` };
  }
  let worst = 0, at = -1;
  for (let i = 0; i < a.length; i++) {
    const d = Math.abs(a[i] - b[i]);
    if (d > worst) { worst = d; at = i; }
  }
  return { diff: worst, at };
}

function chainIdsFrom(chains) {
  const seen = new Map();
  return Int32Array.from(chains, (c) => {
    if (!seen.has(c)) seen.set(c, seen.size);
    return seen.get(c);
  });
}

let failures = 0;
function check(label, ok, detail) {
  const mark = ok ? 'PASS' : 'FAIL';
  if (!ok) failures++;
  console.log(`  [${mark}] ${label}${detail ? '  ' + detail : ''}`);
}

function main() {
  if (!fs.existsSync(FIXTURE)) {
    console.error(
      'crossval.json missing. Build it first:\n' +
      '  python3 sidecar/export_web.py'
    );
    process.exit(2);
  }

  const fx = JSON.parse(fs.readFileSync(FIXTURE, 'utf8'));
  const engine = loadEngine();

  console.log(
    `cross-engine check: ${fx.source}, ${fx.input.one.length} residues, ` +
    `potential ${fx.potential.name}` +
    (fx.potential.published ? '' : ' (STAND-IN)')
  );

  const alphabet = fx.potential.alphabet;
  const aaIndex = new Map(alphabet.map((a, k) => [a, k]));
  const site = {
    x: Float64Array.from(fx.input.x),
    y: Float64Array.from(fx.input.y),
    z: Float64Array.from(fx.input.z),
    aa: Int32Array.from(fx.input.one, (c) => aaIndex.get(c) ?? 0),
    chainIds: chainIdsFrom(fx.input.chain)
  };

  // --- PRNG first: everything downstream depends on it -------------------
  const rand = engine.mulberry32(fx.config.seed);
  const first = [];
  for (let k = 0; k < 8; k++) first.push(rand());
  const inRange = first.every((v) => v >= 0 && v < 1);
  check('PRNG draws land in [0, 1)', inRange, `first=${first[0].toFixed(12)}`);

  // --- contact list ------------------------------------------------------
  const conf = engine.compute(site, fx.potential, {
    nDecoys: fx.config.n_decoys,
    seed: fx.config.seed,
    minSepContact: fx.config.min_sep_contact,
    contactCutoff: fx.config.contact_cutoff,
    flavour: 'configurational'
  });

  check(
    'contact count matches',
    conf.i.length === fx.expected.i.length,
    `js=${conf.i.length} py=${fx.expected.i.length}`
  );

  let orderOk = conf.i.length === fx.expected.i.length;
  if (orderOk) {
    for (let k = 0; k < conf.i.length; k++) {
      if (conf.i[k] !== fx.expected.i[k] || conf.j[k] !== fx.expected.j[k]) {
        orderOk = false;
        console.log(`        first mismatch at ${k}: js=(${conf.i[k]},${conf.j[k]}) py=(${fx.expected.i[k]},${fx.expected.j[k]})`);
        break;
      }
    }
  }
  check('contact order matches', orderOk);

  // --- numerical agreement ----------------------------------------------
  const dens = maxAbsDiff(Array.from(conf.density), fx.expected.density);
  check('local density', dens.diff < TOL_DENSITY, `max|d|=${dens.diff.toExponential(3)}`);

  const en = maxAbsDiff(Array.from(conf.energy), fx.expected.energy);
  check('native contact energy', en.diff < TOL_ENERGY, `max|d|=${en.diff.toExponential(3)}`);

  const fc = maxAbsDiff(Array.from(conf.frustration), fx.expected.configurational);
  check(
    'configurational frustration',
    fc.diff < TOL_FRUSTRATION,
    `max|d|=${fc.diff.toExponential(3)}`
  );

  const mut = engine.compute(site, fx.potential, {
    nDecoys: fx.config.n_decoys,
    seed: fx.config.seed,
    minSepContact: fx.config.min_sep_contact,
    contactCutoff: fx.config.contact_cutoff,
    flavour: 'mutational'
  });
  const fm = maxAbsDiff(Array.from(mut.frustration), fx.expected.mutational);
  check(
    'mutational frustration',
    fm.diff < TOL_FRUSTRATION,
    `max|d|=${fm.diff.toExponential(3)}`
  );

  // --- determinism -------------------------------------------------------
  const again = engine.compute(site, fx.potential, {
    nDecoys: fx.config.n_decoys,
    seed: fx.config.seed,
    minSepContact: fx.config.min_sep_contact,
    contactCutoff: fx.config.contact_cutoff,
    flavour: 'configurational'
  });
  const rep = maxAbsDiff(Array.from(conf.frustration), Array.from(again.frustration));
  check('repeat run is identical', rep.diff === 0, `max|d|=${rep.diff}`);

  console.log(
    failures === 0
      ? '\nall checks passed'
      : `\n${failures} check(s) failed`
  );
  process.exit(failures === 0 ? 0 : 1);
}

main();
