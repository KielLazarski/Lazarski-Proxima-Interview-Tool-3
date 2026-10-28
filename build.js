#!/usr/bin/env node
/* GlueForge build.
 *
 * Assembles src/ into one HTML file. There are two shapes of output and the
 * difference is only whether the validation cache travels inside it:
 *
 *   node build.js                 -> dist/glueforge.html   (served; ~200 kB)
 *   node build.js --standalone    -> dist/glueforge-standalone.html (~4 MB)
 *   node build.js --out=PATH      -> anywhere
 *
 * The served page is what glueforge.py hands to the browser. It fetches the
 * validation set from /api/validation the first time someone asks for it,
 * which keeps the landing screen instant instead of making every start pay
 * for four megabytes nobody has asked to see yet.
 *
 * The standalone page is the export: data inlined, no server needed, opens by
 * double-clicking on a machine with nothing installed. That is the copy you
 * send to someone else.
 *
 * No bundler and no dependencies, on purpose.
 */

'use strict';

const fs = require('fs');
const path = require('path');

const ROOT = __dirname;
const SRC = path.join(ROOT, 'src');
const CACHE = path.join(ROOT, 'cache', 'web');

function arg(name, fallback) {
  const hit = process.argv.find((a) => a.startsWith('--' + name + '='));
  return hit ? hit.split('=').slice(1).join('=') : fallback;
}
const FLAG = (name) => process.argv.includes('--' + name);

function readManifest() {
  const p = path.join(SRC, 'manifest.json');
  if (!fs.existsSync(p)) throw new Error('src/manifest.json is missing');
  return JSON.parse(fs.readFileSync(p, 'utf8'));
}

/** Inline the cache as one object, for the standalone build only. */
function collectData() {
  if (!fs.existsSync(path.join(CACHE, 'index.json'))) {
    console.warn(
      '! cache/web/index.json missing — the standalone page would have no\n' +
      '  validation data in it. Run: python3 sidecar/export_web.py'
    );
    return null;
  }

  const index = JSON.parse(fs.readFileSync(path.join(CACHE, 'index.json'), 'utf8'));
  const systems = {};
  let bytes = 0;

  index.systems.forEach((s) => {
    const f = path.join(CACHE, s.entry + '.json');
    if (!fs.existsSync(f)) {
      console.warn('  ! missing ' + s.entry + '.json, skipping');
      return;
    }
    const raw = fs.readFileSync(f, 'utf8');
    bytes += raw.length;
    systems[s.entry] = JSON.parse(raw);
  });

  // The potential travels with crossval.json; it is the same object the
  // sidecar used, so the page and the reference cannot drift apart.
  let potential = null;
  const cv = path.join(CACHE, 'crossval.json');
  if (fs.existsSync(cv)) {
    potential = JSON.parse(fs.readFileSync(cv, 'utf8')).potential;
  }

  return { index, systems, potential, _bytes: bytes };
}

function main() {
  const standalone = FLAG('standalone');
  const manifest = readManifest();
  let shell = fs.readFileSync(path.join(SRC, 'index.html'), 'utf8');

  const styles = manifest.styles
    .map((rel) => fs.readFileSync(path.join(SRC, rel), 'utf8'))
    .join('\n');

  const missing = manifest.modules.filter(
    (rel) => !fs.existsSync(path.join(SRC, rel))
  );
  if (missing.length) {
    console.error('modules listed in the manifest but not on disk:');
    missing.forEach((m) => console.error('  ' + m));
    process.exit(1);
  }

  const scripts = manifest.modules
    .map((rel) => '/* ===== ' + rel + ' ===== */\n' +
         fs.readFileSync(path.join(SRC, rel), 'utf8'))
    .join('\n\n');

  const data = standalone ? collectData() : null;
  let dataScript = 'globalThis.GF = globalThis.GF || {};\nGF.DATA = ' +
    (data ? JSON.stringify(data) : 'null') + ';';
  // Guard against a closing script tag appearing inside the JSON payload.
  dataScript = dataScript.replace(/<\/script/gi, '<\\/script');

  shell = shell
    .replace('/*__STYLES__*/', () => styles)
    .replace('/*__DATA__*/', () => dataScript)
    // When this page was built, so it can tell whether the server it is
    // talking to started before or after it. Unzipping a new version over a
    // running server is the obvious thing to do and leaves a new page
    // talking to an old process, which fails in ways that look like
    // anything but what they are.
    .replace('/*__SCRIPTS__*/', () =>
      'window.GF=window.GF||{};window.GF.BUILT=' + Math.floor(Date.now() / 1000) + ';\n'
      + scripts);

  const fallback = standalone
    ? path.join('dist', 'glueforge-standalone.html')
    : path.join('dist', 'glueforge.html');
  const out = path.resolve(ROOT, arg('out', fallback));
  fs.mkdirSync(path.dirname(out), { recursive: true });
  fs.writeFileSync(out, shell);

  const kb = (fs.statSync(out).size / 1024).toFixed(0);
  console.log(
    'built ' + path.relative(ROOT, out) + '  (' + kb + ' kB, ' +
    (standalone ? 'standalone' : 'served') + ')'
  );
  if (standalone) {
    if (data && data.index) {
      console.log(
        '  ' + Object.keys(data.systems).length + ' systems inlined, ' +
        'potential ' + data.index.potential +
        (data.index.potential_published ? '' : ' (STAND-IN)')
      );
    }
  } else {
    console.log('  validation data is fetched from the server on demand');
  }
}

main();
