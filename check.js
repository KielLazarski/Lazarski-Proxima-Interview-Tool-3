#!/usr/bin/env node
/* Fast sanity check before a build. Catches the cheap mistakes: a module in
 * the manifest that is not on disk, a module on disk that nobody loads, a
 * syntax error, a cache whose index disagrees with its files. */
'use strict';
const fs = require('fs'), path = require('path'), vm = require('vm');
const ROOT = __dirname, SRC = path.join(ROOT, 'src'), CACHE = path.join(ROOT, 'cache', 'web');
let bad = 0;
const fail = (m) => { console.error('  ! ' + m); bad++; };
const ok = (m) => console.log('  . ' + m);

const manifest = JSON.parse(fs.readFileSync(path.join(SRC, 'manifest.json'), 'utf8'));
manifest.modules.forEach((m) => {
  const p = path.join(SRC, m);
  if (!fs.existsSync(p)) return fail('listed but missing: ' + m);
  try { new vm.Script(fs.readFileSync(p, 'utf8'), { filename: m }); }
  catch (e) { return fail(m + ': ' + e.message); }
});
ok(manifest.modules.length + ' modules parse');

const onDisk = fs.readdirSync(path.join(SRC, 'js')).filter((f) => f.endsWith('.js'));
onDisk.forEach((f) => {
  if (!manifest.modules.includes('js/' + f)) fail('on disk but not in the manifest: js/' + f);
});
ok('manifest and src/js agree');

manifest.styles.forEach((s) => {
  if (!fs.existsSync(path.join(SRC, s))) fail('missing stylesheet: ' + s);
});

const ix = path.join(CACHE, 'index.json');
if (!fs.existsSync(ix)) {
  console.log('  . no cache yet (run python3 sidecar/export_web.py)');
} else {
  const index = JSON.parse(fs.readFileSync(ix, 'utf8'));
  let missing = 0;
  index.systems.forEach((s) => {
    if (!fs.existsSync(path.join(CACHE, s.entry + '.json'))) missing++;
  });
  if (missing) fail(missing + ' system file(s) named in the index but absent');
  else ok(index.systems.length + ' cached systems present');
  if (!index.potential_published) {
    console.log('  . note: cache was built with the stand-in potential, not AWSEM');
  }
}
console.log(bad ? '\n' + bad + ' problem(s)' : '\nok');
process.exit(bad ? 1 : 0);
