/* GlueForge — the landing page.
 *
 * This replaces the terminal menu. It is the first thing you see, and its job
 * is to answer three questions without being asked: what is this, what is
 * missing, and what can I do right now.
 *
 * The design rule here is that nothing is hidden behind a disabled button
 * with no explanation. Every stage that is not available says which package
 * it is waiting for and what that package is for, and the install button for
 * it is right there in the same row.
 */

(function (GF) {
  'use strict';

  const el = GF.el;
  const api = GF.api;

  /* ---- what the tool is, in the fewest words that are still true ------- */

  function intro() {
    return el('div', { class: 'page' }, [
      el('h1', { text: 'GlueForge' }),
      el('p', { class: 'lede' }, [
        'Two proteins and one small molecule. The molecule is already bound ' +
        'to the first protein; the question is whether a second protein can ' +
        'be brought to it, and what chemistry would have to be added to the ' +
        'molecule to hold them together.'
      ]),
      el('p', { class: 'sub' }, [
        'The part of the molecule still exposed once it sits in its pocket is ' +
        'the only surface a second protein can reach. GlueForge measures that ' +
        'face, uses it to aim a docking run, judges whether the resulting ' +
        'interface has a pocket worth designing into, and only then grows ' +
        'chemistry out from the molecule into it.'
      ])
    ]);
  }

  /* ---- reading what pip actually said ----------------------------------- */

  /** Turn a pip failure into the thing that would fix it.
   *
   *  "If this is a permissions error, run it yourself" was the old message
   *  and it is close to useless: it guesses at one cause and then suggests
   *  re-running the command that just failed. pip says why it failed; this
   *  reads that and names the cause.
   *
   *  The common one on Windows is a package with no wheel for the running
   *  Python, which makes pip fall back to compiling from source and fail for
   *  want of a C compiler. The fix is almost never to install a compiler.
   */
  function diagnose(text) {
    const py = (api.state.env && api.state.env.python) || {};
    const ver = py.version || '';
    // "Python ." is what this prints if the status has not arrived, which
    // looks like a broken page rather than a diagnosis.
    const short = ver ? ver.split('.').slice(0, 2).join('.') : '';
    const named = short ? 'Python ' + short : 'this Python';

    const noWheel = /Failed building wheel|Failed to build installable wheels|error: Microsoft Visual C\+\+|building .* extension/i.test(text);
    const names = [];
    const m = text.match(/Failed to build(?: installable wheels for some [^(]*\()?\s*([^)\n]+)/i);
    if (m) {
      m[1].split(/[,\s]+/).forEach((n) => {
        const t = n.replace(/[^A-Za-z0-9._-]/g, '');
        if (t && !names.includes(t)) names.push(t);
      });
    }

    if (noWheel) {
      // LightDock is the case this was written for, and it is the awkward
      // one. It publishes no Windows wheel for any Python -- its own PyPI
      // page says Windows is not officially supported -- and its fastdfire
      // scoring function is compiled C. So "use an older Python" is not a
      // fix for it at any version, and the first version of this message
      // said 3.12, which is wrong twice over: freesasa's newest Windows
      // wheel is for 3.11. Checked against PyPI rather than guessed.
      const windows = (py.platform || '').toLowerCase().indexOf('win') === 0;
      const dock = names.indexOf('lightdock') >= 0;

      const node = el('div', { class: 'callout warn diag' }, [
        el('strong', { text: 'Nothing prebuilt for ' + named + '. ' }),
        (names.length
          ? names.join(', ') + ' ship compiled code. '
          : 'This package ships compiled code. ') +
        'pip found no ready-made build, fell back to compiling from source, ' +
        'and there is no C compiler here.'
      ]);

      if (windows && dock) {
        node.appendChild(el('p', { class: 'sub' }, [
          'A compiler alone will not finish this. LightDock publishes no ' +
          'Windows build for any Python version and does not list Windows ' +
          'as supported, so on Windows the reliable route is ',
          el('strong', { text: 'WSL' }),
          ', where it is built the way its authors build it:'
        ]));
        node.appendChild(el('pre', { class: 'cmd', text:
          'wsl --install -d Ubuntu          (once, then restart)\n' +
          'sudo apt update\n' +
          'sudo apt install -y python3-venv python3-dev build-essential\n' +
          'cd /mnt/c/path/to/GlueForge\n' +
          'python3 -m venv .venv && . .venv/bin/activate\n' +
          'python glueforge.py --no-browser' }));
        node.appendChild(el('p', { class: 'fnote', text:
          'WSL forwards localhost, so the address it prints opens in your ' +
          'normal Windows browser. GlueForge itself is pure Python plus ' +
          'numpy and scipy, so it does not care which side it runs on, and ' +
          'your files stay where they are under /mnt/c.' }));
        node.appendChild(el('p', { class: 'sub' }, [
          el('strong', { text: 'Or skip it. ' }),
          'Docking is the only stage that needs LightDock. Assessing a ' +
          'complex you already have, and growing chemistry into it, need ' +
          'only numpy, scipy and RDKit — all of which have Windows builds ' +
          'for the Python you are running.'
        ]));
      } else if (windows) {
        node.appendChild(el('p', { class: 'sub', text:
          'Find out whether a build exists for any Python before installing ' +
          'a compiler. Run these one at a time — pip stops at the first ' +
          'failure, so a single line only reports on the first package:' }));
        node.appendChild(el('pre', { class: 'cmd', text:
          (names.length ? names : ['<package>'])
            .map((n) => 'py -m pip download --no-deps --only-binary=:all: ' + n)
            .join('\n') }));
        node.appendChild(el('p', { class: 'sub', text:
          'If one of them downloads for an older Python, run GlueForge on ' +
          'that Python instead. From this folder, with the version that has ' +
          'a build:' }));
        node.appendChild(el('pre', { class: 'cmd', text:
          'py -3.11 -m venv .venv\n' +
          '.venv\\Scripts\\activate\n' +
          'python glueforge.py' }));
        node.appendChild(el('p', { class: 'fnote', text:
          'If nothing downloads at any version, the package has no Windows ' +
          'build at all and WSL is the shorter road than a compiler.' }));
      } else {
        node.appendChild(el('p', { class: 'sub', text:
          'On Linux these build from source in seconds once a compiler is ' +
          'present, which is one command:' }));
        node.appendChild(el('pre', { class: 'cmd', text:
          'sudo apt install -y python3-dev build-essential' }));
      }

      node.appendChild(el('p', { class: 'fnote', text:
        'GlueForge installs into whichever Python started it, so once you ' +
        'are in the right environment the buttons on this page work inside ' +
        'it. Nothing else needs changing.' }));
      return node;
    }

    if (/Permission denied|Access is denied|EnvironmentError/i.test(text)) {
      return el('div', { class: 'callout warn diag' }, [
        el('strong', { text: 'Permissions. ' }),
        'pip could not write where it wanted to. Running GlueForge from a ' +
        'virtual environment avoids this entirely:',
        el('pre', { class: 'cmd', text:
          'py -m venv .venv\n.venv\\Scripts\\activate\npython glueforge.py' })
      ]);
    }

    if (/No matching distribution|Could not find a version/i.test(text)) {
      return el('div', { class: 'callout warn diag' }, [
        el('strong', { text: 'Nothing published for this Python. ' }),
        'pip reached the index and found no release matching ' + named +
        '. An older Python, 3.12 for preference, is the fix.'
      ]);
    }

    if (/Temporary failure|getaddrinfo|ProxyError|SSLError|Connection/i.test(text)) {
      return el('div', { class: 'callout warn diag' }, [
        el('strong', { text: 'Could not reach the package index. ' }),
        'This looks like a network or proxy problem rather than anything ' +
        'about the package. If you are behind a corporate proxy, pip needs ' +
        'to be told about it.'
      ]);
    }
    return null;
  }

  /* ---- dependency rows ------------------------------------------------- */

  function packageRow(p, onChanged) {
    const row = el('div', { class: 'dep' + (p.present ? ' on' : '') });

    row.appendChild(el('div', { class: 'dep-dot' }, [
      el('span', { class: 'dot' + (p.present ? ' ok' : (p.required ? ' bad' : ' warn')) })
    ]));

    const main = el('div', { class: 'dep-main' }, [
      el('div', { class: 'dep-name' }, [
        el('code', { text: p.key }),
        p.version
          ? el('span', { class: 'dep-ver mono', text: p.version })
          : el('span', {
              class: 'dep-tag',
              text: p.required ? 'required' : 'needed for ' + p.stage
            })
      ]),
      el('div', { class: 'dep-why', text: p.why })
    ]);
    row.appendChild(main);

    const act = el('div', { class: 'dep-act' });
    if (p.present) {
      act.appendChild(el('span', { class: 'dep-ok', text: 'installed' }));
    } else {
      const log = el('pre', { class: 'dep-log' });
      const btn = el('button', { class: 'btn', text: 'Install' });
      btn.addEventListener('click', () => {
        btn.disabled = true;
        btn.textContent = 'Installing…';
        row.classList.add('busy');
        log.classList.add('open');
        log.textContent = '';
        api.install(p.key, (line) => {
          log.textContent += line + '\n';
          log.scrollTop = log.scrollHeight;
        }, (done) => {
          row.classList.remove('busy');
          if (done && done.code === 0) {
            btn.textContent = 'Installed';
            // The row only turns green once the server confirms the module
            // actually imports. pip exiting 0 is not the same thing: a wheel
            // can land for the wrong interpreter and import will still fail.
            onChanged();
          } else {
            btn.disabled = false;
            btn.textContent = 'Try again';
            log.textContent +=
              '\n--- install failed' +
              (done && done.code !== undefined ? ' (exit ' + done.code + ')' : '') +
              ' ---\n';
            log.scrollTop = log.scrollHeight;
            // Telling someone to re-run the command that just failed is not
            // advice. Read what pip actually said and name the cause.
            const why = diagnose(log.textContent);
            if (why) main.insertBefore(why, log);
          }
        });
      });
      act.appendChild(btn);
      main.appendChild(log);
    }
    row.appendChild(act);
    return row;
  }

  /* ---- stage buttons --------------------------------------------------- */

  const STAGE_ROUTE = {
    assess: '#/new',
    dock: '#/new',
    chemistry: '#/new',
    validation: '#/validation'
  };

  function stageCard(st) {
    const ready = st.ready;
    const node = el(ready ? 'a' : 'div', {
      class: 'stage' + (ready ? '' : ' off'),
      href: ready ? (STAGE_ROUTE[st.key] || '#/') : null
    }, [
      el('div', { class: 'stage-label', text: st.label }),
      el('div', { class: 'stage-note', text: st.note })
    ]);
    if (!ready) {
      node.appendChild(el('div', { class: 'stage-blocked' }, [
        'waiting on ',
        el('code', { text: st.missing.join(', ') })
      ]));
    } else if (st.better_with_missing && st.better_with_missing.length) {
      // A stage that runs, with an alternative available. Said differently
      // from blocked on purpose: docking works without LightDock, so
      // showing the two the same way sends people after an install they do
      // not need. And not called "stronger" -- the two engines have never
      // been compared on this benchmark, so that word would be a number
      // nobody measured.
      node.appendChild(el('div', { class: 'stage-extra' }, [
        'runs on the built-in engine · ',
        el('code', { text: st.better_with_missing.join(', ') }),
        ' is a second engine, not compared with it yet'
      ]));
    }
    return node;
  }

  /* ---- recent runs ----------------------------------------------------- */

  function runList(runs) {
    if (!runs || !runs.length) return null;
    const tbl = el('table', { class: 'grid runs' });
    tbl.appendChild(el('thead', {}, [
      el('tr', {}, [
        el('th', { text: 'Run' }),
        el('th', { text: 'Protein A' }),
        el('th', { text: 'Protein B' }),
        el('th', { text: 'Stage' })
      ])
    ]));
    const body = el('tbody');
    runs.forEach((r) => {
      const a = (r.slots && r.slots.A) || {};
      const b = (r.slots && r.slots.B) || {};
      body.appendChild(el('tr', {}, [
        el('td', {}, [
          el('a', {
            href: '#/run/' + r.id,
            class: 'mono',
            text: r.label || r.id
          })
        ]),
        el('td', { class: 'mono', text: a.entry || a.name || '—' }),
        el('td', { class: 'mono', text: b.entry || b.name || '—' }),
        el('td', { text: r.stage || 'upload' })
      ]));
    });
    tbl.appendChild(body);
    return el('section', { class: 'panel' }, [
      el('h2', { text: 'Earlier runs' }),
      el('p', { class: 'sub' }, [
        'Every run is a folder under ',
        el('code', { text: 'runs/' }),
        ', holding the structures you gave it, the poses, the scores and the ' +
        'settings that produced them. Nothing is kept only in the browser.'
      ]),
      tbl
    ]);
  }

  /* ---- the standalone case --------------------------------------------- */

  function offlineNotice() {
    return el('section', { class: 'panel warn' }, [
      el('h2', { text: 'This is an exported copy' }),
      el('p', {}, [
        'It holds a finished analysis and the validation set, and it will ' +
        'keep working on any machine with no Python and nothing installed. ' +
        'What it cannot do is run anything new: no uploads, no docking, no ' +
        'chemistry.'
      ]),
      el('p', { class: 'note' }, [
        'To run something new, start the tool itself — ',
        el('code', { text: 'py glueforge.py' }),
        ' in the GlueForge folder — and this page will come back with its ' +
        'buttons live.'
      ])
    ]);
  }

  /* ---- assembly -------------------------------------------------------- */

  function landing() {
    const host = el('div', {});
    host.appendChild(intro());

    if (!api.state.live) {
      host.appendChild(offlineNotice());
      const stages = el('div', { class: 'stages' });
      stages.appendChild(stageCard({
        key: 'validation', label: 'Browse validation', ready: true,
        note: 'The 14 benchmark systems and the error budget.', missing: []
      }));
      host.appendChild(el('section', { class: 'panel' }, [stages]));
      return host;
    }

    const env = api.state.env || {};

    // What you can do, first and large. This used to come after the
    // dependency list, which put a wall of install buttons between someone
    // opening the tool and the thing they opened it to do -- and most of
    // the time nothing needs installing at all.
    const stages = el('div', { class: 'stages' });
    (env.stages || []).forEach((st) => stages.appendChild(stageCard(st)));
    host.appendChild(el('section', { class: 'panel' }, [
      el('h2', { text: 'What you can do' }),
      stages
    ]));

    // The setup line. It only becomes loud when something is genuinely
    // blocked: a stage that cannot run is worth interrupting for, and a
    // stage that would merely run better is not.
    const blocked = (env.stages || []).filter((st) => !st.ready);
    const missing = [];
    blocked.forEach((st) => (st.missing || []).forEach((m) => {
      if (missing.indexOf(m) < 0) missing.push(m);
    }));

    if (missing.length) {
      host.appendChild(el('div', { class: 'callout warn setup-line' }, [
        el('strong', { text: missing.join(' and ') + ' not installed. ' }),
        'That leaves ' +
        blocked.map((s2) => s2.label.toLowerCase()).join(' and ') +
        ' unavailable. ',
        el('a', { href: '#/setup', text: 'Install from Setup' }),
        '.'
      ]));
    } else {
      host.appendChild(el('p', { class: 'fnote setup-line' }, [
        'Everything these need is installed. ',
        el('a', { href: '#/setup', text: 'Setup' }),
        ' has the package list, the Python it is running on, and the '
        + 'installer output if you ever need it.'
      ]));
    }

    const runs = runList(env.runs);
    if (runs) host.appendChild(runs);

    return host;
  }

  /* ---- setup, on its own page ------------------------------------------ */

  /** The dependency list, moved off the front.
   *
   *  It belongs somewhere -- an install that fails needs its output read,
   *  and the diagnosis of a failed install is the most useful text in this
   *  whole interface to whoever hits it. It just does not belong in front
   *  of someone who came here to dock two proteins and whose installation
   *  is already fine.
   */
  function setup() {
    const host = el('div', {});
    host.appendChild(el('div', { class: 'page' }, [
      el('h1', { text: 'Setup' }),
      el('p', { class: 'lede' }, [
        'What GlueForge runs on, and what it is missing. Everything here is '
        + 'optional except numpy and scipy: the stages that need a package '
        + 'say so, and the ones that do not are unaffected by it.'
      ])
    ]));

    if (!api.state.live) {
      host.appendChild(offlineNotice());
      return host;
    }

    const env = api.state.env || {};
    const redraw = async () => {
      try {
        await api.refresh();
      } catch (err) { /* keep whatever we last knew */ }
      GF.route();
    };

    const deps = el('section', { class: 'panel' }, [
      el('h2', { text: 'Packages' }),
      el('p', { class: 'sub' }, [
        'The output appears in the page as it installs, so a slow install '
        + 'is distinguishable from a stuck one. If one fails, what it '
        + 'printed is read here and turned into the thing that would fix it.'
      ])
    ]);
    const list = el('div', { class: 'deps' });
    (env.packages || []).forEach((p2) => list.appendChild(packageRow(p2, redraw)));
    deps.appendChild(list);
    host.appendChild(deps);

    const st = el('section', { class: 'panel' }, [
      el('h2', { text: 'What that leaves available' })
    ]);
    const stages = el('div', { class: 'stages' });
    (env.stages || []).forEach((x) => stages.appendChild(stageCard(x)));
    st.appendChild(stages);
    host.appendChild(st);

    const py = env.python || {};
    host.appendChild(el('section', { class: 'panel' }, [
      el('h2', { text: 'The Python this is running on' }),
      el('p', { class: 'sub', text:
        'GlueForge installs into whichever Python started it, so this is the '
        + 'one the buttons above are installing into.' }),
      el('p', { class: 'fnote mono' }, [
        'Python ' + (py.version || '?') + ' · ' + (py.platform || '') +
        ' · ' + (py.machine || '')
      ]),
      el('p', { class: 'fnote mono', text: py.executable || '' })
    ]));

    return host;
  }

  GF.pages = GF.pages || {};
  GF.pages.landing = landing;
  GF.pages.setup = setup;
  // Exported so the diagnosis can be tested against real pip output rather
  // than only seen by whoever happens to hit the failure.
  GF.landing = { diagnose: diagnose };
})(globalThis.GF || (globalThis.GF = {}));
