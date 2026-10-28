/* GlueForge — boot and routing.
 *
 * Hash routing, so an exported copy works from file:// with no server.
 *
 * Two kinds of screen live here and they have different needs. The landing
 * and upload screens need a server and no data; the validation screens need
 * the 14-system cache and no server. A served page therefore boots instantly
 * on the landing screen and only pays for the four megabytes of validation
 * data if someone asks for it. An exported copy has that data inlined and no
 * server, and says so.
 *
 * The old hashes (#/benchmark, #/system/6TD3, #/cavity/6TD3, #/method) are
 * kept exactly as they were, because they are linked from inside those pages
 * and from the help text.
 */

(function (GF) {
  'use strict';

  const el = GF.el;

  /** Routes that cannot render without the validation cache. */
  const NEEDS_DATA = new Set(['benchmark', 'system', 'cavity', 'method', 'validation']);

  function parse(hash) {
    const parts = (hash || '').replace(/^#\/?/, '').split('/').filter(Boolean);
    return { head: parts[0] || '', rest: parts.slice(1) };
  }

  /* ---- chrome ---------------------------------------------------------- */

  function navbar(active) {
    const nav = el('nav', { class: 'nav' });

    nav.appendChild(el('a', { class: 'brand', href: '#/' }, [
      el('span', { class: 'mark' }),
      el('span', { text: 'GlueForge' })
    ]));

    const links = el('div', { class: 'links' });
    const link = (href, label, key) => links.appendChild(el('a', {
      href: href,
      class: 'link' + (active === key ? ' on' : ''),
      text: label
    }));

    if (GF.EXPORT && GF.EXPORT.run) {
      link('#/exported', 'This result', 'exported');
      link('#/calibration', 'Thresholds', 'calibration');
    } else {
      link('#/', 'Start', 'landing');
      if (GF.api.state.live) link('#/new', 'Screen a pair', 'upload');
      // Last in the row and after the work, because it is where you go when
      // something is missing rather than where you start.
      if (GF.api.state.live) link('#/setup', 'Setup', 'setup');
    }
    link('#/benchmark', 'Validation', 'benchmark');

    // The system picker only makes sense once the cache is loaded, and only
    // on the screens it can navigate within.
    if (GF.api.haveValidation() && (active === 'benchmark' ||
        active === 'system' || active === 'cavity')) {
      const sel = el('select', { class: 'sel jump' });
      sel.appendChild(el('option', { value: '', text: 'Open a system…' }));
      GF.data.groups().forEach((g) => {
        const og = document.createElement('optgroup');
        og.label = g.label;
        g.rows.forEach((r) => og.appendChild(el('option', {
          value: r.entry, text: r.entry + '  ' + r.ligand
        })));
        sel.appendChild(og);
      });
      sel.addEventListener('change', () => {
        if (!sel.value) return;
        location.hash = (active === 'cavity' ? '#/cavity/' : '#/system/') + sel.value;
      });
      links.appendChild(sel);
    }
    nav.appendChild(links);

    nav.appendChild(statusChip());
    return nav;
  }

  /** The one-glance state: is a server there, and is the potential published. */
  function statusChip() {
    const box = el('div', { class: 'status' });
    const live = GF.api.state.live;

    box.appendChild(el('span', {
      class: 'dot' + (live ? ' ok' : ' warn'),
      title: live
        ? 'Connected to the local GlueForge server.'
        : 'No server: this is an exported copy, so nothing new can be run.'
    }));
    box.appendChild(el('span', {
      class: 'mono',
      text: live ? 'local' : 'exported'
    }));

    const ix = GF.api.haveValidation() ? GF.data.index() : null;
    if (ix) {
      box.appendChild(el('span', { class: 'status-sep', text: '·' }));
      box.appendChild(el('span', {
        class: 'dot' + (ix.potential_published ? ' ok' : ' warn'),
        title: ix.potential_provenance || ''
      }));
      box.appendChild(el('span', {
        class: 'mono',
        text: ix.potential_published ? ix.potential : ix.potential + ' · stand-in'
      }));
    }
    return box;
  }

  function waiting(message) {
    return el('div', { class: 'page' }, [
      el('h1', { text: 'One moment' }),
      el('p', { class: 'lede', text: message })
    ]);
  }

  function problem(title, detail, extra) {
    return el('div', { class: 'page' }, [
      el('h1', { text: title }),
      el('p', { class: 'lede', text: detail }),
      extra || null
    ]);
  }

  /* ---- routing --------------------------------------------------------- */

  let loadingData = false;

  function route() {
    const app = document.getElementById('app');
    const { head, rest } = parse(location.hash);
    const key = head || 'landing';

    // Validation screens need the cache. Fetch it once, then come back here.
    if (NEEDS_DATA.has(key) && !GF.api.haveValidation()) {
      app.innerHTML = '';
      app.appendChild(navbar(key));
      if (loadingData) {
        app.appendChild(waiting('Loading the validation set.'));
        return;
      }
      loadingData = true;
      app.appendChild(waiting(
        'Loading the 14 benchmark systems. This happens once per session.'
      ));
      GF.api.loadValidation().then(() => {
        loadingData = false;
        route();
      }).catch((err) => {
        loadingData = false;
        app.innerHTML = '';
        app.appendChild(navbar(key));
        app.appendChild(problem(
          'The validation set is not available',
          err.message,
          el('p', { class: 'note' }, [
            'Rebuild it with ',
            el('code', { text: 'python3 sidecar/export_web.py' }),
            '.'
          ])
        ));
      });
      return;
    }

    let view, active = key;
    try {
      if (key === 'landing') {
        view = GF.pages.landing();
      } else if (key === 'setup') {
        view = GF.pages.setup();
      } else if (key === 'new') {
        active = 'upload';
        view = requireServer() || GF.pages.upload(null);
      } else if (key === 'run') {
        active = 'upload';
        view = requireServer() || GF.pages.upload(rest[0]);
      } else if (key === 'assess' && rest[0]) {
        active = 'upload';
        view = requireServer() || GF.pages.assess(rest[0]);
      } else if (key === 'dock' && rest[0]) {
        active = 'upload';
        view = requireServer() || GF.pages.dock(rest[0]);
      } else if (key === 'calibration') {
        active = 'benchmark';
        view = GF.EXPORT && GF.EXPORT.calibration
          ? GF.pages.calibration(rest[0] || null)
          : (requireServer() || GF.pages.calibration(rest[0] || null));
      } else if (key === 'chemistry' && rest[0]) {
        active = 'upload';
        view = requireServer() || GF.pages.chemistry(rest[0]);
      } else if (key === 'exported') {
        active = 'exported';
        view = GF.pages.exported();
      } else if (key === 'cavity' && rest[0]) {
        view = GF.pages.cavity(rest[0].toUpperCase());
      } else if (key === 'system' && rest[0]) {
        view = GF.pages.system(rest[0].toUpperCase());
      } else if (key === 'method') {
        view = GF.pages.method();
      } else if (key === 'benchmark' || key === 'validation') {
        active = 'benchmark';
        view = GF.pages.benchmark();
      } else {
        view = problem(
          'No such screen',
          'Nothing is routed at ' + location.hash + '.',
          el('a', { class: 'btn', href: '#/', text: 'Back to the start' })
        );
      }
    } catch (err) {
      view = el('div', { class: 'page' }, [
        el('h1', { text: 'Something went wrong rendering this screen' }),
        el('pre', { class: 'err', text: String((err && err.stack) || err) })
      ]);
      if (window.console) console.error(err);
    }

    app.innerHTML = '';
    app.appendChild(navbar(active));
    // Above everything, because when it applies it explains most of what
    // else is about to look broken.
    const stale = staleBanner();
    if (stale) app.appendChild(stale);
    app.appendChild(GF.wrapTables(view));
    window.scrollTo(0, 0);
  }

  /** The page is newer than the process answering it.
   *
   *  Unzipping a new version over a running one is the obvious thing to do
   *  and leaves exactly this: new pages, old routes. The symptoms are
   *  scattered and none of them name the cause -- a panel that will not
   *  load, a button that 404s -- so it is said once, at the top, in terms
   *  of the thing to actually do.
   */
  function staleBanner() {
    if (!GF.api.serverIsStale || !GF.api.serverIsStale()) return null;
    return el('div', { class: 'callout warn stale-banner' }, [
      el('strong', { text: 'The running server is older than this page. ' }),
      'GlueForge was updated on disk while the server was running, so it is '
      + 'still serving the code it started with. Stop it in the window it '
      + 'is running in and start it again with ',
      el('code', { text: 'py glueforge.py' }),
      '. Anything that fails to load until then is a symptom of this.'
    ]);
  }

  /** Screens that need the server, when there is none. */
  function requireServer() {
    if (GF.api.state.live) return null;
    return problem(
      'This copy cannot run anything new',
      'Uploading, docking and chemistry all need the local server. This is ' +
      'an exported page, so it has the results it was exported with and ' +
      'nothing else.',
      el('p', { class: 'note' }, [
        'Start it by double-clicking ',
        el('code', { text: 'GlueForge.bat' }),
        ' in the GlueForge folder, or by running ',
        el('code', { text: 'py glueforge.py' }),
        ' there. Leave that window open — it is the server, and closing it '
        + 'is what makes this page stop answering.'
      ])
    );
  }

  /* ---- boot ------------------------------------------------------------ */

  async function boot() {
    const app = document.getElementById('app');
    app.innerHTML = '';
    app.appendChild(waiting('Starting up.'));

    await GF.api.probe();

    // An exported assessment carries its own run and calibration. Open it
    // where it was left rather than on a landing page whose buttons do
    // nothing here.
    if (GF.EXPORT && GF.EXPORT.run && !location.hash) {
      location.hash = '#/exported';
    }

    window.addEventListener('hashchange', route);
    route();
  }

  GF.route = route;
  GF.boot = boot;

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', boot);
  } else {
    boot();
  }
})(globalThis.GF || (globalThis.GF = {}));
