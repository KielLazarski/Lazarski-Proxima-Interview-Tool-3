/* GlueForge — talking to the local server.
 *
 * The page runs in two modes and has to know which. Served by glueforge.py it
 * can install packages, parse uploads and dock; opened as an exported file it
 * can only show what was baked into it. Everything here is written so that a
 * standalone copy degrades to "that button is not available" rather than to a
 * silent failure or a spinner that never stops.
 */

(function (GF) {
  'use strict';

  const state = {
    live: false,          // a server answered
    checked: false,
    base: '',
    env: null,            // last /api/status payload
    error: null
  };

  /** fetch + JSON + a usable error message.
   *
   *  The server reports its own problems as {error: "..."} with a real status
   *  code, so the message a user sees is the one the server wrote rather than
   *  "Failed to fetch".
   */
  async function call(path, opts) {
    const o = Object.assign({}, opts);
    let res;
    try {
      res = await fetch(state.base + path, o);
    } catch (err) {
      // The server has gone since the last successful call. Recorded, so
      // that whatever is on screen can show the whole page is offline
      // rather than reporting it as one panel failing to load -- which is
      // what it looks like from inside the panel and is misleading about
      // what is actually wrong.
      state.live = false;
      state.error = 'the server stopped answering';
      throw new ApiError(
        0,
        'The GlueForge server is not answering. It runs in a terminal ' +
        'window and that window has to stay open — if you closed it, ' +
        'start it again by double-clicking GlueForge.bat in the GlueForge ' +
        'folder.'
      );
    }
    const text = await res.text();
    let body = null;
    if (text) {
      try { body = JSON.parse(text); } catch (e) { body = null; }
    }
    if (!res.ok) {
      const msg = (body && body.error) || ('HTTP ' + res.status);
      throw new ApiError(res.status, msg, body || {});
    }
    return body;
  }

  function ApiError(status, message, extra) {
    this.name = 'ApiError';
    this.status = status;
    this.message = message;
    Object.assign(this, extra || {});
  }
  ApiError.prototype = Object.create(Error.prototype);
  ApiError.prototype.constructor = ApiError;

  const json = (path, obj) => call(path, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(obj || {})
  });

  /* ---- environment ----------------------------------------------------- */

  /** Is a server there? Called once at boot; cheap and never throws. */
  async function probe() {
    if (location.protocol === 'file:') {
      state.checked = true;
      state.live = false;
      return state;
    }
    try {
      state.env = await call('/api/status');
      state.live = true;
      state.error = null;
    } catch (err) {
      state.live = false;
      state.error = err.message;
    }
    state.checked = true;
    return state;
  }

  /** Is this page newer than the server answering it?
   *
   *  True when a new version was unzipped over a running one. The page on
   *  disk is replaced immediately; the Python process keeps running the
   *  code it loaded at startup. Everything then half-works: the pages are
   *  new, the routes they call may not exist, and the errors describe
   *  symptoms rather than the cause.
   */
  function serverIsStale() {
    const built = globalThis.GF && globalThis.GF.BUILT;
    const started = state.env && state.env.started;
    if (!built || !started) return false;
    // A minute of slack, because the two clocks are the same clock but the
    // build and the start are not simultaneous even in the normal case.
    return built > started + 60;
  }

  async function refresh() {
    if (!state.live) return state.env;
    state.env = await call('/api/status');
    return state.env;
  }

  /** Does this stage have everything it needs? */
  function stageReady(key) {
    const env = state.env;
    if (!env) return false;
    const st = (env.stages || []).find((s) => s.key === key);
    return !!(st && st.ready);
  }

  function stage(key) {
    const env = state.env;
    if (!env) return null;
    return (env.stages || []).find((s) => s.key === key) || null;
  }

  function pkg(key) {
    const env = state.env;
    if (!env) return null;
    return (env.packages || []).find((p) => p.key === key) || null;
  }

  /* ---- installing ------------------------------------------------------ */

  /** Start an install and stream its output.
   *
   *  `onLine` gets each line as it arrives; `onDone` gets {code, state}. The
   *  stream is server-sent events, with a polling fallback, because an
   *  install that takes four minutes with no output looks identical to one
   *  that has hung — and that is the single most common thing to go wrong on
   *  a fresh machine.
   */
  async function install(key, onLine, onDone) {
    const started = await json('/api/install', { package: key });
    return follow(started.job, onLine, onDone);
  }

  function follow(jobId, onLine, onDone) {
    let settled = false;
    const finish = (payload) => {
      if (settled) return;
      settled = true;
      if (onDone) onDone(payload);
    };

    if (typeof EventSource === 'function') {
      const es = new EventSource(state.base + '/api/job/' + jobId + '/stream');
      es.onmessage = (ev) => {
        try {
          const d = JSON.parse(ev.data);
          if (d.line !== undefined && onLine) onLine(d.line);
        } catch (e) { /* a malformed frame is not worth killing the stream */ }
      };
      es.addEventListener('done', (ev) => {
        let d = {};
        try { d = JSON.parse(ev.data); } catch (e) { /* ignore */ }
        es.close();
        finish(d);
      });
      es.onerror = () => {
        // The connection dropped. Fall back to polling rather than leaving
        // the button spinning: the install itself is still running server-side.
        es.close();
        if (!settled) poll(jobId, onLine, finish);
      };
      return { cancel: () => { es.close(); settled = true; } };
    }

    poll(jobId, onLine, finish);
    return { cancel: () => { settled = true; } };
  }

  function poll(jobId, onLine, finish) {
    let seen = 0;
    const tick = async () => {
      let snap;
      try {
        snap = await call('/api/job/' + jobId);
      } catch (err) {
        return finish({ code: -1, state: 'failed', error: err.message });
      }
      (snap.lines || []).slice(seen).forEach((l) => onLine && onLine(l));
      seen = (snap.lines || []).length;
      if (snap.state === 'running') {
        setTimeout(tick, 900);
      } else {
        finish({ code: snap.code, state: snap.state });
      }
    };
    tick();
  }

  /* ---- runs ------------------------------------------------------------ */

  const runs = {
    create: (label) => json('/api/run', { label: label || null }),
    get: (id) => call('/api/run/' + encodeURIComponent(id)),
    list: () => call('/api/runs'),

    /** Raw bytes rather than multipart: the stdlib has no multipart parser
     *  worth trusting with 20 MB, and there is nothing to send alongside. */
    upload(runId, slot, file) {
      const q = '?slot=' + encodeURIComponent(slot) +
                '&name=' + encodeURIComponent(file.name);
      return call('/api/run/' + encodeURIComponent(runId) + '/upload' + q, {
        method: 'POST',
        headers: { 'Content-Type': 'application/octet-stream' },
        body: file
      });
    },

    face: (runId, body) =>
      json('/api/run/' + encodeURIComponent(runId) + '/face', body)
  };

  /* ---- validation data ------------------------------------------------- */

  /** The 14-system cache.
   *
   *  A standalone export has it inlined as GF.DATA already. A served page
   *  fetches it the first time a validation screen is opened, which keeps the
   *  landing page instant instead of making it wait on four megabytes.
   */
  let validationPromise = null;

  function haveValidation() {
    return !!(GF.DATA && GF.DATA.index);
  }

  function loadValidation() {
    if (haveValidation()) return Promise.resolve(GF.DATA);
    if (validationPromise) return validationPromise;
    if (!state.live) {
      return Promise.reject(new ApiError(
        0,
        'This is an exported copy with no validation data in it, and there ' +
        'is no server to fetch it from.'
      ));
    }
    validationPromise = call('/api/validation').then((d) => {
      GF.DATA = d;
      return d;
    }).catch((err) => {
      validationPromise = null;
      throw err;
    });
    return validationPromise;
  }

  GF.api = {
    state, probe, refresh, call, json, serverIsStale,
    stageReady, stage, pkg,
    install, follow,
    runs,
    haveValidation, loadValidation,
    ApiError
  };
})(globalThis.GF || (globalThis.GF = {}));
