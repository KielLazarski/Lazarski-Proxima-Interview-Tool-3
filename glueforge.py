#!/usr/bin/env python3
"""GlueForge.

Start this and a page opens in your browser. Everything else happens there:
checking what is installed, installing what is missing, uploading two
structures, docking them, assessing the result, and growing chemistry into it.

    py glueforge.py            Windows
    python3 glueforge.py       macOS and Linux

This file deliberately imports nothing beyond the standard library at module
scope. numpy, scipy, LightDock and RDKit are all installed *from* the landing
page, which means the server has to be able to run before any of them exist.
Every import of those lives inside the function that needs it.
"""

from __future__ import annotations

import http.server
import importlib
import importlib.util
import json
import mimetypes
import os
import platform
import queue
import re
import shutil
import socket
import socketserver
import subprocess
import sys
import threading
import time
import traceback
import uuid
import webbrowser
from pathlib import Path
from urllib.parse import parse_qs, unquote, urlparse

HERE = Path(__file__).resolve().parent
SIDECAR = HERE / "sidecar"
RUNS = HERE / "runs"
CACHE = HERE / "cache"
WEB_CACHE = CACHE / "web"

VERSION = "2.0.0-dev"

#: When this process started. Reported so the page can tell whether it is
#: talking to a server older than itself -- which is what happens when a new
#: version is unzipped over a running one, and it fails in ways that look
#: like anything except what they are.
STARTED = time.time()
DEFAULT_PORT = 8731
MAX_UPLOAD = 192 * 1024 * 1024  # a 200k-atom mmCIF is about 20 MB

# Where the page comes from. dist/ is what build.js writes; the copy at the
# root is what ships in the zip. Either works.
PAGE_CANDIDATES = (HERE / "dist" / "glueforge.html", HERE / "glueforge.html")

#: The build with the validation cache inlined. An exported assessment is
#: built from this so that "why this cut-off" still works on a machine with
#: no server, which is the whole point of exporting it.
STANDALONE_CANDIDATES = (
    HERE / "dist" / "glueforge-standalone.html",
    HERE / "glueforge-standalone.html",
)


# ---------------------------------------------------------------------------
# Dependencies
# ---------------------------------------------------------------------------

#: Each row is one line on the landing page. `stage` is what you lose without
#: it, which is the only thing that makes "optional" a useful word.
PACKAGES = (
    {
        "key": "numpy",
        "module": "numpy",
        "pip": "numpy",
        "required": True,
        "stage": "everything",
        "why": "All the geometry. Nothing analytical runs without it.",
    },
    {
        "key": "scipy",
        "module": "scipy",
        "pip": "scipy",
        "required": True,
        "stage": "everything",
        "why": "Cavity detection and the statistics behind every threshold.",
    },
    {
        "key": "lightdock",
        "module": "lightdock",
        "pip": "lightdock",
        "required": False,
        "stage": "docking",
        "why": (
            "The better of the two docking engines. GlueForge has its own, "
            "which needs nothing and handles small partners well but "
            "degrades on large ones; LightDock handles both. It does not "
            "install on Windows — see the README for the WSL route."
        ),
    },
    {
        "key": "rdkit",
        "module": "rdkit",
        "pip": "rdkit",
        "required": False,
        "stage": "chemistry",
        "why": (
            "Building molecules out from the ligand. Without it you can "
            "still dock, assess, and see where growth would go."
        ),
    },
)

STAGES = (
    {
        "key": "assess",
        "label": "Assess a complex",
        "needs": ["numpy", "scipy"],
        "note": "Upload a ternary complex and get the gate verdict.",
    },
    {
        "key": "dock",
        "label": "Dock a new pair",
        "needs": ["numpy", "scipy"],
        "better_with": ["lightdock"],
        "note": "Upload two separate proteins and place them.",
    },
    {
        "key": "chemistry",
        "label": "Grow chemistry",
        "needs": ["numpy", "scipy", "rdkit"],
        "note": "Build molecules from the ligand into the partner.",
    },
    {
        "key": "validation",
        "label": "Browse validation",
        "needs": [],
        "note": "The 14 benchmark systems. Works with nothing installed.",
    },
)


def module_version(mod_name: str) -> str | None:
    """Import just far enough to read a version string, or give up quietly."""
    try:
        mod = importlib.import_module(mod_name)
    except Exception:  # noqa: BLE001 - a broken install must not kill status
        return None
    for attr in ("__version__", "version", "VERSION"):
        v = getattr(mod, attr, None)
        if isinstance(v, str):
            return v
        if v is not None and hasattr(v, "__str__") and not callable(v):
            return str(v)
    return "installed"


def have(mod_name: str) -> bool:
    """Present *and* importable.

    find_spec alone is not enough: a half-finished pip install leaves a spec
    behind that blows up on import, and the landing page would then show a
    green row for something that cannot be used.
    """
    try:
        if importlib.util.find_spec(mod_name) is None:
            return False
    except Exception:  # noqa: BLE001 - a namespace clash raises here
        return False
    try:
        importlib.import_module(mod_name)
        return True
    except Exception:  # noqa: BLE001
        return False


def environment() -> dict:
    packages = []
    for spec in PACKAGES:
        present = have(spec["module"])
        packages.append({
            **{k: spec[k] for k in ("key", "pip", "required", "stage", "why")},
            "present": present,
            "version": module_version(spec["module"]) if present else None,
        })
    by_key = {p["key"]: p["present"] for p in packages}

    stages = []
    for st in STAGES:
        missing = [n for n in st["needs"] if not by_key.get(n)]
        # Separate from `missing` on purpose. A stage that cannot run and a
        # stage that would run better are different things to tell someone,
        # and showing them the same way sends people after an install they
        # do not need.
        better = [n for n in st.get("better_with", ()) if not by_key.get(n)]
        stages.append({**st, "ready": not missing, "missing": missing,
                       "better_with_missing": better})

    return {
        "version": VERSION,
        "started": STARTED,
        "python": {
            "version": "%d.%d.%d" % sys.version_info[:3],
            "executable": sys.executable,
            "platform": platform.system(),
            "machine": platform.machine(),
            "ok": sys.version_info >= (3, 9),
        },
        "packages": packages,
        "stages": stages,
        "validation": WEB_CACHE.joinpath("index.json").exists(),
        "root": str(HERE),
    }


# ---------------------------------------------------------------------------
# Background jobs
# ---------------------------------------------------------------------------


class Job:
    """A subprocess whose output the page watches while it runs.

    The point of streaming rather than waiting is that pip on a cold cache can
    sit silent for two minutes building a wheel, and a button that has been
    spinning for two minutes is indistinguishable from a button that has hung.
    """

    def __init__(self, label: str, argv: list[str]) -> None:
        self.id = uuid.uuid4().hex[:12]
        self.label = label
        self.argv = argv
        self.lines: list[str] = []
        self.state = "running"
        self.code: int | None = None
        self.started = time.time()
        self.finished: float | None = None
        self._lock = threading.Lock()
        self._watchers: list[queue.Queue] = []

    # -- writing --

    def emit(self, line: str) -> None:
        with self._lock:
            self.lines.append(line)
            # A runaway build log should not become a memory leak.
            if len(self.lines) > 4000:
                del self.lines[:1000]
                self.lines.insert(0, "… earlier output trimmed …")
            watchers = list(self._watchers)
        for w in watchers:
            try:
                w.put_nowait(line)
            except queue.Full:
                pass

    def finish(self, code: int) -> None:
        with self._lock:
            self.code = code
            self.state = "done" if code == 0 else "failed"
            self.finished = time.time()
            watchers = list(self._watchers)
        for w in watchers:
            try:
                w.put_nowait(None)
            except queue.Full:
                pass

    # -- reading --

    def watch(self) -> queue.Queue:
        q: queue.Queue = queue.Queue(maxsize=2048)
        with self._lock:
            for line in self.lines:
                try:
                    q.put_nowait(line)
                except queue.Full:
                    break
            if self.state != "running":
                q.put_nowait(None)
            else:
                self._watchers.append(q)
        return q

    def unwatch(self, q: queue.Queue) -> None:
        with self._lock:
            if q in self._watchers:
                self._watchers.remove(q)

    def snapshot(self) -> dict:
        with self._lock:
            return {
                "id": self.id,
                "label": self.label,
                "state": self.state,
                "code": self.code,
                "lines": list(self.lines),
                "elapsed": round((self.finished or time.time()) - self.started, 1),
            }


class Jobs:
    def __init__(self) -> None:
        self._jobs: dict[str, Job] = {}
        self._lock = threading.Lock()

    def get(self, job_id: str) -> Job | None:
        with self._lock:
            return self._jobs.get(job_id)

    def start(self, label: str, argv: list[str], cwd: Path | None = None) -> Job:
        job = Job(label, argv)
        with self._lock:
            self._jobs[job.id] = job
            # Keep the last 40; the page only ever looks at the live one.
            if len(self._jobs) > 40:
                oldest = sorted(self._jobs.values(), key=lambda j: j.started)
                for j in oldest[:10]:
                    if j.state != "running":
                        self._jobs.pop(j.id, None)

        def run() -> None:
            job.emit("$ " + " ".join(argv))
            try:
                env = dict(os.environ)
                env["PYTHONUNBUFFERED"] = "1"
                proc = subprocess.Popen(
                    argv,
                    cwd=str(cwd or HERE),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    bufsize=1,
                    env=env,
                    # Stop a console window flashing up on Windows.
                    creationflags=(
                        getattr(subprocess, "CREATE_NO_WINDOW", 0)
                        if os.name == "nt" else 0
                    ),
                )
            except OSError as exc:
                job.emit("could not start: %s" % exc)
                job.finish(127)
                return

            assert proc.stdout is not None
            for raw in proc.stdout:
                job.emit(raw.rstrip("\n"))
            proc.wait()
            job.finish(proc.returncode)

        threading.Thread(target=run, daemon=True, name="job-" + job.id).start()
        return job


JOBS = Jobs()


def pip_install(pip_name: str) -> Job:
    argv = [sys.executable, "-m", "pip", "install", "--upgrade", pip_name]
    return JOBS.start("install " + pip_name, argv)


# ---------------------------------------------------------------------------
# Runs
# ---------------------------------------------------------------------------

SAFE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,63}$")


class BadRequest(ValueError):
    """Input that never had a chance of being valid.

    Distinct from a plain ValueError so the dispatcher can answer 400
    without also turning genuine internal errors into client errors. A run
    id with a slash in it reached the generic handler before this existed
    and came back 500 with a traceback in the log -- which reads as the
    server breaking rather than the request being refused, and buries real
    failures among the probes.
    """


def new_run_id() -> str:
    return time.strftime("%Y%m%d-%H%M%S") + "-" + uuid.uuid4().hex[:4]


def run_dir(run_id: str) -> Path:
    if not SAFE_ID.match(run_id):
        raise BadRequest("bad run id")
    return RUNS / run_id


def load_run(run_id: str) -> dict:
    p = run_dir(run_id) / "run.json"
    if not p.exists():
        raise FileNotFoundError(run_id)
    return json.loads(p.read_text(encoding="utf-8"))


def save_run(state: dict) -> None:
    d = run_dir(state["id"])
    d.mkdir(parents=True, exist_ok=True)
    state["updated"] = time.time()
    tmp = d / "run.json.tmp"
    tmp.write_text(json.dumps(state, indent=2), encoding="utf-8")
    tmp.replace(d / "run.json")


def list_runs(limit: int = 50) -> list[dict]:
    if not RUNS.exists():
        return []
    out = []
    for d in sorted(RUNS.iterdir(), reverse=True):
        f = d / "run.json"
        if not f.is_file():
            continue
        try:
            st = json.loads(f.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        out.append({
            "id": st.get("id", d.name),
            "created": st.get("created"),
            "updated": st.get("updated"),
            "stage": st.get("stage", "upload"),
            "label": st.get("label") or d.name,
            "slots": {
                k: {"name": v.get("name"), "entry": v.get("entry")}
                for k, v in (st.get("slots") or {}).items()
            },
        })
        if len(out) >= limit:
            break
    return out


# ---------------------------------------------------------------------------
# Structure parsing (needs numpy, so everything here imports lazily)
# ---------------------------------------------------------------------------


def _sidecar_on_path() -> None:
    s = str(SIDECAR)
    if s not in sys.path:
        sys.path.insert(0, s)


class NeedsPackages(Exception):
    def __init__(self, missing: list[str]) -> None:
        super().__init__("missing: " + ", ".join(missing))
        self.missing = missing


def require(*mods: str) -> None:
    missing = [m for m in mods if not have(m)]
    if missing:
        raise NeedsPackages(missing)


def summarise_structure(path: Path, display_name: str | None = None) -> dict:
    """What the upload screen shows the moment a file lands.

    Parsing here rather than at dock time is the whole point: a mis-read file
    should surface in a second, not twenty minutes into a run.

    `display_name` is the name the file had before it was stored under the
    slot letter. A PDB without a HEADER record has no entry id of its own, so
    the reader falls back to the filename stem -- which, after we have renamed
    the file to `A.pdb`, would show the user an entry called "A".
    """
    require("numpy")
    _sidecar_on_path()
    from glueforge.cif import read_any

    st = read_any(path)
    entry = st.entry_id
    if display_name and (not entry or entry == path.stem):
        entry = Path(display_name).stem[:32] or entry
    sizes = st.chain_sizes()

    # How close each chain comes to the biggest ligand. A deposited
    # structure often carries proteins that have nothing to do with the
    # site -- CRBN arrives bound to DDB1, which in 4CI2 is 34 A from the
    # glue and touches none of it -- and the only way to decide whether a
    # chain belongs in the run is to see where it sits relative to the
    # molecule. Without this the chain list is a row of letters.
    import numpy as _np
    biggest = st.largest_ligand()
    near = {}
    if biggest is not None:
        from scipy.spatial import cKDTree as _KD
        L = _np.array([(a.x, a.y, a.z) for a in biggest.atoms
                       if a.element.upper() != "H"])
        if len(L):
            tree = _KD(L)
            for cid in sizes:
                X = _np.array([(a.x, a.y, a.z) for r in st.residues
                               if r.chain == cid for a in r.atoms
                               if a.element.upper() != "H"])
                if len(X):
                    d = tree.query(X, k=1)[0]
                    near[cid] = (round(float(d.min()), 1),
                                 int((d <= 5.0).sum()))

    chains = []
    for cid, n in sorted(sizes.items(), key=lambda kv: -kv[1]):
        d, c = near.get(cid, (None, 0))
        chains.append({
            "id": cid,
            "residues": n,
            "name": st.name_of(cid),
            "sequence": st.sequence(cid)[:1200],
            "ligand_distance": d,
            "ligand_contacts": c,
        })

    ligands = []
    for lig in st.ligands:
        ligands.append({
            "comp_id": lig.comp_id,
            "chain": lig.chain,
            "seq_id": lig.seq_id,
            "heavy_atoms": lig.heavy_atom_count,
            "centroid": [round(float(v), 2) for v in lig.centroid()],
        })
    ligands.sort(key=lambda d: -d["heavy_atoms"])

    return {
        "entry": entry,
        "format": path.suffix.lower().lstrip("."),
        "chains": chains,
        "n_residues": sum(sizes.values()),
        "n_chains": len(sizes),
        "ligands": ligands,
        "waters_removed": st.waters_removed,
    }


def presented_face_summary(
    path: Path,
    ligand_comp: str,
    ligand_chain: str | None = None,
    ligand_seq: int | None = None,
    anchor_chains: list[str] | None = None,
) -> dict:
    require("numpy", "scipy")
    _sidecar_on_path()
    from glueforge.cif import read_any
    from glueforge.face import default_anchor_chains, presented_face
    from glueforge.restraints import restraints_from_face

    st = read_any(path)
    lig = None
    for cand in st.ligands:
        if cand.comp_id != ligand_comp:
            continue
        if ligand_chain and cand.chain != ligand_chain:
            continue
        if ligand_seq is not None and cand.seq_id != ligand_seq:
            continue
        lig = cand
        break
    if lig is None:
        raise ValueError("no ligand %s in this structure" % ligand_comp)

    # When a file holds the whole ternary complex, measuring the face across
    # every chain the ligand touches is the thing the method exists to avoid:
    # the partner is sitting on the very surface being measured, so the number
    # comes back tiny and means nothing. Default to the chain that buries the
    # ligand most, which is the anchor, and let the caller override.
    contacting = default_anchor_chains(st, lig)
    if anchor_chains is None and len(contacting) > 1:
        anchor_chains = [contacting[0]]

    face = presented_face(st, lig, anchor_chains=anchor_chains)
    out = face.as_dict()
    # The restraint split travels with the face, because what the upload
    # screen should show is not "26 nearby residues" but the two groups a
    # docking run treats differently.
    out["restraints"] = restraints_from_face(face).as_dict()
    out["contacting_chains"] = contacting
    return out


def assess_complex(
    path: Path,
    ligand_comp: str,
    anchor_chains: list[str],
    partner_chains: list[str],
    ligand_chain: str | None = None,
    ligand_seq: int | None = None,
) -> dict:
    """Measure the interface of an assembled complex and judge it.

    Works on anything where both proteins are already positioned: a crystal
    structure, a predicted complex, or a docking pose. It cannot position
    them for you -- that is what the docking stage is for -- and a predicted
    complex is as good as a solved one here, its accuracy just becomes the
    error budget.
    """
    require("numpy", "scipy")
    _sidecar_on_path()
    from glueforge.cif import read_any
    from glueforge.face import presented_face
    from glueforge.cavity import detect_cavity
    from glueforge.gate import assess as judge
    from glueforge.gate import cavity_metrics
    from glueforge.growth import growth_points, summarise
    from glueforge.pose import PoseScorer, coords, heavy

    st = read_any(path)
    lig = None
    for cand in st.ligands:
        if cand.comp_id != ligand_comp:
            continue
        if ligand_chain and cand.chain != ligand_chain:
            continue
        if ligand_seq is not None and cand.seq_id != ligand_seq:
            continue
        lig = cand
        break
    if lig is None:
        raise ValueError("no ligand %s in this structure" % ligand_comp)

    anchor = [a for r in st.residues if r.chain in set(anchor_chains)
              for a in r.atoms]
    partner = [a for r in st.residues if r.chain in set(partner_chains)
               for a in r.atoms]
    if not anchor:
        raise ValueError("no residues in anchor chain(s) %s"
                         % ", ".join(anchor_chains))
    if not partner:
        raise ValueError("no residues in partner chain(s) %s"
                         % ", ".join(partner_chains))

    # The face is measured against the anchor ALONE even though the partner
    # is right here, because that is the surface the ligand presents and the
    # partner covers. Measuring it in the assembled complex would report the
    # sliver the partner left uncovered, which is a different quantity.
    face = presented_face(st, lig, anchor_chains=anchor_chains)

    lig_xyz = coords(heavy(lig.atoms))
    metrics = cavity_metrics(heavy(anchor), heavy(partner), lig_xyz)
    scorer = PoseScorer(anchor, lig.atoms, face)
    score = scorer.score(partner, "as given")
    metrics["face_burial"] = score.face_burial

    verdict = judge(metrics, face=face.as_dict())

    # Growth points come from the same cavity the verdict was built on, so
    # they are computed here rather than recomputed later from a cavity that
    # might not be the same one.
    anchor_h = heavy(anchor)
    partner_h = heavy(partner)
    cav = detect_cavity(
        anchor_h + partner_h, lig_xyz.mean(axis=0),
        set(range(len(anchor_h))),
    )
    partner_res = [r for r in st.residues if r.chain in set(partner_chains)]
    points = growth_points(lig, anchor, partner, face, cav, partner_res)

    return {
        "face": face.as_dict(),
        "metrics": {k: round(float(v), 4) for k, v in metrics.items()},
        "pose": score.as_dict(),
        "verdict": verdict.as_dict(),
        "growth": {
            "points": [g.as_dict() for g in points],
            "summary": summarise(points),
        },
        "anchor_chains": list(anchor_chains),
        "partner_chains": list(partner_chains),
    }


# ---------------------------------------------------------------------------
# API
# ---------------------------------------------------------------------------


class ApiError(Exception):
    def __init__(self, status: int, message: str, **extra) -> None:
        super().__init__(message)
        self.status = status
        self.message = message
        self.extra = extra


class Api:
    """Routes, kept as plain methods so the dispatch table reads as a list."""

    # -- environment --

    def get_status(self, q: dict) -> dict:
        env = environment()
        env["runs"] = list_runs(12)
        return env

    def post_install(self, q: dict, body: bytes) -> dict:
        payload = json_body(body)
        key = payload.get("package")
        spec = next((p for p in PACKAGES if p["key"] == key), None)
        if spec is None:
            raise ApiError(400, "unknown package %r" % (key,))
        job = pip_install(spec["pip"])
        return {"job": job.id, "package": key}

    def get_job(self, q: dict, job_id: str) -> dict:
        job = JOBS.get(job_id)
        if job is None:
            raise ApiError(404, "no such job")
        return job.snapshot()

    # -- runs --

    def get_runs(self, q: dict) -> dict:
        return {"runs": list_runs()}

    def post_run(self, q: dict, body: bytes) -> dict:
        payload = json_body(body)
        state = {
            "id": new_run_id(),
            "created": time.time(),
            "stage": "upload",
            "label": (payload.get("label") or "").strip()[:80] or None,
            "slots": {},
            "notes": [],
        }
        save_run(state)
        return state

    def get_run(self, q: dict, run_id: str) -> dict:
        try:
            return load_run(run_id)
        except FileNotFoundError:
            raise ApiError(404, "no such run") from None

    def post_upload(self, q: dict, body: bytes, run_id: str) -> dict:
        """Raw bytes in, structure summary out.

        Raw rather than multipart on purpose: multipart parsing in the stdlib
        is a liability for 20 MB of mmCIF, and the page has nothing else to
        send alongside the file.
        """
        slot = (q.get("slot", ["A"])[0] or "A").upper()
        if slot not in ("A", "B"):
            raise ApiError(400, "slot must be A or B")
        name = os.path.basename(unquote(q.get("name", ["structure.pdb"])[0]))
        ext = Path(name).suffix.lower()
        if ext not in (".pdb", ".cif", ".mmcif", ".ent"):
            raise ApiError(
                400,
                "%s is not a structure file. Give a .pdb or .cif." % (name,),
            )
        if not body:
            raise ApiError(400, "empty file")

        try:
            state = load_run(run_id)
        except FileNotFoundError:
            raise ApiError(404, "no such run") from None

        d = run_dir(run_id)
        d.mkdir(parents=True, exist_ok=True)
        dest = d / (slot + ext)
        for stale in d.glob(slot + ".*"):
            if stale.name != dest.name:
                stale.unlink()
        dest.write_bytes(body)

        try:
            summary = summarise_structure(dest, display_name=name)
        except NeedsPackages as exc:
            raise ApiError(
                409, "install %s first" % ", ".join(exc.missing),
                missing=exc.missing,
            ) from None
        except Exception as exc:  # noqa: BLE001 - a bad file is a user error
            raise ApiError(
                422,
                "could not read %s: %s" % (name, exc),
            ) from None

        if summary["n_residues"] == 0:
            raise ApiError(
                422,
                "%s parsed, but holds no amino-acid residues. A nucleic acid "
                "or a ligand-only file will not dock." % name,
            )

        slot_state = {
            "name": name,
            "file": dest.name,
            "bytes": len(body),
            "uploaded": time.time(),
            **summary,
        }
        state.setdefault("slots", {})[slot] = slot_state
        save_run(state)
        return {"slot": slot, "structure": slot_state, "run": state["id"]}

    def get_examples(self, q: dict) -> dict:
        """The worked examples that ship with the tool.

        Prepared files rather than PDB entries: each is the chains that
        matter, already split into protein A with its ligand and protein B
        on its own. The numbers in each description were produced by
        `sidecar/make_examples.py` running the same code the tool runs, so
        an example cannot describe itself in terms its own measurement would
        contradict.
        """
        index = CACHE / "examples" / "index.json"
        if not index.exists():
            return {"examples": [], "note": (
                "No examples are built. Run sidecar/make_examples.py, which "
                "needs the benchmark structures in cache/raw."
            )}
        try:
            data = json.loads(index.read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            raise ApiError(500, "could not read the examples: %s" % exc) from None
        return data

    def post_example(self, q: dict, body: bytes, run_id: str) -> dict:
        """Load one example's files into a run.

        Copied in rather than referenced, so an example behaves exactly like
        an upload from here on -- the same chain selection, the same
        rewriting, the same everything. An example that took a different
        path through the tool would be demonstrating a path nobody else
        uses.
        """
        payload = json_body(body)
        want = str(payload.get("example") or "")
        if not SAFE_ID.match(want):
            raise ApiError(400, "bad example id")

        try:
            state = load_run(run_id)
        except FileNotFoundError:
            raise ApiError(404, "no such run") from None

        src = CACHE / "examples" / want
        meta_path = src / "example.json"
        if not meta_path.exists():
            raise ApiError(404, "no example %r" % want)
        meta = json.loads(meta_path.read_text(encoding="utf-8"))

        d = run_dir(run_id)
        d.mkdir(parents=True, exist_ok=True)
        slots = {}
        for slot_id in ("A", "B"):
            name = (meta.get("files") or {}).get(slot_id)
            if not name or not (src / name).exists():
                continue
            for stale in d.glob(slot_id + ".*"):
                stale.unlink()
            dest = d / (slot_id + ".pdb")
            shutil.copyfile(src / name, dest)
            try:
                summary = summarise_structure(dest, display_name="%s %s" % (
                    meta["entry"], "A" if slot_id == "A" else "B"))
            except Exception as exc:  # noqa: BLE001
                raise ApiError(
                    422, "example %s would not parse: %s" % (want, exc)
                ) from None
            slots[slot_id] = {
                "name": "%s — %s" % (meta["entry"], meta.get("title", "")),
                "file": dest.name,
                "bytes": dest.stat().st_size,
                "uploaded": time.time(),
                "from_example": want,
                **summary,
            }

        if "A" not in slots:
            raise ApiError(500, "example %s has no protein A" % want)

        state["slots"] = slots
        # An example arrives with its ligand already known, which is the one
        # thing it can carry that an upload cannot. Everything downstream is
        # then identical to a file the user chose themselves.
        state["ligand"] = {
            "comp_id": meta.get("ligand"),
            "chain": meta.get("anchor_chain"),
            "seq_id": None,
        }
        # The anchor is known for an example, and saying so means the gate
        # can be reached in one more click rather than after a decision the
        # example was chosen to have already made.
        state["anchor_chains"] = [meta["anchor_chain"]]
        for k in ("face", "assessment", "anchor_chains", "dock", "growth"):
            state.pop(k, None)
        state["example"] = {"id": want, "entry": meta.get("entry"),
                            "title": meta.get("title")}
        state["stage"] = "uploaded"
        save_run(state)
        return {"run": state["id"], "example": meta, "slots": slots}

    def post_chains(self, q: dict, body: bytes, run_id: str) -> dict:
        """Keep only some of a structure's chains.

        A deposited structure is whatever crystallised, not what the question
        is about. CRBN comes bound to DDB1; a ternary complex comes with the
        substrate already in place. Both are reasons someone cannot start:
        the extra protein either sits in the way of docking a new one, or
        makes the tool think the complex is finished.

        So the chains are a choice, and it is made by rewriting protein A to
        the chosen ones rather than by carrying a filter through every stage
        after it. `dock.py` is a separate process that reads this file by
        path, as does the exporter and the viewer; one filtered file on disk
        is correct everywhere by construction, where a flag would have to be
        threaded through each of them and would be wrong in whichever one it
        was forgotten.

        The upload is never overwritten, so the choice can be changed or
        undone.
        """
        payload = json_body(body)
        try:
            state = load_run(run_id)
        except FileNotFoundError:
            raise ApiError(404, "no such run") from None

        slot_id = (payload.get("slot") or "A").upper()
        slot = (state.get("slots") or {}).get(slot_id)
        if not slot:
            raise ApiError(400, "nothing uploaded in slot %s" % slot_id)

        keep = payload.get("keep")
        if keep is not None and not isinstance(keep, list):
            raise ApiError(400, "keep must be a list of chain ids")

        d = run_dir(run_id)
        original = slot.get("original") or slot["file"]
        src = d / original
        if not src.exists():
            raise ApiError(410, "the uploaded file is no longer on disk")

        _sidecar_on_path()
        from glueforge.cif import read_any
        from glueforge.restraints import write_pdb

        st = read_any(src)
        available = set(st.chain_sizes())

        # No selection, or all of them: go back to the file as uploaded.
        if keep is None or set(keep) >= available:
            slot["file"] = original
            slot.pop("kept_chains", None)
        else:
            unknown = [c for c in keep if c not in available]
            if unknown:
                raise ApiError(400, "no chain %s in this file" % ", ".join(unknown))
            if not keep:
                raise ApiError(400, "keep at least one chain")

            wanted = set(keep)
            atoms = [a for r in st.residues if r.chain in wanted for a in r.atoms]
            if not atoms:
                raise ApiError(400, "those chains hold no residues")
            # Ligands on a kept chain come with it. A chain kept for its
            # protein but stripped of the molecule bound to it would be a
            # silent way to lose the whole point of the run.
            for lig in st.ligands:
                if lig.chain in wanted:
                    atoms += list(lig.atoms)

            dest = d / (slot_id + "-chains.pdb")
            try:
                write_pdb(atoms, dest)
            except Exception as exc:  # noqa: BLE001
                raise ApiError(
                    422, "could not write the selected chains: %s" % exc
                ) from None
            slot["original"] = original
            slot["file"] = dest.name
            slot["kept_chains"] = sorted(wanted)
            slot["dropped_chains"] = sorted(available - wanted)

            # The chains taken out of protein A are usually the obvious
            # protein B -- that is what dropping them was for. Offered
            # rather than assumed, because "remove this from the anchor"
            # and "dock this against the anchor" are different intentions
            # and only the user knows which one this was.
            if payload.get("dropped_to_b"):
                spare = [a for r in st.residues
                         if r.chain in (available - wanted) for a in r.atoms]
                if not spare:
                    raise ApiError(400, "the dropped chains hold no residues")
                bdest = d / "B-from-A.pdb"
                write_pdb(spare, bdest)
                try:
                    bsum = summarise_structure(
                        bdest, display_name="dropped from protein A")
                except Exception as exc:  # noqa: BLE001
                    raise ApiError(
                        422, "could not re-read the dropped chains: %s" % exc
                    ) from None
                state["slots"]["B"] = {
                    "name": "chain %s, dropped from protein A" % ", ".join(
                        sorted(available - wanted)),
                    "file": bdest.name,
                    "bytes": bdest.stat().st_size,
                    "uploaded": time.time(),
                    "from_slot_a": True,
                    **bsum,
                }

        try:
            summary = summarise_structure(d / slot["file"],
                                          display_name=slot.get("name"))
        except Exception as exc:  # noqa: BLE001
            raise ApiError(422, "could not re-read the selection: %s" % exc) from None
        slot.update(summary)

        # The face, the assessment and anything docked were measured against
        # a different structure. Keeping them beside a changed selection is
        # how a number gets quoted for a complex it was never computed on.
        dropped = [k for k in ("face", "assessment", "anchor_chains", "dock",
                               "growth", "assessed_pose") if k in state]
        for k in dropped:
            state.pop(k, None)
        state["stage"] = "uploaded"
        save_run(state)
        return {"slot": slot_id, "structure": slot, "run": state["id"],
                "invalidated": dropped}

    def post_face(self, q: dict, body: bytes, run_id: str) -> dict:
        payload = json_body(body)
        try:
            state = load_run(run_id)
        except FileNotFoundError:
            raise ApiError(404, "no such run") from None

        slot = state.get("slots", {}).get("A")
        if not slot:
            raise ApiError(400, "upload protein A first")

        comp = payload.get("ligand")
        if not comp:
            raise ApiError(400, "pick a ligand")
        anchor = payload.get("anchor_chains") or None

        try:
            face = presented_face_summary(
                run_dir(run_id) / slot["file"],
                comp,
                ligand_chain=payload.get("ligand_chain"),
                ligand_seq=payload.get("ligand_seq"),
                anchor_chains=anchor,
            )
        except NeedsPackages as exc:
            raise ApiError(
                409, "install %s first" % ", ".join(exc.missing),
                missing=exc.missing,
            ) from None
        except ValueError as exc:
            raise ApiError(400, str(exc)) from None

        state["ligand"] = {
            "comp_id": comp,
            "chain": payload.get("ligand_chain"),
            "seq_id": payload.get("ligand_seq"),
        }
        state["anchor_chains"] = face["anchor_chains"]
        state["face"] = face
        state["stage"] = "face"
        save_run(state)
        return {"face": face, "run": state["id"]}

    def _assemble(self, run_id: str, state: dict, pose):
        """The run as one complex: the anchor, the ligand, and the partner
        wherever it currently is.

        Shared by the viewer and the gate, because they have to be looking
        at the same thing. A pose drawn on screen and a pose sent to be
        judged that were assembled by two different pieces of code is the
        kind of disagreement nothing on screen could show.

        Returns the residues, the ligand, and a label for what this is.
        """
        slots = state.get("slots") or {}
        if "A" not in slots:
            raise ApiError(400, "nothing uploaded yet")
        lig_info = state.get("ligand") or {}
        if not lig_info.get("comp_id"):
            raise ApiError(400, "pick a ligand first")

        from glueforge.cif import Residue, read_any

        d = run_dir(run_id)
        st = read_any(d / slots["A"]["file"])
        lig = None
        for cand in st.ligands:
            if cand.comp_id != lig_info["comp_id"]:
                continue
            if lig_info.get("chain") and cand.chain != lig_info["chain"]:
                continue
            lig = cand
            break
        if lig is None:
            raise ApiError(400, "the chosen ligand is not in protein A")

        anchor_chains = set(state.get("anchor_chains") or [])
        if not anchor_chains:
            raise ApiError(400, "no anchor chain has been chosen")
        anchor_res = [r for r in st.residues if r.chain in anchor_chains]
        if not anchor_res:
            raise ApiError(400, "the anchor chain is not in protein A")

        label, partner_res = "", []
        if pose:
            # The pose name arrives from the page, so it is checked against
            # the same pattern run ids are rather than interpolated into a
            # path and hoped for.
            if not SAFE_ID.match(str(pose)):
                raise ApiError(400, "bad pose name")
            p = d / "dock" / "poses" / ("%s.pdb" % pose)
            if not p.exists():
                raise ApiError(404, "no pose %r in this run" % pose)
            partner_res = list(read_any(p).residues)
            label = str(pose)
        else:
            partner_res = [r for r in st.residues
                           if r.chain not in anchor_chains]
            label = "as uploaded"
            if not partner_res and "B" in slots:
                partner_res = list(
                    read_any(d / slots["B"]["file"]).residues)
                label = "before docking"

        # The two came from different files and may carry the same chain
        # letter, which would merge them into one protein on screen and in
        # every interface calculation. Renamed here rather than at write
        # time, because the pose files are also opened in PyMOL, where the
        # original chain naming is what the user expects to see.
        taken = {r.chain for r in anchor_res}
        if any(r.chain in taken for r in partner_res):
            free = next(c for c in "BCDEFGHIJKLMNOPQRSTUVWXYZ0"
                        if c not in taken)
            partner_res = [
                Residue(r.comp_id, free, r.seq_id, r.insertion, r.one,
                        [type(a)(a.serial, a.name, a.element, a.comp_id, free,
                                 a.seq_id, a.x, a.y, a.z, a.occupancy,
                                 a.hetatm, a.insertion) for a in r.atoms])
                for r in partner_res
            ]

        if not partner_res:
            raise ApiError(
                400,
                "this run has no second protein in place yet. Dock one, or "
                "upload a complex that already holds both."
            )
        return anchor_res, partner_res, lig, label

    def get_view(self, q: dict, run_id: str) -> dict:
        """The run as something the page's 3D viewer can draw.

        Two shapes of run reach here. One has a solved complex in slot A and
        was assessed: both proteins are already in place and the view is of
        the file as given. The other was docked, so protein B's position is
        a result rather than an input, and the view is of one chosen pose --
        which means the partner's coordinates come from the pose file, not
        from the upload it started as.

        Built on demand rather than at dock time. A run can return hundreds
        of poses and almost all of them are never looked at; building every
        view up front would cost more than the docking did.
        """
        try:
            state = load_run(run_id)
        except FileNotFoundError:
            raise ApiError(404, "no such run") from None

        slots = state.get("slots") or {}
        if "A" not in slots:
            raise ApiError(400, "nothing uploaded yet")
        lig_info = state.get("ligand") or {}
        if not lig_info.get("comp_id"):
            raise ApiError(400, "pick a ligand first")

        pose = q.get("pose", [None])
        pose = pose[0] if isinstance(pose, list) else pose

        _sidecar_on_path()
        from glueforge.viewblob import build

        parts = self._assemble(run_id, state, pose)
        anchor_res, partner_res, lig, label = parts

        try:
            blob = build(
                anchor_res + partner_res, lig,
                anchor_chain=anchor_res[0].chain,
                partner_chain=partner_res[0].chain,
                anchor_name=(state.get("slots") or {}).get("A", {}).get("name")
                or "Protein A",
                partner_name=(state.get("slots") or {}).get("B", {}).get("name")
                or "Protein B",
                entry=run_id, label=label,
            )
        except ValueError as exc:
            raise ApiError(400, str(exc)) from None
        return blob

    def get_presets(self, q: dict) -> dict:
        """The three speed/accuracy points, with whatever is known about them.

        `recovery` comes back null until run_recovery.py has been run on a
        machine with LightDock. A preset showing a made-up recovery rate
        beside a real run time is worse than one that admits it does not
        know, because the made-up one gets quoted.
        """
        _sidecar_on_path()
        from glueforge.restraints import PRESETS, DEFAULT_PRESET, missing_tools
        measured = {}
        rec = CACHE / "recovery" / "recovery.json"
        if rec.exists():
            try:
                data = json.loads(rec.read_text(encoding="utf-8"))
                for key, block in (data.get("presets") or {}).items():
                    good = [r for r in block.get("systems", [])
                            if "error" not in r]
                    if good:
                        measured[key] = {
                            "recovered": sum(1 for r in good if r["recovered"]),
                            "designable": sum(1 for r in good if r["designable"]),
                            "n": len(good),
                            "cutoff": data.get("cutoff"),
                            "machine": (data.get("machine") or {}).get("platform"),
                        }
            except (OSError, ValueError):
                measured = {}

        # The built-in engine has its own presets and its own measurements,
        # from run_builtin_recovery.py. They are reported under their own key
        # rather than merged: the two engines search differently and a
        # recovery rate measured on one says nothing about the other.
        from glueforge.fftdock import PRESETS as FFT_PRESETS
        builtin_measured = None
        bfile = CACHE / "recovery" / "builtin.json"
        if bfile.exists():
            try:
                data = json.loads(bfile.read_text(encoding="utf-8"))
                good = [r for r in data.get("systems", []) if "error" not in r]
                if good:
                    builtin_measured = {
                        "recovered": sum(1 for r in good if r["recovered"]),
                        "designable": sum(1 for r in good if r["designable"]),
                        "n": len(good),
                        "cutoff": data.get("cutoff"),
                        "preset": (data.get("settings") or {}).get("key"),
                        "bound_bound": bool(data.get("bound_bound")),
                        "largest_partner": max(
                            (r.get("partner_atoms") or 0) for r in good),
                    }
            except (OSError, ValueError):
                builtin_measured = None

        engines = {
            "lightdock": {
                "available": not missing_tools(),
                "presets": [{**PRESETS[k].as_dict(),
                             "measured": measured.get(k)}
                            for k in ("quick", "standard", "thorough")],
            },
            "builtin": {
                "available": True,
                "presets": [FFT_PRESETS[k].as_dict()
                            for k in ("quick", "standard", "thorough")],
                "measured": builtin_measured,
            },
        }

        out = []
        for key in ("quick", "standard", "thorough"):
            d = PRESETS[key].as_dict()
            d["measured"] = measured.get(key)
            out.append(d)
        return {"presets": out, "default": DEFAULT_PRESET,
                "calibrated": bool(measured),
                "engines": engines,
                "engine_default": "lightdock" if not missing_tools()
                                  else "builtin"}

    def post_dock(self, q: dict, body: bytes, run_id: str) -> dict:
        payload = json_body(body)
        try:
            state = load_run(run_id)
        except FileNotFoundError:
            raise ApiError(404, "no such run") from None

        slots = state.get("slots") or {}
        if "A" not in slots or "B" not in slots:
            raise ApiError(400, "upload both proteins first")
        if not state.get("face"):
            raise ApiError(400, "pick a ligand and measure the face first")

        _sidecar_on_path()
        from glueforge.restraints import PRESETS, missing_tools

        # LightDock when it is there, the built-in engine when it is not.
        # Docking used to refuse outright here, which on Windows meant it
        # refused permanently -- LightDock publishes no build for it.
        engine = payload.get("engine") or "auto"
        if engine not in ("auto", "lightdock", "builtin"):
            raise ApiError(400, "unknown engine %r" % engine)
        missing = missing_tools()
        if engine == "auto":
            engine = "builtin" if missing else "lightdock"
        elif engine == "lightdock" and missing:
            raise ApiError(
                409,
                "LightDock is not on PATH (%s). The built-in engine can run "
                "this instead." % ", ".join(missing),
                missing=["lightdock"],
            )

        preset = payload.get("preset") or "standard"
        if engine == "builtin":
            from glueforge.fftdock import PRESETS as ENGINE_PRESETS
        else:
            ENGINE_PRESETS = PRESETS
        if preset not in ENGINE_PRESETS:
            raise ApiError(400, "unknown preset %r" % preset)

        d = run_dir(run_id)
        lig = (state.get("ligand") or {}).get("comp_id")
        argv = [
            sys.executable, str(SIDECAR / "dock.py"),
            str(d / slots["A"]["file"]), str(d / slots["B"]["file"]),
            "--preset", preset,
            "--engine", engine,
            "--out", str(d / "dock"),
        ]
        if lig:
            argv += ["--ligand", lig]
        for c in (state.get("anchor_chains") or []):
            argv += ["--anchor", c]

        job = JOBS.start("dock " + run_id, argv)
        state["dock"] = {"job": job.id, "preset": preset, "engine": engine,
                         "started": time.time()}
        state["stage"] = "docking"
        save_run(state)
        return {"job": job.id, "preset": preset, "run": run_id}

    def get_dock(self, q: dict, run_id: str) -> dict:
        """What the docking run produced, once it has."""
        try:
            state = load_run(run_id)
        except FileNotFoundError:
            raise ApiError(404, "no such run") from None

        info = dict(state.get("dock") or {})
        job = JOBS.get(info.get("job", ""))
        info["state"] = job.state if job else "unknown"
        if job:
            info["elapsed"] = job.snapshot()["elapsed"]

        scored = run_dir(run_id) / "dock" / "scored.json"
        if scored.exists():
            try:
                info["result"] = json.loads(scored.read_text(encoding="utf-8"))
            except (OSError, ValueError) as exc:
                info["error"] = "could not read scored.json: %s" % exc
        return info

    def post_assess(self, q: dict, body: bytes, run_id: str) -> dict:
        payload = json_body(body)
        try:
            state = load_run(run_id)
        except FileNotFoundError:
            raise ApiError(404, "no such run") from None

        slot_name = (payload.get("slot") or "A").upper()
        slot = state.get("slots", {}).get(slot_name)
        if not slot:
            raise ApiError(400, "upload protein %s first" % slot_name)

        lig = payload.get("ligand") or (state.get("ligand") or {}).get("comp_id")
        if not lig:
            raise ApiError(400, "pick a ligand first")

        # Judging a docked pose. The partner is not in any uploaded file --
        # its position is the result -- so the complex is assembled here
        # from the anchor and that pose, by the same code the viewer uses,
        # and written out for the gate to read.
        #
        # Without this, docking and the gate were two halves that never met:
        # the dock screen offered "Assess this pose" and the assess screen
        # then said nothing had been assessed, because nothing had.
        pose = payload.get("pose")
        source = run_dir(run_id) / slot["file"]
        anchor = payload.get("anchor_chains") or state.get("anchor_chains")
        partner = payload.get("partner_chains")

        if pose:
            _sidecar_on_path()
            from glueforge.restraints import write_pdb

            anchor_res, partner_res, lig_obj, label = self._assemble(
                run_id, state, pose)
            atoms = [a for r in anchor_res for a in r.atoms]
            atoms += list(lig_obj.atoms)
            atoms += [a for r in partner_res for a in r.atoms]
            source = run_dir(run_id) / "dock" / ("assessed-%s.pdb" % pose)
            source.parent.mkdir(parents=True, exist_ok=True)
            write_pdb(atoms, source)
            anchor = [anchor_res[0].chain]
            partner = [partner_res[0].chain]

        if not anchor or not partner:
            raise ApiError(
                400,
                "say which chains are the anchor and which the partner",
            )

        try:
            out = assess_complex(
                source, lig,
                list(anchor), list(partner),
                ligand_chain=payload.get("ligand_chain")
                or (state.get("ligand") or {}).get("chain"),
                ligand_seq=payload.get("ligand_seq")
                if payload.get("ligand_seq") is not None
                else (state.get("ligand") or {}).get("seq_id"),
            )
        except NeedsPackages as exc:
            raise ApiError(409, "install %s first" % ", ".join(exc.missing),
                           missing=exc.missing) from None
        except ValueError as exc:
            raise ApiError(400, str(exc)) from None

        state["assessment"] = out
        if pose:
            # Which pose this verdict is about. An assessment with no note
            # of what it judged is one that will be read as applying to
            # whatever is on screen later.
            state["assessment"]["pose"] = pose
            state["assessed_pose"] = pose
        state["stage"] = "assessed"
        save_run(state)
        return out

    def export_run(self, run_id: str) -> tuple[bytes, str]:
        """A finished assessment as one file that opens anywhere.

        Built from the standalone page, so the calibration travels with it
        and every threshold in the exported copy is still one click from the
        evidence that set it. A shared result that cannot show its working is
        a screenshot with extra steps.
        """
        try:
            state = load_run(run_id)
        except FileNotFoundError:
            raise ApiError(404, "no such run") from None
        if not state.get("assessment"):
            raise ApiError(400, "this run has no assessment to export")

        page = None
        for p in STANDALONE_CANDIDATES:
            if p.exists():
                page = p.read_text(encoding="utf-8")
                break
        if page is None:
            raise ApiError(
                500,
                "The standalone page has not been built. Run "
                "`node build.js --standalone`.",
            )

        cal = None
        f = CACHE / "gate" / "gate.json"
        if f.exists():
            try:
                cal = json.loads(f.read_text(encoding="utf-8"))
            except (OSError, ValueError):
                cal = None

        payload = json.dumps({
            "run": {
                k: state.get(k) for k in
                ("id", "label", "created", "updated", "stage", "slots",
                 "ligand", "anchor_chains", "face", "assessment", "dock")
            },
            "calibration": cal,
            "exported": time.time(),
            "version": VERSION,
        }, allow_nan=False).replace("</script", "<\\/script")

        inject = (
            "<script>globalThis.GF = globalThis.GF || {};"
            "GF.EXPORT = " + payload + ";</script>\n</body>"
        )
        if "</body>" not in page:
            raise ApiError(500, "the standalone page has no </body> to patch")
        out = page.replace("</body>", inject, 1)

        name = "glueforge-%s.html" % (state.get("id") or run_id)
        return out.encode("utf-8"), name

    def get_calibration(self, q: dict) -> dict:
        """The evidence behind every gate threshold.

        Served whole -- it is 39 kB -- because a threshold you cannot see the
        working for is a threshold nobody can argue with, and the arguable
        ones are the only kind worth having.
        """
        f = CACHE / "gate" / "gate.json"
        if not f.exists():
            raise ApiError(
                404,
                "No gate calibration. Run sidecar/run_gate.py to build it.",
            )
        return json.loads(f.read_text(encoding="utf-8"))

    # -- validation --

    def get_validation(self, q: dict) -> dict:
        """The 14-system cache, assembled the way build.js assembles it."""
        index_path = WEB_CACHE / "index.json"
        if not index_path.exists():
            raise ApiError(
                404,
                "No validation cache. Run sidecar/export_web.py to build it.",
            )
        index = json.loads(index_path.read_text(encoding="utf-8"))
        systems = {}
        for s in index.get("systems", []):
            f = WEB_CACHE / (s["entry"] + ".json")
            if f.exists():
                systems[s["entry"]] = json.loads(f.read_text(encoding="utf-8"))
        potential = None
        cv = WEB_CACHE / "crossval.json"
        if cv.exists():
            potential = json.loads(cv.read_text(encoding="utf-8")).get("potential")
        return {"index": index, "systems": systems, "potential": potential}


API = Api()


def json_body(body: bytes) -> dict:
    if not body:
        return {}
    try:
        out = json.loads(body.decode("utf-8"))
    except (UnicodeDecodeError, ValueError):
        raise ApiError(400, "body is not JSON") from None
    if not isinstance(out, dict):
        raise ApiError(400, "body must be a JSON object")
    return out


# ---------------------------------------------------------------------------
# HTTP
# ---------------------------------------------------------------------------


def page_bytes() -> bytes | None:
    for p in PAGE_CANDIDATES:
        if p.exists():
            return p.read_bytes()
    return None


class Handler(http.server.BaseHTTPRequestHandler):
    server_version = "GlueForge/" + VERSION
    protocol_version = "HTTP/1.1"

    # -- plumbing --

    def log_message(self, fmt: str, *args) -> None:
        # One line per request is noise when the page polls; errors already
        # print their own traceback.
        pass

    def _send(self, status: int, body: bytes, ctype: str,
              extra: dict | None = None) -> None:
        self.send_response(status)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        # The page is same-origin; this just stops a stray extension poking it.
        self.send_header("X-Content-Type-Options", "nosniff")
        for k, v in (extra or {}).items():
            self.send_header(k, v)
        self.end_headers()
        if self.command != "HEAD":
            self.wfile.write(body)

    def _json(self, obj, status: int = 200) -> None:
        self._send(status, json.dumps(obj).encode("utf-8"),
                   "application/json; charset=utf-8")

    def _error(self, status: int, message: str, **extra) -> None:
        self._json({"error": message, **extra}, status)

    def _read_body(self) -> bytes:
        length = int(self.headers.get("Content-Length") or 0)
        if length > MAX_UPLOAD:
            raise ApiError(413, "file is larger than %d MB"
                           % (MAX_UPLOAD // (1024 * 1024)))
        if length <= 0:
            return b""
        chunks, left = [], length
        while left > 0:
            chunk = self.rfile.read(min(left, 1 << 20))
            if not chunk:
                break
            chunks.append(chunk)
            left -= len(chunk)
        return b"".join(chunks)

    # -- verbs --

    def do_HEAD(self) -> None:  # noqa: N802
        self.do_GET()

    def do_GET(self) -> None:  # noqa: N802
        url = urlparse(self.path)
        path = url.path
        q = parse_qs(url.query)

        try:
            if path.startswith("/api/"):
                return self._api_get(path, q)
            if path in ("/", "/index.html"):
                return self._serve_page()
            if path == "/favicon.ico":
                return self._send(204, b"", "image/x-icon")
            if path.startswith("/runfile/"):
                return self._serve_runfile(path)
            return self._error(404, "no route for " + path)
        except ApiError as exc:
            return self._error(exc.status, exc.message, **exc.extra)
        except BadRequest as exc:
            return self._error(400, str(exc))
        except BrokenPipeError:
            return
        except Exception as exc:  # noqa: BLE001 - a handler must not kill the server
            traceback.print_exc()
            return self._error(500, "%s: %s" % (type(exc).__name__, exc))

    def do_POST(self) -> None:  # noqa: N802
        url = urlparse(self.path)
        path, q = url.path, parse_qs(url.query)
        try:
            body = self._read_body()
            if path == "/api/install":
                return self._json(API.post_install(q, body))
            if path == "/api/run":
                return self._json(API.post_run(q, body))
            m = re.fullmatch(r"/api/run/([^/]+)/upload", path)
            if m:
                return self._json(API.post_upload(q, body, m.group(1)))
            m = re.fullmatch(r"/api/run/([^/]+)/example", path)
            if m:
                return self._json(API.post_example(q, body, m.group(1)))
            m = re.fullmatch(r"/api/run/([^/]+)/chains", path)
            if m:
                return self._json(API.post_chains(q, body, m.group(1)))
            m = re.fullmatch(r"/api/run/([^/]+)/face", path)
            if m:
                return self._json(API.post_face(q, body, m.group(1)))
            m = re.fullmatch(r"/api/run/([^/]+)/assess", path)
            if m:
                return self._json(API.post_assess(q, body, m.group(1)))
            m = re.fullmatch(r"/api/run/([^/]+)/dock", path)
            if m:
                return self._json(API.post_dock(q, body, m.group(1)))
            if path == "/api/shutdown":
                threading.Timer(0.3, self.server.shutdown).start()
                return self._json({"stopping": True})
            return self._error(404, "no route for " + path)
        except ApiError as exc:
            return self._error(exc.status, exc.message, **exc.extra)
        except BadRequest as exc:
            return self._error(400, str(exc))
        except BrokenPipeError:
            return
        except Exception as exc:  # noqa: BLE001
            traceback.print_exc()
            return self._error(500, "%s: %s" % (type(exc).__name__, exc))

    # -- GET routes --

    def _api_get(self, path: str, q: dict) -> None:
        if path == "/api/status":
            return self._json(API.get_status(q))
        if path == "/api/runs":
            return self._json(API.get_runs(q))
        if path == "/api/validation":
            return self._json(API.get_validation(q))
        if path == "/api/examples":
            return self._json(API.get_examples(q))
        if path == "/api/presets":
            return self._json(API.get_presets(q))
        if path == "/api/calibration":
            return self._json(API.get_calibration(q))
        m = re.fullmatch(r"/api/run/([^/]+)/dock", path)
        if m:
            return self._json(API.get_dock(q, m.group(1)))
        m = re.fullmatch(r"/api/run/([^/]+)/view", path)
        if m:
            return self._json(API.get_view(q, m.group(1)))
        m = re.fullmatch(r"/api/run/([^/]+)/export", path)
        if m:
            body, name = API.export_run(m.group(1))
            return self._send(
                200, body, "text/html; charset=utf-8",
                {"Content-Disposition": 'attachment; filename="%s"' % name},
            )
        m = re.fullmatch(r"/api/job/([^/]+)/stream", path)
        if m:
            return self._stream_job(m.group(1))
        m = re.fullmatch(r"/api/job/([^/]+)", path)
        if m:
            return self._json(API.get_job(q, m.group(1)))
        m = re.fullmatch(r"/api/run/([^/]+)", path)
        if m:
            return self._json(API.get_run(q, m.group(1)))
        return self._error(404, "no route for " + path)

    def _serve_page(self) -> None:
        body = page_bytes()
        if body is None:
            msg = (
                "<!doctype html><meta charset=utf-8>"
                "<title>GlueForge</title>"
                "<body style='font:16px/1.6 system-ui;max-width:40rem;"
                "margin:4rem auto;padding:0 1rem'>"
                "<h1>The page has not been built</h1>"
                "<p>The server is running, but neither "
                "<code>dist/glueforge.html</code> nor "
                "<code>glueforge.html</code> exists beside "
                "<code>glueforge.py</code>.</p>"
                "<p>Build it with <code>node build.js</code>, or use a "
                "release zip, which ships the page already built.</p>"
            ).encode("utf-8")
            return self._send(500, msg, "text/html; charset=utf-8")
        return self._send(200, body, "text/html; charset=utf-8")

    def _serve_runfile(self, path: str) -> None:
        """Hand back a file from inside a run folder, and nothing else."""
        m = re.fullmatch(r"/runfile/([^/]+)/([^/]+)", path)
        if not m:
            raise ApiError(404, "no such file")
        run_id, name = m.group(1), unquote(m.group(2))
        if not SAFE_ID.match(run_id) or "/" in name or "\\" in name or ".." in name:
            raise ApiError(400, "bad path")
        base = run_dir(run_id).resolve()
        target = (base / name).resolve()
        if base not in target.parents or not target.is_file():
            raise ApiError(404, "no such file")
        ctype = mimetypes.guess_type(name)[0] or "application/octet-stream"
        return self._send(200, target.read_bytes(), ctype)

    def _stream_job(self, job_id: str) -> None:
        """Server-sent events, so pip's output lands in the page as it happens."""
        job = JOBS.get(job_id)
        if job is None:
            raise ApiError(404, "no such job")

        self.send_response(200)
        self.send_header("Content-Type", "text/event-stream; charset=utf-8")
        self.send_header("Cache-Control", "no-store")
        self.send_header("Connection", "close")
        self.end_headers()
        self.close_connection = True

        q = job.watch()
        try:
            while True:
                try:
                    line = q.get(timeout=20)
                except queue.Empty:
                    self.wfile.write(b": keepalive\n\n")
                    self.wfile.flush()
                    continue
                if line is None:
                    payload = json.dumps(
                        {"done": True, "code": job.code, "state": job.state}
                    )
                    self.wfile.write(
                        ("event: done\ndata: %s\n\n" % payload).encode("utf-8")
                    )
                    self.wfile.flush()
                    return
                payload = json.dumps({"line": line})
                self.wfile.write(("data: %s\n\n" % payload).encode("utf-8"))
                self.wfile.flush()
        except (BrokenPipeError, ConnectionResetError, OSError):
            return
        finally:
            job.unwatch(q)


class Server(socketserver.ThreadingMixIn, http.server.HTTPServer):
    daemon_threads = True
    allow_reuse_address = True


# ---------------------------------------------------------------------------
# Startup
# ---------------------------------------------------------------------------


def pick_port(preferred: int) -> int:
    for port in range(preferred, preferred + 40):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            try:
                s.bind(("127.0.0.1", port))
                return port
            except OSError:
                continue
    raise SystemExit(
        "Every port from %d to %d is taken. Pass --port=NNNN to choose one."
        % (preferred, preferred + 39)
    )


def ensure_page() -> bool:
    """Build dist/glueforge.html if it is missing and Node can do it.

    `dist/` is generated, so it is not in the repository, so a clone has no
    page. The server would start anyway and the browser would open on "the
    page has not been built" -- an error screen reached by a launcher that
    had just said everything was fine, which is a bad way to find out that
    one command was missing.

    Building takes about a second and writes nothing outside `dist/`, so
    when Node is here it is done rather than described. When it is not,
    the terminal says so before the browser opens, because that message
    arrives while the person is still looking at the window they started.
    """
    if page_bytes() is not None:
        return True
    builder = HERE / "build.js"
    if not builder.exists():
        return False

    node = shutil.which("node")
    if node:
        print("  No built page yet — building it with Node.", flush=True)
        try:
            r = subprocess.run([node, str(builder)], cwd=str(HERE),
                               capture_output=True, text=True, timeout=120)
        except (OSError, subprocess.SubprocessError) as exc:
            print("  Could not run Node: %s" % exc, flush=True)
        else:
            if page_bytes() is not None:
                print("  Built dist/glueforge.html.\n", flush=True)
                return True
            detail = (r.stderr or r.stdout or "").strip().splitlines()
            print("  `node build.js` did not produce a page.%s"
                  % (("\n  " + detail[-1]) if detail else ""), flush=True)

    print(
        "\n  There is no page to serve yet.\n\n"
        "  This is the source repository, where dist/ is generated rather\n"
        "  than committed. Build it once:\n\n"
        "      node build.js\n\n"
        "  Node is needed only for that one step. If you do not have it,\n"
        "  https://nodejs.org installs it, or use a release zip, which\n"
        "  ships the page already built and needs only Python.\n",
        flush=True,
    )
    return False


def check_layout() -> None:
    """Catch the two ways this fails before the browser is even open."""
    looks_temporary = any(
        part.lower() in ("temp", "tmp") for part in HERE.parts
    ) and HERE.name.lower().startswith(("glueforge",))
    if looks_temporary:
        print(
            "\n  This looks like it is running from inside the zip.\n"
            "  Windows lets you open a zip and double-click what is inside,\n"
            "  but the files are then in a temporary folder that gets wiped.\n\n"
            "  Extract the folder somewhere real first, then start it again.\n",
            flush=True,
        )
    if not SIDECAR.is_dir():
        raise SystemExit(
            "sidecar/ is missing from %s.\nThis file needs the rest of the "
            "folder beside it." % HERE
        )


KNOWN_FLAGS = ("--selftest", "--no-browser", "--port")


def arg_value(name: str, default: str | None = None) -> str | None:
    """A flag's value, written either way.

    `--port=8777` was the only form that worked, and `--port 8777` was
    accepted in silence and ignored, so the server came up on a different
    port than the one asked for and said nothing about it. Both forms now
    work, and anything unrecognised is said out loud rather than dropped --
    a typed flag that has no effect is worse than one that is refused.
    """
    prefix = "--%s=" % name
    args = sys.argv[1:]
    for i, a in enumerate(args):
        if a.startswith(prefix):
            return a[len(prefix):]
        if a == "--%s" % name and i + 1 < len(args):
            nxt = args[i + 1]
            if not nxt.startswith("-"):
                return nxt
    return default


def warn_unknown_flags() -> list[str]:
    """Names anything on the command line that nothing reads."""
    unknown = []
    skip = False
    for a in sys.argv[1:]:
        if skip:
            skip = False
            continue
        if not a.startswith("--"):
            continue
        name = a.split("=", 1)[0]
        if name not in KNOWN_FLAGS:
            unknown.append(a)
        elif "=" not in a and name == "--port":
            skip = True
    if unknown:
        print("  Ignoring %s — nothing reads %s. Known flags: %s"
              % (", ".join(unknown),
                 "it" if len(unknown) == 1 else "them",
                 " ".join(KNOWN_FLAGS)))
    return unknown


def selftest() -> int:
    """Run every check there is and say plainly what works.

    One command because `glueforge.py` is the only file anyone is told about,
    and "is this thing installed properly" is the question people have when
    something has gone wrong and they are not in a mood to hunt for test
    runners.

    Suites that need something absent are skipped rather than failed. A check
    that could not run has not failed, and a fresh install reporting red is
    how people conclude a working tool is broken.
    """
    print("\n  GlueForge %s self-test" % VERSION)
    print("  " + "-" * 56)

    env = environment()
    # Which stages a missing package actually blocks, read from the stage
    # table rather than from the package's own `stage` label. LightDock
    # carries stage "docking" but no longer blocks it -- docking falls back
    # to the built-in engine -- and printing "needed for docking" next to it
    # would be the line that sends someone off to install it.
    blocks = {}
    improves = {}
    for st in env["stages"]:
        for k in st.get("missing", ()):
            blocks.setdefault(k, []).append(st["label"].lower())
        for k in st.get("better_with_missing", ()):
            improves.setdefault(k, []).append(st["label"].lower())

    for pkg in env["packages"]:
        if pkg["version"]:
            note = "  " + pkg["version"]
        elif pkg["required"]:
            note = "   (required)"
        elif blocks.get(pkg["key"]):
            note = "   (needed to %s)" % blocks[pkg["key"]][0]
        elif improves.get(pkg["key"]):
            note = ("   (optional; an alternative for \"%s\", which already "
                    "works without it)" % improves[pkg["key"]][0])
        else:
            note = "   (optional)"
        print("    %-10s %s%s" % (
            pkg["key"], "yes" if pkg["present"] else "no ", note))
    structures = len(list((CACHE / "raw").glob("*.cif"))) if (CACHE / "raw").is_dir() else 0
    print("    %-10s %d" % ("structures", structures))
    if not structures:
        print("               (tests needing them will skip; "
              "sidecar/fetch.py downloads them)")
    print("    %-10s %s" % ("page built", "yes" if page_bytes() else "NO"))
    print()

    have_core = all(have(m) for m in ("numpy", "scipy"))
    suites: list[tuple[str, list[str]]] = []
    if have_core:
        # Found rather than listed. A hardcoded list silently stops running
        # whatever was added after it was written, and a check that quietly
        # skips a suite is worse than no check -- it reports green for work
        # it never looked at.
        for f in sorted((SIDECAR / "tests").glob("test_*.py")):
            suites.append((f.stem, [sys.executable, str(f)]))
    else:
        print("    numpy and scipy are missing, so the engine checks cannot")
        print("    run. Start the tool and install them from the page.\n")

    if node_available():
        for f in sorted((HERE / "tests").glob("*.test.js")):
            suites.append((f.name.replace(".test.js", ""), ["node", str(f)]))
    else:
        print("    node is not installed, so the browser-side checks are")
        print("    skipped. Nothing needs node except rebuilding the page.\n")

    if not suites:
        print("  Nothing could be run.\n")
        return 1

    failed = []
    for name, argv in suites:
        print("    %-16s " % name, end="", flush=True)
        try:
            r = subprocess.run(argv, cwd=str(HERE), capture_output=True,
                               text=True, timeout=3600)
        except (OSError, subprocess.SubprocessError) as exc:
            print("could not run: %s" % exc)
            failed.append(name)
            continue
        tail = [l for l in r.stdout.strip().split("\n") if l.strip()]
        last = tail[-1] if tail else "(no output)"
        if r.returncode == 0:
            print(last)
        else:
            print("FAILED \u2014 " + last)
            failed.append(name)
            for line in tail:
                if "[FAIL]" in line or "[ERR " in line:
                    print("        " + line.strip())

    print()
    if failed:
        print("  %d of %d suites failed: %s\n"
              % (len(failed), len(suites), ", ".join(failed)))
        return 1
    print("  all %d suites passed\n" % len(suites))
    return 0


def node_available() -> bool:
    try:
        subprocess.run(["node", "--version"], capture_output=True, timeout=10)
        return True
    except (OSError, subprocess.SubprocessError):
        return False


def main() -> int:
    check_layout()
    if "--selftest" in sys.argv:
        return selftest()
    RUNS.mkdir(exist_ok=True)
    warn_unknown_flags()

    # Before the port, the server and above all the browser: opening a tab
    # on an error page is a worse way to say this than the terminal the
    # person is already looking at.
    have_page = ensure_page()

    try:
        port = int(arg_value("port", str(DEFAULT_PORT)))
    except (TypeError, ValueError):
        print("  --port needs a number. Using %d." % DEFAULT_PORT)
        port = DEFAULT_PORT
    port = pick_port(port)
    url = "http://127.0.0.1:%d/" % port

    httpd = Server(("127.0.0.1", port), Handler)

    env = environment()
    missing_required = [
        p["key"] for p in env["packages"] if p["required"] and not p["present"]
    ]

    print("\n  GlueForge %s" % VERSION)
    print("  %s" % url)
    if not have_page:
        print("  Serving nothing until the page is built — see above.")
    if missing_required:
        print("  %s not installed yet — the page has a button for that."
              % ", ".join(missing_required))
    print("  Ctrl-C here to stop.\n", flush=True)

    # Not opened when there is no page: a tab on an error screen buries the
    # one line in the terminal that says what to do about it.
    if "--no-browser" not in sys.argv and have_page:
        threading.Timer(0.4, lambda: webbrowser.open(url)).start()

    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\n  Stopped.")
    finally:
        httpd.server_close()
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        sys.exit(0)
