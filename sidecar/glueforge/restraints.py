"""Turning the presented face into a LightDock run.

LightDock takes residue restraints in a plain text file and uses them twice:
it discards the swarms that start nowhere near the restrained region, and it
orients the ligand's starting glowworms toward it. On the published 4G6M
example, restraints cut 400 starting swarms to 84. That reduction *is* the
contribution of this tool -- the presented face is a restraint nobody else
has, because it comes from the ligand rather than from the proteins.

Two decisions are worth stating because they are not LightDock's defaults.

**The ligand is not in the receptor file.** LightDock scores protein against
protein; its `fastdfire` function has no parameters for an arbitrary small
molecule, and `--noh` exists precisely because unparametrised atoms break it.
So the receptor handed to LightDock is protein only, and the consequence has
to be handled on our side: nothing in the docking run stops a pose from
settling into the ligand's own vacated pocket, where it would bury the
presented face perfectly by occupying it and rank first. Every returned pose
is therefore re-scored here against the ligand, and ones that overlap it are
discarded. See `pose.PoseScore.viable`.

**Active and passive are split by distance, not by count.** A rim residue
within van der Waals reach of an exposed ligand atom is certainly lining the
opening; one 7 angstroms away might merely be nearby. Making all of them
active over-constrains the run; making all of them passive barely constrains
it at all. The split is at `ACTIVE_WITHIN`.

Nothing in this module imports LightDock. It writes files and builds command
lines, which means the whole of it can be tested on a machine that does not
have LightDock installed -- which is the machine it was written on.
"""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

from .cif import Atom, Structure
from .face import PresentedFace
from .residues import MODIFIED_AS_POLYMER, NON_LIGAND, one_letter

#: How many of the rim become ACTIVE restraints, taken in order of how far
#: they sit along the face's outward normal.
#:
#: Small on purpose. Only 30-57% of rim residues are touched by the real
#: partner in the benchmark complexes, so a large active set is one the
#: correct answer cannot satisfy -- and LightDock's own filter would then
#: discard it. Six is about where the native pose still clears a sensible
#: fraction while the restraint still prunes swarms.
MAX_ACTIVE = 6

#: Kept for callers that still pass it; selection is by normal projection
#: now, because distance to the ligand turned out to carry no signal about
#: which residues a partner actually reaches.
ACTIVE_WITHIN = 5.0

#: What share of the active restraints a pose must satisfy to be kept.
#:
#: LightDock's tutorials use 0.4. On this benchmark the *crystal* pose
#: satisfies only 30-57% of its own rim, so 0.4 sits right at the median and
#: would throw away about half the correct answers. The real discriminator
#: here is face burial, which is re-scored afterwards; the restraint filter
#: only has to avoid being actively harmful.
FNAT = 0.25

#: LightDock's own executables, in the order a run uses them.
TOOLS = (
    "lightdock3_setup.py",
    "lightdock3.py",
    "lgd_generate_conformations.py",
    "lgd_cluster_bsas.py",
    "lgd_rank.py",
)


# ---------------------------------------------------------------------------
# Presets
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Preset:
    """One point on the speed/accuracy trade.

    `recovery` is left None until it has been measured on the 14 benchmark
    systems. It is not guessed: a preset that shows a made-up recovery rate
    next to a real run time is worse than one that admits it does not know.
    """

    key: str
    label: str
    swarms: int
    glowworms: int
    steps: int
    note: str
    recovery: float | None = None

    @property
    def minutes_estimate(self) -> tuple[int, int]:
        """A rough span, from LightDock's own figure of 1-2 minutes per 10
        steps per swarm, assuming 8 cores and that restraints cut the swarms
        to roughly a fifth."""
        swarms_after = max(1, self.swarms // 5)
        per_swarm = self.steps / 10.0
        lo = swarms_after * per_swarm * 1.0 / 8.0
        hi = swarms_after * per_swarm * 2.0 / 8.0
        return (max(1, int(lo)), max(2, int(hi)))

    def as_dict(self) -> dict:
        lo, hi = self.minutes_estimate
        return {
            "key": self.key,
            "label": self.label,
            "swarms": self.swarms,
            "glowworms": self.glowworms,
            "steps": self.steps,
            "note": self.note,
            "minutes_low": lo,
            "minutes_high": hi,
            "recovery": self.recovery,
        }


PRESETS: dict[str, Preset] = {
    "quick": Preset(
        key="quick", label="Quick", swarms=100, glowworms=40, steps=25,
        note=(
            "For triaging a list of pairs to find which deserve a real run. "
            "Treat a negative from this as 'not yet looked at properly'."
        ),
    ),
    "standard": Preset(
        key="standard", label="Standard", swarms=400, glowworms=100, steps=60,
        note="The default, and what the gate thresholds are calibrated on.",
    ),
    "thorough": Preset(
        key="thorough", label="Thorough", swarms=400, glowworms=200, steps=100,
        note=(
            "LightDock's own published settings. For a pair you are going to "
            "act on."
        ),
    ),
}

DEFAULT_PRESET = "standard"


# ---------------------------------------------------------------------------
# Restraint files
# ---------------------------------------------------------------------------


@dataclass
class RestraintSet:
    active: list[str] = field(default_factory=list)
    passive: list[str] = field(default_factory=list)

    def __len__(self) -> int:
        return len(self.active) + len(self.passive)

    def text(self) -> str:
        """The file LightDock reads.

        One residue per line, `R chain.RESNAME.number`, with a trailing `P`
        for passive. `R` because all of ours are on the receptor: we know
        where on protein A the partner has to land and nothing about where on
        protein B it binds, which is the whole asymmetry of the method.
        """
        lines = [f"R {r}" for r in self.active]
        lines += [f"R {r} P" for r in self.passive]
        return "\n".join(lines) + "\n"

    def as_dict(self) -> dict:
        return {
            "active": list(self.active),
            "passive": list(self.passive),
            "n_active": len(self.active),
            "n_passive": len(self.passive),
        }


def restraints_from_face(
    face: PresentedFace,
    active_within: float = ACTIVE_WITHIN,
    max_active: int = MAX_ACTIVE,
) -> RestraintSet:
    """Split the rim into what a partner must touch and what it may.

    Active restraints are the residues furthest along the face's outward
    normal -- the approach side. Ranking them by distance to the ligand
    instead, which is the obvious thing to do and what this did at first,
    selects residues the real partner never touches: across 324 rim residues
    in the benchmark, distance predicts contact at r = -0.07 (p = 0.22) and
    normal projection at r = +0.57 (p < 1e-4).

    Everything else is passive. That is the honest model of what is known:
    the partner has to cover this face, not touch any particular residue
    around it.
    """
    rim = list(face.rim)
    if not rim:
        return RestraintSet()

    # A face that wraps the ligand has no approach direction, so there is
    # nothing to rank by and no residue that can honestly be required.
    # Everything goes passive: the region still prunes swarms, and nothing
    # unsatisfiable is demanded of the right answer. Forcing an active set
    # here is how 6Q0W ended up with a crystal pose satisfying none of its
    # own restraints.
    if not getattr(face, "directional", True):
        return RestraintSet(active=[], passive=[r.label for r in rim])

    by_approach = sorted(rim, key=lambda r: -r.normal_projection)
    n_active = min(max_active, max(1, len(rim) // 3))
    active = [r.label for r in by_approach[:n_active]]
    chosen = set(active)
    passive = [r.label for r in rim if r.label not in chosen]
    return RestraintSet(active=active, passive=passive)


# ---------------------------------------------------------------------------
# PDB output
# ---------------------------------------------------------------------------


def write_pdb(atoms: list[Atom], path: Path) -> int:
    """Write atoms as a PDB LightDock will accept.

    Column-exact, because PDB is a fixed-column format and a reader that
    tolerates sloppiness is not the one you will meet. Atom names follow the
    convention that a one-to-three character name starts in column 14, which
    is what puts the element in the right place for a reader that infers it
    from the name rather than from columns 77-78.

    **The record type comes from what the residue is, not from how the
    depositor flagged it.** A modified polymer residue -- phosphothreonine,
    selenomethionine, pyroglutamate -- is HETATM in most depositions while
    being part of the backbone. CDK12 has a phosphothreonine in its
    activation loop; writing that as HETATM invites any downstream reader to
    treat it as a ligand and drop it, which leaves a one-residue hole in the
    middle of the protein being docked. So polymer residues are written as
    ATOM whatever their input flag, and only genuine heteroatoms get HETATM.

    This function writes every atom it is given. Deciding what to leave out
    is the caller's job -- a filter here is a filter nobody sees.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with path.open("w", encoding="ascii", errors="replace") as fh:
        for a in atoms:
            polymer = (one_letter(a.comp_id) is not None
                       or a.comp_id.upper() in MODIFIED_AS_POLYMER)
            name = a.name.strip()
            if len(name) >= 4:
                field_name = name[:4]
            else:
                field_name = " " + name.ljust(3)
            fh.write(
                "%-6s%5d %-4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f"
                "          %2s\n" % (
                    "ATOM" if polymer else "HETATM",
                    (n + 1) % 100000,
                    field_name,
                    (a.insertion or " ")[:1] if False else " ",
                    a.comp_id[:3].rjust(3),
                    (a.chain or "A")[:1],
                    a.seq_id,
                    (a.insertion or " ")[:1],
                    a.x, a.y, a.z,
                    a.occupancy if a.occupancy else 1.00,
                    0.00,
                    a.element[:2].upper().rjust(2),
                )
            )
            n += 1
        fh.write("END\n")
    return n


def split_for_docking(
    structure: Structure,
    anchor_chains: list[str],
    partner_chains: list[str],
    out_dir: Path,
    ligand: object | None = None,
    partner_structure: Structure | None = None,
) -> dict:
    """Write the receptor and partner PDBs a LightDock run needs.

    The two proteins may arrive in one file (a ternary complex being taken
    apart to re-dock, which is what the validation harness does) or in two,
    which is the normal case and the one the upload screen produces. Pass
    `partner_structure` for the second; chain names in two separately-solved
    structures frequently collide, and looking for B's chains inside A finds
    either nothing or the wrong thing.

    The receptor is protein only. The ligand small molecule is written
    separately so it can be used to re-score and to rebuild the ternary
    complex afterwards, but it never goes into the docking itself.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    partner_src = partner_structure if partner_structure is not None else structure

    anchor = [a for r in structure.residues if r.chain in set(anchor_chains)
              and r.comp_id not in NON_LIGAND for a in r.atoms]
    partner = [a for r in partner_src.residues if r.chain in set(partner_chains)
               and r.comp_id not in NON_LIGAND for a in r.atoms]

    if not anchor:
        raise ValueError(
            "no anchor residues in chain(s) %s" % ", ".join(anchor_chains))
    if not partner:
        have = ", ".join(sorted(partner_src.chain_sizes())) or "none"
        raise ValueError(
            "no partner residues in chain(s) %s. That structure has: %s"
            % (", ".join(partner_chains), have)
        )

    rec = out_dir / "receptor.pdb"
    lig = out_dir / "partner.pdb"
    n_rec = write_pdb(anchor, rec)
    n_lig = write_pdb(partner, lig)

    out = {
        "receptor": str(rec), "partner": str(lig),
        "receptor_atoms": n_rec, "partner_atoms": n_lig,
        "anchor_chains": list(anchor_chains),
        "partner_chains": list(partner_chains),
    }
    if ligand is not None and getattr(ligand, "atoms", None):
        small = out_dir / "ligand.pdb"
        out["ligand"] = str(small)
        out["ligand_atoms"] = write_pdb(list(ligand.atoms), small)
    return out


# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------


def missing_tools() -> list[str]:
    """Which LightDock executables are not on PATH.

    Reported as a list rather than a boolean because "lightdock is not
    installed" and "lightdock is installed but its scripts are not on your
    PATH" are different problems with different fixes, and on Windows the
    second is the common one.
    """
    return [t for t in TOOLS if shutil.which(t) is None]


def setup_command(
    receptor: Path, partner: Path, restraints: Path, preset: Preset
) -> list[str]:
    return [
        "lightdock3_setup.py", str(receptor), str(partner),
        str(preset.swarms), str(preset.glowworms),
        "-anm", "--noxt", "--noh",
        "-rst", str(restraints),
    ]


def run_command(preset: Preset, cores: int = 0) -> list[str]:
    import os
    cores = cores or max(1, (os.cpu_count() or 2) - 1)
    return [
        "lightdock3.py", "setup.json", str(preset.steps),
        "-s", "fastdfire", "-c", str(cores),
    ]


def conformations_command(
    receptor: Path, partner: Path, preset: Preset
) -> list[str]:
    return [
        "lgd_generate_conformations.py",
        str(receptor), str(partner),
        f"gso_{preset.steps}.out", str(preset.glowworms),
    ]


def cluster_command(preset: Preset) -> list[str]:
    return ["lgd_cluster_bsas.py", f"gso_{preset.steps}.out"]


def rank_command(n_swarms: int, preset: Preset) -> list[str]:
    return ["lgd_rank.py", str(n_swarms), str(preset.steps)]


def count_swarms(work: Path) -> int:
    return len([d for d in work.glob("swarm_*") if d.is_dir()])


def run(cmd: list[str], cwd: Path, log=None, timeout: int = 14400) -> int:
    """Run one LightDock step, streaming its output to `log`.

    Returns the exit code rather than raising, because a LightDock step that
    fails usually says something useful on stdout and the caller wants that
    in front of the user rather than a traceback about a CalledProcessError.
    """
    if log:
        log("$ " + " ".join(cmd))
    try:
        proc = subprocess.Popen(
            cmd, cwd=str(cwd), stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, text=True, bufsize=1,
        )
    except OSError as exc:
        if log:
            log(f"could not start {cmd[0]}: {exc}")
        return 127
    assert proc.stdout is not None
    for line in proc.stdout:
        if log:
            log(line.rstrip("\n"))
    try:
        proc.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        proc.kill()
        if log:
            log(f"{cmd[0]} exceeded {timeout}s and was stopped")
        return 124
    return proc.returncode
