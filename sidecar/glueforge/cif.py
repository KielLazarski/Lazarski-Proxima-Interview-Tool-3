"""A small mmCIF reader, written so the sidecar has no binary dependencies.

It reads only what this pipeline needs: the atom_site loop, the entity and
chem_comp tables, and a handful of header items. That is a deliberate limit,
not an oversight -- a general mmCIF parser is a much larger thing and nothing
downstream of here asks for one.

Handles the parts of the format that actually turn up in PDB entries:
multi-line loops, quoted values, semicolon text fields, '.' and '?' nulls,
and multiple models.
"""

from __future__ import annotations

import gzip
import io
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator, Sequence

from .residues import (
    MODIFIED_AS_POLYMER,
    NON_LIGAND,
    THREE_TO_ONE,
    one_letter,
)

NULLS = {".", "?"}


# ---------------------------------------------------------------------------
# Tokenising
# ---------------------------------------------------------------------------


def _split_tokens(line: str) -> list[str]:
    """Split one mmCIF data line into values, respecting quotes."""
    out: list[str] = []
    i, n = 0, len(line)
    while i < n:
        c = line[i]
        if c in " \t":
            i += 1
            continue
        if c == "#":
            break
        if c in "'\"":
            quote = c
            i += 1
            start = i
            # A quote closes only when followed by whitespace or end of line,
            # which is what lets values such as 5'-ATP survive intact.
            while i < n:
                if line[i] == quote and (i + 1 >= n or line[i + 1] in " \t"):
                    break
                i += 1
            out.append(line[start:i])
            i += 1
            continue
        start = i
        while i < n and line[i] not in " \t":
            i += 1
        out.append(line[start:i])
    return out


def _open_text(path: str | Path) -> io.TextIOBase:
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


# ---------------------------------------------------------------------------
# Loop reader
# ---------------------------------------------------------------------------


def _collect_tags(path: str | Path, prefix: str) -> tuple[list[str], int]:
    """Find a loop's column names and the line index where its data starts."""
    with _open_text(path) as fh:
        tags: list[str] = []
        seen_loop = False
        for idx, raw in enumerate(fh):
            stripped = raw.strip()
            if stripped == "loop_":
                seen_loop = True
                tags = []
                continue
            if seen_loop and stripped.startswith("_"):
                if stripped.startswith(prefix):
                    tags.append(stripped.split()[0][len(prefix):])
                elif tags:
                    seen_loop = False
                    tags = []
                else:
                    seen_loop = False
                continue
            if tags:
                return tags, idx
            seen_loop = False
    return [], -1


def read_loop(path: str | Path, category: str) -> Iterator[dict[str, str]]:
    """Stream the rows of one mmCIF loop. See `iter_loop` for the contract."""
    prefix = category if category.endswith(".") else category + "."
    tags, start = _collect_tags(path, prefix)
    if not tags:
        return

    with _open_text(path) as fh:
        buffer: list[str] = []
        pending: list[str] | None = None
        for idx, raw in enumerate(fh):
            if idx < start:
                continue
            line = raw.rstrip("\n").rstrip("\r")

            if pending is not None:
                if line.startswith(";"):
                    buffer.append("\n".join(pending))
                    pending = None
                else:
                    pending.append(line)
                continue

            stripped = line.strip()
            if not stripped:
                continue
            if stripped == "#":
                break
            if stripped.startswith("loop_") or stripped.startswith("data_"):
                break
            if stripped.startswith("_"):
                break
            if stripped.startswith(";"):
                pending = [stripped[1:]]
                continue

            buffer.extend(_split_tokens(line))
            while len(buffer) >= len(tags):
                yield dict(zip(tags, buffer[: len(tags)]))
                del buffer[: len(tags)]


def read_items(path: str | Path, prefix: str) -> dict[str, str]:
    """Read the non-loop `_category.item value` pairs for one category."""
    out: dict[str, str] = {}
    pref = prefix if prefix.endswith(".") else prefix + "."
    with _open_text(path) as fh:
        pending_key: str | None = None
        pending_val: list[str] | None = None
        for raw in fh:
            line = raw.rstrip("\n")
            if pending_val is not None:
                if line.startswith(";"):
                    out[pending_key] = "\n".join(pending_val).strip()
                    pending_key, pending_val = None, None
                else:
                    pending_val.append(line)
                continue
            stripped = line.strip()
            if not stripped.startswith(pref):
                continue
            parts = _split_tokens(stripped)
            key = parts[0][len(pref):]
            if len(parts) > 1:
                out[key] = parts[1]
            else:
                pending_key, pending_val = key, []
    return out


# ---------------------------------------------------------------------------
# Structure model
# ---------------------------------------------------------------------------


@dataclass(slots=True)
class Atom:
    serial: int
    name: str
    element: str
    comp_id: str
    chain: str
    seq_id: int
    x: float
    y: float
    z: float
    occupancy: float
    hetatm: bool
    insertion: str = ""


@dataclass(slots=True)
class Residue:
    comp_id: str
    chain: str
    seq_id: int
    insertion: str
    one: str
    atoms: list[Atom] = field(default_factory=list)

    @property
    def key(self) -> tuple[str, int, str]:
        return (self.chain, self.seq_id, self.insertion)

    def atom(self, name: str) -> Atom | None:
        for a in self.atoms:
            if a.name == name:
                return a
        return None

    def representative(self) -> Atom:
        """CB where it exists, otherwise CA, otherwise the first atom.

        This is the site the contact map is built on, matching the usual
        coarse-grained convention; glycine has no CB so it falls back to CA.
        """
        return self.atom("CB") or self.atom("CA") or self.atoms[0]

    def centroid(self) -> tuple[float, float, float]:
        n = len(self.atoms)
        return (
            sum(a.x for a in self.atoms) / n,
            sum(a.y for a in self.atoms) / n,
            sum(a.z for a in self.atoms) / n,
        )


@dataclass(slots=True)
class Ligand:
    comp_id: str
    chain: str
    seq_id: int
    atoms: list[Atom] = field(default_factory=list)

    @property
    def heavy_atom_count(self) -> int:
        return sum(1 for a in self.atoms if a.element.upper() != "H")

    def centroid(self) -> tuple[float, float, float]:
        n = len(self.atoms)
        return (
            sum(a.x for a in self.atoms) / n,
            sum(a.y for a in self.atoms) / n,
            sum(a.z for a in self.atoms) / n,
        )


@dataclass
class Structure:
    entry_id: str
    residues: list[Residue]
    ligands: list[Ligand]
    waters_removed: int = 0
    #: auth chain id -> the entity description from the mmCIF header.
    chain_names: dict[str, str] = field(default_factory=dict)
    #: auth chain id -> mmCIF entity id, so copies of one protein group.
    chain_entity: dict[str, str] = field(default_factory=dict)

    def name_of(self, chain: str) -> str:
        return self.chain_names.get(chain, f"chain {chain}")

    def chains(self) -> dict[str, list[Residue]]:
        out: dict[str, list[Residue]] = {}
        for r in self.residues:
            out.setdefault(r.chain, []).append(r)
        return out

    def chain_sizes(self) -> dict[str, int]:
        return {c: len(rs) for c, rs in self.chains().items()}

    def sequence(self, chain: str) -> str:
        return "".join(
            r.one for r in self.residues if r.chain == chain
        )

    def select(self, chains: Sequence[str]) -> list[Residue]:
        want = set(chains)
        return [r for r in self.residues if r.chain in want]

    def largest_ligand(self) -> Ligand | None:
        candidates = [
            lig for lig in self.ligands
            if lig.comp_id.upper() not in NON_LIGAND
            and lig.heavy_atom_count >= 8
        ]
        if not candidates:
            return None
        return max(candidates, key=lambda l: l.heavy_atom_count)

    def ligand_by_id(self, comp_id: str) -> Ligand | None:
        up = comp_id.upper()
        hits = [l for l in self.ligands if l.comp_id.upper() == up]
        if not hits:
            return None
        return max(hits, key=lambda l: l.heavy_atom_count)


def read_pdb(
    path: str | Path,
    model: int = 1,
    keep_waters: bool = False,
) -> Structure:
    """Parse a legacy PDB file.

    Column-delimited, so the field positions below are the format spec rather
    than a guess. Chain names are one character and there are no entity
    descriptions, so chains are named after themselves; everything downstream
    works the same.
    """
    path = Path(path)
    residues: dict[tuple[str, int, str], Residue] = {}
    ligands: dict[tuple[str, int, str], Ligand] = {}
    waters = 0
    seen_altloc: set[tuple[str, int, str, str]] = set()
    current_model = 1

    with _open_text(path) as fh:
        for line in fh:
            tag = line[:6]
            if tag == "MODEL ":
                try:
                    current_model = int(line[10:14])
                except ValueError:
                    current_model = 1
                continue
            if tag == "ENDMDL" and current_model >= model:
                break
            if tag not in ("ATOM  ", "HETATM"):
                continue
            if current_model != model:
                continue

            comp_id = line[17:20].strip().upper()
            if comp_id in ("HOH", "DOD", "WAT"):
                waters += 1
                if not keep_waters:
                    continue

            alt = line[16]
            name = line[12:16].strip()
            chain = line[21].strip() or "A"
            ins = line[26].strip()
            try:
                seq_id = int(line[22:26])
                x, y, z = (float(line[30:38]), float(line[38:46]),
                           float(line[46:54]))
            except ValueError:
                continue

            key4 = (chain, seq_id, ins, name)
            if alt not in (" ", ""):
                if key4 in seen_altloc:
                    continue
                seen_altloc.add(key4)

            element = line[76:78].strip()
            if not element:
                # Older files omit it; the atom name's first letter is the
                # usual fallback, skipping a leading digit.
                element = name.lstrip("0123456789")[:1] or "C"

            try:
                occ = float(line[54:60])
            except ValueError:
                occ = 1.0

            atom = Atom(
                serial=int(line[6:11]) if line[6:11].strip().isdigit() else 0,
                name=name, element=element, comp_id=comp_id, chain=chain,
                seq_id=seq_id, x=x, y=y, z=z, occupancy=occ,
                hetatm=(tag == "HETATM"), insertion=ins,
            )

            key = (chain, seq_id, ins)
            is_polymer = comp_id in THREE_TO_ONE and (
                not atom.hetatm or comp_id in MODIFIED_AS_POLYMER
            )
            if is_polymer:
                res = residues.get(key)
                if res is None:
                    res = Residue(comp_id=comp_id, chain=chain, seq_id=seq_id,
                                  insertion=ins, one=one_letter(comp_id) or "X")
                    residues[key] = res
                res.atoms.append(atom)
            else:
                lig = ligands.get(key)
                if lig is None:
                    lig = Ligand(comp_id=comp_id, chain=chain, seq_id=seq_id)
                    ligands[key] = lig
                lig.atoms.append(atom)

    ordered = [r for r in sorted(residues.values(),
                                 key=lambda r: (r.chain, r.seq_id, r.insertion))
               if r.atoms]
    chains = {r.chain for r in ordered}
    return Structure(
        entry_id=path.stem.upper(),
        residues=ordered,
        ligands=sorted(ligands.values(), key=lambda l: (l.chain, l.seq_id)),
        waters_removed=0 if keep_waters else waters,
        chain_names={c: f"chain {c}" for c in chains},
        chain_entity={c: c for c in chains},
    )


def read_any(
    path: str | Path, model: int = 1, keep_waters: bool = False
) -> Structure:
    """Read mmCIF or PDB, chosen by extension then by sniffing the content."""
    path = Path(path)
    suffixes = [s.lower() for s in path.suffixes]
    if ".pdb" in suffixes or ".ent" in suffixes:
        return read_pdb(path, model, keep_waters)
    if ".cif" in suffixes or ".mmcif" in suffixes:
        return read_structure(path, model, keep_waters)
    with _open_text(path) as fh:
        head = fh.read(4000)
    if "_atom_site." in head or head.lstrip().startswith("data_"):
        return read_structure(path, model, keep_waters)
    return read_pdb(path, model, keep_waters)


def read_structure(
    path: str | Path,
    model: int = 1,
    keep_waters: bool = False,
) -> Structure:
    """Parse one mmCIF file into residues and ligands.

    Keeps the first altloc of any atom, uses author chain and sequence
    numbering where present, and treats modified residues such as MSE as
    polymer rather than as ligands.
    """
    path = Path(path)
    entry = read_items(path, "_entry.").get("id", path.stem.upper())

    residues: dict[tuple[str, int, str], Residue] = {}
    ligands: dict[tuple[str, int, str], Ligand] = {}
    waters = 0
    seen_altloc: set[tuple[str, int, str, str]] = set()
    chain_entity: dict[str, str] = {}

    entity_desc: dict[str, str] = {}
    for row in read_loop(path, "_entity."):
        desc = row.get("pdbx_description", "")
        if desc not in NULLS:
            entity_desc[row.get("id", "")] = desc

    for row in read_loop(path, "_atom_site."):
        if row.get("pdbx_PDB_model_num", "1") not in (str(model), ".", "?"):
            continue

        comp_id = row.get("label_comp_id", "UNK").upper()
        if comp_id in ("HOH", "DOD", "WAT"):
            waters += 1
            if not keep_waters:
                continue

        alt = row.get("label_alt_id", ".")
        chain = row.get("auth_asym_id") or row.get("label_asym_id") or "A"
        if chain in NULLS:
            chain = row.get("label_asym_id", "A")

        raw_seq = row.get("auth_seq_id") or row.get("label_seq_id") or "0"
        if raw_seq in NULLS:
            raw_seq = row.get("label_seq_id", "0")
        if raw_seq in NULLS:
            raw_seq = "0"
        try:
            seq_id = int(raw_seq)
        except ValueError:
            continue

        ins = row.get("pdbx_PDB_ins_code", "")
        if ins in NULLS:
            ins = ""

        name = row.get("label_atom_id", "").strip('"')
        key4 = (chain, seq_id, ins, name)
        if alt not in NULLS:
            if key4 in seen_altloc:
                continue
            seen_altloc.add(key4)

        try:
            atom = Atom(
                serial=int(row.get("id", "0")),
                name=name,
                element=row.get("type_symbol", "C"),
                comp_id=comp_id,
                chain=chain,
                seq_id=seq_id,
                x=float(row["Cartn_x"]),
                y=float(row["Cartn_y"]),
                z=float(row["Cartn_z"]),
                occupancy=float(row.get("occupancy", "1") or 1),
                hetatm=row.get("group_PDB", "ATOM") == "HETATM",
                insertion=ins,
            )
        except (KeyError, ValueError):
            continue

        ent = row.get("label_entity_id", "")
        if ent not in NULLS and chain not in chain_entity:
            chain_entity[chain] = ent

        key = (chain, seq_id, ins)
        is_polymer = comp_id in THREE_TO_ONE and (
            not atom.hetatm or comp_id in MODIFIED_AS_POLYMER
        )

        if is_polymer:
            res = residues.get(key)
            if res is None:
                res = Residue(
                    comp_id=comp_id,
                    chain=chain,
                    seq_id=seq_id,
                    insertion=ins,
                    one=one_letter(comp_id) or "X",
                )
                residues[key] = res
            res.atoms.append(atom)
        else:
            lig = ligands.get(key)
            if lig is None:
                lig = Ligand(comp_id=comp_id, chain=chain, seq_id=seq_id)
                ligands[key] = lig
            lig.atoms.append(atom)

    ordered = sorted(residues.values(), key=lambda r: (r.chain, r.seq_id, r.insertion))
    ordered = [r for r in ordered if r.atoms]

    return Structure(
        entry_id=entry,
        residues=ordered,
        ligands=sorted(ligands.values(), key=lambda l: (l.chain, l.seq_id)),
        waters_removed=0 if keep_waters else waters,
        chain_names={
            c: entity_desc.get(e, f"chain {c}") for c, e in chain_entity.items()
        },
        chain_entity=chain_entity,
    )
