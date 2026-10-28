"""Residue-level reference data.

Everything here is a published scale or a plain chemical fact, held in one
place so the Python reference engine and the JavaScript engine can be checked
against identical inputs. No parameter in this file is fitted or invented.
"""

# ---------------------------------------------------------------------------
# Identity
# ---------------------------------------------------------------------------

#: Canonical ordering used for any residue-pair matrix in this package.
#: Matches the conventional one-letter alphabetical ordering.
AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")

AA_INDEX = {a: i for i, a in enumerate(AA_ORDER)}

THREE_TO_ONE = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
    # Common modified residues mapped to their parent, as the PDB does.
    "MSE": "M",  # selenomethionine
    "SEP": "S",  # phosphoserine
    "TPO": "T",  # phosphothreonine
    "PTR": "Y",  # phosphotyrosine
    "CSO": "C", "CME": "C", "CSD": "C", "OCS": "C",
    "MLY": "K", "KCX": "K", "LLP": "K", "M3L": "K",
    "HYP": "P", "PCA": "Q", "SAH": "M", "ABA": "A",
    "HIC": "H", "NEP": "H", "HIP": "H",
    "ASX": "N", "GLX": "Q", "SEC": "C", "PYL": "K",
}

#: Residues the structure layer treats as polymer even when flagged HETATM.
MODIFIED_AS_POLYMER = {
    k for k in THREE_TO_ONE
    if k not in {
        "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU",
        "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR",
    }
}

# ---------------------------------------------------------------------------
# Scales
# ---------------------------------------------------------------------------

#: Kyte & Doolittle hydropathy, J Mol Biol 157:105-132 (1982).
KD_HYDROPATHY = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5, "M": 1.9, "A": 1.8,
    "G": -0.4, "T": -0.7, "S": -0.8, "W": -0.9, "Y": -1.3, "P": -1.6,
    "H": -3.2, "E": -3.5, "Q": -3.5, "D": -3.5, "N": -3.5, "K": -3.9,
    "R": -4.5,
}

#: Kyte-Doolittle rescaled onto [0, 1]; 1.0 is the most hydrophobic (Ile).
KD_NORMALISED = {
    a: (v - (-4.5)) / (4.5 - (-4.5)) for a, v in KD_HYDROPATHY.items()
}

#: Formal sidechain charge at pH 7. Histidine carries a partial charge
#: because its pKa sits close to physiological pH.
SIDECHAIN_CHARGE = {
    "D": -1.0, "E": -1.0, "K": 1.0, "R": 1.0, "H": 0.1,
    **{a: 0.0 for a in "ACFGILMNPQSTVWY"},
}

#: Sidechain hydrogen-bond capability, as a 0/1 indicator.
HB_DONOR = {a: 0.0 for a in AA_ORDER}
HB_DONOR.update({a: 1.0 for a in "RKWNQHSTY"})

HB_ACCEPTOR = {a: 0.0 for a in AA_ORDER}
HB_ACCEPTOR.update({a: 1.0 for a in "DENQHSTY"})

#: Residue volumes in cubic angstroms, Zamyatnin, Prog Biophys Mol Biol
#: 24:107-123 (1972).
RESIDUE_VOLUME = {
    "A": 88.6, "R": 173.4, "N": 114.1, "D": 111.1, "C": 108.5,
    "Q": 143.8, "E": 138.4, "G": 60.1, "H": 153.2, "I": 166.7,
    "L": 166.7, "K": 168.6, "M": 162.9, "F": 189.9, "P": 112.7,
    "S": 89.0, "T": 116.1, "W": 227.8, "Y": 193.6, "V": 140.0,
}

#: Van der Waals radii in angstroms, by element, for solvent accessibility.
#: Bondi, J Phys Chem 68:441-451 (1964).
VDW_RADII = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80,
    "F": 1.47, "CL": 1.75, "BR": 1.85, "I": 1.98, "SE": 1.90,
    "ZN": 1.39, "MG": 1.73, "NA": 2.27, "K": 2.75, "CA": 2.31,
    "FE": 2.00, "MN": 2.05, "CU": 1.40, "CO": 2.00, "NI": 1.63,
}
DEFAULT_VDW = 1.70

# ---------------------------------------------------------------------------
# Heteroatom classification
# ---------------------------------------------------------------------------

#: Components that are never the ligand of interest: solvent, buffer,
#: cryoprotectant, ions and common crystallisation additives.
SOLVENT_AND_IONS = {
    "HOH", "DOD", "WAT", "H2O",
    "ZN", "MG", "CA", "NA", "K", "CL", "BR", "IOD", "FE", "FE2", "MN",
    "CU", "CU1", "NI", "CO", "CD", "HG", "CS", "RB", "SR", "BA", "AU",
    "SO4", "PO4", "NO3", "ACT", "FMT", "CIT", "TAR", "MES", "TRS", "EPE",
    "GOL", "EDO", "PEG", "PG4", "PGE", "1PE", "2PE", "P6G", "PE4", "MPD",
    "DMS", "IPA", "ETA", "ACY", "BME", "DTT", "TCE", "AZI", "SCN", "CAC",
    "IMD", "BTB", "MLI", "OXL", "SIN", "FLC", "UNX", "UNL",
}

#: Cofactors and nucleotides: present and real, but not the glue.
COFACTORS = {
    "ATP", "ADP", "AMP", "ANP", "ACP", "AGS", "GTP", "GDP", "GNP", "GSP",
    "NAD", "NAP", "NDP", "FAD", "FMN", "SAM", "SAH", "COA", "HEM", "HEC",
    "PLP", "TPP", "BTN", "UDP", "UTP", "CTP", "TTP",
}

#: Glycans, which appear as HETATM but belong to the protein.
GLYCANS = {
    "NAG", "NDG", "BMA", "MAN", "GAL", "GLC", "FUC", "XYP", "SIA", "BGC",
}

NON_LIGAND = SOLVENT_AND_IONS | COFACTORS | GLYCANS


def one_letter(comp_id: str) -> str | None:
    """Return the one-letter code for a component, or None if not a residue."""
    return THREE_TO_ONE.get(comp_id.upper())


def element_radius(element: str) -> float:
    """Van der Waals radius for an element symbol, with a carbon default."""
    return VDW_RADII.get(element.upper().strip(), DEFAULT_VDW)
