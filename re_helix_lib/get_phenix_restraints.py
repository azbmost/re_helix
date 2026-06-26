#!/usr/bin/env python3
"""
get_phenix_restraints.py

Read LINK records from a PDB file and convert them into Phenix
geometry_restraints.edits bond restraints in a .params file. In addition,
for LINK records involving phosphate P--O3' or P--O5' bonds, generate
phosphate-centered angle restraints that preserve the local 3D geometry of
P, OP1/OP2 (or legacy O1P/O2P), O5', and the LINK-derived O3'/O5' atom.

V4 update
---------
- Add explicit restraints for standalone X33 3'-3' linker phosphate residues:
  P-OP1/P-OP2 bond restraints and phosphate-centered angle restraints.
- Parse REMARK 950 RE_SCRIPT JUNCTION records written by the newer reciprocal
  exchange / aligner scripts and write a second <input_stem>_junctions.params
  file.  This file keeps the usual P/OP/O5'/C5'/O3' minimization selection for
  all residues, while adding C4', O4', C3', C2', and O2' only for the parsed
  junction residues.

V3 update
---------
V2 generated every pairwise O-P-O angle around each LINK-associated phosphate.
Phenix already supplies several of those angles from standard nucleotide
monomer geometry and standard previous-residue O3'-P chain-link geometry.
Re-adding them can make phenix.geometry_minimization stop with duplicate-angle
errors. V3 omits those built-in-like angle candidates by default while keeping
the custom LINK-dependent angles needed for nonstandard phosphate geometry.

Usage
-----
  python get_phenix_restraints.py model_with_LINKs.pdb
  python get_phenix_restraints.py model.pdb --output-base model
  python get_phenix_restraints.py --gui

The output file can be passed directly to phenix.geometry_minimization, e.g.:

  phenix.geometry_minimization model_with_LINKs.pdb model_links.params

Notes
-----
- Only LINK records are parsed; existing CONECT records are ignored.
- The numeric distance at the end of the LINK line is *ignored*; instead
  we assign a simple "ideal" distance based on the element types of the
  two linked atoms (P-O, C-O, C-N, etc.), with a configurable sigma.
- Phosphate angle restraints are generated only for phosphate groups whose
  geometry can be inferred from LINK records plus atoms present in the input
  PDB. The output atom names for non-bridging phosphate oxygens follow the
  input PDB convention: OP1/OP2 or O1P/O2P.
- The two special phosphate-linkage cases are handled explicitly and reported
  to the terminal when detected:
    * 5'-to-5': one LINK-derived O5' plus the same-residue O5' bind the P.
    * 3'-to-3': two LINK-derived O3' atoms bind an extra linker phosphate P.
- By default, angle restraints already covered by standard Phenix nucleotide
  geometry are not written. This avoids fatal duplicate-angle errors.
- This is intended for non-standard 3'-3' or 5'-5' phosphate linkages
  and 3'-5' backbone LINKs (e.g. from bowtie_link.py).
- X33 is treated as a standalone HETATM linker residue containing exactly
  P, OP1, and OP2; its internal geometry is restrained explicitly because
  Phenix will not have a standard monomer definition for this custom hetID.
"""

import argparse
import os
import shlex
import sys
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set

try:
    import tkinter as tk
    from tkinter import filedialog, messagebox, scrolledtext
except Exception:
    tk = None
    filedialog = None
    messagebox = None
    scrolledtext = None


class LinkRecord:
    def __init__(
        self,
        atom_name_1: str,
        res_name_1: str,
        chain_id_1: str,
        res_seq_1: int,
        atom_name_2: str,
        res_name_2: str,
        chain_id_2: str,
        res_seq_2: int,
        original_distance: Optional[float] = None,
        raw_line: str = "",
    ) -> None:
        self.atom_name_1 = atom_name_1
        self.res_name_1 = res_name_1
        self.chain_id_1 = chain_id_1
        self.res_seq_1 = res_seq_1

        self.atom_name_2 = atom_name_2
        self.res_name_2 = res_name_2
        self.chain_id_2 = chain_id_2
        self.res_seq_2 = res_seq_2

        self.original_distance = original_distance
        self.raw_line = raw_line

    def key(self) -> Tuple[str, str, int, str, str, int]:
        """
        Canonical key for avoiding duplicate restraints.
        Order is sorted by (chain,resSeq,atom_name) so that A-B and B-A
        are treated as the same link.
        """
        a1 = (self.chain_id_1, self.res_seq_1, self.atom_name_1.strip())
        a2 = (self.chain_id_2, self.res_seq_2, self.atom_name_2.strip())
        return tuple(sorted((a1, a2)))


class AtomRef:
    """Minimal atom identifier used to build atom_selection strings."""

    def __init__(
        self,
        atom_name: str,
        res_name: str,
        chain_id: str,
        res_seq: int,
    ) -> None:
        self.atom_name = atom_name.strip()
        self.res_name = res_name.strip()
        self.chain_id = chain_id.strip()
        self.res_seq = res_seq

    def key(self) -> Tuple[str, int, str]:
        return (self.chain_id, self.res_seq, normalize_atom_name(self.atom_name))

    def selection(self) -> str:
        return (
            f"chain {self.chain_id} and resid {self.res_seq} "
            f"and name {atom_name(self.atom_name)}"
        )

    def label(self) -> str:
        return f"{self.atom_name}:{self.res_name}:{self.chain_id}{self.res_seq}"


AtomIndex = Dict[Tuple[str, int], Dict[str, AtomRef]]

NONBRIDGING_1_NAMES = ("OP1", "O1P")
NONBRIDGING_2_NAMES = ("OP2", "O2P")
NONBRIDGING_NAMES = NONBRIDGING_1_NAMES + NONBRIDGING_2_NAMES
O3_NAMES = ("O3'", "O3*")
O5_NAMES = ("O5'", "O5*")

TOOL_NAME = "Get Phenix Restraints"
VERSION = "V1.0"

DEFAULT_LINKER_RESNAME = "X33"
X33_RESNAME = DEFAULT_LINKER_RESNAME
X33_P_OP_DISTANCE = 1.495
DEFAULT_LINK_DISTANCE_CUTOFF = 6.5
DEFAULT_NONBRIDGING_ANGLE = 119.0
DEFAULT_MIXED_ANGLE = 108.0
DEFAULT_BRIDGING_ANGLE = 103.0
REMARK_PREFIX = "REMARK 950 RE_SCRIPT"
STANDARD_NUCLEIC_ACID_RESNAMES = {"A", "C", "G", "U", "DA", "DC", "DG", "DT"}

JUNCTION_BASE_ATOM_NAMES = (" P ", " OP1", " OP2", " O5'", " C5'", " O3'", " O1P", " O2P")
JUNCTION_EXTRA_SUGAR_ATOM_NAMES = (" C4'", " O4'", " C3'", " C2'", " O2'")


def normalize_linker_resname(resname: str) -> str:
    cleaned = (resname or "").strip().upper()
    if not cleaned:
        raise ValueError("The 3'-to-3' linker residue name cannot be blank.")
    if len(cleaned) > 3:
        raise ValueError("The 3'-to-3' linker residue name must fit the PDB residue-name field (1-3 characters).")
    if any(ch.isspace() for ch in cleaned):
        raise ValueError("The 3'-to-3' linker residue name cannot contain spaces.")
    return cleaned


# ---------------------------------------------------------------------------
# Existing bond-restraint helpers
# ---------------------------------------------------------------------------

def _simple_element_from_pdb_name(atom_name: str) -> Optional[str]:
    """
    Very simple guess of element symbol from a 4-character PDB atom name.

    We just return the first alphabetic character; for typical nucleic-acid
    names like 'P', 'O3'', 'C4*' etc., this gives P, O, C, N, etc., which
    is sufficient for rough bond-length guesses.
    """
    for ch in atom_name.strip():
        if ch.isalpha():
            return ch.upper()
    return None


def guess_ideal_distance(atom_name_1: str, atom_name_2: str) -> float:
    """
    Heuristic ideal distance for a covalent bond, based only on atom types.

    This does *not* attempt precise small-molecule geometry; it just gives
    reasonable values for nucleic-acid P-O and C-O/C-N links. All distances
    are in A.
    """
    e1 = _simple_element_from_pdb_name(atom_name_1)
    e2 = _simple_element_from_pdb_name(atom_name_2)
    if e1 is None or e2 is None:
        return 1.6

    pair = {e1, e2}
    if pair == {"P", "O"}:
        # P-O single bond (bridging or terminal)
        return 1.60
    if pair == {"C", "O"}:
        # generic C-O single
        return 1.43
    if pair == {"C", "N"}:
        # generic C-N single
        return 1.33
    if pair == {"C", "C"}:
        # generic C-C single
        return 1.54
    if pair == {"N", "H"}:
        return 1.01
    if pair == {"O", "H"}:
        return 0.98

    # Fallback
    return 1.6


def parse_link_line(line: str) -> Optional[LinkRecord]:
    """
    Parse one LINK line into a LinkRecord.

    We do NOT rely on exact PDB column positions, because the bowtie
    scripts write LINK records using formatted f-strings; instead we
    split on whitespace and interpret the tokens as:

        LINK atom1 resName1 chain1 resSeq1 atom2 resName2 chain2 resSeq2 [sym1] [sym2] [dist]

    For example:

        LINK         O3'   A F   8                   P   A F   9     1555   1555  3.40

    becomes tokens:
        ['LINK', "O3'", 'A', 'F', '8', 'P', 'A', 'F', '9', '1555', '1555', '3.40']
    """
    if not line.startswith("LINK"):
        return None

    tokens = line.split()
    if len(tokens) < 9:
        # Too short to be a standard LINK record
        return None

    if tokens[0] != "LINK":
        return None

    try:
        atom1 = tokens[1]
        resn1 = tokens[2]
        chain1 = tokens[3]
        resseq1 = int(tokens[4])

        atom2 = tokens[5]
        resn2 = tokens[6]
        chain2 = tokens[7]
        resseq2 = int(tokens[8])
    except (IndexError, ValueError):
        return None

    original_distance: Optional[float] = None
    # The last token is often the distance; try to parse but ignore if it fails.
    try:
        original_distance = float(tokens[-1])
    except ValueError:
        original_distance = None

    return LinkRecord(
        atom_name_1=atom1,
        res_name_1=resn1,
        chain_id_1=chain1,
        res_seq_1=resseq1,
        atom_name_2=atom2,
        res_name_2=resn2,
        chain_id_2=chain2,
        res_seq_2=resseq2,
        original_distance=original_distance,
        raw_line=line.rstrip("\n"),
    )


def read_link_records(pdb_path: Path) -> List[LinkRecord]:
    links: List[LinkRecord] = []
    seen_keys: Set[Tuple[str, str, int, str, str, int]] = set()

    with pdb_path.open("r") as fh:
        for line in fh:
            if not line.startswith("LINK"):
                continue
            rec = parse_link_line(line)
            if rec is None:
                continue
            key = rec.key()
            if key in seen_keys:
                # Skip duplicate LINKs between same atom pair
                continue
            seen_keys.add(key)
            links.append(rec)

    return links


def atom_name(name: str) -> str:
    """
    Quote an atom name for use in a Phenix atom_selection string.

    We keep the original script's atom-name formatting so existing bond
    selections remain unchanged. Names containing apostrophes or asterisks
    are passed through exactly as they appear in the input/LINK records.
    """
    return f'{name.strip()}'


# ---------------------------------------------------------------------------
# PDB atom indexing and phosphate angle helpers
# ---------------------------------------------------------------------------

def normalize_atom_name(name: str) -> str:
    """
    Canonicalize atom names only for internal matching.

    The output always uses the original atom name from the PDB/LINK record.
    We treat O3*/O5* as equivalent to O3'/O5' for older-style files.
    """
    return name.strip().upper().replace("*", "'")


def is_phosphorus_name(name: str) -> bool:
    return normalize_atom_name(name) == "P"


def is_o3_name(name: str) -> bool:
    return normalize_atom_name(name) in {normalize_atom_name(n) for n in O3_NAMES}


def is_o5_name(name: str) -> bool:
    return normalize_atom_name(name) in {normalize_atom_name(n) for n in O5_NAMES}


def read_pdb_atom_index(pdb_path: Path) -> AtomIndex:
    """
    Read ATOM/HETATM records and index them by (chainID, resSeq, atom name).

    This is deliberately lightweight and mirrors the original script's lack of
    insertion-code handling. It is used only to discover same-residue phosphate
    atoms such as OP1/OP2 or O1P/O2P.
    """
    atom_index: AtomIndex = {}
    with pdb_path.open("r") as fh:
        for line in fh:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            if len(line) < 26:
                continue
            try:
                atom = AtomRef(
                    atom_name=line[12:16].strip(),
                    res_name=line[17:20].strip(),
                    chain_id=line[21].strip(),
                    res_seq=int(line[22:26]),
                )
            except ValueError:
                continue
            res_key = (atom.chain_id, atom.res_seq)
            atom_index.setdefault(res_key, {})[normalize_atom_name(atom.atom_name)] = atom
    return atom_index


def read_previous_residue_map(
    pdb_path: Path,
) -> Dict[Tuple[str, int], Optional[Tuple[str, int]]]:
    """
    Return the immediately previous residue in file order for each chain.

    Phenix normally supplies phosphate angles involving OP1/OP2 (or O1P/O2P),
    same-residue O5', and the O3' atom of the previous residue in the same
    chain. Those angles must not be re-added as custom edits, otherwise
    phenix.geometry_minimization can stop with duplicate-angle errors. File
    order is used instead of resSeq - 1 so numbering gaps remain sensible.
    """
    residues_by_chain: Dict[str, List[Tuple[str, int]]] = {}
    seen: Set[Tuple[str, int]] = set()

    with pdb_path.open("r") as fh:
        for line in fh:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            if len(line) < 26:
                continue
            try:
                chain_id = line[21].strip()
                res_seq = int(line[22:26])
            except ValueError:
                continue
            key = (chain_id, res_seq)
            if key in seen:
                continue
            seen.add(key)
            residues_by_chain.setdefault(chain_id, []).append(key)

    previous_by_residue: Dict[Tuple[str, int], Optional[Tuple[str, int]]] = {}
    for residues in residues_by_chain.values():
        previous: Optional[Tuple[str, int]] = None
        for key in residues:
            previous_by_residue[key] = previous
            previous = key

    return previous_by_residue


def is_nonbridging_name(name: str) -> bool:
    return normalize_atom_name(name) in {normalize_atom_name(n) for n in NONBRIDGING_NAMES}


def same_residue(atom_1: AtomRef, atom_2: AtomRef) -> bool:
    return atom_1.chain_id == atom_2.chain_id and atom_1.res_seq == atom_2.res_seq


def ligand_has_phenix_builtin_phosphate_geometry(
    center_p: AtomRef,
    ligand: AtomRef,
    previous_residue_by_key: Dict[Tuple[str, int], Optional[Tuple[str, int]]],
    linker_resname: str = DEFAULT_LINKER_RESNAME,
) -> bool:
    """Return True if Phenix normally already has P-centered angles to ligand.

    X33 is a custom standalone HETATM linker residue, so none of its
    phosphate-centered angles are assumed to be supplied by standard Phenix
    monomer/backbone geometry.

    Built-in phosphate ligands are:
      - same-residue OP1/OP2 or O1P/O2P,
      - same-residue O5'/O5*, and
      - O3'/O3* from the immediately previous residue in the same chain.

    Angles where both ligands are in this built-in set are skipped by default.
    Angles involving at least one non-standard LINK ligand are retained.
    """
    if center_p.res_name.strip().upper() == linker_resname.strip().upper():
        return False

    if same_residue(center_p, ligand):
        if is_nonbridging_name(ligand.atom_name):
            return True
        if is_o5_name(ligand.atom_name):
            return True

    previous_residue = previous_residue_by_key.get((center_p.chain_id, center_p.res_seq))
    if previous_residue is not None and is_o3_name(ligand.atom_name):
        return (ligand.chain_id, ligand.res_seq) == previous_residue

    return False


def find_atom(
    atom_index: AtomIndex,
    chain_id: str,
    res_seq: int,
    candidate_names: Tuple[str, ...],
) -> Optional[AtomRef]:
    res_atoms = atom_index.get((chain_id.strip(), res_seq), {})
    for name in candidate_names:
        atom = res_atoms.get(normalize_atom_name(name))
        if atom is not None:
            return atom
    return None


def resolve_atom_ref(atom_index: AtomIndex, ref: AtomRef) -> AtomRef:
    """Return the PDB atom matching a LINK endpoint when available."""
    atom = find_atom(atom_index, ref.chain_id, ref.res_seq, (ref.atom_name,))
    if atom is not None:
        return atom
    return ref


def link_endpoint_refs(link: LinkRecord) -> Tuple[AtomRef, AtomRef]:
    return (
        AtomRef(link.atom_name_1, link.res_name_1, link.chain_id_1, link.res_seq_1),
        AtomRef(link.atom_name_2, link.res_name_2, link.chain_id_2, link.res_seq_2),
    )


def get_phosphate_link_endpoint(
    link: LinkRecord,
) -> Optional[Tuple[AtomRef, AtomRef, str]]:
    """
    If a LINK connects phosphate P to O3'/O5', return (P_atom, O_atom, O_type).

    O_type is 'O3' or 'O5'. Other LINK types are ignored for angle generation
    but still receive the original bond restraint.
    """
    a1, a2 = link_endpoint_refs(link)

    if is_phosphorus_name(a1.atom_name) and is_o3_name(a2.atom_name):
        return a1, a2, "O3"
    if is_phosphorus_name(a1.atom_name) and is_o5_name(a2.atom_name):
        return a1, a2, "O5"
    if is_phosphorus_name(a2.atom_name) and is_o3_name(a1.atom_name):
        return a2, a1, "O3"
    if is_phosphorus_name(a2.atom_name) and is_o5_name(a1.atom_name):
        return a2, a1, "O5"

    return None


def angle_ideal_for_ligand_pair(
    ligand_role_1: str,
    ligand_role_2: str,
    nonbridging_angle: float,
    mixed_angle: float,
    bridging_angle: float,
) -> float:
    """Choose a simple phosphate O-P-O ideal angle from ligand roles."""
    if ligand_role_1 == "nonbridging" and ligand_role_2 == "nonbridging":
        return nonbridging_angle
    if ligand_role_1 == "bridging" and ligand_role_2 == "bridging":
        return bridging_angle
    return mixed_angle


def build_angle_block(
    ligand_1: AtomRef,
    center_p: AtomRef,
    ligand_2: AtomRef,
    angle_ideal: float,
    angle_sigma: float,
) -> List[str]:
    lines: List[str] = []
    lines.append("    angle {")
    lines.append("      action = add")
    lines.append(f"      atom_selection_1 = {ligand_1.selection()}")
    lines.append(f"      atom_selection_2 = {center_p.selection()}")
    lines.append(f"      atom_selection_3 = {ligand_2.selection()}")
    lines.append(f"      angle_ideal = {angle_ideal:.3f}")
    lines.append(f"      sigma = {angle_sigma:.3f}")
    lines.append("    }")
    return lines


def build_phosphate_angle_lines(
    pdb_path: Path,
    links: List[LinkRecord],
    angle_sigma: float,
    nonbridging_angle: float,
    mixed_angle: float,
    bridging_angle: float,
    skip_phenix_builtin_angles: bool = True,
    preexisting_angle_keys: Optional[Set[Tuple[Tuple[str, int, str], Tuple[Tuple[str, int, str], Tuple[str, int, str]]]]] = None,
    linker_resname: str = DEFAULT_LINKER_RESNAME,
) -> Tuple[List[str], List[str]]:
    """
    Build phosphate-centered angle restraints derived from LINK records.

    Returns:
        angle_lines: .params lines to append inside geometry_restraints.edits
        reports: terminal messages for detected special linkages and warnings
    """
    atom_index = read_pdb_atom_index(pdb_path)
    previous_residue_by_key = read_previous_residue_map(pdb_path)

    center_groups: Dict[Tuple[str, int], Dict[str, object]] = {}
    reports: List[str] = []

    for link_number, link in enumerate(links, start=1):
        endpoint = get_phosphate_link_endpoint(link)
        if endpoint is None:
            continue

        p_from_link, linked_o_from_link, linked_o_type = endpoint
        center_p = resolve_atom_ref(atom_index, p_from_link)
        linked_o = resolve_atom_ref(atom_index, linked_o_from_link)
        center_key = (center_p.chain_id, center_p.res_seq)

        group = center_groups.setdefault(
            center_key,
            {
                "center": center_p,
                "linked": [],
                "linked_keys": set(),
            },
        )
        linked_keys = group["linked_keys"]
        assert isinstance(linked_keys, set)
        if linked_o.key() in linked_keys:
            continue
        linked_keys.add(linked_o.key())
        linked = group["linked"]
        assert isinstance(linked, list)
        linked.append((linked_o, linked_o_type, link_number))

    angle_lines: List[str] = []
    seen_angle_keys: Set[Tuple[Tuple[str, int, str], Tuple[Tuple[str, int, str], Tuple[str, int, str]]]] = set(
        preexisting_angle_keys or set()
    )
    skipped_builtin_angle_count = 0
    generated_angle_count = 0

    if center_groups:
        angle_lines.append("    # Phosphate angle restraints inferred from LINK-associated phosphate groups")
        if skip_phenix_builtin_angles:
            angle_lines.append("    # V3 duplicate-avoidance: standard Phenix phosphate angles are omitted")
        else:
            angle_lines.append("    # Diagnostic mode: standard Phenix phosphate angles are included")
        angle_lines.append("    # Angle ideals: nonbridging-nonbridging, mixed, bridging-bridging")
        angle_lines.append(
            f"    # {nonbridging_angle:.3f}, {mixed_angle:.3f}, {bridging_angle:.3f} degrees; sigma {angle_sigma:.3f} degrees"
        )
        angle_lines.append("")

    for center_key in sorted(center_groups.keys()):
        group = center_groups[center_key]
        center_p = group["center"]
        linked = group["linked"]
        assert isinstance(center_p, AtomRef)
        assert isinstance(linked, list)

        op1 = find_atom(atom_index, center_p.chain_id, center_p.res_seq, NONBRIDGING_1_NAMES)
        op2 = find_atom(atom_index, center_p.chain_id, center_p.res_seq, NONBRIDGING_2_NAMES)
        same_o5 = find_atom(atom_index, center_p.chain_id, center_p.res_seq, O5_NAMES)

        linked_o3 = [(atom, link_number) for atom, o_type, link_number in linked if o_type == "O3"]
        linked_o5 = [(atom, link_number) for atom, o_type, link_number in linked if o_type == "O5"]

        # Report the two requested exception classes as soon as they are recognized,
        # even if the final angle set must later be skipped because atoms are missing.
        if len(linked_o3) >= 2:
            preview_pair = linked_o3[:2]
            reports.append(
                "Detected 3'-to-3' linker phosphate for angle restraints: "
                f"P {center_p.label()} linked to "
                + ", ".join(atom.label() for atom, _ in preview_pair)
                + "."
            )
        elif linked_o5:
            linked_o5_atom, _ = linked_o5[0]
            same_o5_text = same_o5.label() if same_o5 is not None else "not found"
            reports.append(
                "Detected 5'-to-5' phosphate linkage for angle restraints: "
                f"P {center_p.label()} linked to {linked_o5_atom.label()} "
                f"with same-residue O5' {same_o5_text}."
            )

        if op1 is None or op2 is None:
            reports.append(
                "WARNING: skipped phosphate angle restraints for "
                f"{center_p.label()} because OP1/OP2 (or O1P/O2P) were not both found."
            )
            continue

        linkage_kind = "standard 3'-5'"
        # Ligand tuple fields are: AtomRef, ideal-angle role, is_custom_link_ligand.
        # Custom is False for ligands Phenix normally already knows around P.
        selected_ligands: List[Tuple[AtomRef, str, bool]] = [
            (op1, "nonbridging", False),
            (op2, "nonbridging", False),
        ]
        selected_link_numbers: List[int] = []

        def linked_o3_ligand(atom: AtomRef) -> Tuple[AtomRef, str, bool]:
            is_builtin = ligand_has_phenix_builtin_phosphate_geometry(
                center_p=center_p,
                ligand=atom,
                previous_residue_by_key=previous_residue_by_key,
                linker_resname=linker_resname,
            )
            return (atom, "bridging", not is_builtin)

        if len(linked_o3) >= 2:
            # 3'-to-3' special case: an extra linker phosphate has two O3' ligands.
            linkage_kind = "3'-to-3' linker phosphate"
            selected_pair = linked_o3[:2]
            selected_ligands.extend(linked_o3_ligand(atom) for atom, _ in selected_pair)
            selected_link_numbers = [link_number for _, link_number in selected_pair]
            if len(linked_o3) > 2:
                reports.append(
                    "WARNING: more than two O3' LINK ligands were found for "
                    f"{center_p.label()}; only the first two were used for phosphate angles."
                )
            if same_o5 is not None:
                reports.append(
                    "WARNING: 3'-to-3' linker phosphate "
                    f"{center_p.label()} also has a same-residue O5' atom ({same_o5.label()}); "
                    "the O5' atom was not included in the 3'-to-3' angle set."
                )

        elif linked_o5:
            # 5'-to-5' special case: same-residue O5' plus a LINK-derived O5'.
            linkage_kind = "5'-to-5' phosphate"
            linked_o5_atom, linked_o5_link_number = linked_o5[0]
            if same_o5 is None:
                reports.append(
                    "WARNING: skipped 5'-to-5' phosphate angle restraints for "
                    f"{center_p.label()} because the same-residue O5' atom was not found."
                )
                continue
            selected_ligands.append((same_o5, "bridging", False))
            selected_ligands.append((linked_o5_atom, "bridging", True))
            selected_link_numbers = [linked_o5_link_number]
            if len(linked_o5) > 1:
                reports.append(
                    "WARNING: more than one LINK-derived O5' ligand was found for "
                    f"{center_p.label()}; only the first was used for phosphate angles."
                )

        elif len(linked_o3) == 1:
            # Standard or inverted 3'-5' phosphate geometry inferred from a P--O3' LINK.
            linked_o3_atom, linked_o3_link_number = linked_o3[0]
            if same_o5 is None:
                reports.append(
                    "WARNING: skipped phosphate angle restraints for "
                    f"{center_p.label()} because same-residue O5' was not found."
                )
                continue
            selected_ligands.append((same_o5, "bridging", False))
            selected_ligands.append(linked_o3_ligand(linked_o3_atom))
            selected_link_numbers = [linked_o3_link_number]

        else:
            continue

        unique_ligand_keys = {atom.key() for atom, _role, _custom in selected_ligands}
        if len(unique_ligand_keys) != 4:
            reports.append(
                "WARNING: skipped phosphate angle restraints for "
                f"{center_p.label()} because the four ligand atoms were not unique."
            )
            continue

        center_angle_lines: List[str] = []
        skipped_this_center = 0
        for (ligand_1, role_1, custom_1), (ligand_2, role_2, custom_2) in combinations(selected_ligands, 2):
            sorted_ligand_keys = tuple(sorted((ligand_1.key(), ligand_2.key())))
            angle_key = (center_p.key(), sorted_ligand_keys)
            if angle_key in seen_angle_keys:
                continue

            # Central V3 fix: if both ligands are already part of the standard
            # Phenix phosphate environment, do not re-add the angle as custom.
            if skip_phenix_builtin_angles and not (custom_1 or custom_2):
                skipped_builtin_angle_count += 1
                skipped_this_center += 1
                continue

            seen_angle_keys.add(angle_key)
            ideal = angle_ideal_for_ligand_pair(
                role_1,
                role_2,
                nonbridging_angle=nonbridging_angle,
                mixed_angle=mixed_angle,
                bridging_angle=bridging_angle,
            )
            center_angle_lines.extend(
                build_angle_block(
                    ligand_1=ligand_1,
                    center_p=center_p,
                    ligand_2=ligand_2,
                    angle_ideal=ideal,
                    angle_sigma=angle_sigma,
                )
            )
            center_angle_lines.append("")
            generated_angle_count += 1

        if center_angle_lines:
            angle_lines.append(
                "    # Phosphate angles for "
                f"{linkage_kind} centered at {center_p.label()}"
            )
            if selected_link_numbers:
                angle_lines.append(
                    "    # Derived from LINK record(s): "
                    + ", ".join(str(n) for n in selected_link_numbers)
                )
            if skipped_this_center:
                angle_lines.append(
                    f"    # Skipped {skipped_this_center} Phenix-standard phosphate angle candidate(s) for this center"
                )
            angle_lines.extend(center_angle_lines)

    if center_groups:
        if skipped_builtin_angle_count:
            reports.append(
                "Skipped "
                f"{skipped_builtin_angle_count} phosphate angle candidate(s) likely already present "
                "in Phenix standard nucleotide restraints."
            )
        reports.append(f"Generated {generated_angle_count} phosphate angle restraint(s).")

    return angle_lines, reports



# ---------------------------------------------------------------------------
# Standalone linker-phosphate helpers
# ---------------------------------------------------------------------------

def read_x33_residues(
    pdb_path: Path,
    linker_resname: str = DEFAULT_LINKER_RESNAME,
) -> Dict[Tuple[str, int], Dict[str, AtomRef]]:
    """Return standalone linker residue atoms keyed by (chainID, resSeq)."""
    residues: Dict[Tuple[str, int], Dict[str, AtomRef]] = {}
    linker_resname_upper = linker_resname.strip().upper()
    with pdb_path.open("r") as fh:
        for line in fh:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            if len(line) < 26:
                continue
            res_name = line[17:20].strip()
            if res_name.upper() != linker_resname_upper:
                continue
            try:
                atom = AtomRef(
                    atom_name=line[12:16].strip(),
                    res_name=res_name,
                    chain_id=line[21].strip(),
                    res_seq=int(line[22:26]),
                )
            except ValueError:
                continue
            residues.setdefault((atom.chain_id, atom.res_seq), {})[normalize_atom_name(atom.atom_name)] = atom
    return residues


def angle_key_for_atoms(center_p: AtomRef, ligand_1: AtomRef, ligand_2: AtomRef):
    """Return the same canonical angle key used by build_phosphate_angle_lines()."""
    sorted_ligand_keys = tuple(sorted((ligand_1.key(), ligand_2.key())))
    return (center_p.key(), sorted_ligand_keys)


def build_bond_block(atom_1: AtomRef, atom_2: AtomRef, distance_ideal: float, sigma: float) -> List[str]:
    lines: List[str] = []
    lines.append("    bond {")
    lines.append("      action = add")
    lines.append(f"      atom_selection_1 = {atom_1.selection()}")
    lines.append(f"      atom_selection_2 = {atom_2.selection()}")
    lines.append(f"      distance_ideal = {distance_ideal:.3f}")
    lines.append(f"      sigma = {sigma:.3f}")
    lines.append("    }")
    return lines


def build_x33_internal_restraint_lines(
    pdb_path: Path,
    sigma: float,
    angle_sigma: float,
    p_op_distance: float = X33_P_OP_DISTANCE,
    op_p_op_angle: float = DEFAULT_NONBRIDGING_ANGLE,
    linker_resname: str = DEFAULT_LINKER_RESNAME,
) -> Tuple[List[str], List[str], Set[Tuple[Tuple[str, int, str], Tuple[Tuple[str, int, str], Tuple[str, int, str]]]]]:
    """Build internal P-OP1/P-OP2 bond and OP1-P-OP2 angle restraints.

    Returns the generated lines, terminal report messages, and the canonical
    OP1-P-OP2 angle keys already written so LINK-derived angle generation can
    avoid duplicates even in --include-phenix-builtin-angles mode.
    """
    linker_resname = normalize_linker_resname(linker_resname)
    x33_residues = read_x33_residues(pdb_path, linker_resname=linker_resname)
    lines: List[str] = []
    reports: List[str] = []
    written_angle_keys: Set[Tuple[Tuple[str, int, str], Tuple[Tuple[str, int, str], Tuple[str, int, str]]]] = set()

    if not x33_residues:
        return lines, reports, written_angle_keys

    lines.append(f"    # Internal restraints for standalone {linker_resname} 3'-3' linker phosphate residues")
    lines.append(
        f"    # {linker_resname} contains P, OP1, and OP2; P-O3' LINK bonds are handled by LINK-derived restraints"
    )

    complete_count = 0
    for (chain_id, res_seq), atoms in sorted(x33_residues.items(), key=lambda kv: (kv[0][0], kv[0][1])):
        p_atom = atoms.get("P")
        op1 = atoms.get("OP1") or atoms.get("O1P")
        op2 = atoms.get("OP2") or atoms.get("O2P")
        label = f"{linker_resname}:{chain_id or '_'}{res_seq}"

        if p_atom is None or op1 is None or op2 is None:
            missing = []
            if p_atom is None:
                missing.append("P")
            if op1 is None:
                missing.append("OP1/O1P")
            if op2 is None:
                missing.append("OP2/O2P")
            reports.append(
                f"WARNING: skipped incomplete {linker_resname} internal restraints for "
                f"{label}; missing {', '.join(missing)}."
            )
            continue

        lines.append(f"    # {linker_resname} internal geometry for {p_atom.res_name}:{chain_id or '_'}{res_seq}")
        lines.extend(build_bond_block(p_atom, op1, p_op_distance, sigma))
        lines.append("")
        lines.extend(build_bond_block(p_atom, op2, p_op_distance, sigma))
        lines.append("")
        lines.extend(
            build_angle_block(
                ligand_1=op1,
                center_p=p_atom,
                ligand_2=op2,
                angle_ideal=op_p_op_angle,
                angle_sigma=angle_sigma,
            )
        )
        lines.append("")
        written_angle_keys.add(angle_key_for_atoms(p_atom, op1, op2))
        complete_count += 1

    if complete_count:
        reports.append(
            f"Generated internal {linker_resname} restraints for {complete_count} standalone linker residue(s)."
        )
    if lines and lines[-1] != "":
        lines.append("")
    return lines, reports, written_angle_keys


# ---------------------------------------------------------------------------
# Junction REMARK parsing and junction-selection params helpers
# ---------------------------------------------------------------------------

def pdb_input_stem(pdb_path: Path) -> str:
    """Return an input filename stem, treating .pdb.txt as one PDB-style suffix."""
    name = pdb_path.name
    lower = name.lower()
    if lower.endswith(".pdb.txt"):
        return name[:-8]
    if lower.endswith(".pdb"):
        return name[:-4]
    if pdb_path.suffix:
        return name[: -len(pdb_path.suffix)]
    return name


def parse_remark_key_values(line: str) -> Dict[str, str]:
    """Parse simple key=value tokens from a REMARK 950 RE_SCRIPT line."""
    fields: Dict[str, str] = {}
    for token in line.strip().split():
        if "=" not in token:
            continue
        key, value = token.split("=", 1)
        fields[key] = value
    return fields


def parse_residue_label_text(label: str) -> Optional[Tuple[str, int, str]]:
    """Parse labels like A:12:DA or _:12:X33 from RE_SCRIPT REMARK fields."""
    parts = label.strip().split(":")
    if len(parts) < 2:
        return None
    chain_id = parts[0]
    if chain_id == "_":
        chain_id = ""
    try:
        res_seq = int(parts[1])
    except ValueError:
        return None
    res_name = parts[2].strip() if len(parts) >= 3 else ""
    return (chain_id, res_seq, res_name)


def read_junction_residues_from_remarks(pdb_path: Path) -> List[Tuple[str, int, str]]:
    """Read final-output junction residue labels from RE_SCRIPT JUNCTION REMARKs."""
    residues: List[Tuple[str, int, str]] = []
    seen: Set[Tuple[str, int, str]] = set()
    with pdb_path.open("r") as fh:
        for line in fh:
            if not line.startswith(REMARK_PREFIX + " JUNCTION"):
                continue
            fields = parse_remark_key_values(line)
            residue_text = fields.get("residues")
            if not residue_text:
                continue
            for label_text in residue_text.split(","):
                parsed = parse_residue_label_text(label_text)
                if parsed is None:
                    continue
                key = (parsed[0], parsed[1], parsed[2].upper())
                if key in seen:
                    continue
                seen.add(key)
                residues.append(parsed)
    return residues


def phenix_name_selection(pdb_atom_name: str) -> str:
    """Return an exact four-character PDB atom-name selection."""
    return f'name "{pdb_atom_name}"'


def phenix_chain_residue_selection(chain_id: str, res_seq: int) -> str:
    if chain_id:
        return f"chain {chain_id} and resid {res_seq}"
    return f'chain " " and resid {res_seq}'


def residue_has_any_named_atom(atom_index: AtomIndex, chain_id: str, res_seq: int, pdb_atom_names: Tuple[str, ...]) -> bool:
    atoms = atom_index.get((chain_id.strip(), res_seq), {})
    for pdb_atom_name in pdb_atom_names:
        if normalize_atom_name(pdb_atom_name) in atoms:
            return True
    return False


def build_junction_selection_params_text(
    pdb_path: Path,
    junction_residues: List[Tuple[str, int, str]],
    link_distance_cutoff: Optional[float] = DEFAULT_LINK_DISTANCE_CUTOFF,
    linker_resname: str = DEFAULT_LINKER_RESNAME,
) -> Tuple[str, List[str]]:
    """Build a min_P_C5-style selection file with extra sugar atoms at junctions."""
    atom_index = read_pdb_atom_index(pdb_path)
    linker_resname = normalize_linker_resname(linker_resname)
    base_selection = " or ".join(phenix_name_selection(name) for name in JUNCTION_BASE_ATOM_NAMES)
    extra_sugar_selection = " or ".join(
        phenix_name_selection(name) for name in JUNCTION_EXTRA_SUGAR_ATOM_NAMES
    )

    included: List[Tuple[str, int, str]] = []
    seen_keys: Set[Tuple[str, int]] = set()
    for chain_id, res_seq, res_name in junction_residues:
        if (chain_id, res_seq) in seen_keys:
            continue
        seen_keys.add((chain_id, res_seq))
        if res_name.upper() == linker_resname:
            continue
        if not residue_has_any_named_atom(atom_index, chain_id, res_seq, JUNCTION_EXTRA_SUGAR_ATOM_NAMES):
            continue
        included.append((chain_id, res_seq, res_name))

    residue_clauses = [
        f"(({phenix_chain_residue_selection(chain_id, res_seq)}) and ({extra_sugar_selection}))"
        for chain_id, res_seq, _res_name in included
    ]

    selection = f"({base_selection})"
    if residue_clauses:
        selection += " or " + " or ".join(residue_clauses)

    labels_all = [
        f"{chain_id or '_'}:{res_seq}:{res_name or 'UNK'}" for chain_id, res_seq, res_name in junction_residues
    ]
    labels_included = [
        f"{chain_id or '_'}:{res_seq}:{res_name or 'UNK'}" for chain_id, res_seq, res_name in included
    ]

    lines: List[str] = []
    lines.append("# Junction-aware minimization selection generated from RE_SCRIPT JUNCTION REMARKs")
    lines.append(f"# source PDB: {pdb_path.name}")
    lines.append("# Base selection matches min_P_C5.params for all residues.")
    lines.append("# Extra sugar atoms C4', O4', C3', C2', and O2' are added only for junction residues.")
    lines.append("# parsed_junction_residues = " + (",".join(labels_all) if labels_all else "none"))
    lines.append("# sugar_expanded_residues = " + (",".join(labels_included) if labels_included else "none"))
    if link_distance_cutoff is not None:
        lines.append("pdb_interpretation {")
        lines.append(f"    link_distance_cutoff = {link_distance_cutoff:.3f}")
        lines.append("}")
    lines.append(f"selection = {selection}")
    lines.append("")

    reports = [
        f"Parsed {len(junction_residues)} junction residue label(s) from RE_SCRIPT REMARKs.",
        f"Added junction-specific sugar atoms for {len(included)} residue(s).",
    ]
    return "\n".join(lines), reports


def default_junction_params_path(pdb_path: Path) -> Path:
    return pdb_path.with_name(pdb_input_stem(pdb_path) + "_junctions.params")

# ---------------------------------------------------------------------------
# Params construction and command-line interface
# ---------------------------------------------------------------------------

def build_params_text(
    pdb_path: Path,
    links: List[LinkRecord],
    sigma: float = 0.02,
    angle_sigma: float = 3.0,
    nonbridging_angle: float = DEFAULT_NONBRIDGING_ANGLE,
    mixed_angle: float = DEFAULT_MIXED_ANGLE,
    bridging_angle: float = DEFAULT_BRIDGING_ANGLE,
    skip_phenix_builtin_angles: bool = True,
    x33_p_op_distance: float = X33_P_OP_DISTANCE,
    linker_resname: str = DEFAULT_LINKER_RESNAME,
    return_reports: bool = False,
):
    """
    Build the text of a geometry_restraints .params file from a list of links.
    """
    linker_resname = normalize_linker_resname(linker_resname)
    lines: List[str] = []
    has_x33_residues = bool(read_x33_residues(pdb_path, linker_resname=linker_resname))

    lines.append("# geometry_restraints generated from LINK records")
    lines.append(f"# source PDB: {pdb_path.name}")
    lines.append("# Each LINK is converted into a geometry_restraints.edits.bond")
    lines.append("# LINK-derived phosphate groups also receive geometry_restraints.edits.angle blocks")
    if skip_phenix_builtin_angles:
        lines.append("# Standard Phenix phosphate angles are omitted to avoid duplicate-angle errors")
    else:
        lines.append("# Diagnostic mode: standard Phenix phosphate angles are included")
    lines.append("# The numeric distance in the LINK record (if any) is ignored;")
    lines.append("# ideal distances are guessed from atom types (P-O, C-O, C-N, etc.).")
    lines.append("# Phosphate non-bridging oxygen names follow the input PDB: OP1/OP2 or O1P/O2P.")
    if has_x33_residues:
        lines.append(
            f"# Standalone {linker_resname} linker residues receive explicit internal P-OP1/P-OP2 bond and OP1-P-OP2 angle restraints."
        )
    lines.append("")
    lines.append("geometry_restraints {")
    lines.append("  edits {")

    for i, link in enumerate(links, start=1):
        ideal = guess_ideal_distance(link.atom_name_1, link.atom_name_2)

        # Optional comment showing original LINK line and distance (if present)
        lines.append(f"    # LINK {i}: {link.raw_line}")
        if link.original_distance is not None:
            lines.append(
                f"    # original LINK distance (ignored for restraints): {link.original_distance:.3f} A"
            )

        sel1 = (
            f"chain {link.chain_id_1} and resid {link.res_seq_1} "
            f"and name {atom_name(link.atom_name_1)}"
        )
        sel2 = (
            f"chain {link.chain_id_2} and resid {link.res_seq_2} "
            f"and name {atom_name(link.atom_name_2)}"
        )

        lines.append("    bond {")
        lines.append("      action = add")
        lines.append(f"      atom_selection_1 = {sel1}")
        lines.append(f"      atom_selection_2 = {sel2}")
        lines.append(f"      distance_ideal = {ideal:.3f}")
        lines.append(f"      sigma = {sigma:.3f}")
        lines.append("    }")
        lines.append("")

    x33_lines, x33_reports, x33_angle_keys = build_x33_internal_restraint_lines(
        pdb_path=pdb_path,
        sigma=sigma,
        angle_sigma=angle_sigma,
        p_op_distance=x33_p_op_distance,
        op_p_op_angle=nonbridging_angle,
        linker_resname=linker_resname,
    )
    lines.extend(x33_lines)

    angle_lines, angle_reports = build_phosphate_angle_lines(
        pdb_path=pdb_path,
        links=links,
        angle_sigma=angle_sigma,
        nonbridging_angle=nonbridging_angle,
        mixed_angle=mixed_angle,
        bridging_angle=bridging_angle,
        skip_phenix_builtin_angles=skip_phenix_builtin_angles,
        preexisting_angle_keys=x33_angle_keys,
        linker_resname=linker_resname,
    )
    lines.extend(angle_lines)

    lines.append("  }")
    lines.append("}")

    text = "\n".join(lines) + "\n"
    if return_reports:
        return text, x33_reports + angle_reports
    return text


GUIDANCE_TEXT = """Recommended Phenix workflow

*_links.params adds custom geometry restraints from LINK records, including nonstandard P-O bonds and phosphate-centered angle restraints.

*_junctions.params is the minimization atom-selection file. It replaces the older min_P_C5.params workflow for most re_helix outputs by keeping phosphate/backbone atoms movable globally, then adding extra sugar atoms only for residues listed in REMARK 950 RE_SCRIPT JUNCTION lines.

Recommended command:

phenix.geometry_minimization \\
  model_rex.pdb \\
  model_rex_links.params \\
  model_rex_junctions.params \\
  X33_phenix_atomtypes.cif \\
  X33_safe_interpretation.params

Use exactly one movement-selection file. Do not combine min_P_C5.params or min.params with *_junctions.params unless you deliberately want order-dependent selection behavior. If you disable *_junctions.params, use min_P_C5.params or a carefully prepared min.params instead.

The linker atomtypes CIF is needed when the standalone 3'-3' linker residue is nonstandard, such as X33. The safe interpretation params file should be listed last and avoids permissive automatic X33-P links. The default link_distance_cutoff here is 6.5 A.
"""


def default_output_base_path(pdb_path: Path) -> Path:
    return pdb_path.with_name(pdb_input_stem(pdb_path))


def output_path_from_base(output_base: Path, suffix: str) -> Path:
    return Path(str(output_base) + suffix)


def default_links_params_path(pdb_path: Path) -> Path:
    return output_path_from_base(default_output_base_path(pdb_path), "_links.params")


def ensure_parent_directory(path: Path) -> None:
    parent = path.expanduser().resolve().parent
    if not parent.exists():
        raise ValueError(f"Output directory does not exist: {parent}")


def write_text_file(path: Path, text: str) -> None:
    ensure_parent_directory(path)
    with path.open("w") as fh:
        fh.write(text)


def build_linker_atomtypes_cif_text(linker_resname: str) -> str:
    linker_resname = normalize_linker_resname(linker_resname)
    return f"""# {linker_resname}_phenix_atomtypes.cif
# Minimal Phenix/CCP4-style monomer definition for the standalone {linker_resname}
# 3'-3' phosphodiester linker phosphate used by DiLiuLab RE scripts.
#
# Use this file together with the *_links.params file generated by
# get_phenix_restraints.py. The params file supplies the internal P-OP1/P-OP2
# bond/angle restraints and LINK-derived external restraints.
#
# {linker_resname} atoms: P, OP1, OP2 only.

data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
{linker_resname} {linker_resname} 'standalone 3-prime-3-prime linker phosphate' non-polymer 3 3 .
#

data_comp_{linker_resname}
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.partial_charge
{linker_resname} P   P P   1.000
{linker_resname} OP1 O OP -0.500
{linker_resname} OP2 O OP -0.500
#
loop_
_chem_comp_tree.comp_id
_chem_comp_tree.atom_id
_chem_comp_tree.atom_back
_chem_comp_tree.atom_forward
_chem_comp_tree.connect_type
{linker_resname} P   n/a OP1 START
{linker_resname} OP1 P   .   END
{linker_resname} OP2 P   .   .
#
"""


def build_safe_interpretation_params_text(
    linker_resname: str,
    link_distance_cutoff: Optional[float] = DEFAULT_LINK_DISTANCE_CUTOFF,
) -> str:
    linker_resname = normalize_linker_resname(linker_resname)
    cutoff_text = ""
    if link_distance_cutoff is not None:
        cutoff_text = (
            "    # The Phenix default is 3.0 A; this re_helix helper defaults to 6.5 A.\n"
            "    # Avoid 7.0 A unless you deliberately need even more permissive automatic-link searching.\n"
            f"    link_distance_cutoff = {link_distance_cutoff:.3f}\n\n"
        )
    return f"""# {linker_resname}_safe_interpretation.params
#
# Use this file LAST in the phenix.geometry_minimization command, especially
# after any file that changes pdb_interpretation settings.
# It keeps normal polymer interpretation as intact as possible, but prevents
# permissive automatic ligand/small-molecule links involving {linker_resname}-P.
#
# Intended {linker_resname} external bonds should come from explicit *_links.params
# geometry_restraints.edits blocks, not from automatic distance-based linking.

pdb_interpretation {{
{cutoff_text}    automatic_linking {{
        link_all = False
        link_amino_acid_rna_dna = False
        link_ligands = False
        link_small_molecules = False
    }}

    # {linker_resname} is a 3'-3' linker phosphate. Its intended external links are
    # P--O3' to the two adjacent junction residues. Exclude automatic {linker_resname}-P
    # links to O5' or other nearby atoms; explicit *_links.params restraints are not
    # removed by this exclusion.
    exclude_from_automatic_linking {{
        selection_1 = resname {linker_resname} and name " P "
        selection_2 = not (name " O3'" or name " O3*")
    }}
    exclude_from_automatic_linking {{
        selection_1 = not (name " O3'" or name " O3*")
        selection_2 = resname {linker_resname} and name " P "
    }}
}}
"""


def write_linker_support_files(
    output_dir: Path,
    linker_resname: str,
    link_distance_cutoff: Optional[float],
) -> Tuple[List[Path], List[str]]:
    linker_resname = normalize_linker_resname(linker_resname)
    if linker_resname in STANDARD_NUCLEIC_ACID_RESNAMES:
        return [], [
            f"Skipped linker support files for standard residue name {linker_resname}; "
            "Phenix should already know that monomer definition."
        ]

    atomtypes_path = output_dir / f"{linker_resname}_phenix_atomtypes.cif"
    safe_path = output_dir / f"{linker_resname}_safe_interpretation.params"
    write_text_file(atomtypes_path, build_linker_atomtypes_cif_text(linker_resname))
    write_text_file(safe_path, build_safe_interpretation_params_text(linker_resname, link_distance_cutoff))
    return [atomtypes_path, safe_path], [
        f"Wrote linker atomtypes CIF to: {atomtypes_path}",
        f"Wrote linker safe interpretation params to: {safe_path}",
    ]


def format_phenix_command(parts: List[str]) -> str:
    if not parts:
        return ""
    lines = []
    for idx, part in enumerate(parts):
        suffix = " \\" if idx < len(parts) - 1 else ""
        lines.append("  " + shlex.quote(str(part)) + suffix)
    return "\n".join(lines)


def build_recommended_phenix_command(
    pdb_path: Path,
    links_out: Optional[Path],
    junctions_out: Optional[Path],
    support_files: List[Path],
) -> str:
    parts: List[str] = ["phenix.geometry_minimization", str(pdb_path)]
    if links_out is not None:
        parts.append(str(links_out))
    if junctions_out is not None:
        parts.append(str(junctions_out))
    parts.extend(str(path) for path in support_files)
    return format_phenix_command(parts)


def generate_phenix_restraint_outputs(args) -> Dict[str, object]:
    pdb_path = Path(args.pdb_in).expanduser()
    if not pdb_path.is_file():
        raise ValueError(f"Input PDB file not found: {pdb_path}")

    linker_resname = normalize_linker_resname(getattr(args, "linker_resname", DEFAULT_LINKER_RESNAME))
    output_base_value = getattr(args, "output_base", None)
    output_base = Path(output_base_value).expanduser() if output_base_value else default_output_base_path(pdb_path)
    params_out_value = getattr(args, "params_out", None)
    junctions_out_value = getattr(args, "junctions_out", None)

    params_out = Path(params_out_value).expanduser() if params_out_value else output_path_from_base(output_base, "_links.params")
    junctions_out = (
        Path(junctions_out_value).expanduser()
        if junctions_out_value
        else output_path_from_base(output_base, "_junctions.params")
    )
    support_dir_value = getattr(args, "support_out_dir", None)
    support_out_dir = Path(support_dir_value).expanduser() if support_dir_value else output_base.parent

    generate_junctions = bool(
        getattr(args, "generate_junctions_params", not getattr(args, "no_junctions_params", False))
    )
    generate_support = bool(
        getattr(args, "generate_linker_support_files", not getattr(args, "no_linker_support_files", False))
    )
    link_distance_cutoff = getattr(args, "link_distance_cutoff", DEFAULT_LINK_DISTANCE_CUTOFF)

    links = read_link_records(pdb_path)
    linker_residues = read_x33_residues(pdb_path, linker_resname=linker_resname)
    junction_residues = read_junction_residues_from_remarks(pdb_path)

    reports: List[str] = [
        f"Found {len(links)} LINK record(s).",
        f"Found {len(linker_residues)} standalone {linker_resname} linker residue(s).",
        f"Found {len(junction_residues)} RE_SCRIPT JUNCTION residue label(s).",
    ]
    links_written: Optional[Path] = None
    junctions_written: Optional[Path] = None
    support_files: List[Path] = []

    if links or linker_residues:
        params_text, angle_reports = build_params_text(
            pdb_path,
            links,
            sigma=args.sigma,
            angle_sigma=args.angle_sigma,
            nonbridging_angle=args.angle_nonbridging,
            mixed_angle=args.angle_mixed,
            bridging_angle=args.angle_bridging,
            skip_phenix_builtin_angles=not args.include_phenix_builtin_angles,
            x33_p_op_distance=args.x33_p_op_distance,
            linker_resname=linker_resname,
            return_reports=True,
        )
        write_text_file(params_out, params_text)
        links_written = params_out
        reports.extend(angle_reports)
        reports.append(f"Wrote geometry_restraints parameters to: {params_out}")
    else:
        reports.append(f"No LINK records or {linker_resname} residues found; no _links.params file was written.")

    if junction_residues and generate_junctions:
        junction_text, junction_reports = build_junction_selection_params_text(
            pdb_path,
            junction_residues,
            link_distance_cutoff=link_distance_cutoff,
            linker_resname=linker_resname,
        )
        write_text_file(junctions_out, junction_text)
        junctions_written = junctions_out
        reports.extend(junction_reports)
        reports.append(f"Wrote junction-aware minimization selection parameters to: {junctions_out}")
    elif junction_residues and not generate_junctions:
        reports.append(
            "Skipped _junctions.params generation. Use exactly one movement-selection file instead, "
            "for example min_P_C5.params or min.params."
        )
    else:
        reports.append("No RE_SCRIPT JUNCTION REMARKs found; no _junctions.params file was written.")

    if generate_support:
        written_support, support_reports = write_linker_support_files(
            support_out_dir,
            linker_resname=linker_resname,
            link_distance_cutoff=link_distance_cutoff,
        )
        support_files.extend(written_support)
        reports.extend(support_reports)
    else:
        reports.append("Skipped linker support CIF/safe-interpretation file generation.")

    if junctions_written is not None:
        reports.append(
            "Use *_junctions.params as the movement-selection file. Do not also include min_P_C5.params or min.params "
            "unless you deliberately want to override or test selection behavior."
        )
    else:
        reports.append(
            "Because no *_junctions.params file was generated, add exactly one movement-selection file such as "
            "min_P_C5.params or a carefully prepared min.params when running phenix.geometry_minimization."
        )

    phenix_command = build_recommended_phenix_command(
        pdb_path=pdb_path,
        links_out=links_written,
        junctions_out=junctions_written,
        support_files=support_files,
    )

    return {
        "pdb_path": pdb_path,
        "linker_resname": linker_resname,
        "links_out": links_written,
        "junctions_out": junctions_written,
        "support_files": support_files,
        "reports": reports,
        "phenix_command": phenix_command,
    }


def format_generation_summary(result: Dict[str, object]) -> str:
    reports = list(result.get("reports", []))
    command = str(result.get("phenix_command", "") or "")
    lines = ["Generated Phenix restraint helper output:"]
    lines.extend(f"- {msg}" for msg in reports)
    if command:
        lines.append("")
        lines.append("Suggested Phenix command:")
        lines.append(command)
    return "\n".join(lines)


def build_cli_command(program: str, args) -> str:
    parts = [
        sys.executable or "python3",
        program,
        args.pdb_in,
        "--linker-resname",
        args.linker_resname,
        "--link-distance-cutoff",
        str(args.link_distance_cutoff),
    ]
    if getattr(args, "output_base", None):
        parts.extend(["--output-base", args.output_base])
    if getattr(args, "params_out", None):
        parts.extend(["--out", args.params_out])
    if getattr(args, "junctions_out", None):
        parts.extend(["--junctions-out", args.junctions_out])
    if getattr(args, "support_out_dir", None):
        parts.extend(["--support-out-dir", args.support_out_dir])
    if not args.generate_junctions_params:
        parts.append("--no-junctions-params")
    if not args.generate_linker_support_files:
        parts.append("--no-linker-support-files")
    if args.include_phenix_builtin_angles:
        parts.append("--include-phenix-builtin-angles")
    if args.sigma != 0.02:
        parts.extend(["--sigma", str(args.sigma)])
    if args.angle_sigma != 3.0:
        parts.extend(["--angle-sigma", str(args.angle_sigma)])
    if args.angle_nonbridging != DEFAULT_NONBRIDGING_ANGLE:
        parts.extend(["--angle-nonbridging", str(args.angle_nonbridging)])
    if args.angle_mixed != DEFAULT_MIXED_ANGLE:
        parts.extend(["--angle-mixed", str(args.angle_mixed)])
    if args.angle_bridging != DEFAULT_BRIDGING_ANGLE:
        parts.extend(["--angle-bridging", str(args.angle_bridging)])
    if args.x33_p_op_distance != X33_P_OP_DISTANCE:
        parts.extend(["--linker-p-op-distance", str(args.x33_p_op_distance)])
    return " ".join(shlex.quote(str(part)) for part in parts)


def run_gui(initial_path=None) -> int:
    if tk is None:
        sys.stderr.write("ERROR: Tkinter is not available in this Python installation.\n")
        return 1

    root = tk.Tk()
    root.title(f"{TOOL_NAME} {VERSION}")
    root.geometry("980x760")
    root.minsize(900, 680)

    input_var = tk.StringVar(value=initial_path or "")
    output_base_var = tk.StringVar()
    if initial_path:
        output_base_var.set(str(default_output_base_path(Path(initial_path))))
    output_base_state = {"last_default": output_base_var.get(), "custom": False}

    linker_var = tk.StringVar(value=DEFAULT_LINKER_RESNAME)
    cutoff_var = tk.StringVar(value=f"{DEFAULT_LINK_DISTANCE_CUTOFF:.1f}")
    generate_junctions_var = tk.BooleanVar(value=True)
    generate_support_var = tk.BooleanVar(value=True)
    sigma_var = tk.StringVar(value="0.02")
    angle_sigma_var = tk.StringVar(value="3.0")
    p_op_distance_var = tk.StringVar(value=f"{X33_P_OP_DISTANCE:.3f}")
    include_builtin_angles_var = tk.BooleanVar(value=False)

    def derived_output_base(path_text: str) -> str:
        return str(default_output_base_path(Path(path_text))) if path_text else ""

    def set_input_path(path: str) -> None:
        input_var.set(path)
        new_default = derived_output_base(path)
        current_output = output_base_var.get().strip()
        if not output_base_state["custom"] or current_output == output_base_state["last_default"]:
            output_base_var.set(new_default)
            output_base_state["custom"] = False
        output_base_state["last_default"] = new_default

    def mark_output_custom(*_args) -> None:
        current_output = output_base_var.get().strip()
        output_base_state["custom"] = bool(current_output and current_output != output_base_state["last_default"])

    output_base_var.trace_add("write", mark_output_custom)

    def choose_input() -> None:
        path = filedialog.askopenfilename(
            title="Choose input PDB",
            filetypes=[("PDB files", "*.pdb *.ent *.pdb.txt"), ("All files", "*")],
        )
        if path:
            set_input_path(path)

    def choose_output_base() -> None:
        initial_dir = ""
        if input_var.get().strip():
            initial_dir = str(Path(input_var.get().strip()).expanduser().parent)
        path = filedialog.asksaveasfilename(
            title="Choose output base",
            initialdir=initial_dir or None,
            initialfile=Path(output_base_var.get().strip()).name if output_base_var.get().strip() else "",
            filetypes=[("Output base", "*"), ("All files", "*")],
        )
        if path:
            output_base_var.set(path)

    root.columnconfigure(1, weight=1)

    header = tk.Label(root, text=f"{TOOL_NAME} {VERSION}", font=("TkDefaultFont", 15, "bold"))
    header.grid(row=0, column=0, columnspan=3, sticky="w", padx=10, pady=(10, 6))

    tk.Label(root, text="Input PDB:", anchor="e").grid(row=1, column=0, sticky="e", padx=8, pady=4)
    tk.Entry(root, textvariable=input_var, width=68).grid(row=1, column=1, sticky="we", padx=8, pady=4)
    tk.Button(root, text="Browse", command=choose_input).grid(row=1, column=2, sticky="we", padx=8, pady=4)

    tk.Label(root, text="Output base:", anchor="e").grid(row=2, column=0, sticky="e", padx=8, pady=4)
    tk.Entry(root, textvariable=output_base_var, width=68).grid(row=2, column=1, sticky="we", padx=8, pady=4)
    tk.Button(root, text="Browse", command=choose_output_base).grid(row=2, column=2, sticky="we", padx=8, pady=4)

    tk.Label(root, text="3'-3' linker resName:", anchor="e").grid(row=3, column=0, sticky="e", padx=8, pady=4)
    tk.Entry(root, textvariable=linker_var, width=10).grid(row=3, column=1, sticky="w", padx=8, pady=4)

    tk.Label(root, text="link_distance_cutoff:", anchor="e").grid(row=4, column=0, sticky="e", padx=8, pady=4)
    tk.Entry(root, textvariable=cutoff_var, width=10).grid(row=4, column=1, sticky="w", padx=8, pady=4)

    tk.Checkbutton(root, text="Generate *_junctions.params", variable=generate_junctions_var).grid(
        row=5, column=1, sticky="w", padx=8, pady=(8, 2)
    )
    tk.Checkbutton(root, text="Generate linker atomtypes CIF and safe interpretation params", variable=generate_support_var).grid(
        row=6, column=1, sticky="w", padx=8, pady=2
    )

    advanced_visible = tk.BooleanVar(value=False)

    def toggle_advanced() -> None:
        if advanced_visible.get():
            advanced.grid_remove()
            advanced_visible.set(False)
            advanced_button.configure(text="Show Advanced restraint values")
        else:
            advanced.grid()
            advanced_visible.set(True)
            advanced_button.configure(text="Hide Advanced restraint values")

    advanced_button = tk.Button(root, text="Show Advanced restraint values", command=toggle_advanced)
    advanced_button.grid(row=7, column=1, sticky="w", padx=8, pady=(8, 4))

    advanced = tk.LabelFrame(root, text="Advanced restraint values", padx=8, pady=6)
    advanced.grid(row=8, column=0, columnspan=3, sticky="we", padx=8, pady=(2, 4))
    advanced.columnconfigure(1, weight=1)
    tk.Label(advanced, text="Bond sigma:").grid(row=0, column=0, sticky="e", padx=6, pady=3)
    tk.Entry(advanced, textvariable=sigma_var, width=10).grid(row=0, column=1, sticky="w", padx=6, pady=3)
    tk.Label(advanced, text="Angle sigma:").grid(row=1, column=0, sticky="e", padx=6, pady=3)
    tk.Entry(advanced, textvariable=angle_sigma_var, width=10).grid(row=1, column=1, sticky="w", padx=6, pady=3)
    tk.Label(advanced, text="Linker P-OP distance:").grid(row=2, column=0, sticky="e", padx=6, pady=3)
    tk.Entry(advanced, textvariable=p_op_distance_var, width=10).grid(row=2, column=1, sticky="w", padx=6, pady=3)
    tk.Checkbutton(
        advanced,
        text="Include phosphate angles Phenix normally supplies",
        variable=include_builtin_angles_var,
    ).grid(row=3, column=1, sticky="w", padx=6, pady=3)
    advanced.grid_remove()

    note_widget = scrolledtext.ScrolledText(root, wrap="word", height=12)
    note_widget.grid(row=9, column=0, columnspan=3, sticky="nsew", padx=8, pady=(8, 4))
    note_widget.insert("1.0", GUIDANCE_TEXT)
    note_widget.configure(state="disabled")

    log_widget = scrolledtext.ScrolledText(root, wrap="word", height=8)
    log_widget.grid(row=10, column=0, columnspan=3, sticky="nsew", padx=8, pady=(4, 8))
    root.rowconfigure(9, weight=1)
    root.rowconfigure(10, weight=1)

    def append_log(text: str) -> None:
        log_widget.insert("end", text)
        log_widget.see("end")
        log_widget.update_idletasks()

    def run_processing() -> None:
        try:
            input_path = input_var.get().strip()
            if not input_path:
                raise ValueError("Please choose an input PDB file.")
            output_base = output_base_var.get().strip()
            if not output_base:
                output_base = derived_output_base(input_path)
                output_base_var.set(output_base)

            linker_resname = normalize_linker_resname(linker_var.get())
            linker_var.set(linker_resname)
            link_distance_cutoff = float(cutoff_var.get())
            if link_distance_cutoff <= 0:
                raise ValueError("link_distance_cutoff must be positive.")

            gui_args = argparse.Namespace(
                pdb_in=input_path,
                output_base=output_base,
                params_out=None,
                junctions_out=None,
                support_out_dir=None,
                linker_resname=linker_resname,
                link_distance_cutoff=link_distance_cutoff,
                generate_junctions_params=generate_junctions_var.get(),
                generate_linker_support_files=generate_support_var.get(),
                sigma=float(sigma_var.get()),
                angle_sigma=float(angle_sigma_var.get()),
                angle_nonbridging=DEFAULT_NONBRIDGING_ANGLE,
                angle_mixed=DEFAULT_MIXED_ANGLE,
                angle_bridging=DEFAULT_BRIDGING_ANGLE,
                include_phenix_builtin_angles=include_builtin_angles_var.get(),
                x33_p_op_distance=float(p_op_distance_var.get()),
            )

            cli_cmd = build_cli_command(os.path.basename(__file__), gui_args)
            print("Equivalent CLI command:", flush=True)
            print(cli_cmd, flush=True)
            append_log("Equivalent CLI command:\n" + cli_cmd + "\n\n")
            result = generate_phenix_restraint_outputs(gui_args)
            summary = format_generation_summary(result)
            print(summary, flush=True)
            append_log(summary + "\n\n")
            messagebox.showinfo(
                "Done",
                "Generated Phenix restraint files.\n\n"
                "See the run log for output paths and the suggested phenix.geometry_minimization command.",
            )
        except Exception as exc:
            print(f"Error: {exc}", flush=True)
            append_log(f"Error: {exc}\n\n")
            messagebox.showerror("Error", str(exc))

    button_row = tk.Frame(root)
    button_row.grid(row=11, column=0, columnspan=3, sticky="we", padx=8, pady=(0, 12))
    button_row.columnconfigure(0, weight=1)
    tk.Button(button_row, text="Generate Phenix restraints", command=run_processing, height=2).grid(
        row=0, column=0, sticky="we", padx=(0, 6)
    )
    tk.Button(button_row, text="Close", command=root.destroy, height=2).grid(row=0, column=1, sticky="e")

    root.mainloop()
    return 0


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Convert LINK records in a PDB file into Phenix geometry_restraints "
            "bond edits plus phosphate angle edits (.params), including standalone linker "
            "restraints and optional junction-selection params, suitable for "
            "phenix.geometry_minimization."
        )
    )
    parser.add_argument("pdb_in", help="Input PDB file containing LINK records.")
    parser.add_argument("--gui", action="store_true", help="Open the Tk GUI.")
    parser.add_argument("-v", "--version", action="version", version=f"{TOOL_NAME} {VERSION}")
    parser.add_argument(
        "--output-base",
        default=None,
        help="Output prefix. Default: <input_stem>; writes <base>_links.params and <base>_junctions.params.",
    )
    parser.add_argument(
        "-o",
        "--out",
        dest="params_out",
        default=None,
        help="Output _links.params filename (default: <base>_links.params).",
    )
    parser.add_argument(
        "--linker-resname",
        default=DEFAULT_LINKER_RESNAME,
        help=f"Residue name for standalone 3'-3' linker phosphate residues (default: {DEFAULT_LINKER_RESNAME}).",
    )
    parser.add_argument(
        "--link-distance-cutoff",
        type=float,
        default=DEFAULT_LINK_DISTANCE_CUTOFF,
        help=f"pdb_interpretation.link_distance_cutoff for generated params files (default: {DEFAULT_LINK_DISTANCE_CUTOFF:.1f}).",
    )
    parser.add_argument(
        "--sigma",
        type=float,
        default=0.02,
        help="Sigma for all bond restraints (default: 0.02 A).",
    )
    parser.add_argument(
        "--angle-sigma",
        type=float,
        default=3.0,
        help="Sigma for phosphate angle restraints in degrees (default: 3.0).",
    )
    parser.add_argument(
        "--angle-nonbridging",
        type=float,
        default=DEFAULT_NONBRIDGING_ANGLE,
        help=f"Ideal angle for OP1/O1P-P-OP2/O2P in degrees (default: {DEFAULT_NONBRIDGING_ANGLE:.1f}).",
    )
    parser.add_argument(
        "--angle-mixed",
        type=float,
        default=DEFAULT_MIXED_ANGLE,
        help=f"Ideal angle for non-bridging O-P-bridging O in degrees (default: {DEFAULT_MIXED_ANGLE:.1f}).",
    )
    parser.add_argument(
        "--angle-bridging",
        type=float,
        default=DEFAULT_BRIDGING_ANGLE,
        help=f"Ideal angle for bridging O-P-bridging O in degrees (default: {DEFAULT_BRIDGING_ANGLE:.1f}).",
    )
    parser.add_argument(
        "--include-phenix-builtin-angles",
        action="store_true",
        help=(
            "Also write phosphate angles that Phenix normally already has. "
            "This reproduces V2-style output for diagnostics but can trigger duplicate-angle errors."
        ),
    )
    parser.add_argument(
        "--x33-p-op-distance",
        "--linker-p-op-distance",
        type=float,
        dest="x33_p_op_distance",
        default=X33_P_OP_DISTANCE,
        help=f"Ideal linker P-OP1/P-OP2 bond distance in A (default: {X33_P_OP_DISTANCE:.3f}).",
    )
    parser.add_argument(
        "--junctions-out",
        default=None,
        help="Output filename for junction-aware selection params (default: <pdb_stem>_junctions.params).",
    )
    parser.add_argument(
        "--no-junctions-params",
        action="store_true",
        help="Do not write the optional <pdb_stem>_junctions.params file even if JUNCTION REMARKs are present.",
    )
    parser.add_argument(
        "--support-out-dir",
        default=None,
        help="Directory for <resname>_phenix_atomtypes.cif and <resname>_safe_interpretation.params.",
    )
    parser.add_argument(
        "--no-linker-support-files",
        action="store_true",
        help="Do not write linker atomtypes CIF and safe interpretation params files.",
    )
    return parser


def main(argv=None) -> int:
    if argv is None:
        argv = sys.argv[1:]
    if len(argv) == 0:
        return run_gui()
    if "--gui" in argv:
        remaining = [item for item in argv if item != "--gui"]
        initial_path = remaining[0] if remaining else None
        return run_gui(initial_path=initial_path)

    parser = build_arg_parser()
    args = parser.parse_args(argv)
    args.generate_junctions_params = not args.no_junctions_params
    args.generate_linker_support_files = not args.no_linker_support_files

    try:
        result = generate_phenix_restraint_outputs(args)
    except Exception as exc:
        sys.stderr.write(f"ERROR: {exc}\n")
        return 1

    print("Equivalent CLI command:")
    command_args = argparse.Namespace(**vars(args))
    if command_args.output_base is None:
        command_args.output_base = str(default_output_base_path(Path(command_args.pdb_in)))
    print(build_cli_command(os.path.basename(__file__), command_args))
    print(format_generation_summary(result))
    return 0


if __name__ == "__main__":
    sys.exit(main())
