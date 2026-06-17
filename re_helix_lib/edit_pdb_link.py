#!/usr/bin/env python3
"""
edit_pdb_link.py

Helpers for parsing and writing LINK records in PDB files.

This module is designed to complement `edit_pdb_atom.py` without
modifying it.  It provides:

  - class pdb_link_record
  - function file2rec_link(inpfile, atom_rec_list, link_rec_list)
  - function rec2file_link(atom_rec_list, link_rec_list, outfile, reorder_serial=False)
  - function have_link(pdb_path)
  - utilities:
        build_atom_index(atom_rec_list)
        find_atom(atom_index, chainID, resSeq, atom_name)
        update_link_distances_from_atoms(link_rec_list, atom_rec_list)

Typical usage
-------------

    from edit_pdb_atom import file2rec, rec2file
    from edit_pdb_link import (
        pdb_link_record,
        file2rec_link,
        rec2file_link,
        have_link,
        build_atom_index,
        update_link_distances_from_atoms,
    )

    atom_recs = []
    link_recs = []

    with open("input.pdb", "r") as fin:
        file2rec_link(fin, atom_recs, link_recs)

    # ... manipulate atoms and/or links ...

    with open("output.pdb", "w") as fout:
        rec2file_link(atom_recs, link_recs, fout, reorder_serial=False)

Notes
-----
- ATOM / HETATM / TER records are parsed with the *same logic* as in
  edit_pdb_atom.file2rec, so behaviour is consistent.
- LINK records are parsed in a token-based way and stored as
  `pdb_link_record` objects.  Updating any field on a link rebuilds
  the underlying PDB LINK line in a canonical format similar to the
  bowtie/merge scripts you already use.
"""

import math
from typing import Dict, List, Optional, Tuple, Iterable, Set, TextIO

try:
    from . import edit_pdb_atom
except ImportError:  # pragma: no cover - keeps direct script execution working.
    import edit_pdb_atom


# ---------------------------------------------------------------------------
# LINK record class
# ---------------------------------------------------------------------------

class pdb_link_record:
    """
    Representation of a PDB LINK record.

    Fields
    ------
    string : str
        The current PDB line for this LINK (including newline).
        This is rebuilt whenever one of the fields is updated.
    recordName : str
        Always 'LINK'.
    # Endpoint 1
    name1 : str
    resName1 : str
    chainID1 : str
    resSeq1 : int
    iCode1 : str   (insertion code, usually ' ')
    # Endpoint 2
    name2 : str
    resName2 : str
    chainID2 : str
    resSeq2 : int
    iCode2 : str
    # Symmetry / distance
    sym1 : str     (e.g. '1555')
    sym2 : str     (e.g. '1555')
    distance : Optional[float]
        Ideal/observed distance in Å, if present; otherwise None.

    Notes
    -----
    The constructor parses the given LINK line (using whitespace splitting)
    to populate these fields, but the output string is always written in a
    canonical format similar to the format used by your bowtie-link scripts:

        LINK        <name1> <res1> <chain1><resSeq1>                \
                    <name2> <res2> <chain2><resSeq2>     sym1   sym2 dist

    where sym1/sym2 default to '1555' and dist is '     ' if distance is None.
    """

    def __init__(self, string: str) -> None:
        # Store the raw input (mainly for debugging)
        raw = string.rstrip("\n")
        self.recordName = "LINK"

        # Defaults (in case parsing fails partially)
        self.name1 = ""
        self.resName1 = ""
        self.chainID1 = ""
        self.resSeq1 = 0
        self.iCode1 = " "

        self.name2 = ""
        self.resName2 = ""
        self.chainID2 = ""
        self.resSeq2 = 0
        self.iCode2 = " "

        self.sym1 = "1555"
        self.sym2 = "1555"
        self.distance: Optional[float] = None

        # Parse using a token-based approach; this works for standard PDB LINK
        # lines as well as the ones written by MergeChain_* scripts.
        tokens = raw.split()
        if len(tokens) < 9 or tokens[0] != "LINK":
            raise ValueError(f"Not a valid LINK line: {string!r}")

        try:
            # LINK atom1 resName1 chain1 resSeq1 atom2 resName2 chain2 resSeq2 ...
            self.name1 = tokens[1]
            self.resName1 = tokens[2]
            self.chainID1 = tokens[3]
            self.resSeq1 = int(tokens[4])

            self.name2 = tokens[5]
            self.resName2 = tokens[6]
            self.chainID2 = tokens[7]
            self.resSeq2 = int(tokens[8])
        except Exception as e:
            raise ValueError(f"Could not parse LINK fields from: {string!r}") from e

        # Optional symmetry and distance
        if len(tokens) >= 10:
            self.sym1 = tokens[9]
        if len(tokens) >= 11:
            self.sym2 = tokens[10]
        if len(tokens) >= 12:
            try:
                self.distance = float(tokens[11])
            except ValueError:
                self.distance = None

        # Build the canonical PDB line
        self._rebuild_string()

    # ----------------------------
    # Internal: rebuild PDB string
    # ----------------------------

    def _rebuild_string(self) -> None:
        """
        Rebuild self.string from the current fields in a canonical format.

        This does *not* try to preserve the original spacing; instead we
        consistently use the same layout, similar to your existing scripts.
        """
        dist_str = "     "
        if self.distance is not None:
            dist_str = f"{self.distance:5.2f}"

        line = (
            "LINK        "
            f"{self.name1:>4s} {self.resName1:>3s} {self.chainID1:1s}{self.resSeq1:4d}"
            "                "
            f"{self.name2:>4s} {self.resName2:>3s} {self.chainID2:1s}{self.resSeq2:4d}"
            f"     {self.sym1:>4s}   {self.sym2:>4s} {dist_str}"
        )
        # Ensure newline at end
        self.string = line + "\n"

    # ----------------------------
    # Endpoint-1 update helpers
    # ----------------------------

    def update_name1(self, new_name1: str) -> None:
        self.name1 = new_name1
        self._rebuild_string()

    def update_resName1(self, new_resName1: str) -> None:
        self.resName1 = new_resName1
        self._rebuild_string()

    def update_chainID1(self, new_chainID1: str) -> None:
        self.chainID1 = new_chainID1
        self._rebuild_string()

    def update_resSeq1(self, new_resSeq1: int) -> None:
        self.resSeq1 = new_resSeq1
        self._rebuild_string()

    def update_iCode1(self, new_iCode1: str = " ") -> None:
        # Not printed explicitly in our canonical format, but kept for future use.
        self.iCode1 = (new_iCode1 or " ")[0]
        self._rebuild_string()

    # ----------------------------
    # Endpoint-2 update helpers
    # ----------------------------

    def update_name2(self, new_name2: str) -> None:
        self.name2 = new_name2
        self._rebuild_string()

    def update_resName2(self, new_resName2: str) -> None:
        self.resName2 = new_resName2
        self._rebuild_string()

    def update_chainID2(self, new_chainID2: str) -> None:
        self.chainID2 = new_chainID2
        self._rebuild_string()

    def update_resSeq2(self, new_resSeq2: int) -> None:
        self.resSeq2 = new_resSeq2
        self._rebuild_string()

    def update_iCode2(self, new_iCode2: str = " ") -> None:
        self.iCode2 = (new_iCode2 or " ")[0]
        self._rebuild_string()

    # ----------------------------
    # Symmetry / distance
    # ----------------------------

    def update_sym1(self, new_sym1: str) -> None:
        self.sym1 = new_sym1.strip() or "1555"
        self._rebuild_string()

    def update_sym2(self, new_sym2: str) -> None:
        self.sym2 = new_sym2.strip() or "1555"
        self._rebuild_string()

    def update_distance(self, new_distance: Optional[float]) -> None:
        self.distance = new_distance
        self._rebuild_string()

    def update_from_atoms(
        self,
        atom1: "edit_pdb_atom.pdb_atom_record",
        atom2: "edit_pdb_atom.pdb_atom_record",
        recompute_distance: bool = True,
    ) -> None:
        """
        Update both endpoints from two atom records.

        If recompute_distance is True, distance is set from the atom coords.
        """
        # Endpoint 1
        self.name1 = atom1.name.strip()
        self.resName1 = atom1.resName.strip()
        self.chainID1 = atom1.chainID
        self.resSeq1 = atom1.resSeq

        # Endpoint 2
        self.name2 = atom2.name.strip()
        self.resName2 = atom2.resName.strip()
        self.chainID2 = atom2.chainID
        self.resSeq2 = atom2.resSeq

        # Distance
        if recompute_distance:
            dx = atom1.x - atom2.x
            dy = atom1.y - atom2.y
            dz = atom1.z - atom2.z
            self.distance = math.sqrt(dx * dx + dy * dy + dz * dz)

        self._rebuild_string()

    # ----------------------------
    # Convenience / query helpers
    # ----------------------------

    def endpoints(self) -> Tuple[Tuple[str, str, int, str], Tuple[str, str, int, str]]:
        """
        Return endpoints as:
            ((chainID1, resName1, resSeq1, name1),
             (chainID2, resName2, resSeq2, name2))
        """
        return (
            (self.chainID1, self.resName1, self.resSeq1, self.name1),
            (self.chainID2, self.resName2, self.resSeq2, self.name2),
        )

    def involves_chain(self, chain_id: str) -> bool:
        """Return True if either endpoint uses the given chain ID."""
        return self.chainID1 == chain_id or self.chainID2 == chain_id

    def involves_residue(self, chain_id: str, res_seq: int) -> bool:
        """Return True if either endpoint is the given (chain, resSeq)."""
        return (
            (self.chainID1 == chain_id and self.resSeq1 == res_seq)
            or (self.chainID2 == chain_id and self.resSeq2 == res_seq)
        )


# ---------------------------------------------------------------------------
# IO helpers: file2rec_link, rec2file_link, have_link
# ---------------------------------------------------------------------------

def file2rec_link(
    inpfile: TextIO,
    atom_rec_list: List[edit_pdb_atom.pdb_atom_record],
    link_rec_list: List[pdb_link_record],
) -> None:
    """
    Read a PDB file into separate ATOM/HETATM/TER and LINK lists.

    This mirrors edit_pdb_atom.file2rec for atoms/TERs, and additionally
    collects LINK records as pdb_link_record objects.

    Parameters
    ----------
    inpfile : file-like object opened for reading text
    atom_rec_list : list
        Will be filled with pdb_atom_record / pdb_ter_record instances.
    link_rec_list : list
        Will be filled with pdb_link_record instances for each LINK line.
    """
    last_serial = 0
    last_resSeq = 0
    added_ter_serial = 0  # recording added serial for TER having no serial

    for eachline in inpfile:
        if eachline[:4] == "ATOM" or eachline[:6] == "HETATM":
            curr_atom = edit_pdb_atom.pdb_atom_record(eachline)
            curr_atom_serial = curr_atom.serial + added_ter_serial
            curr_atom.update_serial(curr_atom_serial)
            last_serial = curr_atom.serial
            last_resSeq = curr_atom.resSeq
            atom_rec_list.append(curr_atom)
        elif eachline[:3] == "TER":
            curr_atom = edit_pdb_atom.pdb_ter_record(eachline)
            if curr_atom.serial is None:
                curr_atom.update_serial(last_serial + 1)
                added_ter_serial = added_ter_serial + 1
            if curr_atom.resSeq is None:
                curr_atom.update_resSeq(last_resSeq)
            atom_rec_list.append(curr_atom)
        elif eachline.startswith("LINK"):
            link_rec = pdb_link_record(eachline)
            link_rec_list.append(link_rec)
        else:
            # other record types (REMARK, CONECT, END, etc.) are ignored here
            continue


def rec2file_link(
    atom_rec_list: List[edit_pdb_atom.pdb_atom_record],
    link_rec_list: List[pdb_link_record],
    outfile: TextIO,
    reorder_serial: bool = False,
) -> None:
    """
    Write LINK records, then ATOM/HETATM/TER records to a PDB file.

    Parameters
    ----------
    atom_rec_list : list
        pdb_atom_record / pdb_ter_record instances.
    link_rec_list : list
        pdb_link_record instances (their .string field is written directly).
    outfile : file-like object opened for writing text
    reorder_serial : bool
        If True, renumber ATOM/HETATM/TER serials starting from 1.
        LINK lines are not affected by this flag.
    """
    # 1) LINK records first (in the order given)
    for link_rec in link_rec_list:
        outfile.write(link_rec.string)

    # 2) Then atoms/TERs using the original helper
    edit_pdb_atom.rec2file(atom_rec_list, outfile, reorder_serial=reorder_serial)


def have_link(pdb_path: str) -> bool:
    """
    Return True if the given PDB file contains at least one LINK record.

    Parameters
    ----------
    pdb_path : str
        Path to a PDB file.

    Returns
    -------
    bool
        True if any line in the file starts with 'LINK', False otherwise.
    """
    with open(pdb_path, "r") as fh:
        for line in fh:
            if line.startswith("LINK"):
                return True
    return False


# ---------------------------------------------------------------------------
# Utilities for atom / LINK consistency
# ---------------------------------------------------------------------------

def build_atom_index(
    atom_rec_list: Iterable[edit_pdb_atom.pdb_atom_record],
) -> Dict[Tuple[str, int, str], List[edit_pdb_atom.pdb_atom_record]]:
    """
    Build a simple index:
        (chainID, resSeq, atomName) -> list of atom records

    This is useful when you want to look up atoms referenced by LINK
    records, e.g. to recompute P--O distances after a geometry change.
    """
    index: Dict[Tuple[str, int, str], List[edit_pdb_atom.pdb_atom_record]] = {}
    for atom in atom_rec_list:
        if getattr(atom, "recordName", "") not in ("ATOM", "HETATM"):
            continue
        key = (atom.chainID, atom.resSeq, atom.name.strip())
        index.setdefault(key, []).append(atom)
    return index


def find_atom(
    atom_index: Dict[Tuple[str, int, str], List[edit_pdb_atom.pdb_atom_record]],
    chain_id: str,
    res_seq: int,
    atom_name: str,
) -> Optional[edit_pdb_atom.pdb_atom_record]:
    """
    Find a unique atom in the index given (chainID, resSeq, atom_name).

    Returns the first atom if found, or None if no entry exists for the key.
    If multiple atoms share the same (chain,resSeq,name), the first is
    returned; callers that require uniqueness should check the list length
    directly in the index.
    """
    key = (chain_id, res_seq, atom_name.strip())
    atoms = atom_index.get(key)
    if not atoms:
        return None
    return atoms[0]


def update_link_distances_from_atoms(
    link_rec_list: List[pdb_link_record],
    atom_rec_list: List[edit_pdb_atom.pdb_atom_record],
) -> None:
    """
    For each LINK in link_rec_list, find the corresponding atom pair in
    atom_rec_list (by chainID, resSeq, atom name) and recompute the
    distance field from their coordinates.

    LINKs that refer to atoms not present in atom_rec_list are left
    unchanged (no error is raised).
    """
    index = build_atom_index(atom_rec_list)

    for link in link_rec_list:
        key1 = (link.chainID1, link.resSeq1, link.name1.strip())
        key2 = (link.chainID2, link.resSeq2, link.name2.strip())

        atoms1 = index.get(key1)
        atoms2 = index.get(key2)

        if not atoms1 or not atoms2:
            # silently skip if we can't find both atoms
            continue

        a1 = atoms1[0]
        a2 = atoms2[0]

        dx = a1.x - a2.x
        dy = a1.y - a2.y
        dz = a1.z - a2.z
        d = math.sqrt(dx * dx + dy * dy + dz * dz)

        link.update_distance(d)


if __name__ == "__main__":
    # Simple self-test / demo: read a PDB from stdin and echo back LINK+ATOM.
    import sys

    atom_list: List[edit_pdb_atom.pdb_atom_record] = []
    link_list: List[pdb_link_record] = []

    file2rec_link(sys.stdin, atom_list, link_list)

    sys.stderr.write(
        f"[edit_pdb_link] Parsed {len(link_list)} LINK records and "
        f"{len(atom_list)} ATOM/HETATM/TER records from stdin.\n"
    )

    update_link_distances_from_atoms(link_list, atom_list)

    rec2file_link(atom_list, link_list, sys.stdout, reorder_serial=False)
