#!/usr/bin/env python3
"""add_pdb_link_record.py

Build LINK-driven chain topology from selected endpoints in a PDB file.

Supported actions
-----------------
1) Automatic circularization
   - For each selected chain, the script finds:
       nucleic acid: P on the 5' residue and O3' on the 3' residue, or
       peptide: N on the N-terminal residue and C on the C-terminal residue
   - The selected chain is added as a LINK-driven topology edit.

2) Manual LINK creation in the GUI
   - The user selects nucleic-acid P/O3' or peptide N/C mode and specifies the
     residue number and chain ID for each endpoint.
   - Clicking "Add Link" stages a default-checked LINK entry below.
   - Clicking "Run" rebuilds the topology from the checked automatic and
     manual LINKs in one pass.

Behavior
--------
- Selected endpoints are used to split chain segments, merge them into new
  output chains, reassign chain IDs, rewrite TER records, and renumber residues.
- Existing LINK lines are preserved, remapped to rebuilt residue numbering when
  possible, and written before any newly selected LINK edits.
- Nucleic-acid output defaults to "_circ" before the extension; peptide output
  defaults to "_peptide_circ".

Important notes
---------------
- 5'/3' and N/C termini are inferred from residue numbering (min/max resSeq)
  within each chain for automatic circularization.
- O3' atom name may appear as "O3'" or "O3*"; both are supported.
- The output topology assumes the selected P/O3' endpoints are valid chain
  endpoints for the intended segment merge.

GUI
---
- If run with no arguments, or with --gui, a small Tk GUI is launched.
- After selecting a PDB and molecule type, the GUI reports either P-O3' or N-C
  terminal distance for every chain.
- The GUI also includes a molecule-type-aware manual LINK panel and a pending
  LINK list.

License/usage
-------------
This script is intended as a lightweight PDB helper. It does not validate the
chemical correctness of the resulting topology.
"""

from __future__ import annotations

import argparse
import math
import shlex
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

try:
    from re_helix_lib.gui_icon import apply_optional_icon
except ImportError:  # pragma: no cover - direct script execution fallback
    from gui_icon import apply_optional_icon

TOOL_NAME = "Add PDB LINK Record"
VERSION = "1.2"


# ----------------------------
# Small PDB helpers
# ----------------------------


@dataclass(frozen=True)
class AtomInfo:
    """Minimal info extracted from an ATOM/HETATM line."""

    name: str
    resName: str
    chainID: str
    resSeq: int
    x: float
    y: float
    z: float
    recordName: str = "ATOM"
    serial: int = 0
    raw_line: str = ""


@dataclass
class ChainEndInfo:
    chainID: str
    first_resSeq: int
    last_resSeq: int
    first_resName: str
    last_resName: str
    p_atom: Optional[AtomInfo]
    o3_atom: Optional[AtomInfo]
    distance: Optional[float]
    n_atom: Optional[AtomInfo]
    c_atom: Optional[AtomInfo]
    peptide_distance: Optional[float]


@dataclass
class PendingManualLink:
    label: str
    endpoint1_atom: AtomInfo
    endpoint2_atom: AtomInfo
    endpoint1: Tuple[str, int, str]
    endpoint2: Tuple[str, int, str]
    var: object


O3_CANDIDATES: Tuple[str, ...] = ("O3'", "O3*", "O3")


def _safe_int(s: str) -> Optional[int]:
    s = s.strip()
    if not s:
        return None
    try:
        return int(s)
    except ValueError:
        return None


def _safe_float(s: str) -> Optional[float]:
    s = s.strip()
    if not s:
        return None
    try:
        return float(s)
    except ValueError:
        return None


def parse_atom_line(line: str) -> Optional[AtomInfo]:
    """Parse an ATOM/HETATM line; return AtomInfo or None if not parseable."""
    if not (line.startswith("ATOM") or line.startswith("HETATM")):
        return None

    try:
        name = line[12:16].strip()
        resName = line[17:20].strip()
        chainID = line[21:22]  # single character; may be a space
        resSeq = _safe_int(line[22:26])
        x = _safe_float(line[30:38])
        y = _safe_float(line[38:46])
        z = _safe_float(line[46:54])
        serial = _safe_int(line[6:11]) or 0
    except Exception:
        return None

    if resSeq is None or x is None or y is None or z is None:
        return None

    return AtomInfo(
        name=name,
        resName=resName,
        chainID=chainID,
        resSeq=resSeq,
        x=x,
        y=y,
        z=z,
        recordName=line[:6].strip() or "ATOM",
        serial=serial,
        raw_line=line.rstrip("\n"),
    )


def distance(a: AtomInfo, b: AtomInfo) -> float:
    dx = a.x - b.x
    dy = a.y - b.y
    dz = a.z - b.z
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def format_link_line(
    atom1_name: str,
    resName1: str,
    chainID1: str,
    resSeq1: int,
    atom2_name: str,
    resName2: str,
    chainID2: str,
    resSeq2: int,
    dist: Optional[float] = None,
    sym1: str = "1555",
    sym2: str = "1555",
) -> str:
    """Return a canonical LINK record line (with newline)."""
    dist_str = "     "
    if dist is not None:
        dist_str = f"{dist:5.2f}"

    line = (
        "LINK        "
        f"{atom1_name:>4s} {resName1:>3s} {chainID1:1s}{resSeq1:4d}"
        "                "
        f"{atom2_name:>4s} {resName2:>3s} {chainID2:1s}{resSeq2:4d}"
        f"     {sym1:>4s}   {sym2:>4s} {dist_str}"
    )
    return line + "\n"


def _norm_atom_for_link(name: str) -> str:
    """Normalize atom names for duplicate-link detection."""
    n = name.strip().upper()
    if n in {"O3'", "O3*", "O3"}:
        return "O3"
    if n == "P":
        return "P"
    return n


@dataclass(frozen=True)
class ExistingLinkRecord:
    """Parsed LINK record from the input PDB."""

    atom1_name: str
    resName1: str
    chainID1: str
    resSeq1: int
    atom2_name: str
    resName2: str
    chainID2: str
    resSeq2: int
    sym1: str = "1555"
    sym2: str = "1555"
    dist: Optional[float] = None
    raw_line: str = ""


def _link_endpoint_pair_key(
    chainID1: str,
    resSeq1: int,
    atom1_name: str,
    chainID2: str,
    resSeq2: int,
    atom2_name: str,
) -> frozenset[Tuple[str, int, str]]:
    return frozenset(
        (
            _format_endpoint(chainID1, resSeq1, atom1_name),
            _format_endpoint(chainID2, resSeq2, atom2_name),
        )
    )


def parse_existing_link_record_line(line: str) -> Optional[ExistingLinkRecord]:
    """Parse a LINK line while preserving enough fields to rewrite it later."""
    if not line.startswith("LINK"):
        return None

    padded = line.rstrip("\n").ljust(80)

    try:
        atom1 = padded[12:16].strip()
        resName1 = padded[17:20].strip()
        chain1 = padded[21:22]
        resSeq1 = int(padded[22:26])
        atom2 = padded[42:46].strip()
        resName2 = padded[47:50].strip()
        chain2 = padded[51:52]
        resSeq2 = int(padded[52:56])
        sym1 = padded[59:65].strip() or "1555"
        sym2 = padded[66:72].strip() or "1555"
        dist = _safe_float(padded[73:78])
        if atom1 and atom2:
            return ExistingLinkRecord(
                atom1_name=atom1,
                resName1=resName1,
                chainID1=chain1,
                resSeq1=resSeq1,
                atom2_name=atom2,
                resName2=resName2,
                chainID2=chain2,
                resSeq2=resSeq2,
                sym1=sym1,
                sym2=sym2,
                dist=dist,
                raw_line=line if line.endswith("\n") else line + "\n",
            )
    except Exception:
        pass

    toks = line.split()
    if len(toks) < 9 or toks[0] != "LINK":
        return None
    try:
        dist = _safe_float(toks[11]) if len(toks) >= 12 else None
        return ExistingLinkRecord(
            atom1_name=toks[1],
            resName1=toks[2],
            chainID1=toks[3],
            resSeq1=int(toks[4]),
            atom2_name=toks[5],
            resName2=toks[6],
            chainID2=toks[7],
            resSeq2=int(toks[8]),
            sym1=toks[9] if len(toks) >= 10 else "1555",
            sym2=toks[10] if len(toks) >= 11 else "1555",
            dist=dist,
            raw_line=line if line.endswith("\n") else line + "\n",
        )
    except Exception:
        return None


def collect_existing_link_records(
    lines: Sequence[str],
) -> Tuple[List[ExistingLinkRecord], List[str]]:
    """Return parsed and unparsable input LINK records."""
    parsed: List[ExistingLinkRecord] = []
    unparsed: List[str] = []
    for line in lines:
        if not line.startswith("LINK"):
            continue
        rec = parse_existing_link_record_line(line)
        if rec is None:
            unparsed.append(line if line.endswith("\n") else line + "\n")
        else:
            parsed.append(rec)
    return parsed, unparsed


def collect_existing_link_entries(lines: Sequence[str]) -> List[object]:
    """Return input LINK records in file order as parsed records or raw lines."""
    entries: List[object] = []
    for line in lines:
        if not line.startswith("LINK"):
            continue
        rec = parse_existing_link_record_line(line)
        if rec is None:
            entries.append(line if line.endswith("\n") else line + "\n")
        else:
            entries.append(rec)
    return entries


def parse_existing_links(lines: Sequence[str]) -> Set[frozenset[Tuple[str, int, str]]]:
    """Parse existing LINK records into a set of normalized endpoint pairs.

    Each stored item is frozenset({(chain,resSeq,atom), (chain,resSeq,atom)}).
    Atom names are normalized so O3'/O3* are treated as the same.
    """
    links: Set[frozenset[Tuple[str, int, str]]] = set()
    parsed, _unparsed = collect_existing_link_records(lines)
    for rec in parsed:
        links.add(
            _link_endpoint_pair_key(
                rec.chainID1,
                rec.resSeq1,
                rec.atom1_name,
                rec.chainID2,
                rec.resSeq2,
                rec.atom2_name,
            )
        )
    return links


def parse_existing_link_endpoints(lines: Sequence[str]) -> Set[Tuple[str, int, str]]:
    """Return the set of all atom endpoints that already appear in LINK records."""
    used: Set[Tuple[str, int, str]] = set()
    for pair in parse_existing_links(lines):
        for endpoint in pair:
            used.add(endpoint)
    return used


def build_atom_index(atoms: Iterable[AtomInfo]) -> Dict[str, Dict[int, Dict[str, AtomInfo]]]:
    """Index atoms by chain -> residue number -> atom name."""
    out: Dict[str, Dict[int, Dict[str, AtomInfo]]] = {}
    for a in atoms:
        out.setdefault(a.chainID, {}).setdefault(a.resSeq, {})
        out[a.chainID][a.resSeq].setdefault(a.name.strip().upper(), a)
    return out


def build_chain_end_info(atoms: Iterable[AtomInfo]) -> Dict[str, ChainEndInfo]:
    """Identify nucleic-acid and peptide terminal atoms for each chain."""
    chain_res_atoms = build_atom_index(atoms)
    out: Dict[str, ChainEndInfo] = {}

    for chain, res_map in chain_res_atoms.items():
        if not res_map:
            continue

        resSeqs = sorted(res_map.keys())
        first_r = resSeqs[0]
        last_r = resSeqs[-1]

        first_atoms = res_map[first_r]
        last_atoms = res_map[last_r]

        p_atom = first_atoms.get("P")
        n_atom = first_atoms.get("N")
        c_atom = last_atoms.get("C")

        o3_atom: Optional[AtomInfo] = None
        for nm in O3_CANDIDATES:
            key = nm.upper()
            if key in last_atoms:
                o3_atom = last_atoms[key]
                break

        dist: Optional[float] = None
        if p_atom is not None and o3_atom is not None:
            dist = distance(p_atom, o3_atom)

        peptide_dist: Optional[float] = None
        if n_atom is not None and c_atom is not None:
            peptide_dist = distance(n_atom, c_atom)

        first_resName = p_atom.resName if p_atom is not None else next(iter(first_atoms.values())).resName
        last_resName = o3_atom.resName if o3_atom is not None else next(iter(last_atoms.values())).resName

        out[chain] = ChainEndInfo(
            chainID=chain,
            first_resSeq=first_r,
            last_resSeq=last_r,
            first_resName=first_resName,
            last_resName=last_resName,
            p_atom=p_atom,
            o3_atom=o3_atom,
            distance=dist,
            n_atom=n_atom,
            c_atom=c_atom,
            peptide_distance=peptide_dist,
        )

    return out


def default_output_path(inp_path: Path) -> Path:
    """Insert '_circ' before the final suffix (extension)."""
    return inp_path.with_name(f"{inp_path.stem}_circ{inp_path.suffix}")


def default_peptide_output_path(inp_path: Path) -> Path:
    """Insert '_peptide_circ' before the final suffix."""
    return inp_path.with_name(f"{inp_path.stem}_peptide_circ{inp_path.suffix}")


def _find_first_atom_line_index(lines: Sequence[str]) -> int:
    for i, ln in enumerate(lines):
        rec = ln[:6].strip()
        if rec in {"ATOM", "HETATM"}:
            return i
    return 0


def _selected_link_to_line(item: object) -> str:
    if isinstance(item, SelectedLinkSpec):
        return item.to_link_line()
    if isinstance(item, str):
        return item if item.endswith("\n") else item + "\n"
    raise TypeError(f"Expected a LINK string or SelectedLinkSpec, got {type(item).__name__}.")


def _insert_link_lines(lines: List[str], new_links: Sequence[object]) -> List[str]:
    insert_idx = _find_first_atom_line_index(lines)
    out_lines = list(lines[:insert_idx])
    out_lines.extend(_selected_link_to_line(link) for link in new_links)
    out_lines.extend(lines[insert_idx:])
    return out_lines


def _format_endpoint(chain: str, resSeq: int, atom_name: str) -> Tuple[str, int, str]:
    return (chain.strip() or chain, resSeq, _norm_atom_for_link(atom_name))


def _validate_endpoint_not_already_linked(
    endpoint: Tuple[str, int, str],
    used_endpoints: Set[Tuple[str, int, str]],
    label: str,
) -> None:
    if endpoint in used_endpoints:
        chain, resSeq, atom = endpoint
        raise ValueError(
            f"{label} is already used in an existing LINK record: "
            f"atom {atom} on residue {resSeq} in chain {chain!r}."
        )


def _add_automatic_terminal_links(
    inp_path: Path,
    chains: Optional[Sequence[str]] = None,
    out_path: Optional[Path] = None,
    verbose: bool = True,
    link_type: str = "nucleic",
) -> Tuple[Path, Dict[str, ChainEndInfo], List[str]]:
    """Circularize selected nucleic-acid or peptide chains with LINK records.

    Returns:
        out_path, chain_info, new_link_lines
    """
    inp_path = Path(inp_path)
    if link_type not in {"nucleic", "peptide"}:
        raise ValueError(f"Unknown automatic LINK type: {link_type!r}")
    if out_path is None:
        out_path = (
            default_peptide_output_path(inp_path)
            if link_type == "peptide"
            else default_output_path(inp_path)
        )

    lines = inp_path.read_text(errors="replace").splitlines(True)

    atoms: List[AtomInfo] = []
    for ln in lines:
        a = parse_atom_line(ln)
        if a is not None:
            atoms.append(a)

    chain_info = build_chain_end_info(atoms)
    if not chain_info:
        raise ValueError("No ATOM/HETATM records found (or could not parse any).")

    requested: Optional[Set[str]]
    if chains is None or len(chains) == 0:
        requested = None
    else:
        requested = set(chains)

    if requested is not None:
        missing = [c for c in sorted(requested) if c not in chain_info]
        if missing:
            raise ValueError("Requested chain(s) not found in PDB: " + ", ".join(repr(c) for c in missing))

    existing_links = parse_existing_links(lines)
    existing_link_endpoints = parse_existing_link_endpoints(lines)

    new_links: List[str] = []

    for chain in sorted(chain_info.keys()):
        if requested is not None and chain not in requested:
            continue

        info = chain_info[chain]
        if link_type == "peptide":
            endpoint1_atom = info.n_atom
            endpoint2_atom = info.c_atom
            endpoint1_name = "N"
            endpoint2_name = "C"
            endpoint1_label = "N-terminal N"
            endpoint2_label = "C-terminal C"
            link_distance = info.peptide_distance
        else:
            endpoint1_atom = info.p_atom
            endpoint2_atom = info.o3_atom
            endpoint1_name = "P"
            endpoint2_name = "O3'"
            endpoint1_label = "5' P"
            endpoint2_label = "3' O3'"
            link_distance = info.distance

        if endpoint1_atom is None:
            if verbose:
                print(
                    f"[WARN] Chain {chain!r}: missing {endpoint1_label} atom on "
                    f"residue {info.first_resSeq}; skipping."
                )
            continue
        if endpoint2_atom is None:
            if verbose:
                print(
                    f"[WARN] Chain {chain!r}: missing {endpoint2_label} atom on "
                    f"residue {info.last_resSeq}; skipping."
                )
            continue

        e1 = _format_endpoint(chain, info.first_resSeq, endpoint1_name)
        e2 = _format_endpoint(chain, info.last_resSeq, endpoint2_name)

        key = frozenset((e1, e2))
        if key in existing_links:
            if verbose:
                print(
                    f"[INFO] Chain {chain!r}: terminal {endpoint1_name}-{endpoint2_name} "
                    "LINK already present; not adding duplicate."
                )
            continue

        # Prevent duplicate use of either endpoint in another LINK record
        if e1 in existing_link_endpoints:
            raise ValueError(
                f"Chain {chain!r}: the {endpoint1_label} atom on residue "
                f"{info.first_resSeq} is already used in an existing LINK record."
            )
        if e2 in existing_link_endpoints:
            raise ValueError(
                f"Chain {chain!r}: the {endpoint2_label} atom on residue "
                f"{info.last_resSeq} is already used in an existing LINK record."
            )

        link_line = format_link_line(
            atom1_name=endpoint1_atom.name,
            resName1=endpoint1_atom.resName,
            chainID1=endpoint1_atom.chainID,
            resSeq1=endpoint1_atom.resSeq,
            atom2_name=endpoint2_atom.name,
            resName2=endpoint2_atom.resName,
            chainID2=endpoint2_atom.chainID,
            resSeq2=endpoint2_atom.resSeq,
            dist=link_distance,
        )
        new_links.append(link_line)

    out_lines = _insert_link_lines(lines, new_links)
    out_path.write_text("".join(out_lines))

    return out_path, chain_info, new_links


def circularize_pdb(
    inp_path: Path,
    chains: Optional[Sequence[str]] = None,
    out_path: Optional[Path] = None,
    verbose: bool = True,
) -> Tuple[Path, Dict[str, ChainEndInfo], List[str]]:
    """Circularize selected nucleic-acid chains with terminal P-O3' LINKs."""
    return _add_automatic_terminal_links(
        inp_path, chains=chains, out_path=out_path, verbose=verbose, link_type="nucleic"
    )


def cyclize_peptide_pdb(
    inp_path: Path,
    chains: Optional[Sequence[str]] = None,
    out_path: Optional[Path] = None,
    verbose: bool = True,
) -> Tuple[Path, Dict[str, ChainEndInfo], List[str]]:
    """Cyclize selected peptide chains with terminal N-C LINKs."""
    return _add_automatic_terminal_links(
        inp_path, chains=chains, out_path=out_path, verbose=verbose, link_type="peptide"
    )


def _prepare_named_manual_link_spec(
    lines: Sequence[str],
    atoms: Sequence[AtomInfo],
    endpoint1_chain: str,
    endpoint1_resseq: int,
    endpoint1_name: str,
    endpoint2_chain: str,
    endpoint2_resseq: int,
    endpoint2_name: str,
) -> Tuple[AtomInfo, AtomInfo, Tuple[str, int, str], Tuple[str, int, str]]:
    """Validate and build a named manual LINK record without writing output.

    Manual endpoints may be on internal residues; the output stage will split
    segments at the selected endpoints and rebuild chain IDs / TER records.
    """
    atom_index = build_atom_index(atoms)

    endpoint1_chain = endpoint1_chain.strip()
    endpoint2_chain = endpoint2_chain.strip()

    if endpoint1_chain == "":
        raise ValueError(f"Please enter a valid chain ID for the {endpoint1_name} atom.")
    if endpoint2_chain == "":
        raise ValueError(f"Please enter a valid chain ID for the {endpoint2_name} atom.")

    try:
        endpoint1_res_map = atom_index[endpoint1_chain]
    except KeyError:
        raise ValueError(f"Chain {endpoint1_chain!r} was not found in the PDB file.")

    try:
        endpoint2_res_map = atom_index[endpoint2_chain]
    except KeyError:
        raise ValueError(f"Chain {endpoint2_chain!r} was not found in the PDB file.")

    if endpoint1_resseq not in endpoint1_res_map:
        raise ValueError(
            f"Residue {endpoint1_resseq} was not found in chain {endpoint1_chain!r}."
        )
    if endpoint2_resseq not in endpoint2_res_map:
        raise ValueError(
            f"Residue {endpoint2_resseq} was not found in chain {endpoint2_chain!r}."
        )

    def find_atom(residue_atoms: Dict[str, AtomInfo], atom_name: str) -> Optional[AtomInfo]:
        if _norm_atom_for_link(atom_name) == "O3":
            for candidate in O3_CANDIDATES:
                found = residue_atoms.get(candidate.upper())
                if found is not None:
                    return found
            return None
        return residue_atoms.get(atom_name.strip().upper())

    endpoint1_atom = find_atom(endpoint1_res_map[endpoint1_resseq], endpoint1_name)
    if endpoint1_atom is None:
        raise ValueError(
            f"No {endpoint1_name} atom was found on residue {endpoint1_resseq} "
            f"in chain {endpoint1_chain!r}."
        )
    endpoint2_atom = find_atom(endpoint2_res_map[endpoint2_resseq], endpoint2_name)
    if endpoint2_atom is None:
        atom_label = "O3' (or O3*)" if _norm_atom_for_link(endpoint2_name) == "O3" else endpoint2_name
        raise ValueError(
            f"No {atom_label} atom was found on residue {endpoint2_resseq} "
            f"in chain {endpoint2_chain!r}."
        )

    e1 = _format_endpoint(endpoint1_chain, endpoint1_resseq, endpoint1_atom.name)
    e2 = _format_endpoint(endpoint2_chain, endpoint2_resseq, endpoint2_atom.name)
    return endpoint1_atom, endpoint2_atom, e1, e2


def _prepare_manual_link_spec(
    lines: Sequence[str],
    atoms: Sequence[AtomInfo],
    p_chain: str,
    p_resseq: int,
    o3_chain: str,
    o3_resseq: int,
) -> Tuple[AtomInfo, AtomInfo, Tuple[str, int, str], Tuple[str, int, str]]:
    """Backward-compatible P/O3' manual LINK preparation."""
    return _prepare_named_manual_link_spec(
        lines,
        atoms,
        endpoint1_chain=p_chain,
        endpoint1_resseq=p_resseq,
        endpoint1_name="P",
        endpoint2_chain=o3_chain,
        endpoint2_resseq=o3_resseq,
        endpoint2_name="O3'",
    )


def manual_link_record(
    inp_path: Path,
    p_chain: str,
    p_resseq: int,
    o3_chain: str,
    o3_resseq: int,
    out_path: Optional[Path] = None,
    verbose: bool = True,
) -> Tuple[Path, List[str]]:
    """Create one manual LINK record and rebuild chain segments."""
    inp_path = Path(inp_path)
    if out_path is None:
        out_path = default_output_path(inp_path)

    lines = inp_path.read_text(errors="replace").splitlines(True)

    atoms: List[AtomInfo] = []
    for ln in lines:
        a = parse_atom_line(ln)
        if a is not None:
            atoms.append(a)

    p_atom, o3_atom, e1, e2 = _prepare_manual_link_spec(
        lines=lines,
        atoms=atoms,
        p_chain=p_chain,
        p_resseq=p_resseq,
        o3_chain=o3_chain,
        o3_resseq=o3_resseq,
    )

    selected = [
        SelectedLinkSpec(
            label="manual link",
            endpoint1_atom=p_atom,
            endpoint2_atom=o3_atom,
            endpoint1=e1,
            endpoint2=e2,
        )
    ]
    out_text = build_relinked_pdb_text(lines, atoms, selected)
    out_path.write_text(out_text)

    if verbose:
        print(f"Wrote output: {out_path}")
        print("Added 1 LINK record and rebuilt chain topology.")

    return out_path, ["manual LINK"]


def _parse_chain_list(chain_args: Optional[Sequence[str]]) -> List[str]:
    """Parse chain list from CLI (accepts comma-separated and/or space-separated)."""
    if not chain_args:
        return []
    out: List[str] = []
    for token in chain_args:
        if token is None:
            continue
        t = token.strip()
        if not t:
            continue
        parts = [p for p in t.replace(";", ",").split(",") if p != ""]
        for p in parts:
            p = p.strip()
            if not p:
                continue
            if len(p) > 1 and all(ch.strip() for ch in p) and (" " not in p):
                if p.isalnum() and not p[1:].isdigit():
                    out.extend(list(p))
                else:
                    out.append(p)
            else:
                out.append(p)
    seen: Set[str] = set()
    uniq: List[str] = []
    for c in out:
        if c not in seen:
            seen.add(c)
            uniq.append(c)
    return uniq



# ----------------------------
# Topology reconstruction for selected LINKs
# ----------------------------


@dataclass(frozen=True)
class SelectedLinkSpec:
    label: str
    endpoint1_atom: AtomInfo
    endpoint2_atom: AtomInfo
    endpoint1: Tuple[str, int, str]
    endpoint2: Tuple[str, int, str]

    def to_link_line(self) -> str:
        return format_link_line(
            atom1_name=self.endpoint1_atom.name,
            resName1=self.endpoint1_atom.resName,
            chainID1=self.endpoint1[0],
            resSeq1=self.endpoint1[1],
            atom2_name=self.endpoint2_atom.name,
            resName2=self.endpoint2_atom.resName,
            chainID2=self.endpoint2[0],
            resSeq2=self.endpoint2[1],
            dist=distance(self.endpoint1_atom, self.endpoint2_atom),
        )

    def __str__(self) -> str:
        return self.to_link_line().rstrip("\n")


def _find_all_atom_indices(lines: Sequence[str]) -> List[int]:
    return [i for i, ln in enumerate(lines) if ln.startswith("ATOM") or ln.startswith("HETATM")]


def _extract_preserved_non_topology_lines(lines: Sequence[str]) -> Tuple[List[str], List[str]]:
    """Return (prefix, suffix) non-ATOM/HETATM lines, dropping old TER/LINK lines."""
    atom_indices = _find_all_atom_indices(lines)
    if not atom_indices:
        preserved = [ln for ln in lines if not ln.startswith(("TER", "LINK"))]
        return preserved, []

    first_atom = atom_indices[0]
    last_atom = atom_indices[-1]
    prefix = [ln for ln in lines[:first_atom] if not ln.startswith(("TER", "LINK"))]
    suffix = [ln for ln in lines[last_atom + 1 :] if not ln.startswith(("TER", "LINK"))]
    return prefix, suffix


def _build_chain_residue_orders(atoms: Sequence[AtomInfo]) -> Dict[str, List[int]]:
    """Build residue order for each chain based on first appearance in the file."""
    out: Dict[str, List[int]] = {}
    seen: Dict[str, Set[int]] = {}
    for a in atoms:
        if a.recordName not in ("ATOM", "HETATM"):
            continue
        out.setdefault(a.chainID, [])
        seen.setdefault(a.chainID, set())
        if a.resSeq not in seen[a.chainID]:
            seen[a.chainID].add(a.resSeq)
            out[a.chainID].append(a.resSeq)
    return out


def _build_residue_atom_map(atoms: Sequence[AtomInfo]) -> Dict[Tuple[str, int], List[AtomInfo]]:
    out: Dict[Tuple[str, int], List[AtomInfo]] = {}
    for a in atoms:
        if a.recordName not in ("ATOM", "HETATM"):
            continue
        out.setdefault((a.chainID, a.resSeq), []).append(a)
    return out


def _rewrite_atom_line(atom: AtomInfo, new_chain: str, new_resseq: int, new_serial: int) -> str:
    line = (atom.raw_line or "").rstrip("\n")
    if not line:
        # Fallback if raw_line is unavailable
        line = (
            f"{atom.recordName:<6}{new_serial:>5d} "
            f"{atom.name:>4s} {atom.resName:>3s} {new_chain:1s}{new_resseq:>4d}    "
            f"{atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f}"
        )
        return line + "\n"

    line = line.ljust(80)
    chars = list(line)
    rec = atom.recordName if atom.recordName in ("ATOM", "HETATM") else "ATOM"
    rec_field = f"{rec:<6s}"
    chars[0:6] = list(rec_field)
    chars[6:11] = list(f"{new_serial:>5d}")
    chars[21] = new_chain if new_chain else " "
    chars[22:26] = list(f"{new_resseq:>4d}")
    return "".join(chars).rstrip() + "\n"


def _format_ter_line(serial: int, resName: str, chainID: str, resSeq: int) -> str:
    return f"TER   {serial:>5d}      {resName:>3s} {chainID:1s}{resSeq:>4d}\n"


def _remap_existing_link_line(
    rec: ExistingLinkRecord,
    residue_to_new: Dict[Tuple[str, int], Tuple[str, int]],
) -> Tuple[str, frozenset[Tuple[str, int, str]]]:
    key1 = (rec.chainID1, rec.resSeq1)
    key2 = (rec.chainID2, rec.resSeq2)
    if key1 in residue_to_new and key2 in residue_to_new:
        new_chain1, new_resSeq1 = residue_to_new[key1]
        new_chain2, new_resSeq2 = residue_to_new[key2]
        return (
            format_link_line(
                atom1_name=rec.atom1_name,
                resName1=rec.resName1,
                chainID1=new_chain1,
                resSeq1=new_resSeq1,
                atom2_name=rec.atom2_name,
                resName2=rec.resName2,
                chainID2=new_chain2,
                resSeq2=new_resSeq2,
                dist=rec.dist,
                sym1=rec.sym1,
                sym2=rec.sym2,
            ),
            _link_endpoint_pair_key(
                new_chain1,
                new_resSeq1,
                rec.atom1_name,
                new_chain2,
                new_resSeq2,
                rec.atom2_name,
            ),
        )

    line = rec.raw_line if rec.raw_line.endswith("\n") else rec.raw_line + "\n"
    return (
        line,
        _link_endpoint_pair_key(
            rec.chainID1,
            rec.resSeq1,
            rec.atom1_name,
            rec.chainID2,
            rec.resSeq2,
            rec.atom2_name,
        ),
    )


def _component_sort_key(component: List[int], segments: Dict[int, dict], chain_order: Dict[str, int]) -> Tuple[int, int, int]:
    best = None
    for seg_id in component:
        seg = segments[seg_id]
        key = (chain_order.get(seg["chainID"], 9999), seg["residues"][0], seg["residues"][-1])
        if best is None or key < best:
            best = key
    return best or (9999, 9999, 9999)


def _build_segments_from_selected_links(
    atoms: Sequence[AtomInfo],
    selected_links: Sequence[SelectedLinkSpec],
) -> Tuple[Dict[int, dict], Dict[Tuple[str, int], int], Dict[int, int], Dict[int, int]]:
    """Split chains at selected P/O3 endpoints and build segment graph helpers.

    Returns:
        segments: seg_id -> dict(chainID, residues, first_res, last_res)
        residue_to_segment: (chain,resSeq) -> seg_id
        seg_out: seg_id -> next seg_id (or -1)
        seg_in: seg_id -> prev seg_id (or -1)
    """
    chain_res_order = _build_chain_residue_orders(atoms)
    residue_atom_map = _build_residue_atom_map(atoms)

    if not chain_res_order:
        raise ValueError("No ATOM/HETATM records found (or could not parse any).")

    cut_before: Dict[str, Set[int]] = {}
    cut_after: Dict[str, Set[int]] = {}
    for link in selected_links:
        first_chain, first_res, _ = link.endpoint1
        second_chain, second_res, _ = link.endpoint2
        cut_before.setdefault(first_chain, set()).add(first_res)
        cut_after.setdefault(second_chain, set()).add(second_res)

    segments: Dict[int, dict] = {}
    residue_to_segment: Dict[Tuple[str, int], int] = {}
    seg_id = 0

    for chain, residues in chain_res_order.items():
        if not residues:
            continue
        start_idx = 0
        for i, res in enumerate(residues):
            is_last = i == len(residues) - 1
            next_res = residues[i + 1] if i + 1 < len(residues) else None
            boundary = False
            if res in cut_after.get(chain, set()):
                boundary = True
            if next_res is not None and next_res in cut_before.get(chain, set()):
                boundary = True
            if is_last:
                boundary = True

            if boundary:
                seg_res = residues[start_idx : i + 1]
                if not seg_res:
                    continue
                segments[seg_id] = {
                    "chainID": chain,
                    "residues": seg_res,
                    "atoms": {r: residue_atom_map[(chain, r)] for r in seg_res if (chain, r) in residue_atom_map},
                }
                for r in seg_res:
                    residue_to_segment[(chain, r)] = seg_id
                seg_id += 1
                start_idx = i + 1

    seg_out: Dict[int, int] = {sid: -1 for sid in segments}
    seg_in: Dict[int, int] = {sid: -1 for sid in segments}

    for link in selected_links:
        first_chain, first_res, _ = link.endpoint1
        second_chain, second_res, _ = link.endpoint2
        s_src = residue_to_segment.get((second_chain, second_res))
        s_tgt = residue_to_segment.get((first_chain, first_res))
        if s_src is None or s_tgt is None:
            raise ValueError(
                "Could not place selected LINK endpoints into segments: "
                f"{second_res}{second_chain} -> {first_res}{first_chain}."
            )

        if segments[s_src]["residues"][-1] != second_res:
            raise ValueError(
                f"The selected {link.endpoint2_atom.name} endpoint "
                f"{second_res}{second_chain} is not the last residue of its segment."
            )
        if segments[s_tgt]["residues"][0] != first_res:
            raise ValueError(
                f"The selected {link.endpoint1_atom.name} endpoint "
                f"{first_res}{first_chain} is not the first residue of its segment."
            )

        if seg_out[s_src] not in (-1, s_tgt):
            raise ValueError(
                f"Segment ending at {second_res}{second_chain} already connects to another selected LINK."
            )
        if seg_in[s_tgt] not in (-1, s_src):
            raise ValueError(
                f"Segment starting at {first_res}{first_chain} already receives another selected LINK."
            )

        seg_out[s_src] = s_tgt
        seg_in[s_tgt] = s_src

    return segments, residue_to_segment, seg_out, seg_in


def build_relinked_pdb_text(
    lines: Sequence[str],
    atoms: Sequence[AtomInfo],
    selected_links: Sequence[SelectedLinkSpec],
) -> str:
    """Rebuild the PDB topology from selected LINKs, rewriting chain IDs and TER."""
    if not selected_links:
        raise ValueError("Please select at least one LINK to add.")

    existing_link_entries = collect_existing_link_entries(lines)
    prefix_lines, suffix_lines = _extract_preserved_non_topology_lines(lines)
    segments, residue_to_segment, seg_out, seg_in = _build_segments_from_selected_links(atoms, selected_links)

    # Build connected components of segments
    undirected: Dict[int, Set[int]] = {sid: set() for sid in segments}
    for src, dst in seg_out.items():
        if dst != -1:
            undirected[src].add(dst)
            undirected[dst].add(src)

    visited: Set[int] = set()
    components: List[List[int]] = []

    for sid in sorted(segments, key=lambda s: (segments[s]["chainID"], segments[s]["residues"][0], segments[s]["residues"][-1])):
        if sid in visited:
            continue
        stack = [sid]
        comp: List[int] = []
        visited.add(sid)
        while stack:
            cur = stack.pop()
            comp.append(cur)
            for nb in undirected[cur]:
                if nb not in visited:
                    visited.add(nb)
                    stack.append(nb)
        components.append(comp)

    # Order components by the first appearance of their original chain IDs in the file
    chain_order: Dict[str, int] = {}
    for a in atoms:
        if a.recordName not in ("ATOM", "HETATM"):
            continue
        if a.chainID not in chain_order:
            chain_order[a.chainID] = len(chain_order)
    components.sort(key=lambda comp: _component_sort_key(comp, segments, chain_order))

    letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    if len(components) > len(letters):
        raise ValueError(f"Need {len(components)} output chains, but only {len(letters)} chain IDs are available.")

    # Determine output order within each component
    ordered_residue_blocks: List[Tuple[str, List[Tuple[str, int]]]] = []
    residue_to_new: Dict[Tuple[str, int], Tuple[str, int]] = {}
    output_links: List[str] = []
    serial = 1

    for comp_idx, comp in enumerate(components):
        new_chain = letters[comp_idx]
        # Determine directed order of segments in this component
        indeg = {sid: 0 for sid in comp}
        outdeg = {sid: 0 for sid in comp}
        for sid in comp:
            dst = seg_out[sid]
            if dst != -1 and dst in indeg:
                indeg[dst] += 1
                outdeg[sid] += 1

        starts = [sid for sid in comp if indeg[sid] == 0]
        ordered_seg_ids: List[int] = []
        if starts:
            if len(starts) > 1 and len(comp) > 1:
                raise ValueError(
                    "The selected LINK records create a branched component that cannot be serialized as a single chain."
                )
            cur = sorted(starts, key=lambda s: (segments[s]["chainID"], segments[s]["residues"][0], segments[s]["residues"][-1]))[0]
            seen_local: Set[int] = set()
            while cur != -1 and cur not in seen_local:
                ordered_seg_ids.append(cur)
                seen_local.add(cur)
                cur = seg_out[cur] if seg_out[cur] in comp else -1
        else:
            # cycle: choose stable start
            cur = sorted(comp, key=lambda s: (segments[s]["chainID"], segments[s]["residues"][0], segments[s]["residues"][-1]))[0]
            seen_local: Set[int] = set()
            while cur not in seen_local:
                ordered_seg_ids.append(cur)
                seen_local.add(cur)
                cur = seg_out[cur]
                if cur not in comp:
                    break

        if not ordered_seg_ids:
            ordered_seg_ids = sorted(comp)

        ordered_residues: List[Tuple[str, int]] = []
        for sid in ordered_seg_ids:
            seg = segments[sid]
            for res in seg["residues"]:
                ordered_residues.append((seg["chainID"], res))

        ordered_residue_blocks.append((new_chain, ordered_residues))

        # Rebuild atoms for this chain
        new_resseq = 1
        last_res: Optional[Tuple[str, int]] = None
        chain_last_atom: Optional[AtomInfo] = None
        chain_lines: List[str] = []

        for old_chain, old_res in ordered_residues:
            key = (old_chain, old_res)
            residue_to_new[key] = (new_chain, new_resseq)
            for atom in sorted(segments[residue_to_segment[key]]["atoms"][old_res], key=lambda a: a.serial if a.serial else 0):
                chain_lines.append(_rewrite_atom_line(atom, new_chain, new_resseq, serial))
                serial += 1
                chain_last_atom = atom
            new_resseq += 1
            last_res = key

        # TER for this chain (use last residue info)
        if last_res is not None and chain_last_atom is not None:
            orig_atom = chain_last_atom
            ter_line = _format_ter_line(serial, orig_atom.resName, new_chain, new_resseq - 1)
            serial += 1
            chain_lines.append(ter_line)

        # store back in component order as explicit lines
        ordered_residue_blocks[-1] = (new_chain, ordered_residues)
        ordered_residue_blocks.append(("__LINES__", chain_lines))

    # Build LINK lines using the new numbering map. Existing input LINK lines are
    # preserved first, with endpoints remapped whenever both endpoint residues
    # still exist in the rebuilt atom table.
    link_lines: List[str] = []
    written_link_keys: Set[frozenset[Tuple[str, int, str]]] = set()
    for existing in existing_link_entries:
        if isinstance(existing, ExistingLinkRecord):
            line, key = _remap_existing_link_line(existing, residue_to_new)
            if key in written_link_keys:
                continue
            link_lines.append(line)
            written_link_keys.add(key)
        else:
            link_lines.append(str(existing))

    for link in selected_links:
        endpoint1_new = residue_to_new[(link.endpoint1[0], link.endpoint1[1])]
        endpoint2_new = residue_to_new[(link.endpoint2[0], link.endpoint2[1])]
        selected_key = _link_endpoint_pair_key(
            endpoint1_new[0],
            endpoint1_new[1],
            link.endpoint1_atom.name,
            endpoint2_new[0],
            endpoint2_new[1],
            link.endpoint2_atom.name,
        )
        if selected_key in written_link_keys:
            continue
        link_lines.append(
            format_link_line(
                atom1_name=link.endpoint1_atom.name,
                resName1=link.endpoint1_atom.resName,
                chainID1=endpoint1_new[0],
                resSeq1=endpoint1_new[1],
                atom2_name=link.endpoint2_atom.name,
                resName2=link.endpoint2_atom.resName,
                chainID2=endpoint2_new[0],
                resSeq2=endpoint2_new[1],
                dist=distance(link.endpoint1_atom, link.endpoint2_atom),
            )
        )
        written_link_keys.add(selected_key)

    # Flatten output: prefix, LINKs, rebuilt atoms/TERs, suffix
    rebuilt_lines: List[str] = []
    rebuilt_lines.extend(prefix_lines)
    rebuilt_lines.extend(link_lines)
    for tag, payload in ordered_residue_blocks:
        if tag == "__LINES__":
            rebuilt_lines.extend(payload)
    rebuilt_lines.extend(suffix_lines)
    return "".join(rebuilt_lines)

# ----------------------------
# GUI
# ----------------------------


def run_gui(initial_path: Optional[str] = None) -> None:
    try:
        import tkinter as tk
        from tkinter import filedialog, messagebox
        from tkinter import ttk
    except Exception as e:
        print("[ERROR] Tkinter is required for --gui mode but could not be imported.")
        print(f"        {e}")
        sys.exit(1)

    class ScrollableFrame(ttk.Frame):
        def __init__(self, container: tk.Widget):
            super().__init__(container)
            canvas = tk.Canvas(self)
            scrollbar = ttk.Scrollbar(self, orient="vertical", command=canvas.yview)
            self.scrollable_frame = ttk.Frame(canvas)

            self.scrollable_frame.bind(
                "<Configure>",
                lambda e: canvas.configure(scrollregion=canvas.bbox("all")),
            )

            canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")
            canvas.configure(yscrollcommand=scrollbar.set)

            canvas.pack(side="left", fill="both", expand=True)
            scrollbar.pack(side="right", fill="y")

    class App(tk.Tk):
        def __init__(self):
            super().__init__()
            self.title(f"{TOOL_NAME} V{VERSION}")
            apply_optional_icon(self, __file__)
            self.geometry("1000x820")

            self.inp_path_var = tk.StringVar(value=initial_path or "")
            self.out_path_var = tk.StringVar(value="")
            self.auto_link_mode_var = tk.StringVar(value="nucleic")
            self.manual_link_mode_var = tk.StringVar(value="nucleic")

            self.manual_endpoint1_res_var = tk.StringVar(value="")
            self.manual_endpoint1_chain_var = tk.StringVar(value="")
            self.manual_endpoint2_res_var = tk.StringVar(value="")
            self.manual_endpoint2_chain_var = tk.StringVar(value="")
            self.manual_endpoint1_label_var = tk.StringVar(value="P (5') atom of residue")
            self.manual_endpoint2_label_var = tk.StringVar(value="O3' (3') atom of residue")

            self.current_lines: List[str] = []
            self.current_atoms: List[AtomInfo] = []
            self.chain_info: Dict[str, ChainEndInfo] = {}
            self.auto_chain_vars: Dict[str, tk.BooleanVar] = {}
            self.manual_items: List[PendingManualLink] = []

            self._build_widgets()

            if initial_path:
                self._load_file(Path(initial_path))

        def _build_widgets(self) -> None:
            pad = 6

            style = ttk.Style(self)
            style.configure("ToolTitle.TLabel", font=("TkDefaultFont", 12, "bold"))
            style.configure("Tool.TLabelframe.Label", font=("TkDefaultFont", 10, "bold"))

            ttk.Label(
                self,
                text=f"{TOOL_NAME} V{VERSION}",
                style="ToolTitle.TLabel",
            ).pack(anchor="w", padx=pad, pady=(pad, 0))

            row1 = ttk.Frame(self)
            row1.pack(fill="x", padx=pad, pady=pad)

            ttk.Label(row1, text="Input PDB:").pack(side="left")
            ent = ttk.Entry(row1, textvariable=self.inp_path_var, width=82)
            ent.pack(side="left", padx=(pad, pad), fill="x", expand=True)
            ttk.Button(row1, text="Browse...", command=self._browse).pack(side="left")

            row2 = ttk.Frame(self)
            row2.pack(fill="x", padx=pad, pady=(0, pad))
            ttk.Label(row2, text="Output:").pack(side="left")
            out_ent = ttk.Entry(row2, textvariable=self.out_path_var, width=82, state="readonly")
            out_ent.pack(side="left", padx=(pad, pad), fill="x", expand=True)

            automatic_box = ttk.LabelFrame(
                self,
                text="Automatic terminal LINK creation",
                style="Tool.TLabelframe",
            )
            automatic_box.pack(fill="x", padx=pad, pady=(0, pad))

            mode_row = ttk.Frame(automatic_box)
            mode_row.pack(fill="x", padx=pad, pady=(pad, 2))
            ttk.Label(mode_row, text="Molecule type:").pack(side="left")
            ttk.Radiobutton(
                mode_row,
                text="Nucleic acid (5' P to 3' O3')",
                variable=self.auto_link_mode_var,
                value="nucleic",
                command=self._on_auto_mode_changed,
            ).pack(side="left", padx=(8, 12))
            ttk.Radiobutton(
                mode_row,
                text="Peptide (N-terminus N to C-terminus C)",
                variable=self.auto_link_mode_var,
                value="peptide",
                command=self._on_auto_mode_changed,
            ).pack(side="left")

            ttk.Label(
                automatic_box,
                text="Check the chains to cyclize; the terminal-atom distance is shown below.",
            ).pack(anchor="w", padx=pad, pady=(2, 2))

            self.auto_scroll = ScrollableFrame(automatic_box)
            self.auto_scroll.pack(fill="both", expand=False, padx=pad, pady=(0, pad))

            manual_box = ttk.LabelFrame(
                self,
                text="Manual LINK creation (stage a LINK first)",
                style="Tool.TLabelframe",
            )
            manual_box.pack(fill="x", padx=pad, pady=(0, pad))

            manual_mode_row = ttk.Frame(manual_box)
            manual_mode_row.pack(fill="x", padx=pad, pady=(pad, 4))
            ttk.Label(manual_mode_row, text="Link type:").pack(side="left")
            ttk.Radiobutton(
                manual_mode_row,
                text="Nucleic acid (P to O3')",
                variable=self.manual_link_mode_var,
                value="nucleic",
                command=self._on_manual_mode_changed,
            ).pack(side="left", padx=(8, 12))
            ttk.Radiobutton(
                manual_mode_row,
                text="Peptide (N to C)",
                variable=self.manual_link_mode_var,
                value="peptide",
                command=self._on_manual_mode_changed,
            ).pack(side="left")

            grid = ttk.Frame(manual_box)
            grid.pack(fill="x", padx=pad, pady=(0, pad))

            ttk.Label(grid, textvariable=self.manual_endpoint1_label_var).grid(row=0, column=0, sticky="w")
            ttk.Entry(grid, textvariable=self.manual_endpoint1_res_var, width=10).grid(row=0, column=1, sticky="w", padx=(4, 8))
            ttk.Label(grid, text="in chain").grid(row=0, column=2, sticky="w")
            ttk.Entry(grid, textvariable=self.manual_endpoint1_chain_var, width=8).grid(row=0, column=3, sticky="w", padx=(4, 16))

            ttk.Label(grid, textvariable=self.manual_endpoint2_label_var).grid(row=1, column=0, sticky="w", pady=(6, 0))
            ttk.Entry(grid, textvariable=self.manual_endpoint2_res_var, width=10).grid(row=1, column=1, sticky="w", padx=(4, 8), pady=(6, 0))
            ttk.Label(grid, text="in chain").grid(row=1, column=2, sticky="w", pady=(6, 0))
            ttk.Entry(grid, textvariable=self.manual_endpoint2_chain_var, width=8).grid(row=1, column=3, sticky="w", padx=(4, 16), pady=(6, 0))

            ttk.Button(grid, text="Add Link", command=self._add_manual_link).grid(
                row=0, column=4, rowspan=2, sticky="ns", padx=(8, 0), pady=(2, 0)
            )
            grid.columnconfigure(5, weight=1)

            ttk.Label(
                self,
                text="Staged manual LINK records (default checked; uncheck any you do not want to write):",
            ).pack(anchor="w", padx=pad, pady=(2, 2))

            self.manual_scroll = ScrollableFrame(self)
            self.manual_scroll.pack(fill="both", expand=True, padx=pad, pady=(0, pad))

            row3 = ttk.Frame(self)
            row3.pack(fill="x", padx=pad, pady=(0, pad))

            ttk.Button(row3, text="Select all automatic", command=self._select_all_auto).pack(side="left")
            ttk.Button(row3, text="Clear automatic", command=self._clear_all_auto).pack(side="left", padx=(pad, 0))
            ttk.Button(row3, text="Clear manual", command=self._clear_manual).pack(side="left", padx=(pad, 0))
            ttk.Button(row3, text="Run", command=self._run).pack(side="right")
            ttk.Button(row3, text="Quit", command=self.destroy).pack(side="right", padx=(0, pad))

        def _browse(self) -> None:
            path = filedialog.askopenfilename(
                title="Select PDB file",
                filetypes=[("PDB files", "*.pdb"), ("All files", "*.*")],
            )
            if not path:
                return
            self.inp_path_var.set(path)
            self._load_file(Path(path))

        def _load_file(self, path: Path) -> None:
            try:
                self.current_lines = path.read_text(errors="replace").splitlines(True)
                self.current_atoms = [a for a in (parse_atom_line(ln) for ln in self.current_lines) if a is not None]
                self.chain_info = build_chain_end_info(self.current_atoms)
                if not self.chain_info:
                    raise ValueError("No ATOM/HETATM records found (or could not parse any).")

                self.out_path_var.set(str(self._default_output_for_mode(path)))
                self.manual_items = []
                self._populate_auto_chain_checkboxes()
                self._populate_manual_items()
            except Exception as e:
                messagebox.showerror("Error", f"Failed to read/parse file:\n{path}\n\n{e}")
                self.current_lines = []
                self.current_atoms = []
                self.chain_info = {}
                self.auto_chain_vars = {}
                self.manual_items = []
                self.out_path_var.set("")
                self._populate_auto_chain_checkboxes()
                self._populate_manual_items()

        def _default_output_for_mode(self, path: Path) -> Path:
            if self.auto_link_mode_var.get() == "peptide":
                return default_peptide_output_path(path)
            return default_output_path(path)

        def _on_auto_mode_changed(self) -> None:
            self._populate_auto_chain_checkboxes()
            input_text = self.inp_path_var.get().strip()
            if input_text:
                self.out_path_var.set(str(self._default_output_for_mode(Path(input_text))))

        def _on_manual_mode_changed(self) -> None:
            if self.manual_link_mode_var.get() == "peptide":
                self.manual_endpoint1_label_var.set("N atom of residue")
                self.manual_endpoint2_label_var.set("C atom of residue")
            else:
                self.manual_endpoint1_label_var.set("P (5') atom of residue")
                self.manual_endpoint2_label_var.set("O3' (3') atom of residue")

        def _populate_auto_chain_checkboxes(self) -> None:
            for w in self.auto_scroll.scrollable_frame.winfo_children():
                w.destroy()
            self.auto_chain_vars = {}

            if not self.chain_info:
                ttk.Label(self.auto_scroll.scrollable_frame, text="(no chains loaded)").pack(anchor="w")
                return

            for chain in sorted(self.chain_info.keys()):
                info = self.chain_info[chain]
                peptide_mode = self.auto_link_mode_var.get() == "peptide"
                endpoint1_atom = info.n_atom if peptide_mode else info.p_atom
                endpoint2_atom = info.c_atom if peptide_mode else info.o3_atom
                terminal_distance = info.peptide_distance if peptide_mode else info.distance
                var = tk.BooleanVar(value=endpoint1_atom is not None and endpoint2_atom is not None)
                self.auto_chain_vars[chain] = var

                chain_label = chain if chain.strip() else "<blank>"
                dist_str = "N/A" if terminal_distance is None else f"{terminal_distance:7.2f}"

                missing_bits: List[str] = []
                if endpoint1_atom is None:
                    missing_bits.append("missing N(N-term)" if peptide_mode else "missing P(5')")
                if endpoint2_atom is None:
                    missing_bits.append("missing C(C-term)" if peptide_mode else "missing O3'(3')")
                missing = "; ".join(missing_bits)
                if missing:
                    missing = "  [" + missing + "]"

                if peptide_mode:
                    txt = (
                        f"Chain {chain_label:>6s} : N-term {info.first_resName}{info.first_resSeq}  "
                        f"C-term {info.last_resName}{info.last_resSeq}   "
                        f"|N - C| = {dist_str} Å{missing}"
                    )
                else:
                    txt = (
                        f"Chain {chain_label:>6s} : 5' {info.first_resName}{info.first_resSeq}  "
                        f"3' {info.last_resName}{info.last_resSeq}   "
                        f"|P - O3'| = {dist_str} Å{missing}"
                    )
                cb = ttk.Checkbutton(self.auto_scroll.scrollable_frame, text=txt, variable=var)
                cb.pack(anchor="w")

        def _populate_manual_items(self) -> None:
            for w in self.manual_scroll.scrollable_frame.winfo_children():
                w.destroy()

            if not self.manual_items:
                ttk.Label(self.manual_scroll.scrollable_frame, text="(no manual LINK records added yet)").pack(anchor="w")
                return

            for idx, item in enumerate(self.manual_items, start=1):
                cb = ttk.Checkbutton(self.manual_scroll.scrollable_frame, text=f"{idx}. {item.label}", variable=item.var)
                cb.pack(anchor="w")

        def _select_all_auto(self) -> None:
            for v in self.auto_chain_vars.values():
                v.set(True)

        def _clear_all_auto(self) -> None:
            for v in self.auto_chain_vars.values():
                v.set(False)

        def _clear_manual(self) -> None:
            for item in self.manual_items:
                item.var.set(False)

        def _get_input_and_output_paths(self) -> Tuple[Path, Path]:
            inp = self.inp_path_var.get().strip()
            if not inp:
                raise ValueError("Please choose an input PDB file.")
            inp_path = Path(inp)
            out_text = self.out_path_var.get().strip()
            out_path = Path(out_text) if out_text else default_output_path(inp_path)
            return inp_path, out_path

        def _add_manual_link(self) -> None:
            try:
                if not self.current_lines or not self.current_atoms:
                    raise ValueError("Please load a PDB file first.")

                endpoint1_resseq = int(self.manual_endpoint1_res_var.get().strip())
                endpoint1_chain = self.manual_endpoint1_chain_var.get().strip()
                endpoint2_resseq = int(self.manual_endpoint2_res_var.get().strip())
                endpoint2_chain = self.manual_endpoint2_chain_var.get().strip()
                peptide_mode = self.manual_link_mode_var.get() == "peptide"
                endpoint1_name = "N" if peptide_mode else "P"
                endpoint2_name = "C" if peptide_mode else "O3'"

                endpoint1_atom, endpoint2_atom, e1, e2 = _prepare_named_manual_link_spec(
                    lines=self.current_lines,
                    atoms=self.current_atoms,
                    endpoint1_chain=endpoint1_chain,
                    endpoint1_resseq=endpoint1_resseq,
                    endpoint1_name=endpoint1_name,
                    endpoint2_chain=endpoint2_chain,
                    endpoint2_resseq=endpoint2_resseq,
                    endpoint2_name=endpoint2_name,
                )

                label = (
                    f"{endpoint1_atom.name} atom: residue {endpoint1_atom.resSeq} "
                    f"in chain {endpoint1_atom.chainID} -> {endpoint2_atom.name} atom: "
                    f"residue {endpoint2_atom.resSeq} in chain {endpoint2_atom.chainID}"
                )
                self.manual_items.append(
                    PendingManualLink(
                        label=label,
                        endpoint1_atom=endpoint1_atom,
                        endpoint2_atom=endpoint2_atom,
                        endpoint1=e1,
                        endpoint2=e2,
                        var=tk.BooleanVar(value=True),
                    )
                )
                self.manual_endpoint1_res_var.set("")
                self.manual_endpoint1_chain_var.set("")
                self.manual_endpoint2_res_var.set("")
                self.manual_endpoint2_chain_var.set("")
                self._populate_manual_items()
            except Exception as e:
                messagebox.showerror("Error", str(e))

        def _build_selected_links(self) -> List[SelectedLinkSpec]:
            if not self.current_lines or not self.current_atoms:
                raise ValueError("Please load a PDB file first.")

            selected_links: List[SelectedLinkSpec] = []
            selected_endpoints: Set[Tuple[str, int, str]] = set()

            # Automatic terminal links
            for chain in sorted(self.auto_chain_vars.keys()):
                if not self.auto_chain_vars[chain].get():
                    continue

                info = self.chain_info.get(chain)
                if info is None:
                    continue
                peptide_mode = self.auto_link_mode_var.get() == "peptide"
                endpoint1_atom = info.n_atom if peptide_mode else info.p_atom
                endpoint2_atom = info.c_atom if peptide_mode else info.o3_atom
                endpoint1_name = "N" if peptide_mode else "P"
                endpoint2_name = "C" if peptide_mode else "O3'"
                endpoint1_label = "N-terminal N" if peptide_mode else "5' P"
                endpoint2_label = "C-terminal C" if peptide_mode else "3' O3'"
                if endpoint1_atom is None:
                    raise ValueError(
                        f"Chain {chain!r}: missing {endpoint1_label} atom on terminal "
                        f"residue {info.first_resSeq}."
                    )
                if endpoint2_atom is None:
                    raise ValueError(
                        f"Chain {chain!r}: missing {endpoint2_label} atom on terminal "
                        f"residue {info.last_resSeq}."
                    )

                e1 = _format_endpoint(chain, info.first_resSeq, endpoint1_name)
                e2 = _format_endpoint(chain, info.last_resSeq, endpoint2_name)
                if e1 in selected_endpoints or e2 in selected_endpoints:
                    raise ValueError(
                        f"Chain {chain!r}: this LINK conflicts with another selected LINK record. Please uncheck one of them."
                    )

                selected_links.append(
                    SelectedLinkSpec(
                        label=(
                            f"Automatic peptide N-C cyclization: chain {chain}"
                            if peptide_mode
                            else f"Automatic nucleic-acid P-O3' circularization: chain {chain}"
                        ),
                        endpoint1_atom=endpoint1_atom,
                        endpoint2_atom=endpoint2_atom,
                        endpoint1=e1,
                        endpoint2=e2,
                    )
                )
                selected_endpoints.update((e1, e2))

            # Manually added LINK records
            for item in self.manual_items:
                if not item.var.get():
                    continue

                if item.endpoint1 in selected_endpoints or item.endpoint2 in selected_endpoints:
                    raise ValueError(
                        "One of the selected manual LINK records conflicts with another selected LINK record. "
                        "Please uncheck one of them."
                    )

                selected_links.append(
                    SelectedLinkSpec(
                        label=item.label,
                        endpoint1_atom=item.endpoint1_atom,
                        endpoint2_atom=item.endpoint2_atom,
                        endpoint1=item.endpoint1,
                        endpoint2=item.endpoint2,
                    )
                )
                selected_endpoints.update((item.endpoint1, item.endpoint2))

            return selected_links

        def _run(self) -> None:
            try:
                inp_path, out_path = self._get_input_and_output_paths()
                new_links = self._build_selected_links()
                auto_chains = [
                    chain for chain in sorted(self.auto_chain_vars.keys()) if self.auto_chain_vars[chain].get()
                ]
                manual_count = sum(1 for item in self.manual_items if item.var.get())

                if manual_count == 0 and auto_chains and all(chain.strip() for chain in auto_chains):
                    chain_option = (
                        "--peptide-chains"
                        if self.auto_link_mode_var.get() == "peptide"
                        else "--chains"
                    )
                    cli_parts = [
                        sys.executable,
                        Path(__file__).name,
                        str(inp_path),
                        chain_option,
                        *auto_chains,
                        "-o",
                        str(out_path),
                    ]
                    print("Related CLI command for automatic terminal LINK insertion:", flush=True)
                    print(" ".join(shlex.quote(str(part)) for part in cli_parts), flush=True)
                else:
                    print("GUI operation summary:", flush=True)
                    print(f"Input PDB: {inp_path}", flush=True)
                    print(f"Output PDB: {out_path}", flush=True)
                    if manual_count:
                        print(
                            "No single CLI command is available for staged manual GUI LINK records.",
                            flush=True,
                        )
                    elif auto_chains:
                        mode_label = (
                            "Peptide N-C cyclization"
                            if self.auto_link_mode_var.get() == "peptide"
                            else "Nucleic-acid P-O3' circularization"
                        )
                        print(f"{mode_label} chains: " + ", ".join(auto_chains), flush=True)

                print("Selected LINK records:", flush=True)
                for idx, link in enumerate(new_links, start=1):
                    print(f"{idx}. {link.label}: {link.to_link_line().strip()}", flush=True)

                out_text = build_relinked_pdb_text(self.current_lines, self.current_atoms, new_links)
                out_path.write_text(out_text)
                self.out_path_var.set(str(out_path))
                print(f"Wrote: {out_path}", flush=True)
                print(f"Added {len(new_links)} LINK record(s).", flush=True)
                messagebox.showinfo("Done", f"Wrote: {out_path}\n\nAdded {len(new_links)} LINK record(s).")
            except Exception as e:
                print(f"Error: {e}", flush=True)
                messagebox.showerror("Error", f"Failed to write PDB:\n\n{e}")

    app = App()
    app.mainloop()


# ----------------------------
# CLI
# ----------------------------


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="add_pdb_link_record.py",
        description=(
            "Add PDB LINK records for nucleic-acid P-O3' or peptide N-C "
            "cyclization, with GUI support for automatic and manual links."
        ),
    )
    p.add_argument("pdb", nargs="?", help="Input PDB file.")
    automatic_group = p.add_mutually_exclusive_group()
    automatic_group.add_argument(
        "-c",
        "--chains",
        nargs="*",
        default=None,
        help=(
            "Chain IDs to circularize. If omitted, circularize all chains. "
            "Accepts space-separated and/or comma-separated values, e.g. --chains A B or --chains A,B"
        ),
    )
    automatic_group.add_argument(
        "--peptide-chains",
        nargs="*",
        default=None,
        help=(
            "Peptide chain IDs to cyclize by linking the N-terminal N to the "
            "C-terminal C. If supplied without IDs, process all usable chains."
        ),
    )
    p.add_argument("-o", "--output", default=None, help="Output PDB filename (optional).")
    p.add_argument(
        "--gui",
        action="store_true",
        help="Launch the GUI (also the default if you run with no arguments).",
    )
    p.add_argument("-q", "--quiet", action="store_true", help="Suppress console messages.")
    p.add_argument("-v", "--version", action="version", version=f"{TOOL_NAME} V{VERSION}")
    return p


def main(argv: Optional[Sequence[str]] = None) -> None:
    if argv is None:
        argv = sys.argv[1:]

    if len(argv) == 0:
        run_gui()
        return

    parser = build_arg_parser()
    args = parser.parse_args(list(argv))

    if args.gui:
        run_gui(initial_path=args.pdb)
        return

    if not args.pdb:
        parser.print_help(sys.stderr)
        sys.exit(2)

    inp_path = Path(args.pdb)
    peptide_mode = args.peptide_chains is not None
    chain_args = args.peptide_chains if peptide_mode else args.chains
    chains = _parse_chain_list(chain_args)
    out_path = Path(args.output) if args.output else None
    verbose = not args.quiet

    link_function = cyclize_peptide_pdb if peptide_mode else circularize_pdb
    out_path2, chain_info, new_links = link_function(
        inp_path, chains=chains if chains else None, out_path=out_path, verbose=verbose
    )

    if verbose:
        print("\nDetected chain ends:")
        for chain in sorted(chain_info.keys()):
            info = chain_info[chain]
            chain_label = chain if chain.strip() else "<blank>"
            if peptide_mode:
                dist_str = (
                    "N/A" if info.peptide_distance is None else f"{info.peptide_distance:.2f}"
                )
                print(
                    f"  Chain {chain_label}: N-term {info.first_resName}{info.first_resSeq}  "
                    f"C-term {info.last_resName}{info.last_resSeq}  |N-C|={dist_str} Å"
                )
            else:
                dist_str = "N/A" if info.distance is None else f"{info.distance:.2f}"
                print(
                    f"  Chain {chain_label}: 5' {info.first_resName}{info.first_resSeq}  "
                    f"3' {info.last_resName}{info.last_resSeq}  |P-O3'|={dist_str} Å"
                )

        print(f"\nAdded {len(new_links)} LINK record(s).")
        print(f"Wrote output: {out_path2}")


if __name__ == "__main__":
    main()
