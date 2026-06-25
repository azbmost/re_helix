#!/usr/bin/env python3
"""add_pdb_link_record.py

Build LINK-driven chain topology from selected endpoints in a PDB file.

Supported actions
-----------------
1) Automatic circularization
   - For each selected chain, the script finds:
       P atom on the 5'-terminal residue (lowest resSeq)
       O3' (or O3*) atom on the 3'-terminal residue (highest resSeq)
   - The selected chain is added as a LINK-driven topology edit.

2) Manual LINK creation in the GUI
   - The user can specify:
       P atom of residue [resSeq] in chain [chainID]
       O3' atom of residue [resSeq] in chain [chainID]
   - Clicking "Add Link" stages a default-checked LINK entry below.
   - Clicking "Run" rebuilds the topology from the checked automatic and
     manual LINKs in one pass.

Behavior
--------
- Selected endpoints are used to split chain segments, merge them into new
  output chains, reassign chain IDs, rewrite TER records, and renumber residues.
- Existing LINK lines are preserved, remapped to rebuilt residue numbering when
  possible, and written before any newly selected LINK edits.
- Output filename is the input name with "_circ" inserted before the final
  extension, e.g.:
      input.pdb -> input_circ.pdb

Important notes
---------------
- 5' and 3' ends are inferred from residue numbering (min/max resSeq) within
  each chain for automatic circularization.
- O3' atom name may appear as "O3'" or "O3*"; both are supported.
- The output topology assumes the selected P/O3' endpoints are valid chain
  endpoints for the intended segment merge.

GUI
---
- If run with no arguments, or with --gui, a small Tk GUI is launched.
- After selecting a PDB, the GUI reports the distance between:
      P of the 5' residue (lowest resSeq)
      and O3' (or O3*) of the 3' residue (highest resSeq)
  for every chain, to help decide which chain(s) to circularize.
- The GUI also includes a manual LINK creation panel and a pending LINK list.

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

TOOL_NAME = "Add PDB LINK Record"
VERSION = "1.1"


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


@dataclass
class PendingManualLink:
    label: str
    p_atom: AtomInfo
    o3_atom: AtomInfo
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
    """For each chain, identify min/max resSeq and get P(5') and O3'(3') atoms."""
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

        o3_atom: Optional[AtomInfo] = None
        for nm in O3_CANDIDATES:
            key = nm.upper()
            if key in last_atoms:
                o3_atom = last_atoms[key]
                break

        dist: Optional[float] = None
        if p_atom is not None and o3_atom is not None:
            dist = distance(p_atom, o3_atom)

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
        )

    return out


def default_output_path(inp_path: Path) -> Path:
    """Insert '_circ' before the final suffix (extension)."""
    return inp_path.with_name(f"{inp_path.stem}_circ{inp_path.suffix}")


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


def circularize_pdb(
    inp_path: Path,
    chains: Optional[Sequence[str]] = None,
    out_path: Optional[Path] = None,
    verbose: bool = True,
) -> Tuple[Path, Dict[str, ChainEndInfo], List[str]]:
    """Circularize selected chains by inserting LINK records.

    Returns:
        out_path, chain_info, new_link_lines
    """
    inp_path = Path(inp_path)
    if out_path is None:
        out_path = default_output_path(inp_path)

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

        if info.p_atom is None:
            if verbose:
                print(f"[WARN] Chain {chain!r}: missing P atom on 5' residue {info.first_resSeq}; skipping.")
            continue
        if info.o3_atom is None:
            if verbose:
                print(
                    f"[WARN] Chain {chain!r}: missing O3' (or O3*) atom on 3' residue {info.last_resSeq}; skipping."
                )
            continue

        e1 = _format_endpoint(chain, info.first_resSeq, "P")
        e2 = _format_endpoint(chain, info.last_resSeq, "O3'")

        key = frozenset((e1, e2))
        if key in existing_links:
            if verbose:
                print(f"[INFO] Chain {chain!r}: circularization LINK already present; not adding duplicate.")
            continue

        # Prevent duplicate use of either endpoint in another LINK record
        if e1 in existing_link_endpoints:
            raise ValueError(
                f"Chain {chain!r}: the 5' P atom on residue {info.first_resSeq} is already used in an existing LINK record."
            )
        if e2 in existing_link_endpoints:
            raise ValueError(
                f"Chain {chain!r}: the 3' O3' atom on residue {info.last_resSeq} is already used in an existing LINK record."
            )

        link_line = format_link_line(
            atom1_name=info.p_atom.name,
            resName1=info.p_atom.resName,
            chainID1=info.p_atom.chainID,
            resSeq1=info.p_atom.resSeq,
            atom2_name=info.o3_atom.name,
            resName2=info.o3_atom.resName,
            chainID2=info.o3_atom.chainID,
            resSeq2=info.o3_atom.resSeq,
            dist=info.distance,
        )
        new_links.append(link_line)

    out_lines = _insert_link_lines(lines, new_links)
    out_path.write_text("".join(out_lines))

    return out_path, chain_info, new_links


def _prepare_manual_link_spec(
    lines: Sequence[str],
    atoms: Sequence[AtomInfo],
    p_chain: str,
    p_resseq: int,
    o3_chain: str,
    o3_resseq: int,
) -> Tuple[AtomInfo, AtomInfo, Tuple[str, int, str], Tuple[str, int, str]]:
    """Validate and build a manual LINK record spec without writing output.

    This version allows internal residues too; the output stage will split
    segments at the selected endpoints and rebuild chain IDs / TER records.
    """
    atom_index = build_atom_index(atoms)

    p_chain = p_chain.strip()
    o3_chain = o3_chain.strip()

    if p_chain == "":
        raise ValueError("Please enter a valid chain ID for the P atom.")
    if o3_chain == "":
        raise ValueError("Please enter a valid chain ID for the O3' atom.")

    try:
        p_res_map = atom_index[p_chain]
    except KeyError:
        raise ValueError(f"Chain {p_chain!r} was not found in the PDB file.")

    try:
        o3_res_map = atom_index[o3_chain]
    except KeyError:
        raise ValueError(f"Chain {o3_chain!r} was not found in the PDB file.")

    if p_resseq not in p_res_map:
        raise ValueError(f"Residue {p_resseq} was not found in chain {p_chain!r}.")
    if o3_resseq not in o3_res_map:
        raise ValueError(f"Residue {o3_resseq} was not found in chain {o3_chain!r}.")

    p_atom = p_res_map[p_resseq].get("P")
    if p_atom is None:
        raise ValueError(f"No P atom was found on residue {p_resseq} in chain {p_chain!r}.")

    o3_atom = None
    for nm in O3_CANDIDATES:
        key = nm.upper()
        if key in o3_res_map[o3_resseq]:
            o3_atom = o3_res_map[o3_resseq][key]
            break
    if o3_atom is None:
        raise ValueError(f"No O3' (or O3*) atom was found on residue {o3_resseq} in chain {o3_chain!r}.")

    e1 = _format_endpoint(p_chain, p_resseq, "P")
    e2 = _format_endpoint(o3_chain, o3_resseq, "O3'")
    return p_atom, o3_atom, e1, e2


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
            p_atom=p_atom,
            o3_atom=o3_atom,
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
    p_atom: AtomInfo
    o3_atom: AtomInfo
    endpoint1: Tuple[str, int, str]
    endpoint2: Tuple[str, int, str]

    def to_link_line(self) -> str:
        return format_link_line(
            atom1_name=self.p_atom.name,
            resName1=self.p_atom.resName,
            chainID1=self.endpoint1[0],
            resSeq1=self.endpoint1[1],
            atom2_name=self.o3_atom.name,
            resName2=self.o3_atom.resName,
            chainID2=self.endpoint2[0],
            resSeq2=self.endpoint2[1],
            dist=distance(self.p_atom, self.o3_atom),
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
        p_chain, p_res, _ = link.endpoint1
        o_chain, o_res, _ = link.endpoint2
        cut_before.setdefault(p_chain, set()).add(p_res)
        cut_after.setdefault(o_chain, set()).add(o_res)

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
        p_chain, p_res, _ = link.endpoint1
        o_chain, o_res, _ = link.endpoint2
        s_src = residue_to_segment.get((o_chain, o_res))
        s_tgt = residue_to_segment.get((p_chain, p_res))
        if s_src is None or s_tgt is None:
            raise ValueError(
                f"Could not place selected LINK endpoints into segments: {o_res}{o_chain} -> {p_res}{p_chain}."
            )

        if segments[s_src]["residues"][-1] != o_res:
            raise ValueError(
                f"The selected O3' endpoint {o_res}{o_chain} is not the terminal residue of its segment."
            )
        if segments[s_tgt]["residues"][0] != p_res:
            raise ValueError(
                f"The selected P endpoint {p_res}{p_chain} is not the terminal residue of its segment."
            )

        if seg_out[s_src] not in (-1, s_tgt):
            raise ValueError(
                f"Segment starting at {o_res}{o_chain} already connects to another selected LINK."
            )
        if seg_in[s_tgt] not in (-1, s_src):
            raise ValueError(
                f"Segment starting at {p_res}{p_chain} already receives another selected LINK."
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
        p_new = residue_to_new[(link.endpoint1[0], link.endpoint1[1])]
        o_new = residue_to_new[(link.endpoint2[0], link.endpoint2[1])]
        selected_key = _link_endpoint_pair_key(
            p_new[0],
            p_new[1],
            link.p_atom.name,
            o_new[0],
            o_new[1],
            link.o3_atom.name,
        )
        if selected_key in written_link_keys:
            continue
        link_lines.append(
            format_link_line(
                atom1_name=link.p_atom.name,
                resName1=link.p_atom.resName,
                chainID1=p_new[0],
                resSeq1=p_new[1],
                atom2_name=link.o3_atom.name,
                resName2=link.o3_atom.resName,
                chainID2=o_new[0],
                resSeq2=o_new[1],
                dist=distance(link.p_atom, link.o3_atom),
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
            self.geometry("1000x820")

            self.inp_path_var = tk.StringVar(value=initial_path or "")
            self.out_path_var = tk.StringVar(value="")

            self.manual_p_res_var = tk.StringVar(value="")
            self.manual_p_chain_var = tk.StringVar(value="")
            self.manual_o3_res_var = tk.StringVar(value="")
            self.manual_o3_chain_var = tk.StringVar(value="")

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

            ttk.Label(
                self,
                text=(
                    "Automatic circularization: check the chains you want to circularize "
                    "(distance shown is |P(5') - O3'(3')| in Angstrom)."
                ),
            ).pack(anchor="w", padx=pad, pady=(2, 2))

            self.auto_scroll = ScrollableFrame(self)
            self.auto_scroll.pack(fill="both", expand=False, padx=pad, pady=(0, pad))

            manual_box = ttk.LabelFrame(
                self,
                text="Manual LINK creation (stage a LINK first)",
                style="Tool.TLabelframe",
            )
            manual_box.pack(fill="x", padx=pad, pady=(0, pad))

            ttk.Label(
                manual_box,
                text=(
                    "P (5’) atom of residue [  ] in chain [  ]\n"
                    "O3' (3') atom of residue [  ] in chain [  ]"
                ),
            ).pack(anchor="w", padx=pad, pady=(pad, 4))

            grid = ttk.Frame(manual_box)
            grid.pack(fill="x", padx=pad, pady=(0, pad))

            ttk.Label(grid, text="P (5') atom of residue").grid(row=0, column=0, sticky="w")
            ttk.Entry(grid, textvariable=self.manual_p_res_var, width=10).grid(row=0, column=1, sticky="w", padx=(4, 8))
            ttk.Label(grid, text="in chain").grid(row=0, column=2, sticky="w")
            ttk.Entry(grid, textvariable=self.manual_p_chain_var, width=8).grid(row=0, column=3, sticky="w", padx=(4, 16))

            ttk.Label(grid, text="O3' (3') atom of residue").grid(row=1, column=0, sticky="w", pady=(6, 0))
            ttk.Entry(grid, textvariable=self.manual_o3_res_var, width=10).grid(row=1, column=1, sticky="w", padx=(4, 8), pady=(6, 0))
            ttk.Label(grid, text="in chain").grid(row=1, column=2, sticky="w", pady=(6, 0))
            ttk.Entry(grid, textvariable=self.manual_o3_chain_var, width=8).grid(row=1, column=3, sticky="w", padx=(4, 16), pady=(6, 0))

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

                self.out_path_var.set(str(default_output_path(path)))
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

        def _populate_auto_chain_checkboxes(self) -> None:
            for w in self.auto_scroll.scrollable_frame.winfo_children():
                w.destroy()
            self.auto_chain_vars = {}

            if not self.chain_info:
                ttk.Label(self.auto_scroll.scrollable_frame, text="(no chains loaded)").pack(anchor="w")
                return

            for chain in sorted(self.chain_info.keys()):
                info = self.chain_info[chain]
                var = tk.BooleanVar(value=True)
                self.auto_chain_vars[chain] = var

                chain_label = chain if chain.strip() else "<blank>"
                dist_str = "N/A" if info.distance is None else f"{info.distance:7.2f}"

                missing_bits: List[str] = []
                if info.p_atom is None:
                    missing_bits.append("missing P(5')")
                if info.o3_atom is None:
                    missing_bits.append("missing O3'(3')")
                missing = "; ".join(missing_bits)
                if missing:
                    missing = "  [" + missing + "]"

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

                p_resseq = int(self.manual_p_res_var.get().strip())
                p_chain = self.manual_p_chain_var.get().strip()
                o3_resseq = int(self.manual_o3_res_var.get().strip())
                o3_chain = self.manual_o3_chain_var.get().strip()

                p_atom, o3_atom, _e1, _e2 = _prepare_manual_link_spec(
                    lines=self.current_lines,
                    atoms=self.current_atoms,
                    p_chain=p_chain,
                    p_resseq=p_resseq,
                    o3_chain=o3_chain,
                    o3_resseq=o3_resseq,
                )

                label = (
                    f"P atom: residue {p_atom.resSeq} in chain {p_atom.chainID} -> "
                    f"O3' atom: residue {o3_atom.resSeq} in chain {o3_atom.chainID}"
                )
                self.manual_items.append(
                    PendingManualLink(
                        label=label,
                        p_atom=p_atom,
                        o3_atom=o3_atom,
                        endpoint1=_format_endpoint(p_atom.chainID, p_atom.resSeq, "P"),
                        endpoint2=_format_endpoint(o3_atom.chainID, o3_atom.resSeq, "O3'"),
                        var=tk.BooleanVar(value=True),
                    )
                )
                self.manual_p_res_var.set("")
                self.manual_p_chain_var.set("")
                self.manual_o3_res_var.set("")
                self.manual_o3_chain_var.set("")
                self._populate_manual_items()
            except Exception as e:
                messagebox.showerror("Error", str(e))

        def _build_selected_links(self) -> List[SelectedLinkSpec]:
            if not self.current_lines or not self.current_atoms:
                raise ValueError("Please load a PDB file first.")

            selected_links: List[SelectedLinkSpec] = []
            selected_endpoints: Set[Tuple[str, int, str]] = set()

            # Automatic circularization links
            for chain in sorted(self.auto_chain_vars.keys()):
                if not self.auto_chain_vars[chain].get():
                    continue

                info = self.chain_info.get(chain)
                if info is None:
                    continue
                if info.p_atom is None:
                    raise ValueError(
                        f"Chain {chain!r}: missing P atom on the 5' terminal residue {info.first_resSeq}."
                    )
                if info.o3_atom is None:
                    raise ValueError(
                        f"Chain {chain!r}: missing O3' atom on the 3' terminal residue {info.last_resSeq}."
                    )

                e1 = _format_endpoint(chain, info.first_resSeq, "P")
                e2 = _format_endpoint(chain, info.last_resSeq, "O3'")
                if e1 in selected_endpoints or e2 in selected_endpoints:
                    raise ValueError(
                        f"Chain {chain!r}: this LINK conflicts with another selected LINK record. Please uncheck one of them."
                    )

                selected_links.append(
                    SelectedLinkSpec(
                        label=f"Automatic circularization: chain {chain}",
                        p_atom=info.p_atom,
                        o3_atom=info.o3_atom,
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
                        p_atom=item.p_atom,
                        o3_atom=item.o3_atom,
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
                    cli_parts = [
                        sys.executable,
                        Path(__file__).name,
                        str(inp_path),
                        "--chains",
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
                        print("Automatic circularization chains: " + ", ".join(auto_chains), flush=True)

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
            "Add PDB LINK records between P and O3' endpoints, with GUI support "
            "for manual links and automatic circularization."
        ),
    )
    p.add_argument("pdb", nargs="?", help="Input PDB file.")
    p.add_argument(
        "-c",
        "--chains",
        nargs="*",
        default=None,
        help=(
            "Chain IDs to circularize. If omitted, circularize all chains. "
            "Accepts space-separated and/or comma-separated values, e.g. --chains A B or --chains A,B"
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
    chains = _parse_chain_list(args.chains)
    out_path = Path(args.output) if args.output else None
    verbose = not args.quiet

    out_path2, chain_info, new_links = circularize_pdb(
        inp_path,
        chains=chains if chains else None,
        out_path=out_path,
        verbose=verbose,
    )

    if verbose:
        print("\nDetected chain ends:")
        for chain in sorted(chain_info.keys()):
            info = chain_info[chain]
            chain_label = chain if chain.strip() else "<blank>"
            dist_str = "N/A" if info.distance is None else f"{info.distance:.2f}"
            print(
                f"  Chain {chain_label}: 5' {info.first_resName}{info.first_resSeq}  "
                f"3' {info.last_resName}{info.last_resSeq}  |P-O3'|={dist_str} Å"
            )

        print(f"\nAdded {len(new_links)} LINK record(s).")
        print(f"Wrote output: {out_path2}")


if __name__ == "__main__":
    main()
