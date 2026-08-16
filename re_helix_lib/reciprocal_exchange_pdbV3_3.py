#!/usr/bin/env python3
"""reciprocal_exchange_pdbV3_3.py

Reciprocal exchange (double / single) and bowtie exchange for DNA PDBs.

Key behavior (V3.6):
  1) Build the ORIGINAL 5'->3' backbone graph from residue numbering plus
     recognized existing LINK topology for each chain.
  2) Parse ALL exchanges and apply them as EDGE REWIRES on that original graph.
     - This removes the order-dependence bug when mixing kinds (double/single/bowtie)
       as long as the requested cut edges do not overlap.
  3) Bowtie handling:
     - For bowtie pos1/pos2, we cut the two original incoming backbone edges and add:
         * a 3'-3' special edge between the *predecessors* (prev(pos1), prev(pos2))
           that will later be expanded by inserting a phosphate-only residue derived
           from the phosphate group of pos2.
         * a 5'-5' special edge between pos1 (P atom) and pos2 (O5' atom).
     - Works even when pos1 and pos2 end up on the same final strand.
  4) Standalone 3'-3' linker phosphate residues:
     - For each bowtie, we cut the phosphate-group atoms (P + non-bridging O's)
       from pos2 and store them (keyed by the ORIGINAL pos2 label, e.g., 23F).
     - When reconstructing final strand paths, we insert a standalone
       phosphate-only residue between the two residues that form the 3'-3'
       edge. The default is HETATM X33, but callers can supply another residue
       name, or use ATOM DA for Phenix-friendly relaxation. The linker contains
       exactly three atoms: P, OP1, and OP2.
  5) LINK records:
     - Existing re_helix LINK records are used to reconstruct the input
       backbone graph, including inverted P--O3', 5'-5', and standalone
       phosphate-bridge connections. This makes generated output composable:
       it can safely be used as the input to a later reciprocal exchange.
     - We ignore CONECT entirely.
     - We write LINK records for:
         (a) every bowtie 3'-3' linkage: P(phosphate-only residue) -- O3' (each side)
         (b) every bowtie 5'-5' linkage: O5'(pos2) -- P(pos1)
         (c) every inverted backbone step (i.e., when traversing a standard O3'--P bond
             in the P->O3 direction). Natural O3->P steps rely on standard connectivity.
         (d) in c/C mode, the closing edge of every circular component, because
             it crosses the output TER boundary even when its orientation is
             natural O3'--P.

Additionally:
  - Cyclic components are supported. We perform a circular permutation (rotation)
    for output numbering, similar to V2/V3. cir_shift is applied exactly modulo
    the cycle length; zero preserves the canonical input-provenance start.
  - Path/cycle direction preserves the input PDB's serialized residue order.
    Path direction preserves the user-visible direction of its terminal input
    fragments before considering interior continuity. This keeps closely
    related products consistently oriented when their cut positions change.
  - LINK minimization is opt-in. With --min_link_records, predicted topology
    LINK count becomes the primary orientation criterion and may reverse an
    entire output strand. Cycles are scored after their exact cir_shift, with
    the requested open or covalently closed output mode.
  - A plain cir_shift leaves that serialized break open. Appending c/C to the
    shift writes the closing LINK and preserves covalent circularization.
  - For debugging, the script prints the final linking paths (using ORIGINAL labels)
    with arrows:
        '->'  : standard 5'->3' (no LINK required)
        '->>' : a bond that requires a LINK record (special or inverted)

Usage:
  python reciprocal_exchange_pdbV3_3.py input.pdb  9C 23A double  23C 23F B  9A 9F B  -o out.pdb

Notes:
  - Residue tokens are like "23A" (resSeq + chainID).
  - Kind tokens: double / single / bowtie; also accepts D/S/B.

"""

from __future__ import annotations

import argparse
import math
import re
import shlex
import sys
from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple

try:
    from . import edit_pdb_atom
    from . import edit_pdb_link
except ImportError:  # pragma: no cover - keeps direct script execution working.
    import edit_pdb_atom
    import edit_pdb_link


SOFTWARE_NAME = "reciprocal_exchange_pdb"
SOFTWARE_VERSION = "V3.6"
SOFTWARE_DEVELOPER = "DiLiuLab"
REMARK_PREFIX = "REMARK 950 RE_SCRIPT"
X33_HETID = "X33"
LINKER_PHOSPHATE_DA_RESNAME = "DA"
LINKER_PHOSPHATE_DEFAULT_RECORD = "HETATM"
LINKER_PHOSPHATE_HETNAM = "3'-3' PHOSPHODIESTER LINKER PHOSPHATE"
X33_HETNAM = LINKER_PHOSPHATE_HETNAM


@dataclass(frozen=True)
class LinkerPhosphateStyle:
    """Output style for phosphate-only residues inserted at bowtie 3'-3' links."""

    resname: str = X33_HETID
    record_name: str = LINKER_PHOSPHATE_DEFAULT_RECORD


@dataclass(frozen=True)
class CirShiftSpec:
    """Circular-strand serialization shift and optional closure-LINK mode."""

    shift: int = 8
    circularize: bool = False

    def __str__(self) -> str:
        return f"{self.shift}{'c' if self.circularize else ''}"


_CIR_SHIFT_RE = re.compile(r"^([+-]?\d+)([cC]?)$")


def parse_cir_shift(value: object) -> CirShiftSpec:
    """Parse an integer shift with an optional ``c`` cycle-closure suffix."""
    if isinstance(value, CirShiftSpec):
        return value
    if isinstance(value, bool):
        raise ValueError("cir_shift must be an integer optionally followed by c")
    if isinstance(value, int):
        return CirShiftSpec(shift=value, circularize=False)

    text = str(value).strip()
    match = _CIR_SHIFT_RE.fullmatch(text)
    if match is None:
        raise ValueError(
            f"Invalid cir_shift value {value!r}; use an integer such as 10 or "
            "append c/C to close circular strands with LINK records, such as 10c."
        )
    return CirShiftSpec(
        shift=int(match.group(1)),
        circularize=bool(match.group(2)),
    )


def coerce_cir_shift(value: object) -> CirShiftSpec:
    """Compatibility alias for callers that pass legacy integer shifts."""
    return parse_cir_shift(value)


def _normalize_linker_phosphate_resname(resname: Optional[str]) -> str:
    if resname is None or not str(resname).strip():
        return X33_HETID
    normalized = str(resname).strip().upper()
    if len(normalized) > 3:
        raise ValueError(
            f"3'-3' linker phosphate residue name '{resname}' is too long; "
            "PDB residue names must be 1-3 characters."
        )
    if not re.match(r"^[A-Z0-9]{1,3}$", normalized):
        raise ValueError(
            f"3'-3' linker phosphate residue name '{resname}' is invalid; "
            "use 1-3 letters/digits."
        )
    return normalized


def _normalize_linker_phosphate_record(record_name: Optional[str]) -> Optional[str]:
    if record_name is None or not str(record_name).strip():
        return None
    normalized = str(record_name).strip().upper()
    if normalized not in {"ATOM", "HETATM"}:
        raise ValueError("3'-3' linker phosphate record type must be ATOM or HETATM.")
    return normalized


def make_linker_phosphate_style(
    resname: Optional[str] = None,
    record_name: Optional[str] = None,
) -> LinkerPhosphateStyle:
    """Return the configured 3'-3' linker phosphate output style.

    Defaults to HETATM X33.  If the residue name is DA/dA and the record type is
    not explicitly supplied, use ATOM DA so refinement programs can treat the
    phosphate-only residue as a standard deoxyadenosine residue with missing
    atoms rather than a custom ligand.
    """
    normalized_resname = _normalize_linker_phosphate_resname(resname)
    normalized_record = _normalize_linker_phosphate_record(record_name)
    if normalized_record is None:
        normalized_record = "ATOM" if normalized_resname == LINKER_PHOSPHATE_DA_RESNAME else "HETATM"
    return LinkerPhosphateStyle(resname=normalized_resname, record_name=normalized_record)


def _coerce_linker_phosphate_style(
    linker_phosphate_style: Optional[object] = None,
) -> LinkerPhosphateStyle:
    if linker_phosphate_style is None:
        return make_linker_phosphate_style()
    if isinstance(linker_phosphate_style, LinkerPhosphateStyle):
        return linker_phosphate_style
    if isinstance(linker_phosphate_style, str):
        return make_linker_phosphate_style(linker_phosphate_style)
    raise TypeError("linker_phosphate_style must be a LinkerPhosphateStyle, residue-name string, or None")


def coerce_linker_phosphate_style(
    linker_phosphate_style: Optional[object] = None,
) -> LinkerPhosphateStyle:
    """Return a LinkerPhosphateStyle from a style object, residue name, or None."""
    return _coerce_linker_phosphate_style(linker_phosphate_style)


def _is_default_x33_style(linker_phosphate_style: LinkerPhosphateStyle) -> bool:
    return (
        linker_phosphate_style.resname == X33_HETID
        and linker_phosphate_style.record_name == LINKER_PHOSPHATE_DEFAULT_RECORD
    )


# -------------------------- Data structures --------------------------


Label = Tuple[str, int]  # (chainID, resSeq)


@dataclass
class ResidueNode:
    orig_chain_id: str
    orig_res_seq: int
    atoms: List[edit_pdb_atom.pdb_atom_record]

    # Bowtie bookkeeping
    is_phos_bridge: bool = False
    phos_source: Optional[Label] = None  # original pos2 label providing phosphate
    no_phosphate: bool = False           # true for bowtie pos2 residues after cutting P/OP1/OP2

    # Output labels
    new_chain_id: str = ""
    new_res_seq: int = 0

    # Zero-based position within the chain as serialized in the input PDB.
    # This remains after the historical positional fields for compatibility
    # with callers that construct ResidueNode positionally.
    input_chain_rank: int = 0

    def orig_label(self) -> Label:
        return (self.orig_chain_id, self.orig_res_seq)


@dataclass
class Edge:
    # Undirected edge between nodes, but with typed endpoints.
    a: int
    b: int
    kind: str  # 'std', '3to3', '5to5'
    end_a: str  # 'O3', 'P', 'O5'
    end_b: str
    phos_key: Optional[Label] = None  # only for kind=='3to3'

    def endpoints(self, u: int, v: int) -> Tuple[str, str]:
        """Return endpoint labels (end_u, end_v) for traversal u->v."""
        if u == self.a and v == self.b:
            return self.end_a, self.end_b
        if u == self.b and v == self.a:
            return self.end_b, self.end_a
        raise ValueError("Edge.endpoints called with non-incident nodes")


class BackboneGraph:
    """Degree-<=2 graph of residue nodes with typed edges."""

    def __init__(self, n_nodes: int):
        self.n_nodes = n_nodes
        self.edges: Dict[frozenset[int], Edge] = {}
        self.neigh: List[List[int]] = [[] for _ in range(n_nodes)]

    def add_edge(
        self,
        a: int,
        b: int,
        kind: str,
        end_a: str,
        end_b: str,
        phos_key: Optional[Label] = None,
    ) -> None:
        if a == b:
            raise ValueError("Self-edge is not allowed")
        k = frozenset((a, b))
        if k in self.edges:
            raise ValueError(f"Edge already exists between nodes {a} and {b}")
        self.edges[k] = Edge(a=a, b=b, kind=kind, end_a=end_a, end_b=end_b, phos_key=phos_key)
        self.neigh[a].append(b)
        self.neigh[b].append(a)
        if len(self.neigh[a]) > 2 or len(self.neigh[b]) > 2:
            raise ValueError(
                f"Invalid graph (degree>2) after adding edge {a}-{b}; "
                f"degrees are {len(self.neigh[a])}, {len(self.neigh[b])}"
            )

    def remove_edge(self, a: int, b: int) -> None:
        k = frozenset((a, b))
        if k not in self.edges:
            raise ValueError(f"Requested to remove missing edge between nodes {a} and {b}")
        del self.edges[k]
        self.neigh[a].remove(b)
        self.neigh[b].remove(a)

    def get_edge(self, a: int, b: int) -> Edge:
        return self.edges[frozenset((a, b))]


@dataclass
class InputBackboneTopology:
    """Backbone graph reconstructed from coordinate order and input LINKs."""

    graph: BackboneGraph
    orig_prev: List[Optional[int]]
    passthrough_links: List[edit_pdb_link.pdb_link_record]
    recognized_link_count: int = 0
    phosphate_bridge_count: int = 0


# -------------------------- Parsing helpers --------------------------


_RES_RE = re.compile(r"^(\d+)([A-Za-z])$")


def parse_res_label(token: str) -> Label:
    m = _RES_RE.match(token.strip())
    if not m:
        raise ValueError(f"Invalid residue label token: '{token}' (expected like 23A)")
    res_seq = int(m.group(1))
    # PDB chain identifiers are case-sensitive.  re_helix deliberately uses
    # lower-case IDs after exhausting A-Z, so preserve the user's exact ID.
    chain_id = m.group(2)
    return (chain_id, res_seq)


def normalize_kind(token: str) -> str:
    t = token.strip().lower()
    if t in {"double", "d"}:
        return "double"
    if t in {"single", "s"}:
        return "single"
    if t in {"bowtie", "b"}:
        return "bowtie"
    raise ValueError(f"Unrecognized exchange kind token: '{token}'")


def parse_exchange_specs(tokens: List[str]) -> List[Dict[str, object]]:
    """Parse a flat token list into exchange specs.

    Expected pattern: (pos1 pos2 kind) repeated.
    Example:  9C 23A double  23C 23F B  9A 9F B
    """
    specs: List[Dict[str, object]] = []
    i = 0
    while i < len(tokens):
        if i + 2 >= len(tokens):
            raise ValueError("Incomplete exchange specification; expected triples pos1 pos2 kind")
        pos1 = parse_res_label(tokens[i])
        pos2 = parse_res_label(tokens[i + 1])
        kind = normalize_kind(tokens[i + 2])
        specs.append({"pos1": pos1, "pos2": pos2, "kind": kind})
        i += 3
    return specs


# -------------------------- PDB residue parsing --------------------------


def build_residue_nodes(atom_recs: List[edit_pdb_atom.pdb_atom_record]) -> Tuple[List[ResidueNode], Dict[Label, int]]:
    """Group ATOM/HETATM records into residue nodes keyed by (chainID,resSeq)."""
    grouped: Dict[Label, List[edit_pdb_atom.pdb_atom_record]] = {}
    input_chain_rank: Dict[Label, int] = {}
    next_rank_by_chain: Dict[str, int] = {}
    for a in atom_recs:
        if a.recordName not in ("ATOM", "HETATM"):
            continue
        key = (a.chainID.strip() or " ", a.resSeq)
        if key not in grouped:
            input_chain_rank[key] = next_rank_by_chain.get(key[0], 0)
            next_rank_by_chain[key[0]] = input_chain_rank[key] + 1
        grouped.setdefault(key, []).append(a)

    nodes: List[ResidueNode] = []
    label_to_idx: Dict[Label, int] = {}
    for key in sorted(grouped.keys(), key=lambda x: (x[0], x[1])):
        chain_id, res_seq = key
        idx = len(nodes)
        nodes.append(
            ResidueNode(
                orig_chain_id=chain_id,
                orig_res_seq=res_seq,
                atoms=grouped[key],
                input_chain_rank=input_chain_rank[key],
            )
        )
        label_to_idx[key] = idx

    return nodes, label_to_idx


def build_original_prev_next(nodes: List[ResidueNode]) -> Tuple[List[Optional[int]], List[Optional[int]]]:
    """Compute original prev/next along each chain using orig_res_seq sort."""
    by_chain: Dict[str, List[Tuple[int, int]]] = {}  # chain -> [(resSeq, idx), ...]
    for idx, n in enumerate(nodes):
        by_chain.setdefault(n.orig_chain_id, []).append((n.orig_res_seq, idx))
    orig_prev: List[Optional[int]] = [None] * len(nodes)
    orig_next: List[Optional[int]] = [None] * len(nodes)
    for chain_id, lst in by_chain.items():
        lst_sorted = sorted(lst, key=lambda t: t[0])
        for i, (_res, idx) in enumerate(lst_sorted):
            if i > 0:
                orig_prev[idx] = lst_sorted[i - 1][1]
            if i < len(lst_sorted) - 1:
                orig_next[idx] = lst_sorted[i + 1][1]
    return orig_prev, orig_next


# -------------------------- Bowtie phosphate handling --------------------------


_PHOS_ATOM_NAMES = {
    "P",
    "OP1",
    "OP2",
    "O1P",
    "O2P",
    "OP3",
    "O3P",
}


def cut_and_store_bowtie_phosphates(
    nodes: List[ResidueNode],
    label_to_idx: Dict[Label, int],
    bowtie_specs: List[Dict[str, object]],
) -> Dict[Label, List[edit_pdb_atom.pdb_atom_record]]:
    """For each bowtie spec, cut the phosphate group atoms from pos2 and store them.

    Returns:
        phos_store: dict mapping original pos2 label -> list of phosphate atom records

    Also mutates nodes[pos2].atoms and sets nodes[pos2].no_phosphate = True.
    """
    # Build the list of pos2 labels (keys) to be stored.
    phos_store: Dict[Label, List[edit_pdb_atom.pdb_atom_record]] = {}
    for sp in bowtie_specs:
        pos2 = sp["pos2"]  # type: ignore[index]
        assert isinstance(pos2, tuple)
        if pos2 in phos_store:
            raise ValueError(f"Duplicate bowtie pos2 (phosphate donor) residue: {pos2[1]}{pos2[0]}")
        phos_store[pos2] = []

    # Print essential info: keys for all members.
    if phos_store:
        keys_str = ", ".join([f"{r}{c}" for (c, r) in sorted(phos_store.keys(), key=lambda x: (x[0], x[1]))])
        print(f"Bowtie phosphate donor residues (pos2): {keys_str}")
    else:
        print("Bowtie phosphate donor residues (pos2): (none)")

    # Now actually cut atoms.
    for pos2 in list(phos_store.keys()):
        if pos2 not in label_to_idx:
            raise ValueError(f"Bowtie pos2 residue {pos2[1]}{pos2[0]} not found in PDB")
        idx2 = label_to_idx[pos2]
        node2 = nodes[idx2]
        if node2.no_phosphate:
            raise ValueError(f"Residue {pos2[1]}{pos2[0]} already had phosphate removed")

        phos_atoms: List[edit_pdb_atom.pdb_atom_record] = []
        keep_atoms: List[edit_pdb_atom.pdb_atom_record] = []
        for a in node2.atoms:
            an = a.name.strip()
            if an in _PHOS_ATOM_NAMES:
                phos_atoms.append(a)
            else:
                keep_atoms.append(a)

        # Basic sanity: require at least P plus two non-bridging oxygens.
        names = {a.name.strip() for a in phos_atoms}
        if "P" not in names:
            raise ValueError(
                f"Cannot form bowtie phosphate from {pos2[1]}{pos2[0]}: missing atom 'P'"
            )
        # Count non-bridging O among common names.
        non_bridge = [n for n in names if n in {"OP1", "OP2", "O1P", "O2P", "OP3", "O3P"}]
        if len(non_bridge) < 2:
            raise ValueError(
                f"Cannot form bowtie phosphate from {pos2[1]}{pos2[0]}: "
                f"expected >=2 non-bridging O atoms, found {sorted(non_bridge)}"
            )

        node2.atoms = keep_atoms
        node2.no_phosphate = True
        phos_store[pos2] = phos_atoms

        print(
            f"  Cut phosphate group from {pos2[1]}{pos2[0]}: moved {len(phos_atoms)} atoms "
            f"({', '.join(sorted(names))}) into storage"
        )

    return phos_store


def _set_atom_name(atom: edit_pdb_atom.pdb_atom_record, new_name: str) -> None:
    """Update the atom-name field in a pdb_atom_record in-place."""
    atom.name = new_name
    atom.string = atom.string[:12] + f"{new_name:>4s}" + atom.string[16:]


def _canonicalize_linker_phosphate_atoms(
    phos_atoms: List[edit_pdb_atom.pdb_atom_record],
    source_label: Label,
    linker_phosphate_style: Optional[object] = None,
) -> List[edit_pdb_atom.pdb_atom_record]:
    """Return exactly P, OP1, OP2 atoms for a standalone linker phosphate.

    The phosphate donor can use either OP1/OP2 or O1P/O2P naming.  We keep
    the coordinates of the donor phosphate atoms, rename the non-bridging
    oxygens to OP1/OP2, and mark all three records with the configured record
    type and residue name.  Extra phosphate atoms, if any, are intentionally
    not carried into the linker because the inserted residue is defined as a
    three-atom phosphate-only group.
    """
    style = _coerce_linker_phosphate_style(linker_phosphate_style)
    by_name: Dict[str, edit_pdb_atom.pdb_atom_record] = {}
    for atom in phos_atoms:
        by_name.setdefault(atom.name.strip().upper(), atom)

    p_atom = by_name.get("P")
    if p_atom is None:
        raise ValueError(
            f"3'-3' linker phosphate from {source_label[1]}{source_label[0]} is missing atom P"
        )

    preferred_op1 = ["OP1", "O1P", "OP3", "O3P", "OP2", "O2P"]
    preferred_op2 = ["OP2", "O2P", "OP3", "O3P", "OP1", "O1P"]

    op1_atom: Optional[edit_pdb_atom.pdb_atom_record] = None
    for nm in preferred_op1:
        cand = by_name.get(nm)
        if cand is not None and cand is not p_atom:
            op1_atom = cand
            break

    op2_atom: Optional[edit_pdb_atom.pdb_atom_record] = None
    for nm in preferred_op2:
        cand = by_name.get(nm)
        if cand is not None and cand is not p_atom and cand is not op1_atom:
            op2_atom = cand
            break

    if op1_atom is None or op2_atom is None:
        oxygen_atoms = [
            atom for atom in phos_atoms
            if atom is not p_atom and atom.name.strip().upper() in _PHOS_ATOM_NAMES
        ]
        if op1_atom is None and oxygen_atoms:
            op1_atom = oxygen_atoms.pop(0)
        if op2_atom is None:
            for atom in oxygen_atoms:
                if atom is not op1_atom:
                    op2_atom = atom
                    break

    if op1_atom is None or op2_atom is None:
        raise ValueError(
            f"3'-3' linker phosphate from {source_label[1]}{source_label[0]} "
            "needs two non-bridging O atoms"
        )

    canonical = [(p_atom, "P"), (op1_atom, "OP1"), (op2_atom, "OP2")]
    for atom, atom_name in canonical:
        atom.update_recordName(style.record_name)
        atom.update_resName(style.resname)
        _set_atom_name(atom, atom_name)

    return [p_atom, op1_atom, op2_atom]


def _canonicalize_x33_phosphate_atoms(
    phos_atoms: List[edit_pdb_atom.pdb_atom_record],
    source_label: Label,
) -> List[edit_pdb_atom.pdb_atom_record]:
    """Compatibility wrapper for callers that expect the historical X33 style."""
    return _canonicalize_linker_phosphate_atoms(
        phos_atoms,
        source_label,
        make_linker_phosphate_style(X33_HETID, LINKER_PHOSPHATE_DEFAULT_RECORD),
    )


# -------------------------- Graph construction --------------------------


def build_original_graph(n_nodes: int, orig_next: List[Optional[int]]) -> BackboneGraph:
    g = BackboneGraph(n_nodes)
    for i in range(n_nodes):
        j = orig_next[i]
        if j is None:
            continue
        # Standard 3'-5' bond: O3'(i) -- P(j), so store i as O3 end and j as P end.
        g.add_edge(i, j, kind="std", end_a="O3", end_b="P")
    return g


def _normalize_backbone_atom_name(atom_name: str) -> str:
    """Normalize PDB atom aliases used by re_helix backbone LINK records."""
    name = atom_name.strip().upper().replace("*", "'")
    if name in {"O3", "O3'"}:
        return "O3"
    if name in {"O5", "O5'"}:
        return "O5"
    return name


def _input_link_label(chain_id: str, res_seq: int) -> Label:
    return (chain_id.strip() or " ", int(res_seq))


def _classify_backbone_link(
    link: edit_pdb_link.pdb_link_record,
) -> Optional[Tuple[str, str, str]]:
    """Return ``(kind, endpoint1, endpoint2)`` for supported backbone LINKs.

    re_helix writes three backbone-link forms: P--O3' standard/inverted
    phosphodiester links, P--O5' bowtie links, and two P--O3' links incident
    on a phosphate-only bridge residue. Other LINK chemistry is retained as a
    passthrough record but does not participate in reciprocal-exchange routing.
    """
    if (link.sym1.strip() or "1555") != "1555" or (link.sym2.strip() or "1555") != "1555":
        return None

    end1 = _normalize_backbone_atom_name(link.name1)
    end2 = _normalize_backbone_atom_name(link.name2)
    endpoint_set = {end1, end2}
    if endpoint_set == {"P", "O3"}:
        return ("std", end1, end2)
    if endpoint_set == {"P", "O5"}:
        return ("5to5", end1, end2)
    return None


def _is_phosphate_only_node(node: ResidueNode) -> bool:
    names = {_normalize_backbone_atom_name(atom.name) for atom in node.atoms}
    non_bridging = names.intersection({"OP1", "OP2", "O1P", "O2P", "OP3", "O3P"})
    return "P" in names and len(non_bridging) >= 2 and names.issubset(_PHOS_ATOM_NAMES)


def _build_graph_predecessors(
    g: BackboneGraph,
    nodes: List[ResidueNode],
) -> List[Optional[int]]:
    """Return each residue's chemical predecessor through its P endpoint."""
    orig_prev: List[Optional[int]] = [None] * len(nodes)
    for idx, node in enumerate(nodes):
        candidates: List[int] = []
        for neighbor in g.neigh[idx]:
            edge = g.get_edge(idx, neighbor)
            end_here, end_there = edge.endpoints(idx, neighbor)
            if edge.kind == "std" and end_here == "P" and end_there == "O3":
                candidates.append(neighbor)

        # A standalone linker phosphate intentionally has two P--O3' bonds and
        # therefore no unique incoming nucleotide edge. It cannot itself be an
        # exchange position.
        if node.is_phos_bridge:
            continue
        if len(candidates) > 1:
            labels = ", ".join(
                f"{nodes[n].orig_res_seq}{nodes[n].orig_chain_id}" for n in candidates
            )
            raise ValueError(
                f"Residue {node.orig_res_seq}{node.orig_chain_id} has multiple P-side "
                f"backbone predecessors ({labels}); input LINK topology is branched."
            )
        if candidates:
            orig_prev[idx] = candidates[0]
    return orig_prev


def build_input_backbone_topology(
    nodes: List[ResidueNode],
    label_to_idx: Dict[Label, int],
    input_links: Optional[List[edit_pdb_link.pdb_link_record]] = None,
) -> InputBackboneTopology:
    """Overlay supported input LINK records on the implicit chain graph.

    Coordinate order supplies the default O3'(i)--P(i+1) edges. A supported
    LINK between the same residue pair replaces that default edge with its
    explicit endpoint orientation; a non-adjacent supported LINK adds a cycle
    closure. Phosphate-only residues with two P--O3' links are marked as bridge
    nodes so their two links are regenerated after output renumbering.
    """
    links = list(input_links or [])
    _numeric_prev, numeric_next = build_original_prev_next(nodes)
    graph = build_original_graph(len(nodes), numeric_next)

    descriptors: List[Tuple[int, int, str, str, str, edit_pdb_link.pdb_link_record]] = []
    passthrough: List[edit_pdb_link.pdb_link_record] = []
    seen_pairs: Set[frozenset[int]] = set()

    for link in links:
        classified = _classify_backbone_link(link)
        if classified is None:
            passthrough.append(link)
            continue

        label1 = _input_link_label(link.chainID1, link.resSeq1)
        label2 = _input_link_label(link.chainID2, link.resSeq2)
        if label1 not in label_to_idx or label2 not in label_to_idx:
            missing = label1 if label1 not in label_to_idx else label2
            raise ValueError(
                "Backbone LINK endpoint references a residue absent from the PDB: "
                f"{missing[1]}{missing[0]} in {link.string.strip()}"
            )
        idx1 = label_to_idx[label1]
        idx2 = label_to_idx[label2]
        if idx1 == idx2:
            passthrough.append(link)
            continue

        pair_key = frozenset((idx1, idx2))
        if pair_key in seen_pairs:
            raise ValueError(
                "Multiple supported backbone LINK records connect the same residue pair: "
                f"{link.string.strip()}"
            )
        seen_pairs.add(pair_key)
        kind, end1, end2 = classified
        descriptors.append((idx1, idx2, kind, end1, end2, link))

    # Identify existing standalone X33/custom/ATOM-DA phosphate bridges from
    # structure rather than residue name, so every output style is composable.
    p_to_o3_neighbors: Dict[int, Set[int]] = {}
    for idx1, idx2, kind, end1, end2, _link in descriptors:
        if kind != "std":
            continue
        if end1 == "P" and end2 == "O3":
            p_to_o3_neighbors.setdefault(idx1, set()).add(idx2)
        elif end2 == "P" and end1 == "O3":
            p_to_o3_neighbors.setdefault(idx2, set()).add(idx1)

    phosphate_bridge_count = 0
    for idx, neighbors in p_to_o3_neighbors.items():
        if len(neighbors) >= 2 and _is_phosphate_only_node(nodes[idx]):
            nodes[idx].is_phos_bridge = True
            nodes[idx].phos_source = nodes[idx].orig_label()
            phosphate_bridge_count += 1

    # Remove implicit edges replaced by explicit LINK chemistry first, then add
    # all LINK edges. Doing the removals as a separate pass avoids transient
    # degree-three failures when a non-adjacent closure is also present.
    for idx1, idx2, _kind, _end1, _end2, _link in descriptors:
        pair_key = frozenset((idx1, idx2))
        if pair_key in graph.edges:
            graph.remove_edge(idx1, idx2)

    for idx1, idx2, kind, end1, end2, link in descriptors:
        try:
            graph.add_edge(idx1, idx2, kind=kind, end_a=end1, end_b=end2)
        except ValueError as exc:
            raise ValueError(
                f"Invalid backbone topology while applying LINK '{link.string.strip()}': {exc}"
            ) from exc

    return InputBackboneTopology(
        graph=graph,
        orig_prev=_build_graph_predecessors(graph, nodes),
        passthrough_links=passthrough,
        recognized_link_count=len(descriptors),
        phosphate_bridge_count=phosphate_bridge_count,
    )


def apply_exchanges_to_graph(
    g: BackboneGraph,
    nodes: List[ResidueNode],
    label_to_idx: Dict[Label, int],
    orig_prev: List[Optional[int]],
    specs: List[Dict[str, object]],
) -> Tuple[Set[int], int, int, int]:
    """Apply all exchanges as edge rewires on the original graph.

    Returns:
        junction_nodes: set of node indices involved in exchanges (pos1,pos2,prev1,prev2)
        n_double, n_single, n_bowtie
    """
    cut_edges: Set[frozenset[int]] = set()
    add_edges: List[Tuple[int, int, str, str, str, Optional[Label]]] = []
    junction_nodes: Set[int] = set()

    n_double = n_single = n_bowtie = 0

    for sp in specs:
        pos1 = sp["pos1"]  # type: ignore[index]
        pos2 = sp["pos2"]  # type: ignore[index]
        kind = sp["kind"]  # type: ignore[index]
        assert isinstance(pos1, tuple) and isinstance(pos2, tuple)
        assert isinstance(kind, str)

        if pos1 not in label_to_idx:
            raise ValueError(f"Residue {pos1[1]}{pos1[0]} not found in PDB")
        if pos2 not in label_to_idx:
            raise ValueError(f"Residue {pos2[1]}{pos2[0]} not found in PDB")
        idx1 = label_to_idx[pos1]
        idx2 = label_to_idx[pos2]

        u1 = orig_prev[idx1]
        u2 = orig_prev[idx2]
        if u1 is None:
            raise ValueError(
                f"Residue {pos1[1]}{pos1[0]} has no unique P-side backbone predecessor; "
                "cannot cut its incoming edge"
            )
        if u2 is None:
            raise ValueError(
                f"Residue {pos2[1]}{pos2[0]} has no unique P-side backbone predecessor; "
                "cannot cut its incoming edge"
            )

        # Bookkeeping for circular permutation avoidance.
        junction_nodes.update({idx1, idx2, u1, u2})

        # Record cuts (original incoming edges).
        for a, b in ((u1, idx1), (u2, idx2)):
            k = frozenset((a, b))
            if k in cut_edges:
                raise ValueError(
                    "Overlapping exchanges: the same backbone edge is being cut more than once: "
                    f"{nodes[a].orig_res_seq}{nodes[a].orig_chain_id}-{nodes[b].orig_res_seq}{nodes[b].orig_chain_id}"
                )
            cut_edges.add(k)

        if kind == "double":
            n_double += 1
            # Swap incoming edges: u1->idx2 and u2->idx1 (std edges)
            add_edges.append((u1, idx2, "std", "O3", "P", None))
            add_edges.append((u2, idx1, "std", "O3", "P", None))
        elif kind == "single":
            n_single += 1
            # Single: only connect u1->idx2
            add_edges.append((u1, idx2, "std", "O3", "P", None))
        elif kind == "bowtie":
            n_bowtie += 1
            # Bowtie special edges:
            #   - 3'-3' between u1 and u2 (phosphate donor is pos2)
            #   - 5'-5' between idx1 (P) and idx2 (O5)
            add_edges.append((u1, u2, "3to3", "O3", "O3", pos2))
            add_edges.append((idx1, idx2, "5to5", "P", "O5", None))
        else:
            raise ValueError(f"Unsupported exchange kind: {kind}")

    # Apply all cuts (from original graph)
    for k in cut_edges:
        a, b = tuple(k)
        g.remove_edge(a, b)

    # Apply all additions
    for a, b, kind, end_a, end_b, phos_key in add_edges:
        g.add_edge(a, b, kind=kind, end_a=end_a, end_b=end_b, phos_key=phos_key)

    return junction_nodes, n_double, n_single, n_bowtie


# -------------------------- Component traversal & ordering --------------------------


def _collect_component_nodes(g: BackboneGraph, start: int, visited: Set[int]) -> List[int]:
    """Return list of nodes in the connected component containing start."""
    stack = [start]
    comp: List[int] = []
    visited.add(start)
    while stack:
        u = stack.pop()
        comp.append(u)
        for v in g.neigh[u]:
            if v not in visited:
                visited.add(v)
                stack.append(v)
    return comp


def _traverse_path(g: BackboneGraph, start: int) -> List[int]:
    """Traverse a path component from one end to the other."""
    order: List[int] = []
    prev: Optional[int] = None
    cur = start
    while True:
        order.append(cur)
        nxts = [v for v in g.neigh[cur] if v != prev]
        if not nxts:
            break
        # Degree<=2, so at most one next if we avoid prev.
        nxt = nxts[0]
        prev, cur = cur, nxt
    return order


def _traverse_cycle(g: BackboneGraph, start: int, first_step: int) -> List[int]:
    """Traverse a cycle, returning a list with each node exactly once."""
    order = [start]
    prev = start
    cur = first_step
    while cur != start:
        order.append(cur)
        nxts = [v for v in g.neigh[cur] if v != prev]
        if not nxts:
            raise ValueError("Broken cycle traversal (dead end)")
        nxt = nxts[0]
        prev, cur = cur, nxt
    return order


def _iter_component_pairs(
    order: List[int],
    is_cycle: bool = False,
):
    """Yield ``(left, right, is_closing_pair)`` for a component order."""
    for i in range(len(order) - 1):
        yield order[i], order[i + 1], False
    if is_cycle and len(order) > 1:
        yield order[-1], order[0], True


def _inverted_cost(order: List[int], g: BackboneGraph, is_cycle: bool = False) -> int:
    """Count standard-edge LINKs required by a serialized component order.

    Inverted P->O3 steps require LINK records everywhere. A cycle's closing
    edge also requires a LINK even in the natural O3->P direction because it
    crosses the output TER boundary.
    """
    cost = 0
    for a, b, is_closing_pair in _iter_component_pairs(order, is_cycle):
        e = g.get_edge(a, b)
        if e.kind != "std":
            continue
        end_a, end_b = e.endpoints(a, b)
        if (end_a == "P" and end_b == "O3") or is_closing_pair:
            cost += 1
    return cost


def _input_forward_continuity_score(
    order: List[int],
    nodes: List[ResidueNode],
    is_cycle: bool = False,
) -> int:
    """Count adjacencies that retain the input PDB's directed chain order.

    LINK topology describes chemical connectivity, but it must not cause a
    later RE run to serialize an entire strand backwards merely to save a LINK
    record.  For cycles the closing pair is included so this score is invariant
    under circular rotation; an input chain's last residue does not wrap to its
    first because their ranks are not consecutive.
    """
    score = 0
    for a, b, _is_closing_pair in _iter_component_pairs(order, is_cycle):
        node_a = nodes[a]
        node_b = nodes[b]
        if (
            node_a.orig_chain_id == node_b.orig_chain_id
            and node_b.input_chain_rank == node_a.input_chain_rank + 1
        ):
            score += 1
    return score


def _orientation_tiebreak_key(
    order: List[int],
    nodes: List[ResidueNode],
) -> Tuple[Tuple[str, int, int, int], ...]:
    """Return a deterministic provenance key for an already positioned order."""
    return tuple(
        (
            nodes[idx].orig_chain_id,
            nodes[idx].input_chain_rank,
            nodes[idx].orig_res_seq,
            idx,
        )
        for idx in order
    )


def _end_fragment_dir(
    order: List[int],
    nodes: List[ResidueNode],
    from_start: bool,
) -> int:
    """Infer the input direction of the fragment exposed at one path end.

    Only the pair immediately incident on the requested endpoint can establish
    that terminal fragment's direction.  A chain switch, a singleton fragment,
    or a non-consecutive provenance jump has no inferable input direction.
    """
    if len(order) < 2:
        return 0
    pair_index = 0 if from_start else len(order) - 2
    node_a = nodes[order[pair_index]]
    node_b = nodes[order[pair_index + 1]]
    if node_a.orig_chain_id != node_b.orig_chain_id:
        return 0
    delta = node_b.input_chain_rank - node_a.input_chain_rank
    if delta == 1:
        return 1
    if delta == -1:
        return -1
    return 0


def _both_terminal_fragments_backward(
    order: List[int],
    nodes: List[ResidueNode],
) -> bool:
    return (
        _end_fragment_dir(order, nodes, from_start=True) == -1
        and _end_fragment_dir(order, nodes, from_start=False) == -1
    )


def _terminal_direction_score(
    order: List[int],
    nodes: List[ResidueNode],
) -> int:
    """Score the user-visible directions of a path's two terminal fragments.

    Each terminal fragment contributes +1 when it follows input serialization,
    -1 when it runs backwards, and 0 when no direction can be inferred.  A path
    with both exposed ends running forward therefore outranks its globally
    reversed representation even when the latter would preserve more interior
    adjacencies or require fewer LINK records.
    """
    return _end_fragment_dir(order, nodes, from_start=True) + _end_fragment_dir(
        order,
        nodes,
        from_start=False,
    )


def _cycle_inverted_direction_cost(order: List[int], g: BackboneGraph) -> int:
    """Return a rotation-independent inverted-edge cost for a cycle direction."""
    cost = 0
    for a, b, _is_closing_pair in _iter_component_pairs(order, is_cycle=True):
        edge = g.get_edge(a, b)
        if edge.kind != "std":
            continue
        end_a, end_b = edge.endpoints(a, b)
        if end_a == "P" and end_b == "O3":
            cost += 1
    return cost


def _output_link_cost(
    order: List[int],
    g: BackboneGraph,
    nodes: List[ResidueNode],
    is_cycle: bool = False,
) -> int:
    """Count topology LINK records emitted for an already positioned order.

    Unlike ``_inverted_cost``, this models the writer's actual behavior.  It
    therefore counts explicit cycle closures and special backbone chemistry,
    and it omits the exact serialized break of an open cycle.
    """
    cost = 0
    for a, b, is_closing_pair in _iter_component_pairs(order, is_cycle):
        node_a = nodes[a]
        node_b = nodes[b]
        if node_a.is_phos_bridge or node_b.is_phos_bridge:
            cost += 1
            continue

        edge = g.get_edge(a, b)
        if edge.kind == "5to5":
            cost += 1
            continue
        if edge.kind == "3to3":
            # insert_phosphate_nodes expands this edge into two P--O3' LINKs.
            cost += 2
            continue
        if edge.kind != "std":
            raise ValueError(f"Unsupported edge kind while scoring LINK records: {edge.kind}")

        end_a, end_b = edge.endpoints(a, b)
        inverted = end_a == "P" and end_b == "O3"
        if inverted or is_closing_pair:
            cost += 1
    return cost


def orient_path_component(
    order: List[int],
    nodes: List[ResidueNode],
    g: BackboneGraph,
    min_link_records: bool = False,
) -> Tuple[List[int], bool]:
    """Choose orientation for a path component.

    Returns:
        (oriented_order, was_reversed)

    By default, the terminal input fragments determine the user-visible strand
    direction, followed by total input-serialization continuity and a stable
    provenance key.  When ``min_link_records`` is true, the topology LINK count
    emitted by the writer becomes primary and may reverse the entire strand.
    """
    if not order:
        return order, False
    rev = list(reversed(order))
    candidates = ((order, False), (rev, True))
    if min_link_records:
        chosen, was_reversed = min(
            candidates,
            key=lambda item: (
                _output_link_cost(item[0], g, nodes, is_cycle=False),
                -_terminal_direction_score(item[0], nodes),
                -_input_forward_continuity_score(item[0], nodes),
                _orientation_tiebreak_key(item[0], nodes),
            ),
        )
    else:
        chosen, was_reversed = min(
            candidates,
            key=lambda item: (
                -_terminal_direction_score(item[0], nodes),
                -_input_forward_continuity_score(item[0], nodes),
                _orientation_tiebreak_key(item[0], nodes),
            ),
        )
    return chosen, was_reversed


def choose_cycle_orientation(
    order: List[int],
    g: BackboneGraph,
    circularize: bool = False,
) -> List[int]:
    """Compatibility helper: choose the lower-LINK direction for a cycle.

    New component construction uses input-provenance-aware selection in
    ``_choose_and_rotate_cycle_orientation`` below. This wrapper retains the
    historical public call shape and return type for external callers.
    """
    if not order:
        return order
    rev = list(reversed(order))
    return min(
        (order, rev),
        key=lambda candidate: _inverted_cost(candidate, g, is_cycle=circularize),
    )


def _choose_and_rotate_cycle_orientation(
    order: List[int],
    nodes: List[ResidueNode],
    g: BackboneGraph,
    junction_nodes: Set[int],
    cir_shift: int,
    min_link_records: bool = False,
    circularize_cycles: bool = False,
) -> Tuple[List[int], int, bool]:
    """Choose and rotate a cycle while preserving input serialization.

    Input provenance is primary by default, and direction is then rotated by
    the exact requested shift.  With ``min_link_records`` enabled, each
    direction is first rotated by that exact shift and scored according to the
    LINK records the writer will emit for an open or circularized cycle.
    """
    if not order:
        return order, 0, False
    # Keep both directions anchored at the same node so the final deterministic
    # key cannot depend on rotation or neighbor insertion order.
    rev = [order[0]] + list(reversed(order[1:]))
    candidates = ((order, False), (rev, True))
    if min_link_records:
        positioned_candidates = []
        for candidate, was_reversed in candidates:
            rotated, rotation_index = rotate_cycle_away_from_junctions(
                candidate,
                g,
                junction_nodes,
                cir_shift,
            )
            positioned_candidates.append(
                (candidate, rotated, rotation_index, was_reversed)
            )
        chosen, rotated, rotation_index, was_reversed = min(
            positioned_candidates,
            key=lambda item: (
                _output_link_cost(
                    item[1],
                    g,
                    nodes,
                    is_cycle=circularize_cycles,
                ),
                -_input_forward_continuity_score(item[0], nodes, is_cycle=True),
                _orientation_tiebreak_key(item[0], nodes),
            ),
        )
        return rotated, rotation_index, was_reversed
    else:
        chosen, was_reversed = min(
            candidates,
            key=lambda item: (
                -_input_forward_continuity_score(item[0], nodes, is_cycle=True),
                _orientation_tiebreak_key(item[0], nodes),
            ),
        )
    rotated, rotation_index = rotate_cycle_away_from_junctions(
        chosen,
        g,
        junction_nodes,
        cir_shift,
    )
    return rotated, rotation_index, was_reversed


def rotate_cycle_away_from_junctions(
    order: List[int],
    g: BackboneGraph,
    junction_nodes: Set[int],
    cir_shift: int,
) -> Tuple[List[int], int]:
    """Apply the requested exact circular shift to a cycle order.

    The function retains its historical name and arguments for compatibility,
    but no longer advances beyond the requested index to avoid junctions or
    special edges. Consequently, ``cir_shift=0`` performs no permutation and
    signed shifts are applied exactly modulo the cycle length.

    Returns:
        (rotated_order, start_index_in_original_order)
    """
    n = len(order)
    if n == 0:
        return order, 0

    rotation_index = cir_shift % n
    rotated = order[rotation_index:] + order[:rotation_index]
    return rotated, rotation_index


def build_ordered_components(
    g: BackboneGraph,
    nodes: List[ResidueNode],
    junction_nodes: Set[int],
    cir_shift: int,
    circularize_cycles: bool = False,
    min_link_records: bool = False,
) -> List[Dict[str, object]]:
    """Return ordered connected components (paths/cycles) of the backbone graph."""
    visited: Set[int] = set()
    components: List[Dict[str, object]] = []

    for start in range(g.n_nodes):
        if start in visited:
            continue
        comp_nodes = _collect_component_nodes(g, start, visited)
        # Determine if cycle.
        is_cycle = all(len(g.neigh[n]) == 2 for n in comp_nodes)

        if not is_cycle:
            # Choose an end as traversal start (degree<=1). If none, pick min node.
            ends = [n for n in comp_nodes if len(g.neigh[n]) <= 1]
            start_end = min(ends) if ends else min(comp_nodes)
            order = _traverse_path(g, start_end)
            order, was_rev = orient_path_component(
                order,
                nodes,
                g,
                min_link_records=min_link_records,
            )
            components.append({
                "order": order,
                "is_cycle": False,
                "circularized": False,
                "was_reversed": was_rev,
                "rotation": None,
            })
        else:
            # Cycle: anchor traversal to the input serialization direction when
            # that successor edge is still present. Fall back to a deterministic
            # provenance key rather than LINK insertion order.
            start_cycle = min(
                comp_nodes,
                key=lambda idx: (
                    nodes[idx].orig_chain_id,
                    nodes[idx].input_chain_rank,
                    nodes[idx].orig_res_seq,
                    idx,
                ),
            )
            nbs = g.neigh[start_cycle]
            if len(nbs) != 2:
                raise ValueError("Cycle node does not have 2 neighbors")
            start_node = nodes[start_cycle]
            forward_nbs = [
                neighbor
                for neighbor in nbs
                if nodes[neighbor].orig_chain_id == start_node.orig_chain_id
                and nodes[neighbor].input_chain_rank == start_node.input_chain_rank + 1
            ]
            if forward_nbs:
                first_step = min(forward_nbs)
            else:
                first_step = min(
                    nbs,
                    key=lambda idx: (
                        nodes[idx].orig_chain_id,
                        nodes[idx].input_chain_rank,
                        nodes[idx].orig_res_seq,
                        idx,
                    ),
                )
            order = _traverse_cycle(g, start_cycle, first_step)
            # Select direction by input provenance, then rotate consistently
            # with the requested cir_shift and open/closed cycle mode.
            order, rot_idx, was_rev = _choose_and_rotate_cycle_orientation(
                order,
                nodes,
                g,
                junction_nodes,
                cir_shift,
                min_link_records=min_link_records,
                circularize_cycles=circularize_cycles,
            )
            components.append({
                "order": order,
                "is_cycle": True,
                "circularized": circularize_cycles,
                "was_reversed": was_rev,
                "rotation": rot_idx,
            })

    return components


# -------------------------- Phosphate insertion into orders --------------------------


def insert_phosphate_nodes(
    base_order: List[int],
    g: BackboneGraph,
    nodes: List[ResidueNode],
    phos_store: Dict[Label, List[edit_pdb_atom.pdb_atom_record]],
    used_phos: Set[Label],
    linker_phosphate_style: Optional[object] = None,
    is_cycle: bool = False,
) -> List[int]:
    """Expand 3to3 edges by inserting phosphate-only residue nodes."""
    if not base_order:
        return base_order

    style = _coerce_linker_phosphate_style(linker_phosphate_style)
    out: List[int] = []

    def append_bridge_for_edge(a: int, b: int) -> None:
        e = g.get_edge(a, b)
        if e.kind == "3to3":
            if e.phos_key is None:
                raise ValueError("3to3 edge missing phosphate key")
            phos_key = e.phos_key
            if phos_key in used_phos:
                raise ValueError(
                    f"Phosphate donor {phos_key[1]}{phos_key[0]} used more than once"
                )
            if phos_key not in phos_store:
                raise ValueError(
                    f"Phosphate donor {phos_key[1]}{phos_key[0]} not found in store (bowtie parse mismatch?)"
                )
            phos_atoms = phos_store[phos_key]
            if not phos_atoms:
                raise ValueError(
                    f"Phosphate donor {phos_key[1]}{phos_key[0]} has empty atom list (was it cut?)"
                )

            # Create a new standalone phosphate-only linker residue.
            linker_atoms = _canonicalize_linker_phosphate_atoms(phos_atoms, phos_key, style)
            new_idx = len(nodes)
            nodes.append(
                ResidueNode(
                    orig_chain_id=phos_key[0],
                    orig_res_seq=phos_key[1],
                    atoms=linker_atoms,
                    is_phos_bridge=True,
                    phos_source=phos_key,
                    no_phosphate=False,
                )
            )
            used_phos.add(phos_key)
            out.append(new_idx)

    for i in range(len(base_order) - 1):
        a = base_order[i]
        b = base_order[i + 1]
        out.append(a)
        append_bridge_for_edge(a, b)

    out.append(base_order[-1])
    if is_cycle and len(base_order) > 1:
        append_bridge_for_edge(base_order[-1], base_order[0])
    return out


# -------------------------- Debug path printing --------------------------


def _node_label(nodes: List[ResidueNode], idx: int) -> str:
    n = nodes[idx]
    base = f"{n.orig_res_seq}{n.orig_chain_id}"
    if n.is_phos_bridge:
        return f"Phos({base})"
    if n.no_phosphate:
        return f"{base}[-P]"
    return base


def _arrow_for_pair(a: int, b: int, g: BackboneGraph, nodes: List[ResidueNode]) -> str:
    na, nb = nodes[a], nodes[b]
    if na.is_phos_bridge or nb.is_phos_bridge:
        return "->>"
    e = g.get_edge(a, b)
    if e.kind != "std":
        return "->>"
    end_a, end_b = e.endpoints(a, b)
    if end_a == "O3" and end_b == "P":
        return "->"
    return "->>"


def format_path_string(order: List[int], g: BackboneGraph, nodes: List[ResidueNode]) -> str:
    if not order:
        return "(empty)"
    parts: List[str] = [_node_label(nodes, order[0])]
    for i in range(len(order) - 1):
        a, b = order[i], order[i + 1]
        parts.append(_arrow_for_pair(a, b, g, nodes))
        parts.append(_node_label(nodes, b))
    return " ".join(parts)


# -------------------------- LINK record generation --------------------------


def _distance(a: edit_pdb_atom.pdb_atom_record, b: edit_pdb_atom.pdb_atom_record) -> float:
    return math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2)


def _format_link_line(a: edit_pdb_atom.pdb_atom_record, b: edit_pdb_atom.pdb_atom_record) -> str:
    """Create a canonical PDB LINK line connecting atom *a* and atom *b*.

    The correct LINK layout places the two symmetry operator fields (sym1/sym2)
    before the optional distance:

        LINK ...  sym1  sym2  dist

    We use the conventional dummy symmetry operators '1555'/'1555' and write
    the observed distance in Å in the final field.
    """
    dist = _distance(a, b)
    return (
        "LINK        "
        f"{a.name:>4s} {a.resName:>3s} {a.chainID:1s}{a.resSeq:4d}"
        "                "
        f"{b.name:>4s} {b.resName:>3s} {b.chainID:1s}{b.resSeq:4d}"
        f"     1555   1555 {dist:5.2f}\n"
    )


def _find_required_atom(
    atom_index: Dict[Tuple[str, int, str], edit_pdb_atom.pdb_atom_record],
    chain_id: str,
    res_seq: int,
    names: List[str],
) -> Optional[edit_pdb_atom.pdb_atom_record]:
    for nm in names:
        a = edit_pdb_link.find_atom(atom_index, chain_id, res_seq, nm)
        if a is not None:
            return a
    return None


def build_link_records(
    component_orders: List[Dict[str, object]],
    g: BackboneGraph,
    nodes: List[ResidueNode],
    output_atoms: List[edit_pdb_atom.pdb_atom_record],
) -> Tuple[List[edit_pdb_link.pdb_link_record], Dict[str, int]]:
    """Generate LINK records based on final ordered strands."""
    atom_index = edit_pdb_link.build_atom_index(output_atoms)

    links: List[edit_pdb_link.pdb_link_record] = []
    counts = {
        "total": 0,
        "backbone_inverted": 0,
        "bowtie_5to5": 0,
        "bowtie_3to3": 0,
        "circular_closure": 0,
        "other_preserved": 0,
    }

    for comp in component_orders:
        order: List[int] = comp["order"]  # type: ignore[index]
        include_closure = bool(comp.get("is_cycle")) and bool(comp.get("circularized"))
        for a_idx, b_idx, is_closing_pair in _iter_component_pairs(order, include_closure):
            a_node, b_node = nodes[a_idx], nodes[b_idx]

            # Case: phosphate-only residue involved -> bowtie 3'-3'
            if a_node.is_phos_bridge or b_node.is_phos_bridge:
                phos_idx = a_idx if a_node.is_phos_bridge else b_idx
                res_idx = b_idx if a_node.is_phos_bridge else a_idx
                phos_node = nodes[phos_idx]
                res_node = nodes[res_idx]

                p_atom = _find_required_atom(atom_index, phos_node.new_chain_id, phos_node.new_res_seq, ["P"])
                o3_atom = _find_required_atom(
                    atom_index, res_node.new_chain_id, res_node.new_res_seq, ["O3'", "O3*"]
                )
                if p_atom is None or o3_atom is None:
                    print(
                        f"Warning: could not build 3'-3' LINK for adjacency "
                        f"{_node_label(nodes, a_idx)} - {_node_label(nodes, b_idx)} (missing P or O3')",
                        file=sys.stderr,
                    )
                    continue
                links.append(edit_pdb_link.pdb_link_record(_format_link_line(p_atom, o3_atom)))
                counts["bowtie_3to3"] += 1
                if is_closing_pair:
                    counts["circular_closure"] += 1
                continue

            # Both are standard residues.
            e = g.get_edge(a_idx, b_idx)

            if e.kind == "5to5":
                # Always LINK between P and O5'
                end_a, end_b = e.endpoints(a_idx, b_idx)
                if end_a == "P" and end_b == "O5":
                    p_node, o5_node = a_node, b_node
                elif end_a == "O5" and end_b == "P":
                    p_node, o5_node = b_node, a_node
                else:
                    raise ValueError("Malformed 5to5 edge endpoints")

                p_atom = _find_required_atom(atom_index, p_node.new_chain_id, p_node.new_res_seq, ["P"])
                o5_atom = _find_required_atom(atom_index, o5_node.new_chain_id, o5_node.new_res_seq, ["O5'", "O5*"])
                if p_atom is None or o5_atom is None:
                    print(
                        f"Warning: could not build 5'-5' LINK between "
                        f"{p_node.new_res_seq}{p_node.new_chain_id} and {o5_node.new_res_seq}{o5_node.new_chain_id} "
                        f"(missing P or O5')",
                        file=sys.stderr,
                    )
                    continue
                links.append(edit_pdb_link.pdb_link_record(_format_link_line(o5_atom, p_atom)))
                counts["bowtie_5to5"] += 1
                if is_closing_pair:
                    counts["circular_closure"] += 1
                continue

            if e.kind != "std":
                # 3to3 edges should have been expanded by inserting phosphate nodes.
                raise ValueError(
                    f"Unexpected non-std edge kind '{e.kind}' between residues in final order; "
                    f"missing phosphate insertion?"
                )

            # Standard edge: add LINK when traversal is inverted (P->O3). A
            # natural O3->P edge normally relies on implicit chain adjacency,
            # but a cycle's last->first edge crosses TER and must be explicit.
            end_a, end_b = e.endpoints(a_idx, b_idx)
            inverted = end_a == "P" and end_b == "O3"
            natural = end_a == "O3" and end_b == "P"
            if not inverted and not (is_closing_pair and natural):
                continue

            if inverted:
                p_node, o3_node = a_node, b_node
            else:
                p_node, o3_node = b_node, a_node
            p_atom = _find_required_atom(
                atom_index, p_node.new_chain_id, p_node.new_res_seq, ["P"]
            )
            o3_atom = _find_required_atom(
                atom_index, o3_node.new_chain_id, o3_node.new_res_seq, ["O3'", "O3*"]
            )
            if p_atom is None or o3_atom is None:
                link_type = "inverted-backbone" if inverted else "circular-closure"
                print(
                    f"Warning: could not build {link_type} LINK for adjacency "
                    f"{a_node.new_res_seq}{a_node.new_chain_id}->{b_node.new_res_seq}{b_node.new_chain_id} "
                    f"(missing P or O3')",
                    file=sys.stderr,
                )
                continue
            links.append(edit_pdb_link.pdb_link_record(_format_link_line(p_atom, o3_atom)))
            if is_closing_pair:
                counts["circular_closure"] += 1
            if inverted:
                counts["backbone_inverted"] += 1

    counts["total"] = len(links)
    return links, counts


def _atom_name_candidates(atom_name: str) -> List[str]:
    normalized = _normalize_backbone_atom_name(atom_name)
    if normalized == "O3":
        return [atom_name.strip(), "O3'", "O3*", "O3"]
    if normalized == "O5":
        return [atom_name.strip(), "O5'", "O5*", "O5"]
    return [atom_name.strip()]


def remap_passthrough_link_records(
    input_links: List[edit_pdb_link.pdb_link_record],
    nodes: List[ResidueNode],
    label_to_idx: Dict[Label, int],
    output_atoms: List[edit_pdb_atom.pdb_atom_record],
) -> List[edit_pdb_link.pdb_link_record]:
    """Remap non-backbone input LINKs through output residue provenance."""
    if not input_links:
        return []

    atom_index = edit_pdb_link.build_atom_index(output_atoms)
    remapped: List[edit_pdb_link.pdb_link_record] = []

    for source in input_links:
        label1 = _input_link_label(source.chainID1, source.resSeq1)
        label2 = _input_link_label(source.chainID2, source.resSeq2)
        if label1 not in label_to_idx or label2 not in label_to_idx:
            raise ValueError(
                "Cannot remap preserved LINK because an endpoint residue is absent: "
                f"{source.string.strip()}"
            )
        node1 = nodes[label_to_idx[label1]]
        node2 = nodes[label_to_idx[label2]]

        atom1 = _find_required_atom(
            atom_index,
            node1.new_chain_id,
            node1.new_res_seq,
            _atom_name_candidates(source.name1),
        )
        atom2 = _find_required_atom(
            atom_index,
            node2.new_chain_id,
            node2.new_res_seq,
            _atom_name_candidates(source.name2),
        )
        if atom1 is None or atom2 is None:
            raise ValueError(
                "Cannot remap preserved LINK because an endpoint atom was removed: "
                f"{source.string.strip()}"
            )

        cloned = edit_pdb_link.pdb_link_record(source.string)
        identity_symmetry = (
            (source.sym1.strip() or "1555") == "1555"
            and (source.sym2.strip() or "1555") == "1555"
        )
        # For crystallographic mates, the raw ASU coordinates do not include
        # the LINK symmetry transforms; retain the input distance instead.
        cloned.update_from_atoms(atom1, atom2, recompute_distance=identity_symmetry)
        remapped.append(cloned)

    return remapped


def _link_record_key(
    link: edit_pdb_link.pdb_link_record,
) -> frozenset[Tuple[str, int, str, str]]:
    return frozenset(
        (
            (
                link.chainID1.strip() or " ",
                int(link.resSeq1),
                _normalize_backbone_atom_name(link.name1),
                link.sym1.strip() or "1555",
            ),
            (
                link.chainID2.strip() or " ",
                int(link.resSeq2),
                _normalize_backbone_atom_name(link.name2),
                link.sym2.strip() or "1555",
            ),
        )
    )


def merge_link_records(
    primary: List[edit_pdb_link.pdb_link_record],
    additional: List[edit_pdb_link.pdb_link_record],
) -> List[edit_pdb_link.pdb_link_record]:
    """Merge LINK lists in order while suppressing duplicate atom pairs."""
    merged: List[edit_pdb_link.pdb_link_record] = []
    seen: Set[frozenset[Tuple[str, int, str, str]]] = set()
    for link in list(primary) + list(additional):
        key = _link_record_key(link)
        if key in seen:
            continue
        seen.add(key)
        merged.append(link)
    return merged


# -------------------------- Output renumbering --------------------------


def _chain_id_pool() -> List[str]:
    # PDB chainID is 1-character; include letters and digits for more headroom.
    return list("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789")


def assign_new_labels_and_collect_atoms(
    component_orders: List[Dict[str, object]],
    nodes: List[ResidueNode],
) -> Tuple[List[edit_pdb_atom.pdb_atom_record], Dict[Label, Tuple[str, int]]]:
    """Assign new chainID/resSeq for each node in each component and collect output records.

    This is where we actually "materialize" the final strands into PDB records.

    Returns
    -------
    output_atoms
        List of ATOM/HETATM records in final chain order, **with a TER record
        appended after the last residue of each output chain**.
    phos_new_label
        Mapping original phosphate donor label -> (newChainID, newResSeq) for
        each inserted phosphate-only residue.
    """
    chain_pool = _chain_id_pool()
    if len(component_orders) > len(chain_pool):
        raise ValueError(
            f"Too many resulting chains ({len(component_orders)}) for available chain IDs ({len(chain_pool)})"
        )

    output_atoms: List[edit_pdb_atom.pdb_atom_record] = []
    phos_new_label: Dict[Label, Tuple[str, int]] = {}

    for ci, comp in enumerate(component_orders):
        chain_id = chain_pool[ci]
        order: List[int] = comp["order"]  # type: ignore[index]
        res_counter = 1

        last_atom: Optional[edit_pdb_atom.pdb_atom_record] = None

        for idx in order:
            node = nodes[idx]
            node.new_chain_id = chain_id
            node.new_res_seq = res_counter

            for a in node.atoms:
                a.update_chainID(chain_id)
                a.update_resSeq(res_counter)
                output_atoms.append(a)
                last_atom = a

            if node.is_phos_bridge and node.phos_source is not None:
                phos_new_label[node.phos_source] = (chain_id, res_counter)

            res_counter += 1

        # Append a TER record for this chain, to make chain boundaries explicit.
        if last_atom is not None:
            ter_line = (
                f"TER   {0:5d}      {last_atom.resName:>3s} "
                f"{last_atom.chainID:1s}{last_atom.resSeq:4d}\n"
            )
            output_atoms.append(edit_pdb_atom.pdb_ter_record(ter_line))

    return output_atoms, phos_new_label


# -------------------------- Header REMARK / HET helpers --------------------------

def _clean_remark_value(value: object) -> str:
    """Return a compact value for parse-friendly REMARK key=value fields."""
    text = str(value)
    return text.replace("\n", " ").replace("\r", " ").strip()


def _residue_label(chain_id: str, res_seq: int, res_name: str = "") -> str:
    chain = chain_id if chain_id.strip() else "_"
    if res_name:
        return f"{chain}:{int(res_seq)}:{res_name.strip()}"
    return f"{chain}:{int(res_seq)}"


def _orig_label_text(label: Label, nodes: Optional[List[ResidueNode]] = None, idx: Optional[int] = None) -> str:
    if nodes is not None and idx is not None and 0 <= idx < len(nodes) and nodes[idx].atoms:
        return _residue_label(label[0], label[1], nodes[idx].atoms[0].resName)
    return _residue_label(label[0], label[1])


def _node_output_label(nodes: List[ResidueNode], idx: int) -> str:
    node = nodes[idx]
    res_name = node.atoms[0].resName if node.atoms else (X33_HETID if node.is_phos_bridge else "UNK")
    if node.new_chain_id and node.new_res_seq:
        return _residue_label(node.new_chain_id, node.new_res_seq, res_name)
    return _orig_label_text(node.orig_label(), nodes, idx)


def _linker_phosphate_label_from_mapping(
    phos_key: Label,
    phos_new_label: Dict[Label, Tuple[str, int]],
    linker_phosphate_style: Optional[object] = None,
) -> str:
    style = _coerce_linker_phosphate_style(linker_phosphate_style)
    if phos_key not in phos_new_label:
        return _residue_label(phos_key[0], phos_key[1], style.resname)
    ch, rs = phos_new_label[phos_key]
    return _residue_label(ch, rs, style.resname)


def _x33_label_from_mapping(phos_key: Label, phos_new_label: Dict[Label, Tuple[str, int]]) -> str:
    return _linker_phosphate_label_from_mapping(
        phos_key,
        phos_new_label,
        make_linker_phosphate_style(X33_HETID, LINKER_PHOSPHATE_DEFAULT_RECORD),
    )


def _collect_residue_labels_by_chain(
    atom_rec_list: List[edit_pdb_atom.pdb_atom_record],
) -> Dict[str, List[str]]:
    by_chain: Dict[str, List[str]] = {}
    seen: Set[Tuple[str, int]] = set()
    for rec in atom_rec_list:
        if getattr(rec, "recordName", "") not in ("ATOM", "HETATM"):
            continue
        key = (rec.chainID, rec.resSeq)
        if key in seen:
            continue
        seen.add(key)
        by_chain.setdefault(rec.chainID, []).append(_residue_label(rec.chainID, rec.resSeq, rec.resName))
    return by_chain


def _chunked(items: List[str], chunk_size: int = 24) -> List[List[str]]:
    if not items:
        return []
    return [items[i:i + chunk_size] for i in range(0, len(items), chunk_size)]


def build_chain_residue_remark_lines(
    atom_rec_list: List[edit_pdb_atom.pdb_atom_record],
) -> List[str]:
    """Build parse-friendly chain start/end and residue-list REMARK lines."""
    lines: List[str] = []
    by_chain = _collect_residue_labels_by_chain(atom_rec_list)
    for chain_id in sorted(by_chain.keys(), key=lambda c: (c == " ", c)):
        labels = by_chain[chain_id]
        chain = chain_id if chain_id.strip() else "_"
        if not labels:
            continue
        lines.append(
            f"{REMARK_PREFIX} CHAIN_RANGE chain={chain} start={labels[0]} end={labels[-1]} count={len(labels)}"
        )
        chunks = _chunked(labels)
        for part_index, chunk in enumerate(chunks, start=1):
            lines.append(
                f"{REMARK_PREFIX} CHAIN_RESIDUES chain={chain} part={part_index}/{len(chunks)} "
                f"residues={','.join(chunk)}"
            )
    return lines


def build_linker_phosphate_het_records(
    atom_rec_list: List[edit_pdb_atom.pdb_atom_record],
    linker_phosphate_style: Optional[object] = None,
) -> List[str]:
    """Return HET/HETNAM records for standalone HETATM linker phosphates."""
    style = _coerce_linker_phosphate_style(linker_phosphate_style)
    if style.record_name != "HETATM":
        return []

    counts: Dict[Tuple[str, int], int] = {}
    for rec in atom_rec_list:
        if getattr(rec, "recordName", "") != "HETATM":
            continue
        if getattr(rec, "resName", "").strip() != style.resname:
            continue
        counts[(rec.chainID, rec.resSeq)] = counts.get((rec.chainID, rec.resSeq), 0) + 1

    if not counts:
        return []

    lines: List[str] = []
    for (chain_id, res_seq), atom_count in sorted(counts.items(), key=lambda x: (x[0][0], x[0][1])):
        # PDB-style HET line; kept intentionally simple and parseable.
        lines.append(f"HET    {style.resname:>3s}  {chain_id:1s}{res_seq:4d}     {atom_count:3d}\n")
    lines.append(f"HETNAM     {style.resname:>3s} {LINKER_PHOSPHATE_HETNAM}\n")
    return lines


def build_x33_het_records(atom_rec_list: List[edit_pdb_atom.pdb_atom_record]) -> List[str]:
    """Return HET/HETNAM records for all standalone X33 linker residues."""
    return build_linker_phosphate_het_records(
        atom_rec_list,
        make_linker_phosphate_style(X33_HETID, LINKER_PHOSPHATE_DEFAULT_RECORD),
    )


def build_junction_remark_lines(
    specs: List[Dict[str, object]],
    label_to_idx: Dict[Label, int],
    orig_prev: List[Optional[int]],
    nodes: List[ResidueNode],
    phos_new_label: Optional[Dict[Label, Tuple[str, int]]] = None,
    linker_phosphate_style: Optional[object] = None,
) -> List[str]:
    """Build parse-friendly REMARK lines describing junction residues.

    Residue lists use final output labels when available.  For each operation:
      - double: four residues (prev1,pos1,prev2,pos2)
      - single: two linked residues (prev1,pos2), excluding the two nick ends
      - bowtie: one 3to3 line with five residues including the linker phosphate,
        and one 5to5 line with the four nucleotide residues around the original
        cut sites.
    """
    phos_new_label = phos_new_label or {}
    style = _coerce_linker_phosphate_style(linker_phosphate_style)
    linker_role = "x33" if _is_default_x33_style(style) else "linker_phosphate"
    lines: List[str] = []
    for op_index, sp in enumerate(specs, start=1):
        pos1 = sp["pos1"]  # type: ignore[index]
        pos2 = sp["pos2"]  # type: ignore[index]
        kind = str(sp["kind"]).lower()  # type: ignore[index]
        assert isinstance(pos1, tuple) and isinstance(pos2, tuple)
        if pos1 not in label_to_idx or pos2 not in label_to_idx:
            continue
        idx1 = label_to_idx[pos1]
        idx2 = label_to_idx[pos2]
        u1 = orig_prev[idx1]
        u2 = orig_prev[idx2]
        if u1 is None or u2 is None:
            continue

        original_fields = (
            f"original_prev1={_orig_label_text(nodes[u1].orig_label(), nodes, u1)} "
            f"original_pos1={_orig_label_text(pos1, nodes, idx1)} "
            f"original_prev2={_orig_label_text(nodes[u2].orig_label(), nodes, u2)} "
            f"original_pos2={_orig_label_text(pos2, nodes, idx2)}"
        )

        if kind == "double":
            residues = [_node_output_label(nodes, x) for x in (u1, idx1, u2, idx2)]
            lines.append(
                f"{REMARK_PREFIX} JUNCTION op={op_index} kind=double link=double "
                f"roles=prev1,pos1,prev2,pos2 residues={','.join(residues)} {original_fields}"
            )
        elif kind == "single":
            residues = [_node_output_label(nodes, u1), _node_output_label(nodes, idx2)]
            lines.append(
                f"{REMARK_PREFIX} JUNCTION op={op_index} kind=single link=single "
                f"roles=linked_prev1,linked_pos2 residues={','.join(residues)} {original_fields} "
                f"excluded_nick_ends={_node_output_label(nodes, idx1)},{_node_output_label(nodes, u2)}"
            )
        elif kind == "bowtie":
            linker_label = _linker_phosphate_label_from_mapping(pos2, phos_new_label, style)
            residues_3to3 = [
                _node_output_label(nodes, u1),
                _node_output_label(nodes, idx1),
                _node_output_label(nodes, u2),
                _node_output_label(nodes, idx2),
                linker_label,
            ]
            lines.append(
                f"{REMARK_PREFIX} JUNCTION op={op_index} kind=bowtie link=3to3 "
                f"roles=prev1,pos1,prev2,pos2,{linker_role} residues={','.join(residues_3to3)} "
                f"core={_node_output_label(nodes, u1)},{linker_label},{_node_output_label(nodes, u2)} "
                f"{original_fields}"
            )
            residues_5to5 = [_node_output_label(nodes, x) for x in (u1, idx1, u2, idx2)]
            lines.append(
                f"{REMARK_PREFIX} JUNCTION op={op_index} kind=bowtie link=5to5 "
                f"roles=prev1,pos1,prev2,pos2 residues={','.join(residues_5to5)} "
                f"core={_node_output_label(nodes, idx1)},{_node_output_label(nodes, idx2)} "
                f"{original_fields}"
            )
    return lines


def build_special_remark_lines(
    component_orders: Optional[List[Dict[str, object]]] = None,
    nodes: Optional[List[ResidueNode]] = None,
    phos_new_label: Optional[Dict[Label, Tuple[str, int]]] = None,
    link_counts: Optional[Dict[str, int]] = None,
    specs: Optional[List[Dict[str, object]]] = None,
    linker_phosphate_style: Optional[object] = None,
) -> List[str]:
    """Build parse-friendly REMARK lines for notable topology events."""
    style = _coerce_linker_phosphate_style(linker_phosphate_style)
    lines: List[str] = []

    if link_counts is not None:
        total_links = int(
            link_counts.get(
                "total",
                sum(
                    int(link_counts.get(key, 0))
                    for key in (
                        "backbone_inverted",
                        "bowtie_5to5",
                        "bowtie_3to3",
                        "other_preserved",
                    )
                ),
            )
        )
        lines.append(
            f"{REMARK_PREFIX} SPECIAL event=link_records total={total_links} "
            f"inverted_backbone={link_counts.get('backbone_inverted', 0)} "
            f"bowtie_5to5={link_counts.get('bowtie_5to5', 0)} "
            f"bowtie_3to3={link_counts.get('bowtie_3to3', 0)} "
            f"circular_closure={link_counts.get('circular_closure', 0)} "
            f"other_preserved={link_counts.get('other_preserved', 0)}"
        )

    if component_orders is not None and nodes is not None:
        for strand_index, comp in enumerate(component_orders, start=1):
            order: List[int] = comp.get("order", [])  # type: ignore[assignment]
            chain_id = "?"
            if order:
                chain_id = nodes[order[0]].new_chain_id or "?"
            if comp.get("was_reversed"):
                lines.append(
                    f"{REMARK_PREFIX} SPECIAL event=inverted_strand_direction strand={strand_index} chain={chain_id}"
                )
            if comp.get("is_cycle"):
                rot = comp.get("rotation")
                circularized = bool(comp.get("circularized"))
                lines.append(
                    f"{REMARK_PREFIX} SPECIAL event=circular_component strand={strand_index} chain={chain_id} "
                    f"rotation_start_index={rot} "
                    f"closure_link={'written' if circularized else 'omitted'}"
                )

    if phos_new_label:
        for src, (ch, rs) in sorted(phos_new_label.items(), key=lambda kv: (kv[1][0], kv[1][1])):
            if _is_default_x33_style(style):
                lines.append(
                    f"{REMARK_PREFIX} SPECIAL event=standalone_x33 source={_residue_label(src[0], src[1])} "
                    f"residue={_residue_label(ch, rs, X33_HETID)} atoms=P,OP1,OP2"
                )
            else:
                lines.append(
                    f"{REMARK_PREFIX} SPECIAL event=standalone_linker_phosphate "
                    f"source={_residue_label(src[0], src[1])} "
                    f"residue={_residue_label(ch, rs, style.resname)} "
                    f"record={style.record_name} atoms=P,OP1,OP2"
                )

    if specs:
        for op_index, sp in enumerate(specs, start=1):
            if str(sp.get("kind", "")).lower() == "bowtie":
                pos1 = sp.get("pos1")
                pos2 = sp.get("pos2")
                if isinstance(pos1, tuple) and isinstance(pos2, tuple):
                    if _is_default_x33_style(style):
                        source_field = f"x33_source={_residue_label(pos2[0], pos2[1])}"
                    else:
                        source_field = (
                            f"linker_phosphate_source={_residue_label(pos2[0], pos2[1])} "
                            f"linker_phosphate_resname={style.resname} "
                            f"linker_phosphate_record={style.record_name}"
                        )
                    lines.append(
                        f"{REMARK_PREFIX} SPECIAL event=bowtie_junction op={op_index} "
                        f"pos1={_residue_label(pos1[0], pos1[1])} pos2={_residue_label(pos2[0], pos2[1])} "
                        f"{source_field}"
                    )

    return lines


def build_re_script_header_lines(
    software_name: str,
    software_version: str,
    developer: str,
    command: Optional[str],
    output_stage: str,
    atom_rec_list: List[edit_pdb_atom.pdb_atom_record],
    specs: Optional[List[Dict[str, object]]] = None,
    label_to_idx: Optional[Dict[Label, int]] = None,
    orig_prev: Optional[List[Optional[int]]] = None,
    nodes: Optional[List[ResidueNode]] = None,
    phos_new_label: Optional[Dict[Label, Tuple[str, int]]] = None,
    component_orders: Optional[List[Dict[str, object]]] = None,
    link_counts: Optional[Dict[str, int]] = None,
    extra_special_events: Optional[List[str]] = None,
    linker_phosphate_style: Optional[object] = None,
) -> List[str]:
    """Build standardized RE_SCRIPT REMARK 950 header lines."""
    style = _coerce_linker_phosphate_style(linker_phosphate_style)
    lines: List[str] = []
    lines.append(
        f"{REMARK_PREFIX} SOFTWARE name={_clean_remark_value(software_name)} "
        f"version={_clean_remark_value(software_version)} developer={_clean_remark_value(developer)}"
    )
    if command:
        lines.append(f"{REMARK_PREFIX} COMMAND text={_clean_remark_value(command)}")
    lines.append(f"{REMARK_PREFIX} OUTPUT_STAGE name={_clean_remark_value(output_stage)}")
    lines.extend(build_chain_residue_remark_lines(atom_rec_list))

    if specs is not None and label_to_idx is not None and orig_prev is not None and nodes is not None:
        lines.extend(
            build_junction_remark_lines(
                specs,
                label_to_idx,
                orig_prev,
                nodes,
                phos_new_label,
                linker_phosphate_style=style,
            )
        )

    lines.extend(
        build_special_remark_lines(
            component_orders,
            nodes,
            phos_new_label,
            link_counts,
            specs,
            linker_phosphate_style=style,
        )
    )

    if extra_special_events:
        for event in extra_special_events:
            lines.append(f"{REMARK_PREFIX} SPECIAL event={_clean_remark_value(event)}")

    return [line if line.endswith("\n") else line + "\n" for line in lines]


def write_pdb_with_header(
    atom_rec_list: List[edit_pdb_atom.pdb_atom_record],
    link_rec_list: List[edit_pdb_link.pdb_link_record],
    outfile,
    header_lines: Optional[List[str]] = None,
    reorder_serial: bool = False,
    linker_phosphate_style: Optional[object] = None,
) -> None:
    """Write REMARK/HET/HETNAM records, then LINK and ATOM/HETATM/TER records."""
    style = _coerce_linker_phosphate_style(linker_phosphate_style)
    for line in header_lines or []:
        outfile.write(line if line.endswith("\n") else line + "\n")
    for line in build_linker_phosphate_het_records(atom_rec_list, style):
        outfile.write(line)
    edit_pdb_link.rec2file_link(atom_rec_list, link_rec_list, outfile, reorder_serial=reorder_serial)


# -------------------------- Main --------------------------


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Apply reciprocal exchanges (double/single) and bowtie exchanges to a DNA PDB."
    )
    ap.add_argument("pdbfile", help="Input PDB file")
    ap.add_argument(
        "exchanges",
        nargs="+",
        help='Exchange specs as repeated triples: "pos1 pos2 kind" (e.g., 9C 23A double 23C 23F B)',
    )
    ap.add_argument("-o", "--out", default=None, help="Output PDB filename")
    ap.add_argument(
        "--cir_shift",
        type=parse_cir_shift,
        default=CirShiftSpec(),
        help=(
            "Exact signed residue rotation, modulo cycle length, used to choose "
            "the numbering start for cyclic components; zero performs no "
            "circular permutation. "
            "Append c/C (for example, 10c) to write a closing LINK and keep "
            "each cyclic strand covalently circularized; plain 10 leaves the "
            "serialized break open (default: 8)."
        ),
    )
    ap.add_argument(
        "--min_link_records",
        "--min-link-records",
        action="store_true",
        help=(
            "Choose each strand direction primarily to minimize topology LINK "
            "records. Cycles are scored after the exact cir_shift and in the "
            "requested open or closed mode. This can reverse whole output "
            "strands. By default, terminal input-fragment direction is preserved."
        ),
    )
    ap.add_argument(
        "--linker_phosphate_resname",
        "--linker-phosphate-resname",
        default=None,
        help=(
            "Residue name for phosphate-only residues inserted at bowtie 3'-3' linkages "
            "(default: X33). Use DA/dA for regular ATOM DA output."
        ),
    )
    ap.add_argument(
        "--linker_phosphate_record",
        "--linker-phosphate-record",
        choices=["ATOM", "HETATM", "atom", "hetatm"],
        default=None,
        help=(
            "Record type for inserted 3'-3' linker phosphates. Default is HETATM, "
            "except DA/dA defaults to ATOM."
        ),
    )

    args = ap.parse_args()
    cir_shift_spec = coerce_cir_shift(args.cir_shift)
    try:
        linker_phosphate_style = make_linker_phosphate_style(
            args.linker_phosphate_resname,
            args.linker_phosphate_record,
        )
    except ValueError as exc:
        ap.error(str(exc))

    # Parse exchange specs.
    specs = parse_exchange_specs(args.exchanges)
    bowtie_specs = [sp for sp in specs if sp["kind"] == "bowtie"]

    # Read PDB.
    rec_list: list[edit_pdb_atom.pdb_record] = []
    input_link_records: List[edit_pdb_link.pdb_link_record] = []
    with open(args.pdbfile) as fin:
        edit_pdb_link.file2rec_link(fin, rec_list, input_link_records)
    atom_recs = [r for r in rec_list if isinstance(r, edit_pdb_atom.pdb_atom_record)]
    nodes, label_to_idx = build_residue_nodes(atom_recs)
    input_topology = build_input_backbone_topology(nodes, label_to_idx, input_link_records)
    g = input_topology.graph
    orig_prev = input_topology.orig_prev

    if input_link_records:
        print(
            "Input LINK topology: "
            f"{input_topology.recognized_link_count} backbone record(s), "
            f"{input_topology.phosphate_bridge_count} phosphate bridge(s), "
            f"{len(input_topology.passthrough_links)} other record(s)."
        )

    # Store/cut bowtie phosphates (pos2 only).
    phos_store = cut_and_store_bowtie_phosphates(nodes, label_to_idx, bowtie_specs)

    # Apply exchanges to the reconstructed input backbone graph.
    junction_nodes, n_double, n_single, n_bowtie = apply_exchanges_to_graph(
        g, nodes, label_to_idx, orig_prev, specs
    )

    print("\nLINK records will be written for:")
    print("  - bowtie 3'-3': P(phosphate-only residue) -- O3' (each side)")
    print("  - bowtie 5'-5': O5'(pos2) -- P(pos1)")
    print("  - inverted backbone steps: P(res i) -- O3'(res j) when traversing P->O3 along a standard bond")

    print("\nExchange summary:")
    print(f"  double: {n_double}")
    print(f"  single: {n_single}")
    print(f"  bowtie: {n_bowtie}")
    if n_bowtie:
        print(
            "  3'-3' linker phosphate style: "
            f"{linker_phosphate_style.record_name} {linker_phosphate_style.resname}"
        )

    # Build ordered components.
    base_components = build_ordered_components(
        g,
        nodes,
        junction_nodes,
        cir_shift_spec.shift,
        circularize_cycles=cir_shift_spec.circularize,
        min_link_records=args.min_link_records,
    )

    # Insert phosphate-only residues for 3to3 edges.
    used_phos: Set[Label] = set()
    final_components: List[Dict[str, object]] = []
    for ci, comp in enumerate(base_components):
        base_order: List[int] = comp["order"]  # type: ignore[index]
        expanded_order = insert_phosphate_nodes(
            base_order,
            g,
            nodes,
            phos_store,
            used_phos,
            linker_phosphate_style=linker_phosphate_style,
            is_cycle=bool(comp.get("circularized")),
        )
        final_components.append({
            **comp,
            "order": expanded_order,
        })

    # Debug-print linking paths.
    print("\nResulting strands (original labels; '->' normal, '->>' LINK-required):")
    chain_pool = _chain_id_pool()
    for i, comp in enumerate(final_components):
        order: List[int] = comp["order"]  # type: ignore[index]
        is_cycle: bool = comp["is_cycle"]  # type: ignore[index]
        rot = comp.get("rotation")
        was_rev = comp.get("was_reversed")
        new_chain_id = chain_pool[i] if i < len(chain_pool) else "?"
        header = f"  Strand {i+1} (new chain {new_chain_id}): len={len(order)}"
        if is_cycle:
            header += " (cycle)"
        if rot is not None and is_cycle:
            header += f" [rot_start_idx={rot}]"
        if was_rev:
            header += " [reversed]"
        print(header)
        print("    " + format_path_string(order, g, nodes))

    # Assign new labels, collect atoms.
    output_atoms, phos_new_label = assign_new_labels_and_collect_atoms(final_components, nodes)

    # Generate LINK records.
    link_records, link_counts = build_link_records(final_components, g, nodes, output_atoms)
    passthrough_links = remap_passthrough_link_records(
        input_topology.passthrough_links,
        nodes,
        label_to_idx,
        output_atoms,
    )
    topology_link_count = len(link_records)
    link_records = merge_link_records(link_records, passthrough_links)
    link_counts["other_preserved"] = len(link_records) - topology_link_count
    link_counts["total"] = len(link_records)

    command_text = " ".join(shlex.quote(arg) for arg in [sys.executable] + sys.argv)
    header_lines = build_re_script_header_lines(
        software_name=SOFTWARE_NAME,
        software_version=SOFTWARE_VERSION,
        developer=SOFTWARE_DEVELOPER,
        command=command_text,
        output_stage="reciprocal_exchange",
        atom_rec_list=output_atoms,
        specs=specs,
        label_to_idx=label_to_idx,
        orig_prev=orig_prev,
        nodes=nodes,
        phos_new_label=phos_new_label,
        component_orders=final_components,
        link_counts=link_counts,
        extra_special_events=[
            "orientation_policy "
            f"mode={'min_link_records' if args.min_link_records else 'preserve_terminal_direction'} "
            f"min_link_records={1 if args.min_link_records else 0}"
        ],
        linker_phosphate_style=linker_phosphate_style,
    )

    print("\nOutput summary:")
    print(f"  resulting chains: {len(final_components)}")
    print(
        "  orientation policy: "
        + ("minimize LINK records" if args.min_link_records else "preserve terminal direction")
    )
    print(f"  LINK records: {len(link_records)}")
    print(
        f"    inverted-backbone: {link_counts['backbone_inverted']}, "
        f"bowtie 5'-5': {link_counts['bowtie_5to5']}, "
        f"bowtie 3'-3': {link_counts['bowtie_3to3']}, "
        f"cycle closures: {link_counts['circular_closure']}, "
        f"other preserved: {link_counts['other_preserved']}"
    )

    # Report where each 3'-3' phosphate ended up (new labels).
    for src, (ch, rs) in sorted(phos_new_label.items(), key=lambda kv: (kv[1][0], kv[1][1])):
        print(f"  3'-3' phosphate from original {src[1]}{src[0]} is now residue {rs}{ch}.")

    # Output file name.
    outname = args.out
    if outname is None:
        base = args.pdbfile
        if base.lower().endswith(".pdb"):
            base = base[:-4]
        elif base.lower().endswith(".pdb.txt"):
            base = base[:-8]
        outname = base + "_rex.pdb"

    with open(outname, 'w') as fout:
        write_pdb_with_header(
            output_atoms,
            link_records,
            fout,
            header_lines=header_lines,
            reorder_serial=True,
            linker_phosphate_style=linker_phosphate_style,
        )
    print(f"\nWrote output to: {outname}")


if __name__ == "__main__":
    main()
