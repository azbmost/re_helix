#!/usr/bin/env python3
"""reciprocal_exchange_pdbV3_3.py

Reciprocal exchange (double / single) and bowtie exchange for DNA PDBs.

Key behavior (V3.3):
  1) Build the ORIGINAL 5'->3' backbone graph from residue numbering for each chain.
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
  4) Standalone X33 linker residues:
     - For each bowtie, we cut the phosphate-group atoms (P + non-bridging O's)
       from pos2 and store them (keyed by the ORIGINAL pos2 label, e.g., 23F).
     - When reconstructing final strand paths, we insert a standalone HETATM
       linker residue named X33 between the two residues that form the 3'-3'
       edge. X33 contains exactly three atoms: P, OP1, and OP2.
  5) LINK records:
     - We ignore CONECT entirely.
     - We write LINK records for:
         (a) every bowtie 3'-3' linkage: P(phosphate-only residue) -- O3' (each side)
         (b) every bowtie 5'-5' linkage: O5'(pos2) -- P(pos1)
         (c) every inverted backbone step (i.e., when traversing a standard O3'--P bond
             in the P->O3 direction). Natural O3->P steps rely on standard connectivity.

Additionally:
  - Cyclic components are supported. We perform a circular permutation (rotation)
    for output numbering, similar to V2/V3. We attempt to choose the break point
    away from junction residues.
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
SOFTWARE_VERSION = "V3.3"
SOFTWARE_DEVELOPER = "DiLiuLab"
REMARK_PREFIX = "REMARK 950 RE_SCRIPT"
X33_HETID = "X33"
X33_HETNAM = "3'-3' PHOSPHODIESTER LINKER PHOSPHATE"


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


# -------------------------- Parsing helpers --------------------------


_RES_RE = re.compile(r"^(\d+)([A-Za-z])$")


def parse_res_label(token: str) -> Label:
    m = _RES_RE.match(token.strip())
    if not m:
        raise ValueError(f"Invalid residue label token: '{token}' (expected like 23A)")
    res_seq = int(m.group(1))
    chain_id = m.group(2).upper()
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
    for a in atom_recs:
        if a.recordName not in ("ATOM", "HETATM"):
            continue
        key = (a.chainID.strip() or " ", a.resSeq)
        key = (key[0].upper(), key[1])
        grouped.setdefault(key, []).append(a)

    nodes: List[ResidueNode] = []
    label_to_idx: Dict[Label, int] = {}
    for key in sorted(grouped.keys(), key=lambda x: (x[0], x[1])):
        chain_id, res_seq = key
        idx = len(nodes)
        nodes.append(ResidueNode(orig_chain_id=chain_id, orig_res_seq=res_seq, atoms=grouped[key]))
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


def _canonicalize_x33_phosphate_atoms(
    phos_atoms: List[edit_pdb_atom.pdb_atom_record],
    source_label: Label,
) -> List[edit_pdb_atom.pdb_atom_record]:
    """Return exactly P, OP1, OP2 atoms for a standalone X33 linker residue.

    The phosphate donor can use either OP1/OP2 or O1P/O2P naming.  We keep
    the coordinates of the donor phosphate atoms, rename the non-bridging
    oxygens to OP1/OP2, and mark all three records as HETATM X33.  Extra
    phosphate atoms, if any, are intentionally not carried into X33 because
    this custom linker residue is defined as a three-atom HET group.
    """
    by_name: Dict[str, edit_pdb_atom.pdb_atom_record] = {}
    for atom in phos_atoms:
        by_name.setdefault(atom.name.strip().upper(), atom)

    p_atom = by_name.get("P")
    if p_atom is None:
        raise ValueError(f"X33 linker from {source_label[1]}{source_label[0]} is missing atom P")

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
            f"X33 linker from {source_label[1]}{source_label[0]} needs two non-bridging O atoms"
        )

    canonical = [(p_atom, "P"), (op1_atom, "OP1"), (op2_atom, "OP2")]
    for atom, atom_name in canonical:
        atom.update_recordName("HETATM")
        atom.update_resName(X33_HETID)
        _set_atom_name(atom, atom_name)

    return [p_atom, op1_atom, op2_atom]


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
            raise ValueError(f"Residue {pos1[1]}{pos1[0]} has no original predecessor; cannot cut incoming edge")
        if u2 is None:
            raise ValueError(f"Residue {pos2[1]}{pos2[0]} has no original predecessor; cannot cut incoming edge")

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


def _inverted_cost(order: List[int], g: BackboneGraph) -> int:
    """Count how many std edges are traversed in the inverted P->O3 direction."""
    cost = 0
    for i in range(len(order) - 1):
        a, b = order[i], order[i + 1]
        e = g.get_edge(a, b)
        if e.kind != "std":
            continue
        end_a, end_b = e.endpoints(a, b)
        if end_a == "P" and end_b == "O3":
            cost += 1
    return cost


def _end_fragment_dir(order: List[int], nodes: List[ResidueNode], from_start: bool) -> int:
    """Infer direction (+1 nat, -1 inv) of the first/last fragment in a path.

    We look for the first (or last) consecutive pair that belongs to the SAME original chain.
    """
    if len(order) < 2:
        return 0
    if from_start:
        rng = range(len(order) - 1)
    else:
        rng = range(len(order) - 2, -1, -1)
    for i in rng:
        a = order[i]
        b = order[i + 1]
        na, nb = nodes[a], nodes[b]
        if na.orig_chain_id != nb.orig_chain_id:
            continue
        d = nb.orig_res_seq - na.orig_res_seq
        if d > 0:
            return 1
        if d < 0:
            return -1
    return 0


def orient_path_component(order: List[int], nodes: List[ResidueNode], g: BackboneGraph) -> Tuple[List[int], bool]:
    """Choose orientation for a path component.

    Returns:
        (oriented_order, was_reversed)

    Heuristics:
      - choose the direction with fewer inverted std edges.
      - if both ends would be 3'->5' (dir=-1 on both ends), reverse.
    """
    if not order:
        return order, False
    rev = list(reversed(order))
    c1 = _inverted_cost(order, g)
    c2 = _inverted_cost(rev, g)
    chosen = order if c1 <= c2 else rev

    first_dir = _end_fragment_dir(chosen, nodes, from_start=True)
    last_dir = _end_fragment_dir(chosen, nodes, from_start=False)
    if first_dir == -1 and last_dir == -1:
        chosen = list(reversed(chosen))
        return chosen, True

    return chosen, (chosen is rev)


def choose_cycle_orientation(order: List[int], g: BackboneGraph) -> List[int]:
    """Choose cycle direction that minimizes inverted std edges."""
    if not order:
        return order
    rev = list(reversed(order))
    c1 = _inverted_cost(order, g)
    c2 = _inverted_cost(rev, g)
    return order if c1 <= c2 else rev


def rotate_cycle_away_from_junctions(
    order: List[int],
    g: BackboneGraph,
    junction_nodes: Set[int],
    cir_shift: int,
) -> Tuple[List[int], int]:
    """Rotate a cycle order list to choose a start index.

    We try to start about cir_shift residues away from junction nodes.
    Additionally, we prefer that the break edge (last->first) is a STANDARD edge,
    not a 3to3/5to5 edge (to avoid hiding special linkages at the break).

    Returns:
        (rotated_order, start_index_in_original_order)
    """
    n = len(order)
    if n == 0:
        return order, 0

    start0 = cir_shift % n

    def is_good_start(i: int) -> bool:
        first = order[i]
        prev = order[(i - 1) % n]
        if first in junction_nodes or prev in junction_nodes:
            return False
        e = g.get_edge(prev, first)
        if e.kind != "std":
            return False
        return True

    # First pass: strict.
    for off in range(n):
        i = (start0 + off) % n
        if is_good_start(i):
            rotated = order[i:] + order[:i]
            return rotated, i

    # Second pass: ignore junction avoidance, but still avoid breaking at special edges.
    for off in range(n):
        i = (start0 + off) % n
        prev = order[(i - 1) % n]
        first = order[i]
        e = g.get_edge(prev, first)
        if e.kind == "std":
            rotated = order[i:] + order[:i]
            return rotated, i

    # Fallback: just rotate by start0.
    rotated = order[start0:] + order[:start0]
    return rotated, start0


def build_ordered_components(
    g: BackboneGraph,
    nodes: List[ResidueNode],
    junction_nodes: Set[int],
    cir_shift: int,
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
            order, was_rev = orient_path_component(order, nodes, g)
            components.append({
                "order": order,
                "is_cycle": False,
                "was_reversed": was_rev,
                "rotation": None,
            })
        else:
            # Cycle: traverse starting from min node with an arbitrary first neighbor.
            start_cycle = min(comp_nodes)
            nbs = g.neigh[start_cycle]
            if len(nbs) != 2:
                raise ValueError("Cycle node does not have 2 neighbors")
            order = _traverse_cycle(g, start_cycle, nbs[0])
            # Orientation and rotation.
            order = choose_cycle_orientation(order, g)
            order, rot_idx = rotate_cycle_away_from_junctions(order, g, junction_nodes, cir_shift)
            components.append({
                "order": order,
                "is_cycle": True,
                "was_reversed": False,
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
) -> List[int]:
    """Expand 3to3 edges by inserting phosphate-only residue nodes."""
    if not base_order:
        return base_order

    out: List[int] = []
    for i in range(len(base_order) - 1):
        a = base_order[i]
        b = base_order[i + 1]
        out.append(a)
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

            # Create a new standalone X33 HETATM linker residue.
            x33_atoms = _canonicalize_x33_phosphate_atoms(phos_atoms, phos_key)
            new_idx = len(nodes)
            nodes.append(
                ResidueNode(
                    orig_chain_id=phos_key[0],
                    orig_res_seq=phos_key[1],
                    atoms=x33_atoms,
                    is_phos_bridge=True,
                    phos_source=phos_key,
                    no_phosphate=False,
                )
            )
            used_phos.add(phos_key)
            out.append(new_idx)

    out.append(base_order[-1])
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
        "backbone_inverted": 0,
        "bowtie_5to5": 0,
        "bowtie_3to3": 0,
    }

    for comp in component_orders:
        order: List[int] = comp["order"]  # type: ignore[index]
        for i in range(len(order) - 1):
            a_idx, b_idx = order[i], order[i + 1]
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
                continue

            if e.kind != "std":
                # 3to3 edges should have been expanded by inserting phosphate nodes.
                raise ValueError(
                    f"Unexpected non-std edge kind '{e.kind}' between residues in final order; "
                    f"missing phosphate insertion?"
                )

            # Standard edge: add LINK only when traversal is inverted (P->O3).
            end_a, end_b = e.endpoints(a_idx, b_idx)
            if end_a == "P" and end_b == "O3":
                p_atom = _find_required_atom(atom_index, a_node.new_chain_id, a_node.new_res_seq, ["P"])
                o3_atom = _find_required_atom(atom_index, b_node.new_chain_id, b_node.new_res_seq, ["O3'", "O3*"])
                if p_atom is None or o3_atom is None:
                    print(
                        f"Warning: could not build inverted-backbone LINK for adjacency "
                        f"{a_node.new_res_seq}{a_node.new_chain_id}->{b_node.new_res_seq}{b_node.new_chain_id} "
                        f"(missing P or O3')",
                        file=sys.stderr,
                    )
                    continue
                links.append(edit_pdb_link.pdb_link_record(_format_link_line(p_atom, o3_atom)))
                counts["backbone_inverted"] += 1

    return links, counts


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


def _x33_label_from_mapping(phos_key: Label, phos_new_label: Dict[Label, Tuple[str, int]]) -> str:
    if phos_key not in phos_new_label:
        return _residue_label(phos_key[0], phos_key[1], X33_HETID)
    ch, rs = phos_new_label[phos_key]
    return _residue_label(ch, rs, X33_HETID)


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


def build_x33_het_records(atom_rec_list: List[edit_pdb_atom.pdb_atom_record]) -> List[str]:
    """Return HET/HETNAM records for all standalone X33 linker residues."""
    counts: Dict[Tuple[str, int], int] = {}
    for rec in atom_rec_list:
        if getattr(rec, "recordName", "") != "HETATM":
            continue
        if getattr(rec, "resName", "").strip() != X33_HETID:
            continue
        counts[(rec.chainID, rec.resSeq)] = counts.get((rec.chainID, rec.resSeq), 0) + 1

    if not counts:
        return []

    lines: List[str] = []
    for (chain_id, res_seq), atom_count in sorted(counts.items(), key=lambda x: (x[0][0], x[0][1])):
        # PDB-style HET line; kept intentionally simple and parseable.
        lines.append(f"HET    {X33_HETID:>3s}  {chain_id:1s}{res_seq:4d}     {atom_count:3d}\n")
    lines.append(f"HETNAM     {X33_HETID:>3s} {X33_HETNAM}\n")
    return lines


def build_junction_remark_lines(
    specs: List[Dict[str, object]],
    label_to_idx: Dict[Label, int],
    orig_prev: List[Optional[int]],
    nodes: List[ResidueNode],
    phos_new_label: Optional[Dict[Label, Tuple[str, int]]] = None,
) -> List[str]:
    """Build parse-friendly REMARK lines describing junction residues.

    Residue lists use final output labels when available.  For each operation:
      - double: four residues (prev1,pos1,prev2,pos2)
      - single: two linked residues (prev1,pos2), excluding the two nick ends
      - bowtie: one 3to3 line with five residues including X33, and one 5to5
        line with the four nucleotide residues around the original cut sites.
    """
    phos_new_label = phos_new_label or {}
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
            x33_label = _x33_label_from_mapping(pos2, phos_new_label)
            residues_3to3 = [
                _node_output_label(nodes, u1),
                _node_output_label(nodes, idx1),
                _node_output_label(nodes, u2),
                _node_output_label(nodes, idx2),
                x33_label,
            ]
            lines.append(
                f"{REMARK_PREFIX} JUNCTION op={op_index} kind=bowtie link=3to3 "
                f"roles=prev1,pos1,prev2,pos2,x33 residues={','.join(residues_3to3)} "
                f"core={_node_output_label(nodes, u1)},{x33_label},{_node_output_label(nodes, u2)} "
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
) -> List[str]:
    """Build parse-friendly REMARK lines for notable topology events."""
    lines: List[str] = []

    if link_counts is not None:
        total_links = sum(int(v) for v in link_counts.values())
        lines.append(
            f"{REMARK_PREFIX} SPECIAL event=link_records total={total_links} "
            f"inverted_backbone={link_counts.get('backbone_inverted', 0)} "
            f"bowtie_5to5={link_counts.get('bowtie_5to5', 0)} "
            f"bowtie_3to3={link_counts.get('bowtie_3to3', 0)}"
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
                lines.append(
                    f"{REMARK_PREFIX} SPECIAL event=circular_component strand={strand_index} chain={chain_id} "
                    f"rotation_start_index={rot}"
                )

    if phos_new_label:
        for src, (ch, rs) in sorted(phos_new_label.items(), key=lambda kv: (kv[1][0], kv[1][1])):
            lines.append(
                f"{REMARK_PREFIX} SPECIAL event=standalone_x33 source={_residue_label(src[0], src[1])} "
                f"residue={_residue_label(ch, rs, X33_HETID)} atoms=P,OP1,OP2"
            )

    if specs:
        for op_index, sp in enumerate(specs, start=1):
            if str(sp.get("kind", "")).lower() == "bowtie":
                pos1 = sp.get("pos1")
                pos2 = sp.get("pos2")
                if isinstance(pos1, tuple) and isinstance(pos2, tuple):
                    lines.append(
                        f"{REMARK_PREFIX} SPECIAL event=bowtie_junction op={op_index} "
                        f"pos1={_residue_label(pos1[0], pos1[1])} pos2={_residue_label(pos2[0], pos2[1])} "
                        f"x33_source={_residue_label(pos2[0], pos2[1])}"
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
) -> List[str]:
    """Build standardized RE_SCRIPT REMARK 950 header lines."""
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
        lines.extend(build_junction_remark_lines(specs, label_to_idx, orig_prev, nodes, phos_new_label))

    lines.extend(build_special_remark_lines(component_orders, nodes, phos_new_label, link_counts, specs))

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
) -> None:
    """Write REMARK/HET/HETNAM records, then LINK and ATOM/HETATM/TER records."""
    for line in header_lines or []:
        outfile.write(line if line.endswith("\n") else line + "\n")
    for line in build_x33_het_records(atom_rec_list):
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
        type=int,
        default=8,
        help="Shift (in residues) used to choose a numbering start for cyclic components",
    )

    args = ap.parse_args()

    # Parse exchange specs.
    specs = parse_exchange_specs(args.exchanges)
    bowtie_specs = [sp for sp in specs if sp["kind"] == "bowtie"]

    # Read PDB.
    rec_list: list[edit_pdb_atom.pdb_record] = []
    with open(args.pdbfile) as fin:
        edit_pdb_atom.file2rec(fin, rec_list)
    atom_recs = [r for r in rec_list if isinstance(r, edit_pdb_atom.pdb_atom_record)]
    nodes, label_to_idx = build_residue_nodes(atom_recs)
    orig_prev, orig_next = build_original_prev_next(nodes)

    # Store/cut bowtie phosphates (pos2 only).
    phos_store = cut_and_store_bowtie_phosphates(nodes, label_to_idx, bowtie_specs)

    # Build original backbone graph and apply exchanges.
    g = build_original_graph(len(nodes), orig_next)
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

    # Build ordered components.
    base_components = build_ordered_components(g, nodes, junction_nodes, args.cir_shift)

    # Insert phosphate-only residues for 3to3 edges.
    used_phos: Set[Label] = set()
    final_components: List[Dict[str, object]] = []
    for ci, comp in enumerate(base_components):
        base_order: List[int] = comp["order"]  # type: ignore[index]
        expanded_order = insert_phosphate_nodes(base_order, g, nodes, phos_store, used_phos)
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
    )

    print("\nOutput summary:")
    print(f"  resulting chains: {len(final_components)}")
    print(f"  LINK records: {len(link_records)}")
    print(
        f"    inverted-backbone: {link_counts['backbone_inverted']}, "
        f"bowtie 5'-5': {link_counts['bowtie_5to5']}, "
        f"bowtie 3'-3': {link_counts['bowtie_3to3']}"
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
        write_pdb_with_header(output_atoms, link_records, fout, header_lines=header_lines, reorder_serial=True)
    print(f"\nWrote output to: {outname}")


if __name__ == "__main__":
    main()
