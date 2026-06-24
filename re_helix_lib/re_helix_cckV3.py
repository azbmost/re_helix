#!/usr/bin/env python3
"""
re_helix_cckV3.py

Closure-residual-based cyclic-constraint alignment for reciprocal-exchange
helices, inspired by HolT Hunter. V3 keeps the post-alignment
twist_rod / twist_helix diagnostic report in which twist_helix is measured
directly from reciprocal-exchange-site phosphate rotations around the helix
axis.

Compared with re_helix_ccgV3.py, this version:

  * keeps the existing replication, helix detection, tree / pairwise alignment,
    and reciprocal-exchange machinery from the current re_helix.py helper
    module, falling back to legacy local helper modules when available;

  * replaces the cyclic-component weighted geometric least-squares objective as
    the *primary* optimisation target with an explicit closure residual,
    closer in spirit to HolT Hunter's root solve;

  * uses the geometry score (P–P distances + line topology + axis-distance
    mismatch) only as a secondary ranking / filter among closure solutions.

Pipeline for cyclic components when --axis_parallel n
-----------------------------------------------------
1) Build a BFS spanning tree over the helix graph.
2) Perform a base tree alignment in which each tree edge is aligned by
   optimising only (d, theta, phi) with rho fixed at 0.
3) For each tree edge, store the local phi axis and rho axis.
4) Solve a continuous closure problem over per-tree-edge (phi_offset, rho)
   variables:
      - primary target: explicit closure residual on cycle edges,
      - secondary ranking: geometry score on the whole component.
5) By default, run a final geometry-based cyclic polish (from
   re_helix_ccgV3.py) starting from the chosen closure basin. This
   keeps the HolT-Hunter-like closure search, but finishes with the stronger
   global geometric residual that has been working better in practice.

The closure residual uses, for each cycle-edge P-pair:
  * a 3-vector equal to the Cartesian difference between the two transformed
    constrained P atoms, and
  * one scalar axis-distance residual per cycle edge.

If the closure residual dimension matches the number of variables, scipy.root
is used (HolT-Hunter-like). Otherwise, scipy.least_squares is used on the same
closure residual. Candidate roots are then ranked primarily by closure norm and
secondarily by the geometry score.

Outputs
-------
  <base>_aligned_cck.pdb      : aligned structure before reciprocal exchange
  <base>_aligned_cck_rex.pdb  : aligned structure after reciprocal exchange
  <base>_aligned_cck_twist.tsv : post-alignment twist_rod / twist_helix report

Recent fixes / diagnostics
--------------------------
The earlier cck script evaluated closure candidates against a frozen snapshot,
but accidentally applied the chosen transforms only to that detached snapshot
instead of the live rec_list that gets written out. This version keeps the
live-structure application fix and the optional geometry polish.

The twist report is a post-processing diagnostic. It follows the HolT-Hunter
rod-versus-helix comparison only at the diagnostic level: twist_rod is inferred
from the co-perpendicular connector vectors between aligned helix axes, while
twist_helix is measured directly from the reciprocal-exchange-site phosphate
atoms around the helix axis. The script does not perform a separate exhaustive
ideal triangular-rod search like HolT Hunter.
"""

import argparse
import importlib.util
import math
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
from scipy.optimize import least_squares, root

SOFTWARE_NAME = "re_helix_cck"
SOFTWARE_VERSION = "V3"


# ---------------------------------------------------------------------------
# Load project modules (some filenames contain dots)
# ---------------------------------------------------------------------------

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent
for path in (str(PARENT), str(HERE)):
    if path not in sys.path:
        sys.path.insert(0, path)


from re_helix_lib.edit_pdb_atom import file2rec, rec2file, pdb_atom_record  # type: ignore
from re_helix_lib.edit_pdb_link import rec2file_link  # type: ignore


def _load_module(module_name: str, filenames: List[str]):
    for fname in filenames:
        path = HERE / fname
        if not path.exists():
            continue
        spec = importlib.util.spec_from_file_location(module_name, str(path))
        if spec is None or spec.loader is None:
            raise ImportError("Could not create import spec for %s" % path)
        mod = importlib.util.module_from_spec(spec)
        sys.modules[module_name] = mod
        spec.loader.exec_module(mod)  # type: ignore[attr-defined]
        return mod
    raise ImportError("Could not find any of: %s" % ", ".join(filenames))


alignmod = _load_module(
    "re_helix_current_mod",
    ["../re_helix.py"],
)

# Aliases to the latest project helpers.
HelixID = alignmod.HelixID
parse_exchange_specs = alignmod.parse_exchange_specs
parse_helix_definition_tokens = alignmod.parse_helix_definition_tokens
build_nucleic_acid_maps = alignmod.build_nucleic_acid_maps
compute_chain_partner_map = alignmod.compute_chain_partner_map
build_chain_to_helix_from_defs = alignmod.build_chain_to_helix_from_defs
build_helix_pair_graph = alignmod.build_helix_pair_graph
compute_helix_axis = alignmod.compute_helix_axis
apply_reciprocal_exchanges_in_memory = alignmod.apply_reciprocal_exchanges_in_memory
replicate_all_chains = alignmod.replicate_all_chains
helix_id_str = alignmod.helix_id_str
v_add = alignmod.v_add
v_sub = alignmod.v_sub
v_scale = alignmod.v_scale
v_dot = alignmod.v_dot
v_cross = alignmod.v_cross
v_length = alignmod.v_length
v_norm = alignmod.v_norm
rotate_around_line = alignmod.rotate_around_line
align_axes_for_pair = alignmod.align_axes_for_pair
build_pair_objective = alignmod.build_pair_objective
coordinate_descent = alignmod.coordinate_descent
apply_transform_to_helix = alignmod.apply_transform_to_helix
select_pairs_for_alignment = alignmod.select_pairs_for_alignment


def _load_ccg_module():
    """Load the bundled geometry-based cyclic polish module if present.

    The closure-based solver can optionally hand its chosen basin to the
    geometry-based cyclic refiner for a final polish. This keeps the cck
    closure logic but uses the proven ccg residual to finish the fit.
    """
    try:
        return _load_module(
            "re_helix_ccg_mod",
            ["re_helix_ccgV3.py"],
        )
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Graph utilities
# ---------------------------------------------------------------------------


def build_helix_components(adjacency):
    """Return connected components of the helix graph."""
    seen = set()  # type: Set[HelixID]
    comps = []  # type: List[Set[HelixID]]
    for h in sorted(adjacency.keys()):
        if h in seen:
            continue
        comp = set()  # type: Set[HelixID]
        stack = [h]
        seen.add(h)
        while stack:
            x = stack.pop()
            comp.add(x)
            for nb in adjacency.get(x, set()):
                if nb not in seen:
                    seen.add(nb)
                    stack.append(nb)
        comps.append(comp)
    return comps


def component_edges(component, adjacency):
    """Return all undirected edges inside a component."""
    edges = set()  # type: Set[frozenset]
    for h in component:
        for nb in adjacency.get(h, set()):
            if nb not in component:
                continue
            if h < nb:
                edges.add(frozenset((h, nb)))
    return edges


def build_tree_for_component(component, adjacency, root):
    """Build a BFS spanning tree for one component."""
    if root not in component:
        raise ValueError("Root helix is not in the component.")

    order = [root]  # type: List[HelixID]
    parent = {}  # type: Dict[HelixID, HelixID]
    seen = set([root])  # type: Set[HelixID]
    queue = [root]

    while queue:
        h = queue.pop(0)
        for nb in sorted(adjacency.get(h, set())):
            if nb not in component or nb in seen:
                continue
            seen.add(nb)
            parent[nb] = h
            order.append(nb)
            queue.append(nb)

    if seen != component:
        missing = component - seen
        raise RuntimeError(
            "BFS did not cover component; missing helices: "
            + ", ".join(helix_id_str(h) for h in sorted(missing))
        )
    return order, parent


# ---------------------------------------------------------------------------
# Small geometry helpers
# ---------------------------------------------------------------------------


def skew_matrix(v):
    return np.array(
        [
            [0.0, -v[2], v[1]],
            [v[2], 0.0, -v[0]],
            [-v[1], v[0], 0.0],
        ],
        dtype=float,
    )



def axis_angle_to_matrix(axis, angle):
    """Rodrigues rotation matrix around unit axis by angle."""
    axis = np.asarray(axis, dtype=float)
    n = float(np.linalg.norm(axis))
    if n < 1.0e-12 or abs(angle) < 1.0e-12:
        return np.eye(3, dtype=float)
    k = axis / n
    K = skew_matrix(k)
    return np.eye(3, dtype=float) + math.sin(angle) * K + (1.0 - math.cos(angle)) * np.dot(K, K)



def compose_transform(R2, t2, R1, t1):
    """Return transform T = T2 o T1."""
    R = np.dot(R2, R1)
    t = np.dot(R2, t1) + t2
    return R, t



def apply_transform_points(R, t, pts):
    return np.dot(pts, R.T) + t



def apply_transform_point(R, t, p):
    return np.dot(R, p) + t



def apply_transform_dir(R, d):
    return np.dot(R, d)



def rotation_about_line(point, direction, angle):
    """Return (R, t) for rotation around a 3D line."""
    point = np.asarray(point, dtype=float)
    direction = np.asarray(direction, dtype=float)
    dlen = float(np.linalg.norm(direction))
    if dlen < 1.0e-12 or abs(angle) < 1.0e-12:
        return np.eye(3, dtype=float), np.zeros(3, dtype=float)
    u = direction / dlen
    R = axis_angle_to_matrix(u, angle)
    t = point - np.dot(R, point)
    return R, t



def line_distance(p1, d1, p2, d2, eps=1.0e-9):
    """Minimum distance between two lines p1+t d1 and p2+s d2."""
    d1u = d1 / max(float(np.linalg.norm(d1)), eps)
    d2u = d2 / max(float(np.linalg.norm(d2)), eps)
    cross_d = np.cross(d1u, d2u)
    norm_cross = float(np.linalg.norm(cross_d))
    w = p2 - p1
    if norm_cross > eps:
        n_unit = cross_d / norm_cross
        return abs(float(np.dot(w, n_unit)))
    proj = float(np.dot(w, d1u))
    w_perp = w - proj * d1u
    return float(np.linalg.norm(w_perp))



def axis_distance(c1, u1, c2, u2, eps=1.0e-9):
    """Minimum distance between two axes (infinite lines)."""
    u1u = u1 / max(float(np.linalg.norm(u1)), eps)
    u2u = u2 / max(float(np.linalg.norm(u2)), eps)
    w0 = c2 - c1
    cross_u = np.cross(u1u, u2u)
    norm_cross = float(np.linalg.norm(cross_u))
    if norm_cross > eps:
        n_unit = cross_u / norm_cross
        return abs(float(np.dot(w0, n_unit)))
    proj = float(np.dot(w0, u1u))
    w_perp = w0 - proj * u1u
    return float(np.linalg.norm(w_perp))



def line_closure_vector(p1, d1, p2, d2, eps=1.0e-9):
    """Vector between closest points of lines p1+t d1 and p2+s d2.

    If the lines intersect, this is zero. For near-parallel lines, returns the
    perpendicular component from p1 to p2 relative to d1.
    """
    d1u = d1 / max(float(np.linalg.norm(d1)), eps)
    d2u = d2 / max(float(np.linalg.norm(d2)), eps)
    w = p1 - p2
    b = float(np.dot(d1u, d2u))
    d = float(np.dot(d1u, w))
    e = float(np.dot(d2u, w))
    denom = 1.0 - b * b
    if abs(denom) > eps:
        t = (b * e - d) / denom
        s = (e - b * d) / denom
        q1 = p1 + t * d1u
        q2 = p2 + s * d2u
        return q1 - q2
    # nearly parallel fallback
    proj = float(np.dot(w, d1u))
    w_perp = w - proj * d1u
    return w_perp




def unpack_exchange_spec_compat(spec):
    """Return (pos1, pos2, kind, rho_deg_or_None) for V2/V3-style specs.

    Older helper versions use (pos1, pos2, kind). Newer versions can
    carry an alignment-only rho angle as (pos1, pos2, kind, rho_deg). This helper
    keeps the cck layer compatible with either form.
    """
    if hasattr(alignmod, "unpack_exchange_spec"):
        return alignmod.unpack_exchange_spec(spec)
    if len(spec) == 3:
        pos1, pos2, kind = spec
        return pos1, pos2, kind, None
    if len(spec) == 4:
        pos1, pos2, third, fourth = spec
        if isinstance(third, str):
            return pos1, pos2, third, fourth
        return pos1, pos2, fourth, third
    raise ValueError("Invalid exchange spec record with %d fields: %r" % (len(spec), spec))


# ---------------------------------------------------------------------------
# Symmetry detection / enforcement helpers
# ---------------------------------------------------------------------------


def _kabsch_fit(P, Q):
    """Return (R, t) mapping row-vector points P -> Q as Q ~= P R^T + t."""
    P = np.asarray(P, dtype=float)
    Q = np.asarray(Q, dtype=float)
    if P.shape != Q.shape or P.ndim != 2 or P.shape[1] != 3 or P.shape[0] == 0:
        raise ValueError("Kabsch fit requires matching Nx3 point arrays.")

    Pc = np.mean(P, axis=0)
    Qc = np.mean(Q, axis=0)
    X = P - Pc
    Y = Q - Qc

    H = np.dot(X.T, Y)
    U, _S, Vt = np.linalg.svd(H)
    R = np.dot(Vt.T, U.T)
    if float(np.linalg.det(R)) < 0.0:
        Vt[-1, :] *= -1.0
        R = np.dot(Vt.T, U.T)
    t = Qc - np.dot(Pc, R.T)
    return R, t


def _rotation_axis_angle(R, eps=1.0e-10):
    """Return (unit_axis, angle in [0, pi]) for a rotation matrix."""
    R = np.asarray(R, dtype=float)
    tr = float(np.trace(R))
    cosang = max(-1.0, min(1.0, 0.5 * (tr - 1.0)))
    angle = math.acos(cosang)
    if abs(angle) < eps:
        return np.array([0.0, 0.0, 1.0], dtype=float), 0.0

    vec = np.array(
        [
            R[2, 1] - R[1, 2],
            R[0, 2] - R[2, 0],
            R[1, 0] - R[0, 1],
        ],
        dtype=float,
    )
    nvec = float(np.linalg.norm(vec))
    if nvec > eps:
        return vec / nvec, angle

    # Near-pi fallback: use the eigenvector with eigenvalue 1.
    vals, vecs = np.linalg.eig(R)
    best_idx = int(np.argmin(np.abs(vals - 1.0)))
    axis = np.real(vecs[:, best_idx])
    n_axis = float(np.linalg.norm(axis))
    if n_axis < eps:
        return np.array([0.0, 0.0, 1.0], dtype=float), angle
    return axis / n_axis, angle


def _cycle_orders_for_component(component, adjacency, start):
    """Return the two possible simple-cycle traversals starting at start."""
    comp = set(component)
    if start not in comp:
        return []

    neigh = {}
    for h in comp:
        nbs = sorted(nb for nb in adjacency.get(h, set()) if nb in comp)
        if len(nbs) != 2:
            return []
        neigh[h] = nbs

    orders = []
    seen = set()
    for first in neigh[start]:
        order = [start]
        prev = start
        curr = first
        ok = True
        while True:
            if curr == start:
                break
            if curr in order:
                ok = False
                break
            order.append(curr)
            next_candidates = [nb for nb in neigh[curr] if nb != prev]
            if len(next_candidates) != 1:
                ok = False
                break
            prev, curr = curr, next_candidates[0]
            if len(order) > len(comp):
                ok = False
                break

        if ok and len(order) == len(comp) and start in neigh[order[-1]]:
            key = tuple(order)
            if key not in seen:
                seen.add(key)
                orders.append(order)
    return orders


def _directed_edge_local_pattern(exchange_specs, chain_to_helix, h_from, h_to):
    """Normalise one directed helix-edge pattern into local chain indices."""
    idx_from = dict((ch, i) for i, ch in enumerate(tuple(sorted(h_from))))
    idx_to = dict((ch, i) for i, ch in enumerate(tuple(sorted(h_to))))
    pattern = []

    for spec in exchange_specs:
        (c1, r1), (c2, r2), kind, _rho_deg = unpack_exchange_spec_compat(spec)
        hh1 = chain_to_helix.get(c1)
        hh2 = chain_to_helix.get(c2)
        if hh1 == h_from and hh2 == h_to:
            pattern.append((idx_from[c1], int(r1), idx_to[c2], int(r2), kind))
        elif hh1 == h_to and hh2 == h_from:
            pattern.append((idx_from[c2], int(r2), idx_to[c1], int(r1), kind))

    pattern.sort()
    return pattern


def detect_symmetric_cycle_component(
    exchange_specs,
    chain_to_helix,
    component,
    adjacency,
    preferred_root=None,
):
    """Detect a simple replicated cyclic pattern such as C2-F21, F2-I21, I2-C21.

    The pattern is recognised when the connected component is a simple cycle and
    every directed helix edge carries the same local exchange pattern after the
    chain IDs are normalised to local indices within each helix group.
    """
    component = set(component)
    if len(component) < 3:
        return None

    edges = component_edges(component, adjacency)
    if len(edges) != len(component):
        return None

    helix_sizes = set(len(tuple(sorted(h))) for h in component)
    if len(helix_sizes) != 1:
        return None

    root = preferred_root if (preferred_root is not None and preferred_root in component) else sorted(component)[0]

    for order in _cycle_orders_for_component(component, adjacency, root):
        base_pattern = _directed_edge_local_pattern(exchange_specs, chain_to_helix, order[0], order[1])
        if not base_pattern:
            continue
        ok = True
        for i in range(len(order)):
            h1 = order[i]
            h2 = order[(i + 1) % len(order)]
            if _directed_edge_local_pattern(exchange_specs, chain_to_helix, h1, h2) != base_pattern:
                ok = False
                break
        if ok:
            return {
                "order": order,
                "edge_pattern": base_pattern,
            }

    return None


def _collect_group_atom_arrays(rec_list, helix_id):
    """Return (keys, atoms, xyz) for one helix, keyed by local-chain index."""
    helix_key = tuple(sorted(helix_id))
    chain_index = dict((ch, i) for i, ch in enumerate(helix_key))
    items = []

    for atom in rec_list:
        if atom.recordName not in ("ATOM", "HETATM"):
            continue
        if atom.chainID not in chain_index:
            continue

        key = (
            chain_index[atom.chainID],
            atom.resSeq,
            atom.string[26],
            atom.name,
            atom.resName,
            atom.string[16],
            atom.recordName,
        )
        items.append((key, atom, np.array([atom.x, atom.y, atom.z], dtype=float)))

    items.sort(key=lambda row: row[0])
    keys = [row[0] for row in items]
    atoms = [row[1] for row in items]
    xyz = np.array([row[2] for row in items], dtype=float)
    return keys, atoms, xyz


def has_detected_symmetric_cycle(
    rec_list,
    exchange_specs,
    explicit_helices=None,
    fix_chain=None,
):
    """Return True if the input specifications indicate a simple symmetric cycle."""
    _, _, chain_to_P_atoms = build_nucleic_acid_maps(rec_list)
    if explicit_helices:
        chain_to_helix = build_chain_to_helix_from_defs(explicit_helices, chain_to_P_atoms)
    else:
        chain_to_helix = compute_chain_partner_map(chain_to_P_atoms)

    fix_helix = None
    if fix_chain is not None and fix_chain in chain_to_helix:
        fix_helix = chain_to_helix[fix_chain]

    _helix_pair_data, adjacency = build_helix_pair_graph(exchange_specs, chain_to_helix)
    if not adjacency:
        return False

    for comp in build_helix_components(adjacency):
        detected = detect_symmetric_cycle_component(
            exchange_specs,
            chain_to_helix,
            comp,
            adjacency,
            preferred_root=fix_helix,
        )
        if detected is not None:
            return True
    return False


def enforce_detected_cyclic_symmetry(
    rec_list,
    exchange_specs,
    explicit_helices=None,
    fix_chain=None,
    max_allowed_p_mismatch=1.0,
):
    """Project simple replicated cyclic components onto exact C_n symmetry.

    This is intentionally conservative: it only triggers for components whose
    exchange specifications themselves indicate a repeated cyclic pattern after
    local chain-index normalisation across replicated helix groups. When it
    triggers, the first/root helix is kept fixed and the remaining helix groups
    are rebuilt by repeated application of one common exact C_n transform.

    The projection is accepted only if the projected symmetric model still keeps
    the P atoms of the repeated terminal pattern close together (within
    max_allowed_p_mismatch, in Å). This avoids forcing an exact symmetry that
    would destroy the intended reciprocal-exchange geometry.
    """
    _, _, chain_to_P_atoms = build_nucleic_acid_maps(rec_list)
    if explicit_helices:
        chain_to_helix = build_chain_to_helix_from_defs(explicit_helices, chain_to_P_atoms)
    else:
        chain_to_helix = compute_chain_partner_map(chain_to_P_atoms)

    fix_helix = None
    if fix_chain is not None and fix_chain in chain_to_helix:
        fix_helix = chain_to_helix[fix_chain]

    _helix_pair_data, adjacency = build_helix_pair_graph(exchange_specs, chain_to_helix)
    if not adjacency:
        return

    components = build_helix_components(adjacency)
    for comp in components:
        detected = detect_symmetric_cycle_component(
            exchange_specs,
            chain_to_helix,
            comp,
            adjacency,
            preferred_root=fix_helix,
        )
        if detected is None:
            continue

        order = detected["order"]
        n_fold = len(order)
        if n_fold < 3:
            continue

        group_keys = []
        atoms_by_helix = {}
        xyz_by_helix = {}
        can_apply = True

        for h in order:
            keys, atoms, xyz = _collect_group_atom_arrays(rec_list, h)
            if xyz.size == 0:
                can_apply = False
                break
            group_keys.append(keys)
            atoms_by_helix[h] = atoms
            xyz_by_helix[h] = xyz

        if not can_apply:
            continue
        base_keys = group_keys[0]
        if any(keys != base_keys for keys in group_keys[1:]):
            sys.stderr.write(
                "[re_helix_cck] Symmetry pattern detected, but atom correspondence differs across replicated helices; skipping exact symmetry projection.\n"
            )
            continue

        # Best common transform from each helix to the next one in the cycle.
        P_blocks = []
        Q_blocks = []
        for i in range(n_fold):
            h1 = order[i]
            h2 = order[(i + 1) % n_fold]
            P_blocks.append(xyz_by_helix[h1])
            Q_blocks.append(xyz_by_helix[h2])
        P_concat = np.vstack(P_blocks)
        Q_concat = np.vstack(Q_blocks)

        R_fit, _t_fit = _kabsch_fit(P_concat, Q_concat)
        axis_dir, _angle_fit = _rotation_axis_angle(R_fit)
        if float(np.linalg.norm(axis_dir)) < 1.0e-10:
            axis_dir = np.array([0.0, 0.0, 1.0], dtype=float)
        axis_dir = axis_dir / float(np.linalg.norm(axis_dir))

        exact_angle = (2.0 * math.pi) / float(n_fold)
        R_exact = axis_angle_to_matrix(axis_dir, exact_angle)

        # Find an axis point p0 so that q ~= (p - p0) R^T + p0.
        mean_b = np.mean(Q_concat - np.dot(P_concat, R_exact.T), axis=0)
        A = np.eye(3, dtype=float) - R_exact
        p0, _resid, _rank, _sing = np.linalg.lstsq(A, mean_b, rcond=None)

        root_h = order[0]
        root_xyz = xyz_by_helix[root_h]
        new_xyz_by_helix = {root_h: root_xyz.copy()}
        for i in range(1, n_fold):
            R_pow = np.linalg.matrix_power(R_exact, i)
            new_xyz_by_helix[order[i]] = np.dot(root_xyz - p0, R_pow.T) + p0

        # Validate that the exact-symmetry projection still respects the repeated
        # terminal P-atom contacts implied by the symmetric exchange pattern.
        p_lookup = {}
        for h in order:
            lookup = {}
            atoms = atoms_by_helix[h]
            xyz_new = new_xyz_by_helix[h]
            for atom, xyz in zip(atoms, xyz_new):
                if atom.name == 'P':
                    local_idx = tuple(sorted(h)).index(atom.chainID)
                    lookup[(local_idx, atom.resSeq)] = xyz
            p_lookup[h] = lookup

        p_mismatches = []
        for i in range(n_fold):
            h1 = order[i]
            h2 = order[(i + 1) % n_fold]
            for idx1, res1, idx2, res2, _kind in detected["edge_pattern"]:
                q1 = p_lookup[h1].get((idx1, res1))
                q2 = p_lookup[h2].get((idx2, res2))
                if q1 is None or q2 is None:
                    continue
                p_mismatches.append(float(np.linalg.norm(q1 - q2)))

        if p_mismatches:
            mean_p_mismatch = float(np.mean(p_mismatches))
            max_p_mismatch = float(np.max(p_mismatches))
            if max_p_mismatch > float(max_allowed_p_mismatch):
                sys.stderr.write(
                    "[re_helix_cck] Symmetry pattern detected, but exact C_%d projection would worsen repeated P-atom contacts (mean=%.3f Å, max=%.3f Å); skipping projection.\n"
                    % (n_fold, mean_p_mismatch, max_p_mismatch)
                )
                continue

        for h in order:
            atoms = atoms_by_helix[h]
            xyz_new = new_xyz_by_helix[h]
            for atom, xyz in zip(atoms, xyz_new):
                atom.update_xyz(float(xyz[0]), float(xyz[1]), float(xyz[2]))

        sys.stderr.write(
            "[re_helix_cck] Applied exact C_%d symmetry projection to cyclic component: %s\n"
            % (n_fold, ", ".join(helix_id_str(h) for h in order))
        )


# ---------------------------------------------------------------------------
# Pairwise tree alignment helpers
# ---------------------------------------------------------------------------


def pairwise_align_component_sequential(
    rec_list,
    component,
    helix_pair_data,
    chain_to_helix,
    axis_dist,
    axis_parallel_flag,
    residue_to_P_atom,
    chain_to_P_atoms,
    order,
    parent,
):
    """Tree-structured pairwise alignment for acyclic components or axis_parallel=y."""
    if not order:
        return

    axis_parallel = axis_parallel_flag

    for moving in order[1:]:
        fixed = parent[moving]
        key = frozenset((fixed, moving))
        entry = helix_pair_data.get(key)
        if entry is None:
            raise RuntimeError(
                "Missing helix_pair_data for edge %s--%s" %
                (helix_id_str(fixed), helix_id_str(moving))
            )

        h1 = entry["helix1"]
        h2 = entry["helix2"]
        all_pairs = entry["pairs"]

        if fixed == h1 and moving == h2:
            pairs_fixed_moving = all_pairs
        elif fixed == h2 and moving == h1:
            pairs_fixed_moving = [(b, a) for (a, b) in all_pairs]
        else:
            raise RuntimeError("Inconsistent helix orientation.")

        chains_fixed_for_axis = set()
        chains_moving_for_axis = set()
        for (c1, _r1), (c2, _r2) in pairs_fixed_moving:
            chains_fixed_for_axis.add(c1)
            chains_moving_for_axis.add(c2)

        axis_dir, center_fixed, center_moving = align_axes_for_pair(
            rec_list,
            fixed,
            moving,
            chain_to_P_atoms,
            axis_dist,
            subset_fixed=chains_fixed_for_axis,
            subset_moving=chains_moving_for_axis,
        )

        _, residue_to_P_atom, chain_to_P_atoms = build_nucleic_acid_maps(rec_list)

        effective_pairs = select_pairs_for_alignment(
            pairs_fixed_moving,
            residue_to_P_atom,
            axis_dir,
            center_fixed,
        )
        if not effective_pairs:
            continue

        objective, p1_list, p2_list, d0, anchor1 = build_pair_objective(
            effective_pairs,
            residue_to_P_atom,
            axis_dir,
            center_fixed,
            center_moving,
            axis_parallel,
        )
        if not p1_list:
            continue

        if axis_parallel:
            x0 = [d0, 0.0, 0.0]
            steps0 = [3.4, math.pi / 2.0, math.pi / 2.0]
            angle_indices = set([1, 2])
        else:
            x0 = [d0, 0.0, 0.0, 0.0]
            steps0 = [3.4, math.pi / 2.0, math.pi / 2.0, math.pi / 2.0]
            angle_indices = set([1, 2, 3])

        best_params, best_val = coordinate_descent(
            objective,
            x0,
            steps0,
            angle_indices=angle_indices,
        )

        if axis_parallel:
            d_opt, theta_opt, phi_opt = best_params
            rho_opt = 0.0
        else:
            d_opt, theta_opt, phi_opt, rho_opt = best_params

        sys.stderr.write(
            "[re_helix_cck]   Pairwise optimised: d = %.3f Å, theta = %.2f°, phi = %.2f°" %
            (d_opt, theta_opt * 180.0 / math.pi, phi_opt * 180.0 / math.pi)
        )
        if not axis_parallel:
            sys.stderr.write(", rho = %.2f°" % (rho_opt * 180.0 / math.pi))
        sys.stderr.write("; sum(dist^2) = %.3f.\n" % best_val)

        apply_transform_to_helix(
            rec_list,
            moving,
            axis_dir,
            center_fixed,
            center_moving,
            d_opt,
            theta_opt,
            phi_opt,
            rho_opt,
            axis_parallel,
            anchor1,
        )

        _, residue_to_P_atom, chain_to_P_atoms = build_nucleic_acid_maps(rec_list)



def compute_rho_axis(axis_dir, axis1_point, axis2_point, anchor1, d, theta, phi):
    """Compute the rho rotation axis used by build_pair_objective/apply_transform_to_helix."""
    u = v_norm(axis_dir)
    C2_phi = rotate_around_line(axis2_point, axis1_point, u, phi)
    C2_phi_d = v_add(C2_phi, v_scale(u, d))

    delta_a2 = v_sub(anchor1, C2_phi_d)
    t = v_dot(delta_a2, u)
    C2_base = v_add(C2_phi_d, v_scale(u, t))

    r_perp = v_sub(C2_base, anchor1)
    r_len = v_length(r_perp)
    if r_len < 1.0e-6:
        if abs(u[0]) < 0.9:
            base_vec = (1.0, 0.0, 0.0)
        else:
            base_vec = (0.0, 1.0, 0.0)
        rho_axis_dir = v_norm(v_cross(u, base_vec))
    else:
        rho_axis_dir = v_scale(r_perp, 1.0 / r_len)

    return anchor1, rho_axis_dir



def pairwise_align_component_base_for_closure(
    rec_list,
    component,
    helix_pair_data,
    chain_to_helix,
    axis_dist,
    order,
    parent,
):
    """Base tree alignment for cyclic components: optimise only d/theta/phi with rho=0.

    Returns per-edge information for the later closure solve.
    """
    edge_info = {}  # type: Dict[Tuple[HelixID, HelixID], Dict[str, object]]

    if not order:
        return edge_info

    _, residue_to_P_atom, chain_to_P_atoms = build_nucleic_acid_maps(rec_list)

    for child in order[1:]:
        parent_h = parent[child]
        key = frozenset((parent_h, child))
        entry = helix_pair_data.get(key)
        if entry is None:
            raise RuntimeError(
                "Missing helix_pair_data for edge %s--%s" %
                (helix_id_str(parent_h), helix_id_str(child))
            )

        h1 = entry["helix1"]
        h2 = entry["helix2"]
        all_pairs = entry["pairs"]

        if parent_h == h1 and child == h2:
            pairs_fixed_moving = all_pairs
        elif parent_h == h2 and child == h1:
            pairs_fixed_moving = [(b, a) for (a, b) in all_pairs]
        else:
            raise RuntimeError("Inconsistent helix orientation.")

        chains_fixed_for_axis = set()
        chains_moving_for_axis = set()
        for (c1, _r1), (c2, _r2) in pairs_fixed_moving:
            chains_fixed_for_axis.add(c1)
            chains_moving_for_axis.add(c2)

        axis_dir, center_fixed, center_moving = align_axes_for_pair(
            rec_list,
            parent_h,
            child,
            chain_to_P_atoms,
            axis_dist,
            subset_fixed=chains_fixed_for_axis,
            subset_moving=chains_moving_for_axis,
        )

        _, residue_to_P_atom, chain_to_P_atoms = build_nucleic_acid_maps(rec_list)

        effective_pairs = select_pairs_for_alignment(
            pairs_fixed_moving,
            residue_to_P_atom,
            axis_dir,
            center_fixed,
        )
        if not effective_pairs:
            continue

        objective4, p1_list, p2_list, d0, anchor1 = build_pair_objective(
            effective_pairs,
            residue_to_P_atom,
            axis_dir,
            center_fixed,
            center_moving,
            axis_parallel=False,
        )
        if not p1_list:
            continue

        def objective3(params3):
            d, theta, phi = params3
            return objective4([d, theta, phi, 0.0])

        x0 = [d0, 0.0, 0.0]
        steps0 = [3.4, math.pi / 2.0, math.pi / 2.0]
        angle_indices = set([1, 2])

        best3, best_val = coordinate_descent(
            objective3,
            x0,
            steps0,
            angle_indices=angle_indices,
        )
        d_opt, theta_opt, phi_opt = best3

        sys.stderr.write(
            "[re_helix_cck]   Base closure edge %s->%s: d = %.3f Å, theta = %.2f°, phi = %.2f°; sum(dist^2) = %.3f.\n"
            % (
                helix_id_str(parent_h),
                helix_id_str(child),
                d_opt,
                theta_opt * 180.0 / math.pi,
                phi_opt * 180.0 / math.pi,
                best_val,
            )
        )

        rho_axis_point, rho_axis_dir = compute_rho_axis(
            axis_dir,
            center_fixed,
            center_moving,
            anchor1,
            d_opt,
            theta_opt,
            phi_opt,
        )

        apply_transform_to_helix(
            rec_list,
            child,
            axis_dir,
            center_fixed,
            center_moving,
            d_opt,
            theta_opt,
            phi_opt,
            0.0,
            False,
            anchor1,
        )

        _, residue_to_P_atom, chain_to_P_atoms = build_nucleic_acid_maps(rec_list)

        edge_info[(parent_h, child)] = {
            "parent": parent_h,
            "child": child,
            "axis_dir": np.array(v_norm(axis_dir), dtype=float),
            "axis1_point": np.array(center_fixed, dtype=float),
            "axis2_point": np.array(center_moving, dtype=float),
            "anchor1": np.array(anchor1, dtype=float),
            "d": float(d_opt),
            "theta": float(theta_opt),
            "phi": float(phi_opt),
            "rho_axis_point": np.array(rho_axis_point, dtype=float),
            "rho_axis_dir": np.array(v_norm(rho_axis_dir), dtype=float),
        }

    return edge_info


# ---------------------------------------------------------------------------
# Solver data and transform evaluation
# ---------------------------------------------------------------------------


def build_component_solver_data(
    rec_list,
    component,
    helix_pair_data,
    chain_to_helix,
    cycle_edges,
    target_rec_list=None,
):
    """Precompute base coordinates, axes, and P-pair data from the current base geometry.

    rec_list is the frozen base geometry used for evaluating transforms.
    target_rec_list, if provided, supplies the atom objects that should be
    updated when the chosen solution is applied. This avoids accidentally
    applying the final transforms only to a detached snapshot.
    """
    atoms_by_helix = {}  # type: Dict[HelixID, List[pdb_atom_record]]
    base_coords_by_helix = {}  # type: Dict[HelixID, np.ndarray]

    if target_rec_list is None:
        target_rec_list = rec_list

    for h in component:
        base_atoms = []  # type: List[pdb_atom_record]
        target_atoms = []  # type: List[pdb_atom_record]
        for atom in rec_list:
            if atom.recordName not in ("ATOM", "HETATM"):
                continue
            if atom.chainID in h:
                base_atoms.append(atom)
        for atom in target_rec_list:
            if atom.recordName not in ("ATOM", "HETATM"):
                continue
            if atom.chainID in h:
                target_atoms.append(atom)
        atoms_by_helix[h] = target_atoms
        base_coords_by_helix[h] = np.array([[a.x, a.y, a.z] for a in base_atoms], dtype=float)

    _, residue_to_P_atom, chain_to_P_atoms = build_nucleic_acid_maps(rec_list)

    base_axis = {}  # type: Dict[HelixID, Tuple[np.ndarray, np.ndarray]]
    for h in component:
        axis_dir, center = compute_helix_axis(chain_to_P_atoms, h)
        base_axis[h] = (
            np.array(v_norm(axis_dir), dtype=float),
            np.array(center, dtype=float),
        )

    p_pairs = []  # type: List[Tuple[HelixID, HelixID, np.ndarray, np.ndarray, bool]]
    axis_pairs = []  # type: List[Tuple[HelixID, HelixID]]
    axis_pair_seen = set()  # type: Set[Tuple[HelixID, HelixID]]
    cycle_pairs = []  # type: List[Tuple[HelixID, HelixID, np.ndarray, np.ndarray]]

    for _key, entry in helix_pair_data.items():
        h1 = entry["helix1"]
        h2 = entry["helix2"]
        if h1 not in component or h2 not in component:
            continue
        if (h1, h2) not in axis_pair_seen:
            axis_pairs.append((h1, h2))
            axis_pair_seen.add((h1, h2))

        is_cycle = frozenset((h1, h2)) in cycle_edges
        for (c1, r1), (c2, r2) in entry["pairs"]:
            a1 = residue_to_P_atom.get((c1, r1))
            a2 = residue_to_P_atom.get((c2, r2))
            if a1 is None or a2 is None:
                continue
            p1 = np.array([a1.x, a1.y, a1.z], dtype=float)
            p2 = np.array([a2.x, a2.y, a2.z], dtype=float)
            p_pairs.append((h1, h2, p1, p2, is_cycle))
            if is_cycle:
                cycle_pairs.append((h1, h2, p1, p2))

    return {
        "atoms_by_helix": atoms_by_helix,
        "base_coords_by_helix": base_coords_by_helix,
        "base_axis": base_axis,
        "p_pairs": p_pairs,
        "cycle_pairs": cycle_pairs,
        "axis_pairs": axis_pairs,
    }



def compute_helix_transforms(
    order,
    parent,
    tree_edges,
    edge_info,
    phi_offsets,
    rho_params,
):
    """Compute final per-helix transforms from base geometry for given phi/rho vars.

    Each tree edge contributes an incremental phi rotation around the parent axis
    and an incremental rho rotation around the rho axis. These are local to the
    current parent/subtree frame and are composed recursively down the tree.
    """
    M = len(tree_edges)
    if len(phi_offsets) != M or len(rho_params) != M:
        raise ValueError("phi_offsets and rho_params must match tree_edges length.")

    edge_to_index = dict((edge, i) for i, edge in enumerate(tree_edges))

    transforms = {}  # type: Dict[HelixID, Tuple[np.ndarray, np.ndarray]]
    if not order:
        return transforms
    root = order[0]
    transforms[root] = (np.eye(3, dtype=float), np.zeros(3, dtype=float))

    for child in order[1:]:
        parent_h = parent[child]
        edge = (parent_h, child)
        idx = edge_to_index[edge]
        phi = float(phi_offsets[idx])
        rho = float(rho_params[idx])
        info = edge_info[edge]

        R_parent, t_parent = transforms[parent_h]

        # Parent-axis line for incremental phi, transformed by ancestor transforms.
        phi_axis_point = apply_transform_point(R_parent, t_parent, info["axis1_point"])  # type: ignore[arg-type]
        phi_axis_dir = apply_transform_dir(R_parent, info["axis_dir"])  # type: ignore[arg-type]
        R_phi, t_phi = rotation_about_line(phi_axis_point, phi_axis_dir, phi)
        R_tmp, t_tmp = compose_transform(R_phi, t_phi, R_parent, t_parent)

        # Rho-axis line, transformed by the state after phi.
        rho_axis_point = apply_transform_point(R_tmp, t_tmp, info["rho_axis_point"])  # type: ignore[arg-type]
        rho_axis_dir = apply_transform_dir(R_tmp, info["rho_axis_dir"])  # type: ignore[arg-type]
        R_rho, t_rho = rotation_about_line(rho_axis_point, rho_axis_dir, rho)
        R_child, t_child = compose_transform(R_rho, t_rho, R_tmp, t_tmp)

        transforms[child] = (R_child, t_child)

    return transforms



def evaluate_geometry_score_from_transforms(
    transforms,
    solver_data,
    component,
    root,
    axis_dist,
    w_pp,
    w_line,
    w_axis,
):
    """Secondary geometry score used only for ranking/filtering roots."""
    base_axis = solver_data["base_axis"]  # type: ignore[assignment]
    p_pairs = solver_data["p_pairs"]  # type: ignore[assignment]
    axis_pairs = solver_data["axis_pairs"]  # type: ignore[assignment]

    total_pp = 0.0
    total_line = 0.0
    total_axis = 0.0

    # transformed axes
    t_axes = {}  # type: Dict[HelixID, Tuple[np.ndarray, np.ndarray]]
    for h in component:
        R, t = transforms[h]
        base_u, base_c = base_axis[h]
        c_t = apply_transform_point(R, t, base_c)
        u_t = apply_transform_dir(R, base_u)
        t_axes[h] = (c_t, u_t)

    # P–P and line terms
    for h1, h2, p1, p2, is_cycle in p_pairs:
        R1, t1 = transforms[h1]
        R2, t2 = transforms[h2]
        q1 = apply_transform_point(R1, t1, p1)
        q2 = apply_transform_point(R2, t2, p2)

        if w_pp != 0.0:
            dq = q1 - q2
            total_pp += float(np.dot(dq, dq))

        if is_cycle and w_line != 0.0:
            c1, u1 = t_axes[h1]
            c2, u2 = t_axes[h2]

            delta1 = q1 - c1
            s1 = float(np.dot(delta1, u1))
            A1 = c1 + s1 * u1
            r1 = q1 - A1
            nr1 = float(np.linalg.norm(r1))
            if nr1 < 1.0e-9:
                continue
            r1u = r1 / nr1

            delta2 = q2 - c2
            s2 = float(np.dot(delta2, u2))
            A2 = c2 + s2 * u2
            r2 = q2 - A2
            nr2 = float(np.linalg.norm(r2))
            if nr2 < 1.0e-9:
                continue
            r2u = r2 / nr2

            d_line = line_distance(q1, r1u, q2, r2u)
            total_line += d_line * d_line

    # axis-distance term
    if w_axis != 0.0:
        for h1, h2 in axis_pairs:
            c1, u1 = t_axes[h1]
            c2, u2 = t_axes[h2]
            d_ax = axis_distance(c1, u1, c2, u2)
            diff = d_ax - axis_dist
            total_axis += diff * diff

    return w_pp * total_pp + w_line * total_line + w_axis * total_axis



def evaluate_closure_residual_from_transforms(
    transforms,
    solver_data,
    component,
    root,
    axis_dist,
):
    """Explicit closure residual used as the primary solve target.

    In the HolT-Hunter spirit, this is an explicit geometric closure vector,
    not the weighted geometry score. For each cycle-edge P pair we use the
    Cartesian closure vector (q1 - q2), and for each cycle edge we add one
    scalar axis-distance closure residual (d_axis - axis_dist).

    This makes the primary solve target analogous to HolT Hunter's ec-target
    closure equation, while the weighted geometry score is used only as a
    secondary ranking/filter.
    """
    base_axis = solver_data["base_axis"]  # type: ignore[assignment]
    cycle_pairs = solver_data["cycle_pairs"]  # type: ignore[assignment]

    # transformed axes
    t_axes = {}  # type: Dict[HelixID, Tuple[np.ndarray, np.ndarray]]
    for h in component:
        R, t = transforms[h]
        base_u, base_c = base_axis[h]
        c_t = apply_transform_point(R, t, base_c)
        u_t = apply_transform_dir(R, base_u)
        t_axes[h] = (c_t, u_t)

    residuals = []  # type: List[float]
    processed_cycle_edges = set()  # type: Set[frozenset]

    for h1, h2, p1, p2 in cycle_pairs:
        R1, t1 = transforms[h1]
        R2, t2 = transforms[h2]
        q1 = apply_transform_point(R1, t1, p1)
        q2 = apply_transform_point(R2, t2, p2)

        # Explicit point-closure residual for the skipped cyclic constraint.
        residuals.extend((q1 - q2).tolist())

        ekey = frozenset((h1, h2))
        if ekey not in processed_cycle_edges:
            processed_cycle_edges.add(ekey)
            c1, u1 = t_axes[h1]
            c2, u2 = t_axes[h2]
            d_ax = axis_distance(c1, u1, c2, u2)
            residuals.append(d_ax - axis_dist)

    if not residuals:
        return np.array([0.0], dtype=float)
    return np.array(residuals, dtype=float)


def apply_component_solution_in_place(
    rec_list,
    solver_data,
    component,
    transforms,
):
    """Apply final per-helix transforms to the real rec_list."""
    atoms_by_helix = solver_data["atoms_by_helix"]  # type: ignore[assignment]
    base_coords_by_helix = solver_data["base_coords_by_helix"]  # type: ignore[assignment]

    for h in component:
        atoms = atoms_by_helix[h]
        base_xyz = base_coords_by_helix[h]
        R, t = transforms[h]
        new_xyz = np.dot(base_xyz, R.T) + t
        for atom, xyz in zip(atoms, new_xyz):
            atom.update_xyz(float(xyz[0]), float(xyz[1]), float(xyz[2]))


# ---------------------------------------------------------------------------
# Post-alignment twist_rod / twist_helix diagnostics
# ---------------------------------------------------------------------------


def _mod360(angle_deg):
    """Return angle in [0, 360)."""
    x = math.fmod(float(angle_deg), 360.0)
    if x < 0.0:
        x += 360.0
    if abs(x) < 1.0e-9 or abs(x - 360.0) < 1.0e-9:
        return 0.0
    return x



def _wrap_signed_180(angle_deg):
    """Return signed angle difference in (-180, 180]."""
    x = _mod360(float(angle_deg) + 180.0) - 180.0
    if x <= -180.0:
        x += 360.0
    return x



def _unit_np(v, eps=1.0e-9):
    """Return a NumPy unit vector, or None for near-zero input."""
    n = float(np.linalg.norm(v))
    if n < eps:
        return None
    return np.asarray(v, dtype=float) / n



def _project_to_axis(point, axis_center, axis_dir):
    """Project a point onto an infinite helix-axis line."""
    p = np.asarray(point, dtype=float)
    c = np.asarray(axis_center, dtype=float)
    u = np.asarray(axis_dir, dtype=float)
    return c + float(np.dot(p - c, u)) * u



def _perp_to_axis(vec, axis_dir):
    """Remove the component of vec parallel to axis_dir."""
    v = np.asarray(vec, dtype=float)
    u = np.asarray(axis_dir, dtype=float)
    return v - float(np.dot(v, u)) * u



def _radial_unit_from_axis(point, axis_center, axis_dir, eps=1.0e-9):
    """Unit vector from a helix axis to a point, projected perpendicular to axis."""
    p = np.asarray(point, dtype=float)
    a = _project_to_axis(p, axis_center, axis_dir)
    return _unit_np(p - a, eps=eps)



def _signed_angle_about_axis_deg(v_from, v_to, axis_dir, eps=1.0e-9):
    """Signed angle from v_from to v_to around axis_dir, in degrees."""
    a = _unit_np(v_from, eps=eps)
    b = _unit_np(v_to, eps=eps)
    u = _unit_np(axis_dir, eps=eps)
    if a is None or b is None or u is None:
        return None
    sin_term = float(np.dot(u, np.cross(a, b)))
    cos_term = float(np.dot(a, b))
    return math.degrees(math.atan2(sin_term, cos_term))



def _res_label(chain_id, res_seq):
    return "%d%s" % (int(res_seq), str(chain_id))



def _choose_reference_chain_for_helix(chain_to_P_atoms, helix_id):
    """Choose the chain with the most P atoms as this helix's register reference."""
    best_chain = None
    best_count = -1
    for ch in sorted(helix_id):
        atoms = [a for a in chain_to_P_atoms.get(ch, []) if getattr(a, "recordName", "ATOM") in ("ATOM", "HETATM")]
        if len(atoms) > best_count:
            best_chain = ch
            best_count = len(atoms)
    return best_chain



def _chain_p_atoms_sorted_by_resseq(chain_to_P_atoms, chain_id):
    return sorted(
        [a for a in chain_to_P_atoms.get(chain_id, []) if getattr(a, "recordName", "ATOM") in ("ATOM", "HETATM")],
        key=lambda a: (int(a.resSeq), int(a.serial)),
    )



def _chain_p_atoms_sorted_by_axis(chain_to_P_atoms, chain_id, axis_center, axis_dir):
    """Return P atoms on one chain sorted by coordinate along a helix axis."""
    c = np.asarray(axis_center, dtype=float)
    u = np.asarray(axis_dir, dtype=float)
    atoms = [a for a in chain_to_P_atoms.get(chain_id, []) if getattr(a, "recordName", "ATOM") in ("ATOM", "HETATM")]

    def key(atom):
        p = np.array([atom.x, atom.y, atom.z], dtype=float)
        s = float(np.dot(p - c, u))
        return (s, int(atom.resSeq), int(atom.serial))

    return sorted(atoms, key=key)



def _orient_axis_and_calibrate_twist(chain_to_P_atoms, helix_id, axis_dir, axis_center, helical_repeat):
    """Orient one helix axis and measure its native phosphate twist direction.

    The PCA axis sign is arbitrary. We orient +axis to follow increasing residue
    numbers on a reference chain when possible. We then measure the observed
    signed phosphate rotation between neighbouring P atoms on that reference
    chain. The sign of that observed rotation defines the native-positive twist
    convention used by the twist report.

    ``helical_repeat`` is used only as a fallback when the observed phosphate
    step cannot be measured from the input coordinates.
    """
    u = np.array(v_norm(tuple(axis_dir)), dtype=float)
    c = np.asarray(axis_center, dtype=float)
    ref_chain = _choose_reference_chain_for_helix(chain_to_P_atoms, helix_id)

    if ref_chain is not None:
        ref_atoms = _chain_p_atoms_sorted_by_resseq(chain_to_P_atoms, ref_chain)
        if len(ref_atoms) >= 2:
            p_first = np.array([ref_atoms[0].x, ref_atoms[0].y, ref_atoms[0].z], dtype=float)
            p_last = np.array([ref_atoms[-1].x, ref_atoms[-1].y, ref_atoms[-1].z], dtype=float)
            if float(np.dot(u, p_last - p_first)) < 0.0:
                u = -u

    rises = []
    twist_steps = []
    if ref_chain is not None:
        ref_atoms = _chain_p_atoms_sorted_by_resseq(chain_to_P_atoms, ref_chain)
        for a1, a2 in zip(ref_atoms[:-1], ref_atoms[1:]):
            delta_res = abs(int(a2.resSeq) - int(a1.resSeq))
            if delta_res <= 0:
                continue
            p1 = np.array([a1.x, a1.y, a1.z], dtype=float)
            p2 = np.array([a2.x, a2.y, a2.z], dtype=float)
            ds = abs(float(np.dot(p2 - p1, u)))
            if ds > 1.0e-6:
                rises.append(ds / float(delta_res))
            r1 = _radial_unit_from_axis(p1, c, u)
            r2 = _radial_unit_from_axis(p2, c, u)
            if r1 is not None and r2 is not None:
                ang = _signed_angle_about_axis_deg(r1, r2, u)
                if ang is not None and abs(ang) > 1.0e-6:
                    twist_steps.append(float(ang) / float(delta_res))

    if rises:
        rise_per_bp = float(np.median(np.array(rises, dtype=float)))
    else:
        rise_per_bp = 3.4

    if twist_steps:
        median_step = float(np.median(np.array(twist_steps, dtype=float)))
        twist_sign = 1.0 if median_step >= 0.0 else -1.0
        native_step_abs = abs(median_step)
    else:
        median_step = 360.0 / float(helical_repeat)
        twist_sign = 1.0
        native_step_abs = abs(median_step)

    return u, c, ref_chain, rise_per_bp, twist_sign, median_step, native_step_abs



def _axis_coperpendicular_unit_to_neighbor(
    own_axis_center,
    own_axis_dir,
    neighbor_axis_center,
    neighbor_axis_dir,
    eps=1.0e-9,
):
    """Return the co-perpendicular connector from one helix axis to a neighbour.

    For non-parallel axes, this is the shortest connector between the two
    infinite axis lines, which is perpendicular to both axes. For nearly
    parallel axes, the shortest connector is the center-to-center vector after
    removing the component along the helix axis. The returned direction always
    points from this helix axis toward the neighbouring helix axis.
    """
    c1 = np.asarray(own_axis_center, dtype=float)
    u1 = _unit_np(own_axis_dir, eps=eps)
    c2 = np.asarray(neighbor_axis_center, dtype=float)
    u2 = _unit_np(neighbor_axis_dir, eps=eps)
    if u1 is None or u2 is None:
        return None, "unavailable"

    w0 = c1 - c2
    b = float(np.dot(u1, u2))
    d = float(np.dot(u1, w0))
    e = float(np.dot(u2, w0))
    denom = 1.0 - b * b

    if abs(denom) > eps:
        s = (b * e - d) / denom
        t = (e - b * d) / denom
        p_own = c1 + s * u1
        p_nei = c2 + t * u2
        vec = p_nei - p_own
        source = "axis_coperpendicular_connector"
    else:
        vec = c2 - c1
        source = "parallel_axis_coperpendicular_connector"

    vec = _perp_to_axis(vec, u1)
    unit = _unit_np(vec, eps=eps)

    if unit is None:
        return None, "unavailable"

    return unit, source



def _estimate_bp_span_for_twist(site_a, site_b, axial_span, rise_per_bp):
    """Estimate base-step span between two exchange sites on the same helix."""
    if site_a["chain"] == site_b["chain"]:
        delta_res = abs(int(site_b["resSeq"]) - int(site_a["resSeq"]))
        if delta_res > 0:
            return float(delta_res), "resSeq_same_chain"
    if rise_per_bp > 1.0e-9:
        return abs(float(axial_span)) / float(rise_per_bp), "axis_span_over_observed_rise"
    return float("nan"), "unavailable"



def _choose_unwrapped_total_from_mod(mod_angle_deg, target_total_deg):
    """Choose mod_angle + 360*k closest to a positive target total twist."""
    mod_angle = _mod360(mod_angle_deg)
    if not math.isfinite(target_total_deg) or target_total_deg <= 1.0e-9:
        return mod_angle

    k0 = int(round((float(target_total_deg) - mod_angle) / 360.0))
    candidates = []
    for k in range(max(0, k0 - 3), max(0, k0 + 4)):
        val = mod_angle + 360.0 * float(k)
        if val > 1.0e-9:
            candidates.append(val)
    if not candidates:
        candidates.append(360.0 if mod_angle < 1.0e-9 else mod_angle)
    return min(candidates, key=lambda v: abs(v - float(target_total_deg)))



def _build_chain_twist_path(chain_to_P_atoms, chain_id, axis_center, axis_dir, native_sign):
    """Build cumulative observed P-atom twist for one chain along the axis.

    The cumulative values are derived only from observed phosphate rotations
    around the helix axis; no ideal helical-repeat model is used here.
    """
    atoms = _chain_p_atoms_sorted_by_axis(chain_to_P_atoms, chain_id, axis_center, axis_dir)
    entries = []
    c = np.asarray(axis_center, dtype=float)
    u = np.asarray(axis_dir, dtype=float)
    for atom in atoms:
        p = np.array([atom.x, atom.y, atom.z], dtype=float)
        r = _radial_unit_from_axis(p, c, u)
        if r is None:
            continue
        s = float(np.dot(p - c, u))
        entries.append({"atom": atom, "resSeq": int(atom.resSeq), "point": p, "radial": r, "s": s, "cum": 0.0})

    if len(entries) < 2:
        return entries, {}

    cum = 0.0
    entries[0]["cum"] = cum
    for idx in range(len(entries) - 1):
        raw = _signed_angle_about_axis_deg(entries[idx]["radial"], entries[idx + 1]["radial"], u)
        if raw is None:
            inc = 0.0
        else:
            inc = float(native_sign) * float(raw)
            # Adjacent phosphate steps should be the small helical step, not the
            # complementary 360-step. Keep each increment in a near-principal range.
            while inc <= -180.0:
                inc += 360.0
            while inc > 180.0:
                inc -= 360.0
        cum += inc
        entries[idx + 1]["cum"] = cum

    by_res = {}
    for item in entries:
        by_res[int(item["resSeq"])] = item
    return entries, by_res



def _estimate_total_twist_helix_deg(
    site_a,
    site_b,
    axis_rec,
    chain_to_P_atoms,
    twist_helix_mod,
    bp_span,
    helical_repeat,
):
    """Estimate total observed helical twist for strain denominator.

    The modulo twist_helix is measured directly from the two RE-site phosphate
    radial vectors. For the denominator, we need the total twist, including full
    turns. When both RE sites are on one chain, we sum observed inter-phosphate
    rotations along that chain. Otherwise, we infer the full-turn count by using
    the observed local phosphate twist per base pair, with helical_repeat only
    as a final fallback.
    """
    native_sign = float(axis_rec["twist_sign"])
    u = axis_rec["u"]
    c = axis_rec["c"]

    if site_a["chain"] == site_b["chain"]:
        chain_id = site_a["chain"]
        cache = axis_rec.setdefault("chain_twist_cache", {})
        if chain_id not in cache:
            cache[chain_id] = _build_chain_twist_path(chain_to_P_atoms, chain_id, c, u, native_sign)
        _entries, by_res = cache[chain_id]
        item_a = by_res.get(int(site_a["resSeq"]))
        item_b = by_res.get(int(site_b["resSeq"]))
        if item_a is not None and item_b is not None:
            observed_total = abs(float(item_b["cum"]) - float(item_a["cum"]))
            if observed_total > 1.0e-9:
                total = _choose_unwrapped_total_from_mod(twist_helix_mod, observed_total)
                return total, "same_chain_cumulative_observed_P_rotation"

    native_step_abs = abs(float(axis_rec.get("native_step_abs", float("nan"))))
    if math.isfinite(bp_span) and bp_span > 1.0e-9 and native_step_abs > 1.0e-9:
        target_total = float(bp_span) * native_step_abs
        total = _choose_unwrapped_total_from_mod(twist_helix_mod, target_total)
        return total, "endpoint_P_rotation_mod_plus_axis_span_observed_step"

    if math.isfinite(bp_span) and bp_span > 1.0e-9 and helical_repeat > 0.0:
        target_total = float(bp_span) * 360.0 / float(helical_repeat)
        total = _choose_unwrapped_total_from_mod(twist_helix_mod, target_total)
        return total, "endpoint_P_rotation_mod_plus_helical_repeat_fallback"

    total = twist_helix_mod if twist_helix_mod > 1.0e-9 else float("nan")
    return total, "endpoint_P_rotation_mod_only"



def compute_twist_rod_helix_rows(
    rec_list,
    exchange_specs,
    explicit_helices=None,
    helical_repeat=10.5,
    pairing="consecutive",
):
    """Compute post-alignment twist_rod/twist_helix diagnostics.

    For each helix with two or more reciprocal-exchange P sites, this reports
    segments between exchange sites.

    * twist_rod is the PDB-native rod-required phase: the rotation, around the
      helix axis, between the co-perpendicular connector vectors from this
      helix axis to the neighbouring helix axes at the two exchange sites.
    * twist_helix is measured directly from the rotation of the two
      reciprocal-exchange-site phosphate radial vectors around the same helix
      axis.

    Both are reported modulo 360. The signed difference is
    wrap180(twist_rod_mod - twist_helix_mod), matching the HolT-Hunter strain
    convention that asks how much the helix must twist to fit the rod geometry:
    positive means overtwist is needed, negative means undertwist. Strain is
    difference_deg / total_twist_helix_deg.
    """
    notes = []  # type: List[str]
    rows = []  # type: List[Dict[str, object]]

    if helical_repeat <= 0.0:
        raise ValueError("helical_repeat must be positive.")

    _, residue_to_P_atom, chain_to_P_atoms = build_nucleic_acid_maps(rec_list)
    if explicit_helices:
        chain_to_helix = build_chain_to_helix_from_defs(explicit_helices, chain_to_P_atoms)
    else:
        chain_to_helix = compute_chain_partner_map(chain_to_P_atoms)

    helix_pair_data, adjacency = build_helix_pair_graph(exchange_specs, chain_to_helix)
    if not helix_pair_data:
        notes.append("No inter-helix reciprocal-exchange pairs were found; twist report is empty.")
        return rows, notes

    components = build_helix_components(adjacency)
    component_id_by_helix = {}
    for comp_idx, comp in enumerate(components, start=1):
        for h in comp:
            component_id_by_helix[h] = comp_idx

    axes = {}
    for h in sorted(component_id_by_helix.keys()):
        axis_dir, center = compute_helix_axis(chain_to_P_atoms, h)
        u, c, ref_chain, rise_per_bp, twist_sign, median_step, native_step_abs = _orient_axis_and_calibrate_twist(
            chain_to_P_atoms,
            h,
            axis_dir,
            center,
            helical_repeat,
        )
        axes[h] = {
            "u": u,
            "c": c,
            "ref_chain": ref_chain,
            "rise_per_bp": rise_per_bp,
            "twist_sign": twist_sign,
            "median_step": median_step,
            "native_step_abs": native_step_abs,
            "chain_twist_cache": {},
        }

    sites_by_helix = {}  # type: Dict[HelixID, List[Dict[str, object]]]

    for spec_idx, spec in enumerate(exchange_specs, start=1):
        (c1, r1), (c2, r2), kind, _rho_deg = unpack_exchange_spec_compat(spec)
        h1 = chain_to_helix.get(c1)
        h2 = chain_to_helix.get(c2)
        if h1 is None or h2 is None:
            notes.append("Skipping spec %d because one chain is not assigned to a helix." % spec_idx)
            continue
        if h1 == h2:
            continue

        endpoints = [
            (h1, h2, c1, r1, c2, r2),
            (h2, h1, c2, r2, c1, r1),
        ]
        for h, neighbor_h, c_self, r_self, c_partner, r_partner in endpoints:
            atom = residue_to_P_atom.get((c_self, r_self))
            partner_atom = residue_to_P_atom.get((c_partner, r_partner))
            if atom is None:
                notes.append("Skipping %s in spec %d because its P atom was not found." % (_res_label(c_self, r_self), spec_idx))
                continue
            if partner_atom is None:
                notes.append("Skipping partner %s in spec %d because its P atom was not found." % (_res_label(c_partner, r_partner), spec_idx))
                continue
            axis_rec = axes.get(h)
            neighbor_axis_rec = axes.get(neighbor_h)
            if axis_rec is None or neighbor_axis_rec is None:
                continue
            p = np.array([atom.x, atom.y, atom.z], dtype=float)
            p_partner = np.array([partner_atom.x, partner_atom.y, partner_atom.z], dtype=float)
            u = axis_rec["u"]
            c = axis_rec["c"]
            self_radial = _radial_unit_from_axis(p, c, u)
            if self_radial is None:
                notes.append("Skipping %s in spec %d because its phosphate radial vector is near zero." % (_res_label(c_self, r_self), spec_idx))
                continue

            rod_radial, rod_source = _axis_coperpendicular_unit_to_neighbor(
                c,
                u,
                neighbor_axis_rec["c"],
                neighbor_axis_rec["u"],
            )
            if rod_radial is None:
                notes.append("Skipping %s in spec %d because no axis co-perpendicular connector vector could be defined." % (_res_label(c_self, r_self), spec_idx))
                continue

            s_coord = float(np.dot(p - c, u))
            site = {
                "component": component_id_by_helix.get(h, 0),
                "helix": h,
                "neighbor": neighbor_h,
                "chain": c_self,
                "resSeq": int(r_self),
                "partner_chain": c_partner,
                "partner_resSeq": int(r_partner),
                "kind": str(kind),
                "spec_index": int(spec_idx),
                "point": p,
                "partner_point": p_partner,
                "s": s_coord,
                "self_radial": self_radial,
                "rod_radial": rod_radial,
                "rod_vector_source": rod_source,
                "site_label": _res_label(c_self, r_self),
                "partner_label": _res_label(c_partner, r_partner),
            }
            sites_by_helix.setdefault(h, []).append(site)

    for h in sorted(sites_by_helix.keys()):
        sites = sites_by_helix[h]
        sites.sort(key=lambda x: (float(x["s"]), str(x["chain"]), int(x["resSeq"]), int(x["spec_index"])))
        if len(sites) < 2:
            continue

        if pairing == "all":
            pair_iter = []
            for i in range(len(sites) - 1):
                for j in range(i + 1, len(sites)):
                    pair_iter.append((sites[i], sites[j]))
        else:
            pair_iter = list(zip(sites[:-1], sites[1:]))

        axis_rec = axes[h]
        u = axis_rec["u"]
        rise_per_bp = float(axis_rec["rise_per_bp"])
        twist_sign = float(axis_rec["twist_sign"])
        ref_chain = axis_rec["ref_chain"]
        median_step = float(axis_rec["median_step"])
        native_step_abs = float(axis_rec["native_step_abs"])

        for site_a, site_b in pair_iter:
            axial_span = float(site_b["s"]) - float(site_a["s"])
            if abs(axial_span) < 1.0e-9 and site_a["chain"] == site_b["chain"] and site_a["resSeq"] == site_b["resSeq"]:
                continue

            bp_span, bp_span_source = _estimate_bp_span_for_twist(site_a, site_b, axial_span, rise_per_bp)

            rod_axis_angle = _signed_angle_about_axis_deg(site_a["rod_radial"], site_b["rod_radial"], u)
            helix_axis_angle = _signed_angle_about_axis_deg(site_a["self_radial"], site_b["self_radial"], u)
            if rod_axis_angle is None or helix_axis_angle is None:
                continue

            # Convert both angles into the observed native-positive helix twist
            # convention before taking modulo 360.
            twist_rod_mod = _mod360(twist_sign * float(rod_axis_angle))
            twist_helix_endpoint_mod = _mod360(twist_sign * float(helix_axis_angle))

            total_twist_helix, total_twist_source = _estimate_total_twist_helix_deg(
                site_a,
                site_b,
                axis_rec,
                chain_to_P_atoms,
                twist_helix_endpoint_mod,
                bp_span,
                helical_repeat,
            )
            twist_helix_mod = _mod360(twist_helix_endpoint_mod)

            difference = _wrap_signed_180(twist_rod_mod - twist_helix_mod)
            strain = difference / total_twist_helix if math.isfinite(total_twist_helix) and abs(total_twist_helix) > 1.0e-12 else float("nan")

            if difference > 1.0e-6:
                twist_call = "overtwist"
            elif difference < -1.0e-6:
                twist_call = "undertwist"
            else:
                twist_call = "matched"

            rows.append(
                {
                    "component": int(component_id_by_helix.get(h, 0)),
                    "helix": helix_id_str(h),
                    "reference_chain": ref_chain if ref_chain is not None else "",
                    "site_1": site_a["site_label"],
                    "site_2": site_b["site_label"],
                    "neighbor_1": helix_id_str(site_a["neighbor"]),
                    "neighbor_2": helix_id_str(site_b["neighbor"]),
                    "partner_1": site_a["partner_label"],
                    "partner_2": site_b["partner_label"],
                    "kind_1": site_a["kind"],
                    "kind_2": site_b["kind"],
                    "spec_index_1": int(site_a["spec_index"]),
                    "spec_index_2": int(site_b["spec_index"]),
                    "axis_s_1_A": float(site_a["s"]),
                    "axis_s_2_A": float(site_b["s"]),
                    "axis_span_A": abs(float(axial_span)),
                    "estimated_rise_A_per_bp": rise_per_bp,
                    "bp_span": float(bp_span) if math.isfinite(bp_span) else float("nan"),
                    "bp_span_source": bp_span_source,
                    "fallback_helical_repeat_bp_per_turn": float(helical_repeat),
                    "observed_native_step_deg_per_bp": median_step,
                    "observed_native_step_abs_deg_per_bp": native_step_abs,
                    "twist_rod_mod_deg": twist_rod_mod,
                    "twist_helix_mod_deg": twist_helix_mod,
                    "total_twist_helix_deg": total_twist_helix,
                    "total_twist_helix_source": total_twist_source,
                    "difference_deg": difference,
                    "twist_call": twist_call,
                    "strain": strain,
                    "raw_rod_axis_angle_deg": float(rod_axis_angle),
                    "raw_helix_axis_angle_deg": float(helix_axis_angle),
                    "native_twist_sign": twist_sign,
                    "rod_vector_source_1": site_a["rod_vector_source"],
                    "rod_vector_source_2": site_b["rod_vector_source"],
                }
            )

    rows.sort(key=lambda row: (int(row["component"]), str(row["helix"]), float(row["axis_s_1_A"]), float(row["axis_s_2_A"])))
    if not rows:
        notes.append("No helix had at least two usable reciprocal-exchange P sites for twist reporting.")
    return rows, notes



def _format_tsv_value(value):
    if isinstance(value, float):
        if math.isnan(value):
            return "nan"
        if math.isinf(value):
            return "inf" if value > 0.0 else "-inf"
        return "%.6f" % value
    if isinstance(value, int):
        return str(value)
    if value is None:
        return ""
    return str(value).replace("\t", " ")



def write_twist_report_tsv(path, rows, notes, pairing):
    """Write twist diagnostics to a tab-separated text file."""
    columns = [
        "component",
        "helix",
        "reference_chain",
        "site_1",
        "site_2",
        "neighbor_1",
        "neighbor_2",
        "partner_1",
        "partner_2",
        "kind_1",
        "kind_2",
        "spec_index_1",
        "spec_index_2",
        "axis_s_1_A",
        "axis_s_2_A",
        "axis_span_A",
        "estimated_rise_A_per_bp",
        "bp_span",
        "bp_span_source",
        "fallback_helical_repeat_bp_per_turn",
        "observed_native_step_deg_per_bp",
        "observed_native_step_abs_deg_per_bp",
        "twist_rod_mod_deg",
        "twist_helix_mod_deg",
        "total_twist_helix_deg",
        "total_twist_helix_source",
        "difference_deg",
        "twist_call",
        "strain",
        "raw_rod_axis_angle_deg",
        "raw_helix_axis_angle_deg",
        "native_twist_sign",
        "rod_vector_source_1",
        "rod_vector_source_2",
    ]
    with open(path, "w") as fout:
        fout.write("# re_helix_cckV3 twist diagnostics\n")
        fout.write("# pairing=%s\n" % pairing)
        fout.write("# twist_rod_mod_deg = rod-required phase between axis-to-axis co-perpendicular connector vectors around the helix axis, modulo 360.\n")
        fout.write("# twist_helix_mod_deg = direct phase between the two reciprocal-exchange-site P-atom radial vectors around the helix axis, modulo 360.\n")
        fout.write("# difference_deg = signed wrap180(twist_rod_mod_deg - twist_helix_mod_deg); positive=overtwist needed, negative=undertwist needed.\n")
        fout.write("# strain = difference_deg / total_twist_helix_deg. total_twist_helix_deg includes full turns when they can be inferred from observed P-atom rotations.\n")
        fout.write("# fallback_helical_repeat_bp_per_turn is used only when observed phosphate-step information is unavailable.\n")
        for note in notes:
            fout.write("# NOTE: %s\n" % str(note).replace("\n", " "))
        fout.write("\t".join(columns) + "\n")
        for row in rows:
            fout.write("\t".join(_format_tsv_value(row.get(col, "")) for col in columns) + "\n")



# ---------------------------------------------------------------------------
# HolT-Hunter-like multi-start closure solver
# ---------------------------------------------------------------------------


def make_closure_seed_bank(num_tree_edges, max_root_attempts):
    """Deterministic multi-start seeds for (phi_offsets, rho) variables.

    This is intentionally not an exhaustive angle screen. It is a small,
    deterministic bank of root starts, in the spirit of HolT Hunter's multiple
    FindRoot seeds.
    """
    M = int(num_tree_edges)
    if M <= 0:
        return [np.zeros(0, dtype=float)]

    seeds = []  # type: List[np.ndarray]

    # sign patterns for rho across edges (limited but deterministic)
    patterns = []  # type: List[List[float]]
    max_patterns = min(1 << M, 8)
    for k in range(max_patterns):
        pattern = []
        for i in range(M):
            sign = 1.0 if ((k >> i) & 1) == 0 else -1.0
            pattern.append(sign)
        patterns.append(pattern)

    rho_bank = [
        0.0,
        math.pi / 2.0,
        -math.pi / 2.0,
        math.pi,
        -math.pi,
        2.0 * math.pi / 3.0,
        -2.0 * math.pi / 3.0,
        math.pi / 3.0,
        -math.pi / 3.0,
    ]
    phi_bank = [
        0.0,
        math.pi / 6.0,
        -math.pi / 6.0,
        math.pi / 4.0,
        -math.pi / 4.0,
        math.pi / 2.0,
        -math.pi / 2.0,
    ]

    # zero seed first
    seeds.append(np.zeros(2 * M, dtype=float))

    # rho-biased seeds
    for amp in rho_bank[1:]:
        for pat in patterns:
            phi = [0.0] * M
            rho = [pat[i] * amp for i in range(M)]
            seeds.append(np.array(phi + rho, dtype=float))
            if len(seeds) >= max_root_attempts:
                return seeds[:max_root_attempts]

    # mixed phi/rho seeds
    for amp_phi in phi_bank[1:]:
        for amp_rho in (math.pi / 2.0, -math.pi / 2.0):
            for pat in patterns:
                phi = [pat[i] * amp_phi for i in range(M)]
                rho = [pat[(i + 1) % len(pat)] * amp_rho for i in range(M)]
                seeds.append(np.array(phi + rho, dtype=float))
                if len(seeds) >= max_root_attempts:
                    return seeds[:max_root_attempts]

    return seeds[:max_root_attempts]



def refine_cyclic_components_closure(
    rec_list,
    exchange_specs,
    axis_dist,
    axis_parallel_flag,
    explicit_helices=None,
    fix_chain=None,
    w_pp=1.0,
    w_line=1.0,
    w_axis=1.0e4,
    max_root_attempts=20,
    root_maxfev=800,
    closure_rel_tol=0.10,
    closure_abs_tol=1.0e-3,
):
    """Refine cyclic components by solving explicit closure residuals.

    This is the HolT-Hunter-like stage:
      - primary solve target: closure residual on cycle edges,
      - secondary ranking/filter: geometry score.

    For square systems (len(residual)==len(params)) we use scipy.optimize.root.
    Otherwise we use least_squares on the same closure residual.
    """
    # Current helix assignment from the current structure.
    _, _, chain_to_P_atoms = build_nucleic_acid_maps(rec_list)
    if explicit_helices:
        chain_to_helix = build_chain_to_helix_from_defs(explicit_helices, chain_to_P_atoms)
    else:
        chain_to_helix = compute_chain_partner_map(chain_to_P_atoms)

    fix_helix = None  # type: Optional[HelixID]
    if fix_chain is not None and fix_chain in chain_to_helix:
        fix_helix = chain_to_helix[fix_chain]

    helix_pair_data, adjacency = build_helix_pair_graph(exchange_specs, chain_to_helix)
    if not helix_pair_data:
        return

    components = build_helix_components(adjacency)

    for comp in components:
        edges = component_edges(comp, adjacency)
        if len(edges) <= max(0, len(comp) - 1):
            continue  # acyclic

        root_h = fix_helix if (fix_helix is not None and fix_helix in comp) else sorted(comp)[0]
        order, parent = build_tree_for_component(comp, adjacency, root_h)
        tree_edge_keys = set(frozenset((parent[ch], ch)) for ch in order[1:])
        cycle_edges = set(e for e in edges if e not in tree_edge_keys)

        def _edge_label(edge):
            h1, h2 = sorted(list(edge))
            return "%s--%s" % (helix_id_str(h1), helix_id_str(h2))

        sys.stderr.write(
            "[re_helix_cck] Closure refinement for cyclic component: "
            + ", ".join(helix_id_str(h) for h in sorted(comp))
            + "; cycle edges="
            + (", ".join(sorted(_edge_label(edge) for edge in cycle_edges)) if cycle_edges else "(none)")
            + "\n"
        )

        # Base DTP alignment with rho=0
        edge_info = pairwise_align_component_base_for_closure(
            rec_list,
            comp,
            helix_pair_data,
            chain_to_helix,
            axis_dist,
            order,
            parent,
        )

        tree_edges = []  # type: List[Tuple[HelixID, HelixID]]
        for child_h in order[1:]:
            edge = (parent[child_h], child_h)
            if edge in edge_info:
                tree_edges.append(edge)

        if not tree_edges:
            sys.stderr.write(
                "[re_helix_cck]   Warning: no valid tree edges for closure solve; keeping base geometry.\n"
            )
            continue

        # Base snapshot after DTP alignment.
        rec_base = [type(a)(a.string) for a in rec_list]

        # Solver data on the base geometry.
        solver_data = build_component_solver_data(
            rec_base,
            comp,
            helix_pair_data,
            chain_to_helix,
            cycle_edges,
            target_rec_list=rec_list,
        )

        M = len(tree_edges)
        n_params = 2 * M

        def transforms_from_params(params):
            phi_vec = params[:M]
            rho_vec = params[M:]
            return compute_helix_transforms(order, parent, tree_edges, edge_info, phi_vec, rho_vec)

        def closure_residual(params):
            tr = transforms_from_params(params)
            return evaluate_closure_residual_from_transforms(
                tr,
                solver_data,
                comp,
                root_h,
                axis_dist,
            )

        def geometry_score(params):
            tr = transforms_from_params(params)
            return evaluate_geometry_score_from_transforms(
                tr,
                solver_data,
                comp,
                root_h,
                axis_dist,
                w_pp,
                w_line,
                w_axis,
            )

        # Root-like deterministic seed bank.
        seeds = make_closure_seed_bank(M, max_root_attempts)
        candidates = []  # type: List[Tuple[float, float, np.ndarray, bool, int, str]]

        # decide whether the explicit closure residual is square
        res0 = closure_residual(seeds[0])
        n_resid = int(len(res0))
        use_root = (n_resid == n_params and n_params > 0)
        sys.stderr.write(
            "[re_helix_cck]   closure system: vars=%d, residuals=%d, solver=%s\n"
            % (n_params, n_resid, "root" if use_root else "least_squares")
        )

        for attempt_idx, x0 in enumerate(seeds, start=1):
            try:
                if use_root:
                    sol = root(
                        closure_residual,
                        x0,
                        method="hybr",
                        options={"maxfev": int(root_maxfev)},
                    )
                    x = np.array(sol.x, dtype=float)
                    success = bool(sol.success)
                    msg = getattr(sol, "message", "")
                    nfev = int(getattr(sol, "nfev", -1))
                    # HolT-Hunter-like root first, then a closure-only polish.
                    try:
                        sol_polish = least_squares(
                            closure_residual,
                            x,
                            method="trf",
                            max_nfev=max(50, int(root_maxfev // 2)),
                            xtol=1.0e-10,
                            ftol=1.0e-10,
                            gtol=1.0e-10,
                        )
                        x = np.array(sol_polish.x, dtype=float)
                        success = bool(success or sol_polish.success)
                        nfev += int(getattr(sol_polish, "nfev", 0))
                    except Exception:
                        pass
                else:
                    sol2 = least_squares(
                        closure_residual,
                        x0,
                        method="trf",
                        max_nfev=int(root_maxfev),
                        xtol=1.0e-10,
                        ftol=1.0e-10,
                        gtol=1.0e-10,
                    )
                    x = np.array(sol2.x, dtype=float)
                    success = bool(sol2.success)
                    msg = getattr(sol2, "message", "")
                    nfev = int(getattr(sol2, "nfev", -1))
            except Exception as exc:
                sys.stderr.write(
                    "[re_helix_cck]   attempt %d failed with exception: %s\n"
                    % (attempt_idx, exc)
                )
                continue

            r = closure_residual(x)
            closure_norm = float(np.linalg.norm(r))
            gscore = float(geometry_score(x))

            sys.stderr.write(
                "[re_helix_cck]   attempt %d: closure=%.6f, geom=%.6f, success=%s, nfev=%d\n"
                % (attempt_idx, closure_norm, gscore, str(success), nfev)
            )
            candidates.append((closure_norm, gscore, x, success, nfev, str(msg)))

        if not candidates:
            sys.stderr.write(
                "[re_helix_cck]   Warning: no closure candidates found; keeping base geometry.\n"
            )
            continue

        # Primary ranking = closure residual norm.
        candidates.sort(key=lambda row: (row[0], row[1]))
        best_closure = candidates[0][0]
        closure_cut = max(best_closure + closure_abs_tol, best_closure * (1.0 + closure_rel_tol))
        near_candidates = [row for row in candidates if row[0] <= closure_cut]
        if not near_candidates:
            near_candidates = [candidates[0]]

        # Secondary filter/ranking = geometry score within the near-best closure set.
        near_candidates.sort(key=lambda row: (row[1], row[0]))
        chosen_closure, chosen_geom, chosen_x, chosen_success, chosen_nfev, chosen_msg = near_candidates[0]

        sys.stderr.write(
            "[re_helix_cck]   selected solution: closure=%.6f, geom=%.6f, success=%s\n"
            % (chosen_closure, chosen_geom, str(chosen_success))
        )

        best_transforms = transforms_from_params(chosen_x)
        apply_component_solution_in_place(rec_list, solver_data, comp, best_transforms)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Align nucleic-acid helices with a HolT-Hunter-like closure solve for "
            "cyclic components, then apply reciprocal exchanges."
        )
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%s %s" % (SOFTWARE_NAME, SOFTWARE_VERSION),
        help="Show the script version and exit.",
    )
    parser.add_argument("pdb_in", help="Input PDB file.")
    parser.add_argument(
        "ops",
        nargs="+",
        help=(
            "Tokens including optional helix definitions and exchange specs. "
            "Examples:\n"
            "  (AB) (CD) 30A 8D d 13B 24C s\n"
            "  (ABMN) 30A 8D d 13B 24C s\n"
            "  30A 8D d 13B 24C s"
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="pdb_out_base",
        default=None,
        help=(
            "Base name for output files (extension optional). "
            "Outputs will be <base>_aligned_cck.pdb and <base>_aligned_cck_rex.pdb. "
            "Default base: input filename without extension."
        ),
    )
    parser.add_argument(
        "--axis_dist",
        type=float,
        default=22.0,
        help="Target distance between helical axes in Å (default: 22.0).",
    )
    parser.add_argument(
        "--axis_parallel",
        choices=["y", "Y", "n", "N"],
        default="y",
        help=(
            "If 'y' (default), helices are treated with the original tree alignment. "
            "If 'n', cyclic components are additionally refined with an explicit "
            "closure solve inspired by HolT Hunter."
        ),
    )
    parser.add_argument(
        "--fix",
        dest="fix_chain",
        default=None,
        help=(
            "Chain ID whose helix group should remain fixed during alignment. "
            "For example, '--fix A' keeps the helix group containing chain A "
            "as the immobile reference within its connected component."
        ),
    )
    parser.add_argument(
        "--replicate",
        action="store_true",
        help=(
            "Replicate the entire set of chains (same semantics as in re_helix.py). "
            "If the input PDB appears to contain exactly one helix component, "
            "replication is also enabled automatically."
        ),
    )
    parser.add_argument(
        "--re",
        action="store_true",
        help=(
            "Legacy flag: reciprocal exchanges are always applied and an "
            "*_aligned_cck_rex.pdb file is always written."
        ),
    )
    parser.add_argument(
        "--cir_shift",
        type=int,
        default=8,
        help="Residue shift for circular strands during reciprocal exchange (default: 8).",
    )
    parser.add_argument(
        "--w_pp",
        type=float,
        default=1.0,
        help="Weight for P–P distance score used as the secondary ranking/filter (default: 1.0).",
    )
    parser.add_argument(
        "--w_line",
        type=float,
        default=1.0,
        help=(
            "Weight for line-topology score on cycle edges, used only as the secondary "
            "ranking/filter (default: 1.0)."
        ),
    )
    parser.add_argument(
        "--w_axis",
        type=float,
        default=1.0e4,
        help=(
            "Weight for (d_axis - axis_dist)^2 in the secondary geometry score "
            "(default: 1.0e4)."
        ),
    )
    parser.add_argument(
        "--max_root_attempts",
        type=int,
        default=20,
        help=(
            "Number of deterministic closure-root starts per cyclic component "
            "(default: 20), mirroring HolT Hunter's multiple FindRoot seeds."
        ),
    )
    parser.add_argument(
        "--root_maxfev",
        type=int,
        default=800,
        help="Maximum function evaluations per closure solve attempt (default: 800).",
    )
    parser.add_argument(
        "--closure_rel_tol",
        type=float,
        default=0.10,
        help=(
            "Relative tolerance for defining the near-best closure set before applying "
            "the secondary geometry ranking/filter (default: 0.10)."
        ),
    )
    parser.add_argument(
        "--closure_abs_tol",
        type=float,
        default=1.0e-3,
        help=(
            "Absolute tolerance for defining the near-best closure set before applying "
            "the secondary geometry ranking/filter (default: 1e-3)."
        ),
    )
    parser.add_argument(
        "--geom_polish",
        choices=["y", "Y", "n", "N"],
        default="y",
        help=(
            "If 'y' (default), run a final geometry-based cyclic polish using "
            "re_helix_ccgV3.py after the closure solve. This fixes the "
            "live-structure application bug from the earlier cck script and "
            "usually brings the result close to the better ccg geometry."
        ),
    )
    parser.add_argument(
        "--geom_attempts",
        type=int,
        default=4,
        help=(
            "Number of deterministic geometry-polish starts per cyclic component "
            "(default: 4). Ignored if --geom_polish n."
        ),
    )
    parser.add_argument(
        "--geom_max_nfev",
        type=int,
        default=1200,
        help=(
            "Maximum function evaluations in the coarse geometry-polish stage "
            "(default: 1200). Ignored if --geom_polish n."
        ),
    )
    parser.add_argument(
        "--twist_report",
        choices=["y", "Y", "n", "N"],
        default="y",
        help=(
            "If 'y' (default), write a post-alignment twist_rod/twist_helix "
            "diagnostic TSV file before reciprocal exchange is applied."
        ),
    )
    parser.add_argument(
        "--twist_report_file",
        default=None,
        help=(
            "Optional path for the twist diagnostic TSV file. Default: "
            "<base>_aligned_cck_twist.tsv."
        ),
    )
    parser.add_argument(
        "--twist_helical_repeat",
        type=float,
        default=10.5,
        help=(
            "Fallback helical repeat in bp/turn. V3 measures twist_helix "
            "directly from RE-site phosphate rotations; this value is used only "
            "when observed phosphate-step information is unavailable (default: 10.5)."
        ),
    )
    parser.add_argument(
        "--twist_pairing",
        choices=["consecutive", "all"],
        default="consecutive",
        help=(
            "Which exchange-site pairs to report on each helix: consecutive sites "
            "along the helix axis (default), or all pairwise combinations."
        ),
    )

    args = parser.parse_args()

    alignmod.reset_helix_axis_overrides()

    # Determine base name for outputs
    if args.pdb_out_base is None:
        pdb_in_lower = args.pdb_in.lower()
        if pdb_in_lower.endswith('.pdb.txt'):
            base = args.pdb_in[:-8]
        elif pdb_in_lower.endswith('.pdb'):
            base = args.pdb_in[:-4]
        elif pdb_in_lower.endswith('.txt'):
            base = args.pdb_in[:-4]
        else:
            base = args.pdb_in
    else:
        base = args.pdb_out_base
        base_lower = base.lower()
        if base_lower.endswith('.pdb.txt'):
            base = base[:-8]
        elif base_lower.endswith('.pdb'):
            base = base[:-4]
        elif base_lower.endswith('.txt'):
            base = base[:-4]

    pdb_out_aligned = base + '_aligned_cck.pdb'
    pdb_out_aligned_rex = base + '_aligned_cck_rex.pdb'

    # Flatten ops tokens (support quoted strings)
    if len(args.ops) == 1 and (' ' in args.ops[0] or '\t' in args.ops[0]):
        tokens = args.ops[0].split()
    else:
        tokens = args.ops

    helix_defs, exch_tokens = parse_helix_definition_tokens(tokens)

    try:
        exchange_specs = parse_exchange_specs(exch_tokens)
    except ValueError as exc:
        print("Error parsing exchange specifications: %s" % exc, file=sys.stderr)
        sys.exit(1)

    rec_list = []  # type: List[pdb_atom_record]
    try:
        with open(args.pdb_in, 'r') as fin:
            file2rec(fin, rec_list)
    except OSError as exc:
        print("Error reading PDB file '%s': %s" % (args.pdb_in, exc), file=sys.stderr)
        sys.exit(1)

    if not rec_list:
        print("No ATOM/HETATM/TER records found in input PDB.", file=sys.stderr)
        sys.exit(1)

    # Minimal maps for replication decision
    _, _, chain_to_P_atoms0 = build_nucleic_acid_maps(rec_list)
    if not chain_to_P_atoms0:
        print("No P atoms found in input PDB; alignment cannot proceed.", file=sys.stderr)
        sys.exit(1)

    # Determine how many helices the input PDB appears to contain
    if helix_defs:
        helices0 = sorted(set(tuple(sorted(h)) for h in helix_defs))
    else:
        chain_to_helix0 = compute_chain_partner_map(chain_to_P_atoms0)
        helices0 = sorted(set(chain_to_helix0[ch] for ch in chain_to_helix0))

    auto_single_helix = (len(helices0) == 1)
    replicate_active = args.replicate or auto_single_helix

    if replicate_active:
        try:
            exchange_specs, helix_defs_repl = replicate_all_chains(
                rec_list,
                helices0,
                exchange_specs,
                helix_defs if helix_defs else None,
            )
        except Exception as exc:
            print("Error during helix replication: %s" % exc, file=sys.stderr)
            sys.exit(1)
        helix_defs_for_align = helix_defs_repl
    else:
        helix_defs_for_align = helix_defs if helix_defs else None

    axis_parallel_flag = args.axis_parallel.lower() == 'y'

    # Detect whether the user-supplied terminal pattern itself indicates a
    # replicated symmetric cycle (for example C2-F21, F2-I21, I2-C21). If so,
    # we later make sure a geometry-based cyclic polish is available before the
    # exact C_n symmetry projection.
    try:
        symmetric_cycle_requested = has_detected_symmetric_cycle(
            rec_list,
            exchange_specs,
            explicit_helices=helix_defs_for_align,
            fix_chain=args.fix_chain,
        )
    except Exception as exc:
        print("Error while detecting symmetric terminal cycles: %s" % exc, file=sys.stderr)
        sys.exit(1)

    # Baseline alignment strategy:
    #   - For axis_parallel=y, use the original align_helices_for_exchanges directly.
    #   - For axis_parallel=n, handle components ourselves so that cyclic components
    #     can use the HolT-Hunter-like closure solve after a DTP base alignment.
    if axis_parallel_flag:
        try:
            alignmod.align_helices_for_exchanges(
                rec_list,
                exchange_specs,
                axis_dist=args.axis_dist,
                axis_parallel_flag=axis_parallel_flag,
                explicit_helices=helix_defs_for_align,
                fix_chain=args.fix_chain,
            )
        except Exception as exc:
            print("Error during baseline helix alignment: %s" % exc, file=sys.stderr)
            sys.exit(1)
    else:
        # Manual component-wise handling
        _, _, chain_to_P_atoms = build_nucleic_acid_maps(rec_list)
        if helix_defs_for_align:
            chain_to_helix = build_chain_to_helix_from_defs(helix_defs_for_align, chain_to_P_atoms)
            sys.stderr.write(
                "[re_helix_cck] Using user-defined helices: "
                + ", ".join(helix_id_str(tuple(sorted(h))) for h in helix_defs_for_align)
                + "\n"
            )
        else:
            chain_to_helix = compute_chain_partner_map(chain_to_P_atoms)
            seen_helices = sorted(set(chain_to_helix[ch] for ch in chain_to_helix))
            sys.stderr.write(
                "[re_helix_cck] Using automatic helix detection. Helices: "
                + ", ".join(helix_id_str(h) for h in seen_helices)
                + "\n"
            )

        fix_helix = None  # type: Optional[HelixID]
        if args.fix_chain is not None:
            if args.fix_chain not in chain_to_helix:
                print(
                    "Requested fixed chain '%s' is not part of any helix." % args.fix_chain,
                    file=sys.stderr,
                )
                sys.exit(1)
            fix_helix = chain_to_helix[args.fix_chain]
            sys.stderr.write(
                "[re_helix_cck] Helix %s containing chain '%s' will be used as root (fixed).\n"
                % (helix_id_str(fix_helix), args.fix_chain)
            )

        helix_pair_data, adjacency = build_helix_pair_graph(exchange_specs, chain_to_helix)
        if helix_pair_data:
            components = build_helix_components(adjacency)
            for comp in components:
                if not comp:
                    continue
                root_h = fix_helix if (fix_helix is not None and fix_helix in comp) else sorted(comp)[0]
                order, parent = build_tree_for_component(comp, adjacency, root_h)
                _, residue_to_P_atom, chain_to_P_atoms = build_nucleic_acid_maps(rec_list)

                edges = component_edges(comp, adjacency)
                has_cycle = len(edges) > max(0, len(comp) - 1)
                sys.stderr.write(
                    "[re_helix_cck] Baseline handling component: %s; %s\n"
                    % (
                        ", ".join(helix_id_str(h) for h in sorted(comp)),
                        "cyclic" if has_cycle else "acyclic",
                    )
                )

                if has_cycle:
                    # DTP base for closure solver
                    pairwise_align_component_base_for_closure(
                        rec_list,
                        comp,
                        helix_pair_data,
                        chain_to_helix,
                        args.axis_dist,
                        order,
                        parent,
                    )
                else:
                    pairwise_align_component_sequential(
                        rec_list,
                        comp,
                        helix_pair_data,
                        chain_to_helix,
                        args.axis_dist,
                        axis_parallel_flag=False,
                        residue_to_P_atom=residue_to_P_atom,
                        chain_to_P_atoms=chain_to_P_atoms,
                        order=order,
                        parent=parent,
                    )

    # HolT-Hunter-like closure refinement on cyclic components (only axis_parallel=n)
    if not axis_parallel_flag:
        try:
            refine_cyclic_components_closure(
                rec_list,
                exchange_specs,
                axis_dist=args.axis_dist,
                axis_parallel_flag=axis_parallel_flag,
                explicit_helices=helix_defs_for_align,
                fix_chain=args.fix_chain,
                w_pp=args.w_pp,
                w_line=args.w_line,
                w_axis=args.w_axis,
                max_root_attempts=args.max_root_attempts,
                root_maxfev=args.root_maxfev,
                closure_rel_tol=args.closure_rel_tol,
                closure_abs_tol=args.closure_abs_tol,
            )
        except Exception as exc:
            print("Error during closure-based cyclic refinement: %s" % exc, file=sys.stderr)
            sys.exit(1)

    # Optional geometry-based cyclic polish, seeded by the closure result.
    run_geom_polish = (not axis_parallel_flag) and (
        (args.geom_polish.lower() == 'y') or symmetric_cycle_requested
    )
    if run_geom_polish:
        ccgmod = _load_ccg_module()
        if ccgmod is None:
            sys.stderr.write(
                "[re_helix_cck] Warning: re_helix_ccgV3.py was not found "
                "or could not be loaded; skipping final geometry polish.\n"
            )
        else:
            try:
                if hasattr(ccgmod, "alignmod") and hasattr(ccgmod.alignmod, "set_helix_axis_overrides"):
                    ccgmod.alignmod.set_helix_axis_overrides(alignmod.get_helix_axis_overrides())
                if symmetric_cycle_requested and args.geom_polish.lower() != 'y':
                    sys.stderr.write(
                        "[re_helix_cck] Symmetric terminal cycle detected; running geometry-based cyclic polish even though --geom_polish n was requested.\n"
                    )
                sys.stderr.write(
                    "[re_helix_cck] Running final geometry-based cyclic polish from re_helix_ccgV3.py.\n"
                )
                ccgmod.refine_cyclic_components_geometry(
                    rec_list,
                    exchange_specs,
                    axis_dist=args.axis_dist,
                    axis_parallel_flag=axis_parallel_flag,
                    explicit_helices=helix_defs_for_align,
                    fix_chain=args.fix_chain,
                    w_pp=args.w_pp,
                    w_line=args.w_line,
                    w_axis=args.w_axis,
                    geom_attempts=args.geom_attempts,
                    geom_max_nfev=args.geom_max_nfev,
                )
            except Exception as exc:
                print(
                    "Error during final geometry-based cyclic polish: %s" % exc,
                    file=sys.stderr,
                )
                sys.exit(1)

    # If the terminal exchange pattern itself indicates a replicated symmetric
    # cycle (for example C2-F21, F2-I21, I2-C21), project the aligned model
    # onto exact C_n symmetry before writing the pre-RE structure.
    try:
        enforce_detected_cyclic_symmetry(
            rec_list,
            exchange_specs,
            explicit_helices=helix_defs_for_align,
            fix_chain=args.fix_chain,
        )
    except Exception as exc:
        print(
            "Error during exact symmetry projection: %s" % exc,
            file=sys.stderr,
        )
        sys.exit(1)

    # Post-alignment twist_rod / twist_helix diagnostics are calculated before
    # reciprocal exchange changes chain IDs and residue numbering.
    if args.twist_report.lower() == 'y':
        twist_report_path = args.twist_report_file
        if twist_report_path is None:
            twist_report_path = base + '_aligned_cck_twist.tsv'
        try:
            twist_rows, twist_notes = compute_twist_rod_helix_rows(
                rec_list,
                exchange_specs,
                explicit_helices=helix_defs_for_align,
                helical_repeat=args.twist_helical_repeat,
                pairing=args.twist_pairing,
            )
            write_twist_report_tsv(
                twist_report_path,
                twist_rows,
                twist_notes,
                args.twist_pairing,
            )
        except Exception as exc:
            print("Error writing twist_rod/twist_helix report: %s" % exc, file=sys.stderr)
            sys.exit(1)
        print(
            "Wrote twist_rod/twist_helix report with %d segment row(s) to '%s'."
            % (len(twist_rows), twist_report_path)
        )

    # Write aligned PDB (pre-RE)
    try:
        with open(pdb_out_aligned, 'w') as fout:
            rec2file(rec_list, fout, reorder_serial=True)
    except OSError as exc:
        print("Error writing aligned PDB '%s': %s" % (pdb_out_aligned, exc), file=sys.stderr)
        sys.exit(1)

    print("Wrote aligned PDB (pre-RE, closure-based cyclic refinement) to '%s'." % pdb_out_aligned)

    # Apply reciprocal exchange and write final PDB
    try:
        rec_rex = apply_reciprocal_exchanges_in_memory(
            rec_list,
            exchange_specs,
            cir_shift=args.cir_shift,
        )
    except Exception as exc:
        print("Error during reciprocal exchange stage: %s" % exc, file=sys.stderr)
        sys.exit(1)

    try:
        with open(pdb_out_aligned_rex, 'w') as fout:
            if isinstance(rec_rex, tuple) and len(rec_rex) == 2:
                atom_recs, link_recs = rec_rex
                rec2file_link(atom_recs, link_recs, fout, reorder_serial=True)
            else:
                rec2file(rec_rex, fout, reorder_serial=True)
    except OSError as exc:
        print("Error writing aligned+RE PDB '%s': %s" % (pdb_out_aligned_rex, exc), file=sys.stderr)
        sys.exit(1)

    print("Wrote aligned+RE PDB (closure-based cyclic refinement) to '%s'." % pdb_out_aligned_rex)


if __name__ == '__main__':
    main()
