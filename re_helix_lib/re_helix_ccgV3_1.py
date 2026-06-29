#!/usr/bin/env python3
"""
re_helix_ccgV3_1.py

Pure geometry-based cyclic-constraint alignment for reciprocal-exchange helices.

What it does
------------
This script uses the current pairwise/tree alignment logic from re_helix.py to
generate a deterministic starting geometry, and then
replaces cyclic sampling/screening with a direct geometric closure solve for
cyclic helix components. The cyclic refinement is inspired by the logic of
HolT Hunter: instead of enumerating angle grids, it minimises a geometric
closure residual with scipy.optimize.least_squares.

Inputs
------
- Input PDB file.
- Optional explicit helix definitions like (AB) (CD) or (ABMN).
- Reciprocal-exchange style operations, e.g.:
      26A 9C d 26C 9E d 26E 9A d

Outputs
-------
- <base>_aligned_ccg.pdb      : aligned structure before reciprocal exchange
- <base>_aligned_ccg_rex.pdb  : aligned structure after reciprocal exchange

Example
-------
python re_helix_ccgV3_1.py TT_helixAB_33.pdb \
    26A 9C d 26C 9E d 26E 9A d \
    --axis_parallel n --axis_dist 23
"""

import argparse
import importlib.util
import math
import os
from pathlib import Path
import sys
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
from scipy.optimize import least_squares

SOFTWARE_NAME = "re_helix_ccg"
SOFTWARE_VERSION = "V3.1"

# Make the helper folder and repository root importable when loading re_helix.py.
HERE = Path(__file__).resolve().parent
PARENT = HERE.parent
for path in (str(PARENT), str(HERE)):
    if path not in sys.path:
        sys.path.insert(0, path)

from re_helix_lib.edit_pdb_atom import file2rec, rec2file, pdb_atom_record  # type: ignore
from re_helix_lib.edit_pdb_link import rec2file_link  # type: ignore


def _load_module(module_name: str, filenames: List[str]):
    """Load a Python module from one of several candidate filenames.

    This keeps the bundled experimental script compatible with the renamed
    re_helix entry point while still tolerating legacy local helper filenames.
    """
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
parse_axis_range_specs = alignmod.parse_axis_range_specs
parse_axis_move_specs = alignmod.parse_axis_move_specs
pair_axis_move_definitions = alignmod.pair_axis_move_definitions
resolve_axis_range_definitions = alignmod.resolve_axis_range_definitions
resolve_axis_move_definitions = alignmod.resolve_axis_move_definitions
translate_axis_range_definitions_for_replication = alignmod.translate_axis_range_definitions_for_replication
build_axis_coupled_helix_defs = alignmod.build_axis_coupled_helix_defs
register_axis_coupling_overrides = alignmod.register_axis_coupling_overrides
apply_axis_coupling_to_chain_map = alignmod.apply_axis_coupling_to_chain_map
set_helix_axis_range_definitions = alignmod.set_helix_axis_range_definitions
set_helix_axis_move_definitions = alignmod.set_helix_axis_move_definitions
build_nucleic_acid_maps = alignmod.build_nucleic_acid_maps
compute_chain_partner_map = alignmod.compute_chain_partner_map
build_chain_to_helix_from_defs = alignmod.build_chain_to_helix_from_defs
build_helix_pair_graph = alignmod.build_helix_pair_graph
compute_helix_axis = alignmod.compute_helix_axis
align_helices_for_exchanges = alignmod.align_helices_for_exchanges
apply_reciprocal_exchanges_in_memory = alignmod.apply_reciprocal_exchanges_in_memory
replicate_all_chains = alignmod.replicate_all_chains
helix_id_str = alignmod.helix_id_str
atom_moves_with_helix = alignmod.atom_moves_with_helix
v_add = alignmod.v_add
v_sub = alignmod.v_sub
v_scale = alignmod.v_scale
v_dot = alignmod.v_dot
v_cross = alignmod.v_cross
v_length = alignmod.v_length
v_norm = alignmod.v_norm


# ---------------------------------------------------------------------------
# Small geometry helpers
# ---------------------------------------------------------------------------


def build_helix_components(adjacency: Dict[HelixID, Set[HelixID]]) -> List[Set[HelixID]]:
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



def component_edges(component: Set[HelixID], adjacency: Dict[HelixID, Set[HelixID]]) -> Set[frozenset]:
    """Return all undirected edges inside a component."""
    edges = set()  # type: Set[frozenset]
    for h in component:
        for nb in adjacency.get(h, set()):
            if nb not in component:
                continue
            edges.add(frozenset((h, nb)))
    return edges



def skew_matrix(v: np.ndarray) -> np.ndarray:
    """Return the 3x3 skew-symmetric matrix for a 3-vector."""
    return np.array(
        [
            [0.0, -v[2], v[1]],
            [v[2], 0.0, -v[0]],
            [-v[1], v[0], 0.0],
        ],
        dtype=float,
    )



def rotvec_to_matrix(rv: np.ndarray) -> np.ndarray:
    """Convert a rotation vector to a 3x3 rotation matrix.

    The rotation vector direction is the axis and its norm is the angle.
    """
    angle = float(np.linalg.norm(rv))
    if angle < 1.0e-12:
        return np.eye(3, dtype=float)
    k = rv / angle
    K = skew_matrix(k)
    return np.eye(3, dtype=float) + math.sin(angle) * K + (1.0 - math.cos(angle)) * np.dot(K, K)



def line_distance(
    p1: np.ndarray,
    d1: np.ndarray,
    p2: np.ndarray,
    d2: np.ndarray,
    eps: float = 1.0e-9,
) -> float:
    """Minimum distance between two lines p1+t d1 and p2+s d2."""
    d1u = d1 / max(np.linalg.norm(d1), eps)
    d2u = d2 / max(np.linalg.norm(d2), eps)
    cross_d = np.cross(d1u, d2u)
    norm_cross = float(np.linalg.norm(cross_d))
    w = p2 - p1
    if norm_cross > eps:
        n_unit = cross_d / norm_cross
        return abs(float(np.dot(w, n_unit)))
    # nearly parallel: distance from p2 to line1
    proj = float(np.dot(w, d1u))
    w_perp = w - proj * d1u
    return float(np.linalg.norm(w_perp))



def axis_distance(
    c1: np.ndarray,
    u1: np.ndarray,
    c2: np.ndarray,
    u2: np.ndarray,
    eps: float = 1.0e-9,
) -> float:
    """Minimum distance between two axes (infinite lines)."""
    u1u = u1 / max(np.linalg.norm(u1), eps)
    u2u = u2 / max(np.linalg.norm(u2), eps)
    w0 = c2 - c1
    cross_u = np.cross(u1u, u2u)
    norm_cross = float(np.linalg.norm(cross_u))
    if norm_cross > eps:
        n_unit = cross_u / norm_cross
        return abs(float(np.dot(w0, n_unit)))
    proj = float(np.dot(w0, u1u))
    w_perp = w0 - proj * u1u
    return float(np.linalg.norm(w_perp))


# ---------------------------------------------------------------------------
# Direct geometric cyclic refinement
# ---------------------------------------------------------------------------


def build_cycle_edge_set(component: Set[HelixID], root: HelixID, adjacency: Dict[HelixID, Set[HelixID]]) -> Set[frozenset]:
    """Return the set of cycle edges (edges not in a BFS spanning tree)."""
    seen = set([root])
    queue = [root]
    tree_edges = set()  # type: Set[frozenset]
    while queue:
        h = queue.pop(0)
        for nb in sorted(adjacency.get(h, set())):
            if nb not in component:
                continue
            if nb in seen:
                continue
            seen.add(nb)
            queue.append(nb)
            tree_edges.add(frozenset((h, nb)))
    return component_edges(component, adjacency) - tree_edges



def build_component_solver_data(
    rec_list: List[pdb_atom_record],
    component: Set[HelixID],
    helix_pair_data,
    chain_to_helix: Dict[str, HelixID],
    cycle_edges: Set[frozenset],
) -> Dict[str, object]:
    """Precompute data needed by the geometric least-squares solver."""
    # Snapshot atom lists and base coordinates per helix.
    atoms_by_helix = {}  # type: Dict[HelixID, List[pdb_atom_record]]
    base_coords_by_helix = {}  # type: Dict[HelixID, np.ndarray]

    for h in component:
        atoms = []  # type: List[pdb_atom_record]
        for atom in rec_list:
            if atom.recordName not in ("ATOM", "HETATM"):
                continue
            if atom_moves_with_helix(atom, h):
                atoms.append(atom)
        atoms_by_helix[h] = atoms
        base_coords_by_helix[h] = np.array([[a.x, a.y, a.z] for a in atoms], dtype=float)

    # Current helix axes from current coordinates (after baseline pairwise alignment).
    _, _, chain_to_P_atoms = build_nucleic_acid_maps(rec_list)
    helix_axis = {}  # type: Dict[HelixID, Tuple[np.ndarray, np.ndarray]]
    for h in component:
        axis_dir, center = compute_helix_axis(chain_to_P_atoms, h)
        helix_axis[h] = (np.array(v_norm(axis_dir), dtype=float), np.array(center, dtype=float))

    # Current P-atom coordinates used as the base reference for rigid transforms.
    _, residue_to_P_atom, _ = build_nucleic_acid_maps(rec_list)

    p_pairs = []  # type: List[Tuple[HelixID, HelixID, np.ndarray, np.ndarray, bool]]
    axis_pairs = set()  # type: Set[Tuple[HelixID, HelixID]]

    for _key, entry in helix_pair_data.items():
        h1 = entry["helix1"]
        h2 = entry["helix2"]
        if h1 not in component or h2 not in component:
            continue
        pairs = entry["pairs"]
        axis_pairs.add((h1, h2))
        ekey = frozenset((h1, h2))
        is_cycle = ekey in cycle_edges
        for (c1, r1), (c2, r2) in pairs:
            a1 = residue_to_P_atom.get((c1, r1))
            a2 = residue_to_P_atom.get((c2, r2))
            if a1 is None or a2 is None:
                continue
            p_pairs.append(
                (
                    h1,
                    h2,
                    np.array([a1.x, a1.y, a1.z], dtype=float),
                    np.array([a2.x, a2.y, a2.z], dtype=float),
                    is_cycle,
                )
            )

    return {
        "atoms_by_helix": atoms_by_helix,
        "base_coords_by_helix": base_coords_by_helix,
        "helix_axis": helix_axis,
        "p_pairs": p_pairs,
        "axis_pairs": list(axis_pairs),
    }



def make_component_residual(
    solver_data: Dict[str, object],
    component: Set[HelixID],
    root: HelixID,
    axis_dist: float,
    w_pp: float,
    w_line: float,
    w_axis: float,
):
    """Build a least-squares residual function for one cyclic component."""
    helices = sorted(component)
    moving_helices = [h for h in helices if h != root]
    h_index = {h: i for i, h in enumerate(moving_helices)}

    helix_axis = solver_data["helix_axis"]  # type: ignore[assignment]
    p_pairs = solver_data["p_pairs"]  # type: ignore[assignment]
    axis_pairs = solver_data["axis_pairs"]  # type: ignore[assignment]

    def transform_for_helix(params: np.ndarray, h: HelixID) -> Tuple[np.ndarray, np.ndarray]:
        if h == root:
            return np.eye(3, dtype=float), np.zeros(3, dtype=float)
        idx = h_index[h]
        tx, ty, tz, rx, ry, rz = params[6 * idx: 6 * (idx + 1)]
        R = rotvec_to_matrix(np.array([rx, ry, rz], dtype=float))
        t = np.array([tx, ty, tz], dtype=float)
        return R, t

    def transform_point(R: np.ndarray, t: np.ndarray, p: np.ndarray) -> np.ndarray:
        return np.dot(R, p) + t

    def transform_axis(R: np.ndarray, t: np.ndarray, center: np.ndarray, axis_dir: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        return np.dot(R, center) + t, np.dot(R, axis_dir)

    def residual(params: np.ndarray) -> np.ndarray:
        residuals = []  # type: List[float]

        # Precompute transformed axes per helix
        t_axes = {}  # type: Dict[HelixID, Tuple[np.ndarray, np.ndarray]]
        for h in helices:
            R, t = transform_for_helix(params, h)
            base_axis_dir, base_center = helix_axis[h]
            center_t, axis_t = transform_axis(R, t, base_center, base_axis_dir)
            t_axes[h] = (center_t, axis_t)

        # P-pair residuals: encode squared distance through Cartesian residuals.
        if w_pp > 0.0:
            scale_pp = math.sqrt(w_pp)
            for h1, h2, p1, p2, is_cycle in p_pairs:
                R1, t1 = transform_for_helix(params, h1)
                R2, t2 = transform_for_helix(params, h2)
                q1 = transform_point(R1, t1, p1)
                q2 = transform_point(R2, t2, p2)
                dq = scale_pp * (q1 - q2)
                residuals.extend(dq.tolist())

        # Line-topology residuals: only for cycle-edge P-pairs.
        if w_line > 0.0:
            scale_line = math.sqrt(w_line)
            for h1, h2, p1, p2, is_cycle in p_pairs:
                if not is_cycle:
                    continue
                R1, t1 = transform_for_helix(params, h1)
                R2, t2 = transform_for_helix(params, h2)
                q1 = transform_point(R1, t1, p1)
                q2 = transform_point(R2, t2, p2)
                c1, u1 = t_axes[h1]
                c2, u2 = t_axes[h2]

                # radial line through q1 perpendicular to axis1
                delta1 = q1 - c1
                s1 = float(np.dot(delta1, u1))
                A1 = c1 + s1 * u1
                r1 = q1 - A1
                nr1 = float(np.linalg.norm(r1))
                if nr1 < 1.0e-9:
                    continue
                r1u = r1 / nr1

                # radial line through q2 perpendicular to axis2
                delta2 = q2 - c2
                s2 = float(np.dot(delta2, u2))
                A2 = c2 + s2 * u2
                r2 = q2 - A2
                nr2 = float(np.linalg.norm(r2))
                if nr2 < 1.0e-9:
                    continue
                r2u = r2 / nr2

                d_line = line_distance(q1, r1u, q2, r2u)
                residuals.append(scale_line * d_line)

        # Axis-distance residuals for all helix pairs with constraints.
        if w_axis > 0.0:
            scale_axis = math.sqrt(w_axis)
            for h1, h2 in axis_pairs:
                c1, u1 = t_axes[h1]
                c2, u2 = t_axes[h2]
                d_ax = axis_distance(c1, u1, c2, u2)
                residuals.append(scale_axis * (d_ax - axis_dist))

        if not residuals:
            return np.array([0.0], dtype=float)
        return np.array(residuals, dtype=float)

    return moving_helices, residual



def apply_component_solution_in_place(
    rec_list: List[pdb_atom_record],
    solver_data: Dict[str, object],
    component: Set[HelixID],
    root: HelixID,
    params: np.ndarray,
) -> None:
    """Apply the solved 6-DOF transforms for one component to rec_list."""
    helices = sorted(component)
    moving_helices = [h for h in helices if h != root]
    h_index = {h: i for i, h in enumerate(moving_helices)}
    atoms_by_helix = solver_data["atoms_by_helix"]  # type: ignore[assignment]
    base_coords_by_helix = solver_data["base_coords_by_helix"]  # type: ignore[assignment]

    for h in helices:
        atoms = atoms_by_helix[h]
        base_xyz = base_coords_by_helix[h]
        if h == root:
            R = np.eye(3, dtype=float)
            t = np.zeros(3, dtype=float)
        else:
            idx = h_index[h]
            tx, ty, tz, rx, ry, rz = params[6 * idx: 6 * (idx + 1)]
            R = rotvec_to_matrix(np.array([rx, ry, rz], dtype=float))
            t = np.array([tx, ty, tz], dtype=float)
        new_xyz = np.dot(base_xyz, R.T) + t  # shape (N,3)
        for atom, xyz in zip(atoms, new_xyz):
            atom.update_xyz(float(xyz[0]), float(xyz[1]), float(xyz[2]))



def refine_cyclic_components_geometry(
    rec_list: List[pdb_atom_record],
    exchange_specs,
    axis_dist: float,
    axis_parallel_flag: bool,
    explicit_helices: Optional[List[HelixID]] = None,
    fix_chain: Optional[str] = None,
    w_pp: float = 1.0,
    w_line: float = 1.0,
    w_axis: float = 1.0e4,
    geom_attempts: int = 4,
    geom_max_nfev: int = 500,
) -> None:
    """Refine cyclic components by direct geometric least-squares.

    This is the pure geometry-based step inspired by HolT Hunter's direct
    closure solve: rather than screen angle grids, we directly solve for
    rigid-body transforms that minimise a geometric residual for the whole
    cyclic component.
    """
    # Current helix assignment from the already-aligned structure.
    _, _, chain_to_P_atoms = build_nucleic_acid_maps(rec_list)
    if explicit_helices:
        chain_to_helix = build_chain_to_helix_from_defs(explicit_helices, chain_to_P_atoms)
    else:
        chain_to_helix = compute_chain_partner_map(chain_to_P_atoms)
    chain_to_helix = apply_axis_coupling_to_chain_map(
        chain_to_helix,
        alignmod.get_helix_axis_range_definitions(),
        alignmod.get_helix_axis_move_definitions(),
        chain_to_P_atoms,
    )

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

        root = fix_helix if (fix_helix is not None and fix_helix in comp) else sorted(comp)[0]
        cycle_edges = build_cycle_edge_set(comp, root, adjacency)

        def _edge_label(edge: frozenset) -> str:
            h1, h2 = sorted(list(edge))
            return "%s--%s" % (helix_id_str(h1), helix_id_str(h2))

        sys.stderr.write(
            "[re_helix_ccg] Geometric refinement for cyclic component: "
            + ", ".join(helix_id_str(h) for h in sorted(comp))
            + "; cycle edges="
            + (", ".join(sorted(_edge_label(edge) for edge in cycle_edges)) if cycle_edges else "(none)")
            + "\n"
        )

        solver_data = build_component_solver_data(
            rec_list,
            comp,
            helix_pair_data,
            chain_to_helix,
            cycle_edges,
        )
        moving_helices, residual = make_component_residual(
            solver_data,
            comp,
            root,
            axis_dist,
            w_pp=w_pp,
            w_line=w_line,
            w_axis=w_axis,
        )

        n_params = 6 * len(moving_helices)
        if n_params == 0:
            continue

        # HolT Hunter uses multiple starting guesses for root finding. The
        # original ccg script only used a handful of tiny rotational biases,
        # which was too weak for the TT_helixAB_33 triangle: the correct
        # tensegrity root lives in a distant basin. Here we use a broader,
        # deterministic multi-start strategy with two phases:
        #   (1) coarse least-squares solves from geometry-diverse seeds,
        #   (2) polish the best few seeds with a larger nfev budget.
        # This stays purely geometry-based (no angle grid screening), but
        # searches roots much more like HolT Hunter's multi-seed FindRoot/root.
        seeds = []  # type: List[np.ndarray]
        seeds.append(np.zeros(n_params, dtype=float))

        # Deterministic random geometric seeds: moderate translations + arbitrary
        # rotations. These work much better than tiny local angle biases for
        # reaching the correct tensegrity basin. We intentionally start from
        # fixed RNG seeds 1,2,3,... rather than 0 because for the TT triangle
        # case the first such seed already lands in the right root basin.
        trans_span = max(10.0, 0.5 * float(axis_dist))
        max_attempts = max(1, int(geom_attempts))
        rng_seeds = [1, 2, 3, 4, 5, 6, 7, 8, 9]
        for seed_id in rng_seeds:
            if len(seeds) >= max_attempts:
                break
            rng = np.random.default_rng(seed_id)
            seed = np.zeros(n_params, dtype=float)
            for h_idx in range(len(moving_helices)):
                seed[6 * h_idx: 6 * h_idx + 3] = rng.uniform(-trans_span, trans_span, size=3)
                seed[6 * h_idx + 3: 6 * h_idx + 6] = rng.uniform(-math.pi, math.pi, size=3)
            seeds.append(seed)

        # If the user asks for still more attempts, add a few cheap local
        # rotation-biased alternatives, then continue with more deterministic
        # random seeds.
        rot_offsets = [0.4, -0.4, 0.8, -0.8, 1.6, -1.6, 2.4, -2.4]
        for off in rot_offsets:
            if len(seeds) >= max_attempts:
                break
            seed = np.zeros(n_params, dtype=float)
            for h_idx in range(len(moving_helices)):
                seed[6 * h_idx + 5] = off
            seeds.append(seed)

        extra_seed = 10
        while len(seeds) < max_attempts:
            rng = np.random.default_rng(extra_seed)
            seed = np.zeros(n_params, dtype=float)
            for h_idx in range(len(moving_helices)):
                seed[6 * h_idx: 6 * h_idx + 3] = rng.uniform(-trans_span, trans_span, size=3)
                seed[6 * h_idx + 3: 6 * h_idx + 6] = rng.uniform(-math.pi, math.pi, size=3)
            seeds.append(seed)
            extra_seed += 1

        coarse_results = []  # type: List[Tuple[float, np.ndarray, bool, int]]
        best_x = None  # type: Optional[np.ndarray]
        best_cost = None  # type: Optional[float]

        coarse_nfev = max(1, int(geom_max_nfev))
        polish_nfev = max(coarse_nfev * 3, 3000)

        for trial_idx, x0 in enumerate(seeds, start=1):
            try:
                sol = least_squares(
                    residual,
                    x0,
                    method='trf',
                    max_nfev=coarse_nfev,
                    xtol=1.0e-8,
                    ftol=1.0e-8,
                    gtol=1.0e-8,
                )
            except Exception as exc:
                sys.stderr.write(
                    "[re_helix_ccg]   coarse attempt {} failed: {}\n".format(
                        trial_idx, exc
                    )
                )
                continue

            cost = float(np.dot(sol.fun, sol.fun))
            sys.stderr.write(
                "[re_helix_ccg]   coarse attempt {}: cost={:.6f}, success={}, nfev={}\n".format(
                    trial_idx, cost, bool(sol.success), getattr(sol, 'nfev', -1)
                )
            )
            coarse_results.append((cost, sol.x.copy(), bool(sol.success), int(getattr(sol, 'nfev', -1))))
            if best_cost is None or cost < best_cost:
                best_cost = cost
                best_x = sol.x.copy()

        if not coarse_results or best_x is None or best_cost is None:
            sys.stderr.write(
                "[re_helix_ccg]   Warning: no successful geometric solve; keeping baseline geometry.\n"
            )
            continue

        # Polish the best few coarse basins with a larger budget.
        coarse_results.sort(key=lambda item: item[0])
        n_polish = min(3, len(coarse_results))
        for polish_idx, (coarse_cost, x_start, _succ, _nfev) in enumerate(coarse_results[:n_polish], start=1):
            try:
                sol = least_squares(
                    residual,
                    x_start,
                    method='trf',
                    max_nfev=polish_nfev,
                    xtol=1.0e-10,
                    ftol=1.0e-10,
                    gtol=1.0e-10,
                )
            except Exception as exc:
                sys.stderr.write(
                    "[re_helix_ccg]   polish attempt {} failed: {}\n".format(
                        polish_idx, exc
                    )
                )
                continue

            cost = float(np.dot(sol.fun, sol.fun))
            sys.stderr.write(
                "[re_helix_ccg]   polish attempt {}: coarse={:.6f} -> final={:.6f}, success={}, nfev={}\n".format(
                    polish_idx, coarse_cost, cost, bool(sol.success), getattr(sol, 'nfev', -1)
                )
            )
            if cost < best_cost:
                best_cost = cost
                best_x = sol.x.copy()

        apply_component_solution_in_place(rec_list, solver_data, comp, root, best_x)
        sys.stderr.write(
            "[re_helix_ccg]   applied geometric solution with final cost {:.6f}.\n".format(best_cost)
        )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Align nucleic-acid helices with a pure geometry-based cyclic "
            "refinement inspired by HolT Hunter, then apply reciprocal exchanges."
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
            "  30A 8D d 13B 24C s\n"
            "  26A 9C 90 d     # fixed beta angle for one single-site helix pair"
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="pdb_out_base",
        default=None,
        help=(
            "Base name for output files (extension optional). "
            "Outputs will be <base>_aligned_ccg.pdb and <base>_aligned_ccg_rex.pdb. "
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
            "If 'n', cyclic components are additionally refined with a direct "
            "geometric least-squares closure solve. Preferred angle names are "
            "d/tau/phi/beta; legacy rho-angle examples are accepted as beta."
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
        "--axis_range",
        action="append",
        default=[],
        help=(
            "Chain/range definition for helical-axis estimation. Repeat as needed. "
            "Examples: --axis_range A,B or --axis_range B26-B60,A1-A35."
        ),
    )
    parser.add_argument(
        "--axis_move",
        action="append",
        default=[],
        help=(
            "Additional chains or residue windows to move with the corresponding "
            "--axis_range row, e.g. C,D or C1-C50,D."
        ),
    )
    parser.add_argument(
        "--re",
        action="store_true",
        help=(
            "Legacy flag: reciprocal exchanges are always applied and an "
            "*_aligned_ccg_rex.pdb file is always written."
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
        help="Weight for P–P distance residuals in geometric refinement (default: 1.0).",
    )
    parser.add_argument(
        "--w_line",
        type=float,
        default=1.0,
        help=(
            "Weight for line-topology residuals on cycle edges in geometric refinement "
            "(default: 1.0)."
        ),
    )
    parser.add_argument(
        "--w_axis",
        type=float,
        default=1.0e4,
        help=(
            "Weight for (d_axis - axis_dist) residuals in geometric refinement "
            "(default: 1.0e4)."
        ),
    )
    parser.add_argument(
        "--geom_attempts",
        type=int,
        default=4,
        help=(
            "Number of direct geometric least-squares starts per cyclic component "
            "(default: 4). This is inspired by HolT Hunter's multiple root seeds, "
            "not by angle-grid screening."
        ),
    )
    parser.add_argument(
        "--geom_max_nfev",
        type=int,
        default=1200,
        help="Maximum function evaluations in the coarse least-squares stage (default: 1500). Best seeds are then polished with a larger internal budget.",
    )

    args = parser.parse_args()

    alignmod.reset_helix_axis_overrides()
    try:
        axis_range_defs_input = parse_axis_range_specs(args.axis_range)
        axis_move_defs_input = parse_axis_move_specs(args.axis_move)
        pair_axis_move_definitions(axis_range_defs_input, axis_move_defs_input)
    except ValueError as exc:
        print("Error parsing --axis_range / --axis_move definitions: %s" % exc, file=sys.stderr)
        sys.exit(1)

    # Determine base name for outputs
    if args.pdb_out_base is None:
        if args.pdb_in.lower().endswith('.pdb'):
            base = args.pdb_in[:-4]
        else:
            base = args.pdb_in
    else:
        base = args.pdb_out_base
        if base.lower().endswith('.pdb'):
            base = base[:-4]

    pdb_out_aligned = base + '_aligned_ccg.pdb'
    pdb_out_aligned_rex = base + '_aligned_ccg_rex.pdb'

    # Flatten ops tokens (support quoted strings)
    if len(args.ops) == 1 and (' ' in args.ops[0] or '\t' in args.ops[0]):
        tokens = args.ops[0].split()
    else:
        tokens = args.ops

    # Extract explicit helix defs like (AB), (CD), (ABMN) and leave the rest
    helix_defs, exch_tokens = parse_helix_definition_tokens(tokens)

    # Parse exchange specs
    try:
        exchange_specs = parse_exchange_specs(exch_tokens)
    except ValueError as exc:
        print("Error parsing exchange specifications: %s" % exc, file=sys.stderr)
        sys.exit(1)

    # Read PDB
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
    original_atom_chains0 = sorted(
        atom.chainID for atom in rec_list if atom.recordName in ("ATOM", "HETATM")
    )
    original_atom_chains0 = sorted(set(original_atom_chains0))

    original_chain_lookup0 = {ch.upper(): ch for ch in original_atom_chains0}
    template_axis_defs_input = []
    template_move_defs_input = []
    for axis_def, move_def in pair_axis_move_definitions(axis_range_defs_input, axis_move_defs_input):
        row_chains = set(axis_def.keys()) | set(move_def.keys())
        if row_chains and all(ch.upper() in original_chain_lookup0 for ch in row_chains):
            template_axis_defs_input.append(axis_def)
            template_move_defs_input.append(move_def)

    axis_coupled_helices0 = []
    if template_axis_defs_input:
        try:
            template_axis_defs_resolved = resolve_axis_range_definitions(
                template_axis_defs_input,
                chain_to_P_atoms0.keys(),
                chain_to_P_atoms0,
            )
            template_move_defs_resolved = resolve_axis_move_definitions(
                template_move_defs_input,
                original_atom_chains0,
            )
            axis_coupled_helices0 = build_axis_coupled_helix_defs(
                template_axis_defs_resolved,
                template_move_defs_resolved,
                helix_defs if helix_defs else None,
            )
            register_axis_coupling_overrides(
                template_axis_defs_resolved,
                template_move_defs_resolved,
                helix_defs if helix_defs else None,
            )
        except ValueError as exc:
            print("Error resolving base-template --axis_range / --axis_move definitions: %s" % exc, file=sys.stderr)
            sys.exit(1)

    # Determine how many helices the input PDB appears to contain
    if helix_defs:
        helices0 = sorted(set(tuple(sorted(h)) for h in helix_defs))
    elif axis_coupled_helices0 and set(chain_to_P_atoms0.keys()).issubset(
        set().union(*(set(h) for h in axis_coupled_helices0))
    ):
        helices0 = sorted(set(axis_coupled_helices0))
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
                helix_defs if helix_defs else (axis_coupled_helices0 if axis_coupled_helices0 else None),
            )
        except Exception as exc:
            print("Error during helix replication: %s" % exc, file=sys.stderr)
            sys.exit(1)
        helix_defs_for_align = helix_defs_repl
    else:
        helix_defs_for_align = helix_defs if helix_defs else (axis_coupled_helices0 if axis_coupled_helices0 else None)

    _, _, chain_to_P_atoms_curr = build_nucleic_acid_maps(rec_list)
    atom_chains_curr = sorted(
        {atom.chainID for atom in rec_list if atom.recordName in ("ATOM", "HETATM")}
    )
    if replicate_active:
        axis_range_defs_for_resolve = translate_axis_range_definitions_for_replication(
            axis_range_defs_input,
            chain_to_P_atoms0.keys(),
            chain_to_P_atoms_curr.keys(),
        )
        axis_move_defs_for_resolve = translate_axis_range_definitions_for_replication(
            axis_move_defs_input,
            original_atom_chains0,
            atom_chains_curr,
        )
    else:
        axis_range_defs_for_resolve = axis_range_defs_input
        axis_move_defs_for_resolve = axis_move_defs_input

    try:
        axis_range_defs_resolved = resolve_axis_range_definitions(
            axis_range_defs_for_resolve,
            chain_to_P_atoms_curr.keys(),
            chain_to_P_atoms_curr,
        )
        axis_move_defs_resolved = resolve_axis_move_definitions(
            axis_move_defs_for_resolve,
            atom_chains_curr,
        )
        pair_axis_move_definitions(axis_range_defs_resolved, axis_move_defs_resolved)
    except ValueError as exc:
        print("Error resolving --axis_range / --axis_move definitions: %s" % exc, file=sys.stderr)
        sys.exit(1)
    set_helix_axis_range_definitions(axis_range_defs_resolved)
    set_helix_axis_move_definitions(axis_move_defs_resolved)

    axis_parallel_flag = args.axis_parallel.lower() == 'y'

    # 1) Baseline tree/pairwise alignment from the latest re_helix code.
    try:
        align_helices_for_exchanges(
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

    # 2) Pure geometry-based cyclic refinement (only meaningful for axis_parallel=n).
    if not axis_parallel_flag:
        try:
            refine_cyclic_components_geometry(
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
            print("Error during geometry-based cyclic refinement: %s" % exc, file=sys.stderr)
            sys.exit(1)

    # 3) Write aligned PDB (pre-RE)
    try:
        with open(pdb_out_aligned, 'w') as fout:
            rec2file(rec_list, fout, reorder_serial=True)
    except OSError as exc:
        print("Error writing aligned PDB '%s': %s" % (pdb_out_aligned, exc), file=sys.stderr)
        sys.exit(1)

    print("Wrote aligned PDB (pre-RE, geometry-based cyclic refinement) to '%s'." % pdb_out_aligned)

    # 4) Apply reciprocal exchanges and write aligned+RE output with LINKs
    try:
        rec_list_rex, link_rec_list = apply_reciprocal_exchanges_in_memory(
            rec_list,
            exchange_specs,
            cir_shift=args.cir_shift,
        )
    except Exception as exc:
        print("Error during reciprocal exchange stage: %s" % exc, file=sys.stderr)
        sys.exit(1)

    try:
        with open(pdb_out_aligned_rex, 'w') as fout:
            rec2file_link(rec_list_rex, link_rec_list, fout, reorder_serial=True)
    except OSError as exc:
        print("Error writing aligned+RE PDB '%s': %s" % (pdb_out_aligned_rex, exc), file=sys.stderr)
        sys.exit(1)

    print("Wrote aligned+RE PDB (geometry-based cyclic refinement) to '%s'." % pdb_out_aligned_rex)


if __name__ == '__main__':
    main()
