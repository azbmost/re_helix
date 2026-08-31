#!/usr/bin/env python3
"""
bend_helixV2_5.py

Bend a two-chain nucleic-acid helix at a user-selected phosphorus site.

Usage modes
-----------
1) Positional CLI arguments (backward-compatible):
       bend_helixV2_5.py input.pdb A36 0 30
       # tau defaults to 0 degrees in positional mode.
       # align defaults to y unless you set --align n.

2) Named CLI arguments:
       bend_helixV2_5.py --input A60-heli.pdb --pivot A36 --phi 0 --beta 30 --tau 10 --align y --origin y
       bend_helixV2_5.py --input A60-heli.pdb --pivot A36 --phi 0 --beta 30 --axis_range A1-A35,B60-B26 -o bent.pdb

3) GUI mode:
       bend_helixV2_5.py
       bend_helixV2_5.py --gui

Geometry implemented here
-------------------------
- The input duplex is treated as two paired strands.
- By default, the helix axis is estimated from paired P-atom midpoints.
- With --axis_range, a local helical axis can instead be estimated from user-defined
  residue windows, which is useful for already-bent helices. Use one range per
  line in the GUI or repeat --axis_range on the CLI, e.g. A1-A35,B60-B26.
- The chosen residue defines the *border* base pair between the two rigid pieces.
  Piece 2 contains the chosen base pair and every base pair after it along the
  duplex axis; piece 1 contains the remaining base pairs.
- At the axial height of the chosen P atom, we project that atom to the helix
  axis, creating a circle in the plane perpendicular to the axis.
- phi rotates the chosen P position around the axis on that circle.
- The hinge is the tangent line to that circle at the rotated point.
- beta rotates movable piece 2 rigidly about that hinge. Positive beta bends
  piece 2 away from the helix axis in the radial direction defined by phi.
- tau adds an extra twist of movable piece 2 about its own bent helical axis
  relative to fixed piece 1. Positive tau follows the right-hand rule about the bent axis;
  negative tau is left-handed.
- With --align y (default), piece 2 is translated after the bend/twist so that
  the pivot residue's P atom returns to its original pre-bend position.
- With --align n, the result matches the V2.1 bend/twist behaviour.
- With --sep y, piece 2 is additionally written under new chain IDs so that
  piece 1 and piece 2 are separated in the final PDB.
- Output filenames are written as *_PxByTz.pdb, or *_PxByTz_sep.pdb when
  --sep y is used, unless -o/--output is provided.
- With --origin y, an additional <main-output>-ori.pdb file is written that
  contains the original full helix and the same rigid full-helix transform used
  for piece 2, under sequential chain IDs.

Assumptions
-----------
- The input has exactly two chains containing P atoms. Straight helices work
  automatically; already-bent helices should use --axis_range to define a local axis.
- The two strands have the same number of P-bearing residues.
- Splitting is determined from axial pairing, so selecting either residue of the
  same base pair gives the same split.
"""

import argparse
import itertools
import math
import os
import re
import shlex
import string
import sys
from collections import OrderedDict
from dataclasses import dataclass, field
from typing import Callable, Dict, FrozenSet, Iterable, List, Mapping, Optional, Sequence, Tuple

try:
    from re_helix_lib.gui_icon import apply_optional_icon
except ImportError:  # pragma: no cover - direct script execution fallback
    from gui_icon import apply_optional_icon

EPS = 1.0e-8
CHAIN_ID_CANDIDATES = string.ascii_uppercase + string.ascii_lowercase + string.digits
Point3D = Tuple[float, float, float]
TOOL_NAME = "Bend Helix"
VERSION = "V2.5"
APP_TITLE = f"{TOOL_NAME} {VERSION}"
SCREEN_REFINEMENT_TOLERANCE_DEG = 1.0e-3
SCREEN_REFINEMENT_SEED_LIMIT = 8
SCREEN_REFINEMENT_MAX_ROUNDS = 80
DEFAULT_SCREEN_STEP_DEG = 6.0
DEFAULT_SCREEN_RANGES: Mapping[str, Tuple[float, float]] = {
    "phi": (-90.0, 90.0),
    "beta": (-180.0, 180.0),
    "tau": (-180.0, 180.0),
}

SCREENING_GUI_HELP: Mapping[str, str] = {
    "grid_from": (
        "From angle (degrees)\n\n"
        "The first coarse-grid angle value, in degrees. From may be smaller "
        "or larger than To, so both ascending and descending grids are allowed.\n\n"
        "Example: From = -30, To = 30, Step = 10 starts at -30 degrees."
    ),
    "grid_to": (
        "To angle (degrees)\n\n"
        "The final coarse-grid angle value, in degrees. To is always included, "
        "even when Step does not land on it exactly.\n\n"
        "Example: From = 0, To = 10, Step = 4 tests 0, 4, 8, and 10 degrees."
    ),
    "grid_step": (
        "Step (degrees)\n\n"
        "A positive coarse-grid spacing in degrees. Bend Helix first tests the inclusive "
        "grid, then efficiently searches between nearby grid values around the strongest "
        "coarse candidates. The adaptive refinement halves its spacing until it reaches "
        "0.001 degree. The default coarse Step is 6 degrees. With two screened angles, "
        "the coarse search uses their Cartesian "
        "product and refinement adjusts both angles.\n\n"
        "Example: From = 0, To = 10, Step = 4 starts with 0, 4, 8, and 10, then can "
        "select an in-between value such as 6.35 degrees."
    ),
    "mode": (
        "Screening mode\n\n"
        "Screening for distance finds the coarse or refined candidate whose two endpoints "
        "are closest to a requested separation. Screening for rotation finds the candidate whose "
        "signed endpoint-to-endpoint angle around an axis is closest to a requested angle."
    ),
    "target": (
        "Target value\n\n"
        "For distance mode, enter a nonnegative distance in angstroms. For rotation mode, "
        "enter a signed angle in degrees. The target need not be achievable: Bend Helix "
        "keeps the tested coarse or refined candidate with the smallest error.\n\n"
        "Examples: 12.5 angstroms; -45 degrees."
    ),
    "endpoint1_atom": (
        "Endpoint 1 overlay atom\n\n"
        "Enter an atom from the origin-overlay PDB as CHAIN:RESIDUE:ATOM. The chain ID, "
        "residue number, and atom name must match that overlay, not the input PDB's chain "
        "namespace. Endpoint order determines the sign in rotation mode.\n\n"
        "Examples: A:36:P or C:36:O5'."
    ),
    "endpoint2_source": (
        "Endpoint 2 source\n\n"
        "Choose Overlay atom to use another atom in the origin overlay, or XYZ point to "
        "use a fixed coordinate. Rotation mode also allows Phi-corrected pivot P, whose "
        "position is recomputed for every candidate phi value."
    ),
    "endpoint2_atom": (
        "Endpoint 2 overlay atom\n\n"
        "Enter the second origin-overlay atom as CHAIN:RESIDUE:ATOM. Its identity is "
        "resolved in the candidate model on every coarse or refinement test.\n\n"
        "Examples: B:24:P or D:24:O3'."
    ),
    "endpoint2_xyz": (
        "Endpoint 2 XYZ\n\n"
        "Enter a fixed Cartesian point as x, y, and z coordinates in angstroms. The point "
        "does not move as the screened model changes.\n\n"
        "Example: x = 15.0, y = -2.5, z = 8.0."
    ),
    "endpoint2_pivot": (
        "Phi-corrected pivot P\n\n"
        "Use the selected pivot residue's P position after applying the current candidate's "
        "phi correction. This endpoint is candidate-dependent and is available only for "
        "rotation screening."
    ),
    "axis_source": (
        "Rotation axis source\n\n"
        "Geometric definition lets you provide an axis point and direction below. Local "
        "axis range(s) uses the helix axis derived from the range selected in the main "
        "Bend Helix window."
    ),
    "axis_point_source": (
        "Axis point source\n\n"
        "Choose XYZ point for a fixed point on the rotation axis, or Overlay atom to put "
        "the axis through an origin-overlay atom."
    ),
    "axis_point_xyz": (
        "Axis point XYZ\n\n"
        "Enter any fixed Cartesian point on the rotation axis in angstroms.\n\n"
        "Example: x = 0, y = 0, z = 12.5."
    ),
    "axis_point_atom": (
        "Axis point overlay atom\n\n"
        "Enter an origin-overlay atom as CHAIN:RESIDUE:ATOM. The candidate-dependent "
        "position of that atom is used as a point on the rotation axis.\n\n"
        "Example: A:36:P."
    ),
    "axis_vector_source": (
        "Axis vector source\n\n"
        "Choose how to define the positive axis direction: a direct vector, the direction "
        "from point 1 to point 2, the direction from overlay atom 1 to overlay atom 2, or "
        "the right-hand normal (vector 1 crossed with vector 2)."
    ),
    "direct_vector": (
        "Direct axis vector\n\n"
        "Enter the x, y, and z components of a nonzero vector. Its magnitude does not "
        "matter; its sign sets the positive direction used by the right-hand rule.\n\n"
        "Example: x = 0, y = 0, z = 1."
    ),
    "two_xyz_points": (
        "Axis direction from two XYZ points\n\n"
        "Enter two distinct fixed Cartesian points. The positive axis direction runs from "
        "Point 1 toward Point 2.\n\n"
        "Example: Point 1 = (0, 0, 0), Point 2 = (0, 0, 10)."
    ),
    "two_overlay_atoms": (
        "Axis direction from two overlay atoms\n\n"
        "Enter two distinct origin-overlay atoms as CHAIN:RESIDUE:ATOM. The positive axis "
        "direction runs from Overlay atom 1 toward Overlay atom 2.\n\n"
        "Example: A:10:P to A:30:P."
    ),
    "normal_vectors": (
        "Axis normal to two vectors\n\n"
        "Enter two nonparallel vectors. The positive axis direction is vector 1 crossed "
        "with vector 2, following the right-hand rule.\n\n"
        "Example: vector 1 = (1, 0, 0), vector 2 = (0, 1, 0) gives +z."
    ),
}


@dataclass
class AtomRecord:
    line: str
    record_name: str
    atom_name: str
    chain_id: str
    res_seq: int
    x: float
    y: float
    z: float

    def coord(self) -> Tuple[float, float, float]:
        return (self.x, self.y, self.z)

    def set_coord(self, xyz: Tuple[float, float, float]) -> None:
        self.x, self.y, self.z = xyz

    def set_chain_id(self, new_chain_id: str) -> None:
        if len(new_chain_id) != 1:
            raise ValueError(f"Chain ID must be exactly one character, got {new_chain_id!r}.")
        self.chain_id = new_chain_id

    def to_line(self) -> str:
        return (
            f"{self.line[:21]}{self.chain_id}{self.line[22:30]}"
            f"{self.x:8.3f}{self.y:8.3f}{self.z:8.3f}{self.line[54:]}"
        )


@dataclass
class RawRecord:
    line: str


@dataclass
class Residue:
    chain_id: str
    res_seq: int
    first_seen: int
    atoms: List[AtomRecord] = field(default_factory=list)
    p_atom: AtomRecord = None  # type: ignore[assignment]

    def p_coord(self) -> Tuple[float, float, float]:
        if self.p_atom is None:
            raise ValueError(f"Residue {self.chain_id}{self.res_seq} does not contain a P atom.")
        return self.p_atom.coord()


@dataclass
class AxisRangeTerm:
    chain_id: str
    start_res: int
    end_res: int
    lo: int
    hi: int


@dataclass
class AxisRangeDefinition:
    original_text: str
    terms: List[AxisRangeTerm]

    def ranges(self) -> Dict[str, Tuple[int, int]]:
        return {term.chain_id: (term.lo, term.hi) for term in self.terms}


@dataclass(frozen=True)
class BendGeometryPreparation:
    """Immutable input geometry shared by all bend-angle candidates."""

    selected_key: Tuple[str, int]
    selected_p: Point3D
    axis_point: Point3D
    axis_dir: Point3D
    axis_foot: Point3D
    radial: Point3D
    radius: float
    pair_idx: int
    pair_keys: Tuple[Tuple[Tuple[str, int], Tuple[str, int]], ...]
    piece2_keys: FrozenSet[Tuple[str, int]]
    duplex_chains: Tuple[str, ...]
    axis_source: str
    axis_range_used: Optional[str]


@dataclass(frozen=True)
class BendTransform:
    """Pure rigid transform derived from one phi/beta/tau candidate."""

    phi_deg: float
    beta_deg: float
    tau_deg: float
    hinge_point: Point3D
    hinge_dir: Point3D
    twist_axis_point_pre_align: Point3D
    twist_axis_point: Point3D
    twist_axis_dir: Point3D
    align_mode: str
    align_translation: Point3D

    def transform_coord(self, coord: Point3D) -> Point3D:
        """Apply the full-helix transform, including optional pivot alignment."""
        beta_rad = math.radians(self.beta_deg)
        tau_rad = math.radians(self.tau_deg)
        transformed = rotate_point_about_line(
            coord, self.hinge_point, self.hinge_dir, beta_rad
        )
        if abs(tau_rad) > 0.0:
            transformed = rotate_point_about_line(
                transformed,
                self.twist_axis_point_pre_align,
                self.twist_axis_dir,
                tau_rad,
            )
        if self.align_mode == "y":
            transformed = v_add(transformed, self.align_translation)
        return transformed


@dataclass(frozen=True)
class OverlayAtomSelector:
    """Atom identity expressed using origin-overlay chain IDs."""

    chain_id: str
    res_seq: int
    atom_name: str


@dataclass(frozen=True)
class ScreeningPoint:
    """A screening point: overlay atom, XYZ point, or phi-corrected pivot."""

    kind: str
    value: object = None


@dataclass(frozen=True)
class ScreeningAxis:
    """Local or geometrically defined axis for a rotation screening target."""

    source: str
    point: Optional[ScreeningPoint] = None
    vector_source: str = "direct_vector"
    vector: object = None
    point1: Optional[ScreeningPoint] = None
    point2: Optional[ScreeningPoint] = None
    normal_vector1: object = None
    normal_vector2: object = None


@dataclass(frozen=True)
class ScreeningRequest:
    """Distance or signed-rotation value to approach during angle screening."""

    mode: str
    target: float
    point1: ScreeningPoint
    point2: ScreeningPoint
    axis: Optional[ScreeningAxis] = None


@dataclass(frozen=True)
class ScreenAngleRange:
    """Inclusive grid bounds for one screened angle."""

    name: str
    start: float
    stop: float
    step: float


@dataclass(frozen=True)
class ScreeningAtom:
    """Immutable source-PDB atom snapshot used by the screening engine."""

    source_chain_id: str
    res_seq: int
    atom_name: str
    coord: Point3D


@dataclass(frozen=True)
class ScreeningContext:
    """Immutable bend preparation and origin-overlay atom namespace."""

    input_pdb: str
    preparation: BendGeometryPreparation
    atoms: Tuple[ScreeningAtom, ...]
    model1_chain_map_items: Tuple[Tuple[str, str], ...]
    model2_chain_map_items: Tuple[Tuple[str, str], ...]

    @property
    def origin_chain_map_model1(self) -> Dict[str, str]:
        return dict(self.model1_chain_map_items)

    @property
    def origin_chain_map_model2(self) -> Dict[str, str]:
        return dict(self.model2_chain_map_items)

    @property
    def overlay_chain_map(self) -> Dict[str, Tuple[str, bool]]:
        """Map overlay chain ID to ``(source_chain_id, is_transformed)``."""
        result = {
            overlay: (source, False)
            for source, overlay in self.model1_chain_map_items
        }
        result.update(
            {
                overlay: (source, True)
                for source, overlay in self.model2_chain_map_items
            }
        )
        return result


@dataclass(frozen=True)
class ScreeningResult:
    """Best coarse-to-fine candidate and its achieved target metric."""

    phi_deg: float
    beta_deg: float
    tau_deg: float
    achieved_value: float
    error: float
    candidate_count: int
    refinement_candidate_count: int = 0

    @property
    def angles(self) -> Dict[str, float]:
        return {
            "phi": self.phi_deg,
            "beta": self.beta_deg,
            "tau": self.tau_deg,
        }

    @property
    def observed_value(self) -> float:
        """Alias for callers that describe the achieved metric as observed."""
        return self.achieved_value

    @property
    def target_error(self) -> float:
        return self.error

    @property
    def evaluated_candidate_count(self) -> int:
        """Total unique coarse and refinement candidates evaluated."""
        return self.candidate_count + self.refinement_candidate_count


# ---------------------------------------------------------------------------
# Small vector helpers
# ---------------------------------------------------------------------------


def v_add(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> Tuple[float, float, float]:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])



def v_sub(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> Tuple[float, float, float]:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])



def v_scale(a: Tuple[float, float, float], s: float) -> Tuple[float, float, float]:
    return (a[0] * s, a[1] * s, a[2] * s)



def v_dot(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]



def v_cross(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> Tuple[float, float, float]:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )



def v_len(a: Tuple[float, float, float]) -> float:
    return math.sqrt(v_dot(a, a))



def v_norm(a: Tuple[float, float, float]) -> Tuple[float, float, float]:
    l = v_len(a)
    if l < EPS:
        raise ValueError("Encountered a near-zero vector during normalisation.")
    return (a[0] / l, a[1] / l, a[2] / l)



def centroid(points: List[Tuple[float, float, float]]) -> Tuple[float, float, float]:
    if not points:
        raise ValueError("Cannot compute centroid of an empty point set.")
    sx = sy = sz = 0.0
    for x, y, z in points:
        sx += x
        sy += y
        sz += z
    n = float(len(points))
    return (sx / n, sy / n, sz / n)



def mat_vec_mul(m: List[List[float]], v: Tuple[float, float, float]) -> Tuple[float, float, float]:
    return (
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    )



def principal_axis(points: List[Tuple[float, float, float]]) -> Tuple[Tuple[float, float, float], Tuple[float, float, float]]:
    """Return (centroid, unit principal axis) from a 3x3 covariance matrix."""
    if len(points) < 2:
        raise ValueError("At least two points are required to estimate a line axis.")

    c = centroid(points)
    cov = [[0.0, 0.0, 0.0] for _ in range(3)]
    for p in points:
        d = v_sub(p, c)
        cov[0][0] += d[0] * d[0]
        cov[0][1] += d[0] * d[1]
        cov[0][2] += d[0] * d[2]
        cov[1][0] += d[1] * d[0]
        cov[1][1] += d[1] * d[1]
        cov[1][2] += d[1] * d[2]
        cov[2][0] += d[2] * d[0]
        cov[2][1] += d[2] * d[1]
        cov[2][2] += d[2] * d[2]

    v = v_norm((1.0, 1.0, 1.0))
    for _ in range(64):
        mv = mat_vec_mul(cov, v)
        if v_len(mv) < EPS:
            break
        v = v_norm(mv)

    if v_len(v) < EPS:
        fallback = v_sub(points[-1], points[0])
        if v_len(fallback) < EPS:
            raise ValueError("Failed to estimate a stable helix axis.")
        v = v_norm(fallback)

    return c, v



def rotate_vector(vec: Tuple[float, float, float], axis: Tuple[float, float, float], angle_rad: float) -> Tuple[float, float, float]:
    """Rodrigues rotation of a vector around a unit axis through the origin."""
    k = v_norm(axis)
    c = math.cos(angle_rad)
    s = math.sin(angle_rad)
    term1 = v_scale(vec, c)
    term2 = v_scale(v_cross(k, vec), s)
    term3 = v_scale(k, v_dot(k, vec) * (1.0 - c))
    return v_add(v_add(term1, term2), term3)



def rotate_point_about_line(
    point: Tuple[float, float, float],
    line_point: Tuple[float, float, float],
    line_dir: Tuple[float, float, float],
    angle_rad: float,
) -> Tuple[float, float, float]:
    rel = v_sub(point, line_point)
    rel_rot = rotate_vector(rel, line_dir, angle_rad)
    return v_add(line_point, rel_rot)



def project_point_to_line(
    point: Tuple[float, float, float],
    line_point: Tuple[float, float, float],
    line_dir: Tuple[float, float, float],
) -> Tuple[float, float, float]:
    u = v_norm(line_dir)
    t = v_dot(v_sub(point, line_point), u)
    return v_add(line_point, v_scale(u, t))


# ---------------------------------------------------------------------------
# PDB parsing / formatting
# ---------------------------------------------------------------------------


def parse_atom_record(line: str) -> AtomRecord:
    return AtomRecord(
        line=line,
        record_name=line[:6].strip(),
        atom_name=line[12:16].strip(),
        chain_id=line[21],
        res_seq=int(line[22:26]),
        x=float(line[30:38]),
        y=float(line[38:46]),
        z=float(line[46:54]),
    )



def read_pdb(path: str):
    records = []
    residues: "OrderedDict[Tuple[str, int], Residue]" = OrderedDict()

    with open(path, "r") as handle:
        for line_no, line in enumerate(handle):
            rec_name = line[:6].strip()
            if rec_name in ("ATOM", "HETATM") and len(line) >= 54:
                atom = parse_atom_record(line)
                records.append(atom)
                key = (atom.chain_id, atom.res_seq)
                if key not in residues:
                    residues[key] = Residue(chain_id=atom.chain_id, res_seq=atom.res_seq, first_seen=line_no)
                residues[key].atoms.append(atom)
                if atom.atom_name == "P" and residues[key].p_atom is None:
                    residues[key].p_atom = atom
            else:
                records.append(RawRecord(line=line))

    return records, residues



def parse_residue_token(token: str) -> Tuple[str, int]:
    t = token.strip().replace(".", "")
    if len(t) < 2:
        raise ValueError(f"Invalid residue token '{token}'.")

    if t[0].isalpha() and t[1:].isdigit():
        return (t[0], int(t[1:]))
    if t[-1].isalpha() and t[:-1].isdigit():
        return (t[-1], int(t[:-1]))

    raise ValueError(
        f"Invalid residue token '{token}': expected A36, 36A, A.36, or 36.A."
    )



def resolve_selected_key(raw_key: Tuple[str, int], residues: "OrderedDict[Tuple[str, int], Residue]") -> Tuple[str, int]:
    if raw_key in residues:
        return raw_key

    chain_id, res_seq = raw_key
    alternatives = [(chain_id.swapcase(), res_seq), (chain_id.upper(), res_seq), (chain_id.lower(), res_seq)]
    for alt in alternatives:
        if alt in residues:
            return alt

    raise ValueError(f"Residue {chain_id}{res_seq} was not found in the input PDB.")


def parse_axis_range_term(term: str) -> AxisRangeTerm:
    token = term.strip().replace(".", "")
    if not token:
        raise ValueError("Empty axis-range term.")
    if "-" not in token:
        raise ValueError(
            f"Invalid axis-range term '{term}': expected a residue range like 'A1-A35'."
        )

    left, right = token.split("-", 1)
    chain1, res1 = parse_residue_token(left)

    right = right.strip()
    if not right:
        raise ValueError(f"Invalid axis-range term '{term}': missing end residue after '-'.")

    if any(ch.isalpha() for ch in right):
        chain2, res2 = parse_residue_token(right)
    else:
        chain2, res2 = chain1, int(right)

    if chain1.upper() != chain2.upper():
        raise ValueError(
            f"Axis-range term '{term}' must stay on one chain; got '{chain1}' and '{chain2}'."
        )

    return AxisRangeTerm(
        chain_id=chain1,
        start_res=res1,
        end_res=res2,
        lo=min(res1, res2),
        hi=max(res1, res2),
    )


def parse_axis_range_spec(spec: str) -> AxisRangeDefinition:
    text = spec.strip()
    parts = [part.strip() for part in text.split(",") if part.strip()]
    if not parts:
        raise ValueError("Axis-range specification cannot be empty.")

    terms: List[AxisRangeTerm] = []
    seen: Dict[str, str] = {}
    for part in parts:
        term = parse_axis_range_term(part)
        key = term.chain_id.upper()
        if key in seen:
            raise ValueError(
                f"Axis-range specification '{spec}' defines chain '{term.chain_id}' more than once."
            )
        seen[key] = term.chain_id
        terms.append(term)

    return AxisRangeDefinition(original_text=text, terms=terms)


def parse_axis_range_specs(specs: Optional[Iterable[str]]) -> List[AxisRangeDefinition]:
    parsed: List[AxisRangeDefinition] = []
    if not specs:
        return parsed
    for spec in specs:
        spec = str(spec).strip()
        if not spec:
            continue
        # GUI users can put one --axis_range per line; semicolons are also accepted.
        for chunk in spec.replace(";", "\n").splitlines():
            chunk = chunk.strip()
            if chunk:
                parsed.append(parse_axis_range_spec(chunk))
    return parsed


def format_axis_range_spec(axis_def: AxisRangeDefinition) -> str:
    parts = []
    for term in axis_def.terms:
        parts.append(f"{term.chain_id}{term.start_res}-{term.chain_id}{term.end_res}")
    return ",".join(parts)


def resolve_axis_range_definitions(
    defs: List[AxisRangeDefinition],
    available_chains: Iterable[str],
) -> List[AxisRangeDefinition]:
    lookup = {str(ch).upper(): str(ch) for ch in available_chains}
    resolved: List[AxisRangeDefinition] = []

    for axis_def in defs:
        terms: List[AxisRangeTerm] = []
        seen: Dict[str, str] = {}
        for term in axis_def.terms:
            key = term.chain_id.upper()
            if key not in lookup:
                raise ValueError(
                    f"Axis-range definition {format_axis_range_spec(axis_def)} refers to chain "
                    f"'{term.chain_id}', but that chain is not present in the input PDB."
                )
            actual_chain = lookup[key]
            if actual_chain in seen:
                raise ValueError(
                    f"Axis-range definition {format_axis_range_spec(axis_def)} resolves chain "
                    f"'{term.chain_id}' to '{actual_chain}' more than once."
                )
            seen[actual_chain] = actual_chain
            terms.append(
                AxisRangeTerm(
                    chain_id=actual_chain,
                    start_res=term.start_res,
                    end_res=term.end_res,
                    lo=min(term.start_res, term.end_res),
                    hi=max(term.start_res, term.end_res),
                )
            )
        resolved.append(AxisRangeDefinition(original_text=axis_def.original_text, terms=terms))

    return resolved


def axis_range_contains_residue(axis_def: AxisRangeDefinition, key: Tuple[str, int]) -> bool:
    chain_id, res_seq = key
    for term in axis_def.terms:
        if term.chain_id == chain_id and term.lo <= res_seq <= term.hi:
            return True
    return False


def selected_axis_range_definition(
    defs: List[AxisRangeDefinition],
    selected_key: Tuple[str, int],
) -> Optional[AxisRangeDefinition]:
    if not defs:
        return None

    matching = [axis_def for axis_def in defs if axis_range_contains_residue(axis_def, selected_key)]
    if len(matching) == 1:
        return matching[0]
    if len(matching) > 1:
        texts = "; ".join(format_axis_range_spec(axis_def) for axis_def in matching)
        raise ValueError(
            f"Multiple --axis_range definitions contain pivot residue {selected_key[0]}{selected_key[1]}: "
            f"{texts}. Please make the ranges non-overlapping or provide only one matching range."
        )

    if len(defs) == 1:
        return defs[0]

    raise ValueError(
        "More than one --axis_range was supplied, but none contains the pivot residue "
        f"{selected_key[0]}{selected_key[1]}. Include the pivot in the intended local-axis range, "
        "or supply only one --axis_range to force that local axis."
    )


def collect_axis_range_residues(
    chain_to_residues: Dict[str, List[Residue]],
    axis_def: AxisRangeDefinition,
) -> Dict[str, List[Residue]]:
    range_map = axis_def.ranges()
    selected: Dict[str, List[Residue]] = OrderedDict()
    for chain_id, residues_for_chain in chain_to_residues.items():
        if chain_id not in range_map:
            continue
        lo, hi = range_map[chain_id]
        curr = [res for res in residues_for_chain if lo <= res.res_seq <= hi and res.p_atom is not None]
        if not curr:
            raise ValueError(
                f"Axis-range definition {format_axis_range_spec(axis_def)} selects no P atoms on chain "
                f"'{chain_id}' in residue range {lo}-{hi}."
            )
        selected[chain_id] = curr
    return selected


def orientation_vector_from_axis_range(
    axis_def: AxisRangeDefinition,
    residues: "OrderedDict[Tuple[str, int], Residue]",
    range_residues: Dict[str, List[Residue]],
) -> Optional[Tuple[float, float, float]]:
    for term in axis_def.terms:
        start = residues.get((term.chain_id, term.start_res))
        end = residues.get((term.chain_id, term.end_res))
        if start is not None and end is not None and start.p_atom is not None and end.p_atom is not None:
            vec = v_sub(end.p_coord(), start.p_coord())
            if v_len(vec) >= EPS:
                return vec

    for term in axis_def.terms:
        curr = list(range_residues.get(term.chain_id, []))
        if len(curr) < 2:
            continue
        curr.sort(key=lambda r: (r.res_seq, r.first_seen), reverse=(term.end_res < term.start_res))
        vec = v_sub(curr[-1].p_coord(), curr[0].p_coord())
        if v_len(vec) >= EPS:
            return vec

    return None


def estimate_axis_from_manual_range(
    residues: "OrderedDict[Tuple[str, int], Residue]",
    chain_to_residues: Dict[str, List[Residue]],
    chains: List[str],
    axis_def: AxisRangeDefinition,
) -> Tuple[Tuple[float, float, float], Tuple[float, float, float]]:
    range_residues = collect_axis_range_residues(chain_to_residues, axis_def)
    coords = [res.p_coord() for chain_res in range_residues.values() for res in chain_res]
    if len(coords) < 2:
        raise ValueError(
            f"Axis-range definition {format_axis_range_spec(axis_def)} selects fewer than two P atoms."
        )

    axis_point, axis_dir = principal_axis(coords)
    orientation_vec = orientation_vector_from_axis_range(axis_def, residues, range_residues)
    if orientation_vec is not None and v_dot(orientation_vec, axis_dir) < 0.0:
        axis_dir = v_scale(axis_dir, -1.0)

    # If both duplex chains are covered with matching numbers of P atoms, refine
    # the axis from paired P-atom midpoints, preserving the user-defined direction.
    if all(chain in range_residues for chain in chains):
        n0 = len(range_residues[chains[0]])
        n1 = len(range_residues[chains[1]])
        if n0 == n1 and n0 >= 2:
            range_pair_list = pair_residues_by_axis(range_residues, chains, axis_point, axis_dir)
            axis_point_mid, axis_dir_mid = principal_axis([
                v_scale(v_add(r1.p_coord(), r2.p_coord()), 0.5) for r1, r2 in range_pair_list
            ])
            if orientation_vec is not None and v_dot(orientation_vec, axis_dir_mid) < 0.0:
                axis_dir_mid = v_scale(axis_dir_mid, -1.0)
            axis_point, axis_dir = axis_point_mid, axis_dir_mid

    return axis_point, axis_dir



def normalize_sep(value: str) -> str:
    s = value.strip().lower()
    if s not in ("y", "n"):
        raise ValueError(f"Invalid --sep value '{value}': expected 'y' or 'n'.")
    return s



def normalize_align(value: str) -> str:
    s = value.strip().lower()
    if s not in ("y", "n"):
        raise ValueError(f"Invalid --align value '{value}': expected 'y' or 'n'.")
    return s



def normalize_origin(value: str) -> str:
    s = value.strip().lower()
    if s not in ("y", "n"):
        raise ValueError(f"Invalid --origin value '{value}': expected 'y' or 'n'.")
    return s



def format_angle_for_filename(value: float) -> str:
    if abs(value - round(value)) < 1.0e-8:
        s = str(int(round(value)))
    else:
        s = f"{value:.3f}".rstrip("0").rstrip(".")
    s = s.replace("-", "m").replace(".", "p")
    return s



def format_float_for_cli(value: float) -> str:
    if abs(value - round(value)) < 1.0e-8:
        return str(int(round(value)))
    return f"{value:.15g}"


def build_equivalent_cli_command(
    input_pdb: str,
    pivot_residue: str,
    phi_deg: float,
    beta_deg: float,
    tau_deg: float,
    sep_mode: str,
    align_mode: str,
    origin_mode: str,
    output_pdb: Optional[str] = None,
    axis_range_specs: Optional[Iterable[str]] = None,
) -> str:
    script_name = os.path.basename(__file__) if "__file__" in globals() else "bend_helixV2_5.py"
    parts = [
        "python",
        script_name,
        "--input",
        input_pdb,
        "--pivot",
        pivot_residue,
        "--phi",
        format_float_for_cli(phi_deg),
        "--beta",
        format_float_for_cli(beta_deg),
        "--tau",
        format_float_for_cli(tau_deg),
        "--sep",
        sep_mode,
        "--align",
        align_mode,
        "--origin",
        origin_mode,
    ]
    if output_pdb:
        parts.extend(["-o", output_pdb])
    for spec in axis_range_specs or []:
        spec = str(spec).strip()
        if spec:
            parts.extend(["--axis_range", spec])
    return " ".join(shlex.quote(part) for part in parts)

def make_output_name(
    inp: str,
    phi_deg: float,
    beta_deg: float,
    tau_deg: float,
    sep_mode: str = "n",
    screen_mode: bool = False,
) -> str:
    stem, ext = os.path.splitext(inp)
    if not ext:
        ext = ".pdb"
    suffix = (
        f"_P{format_angle_for_filename(phi_deg)}"
        f"B{format_angle_for_filename(beta_deg)}"
        f"T{format_angle_for_filename(tau_deg)}"
    )
    if screen_mode:
        suffix += "_scr"
    if sep_mode == "y":
        suffix += "_sep"
    return f"{stem}{suffix}{ext}"



def make_origin_output_name(main_out_path: str) -> str:
    stem, ext = os.path.splitext(main_out_path)
    if not ext:
        ext = ".pdb"
    return f"{stem}-ori{ext}"


def normalize_output_path(output_pdb: Optional[str]) -> Optional[str]:
    if output_pdb is None:
        return None
    out = output_pdb.strip()
    if not out:
        return None
    stem, ext = os.path.splitext(out)
    if not ext:
        out = stem + ".pdb"
    return out


# ---------------------------------------------------------------------------
# Duplex construction
# ---------------------------------------------------------------------------


def build_duplex_data(
    residues: "OrderedDict[Tuple[str, int], Residue]",
    selected_key: Optional[Tuple[str, int]] = None,
    axis_range_defs: Optional[List[AxisRangeDefinition]] = None,
):
    chain_to_residues: Dict[str, List[Residue]] = OrderedDict()
    for residue in residues.values():
        if residue.p_atom is None:
            continue
        chain_to_residues.setdefault(residue.chain_id, []).append(residue)

    if len(chain_to_residues) != 2:
        chains = ", ".join(repr(c) for c in chain_to_residues.keys())
        raise ValueError(
            f"{APP_TITLE} currently expects exactly two chains that contain P atoms; "
            f"found {len(chain_to_residues)} ({chains})."
        )

    chains = list(chain_to_residues.keys())
    ref_chain = chains[0]
    ref_res_by_seq = sorted(chain_to_residues[ref_chain], key=lambda r: (r.res_seq, r.first_seen))
    if len(ref_res_by_seq) < 2:
        raise ValueError(f"Reference chain {ref_chain!r} does not have enough P atoms.")

    axis_range_used = None
    if axis_range_defs:
        if selected_key is None:
            raise ValueError("Internal error: selected_key is required when --axis_range is used.")
        axis_range_used = selected_axis_range_definition(axis_range_defs, selected_key)

    if axis_range_used is None:
        # Provisional axis from all P atoms.
        all_p = [res.p_coord() for chain in chains for res in chain_to_residues[chain]]
        axis_point, axis_dir = principal_axis(all_p)

        # Orient the axis so that the first chain encountered in the file increases
        # in projection as its residue numbers increase.
        if v_dot(v_sub(ref_res_by_seq[-1].p_coord(), ref_res_by_seq[0].p_coord()), axis_dir) < 0.0:
            axis_dir = v_scale(axis_dir, -1.0)

        # Pair residues by axial order, refine the axis from paired P-midpoints, then pair again.
        pair_list = pair_residues_by_axis(chain_to_residues, chains, axis_point, axis_dir)
        axis_point, axis_dir = principal_axis([
            v_scale(v_add(r1.p_coord(), r2.p_coord()), 0.5) for r1, r2 in pair_list
        ])
        if v_dot(v_sub(ref_res_by_seq[-1].p_coord(), ref_res_by_seq[0].p_coord()), axis_dir) < 0.0:
            axis_dir = v_scale(axis_dir, -1.0)
        pair_list = pair_residues_by_axis(chain_to_residues, chains, axis_point, axis_dir)
        axis_source = "automatic whole-duplex axis"
    else:
        axis_point, axis_dir = estimate_axis_from_manual_range(
            residues=residues,
            chain_to_residues=chain_to_residues,
            chains=chains,
            axis_def=axis_range_used,
        )
        pair_list = pair_residues_by_axis(chain_to_residues, chains, axis_point, axis_dir)
        axis_source = "manual --axis_range " + format_axis_range_spec(axis_range_used)

    return {
        "chains": chains,
        "ref_chain": ref_chain,
        "chain_to_residues": chain_to_residues,
        "axis_point": axis_point,
        "axis_dir": axis_dir,
        "pairs": pair_list,
        "axis_range_used": axis_range_used,
        "axis_source": axis_source,
    }


def pair_residues_by_axis(
    chain_to_residues: Dict[str, List[Residue]],
    chains: List[str],
    axis_point: Tuple[float, float, float],
    axis_dir: Tuple[float, float, float],
) -> List[Tuple[Residue, Residue]]:
    ordered = []
    for chain in chains:
        curr = sorted(
            chain_to_residues[chain],
            key=lambda r: v_dot(v_sub(r.p_coord(), axis_point), axis_dir),
        )
        ordered.append(curr)

    n0 = len(ordered[0])
    n1 = len(ordered[1])
    if n0 != n1:
        raise ValueError(
            "The two chains do not have the same number of P-bearing residues "
            f"({chains[0]}: {n0}, {chains[1]}: {n1}); axial base-pairing would be ambiguous."
        )

    return list(zip(ordered[0], ordered[1]))


# ---------------------------------------------------------------------------
# Bending logic
# ---------------------------------------------------------------------------


def find_pair_index(pairs: List[Tuple[Residue, Residue]], target_key: Tuple[str, int]) -> int:
    for i, (r1, r2) in enumerate(pairs):
        if (r1.chain_id, r1.res_seq) == target_key or (r2.chain_id, r2.res_seq) == target_key:
            return i
    raise ValueError(
        f"Selected residue {target_key[0]}{target_key[1]} is not part of the paired duplex axis model."
    )



def choose_new_chain_ids(existing_chain_ids: List[str], n_needed: int) -> List[str]:
    used = set(existing_chain_ids)
    available = [cid for cid in CHAIN_ID_CANDIDATES if cid not in used]
    if len(available) < n_needed:
        raise ValueError(
            "Unable to assign new chain IDs for separated piece 2: "
            f"need {n_needed}, but only {len(available)} unused one-character IDs are available."
        )
    return available[:n_needed]



def separate_piece2_chains(records, piece2_keys: set, duplex_chains: List[str]) -> Dict[str, str]:
    existing_chain_ids = []
    seen_existing = set()
    for rec in records:
        if isinstance(rec, AtomRecord) and rec.chain_id not in seen_existing:
            existing_chain_ids.append(rec.chain_id)
            seen_existing.add(rec.chain_id)

    new_chain_ids = choose_new_chain_ids(existing_chain_ids, len(duplex_chains))
    chain_map = {old: new for old, new in zip(duplex_chains, new_chain_ids)}

    for rec in records:
        if isinstance(rec, AtomRecord):
            key = (rec.chain_id, rec.res_seq)
            if key in piece2_keys:
                rec.set_chain_id(chain_map[rec.chain_id])

    return chain_map


def prepare_bend_geometry(
    residues: "OrderedDict[Tuple[str, int], Residue]",
    selected_key: Tuple[str, int],
    axis_range_defs: Optional[List[AxisRangeDefinition]] = None,
) -> BendGeometryPreparation:
    """Prepare pivot, duplex, and local-axis geometry without changing atoms."""
    duplex = build_duplex_data(
        residues, selected_key=selected_key, axis_range_defs=axis_range_defs
    )
    pairs = duplex["pairs"]

    if selected_key not in residues:
        raise ValueError(
            f"Residue {selected_key[0]}{selected_key[1]} was not found in the input PDB."
        )
    if residues[selected_key].p_atom is None:
        raise ValueError(
            f"Residue {selected_key[0]}{selected_key[1]} does not contain a P atom."
        )

    pair_idx = find_pair_index(pairs, selected_key)
    pair_keys = tuple(
        (
            (first.chain_id, first.res_seq),
            (second.chain_id, second.res_seq),
        )
        for first, second in pairs
    )
    piece2_keys = frozenset(
        key for pair in pair_keys[pair_idx:] for key in pair
    )

    axis_point = duplex["axis_point"]
    axis_dir = duplex["axis_dir"]
    selected_p = residues[selected_key].p_coord()
    axis_foot = project_point_to_line(selected_p, axis_point, axis_dir)
    radial = v_sub(selected_p, axis_foot)
    radius = v_len(radial)
    if radius < EPS:
        raise ValueError(
            f"Selected P atom at {selected_key[0]}{selected_key[1]} lies too close "
            "to the estimated helix axis."
        )

    axis_range_used = duplex.get("axis_range_used")
    return BendGeometryPreparation(
        selected_key=selected_key,
        selected_p=selected_p,
        axis_point=axis_point,
        axis_dir=axis_dir,
        axis_foot=axis_foot,
        radial=radial,
        radius=radius,
        pair_idx=pair_idx,
        pair_keys=pair_keys,
        piece2_keys=piece2_keys,
        duplex_chains=tuple(duplex["chains"]),
        axis_source=str(duplex.get("axis_source") or "automatic whole-duplex axis"),
        axis_range_used=(
            format_axis_range_spec(axis_range_used)
            if axis_range_used is not None
            else None
        ),
    )


def build_bend_transform(
    preparation: BendGeometryPreparation,
    phi_deg: float,
    beta_deg: float,
    tau_deg: float = 0.0,
    align_mode: str = "y",
) -> BendTransform:
    """Build one candidate's full-helix transform without mutating records."""
    phi_deg = float(phi_deg)
    beta_deg = float(beta_deg)
    tau_deg = float(tau_deg)
    align_mode = normalize_align(align_mode)

    phi_rad = math.radians(phi_deg)
    beta_rad = math.radians(beta_deg)
    tau_rad = math.radians(tau_deg)
    radial_phi = rotate_vector(preparation.radial, preparation.axis_dir, phi_rad)
    hinge_point = v_add(preparation.axis_foot, radial_phi)
    hinge_dir = v_cross(preparation.axis_dir, radial_phi)
    if v_len(hinge_dir) < EPS:
        raise ValueError("Failed to construct a non-zero tangent direction for the hinge.")
    hinge_dir = v_norm(hinge_dir)

    twist_axis_point_pre_align = rotate_point_about_line(
        preparation.axis_foot, hinge_point, hinge_dir, beta_rad
    )
    twist_axis_dir = v_norm(
        rotate_vector(preparation.axis_dir, hinge_dir, beta_rad)
    )

    def transform_without_alignment(coord: Point3D) -> Point3D:
        transformed = rotate_point_about_line(
            coord, hinge_point, hinge_dir, beta_rad
        )
        if abs(tau_rad) > 0.0:
            transformed = rotate_point_about_line(
                transformed,
                twist_axis_point_pre_align,
                twist_axis_dir,
                tau_rad,
            )
        return transformed

    align_translation: Point3D = (0.0, 0.0, 0.0)
    if align_mode == "y":
        align_translation = v_sub(
            preparation.selected_p,
            transform_without_alignment(preparation.selected_p),
        )

    return BendTransform(
        phi_deg=phi_deg,
        beta_deg=beta_deg,
        tau_deg=tau_deg,
        hinge_point=hinge_point,
        hinge_dir=hinge_dir,
        twist_axis_point_pre_align=twist_axis_point_pre_align,
        twist_axis_point=v_add(twist_axis_point_pre_align, align_translation),
        twist_axis_dir=twist_axis_dir,
        align_mode=align_mode,
        align_translation=align_translation,
    )



def bend_structure(
    records,
    residues: "OrderedDict[Tuple[str, int], Residue]",
    input_pdb: str,
    selected_key: Tuple[str, int],
    phi_deg: float,
    beta_deg: float,
    tau_deg: float = 0.0,
    sep_mode: str = "n",
    align_mode: str = "y",
    output_pdb: Optional[str] = None,
    axis_range_defs: Optional[List[AxisRangeDefinition]] = None,
) -> Tuple[str, Dict[str, object]]:
    preparation = prepare_bend_geometry(
        residues, selected_key, axis_range_defs=axis_range_defs
    )
    transform = build_bend_transform(
        preparation,
        phi_deg=phi_deg,
        beta_deg=beta_deg,
        tau_deg=tau_deg,
        align_mode=align_mode,
    )

    for rec in records:
        if isinstance(rec, AtomRecord):
            key = (rec.chain_id, rec.res_seq)
            if key in preparation.piece2_keys:
                rec.set_coord(transform.transform_coord(rec.coord()))

    piece2_chain_map = None
    if sep_mode == "y":
        piece2_chain_map = separate_piece2_chains(
            records,
            set(preparation.piece2_keys),
            list(preparation.duplex_chains),
        )

    out_path = normalize_output_path(output_pdb) or make_output_name(input_pdb, phi_deg, beta_deg, tau_deg, sep_mode=sep_mode)
    info = {
        "pair_idx": preparation.pair_idx,
        "n_pairs": len(preparation.pair_keys),
        "piece2_pair_start": preparation.pair_idx + 1,
        "piece2_pair_end": len(preparation.pair_keys),
        "radius": preparation.radius,
        "axis_point": preparation.axis_point,
        "axis_dir": preparation.axis_dir,
        "hinge_point": transform.hinge_point,
        "hinge_dir": transform.hinge_dir,
        "phi_deg": transform.phi_deg,
        "beta_deg": transform.beta_deg,
        "tau_deg": transform.tau_deg,
        "twist_axis_point_pre_align": transform.twist_axis_point_pre_align,
        "twist_axis_point": transform.twist_axis_point,
        "twist_axis_dir": transform.twist_axis_dir,
        "align_mode": transform.align_mode,
        "align_translation": transform.align_translation,
        "pair": preparation.pair_keys[preparation.pair_idx],
        "sep_mode": sep_mode,
        "piece2_chain_map": piece2_chain_map,
        "duplex_chains": preparation.duplex_chains,
        "axis_source": preparation.axis_source,
        "axis_range_used": preparation.axis_range_used,
        "output_pdb": out_path,
    }
    return out_path, info



def rewrite_ter_line(template_line: str, last_atom: Optional[AtomRecord]) -> str:
    if last_atom is None:
        return template_line

    raw = template_line.rstrip("\n")
    if len(raw) < 27:
        raw = raw.ljust(27)

    serial_field = raw[6:11] if len(raw) >= 11 else "     "
    updated = (
        f"{raw[:6]}{serial_field}{raw[11:17]}"
        f"{last_atom.line[17:20]}{raw[20:21]}{last_atom.chain_id}{last_atom.res_seq:4d}"
        f"{raw[26:]}"
    )
    return updated + "\n"



def write_pdb(records, out_path: str, update_ter: bool = False) -> None:
    last_atom: Optional[AtomRecord] = None
    with open(out_path, "w") as out:
        for rec in records:
            if isinstance(rec, AtomRecord):
                out.write(rec.to_line())
                last_atom = rec
            else:
                if update_ter and rec.line[:3] == "TER":
                    out.write(rewrite_ter_line(rec.line, last_atom))
                else:
                    out.write(rec.line)



def clone_atom_record(atom: AtomRecord) -> AtomRecord:
    return AtomRecord(
        line=atom.line,
        record_name=atom.record_name,
        atom_name=atom.atom_name,
        chain_id=atom.chain_id,
        res_seq=atom.res_seq,
        x=atom.x,
        y=atom.y,
        z=atom.z,
    )



def clone_records(records):
    cloned = []
    for rec in records:
        if isinstance(rec, AtomRecord):
            cloned.append(clone_atom_record(rec))
        else:
            cloned.append(RawRecord(line=rec.line))
    return cloned



def atom_to_line_with_serial(atom: AtomRecord, serial: int) -> str:
    return (
        f"{atom.line[:6]}{serial:5d}{atom.line[11:21]}{atom.chain_id}{atom.line[22:30]}"
        f"{atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f}{atom.line[54:]}"
    )



def ter_line_with_serial(serial: int, last_atom: AtomRecord) -> str:
    return f"TER   {serial:5d}      {last_atom.line[17:20]} {last_atom.chain_id}{last_atom.res_seq:4d}\n"



def build_full_helix_transform(info: Dict[str, object]):
    hinge_point = info["hinge_point"]
    hinge_dir = info["hinge_dir"]
    twist_axis_point_pre_align = info["twist_axis_point_pre_align"]
    twist_axis_dir = info["twist_axis_dir"]
    beta_rad = math.radians(float(info["beta_deg"]))
    tau_rad = math.radians(float(info["tau_deg"]))
    align_mode = str(info["align_mode"])
    align_translation = info["align_translation"]

    def transform(coord: Tuple[float, float, float]) -> Tuple[float, float, float]:
        new_coord = rotate_point_about_line(coord, hinge_point, hinge_dir, beta_rad)
        if abs(tau_rad) > 0.0:
            new_coord = rotate_point_about_line(new_coord, twist_axis_point_pre_align, twist_axis_dir, tau_rad)
        if align_mode == "y":
            new_coord = v_add(new_coord, align_translation)
        return new_coord

    return transform


def build_origin_overlay_chain_maps(
    chain_order: Sequence[str],
) -> Tuple[Dict[str, str], Dict[str, str]]:
    """Return the source-to-overlay chain maps used by ``-ori`` output."""
    chains = list(chain_order)
    if not chains:
        raise ValueError("No duplex chains are available for origin-overlay mapping.")
    overlay_ids = choose_new_chain_ids([], len(chains) * 2)
    model1 = {old: overlay_ids[i] for i, old in enumerate(chains)}
    model2 = {
        old: overlay_ids[len(chains) + i] for i, old in enumerate(chains)
    }
    return model1, model2


def _normalized_screening_token(value: str) -> str:
    return value.strip().lower().replace("-", "_").replace(" ", "_")


def _normalized_atom_name(value: str) -> str:
    return value.strip().upper().replace("*", "'")


def _finite_point3d(value: object, label: str) -> Point3D:
    if isinstance(value, str):
        raise ValueError(f"{label} must contain exactly three numeric values.")
    try:
        values = tuple(float(component) for component in value)  # type: ignore[union-attr]
    except (TypeError, ValueError):
        raise ValueError(f"{label} must contain exactly three numeric values.")
    if len(values) != 3:
        raise ValueError(f"{label} must contain exactly three numeric values.")
    if not all(math.isfinite(component) for component in values):
        raise ValueError(f"{label} values must be finite.")
    return values  # type: ignore[return-value]


def _finite_unit_vector(value: object, label: str) -> Point3D:
    vector = _finite_point3d(value, label)
    if v_len(vector) < EPS:
        raise ValueError(f"{label} cannot be zero.")
    return v_norm(vector)


def parse_overlay_atom_selector(value: object) -> OverlayAtomSelector:
    """Parse an atom selector in the chain namespace of the origin overlay."""
    if isinstance(value, OverlayAtomSelector):
        chain_id = value.chain_id
        res_seq = value.res_seq
        atom_name = value.atom_name
    elif isinstance(value, (tuple, list)) and len(value) == 3:
        chain_id, res_seq, atom_name = value
    elif isinstance(value, str):
        text = value.strip()
        separated = text.replace("/", ":").split(":")
        if len(separated) == 3:
            chain_id, residue_text, atom_name = (
                part.strip() for part in separated
            )
            try:
                res_seq = int(residue_text)
            except ValueError:
                raise ValueError(
                    f"Invalid overlay atom selector '{value}'. Use A:30:P, A30:P, or 30A:P."
                )
        elif len(separated) == 2:
            residue_selector, atom_name = (
                part.strip() for part in separated
            )
            match = re.fullmatch(
                r"(?:([A-Za-z0-9])(-?\d+)|(-?\d+)([A-Za-z0-9]))",
                residue_selector,
            )
            if match is None:
                raise ValueError(
                    f"Invalid overlay atom selector '{value}'. Use A:30:P, A30:P, or 30A:P."
                )
            if match.group(1) is not None:
                chain_id = match.group(1)
                res_seq = int(match.group(2))
            else:
                chain_id = match.group(4)
                res_seq = int(match.group(3))
        else:
            raise ValueError(
                f"Invalid overlay atom selector '{value}'. Use A:30:P, A30:P, or 30A:P."
            )
    else:
        raise ValueError(
            "Overlay atom selector must be A:30:P, a compact equivalent, or "
            "(chain, residue, atom)."
        )

    chain_text = str(chain_id).strip()
    atom_text = str(atom_name).strip()
    try:
        residue_number = int(res_seq)
    except (TypeError, ValueError):
        raise ValueError("Overlay atom residue number must be an integer.")
    if len(chain_text) != 1 or not atom_text:
        raise ValueError(
            "Overlay atom selector requires one chain-ID character and a nonblank atom name."
        )
    return OverlayAtomSelector(chain_text, residue_number, atom_text)


def _overlay_source_atom(
    context: ScreeningContext,
    selector: object,
) -> Tuple[ScreeningAtom, bool]:
    parsed = parse_overlay_atom_selector(selector)
    chain_mapping = context.overlay_chain_map
    if parsed.chain_id not in chain_mapping:
        available = ", ".join(sorted(chain_mapping))
        raise ValueError(
            f"Overlay chain '{parsed.chain_id}' is not available; choose one of {available}."
        )
    source_chain, is_transformed = chain_mapping[parsed.chain_id]
    wanted_name = _normalized_atom_name(parsed.atom_name)
    matches = [
        atom
        for atom in context.atoms
        if atom.source_chain_id == source_chain
        and atom.res_seq == parsed.res_seq
        and _normalized_atom_name(atom.atom_name) == wanted_name
    ]
    label = f"{parsed.chain_id}:{parsed.res_seq}:{parsed.atom_name}"
    if not matches:
        raise ValueError(f"Overlay atom selector '{label}' did not match any atom.")
    if len(matches) > 1:
        raise ValueError(
            f"Overlay atom selector '{label}' matched multiple atoms; choose an unambiguous atom."
        )
    return matches[0], is_transformed


def resolve_overlay_atom_coordinate(
    context: ScreeningContext,
    selector: object,
    transform: BendTransform,
) -> Point3D:
    """Resolve an A/B original or C/D transformed origin-overlay atom."""
    atom, is_transformed = _overlay_source_atom(context, selector)
    return transform.transform_coord(atom.coord) if is_transformed else atom.coord


def _compile_screening_point(
    context: ScreeningContext,
    reference: ScreeningPoint,
) -> Callable[[BendTransform], Point3D]:
    if not isinstance(reference, ScreeningPoint):
        raise ValueError("Screening point must be a ScreeningPoint value.")
    kind = _normalized_screening_token(reference.kind)
    if kind in ("overlay_atom", "atom"):
        atom, is_transformed = _overlay_source_atom(context, reference.value)
        if is_transformed:
            return lambda transform: transform.transform_coord(atom.coord)
        return lambda _transform: atom.coord
    if kind in ("xyz", "xyz_point", "point"):
        point = _finite_point3d(reference.value, "XYZ screening point")
        return lambda _transform: point
    if kind in (
        "phi_corrected_pivot",
        "phi_corrected_pivot_p",
        "pivot",
    ):
        return lambda transform: transform.hinge_point
    raise ValueError(
        f"Unknown screening point kind '{reference.kind}'; expected overlay_atom, "
        "xyz, or phi_corrected_pivot."
    )


def resolve_screening_point(
    context: ScreeningContext,
    reference: ScreeningPoint,
    transform: BendTransform,
) -> Point3D:
    """Resolve a point reference for one bend candidate."""
    return _compile_screening_point(context, reference)(transform)


def _axis_pair_references(
    axis: ScreeningAxis,
    expected_kind: str,
) -> Tuple[ScreeningPoint, ScreeningPoint]:
    if axis.point1 is not None and axis.point2 is not None:
        return axis.point1, axis.point2
    raw = axis.vector
    if isinstance(raw, (tuple, list)) and len(raw) == 2:
        first, second = raw
        if expected_kind == "xyz":
            return ScreeningPoint("xyz", first), ScreeningPoint("xyz", second)
        return ScreeningPoint("overlay_atom", first), ScreeningPoint(
            "overlay_atom", second
        )
    if isinstance(raw, (tuple, list)) and len(raw) == 6 and expected_kind == "xyz":
        return ScreeningPoint("xyz", raw[:3]), ScreeningPoint("xyz", raw[3:])
    raise ValueError(
        "Two-point axis vector source requires point1 and point2 (or a matching pair in vector)."
    )


def _compile_screening_axis(
    context: ScreeningContext,
    axis: ScreeningAxis,
) -> Callable[[BendTransform], Tuple[Point3D, Point3D]]:
    if not isinstance(axis, ScreeningAxis):
        raise ValueError("Rotation screening requires a ScreeningAxis value.")
    source = _normalized_screening_token(axis.source)
    if source in ("local", "local_axis", "axis_range", "axis_ranges"):
        point = context.preparation.axis_point
        direction = context.preparation.axis_dir
        return lambda _transform: (point, direction)
    if source not in ("geometric", "geometry"):
        raise ValueError(
            f"Unknown screening axis source '{axis.source}'; expected local_axis or geometric."
        )
    if axis.point is None:
        raise ValueError("A geometric screening axis requires an axis point.")
    point_resolver = _compile_screening_point(context, axis.point)
    vector_source = _normalized_screening_token(axis.vector_source)

    if vector_source in ("direct", "direct_vector", "vector"):
        direction = _finite_unit_vector(axis.vector, "Axis direction vector")
        return lambda transform: (point_resolver(transform), direction)

    if vector_source in ("two_xyz", "two_xyz_points", "xyz_points"):
        first_ref, second_ref = _axis_pair_references(axis, "xyz")
        first_resolver = _compile_screening_point(context, first_ref)
        second_resolver = _compile_screening_point(context, second_ref)
    elif vector_source in ("two_atoms", "two_overlay_atoms", "overlay_atoms"):
        first_ref, second_ref = _axis_pair_references(axis, "overlay_atom")
        first_resolver = _compile_screening_point(context, first_ref)
        second_resolver = _compile_screening_point(context, second_ref)
    elif vector_source in ("normal", "normal_vectors", "normal_to_two_vectors"):
        raw_first = axis.normal_vector1
        raw_second = axis.normal_vector2
        if raw_first is None or raw_second is None:
            raw = axis.vector
            if isinstance(raw, (tuple, list)) and len(raw) == 2:
                raw_first, raw_second = raw
            elif isinstance(raw, (tuple, list)) and len(raw) == 6:
                raw_first, raw_second = raw[:3], raw[3:]
            else:
                raise ValueError("Normal-vector axis source requires two XYZ vectors.")
        first = _finite_point3d(raw_first, "Normal-defining vector 1")
        second = _finite_point3d(raw_second, "Normal-defining vector 2")
        normal = v_cross(first, second)
        if v_len(normal) < EPS:
            raise ValueError("Normal-defining vectors cannot be zero or parallel.")
        direction = v_norm(normal)
        return lambda transform: (point_resolver(transform), direction)
    else:
        raise ValueError(
            f"Unknown geometric axis vector source '{axis.vector_source}'."
        )

    def resolve_from_points(transform: BendTransform) -> Tuple[Point3D, Point3D]:
        first = first_resolver(transform)
        second = second_resolver(transform)
        direction = v_sub(second, first)
        if v_len(direction) < EPS:
            raise ValueError("Axis vector-defining points cannot coincide.")
        return point_resolver(transform), v_norm(direction)

    return resolve_from_points


def resolve_screening_axis(
    context: ScreeningContext,
    axis: ScreeningAxis,
    transform: BendTransform,
) -> Tuple[Point3D, Point3D]:
    """Resolve ``(axis_point, unit_axis_direction)`` for one candidate."""
    return _compile_screening_axis(context, axis)(transform)


def distance_between_points(point1: Point3D, point2: Point3D) -> float:
    """Return the Euclidean distance between two finite XYZ points."""
    first = _finite_point3d(point1, "Distance point 1")
    second = _finite_point3d(point2, "Distance point 2")
    return v_len(v_sub(second, first))


def signed_projected_angle_deg(
    point1: Point3D,
    point2: Point3D,
    axis_point: Point3D,
    axis_dir: Point3D,
) -> float:
    """Return the signed point-1-to-point-2 rotation around an axis."""
    first = _finite_point3d(point1, "Rotation point 1")
    second = _finite_point3d(point2, "Rotation point 2")
    origin = _finite_point3d(axis_point, "Rotation axis point")
    direction = _finite_unit_vector(axis_dir, "Rotation axis direction")
    first_rel = v_sub(first, origin)
    second_rel = v_sub(second, origin)
    first_projected = v_sub(
        first_rel, v_scale(direction, v_dot(first_rel, direction))
    )
    second_projected = v_sub(
        second_rel, v_scale(direction, v_dot(second_rel, direction))
    )
    if v_len(first_projected) < EPS or v_len(second_projected) < EPS:
        raise ValueError(
            "Rotation screening points must not lie on the rotation axis."
        )
    first_unit = v_norm(first_projected)
    second_unit = v_norm(second_projected)
    sine = v_dot(direction, v_cross(first_unit, second_unit))
    cosine = max(-1.0, min(1.0, v_dot(first_unit, second_unit)))
    return math.degrees(math.atan2(sine, cosine))


def wrapped_angle_error_deg(actual: float, target: float) -> float:
    """Return the absolute shortest circular difference in degrees."""
    actual_value = float(actual)
    target_value = float(target)
    if not math.isfinite(actual_value) or not math.isfinite(target_value):
        raise ValueError("Rotation angles must be finite.")
    difference = (actual_value - target_value + 180.0) % 360.0 - 180.0
    return abs(difference)


def _normalize_angle_name(name: str) -> str:
    normalized = _normalized_screening_token(name)
    if normalized.endswith("_deg"):
        normalized = normalized[:-4]
    if normalized not in ("phi", "beta", "tau"):
        raise ValueError(
            f"Unknown screened angle '{name}'; expected phi, beta, or tau."
        )
    return normalized


def _angle_range_value_count(angle_range: ScreenAngleRange) -> int:
    start = float(angle_range.start)
    stop = float(angle_range.stop)
    step = float(angle_range.step)
    if not all(math.isfinite(value) for value in (start, stop, step)):
        raise ValueError("Screen angle range values must be finite.")
    if step <= 0.0:
        raise ValueError("Screen angle range step must be positive.")
    span = abs(stop - start)
    if not math.isfinite(span):
        raise ValueError("Screen angle range span is too large.")
    if span == 0.0:
        return 1
    quotient = span / step
    if not math.isfinite(quotient):
        raise ValueError("Screen angle range contains too many values.")
    if quotient < 1.0:
        return 2
    nearest_integer = round(quotient)
    if math.isclose(quotient, nearest_integer, rel_tol=1.0e-12, abs_tol=1.0e-12):
        return int(nearest_integer) + 1
    return int(math.floor(quotient)) + 2


def inclusive_angle_values(angle_range: ScreenAngleRange) -> Tuple[float, ...]:
    """Expand an inclusive ascending or descending positive-step angle grid."""
    count = _angle_range_value_count(angle_range)
    start = float(angle_range.start)
    stop = float(angle_range.stop)
    step = float(angle_range.step)
    if count == 1:
        return (start,)
    direction = 1.0 if stop > start else -1.0
    values = [start + direction * step * index for index in range(count - 1)]
    return tuple(values) + (stop,)


def validate_screen_angle_ranges(
    ranges: Sequence[ScreenAngleRange],
    candidate_cap: int,
) -> Tuple[
    Tuple[ScreenAngleRange, ...],
    Tuple[Tuple[float, ...], ...],
    int,
]:
    """Validate one/two unique grids and enforce the total candidate cap."""
    supplied = tuple(ranges)
    if len(supplied) not in (1, 2):
        raise ValueError("Angle screening requires exactly one or two screened angles.")
    if isinstance(candidate_cap, bool) or int(candidate_cap) != candidate_cap:
        raise ValueError("Screening candidate cap must be a positive integer.")
    cap = int(candidate_cap)
    if cap <= 0:
        raise ValueError("Screening candidate cap must be a positive integer.")

    normalized: List[ScreenAngleRange] = []
    names = set()
    counts: List[int] = []
    for angle_range in supplied:
        if not isinstance(angle_range, ScreenAngleRange):
            raise ValueError("Each screening grid must be a ScreenAngleRange.")
        name = _normalize_angle_name(angle_range.name)
        if name in names:
            raise ValueError(f"Screened angle '{name}' was supplied more than once.")
        names.add(name)
        normalized_range = ScreenAngleRange(
            name,
            float(angle_range.start),
            float(angle_range.stop),
            float(angle_range.step),
        )
        normalized.append(normalized_range)
        counts.append(_angle_range_value_count(normalized_range))

    candidate_count = math.prod(counts)
    if candidate_count > cap:
        raise ValueError(
            f"Screening grid contains {candidate_count} candidates, exceeding the "
            f"candidate cap of {cap}."
        )
    values = tuple(inclusive_angle_values(item) for item in normalized)
    return tuple(normalized), values, candidate_count


def format_screening_grid_preview(
    ranges: Sequence[ScreenAngleRange],
    max_values_per_angle: int = 24,
) -> str:
    """Format calculated coarse-grid values and their Cartesian-product count."""
    normalized, grid_values, candidate_count = validate_screen_angle_ranges(
        ranges, candidate_cap=250000
    )
    lines = []
    for angle_range, values in zip(normalized, grid_values):
        formatted = [format_float_for_cli(value) for value in values]
        if len(formatted) > max_values_per_angle:
            displayed = formatted[:12] + ["..."] + formatted[-6:]
        else:
            displayed = formatted
        lines.append(
            f"{angle_range.name.capitalize()} grid ({len(values)}): "
            + ", ".join(displayed)
            + " deg"
        )
    lines.append(f"Total coarse candidates: {candidate_count}")
    return "\n".join(lines)


def prepare_screening_context_from_records(
    input_pdb: str,
    records: Iterable[object],
    residues: "OrderedDict[Tuple[str, int], Residue]",
    selected_key: Tuple[str, int],
    axis_range_defs: Optional[List[AxisRangeDefinition]] = None,
) -> ScreeningContext:
    """Snapshot already-parsed records into an immutable screening context."""
    preparation = prepare_bend_geometry(
        residues, selected_key, axis_range_defs=axis_range_defs
    )
    model1, model2 = build_origin_overlay_chain_maps(
        preparation.duplex_chains
    )
    atoms = tuple(
        ScreeningAtom(
            source_chain_id=record.chain_id,
            res_seq=record.res_seq,
            atom_name=record.atom_name,
            coord=record.coord(),
        )
        for record in records
        if isinstance(record, AtomRecord)
        and record.chain_id in preparation.duplex_chains
    )
    return ScreeningContext(
        input_pdb=input_pdb,
        preparation=preparation,
        atoms=atoms,
        model1_chain_map_items=tuple(model1.items()),
        model2_chain_map_items=tuple(model2.items()),
    )


def prepare_screening_context(
    input_pdb: str,
    pivot_residue: str,
    axis_range_specs: Optional[Iterable[str]] = None,
) -> ScreeningContext:
    """Read a PDB and prepare a non-mutating context for repeated candidates."""
    records, residues = read_pdb(input_pdb)
    raw_axis_defs = parse_axis_range_specs(axis_range_specs)
    available_chains: List[str] = []
    seen_chains = set()
    for residue in residues.values():
        if residue.p_atom is not None and residue.chain_id not in seen_chains:
            available_chains.append(residue.chain_id)
            seen_chains.add(residue.chain_id)
    axis_range_defs = resolve_axis_range_definitions(
        raw_axis_defs, available_chains
    )
    selected_key = resolve_selected_key(
        parse_residue_token(pivot_residue), residues
    )
    return prepare_screening_context_from_records(
        input_pdb=input_pdb,
        records=records,
        residues=residues,
        selected_key=selected_key,
        axis_range_defs=axis_range_defs,
    )


def screen_bend_angles(
    context: ScreeningContext,
    fixed_angles: Mapping[str, float],
    ranges: Sequence[ScreenAngleRange],
    request: ScreeningRequest,
    align_mode: str = "y",
    candidate_cap: int = 250000,
) -> ScreeningResult:
    """Screen a coarse grid, then adaptively refine its best local candidates."""
    if not isinstance(context, ScreeningContext):
        raise ValueError("Angle screening requires a ScreeningContext.")
    if not isinstance(request, ScreeningRequest):
        raise ValueError("Angle screening requires a ScreeningRequest.")
    normalized_ranges, grid_values, candidate_count = (
        validate_screen_angle_ranges(ranges, candidate_cap)
    )
    align_mode = normalize_align(align_mode)
    target = float(request.target)
    if not math.isfinite(target):
        raise ValueError("Screening target must be finite.")

    screened_names = {angle_range.name for angle_range in normalized_ranges}
    normalized_fixed: Dict[str, float] = {}
    for raw_name, raw_value in fixed_angles.items():
        name = _normalize_angle_name(raw_name)
        value = float(raw_value)
        if not math.isfinite(value):
            raise ValueError(f"Fixed angle '{name}' must be finite.")
        normalized_fixed[name] = value
    for name in ("phi", "beta", "tau"):
        if name not in screened_names and name not in normalized_fixed:
            raise ValueError(f"Fixed value for unscreened angle '{name}' is required.")

    point1_resolver = _compile_screening_point(context, request.point1)
    point2_resolver = _compile_screening_point(context, request.point2)
    mode = _normalized_screening_token(request.mode)
    if mode in ("distance", "screening_for_distance"):
        if target < 0.0:
            raise ValueError("Distance screening target cannot be negative.")
        axis_resolver = None
    elif mode in ("rotation", "screening_for_rotation"):
        if request.axis is None:
            raise ValueError("Rotation screening requires an axis definition.")
        axis_resolver = _compile_screening_axis(context, request.axis)
    else:
        raise ValueError(
            f"Unknown screening mode '{request.mode}'; expected distance or rotation."
        )

    best_evaluation = None
    first_candidate_error: Optional[Exception] = None

    evaluation_cache = {}

    def evaluate_angles(angles: Mapping[str, float]):
        nonlocal best_evaluation, first_candidate_error
        angle_key = (angles["phi"], angles["beta"], angles["tau"])
        if angle_key in evaluation_cache:
            return evaluation_cache[angle_key]
        try:
            transform = build_bend_transform(
                context.preparation,
                phi_deg=angles["phi"],
                beta_deg=angles["beta"],
                tau_deg=angles["tau"],
                align_mode=align_mode,
            )
            point1 = point1_resolver(transform)
            point2 = point2_resolver(transform)
            if axis_resolver is None:
                achieved = distance_between_points(point1, point2)
                error = abs(achieved - target)
            else:
                axis_point, axis_dir = axis_resolver(transform)
                achieved = signed_projected_angle_deg(
                    point1, point2, axis_point, axis_dir
                )
                error = wrapped_angle_error_deg(achieved, target)
        except (ValueError, OverflowError) as exc:
            if first_candidate_error is None:
                first_candidate_error = exc
            evaluation_cache[angle_key] = None
            return None

        candidate_angles = dict(angles)
        key = (error, angles["phi"], angles["beta"], angles["tau"])
        evaluation = (key, candidate_angles, achieved)
        evaluation_cache[angle_key] = evaluation
        if best_evaluation is None or key < best_evaluation[0]:
            best_evaluation = evaluation
        return evaluation

    strongest_coarse_candidates = []
    for candidate_values in itertools.product(*grid_values):
        angles = dict(normalized_fixed)
        for angle_range, value in zip(normalized_ranges, candidate_values):
            angles[angle_range.name] = value
        evaluation = evaluate_angles(angles)
        if evaluation is None:
            continue
        if (
            len(strongest_coarse_candidates) < SCREEN_REFINEMENT_SEED_LIMIT
            or evaluation[0] < strongest_coarse_candidates[-1][0]
        ):
            strongest_coarse_candidates.append(evaluation)
            strongest_coarse_candidates.sort(key=lambda item: item[0])
            del strongest_coarse_candidates[SCREEN_REFINEMENT_SEED_LIMIT:]

    if best_evaluation is None:
        detail = (
            f" First candidate error: {first_candidate_error}"
            if first_candidate_error is not None
            else ""
        )
        raise ValueError(
            "Screening grid did not contain any geometrically valid candidates."
            + detail
        )

    refinement_bounds = {
        angle_range.name: (
            min(angle_range.start, angle_range.stop),
            max(angle_range.start, angle_range.stop),
        )
        for angle_range in normalized_ranges
    }
    initial_refinement_steps = {
        angle_range.name: min(
            angle_range.step / 2.0,
            abs(angle_range.stop - angle_range.start) / 2.0,
        )
        for angle_range in normalized_ranges
    }
    screened_order = tuple(angle_range.name for angle_range in normalized_ranges)
    for seed_evaluation in strongest_coarse_candidates:
        current_evaluation = seed_evaluation
        refinement_steps = dict(initial_refinement_steps)
        for _round_index in range(SCREEN_REFINEMENT_MAX_ROUNDS):
            if not any(
                refinement_steps[name] > SCREEN_REFINEMENT_TOLERANCE_DEG
                for name in screened_order
            ):
                break

            local_best = current_evaluation
            for offsets in itertools.product((-1.0, 0.0, 1.0), repeat=len(screened_order)):
                if all(offset == 0.0 for offset in offsets):
                    continue
                trial_angles = dict(current_evaluation[1])
                changed = False
                for name, offset in zip(screened_order, offsets):
                    low, high = refinement_bounds[name]
                    trial_value = max(
                        low,
                        min(
                            high,
                            current_evaluation[1][name]
                            + offset * refinement_steps[name],
                        ),
                    )
                    if trial_value != current_evaluation[1][name]:
                        changed = True
                    trial_angles[name] = trial_value
                if not changed:
                    continue
                trial_evaluation = evaluate_angles(trial_angles)
                if (
                    trial_evaluation is not None
                    and trial_evaluation[0] < local_best[0]
                ):
                    local_best = trial_evaluation

            if local_best[0] < current_evaluation[0]:
                current_evaluation = local_best
            else:
                refinement_steps = {
                    name: step / 2.0 for name, step in refinement_steps.items()
                }

    best_key, best_angles, best_achieved = best_evaluation
    refinement_candidate_count = max(0, len(evaluation_cache) - candidate_count)
    return ScreeningResult(
        phi_deg=best_angles["phi"],
        beta_deg=best_angles["beta"],
        tau_deg=best_angles["tau"],
        achieved_value=best_achieved,
        error=best_key[0],
        candidate_count=candidate_count,
        refinement_candidate_count=refinement_candidate_count,
    )



def write_origin_overlay_pdb(
    records,
    out_path: str,
    transform_coord,
    chain_order: List[str],
) -> Dict[str, Dict[str, str]]:
    if not chain_order:
        raise ValueError("No duplex chains are available for --origin output.")

    chain_to_atoms: "OrderedDict[str, List[AtomRecord]]" = OrderedDict((cid, []) for cid in chain_order)
    for rec in records:
        if isinstance(rec, AtomRecord) and rec.chain_id in chain_to_atoms:
            chain_to_atoms[rec.chain_id].append(rec)

    missing = [cid for cid, atoms in chain_to_atoms.items() if not atoms]
    if missing:
        raise ValueError(
            "Unable to build --origin output because no atoms were found for chain(s): "
            + ", ".join(missing)
        )

    chain_map_model1, chain_map_model2 = build_origin_overlay_chain_maps(
        chain_order
    )

    serial = 1
    with open(out_path, "w") as out:
        for chain_map, apply_transform in ((chain_map_model1, False), (chain_map_model2, True)):
            for orig_chain in chain_order:
                last_atom: Optional[AtomRecord] = None
                for atom in chain_to_atoms[orig_chain]:
                    new_atom = clone_atom_record(atom)
                    new_atom.set_chain_id(chain_map[orig_chain])
                    if apply_transform:
                        new_atom.set_coord(transform_coord(atom.coord()))
                    out.write(atom_to_line_with_serial(new_atom, serial))
                    serial += 1
                    last_atom = new_atom
                if last_atom is not None:
                    out.write(ter_line_with_serial(serial, last_atom))
                    serial += 1
        out.write("END\n")

    return {
        "origin_chain_map_model1": chain_map_model1,
        "origin_chain_map_model2": chain_map_model2,
    }


# ---------------------------------------------------------------------------
# CLI / GUI helpers
# ---------------------------------------------------------------------------



def normalize_pivot_for_compare(token: str) -> Tuple[str, int]:
    chain_id, res_seq = parse_residue_token(token)
    return (chain_id.upper(), res_seq)



def merge_cli_value(
    label: str,
    positional,
    optional,
    normalizer=None,
):
    if positional is not None and optional is not None:
        lhs = normalizer(positional) if normalizer is not None else positional
        rhs = normalizer(optional) if normalizer is not None else optional
        if lhs != rhs:
            raise ValueError(
                f"Conflicting values were provided for {label}: {positional!r} and {optional!r}."
            )
        return positional
    if positional is not None:
        return positional
    return optional



def resolve_run_parameters(args) -> Tuple[str, str, float, float, float, str, str, str, Optional[str], List[str]]:
    input_pdb = merge_cli_value(
        label="input PDB",
        positional=args.input_pdb,
        optional=args.input_pdb_opt,
        normalizer=lambda path: os.path.abspath(path),
    )
    pivot_residue = merge_cli_value(
        label="pivot residue",
        positional=args.pivot_residue,
        optional=args.pivot_residue_opt,
        normalizer=normalize_pivot_for_compare,
    )
    phi_deg = merge_cli_value(
        label="phi",
        positional=args.phi_deg,
        optional=args.phi_deg_opt,
    )
    beta_deg = merge_cli_value(
        label="beta",
        positional=args.beta_deg,
        optional=args.beta_deg_opt,
    )
    tau_deg = 0.0 if args.tau_deg_opt is None else float(args.tau_deg_opt)
    sep_mode = normalize_sep(args.sep)
    align_mode = normalize_align(args.align)
    origin_mode = normalize_origin(args.origin)
    output_pdb = normalize_output_path(args.output_pdb)
    axis_range_specs = list(args.axis_ranges or [])

    missing = []
    if input_pdb is None:
        missing.append("input PDB")
    if pivot_residue is None:
        missing.append("pivot residue")
    if phi_deg is None:
        missing.append("phi")
    if beta_deg is None:
        missing.append("beta")

    if missing:
        raise ValueError(
            "Missing required parameter(s): "
            + ", ".join(missing)
            + ". Provide all four values positionally, with --input/--pivot/--phi/--beta, or use --gui."
        )

    return (
        input_pdb,
        pivot_residue,
        float(phi_deg),
        float(beta_deg),
        tau_deg,
        sep_mode,
        align_mode,
        origin_mode,
        output_pdb,
        axis_range_specs,
    )


def run_bending(
    input_pdb: str,
    pivot_residue: str,
    phi_deg: float,
    beta_deg: float,
    tau_deg: float = 0.0,
    sep_mode: str = "n",
    align_mode: str = "y",
    origin_mode: str = "n",
    output_pdb: Optional[str] = None,
    axis_range_specs: Optional[Iterable[str]] = None,
) -> Tuple[str, Dict[str, object]]:
    records, residues = read_pdb(input_pdb)
    origin_source_records = clone_records(records) if origin_mode == "y" else None

    raw_axis_defs = parse_axis_range_specs(axis_range_specs)
    available_chains = []
    seen_chains = set()
    for residue in residues.values():
        if residue.p_atom is not None and residue.chain_id not in seen_chains:
            available_chains.append(residue.chain_id)
            seen_chains.add(residue.chain_id)
    axis_range_defs = resolve_axis_range_definitions(raw_axis_defs, available_chains)

    selected_key = resolve_selected_key(parse_residue_token(pivot_residue), residues)
    out_path, info = bend_structure(
        records=records,
        residues=residues,
        input_pdb=input_pdb,
        selected_key=selected_key,
        phi_deg=phi_deg,
        beta_deg=beta_deg,
        tau_deg=tau_deg,
        sep_mode=sep_mode,
        align_mode=align_mode,
        output_pdb=output_pdb,
        axis_range_defs=axis_range_defs,
    )
    write_pdb(records, out_path, update_ter=(sep_mode == "y"))

    info["origin_mode"] = origin_mode
    info["origin_out_path"] = None
    info["origin_chain_map_model1"] = None
    info["origin_chain_map_model2"] = None
    info["axis_range_specs"] = [format_axis_range_spec(axis_def) for axis_def in axis_range_defs]

    if origin_mode == "y":
        if origin_source_records is None:
            raise ValueError("Internal error: missing source records for --origin output.")
        origin_out_path = make_origin_output_name(out_path)
        origin_info = write_origin_overlay_pdb(
            records=origin_source_records,
            out_path=origin_out_path,
            transform_coord=build_full_helix_transform(info),
            chain_order=list(info["duplex_chains"]),
        )
        info["origin_out_path"] = origin_out_path
        info.update(origin_info)

    return out_path, info


def format_run_summary(out_path: str, info: Dict[str, object]) -> str:
    pair_a, pair_b = info["pair"]
    axis_dir = info["axis_dir"]
    hinge_point = info["hinge_point"]
    hinge_dir = info["hinge_dir"]
    twist_axis_point = info["twist_axis_point"]
    twist_axis_dir = info["twist_axis_dir"]
    align_translation = info["align_translation"]
    piece1_end = int(info["piece2_pair_start"]) - 1

    lines = [
        f"Wrote {out_path}",
        f"Split base pair / pivot border: {pair_a[0]}{pair_a[1]} / {pair_b[0]}{pair_b[1]}",
        f"Piece #1 is fixed: base-pair range 1..{piece1_end} of {info['n_pairs']} in axial order",
        (
            "Piece #2 is movable: base-pair range "
            f"{info['piece2_pair_start']}..{info['piece2_pair_end']} of {info['n_pairs']} in axial order"
        ),
        f"Axis source: {info.get('axis_source') or 'automatic whole-duplex axis'}",
        f"Estimated helix radius at pivot P: {info['radius']:.3f} A",
        f"Axis direction: ({axis_dir[0]:.6f}, {axis_dir[1]:.6f}, {axis_dir[2]:.6f})",
        f"Hinge point: ({hinge_point[0]:.3f}, {hinge_point[1]:.3f}, {hinge_point[2]:.3f})",
        f"Hinge direction: ({hinge_dir[0]:.6f}, {hinge_dir[1]:.6f}, {hinge_dir[2]:.6f})",
        f"Tau twist of movable piece #2: {info['tau_deg']:.6f} deg",
        f"Twist axis point: ({twist_axis_point[0]:.3f}, {twist_axis_point[1]:.3f}, {twist_axis_point[2]:.3f})",
        f"Twist axis direction: ({twist_axis_dir[0]:.6f}, {twist_axis_dir[1]:.6f}, {twist_axis_dir[2]:.6f})",
        f"Align pivot P back to original position (--align): {info['align_mode']}",
    ]

    if info.get("align_mode") == "y":
        lines.append(
            "Alignment translation applied to movable piece #2: "
            f"({align_translation[0]:.3f}, {align_translation[1]:.3f}, {align_translation[2]:.3f})"
        )

    if info.get("sep_mode") == "y":
        chain_map = info.get("piece2_chain_map") or {}
        mapping_text = ", ".join(f"{old}->{new}" for old, new in chain_map.items())
        lines.append(f"Separated piece #2 chain IDs: {mapping_text}")

    if info.get("origin_mode") == "y":
        lines.append(f"Origin comparison PDB (--origin): {info['origin_out_path']}")
        chain_map_model1 = info.get("origin_chain_map_model1") or {}
        chain_map_model2 = info.get("origin_chain_map_model2") or {}
        duplex_chains = info.get("duplex_chains") or ()
        if chain_map_model1 and chain_map_model2 and duplex_chains:
            model1_ids = ", ".join(chain_map_model1[cid] for cid in duplex_chains if cid in chain_map_model1)
            model2_ids = ", ".join(chain_map_model2[cid] for cid in duplex_chains if cid in chain_map_model2)
            lines.append(f"Origin full-helix chain IDs: model 1 [{model1_ids}], model 2 [{model2_ids}]")

    return "\n".join(lines)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Bend a two-chain nucleic-acid helix at a chosen phosphorus site. "
            "Piece #1 stays fixed; piece #2 is moved by bend/twist operations."
        )
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=APP_TITLE,
    )
    parser.add_argument("input_pdb", nargs="?", help="input PDB file")
    parser.add_argument("pivot_residue", nargs="?", help="pivot phosphorus residue, e.g. A36 or 36A")
    parser.add_argument("phi_deg", nargs="?", type=float, help="phi angle in degrees")
    parser.add_argument("beta_deg", nargs="?", type=float, help="beta bend angle in degrees")

    parser.add_argument("--input_pdb", "--input", dest="input_pdb_opt", help="input PDB file")
    parser.add_argument("--pivot_residue", "--pivot", dest="pivot_residue_opt", help="pivot phosphorus residue")
    parser.add_argument("--phi_deg", "--phi", dest="phi_deg_opt", type=float, help="phi angle in degrees")
    parser.add_argument("--beta_deg", "--beta", dest="beta_deg_opt", type=float, help="beta bend angle in degrees")
    parser.add_argument("--tau_deg", "--tau", dest="tau_deg_opt", type=float, help="additional twist angle in degrees (default: 0)")
    parser.add_argument(
        "-o",
        "--output",
        "--output_pdb",
        dest="output_pdb",
        help="optional output PDB filename; if omitted, automatic *_PxByTz.pdb naming is used, with _sep added when --sep y",
    )
    parser.add_argument(
        "--axis_range",
        "--axis-range",
        dest="axis_ranges",
        action="append",
        default=[],
        help=(
            "optional local helical-axis residue range, repeatable; examples: "
            "A1-A35,B60-B26 or A36-A60,B25-B1. If multiple are supplied, the one "
            "containing the pivot residue is used."
        ),
    )
    parser.add_argument(
        "--sep",
        default="n",
        type=normalize_sep,
        choices=("y", "n"),
        help="y: give movable piece #2 new chain IDs; n: keep original chain IDs (default: n)",
    )
    parser.add_argument(
        "--align",
        default="y",
        type=normalize_align,
        choices=("y", "n"),
        help="y: translate movable piece #2 so the pivot P returns to its original position; n: no post-bend translation (default: y)",
    )
    parser.add_argument(
        "--origin",
        default="n",
        type=normalize_origin,
        choices=("y", "n"),
        help="y: also write a -ori PDB containing original and fully transformed helix overlays; n: write only the bent output (default: n)",
    )
    parser.add_argument("--gui", action="store_true", help="launch the graphical interface")
    return parser


def launch_gui(defaults: Optional[Dict[str, str]] = None) -> int:
    defaults = defaults or {}

    try:
        import tkinter as tk
        from tkinter import filedialog, messagebox, scrolledtext, ttk
    except Exception as exc:  # pragma: no cover - environment dependent
        print(f"Error: GUI mode requires tkinter ({exc}).", file=sys.stderr)
        return 1

    try:
        root = tk.Tk()
    except tk.TclError as exc:  # pragma: no cover - environment dependent
        print(f"Error: GUI mode requires a graphical display ({exc}).", file=sys.stderr)
        return 1

    root.title(APP_TITLE)
    apply_optional_icon(root, __file__)
    root.geometry("1020x840")
    root.resizable(True, True)
    root.columnconfigure(1, weight=1)
    root.rowconfigure(12, weight=1)

    input_var = tk.StringVar(value=defaults.get("input_pdb", ""))
    output_var = tk.StringVar(value=defaults.get("output_pdb", ""))
    pivot_var = tk.StringVar(value=defaults.get("pivot_residue", ""))
    phi_var = tk.StringVar(value=defaults.get("phi_deg", ""))
    beta_var = tk.StringVar(value=defaults.get("beta_deg", ""))
    tau_var = tk.StringVar(value=defaults.get("tau_deg", "0"))
    sep_var = tk.StringVar(value=defaults.get("sep", "n"))
    align_var = tk.StringVar(value=defaults.get("align", "y"))
    origin_var = tk.StringVar(value=defaults.get("origin", "n"))

    angle_vars = {"phi": phi_var, "beta": beta_var, "tau": tau_var}
    screen_angle_vars = {
        name: tk.BooleanVar(value=False) for name in ("phi", "beta", "tau")
    }
    screen_range_vars = {
        name: {
            "start": tk.StringVar(
                value=format_float_for_cli(DEFAULT_SCREEN_RANGES[name][0])
            ),
            "stop": tk.StringVar(
                value=format_float_for_cli(DEFAULT_SCREEN_RANGES[name][1])
            ),
            "step": tk.StringVar(value=format_float_for_cli(DEFAULT_SCREEN_STEP_DEG)),
        }
        for name in ("phi", "beta", "tau")
    }
    screen_mode_var = tk.StringVar(value="Screening for distance")
    screen_target_var = tk.StringVar()
    screen_point1_atom_var = tk.StringVar()
    screen_point2_source_var = tk.StringVar(value="Overlay atom")
    screen_point2_atom_var = tk.StringVar()
    screen_point2_xyz_vars = [tk.StringVar() for _ in range(3)]
    screen_axis_source_var = tk.StringVar(value="Geometric definition")
    screen_axis_point_source_var = tk.StringVar(value="XYZ point")
    screen_axis_point_xyz_vars = [tk.StringVar() for _ in range(3)]
    screen_axis_point_atom_var = tk.StringVar()
    screen_axis_vector_source_var = tk.StringVar(value="Direct vector")
    screen_axis_vector_vars = [tk.StringVar() for _ in range(3)]
    screen_axis_two_xyz_vars = [tk.StringVar() for _ in range(6)]
    screen_axis_two_atom_vars = [tk.StringVar() for _ in range(2)]
    screen_axis_normal_vars = [tk.StringVar() for _ in range(6)]
    screen_status_var = tk.StringVar(
        value="Check one or two angles, then configure the screening target."
    )
    screen_dialog_state: Dict[str, object] = {"window": None, "configured": False}

    help_text = {
        "input": (
            "Input PDB file\n\n"
            "Choose a PDB containing a two-chain nucleic-acid helix. The current bending model "
            "uses the paired P atoms from the two chains.\n\n"
            "Example: A60-heli.pdb"
        ),
        "output": (
            "Output PDB file\n\n"
            "Optional. Leave this blank to use automatic naming such as input_P0B30T0.pdb "
            "or input_P0B30T0_sep.pdb. If you provide a filename, that path is used for the "
            "main bent model. Automatic screening names add _scr after the angle values, "
            "for example input_P0B30T0_scr.pdb or input_P0B30T0_scr_sep.pdb."
        ),
        "pivot": (
            "Pivot P residue\n\n"
            "Residue containing the P atom that marks the border between the two pieces. "
            "Piece #1 is fixed. Piece #2 starts at the pivot base pair and is movable.\n\n"
            "Accepted forms: A36, 36A, A.36, 36.A"
        ),
        "phi": (
            "Phi angle (degrees)\n\n"
            "Phi chooses the hinge direction around the local helix axis. phi = 0 uses the "
            "tangent through the selected pivot P position; nonzero phi rotates that point "
            "around the helix axis before building the hinge.\n\n"
            "Positive phi follows the right-hand rule about the positive local helix axis: "
            "point your right thumb along +axis, and your curled fingers show the direction "
            "of increasing phi. Viewed from the +axis end looking toward the -axis end, "
            "positive phi is counterclockwise.\n\n"
            "For the automatic axis, +axis points toward increasing residue numbers on the "
            "first P-containing chain encountered in the PDB (file order, not alphabetical "
            "order). The pivot chain does not change this axis direction. For example, if "
            "chain A is encountered before chain B, increasing A residue numbers define "
            "+axis even when the pivot is on chain B; the B-chain pivot changes only the "
            "phi = 0 radial starting direction. With --axis_range, the start-to-end order "
            "of the first range sets +axis.\n\n"
            "Examples: 0, 90, -45"
        ),
        "beta": (
            "Beta bend angle (degrees)\n\n"
            "Rigid bend angle applied to movable piece #2 relative to fixed piece #1. "
            "Positive beta bends piece #2 away from the helical axis in the phi-defined direction.\n\n"
            "Example: 30"
        ),
        "tau": (
            "Tau twist angle (degrees)\n\n"
            "Additional twist of movable piece #2 around its bent helical axis. Positive tau "
            "is right-handed by the right-hand rule; negative tau is left-handed.\n\n"
            "Examples: 10, -15, 0"
        ),
        "screening": (
            "Angle screening\n\n"
            "Check Screen beside exactly one or two of phi, beta, and tau. The checked "
            "angle fields are replaced by inclusive From/To/Step grids configured in the "
            "Screening to achieve... window; unchecked angles stay fixed at their entered "
            "values.\n\n"
            "From, To, and Step are all in degrees. After testing the coarse grid, Bend "
            "Helix efficiently refines between steps around the strongest candidates down "
            "to 0.001 degree.\n\n"
            "Each candidate is evaluated against an origin-overlay distance or signed-rotation "
            "target. The closest candidate is used even when the exact target is unavailable, "
            "and the winning origin-overlay PDB is written automatically."
        ),
        "axis_range": (
            "Local helix-axis residue range(s)\n\n"
            "Optional. Use this when the full input is already bent and you want the bend "
            "operation to use a local straight segment for the helix axis. Enter one "
            "--axis_range per line. Each line is a comma-separated pair of chain ranges.\n\n"
            "Examples:\n"
            "A1-A35,B60-B26\n"
            "A36-A60,B25-B1\n\n"
            "If more than one line is provided, the line containing the pivot residue is used. "
            "If only one line is provided, it is used even when the pivot lies just outside it. "
            "The start-to-end order of the first range sets the positive axis direction."
        ),
        "sep": (
            "Separate output pieces (--sep)\n\n"
            "n: keep the final model as one helix under the original chain IDs.\n"
            "y: keep piece #1 and piece #2 separated by giving movable piece #2 new chain IDs."
        ),
        "align": (
            "Align pivot after bending (--align)\n\n"
            "y: after bend/twist, translate movable piece #2 so the pivot residue P atom "
            "returns to its original position before bending.\n"
            "n: keep the unaligned placement used by V2.1."
        ),
        "origin": (
            "Write origin overlay (--origin)\n\n"
            "n: write only the main bent PDB output.\n"
            "y: also write a -ori PDB containing the original full helix and the same rigid "
            "full-helix transform used for movable piece #2. The two copies use sequential chain IDs."
        ),
    }

    def show_help(key: str) -> None:
        messagebox.showinfo(f"{APP_TITLE} help", help_text[key], parent=root)

    def help_button(row: int, key: str):
        btn = tk.Button(
            root,
            text="?",
            command=lambda k=key: show_help(k),
            bg="#d9ecff",
            activebackground="#c4e0ff",
            highlightbackground="#d9ecff",
            width=2,
            relief="raised",
        )
        btn.grid(row=row, column=3, sticky="w", padx=(0, 8), pady=4)
        return btn

    ttk.Label(
        root,
        text="Bend helix: Piece #1 is fixed; Piece #2 is movable and receives beta/tau transformations.",
        font=("TkDefaultFont", 10, "bold"),
    ).grid(row=0, column=0, columnspan=4, sticky="w", padx=8, pady=(8, 6))

    ttk.Label(root, text="Input PDB file").grid(row=1, column=0, sticky="w", padx=8, pady=4)
    input_entry = ttk.Entry(root, textvariable=input_var)
    input_entry.grid(row=1, column=1, sticky="ew", padx=8, pady=4)

    def browse_input_file() -> None:
        path = filedialog.askopenfilename(
            title="Select input PDB",
            filetypes=[("PDB files", "*.pdb"), ("Text/PDB files", "*.txt"), ("All files", "*")],
        )
        if path:
            input_var.set(path)

    ttk.Button(root, text="Browse...", command=browse_input_file).grid(row=1, column=2, sticky="ew", padx=(0, 8), pady=4)
    help_button(1, "input")

    ttk.Label(root, text="Output PDB file (optional)").grid(row=2, column=0, sticky="w", padx=8, pady=4)
    ttk.Entry(root, textvariable=output_var).grid(row=2, column=1, sticky="ew", padx=8, pady=4)

    def browse_output_file() -> None:
        path = filedialog.asksaveasfilename(
            title="Select output PDB",
            defaultextension=".pdb",
            filetypes=[("PDB files", "*.pdb"), ("All files", "*")],
        )
        if path:
            output_var.set(path)

    ttk.Button(root, text="Save as...", command=browse_output_file).grid(row=2, column=2, sticky="ew", padx=(0, 8), pady=4)
    help_button(2, "output")

    ttk.Label(root, text="Pivot P residue (piece border)").grid(row=3, column=0, sticky="w", padx=8, pady=4)
    ttk.Entry(root, textvariable=pivot_var).grid(row=3, column=1, sticky="ew", padx=8, pady=4)
    help_button(3, "pivot")

    ttk.Label(root, text="Phi: hinge direction (degrees)").grid(row=4, column=0, sticky="w", padx=8, pady=4)
    phi_entry = ttk.Entry(root, textvariable=phi_var)
    phi_entry.grid(row=4, column=1, sticky="ew", padx=8, pady=4)
    phi_screen_check = ttk.Checkbutton(
        root, text="Screen", variable=screen_angle_vars["phi"]
    )
    phi_screen_check.grid(row=4, column=2, sticky="w", padx=(0, 8), pady=4)
    help_button(4, "phi")

    ttk.Label(root, text="Beta: bend movable piece #2 (degrees)").grid(row=5, column=0, sticky="w", padx=8, pady=4)
    beta_entry = ttk.Entry(root, textvariable=beta_var)
    beta_entry.grid(row=5, column=1, sticky="ew", padx=8, pady=4)
    beta_screen_check = ttk.Checkbutton(
        root, text="Screen", variable=screen_angle_vars["beta"]
    )
    beta_screen_check.grid(row=5, column=2, sticky="w", padx=(0, 8), pady=4)
    help_button(5, "beta")

    ttk.Label(root, text="Tau: twist movable piece #2 (degrees)").grid(row=6, column=0, sticky="w", padx=8, pady=4)
    tau_entry = ttk.Entry(root, textvariable=tau_var)
    tau_entry.grid(row=6, column=1, sticky="ew", padx=8, pady=4)
    tau_screen_check = ttk.Checkbutton(
        root, text="Screen", variable=screen_angle_vars["tau"]
    )
    tau_screen_check.grid(row=6, column=2, sticky="w", padx=(0, 8), pady=4)
    help_button(6, "tau")

    angle_entries = {"phi": phi_entry, "beta": beta_entry, "tau": tau_entry}
    screen_checkbuttons = {
        "phi": phi_screen_check,
        "beta": beta_screen_check,
        "tau": tau_screen_check,
    }

    ttk.Label(root, text="Angle screening target").grid(
        row=7, column=0, sticky="w", padx=8, pady=4
    )
    screening_button = ttk.Button(root, text="Screening to achieve...", state="disabled")
    screening_button.grid(row=7, column=1, sticky="w", padx=8, pady=4)
    ttk.Label(root, textvariable=screen_status_var, wraplength=300).grid(
        row=7, column=2, sticky="w", padx=(0, 8), pady=4
    )
    help_button(7, "screening")

    ttk.Label(root, text="Local axis range(s), one per line").grid(row=8, column=0, sticky="nw", padx=8, pady=4)
    axis_range_text = tk.Text(root, height=3, width=48, wrap="none")
    axis_range_text.grid(row=8, column=1, sticky="ew", padx=8, pady=4)
    axis_default = defaults.get("axis_ranges", "")
    if axis_default:
        axis_range_text.insert("1.0", axis_default)
    ttk.Label(root, text="Example: A1-A35,B60-B26").grid(row=8, column=2, sticky="w", padx=(0, 8), pady=4)
    help_button(8, "axis_range")

    ttk.Label(root, text="Separate piece #2 chain IDs (--sep)").grid(row=9, column=0, sticky="w", padx=8, pady=4)
    sep_box = ttk.Combobox(root, textvariable=sep_var, values=("n", "y"), state="readonly", width=6)
    sep_box.grid(row=9, column=1, sticky="w", padx=8, pady=4)
    sep_box.set(sep_var.get())
    help_button(9, "sep")

    ttk.Label(root, text="Realign pivot P after bending (--align)").grid(row=10, column=0, sticky="w", padx=8, pady=4)
    align_box = ttk.Combobox(root, textvariable=align_var, values=("y", "n"), state="readonly", width=6)
    align_box.grid(row=10, column=1, sticky="w", padx=8, pady=4)
    align_box.set(align_var.get())
    help_button(10, "align")

    ttk.Label(root, text="Write origin overlay PDB (--origin)").grid(row=11, column=0, sticky="w", padx=8, pady=4)
    origin_box = ttk.Combobox(root, textvariable=origin_var, values=("n", "y"), state="readonly", width=6)
    origin_box.grid(row=11, column=1, sticky="w", padx=8, pady=4)
    origin_box.set(origin_var.get())
    help_button(11, "origin")

    result_text = scrolledtext.ScrolledText(root, height=14, width=90, wrap="word")
    result_text.grid(row=12, column=0, columnspan=4, sticky="nsew", padx=8, pady=(8, 4))

    button_frame = ttk.Frame(root)
    button_frame.grid(row=13, column=0, columnspan=4, sticky="e", padx=8, pady=(4, 8))

    def get_axis_range_specs_from_gui() -> List[str]:
        text = axis_range_text.get("1.0", tk.END)
        specs: List[str] = []
        for line in text.replace(";", "\n").splitlines():
            line = line.strip()
            if line:
                specs.append(line)
        return specs

    def selected_screen_angles() -> List[str]:
        return [
            name
            for name in ("phi", "beta", "tau")
            if bool(screen_angle_vars[name].get())
        ]

    def gui_finite_float(raw_value: str, label: str) -> float:
        text = raw_value.strip()
        if not text:
            raise ValueError(f"Please provide {label}.")
        try:
            value = float(text)
        except ValueError:
            raise ValueError(f"{label} must be a numeric value.")
        if not math.isfinite(value):
            raise ValueError(f"{label} must be finite.")
        return value

    def gui_xyz(values, label: str) -> Point3D:
        return tuple(
            gui_finite_float(var.get(), f"{label} {axis_name}")
            for var, axis_name in zip(values, ("x", "y", "z"))
        )  # type: ignore[return-value]

    def gui_overlay_point(raw_value: str, label: str) -> ScreeningPoint:
        selector = raw_value.strip()
        if not selector:
            raise ValueError(f"Please provide {label} as CHAIN:RESIDUE:ATOM.")
        parse_overlay_atom_selector(selector)
        return ScreeningPoint("overlay_atom", selector)

    def screening_ranges_from_gui() -> List[ScreenAngleRange]:
        selected = selected_screen_angles()
        if len(selected) not in (1, 2):
            raise ValueError("Check Screen for exactly one or two angles.")
        ranges: List[ScreenAngleRange] = []
        for name in selected:
            values = screen_range_vars[name]
            ranges.append(
                ScreenAngleRange(
                    name=name,
                    start=gui_finite_float(
                        values["start"].get(), f"{name} From (degrees)"
                    ),
                    stop=gui_finite_float(
                        values["stop"].get(), f"{name} To (degrees)"
                    ),
                    step=gui_finite_float(
                        values["step"].get(), f"{name} Step (degrees)"
                    ),
                )
            )
        validate_screen_angle_ranges(ranges, candidate_cap=250000)
        return ranges

    def screening_request_from_gui() -> ScreeningRequest:
        mode_label = screen_mode_var.get().strip()
        is_rotation = mode_label == "Screening for rotation"
        mode = "rotation" if is_rotation else "distance"
        target_label = "target rotation (degrees)" if is_rotation else "target distance (angstroms)"
        target = gui_finite_float(screen_target_var.get(), target_label)
        if not is_rotation and target < 0.0:
            raise ValueError("Target distance cannot be negative.")

        point1 = gui_overlay_point(
            screen_point1_atom_var.get(), "screening endpoint 1 overlay atom"
        )
        point2_source = screen_point2_source_var.get().strip()
        if point2_source == "Overlay atom":
            point2 = gui_overlay_point(
                screen_point2_atom_var.get(), "screening endpoint 2 overlay atom"
            )
        elif point2_source == "XYZ point":
            point2 = ScreeningPoint(
                "xyz", gui_xyz(screen_point2_xyz_vars, "screening endpoint 2")
            )
        elif point2_source == "Phi-corrected pivot P" and is_rotation:
            point2 = ScreeningPoint("phi_corrected_pivot")
        else:
            raise ValueError("Choose a valid source for screening endpoint 2.")

        if not is_rotation:
            return ScreeningRequest(
                mode=mode,
                target=target,
                point1=point1,
                point2=point2,
            )

        axis_source = screen_axis_source_var.get().strip()
        if axis_source == "Local axis range(s)":
            if not get_axis_range_specs_from_gui():
                raise ValueError(
                    "Local-axis rotation screening requires at least one Local axis "
                    "range in the main Bend Helix window."
                )
            axis = ScreeningAxis(source="local_axis")
        elif axis_source == "Geometric definition":
            if screen_axis_point_source_var.get().strip() == "Overlay atom":
                axis_point = gui_overlay_point(
                    screen_axis_point_atom_var.get(), "rotation-axis point overlay atom"
                )
            else:
                axis_point = ScreeningPoint(
                    "xyz", gui_xyz(screen_axis_point_xyz_vars, "rotation-axis point")
                )

            vector_source = screen_axis_vector_source_var.get().strip()
            if vector_source == "Direct vector":
                axis = ScreeningAxis(
                    source="geometric",
                    point=axis_point,
                    vector_source="direct_vector",
                    vector=gui_xyz(screen_axis_vector_vars, "rotation-axis vector"),
                )
            elif vector_source == "Two XYZ points":
                axis = ScreeningAxis(
                    source="geometric",
                    point=axis_point,
                    vector_source="two_xyz_points",
                    point1=ScreeningPoint(
                        "xyz", gui_xyz(screen_axis_two_xyz_vars[:3], "axis-vector point 1")
                    ),
                    point2=ScreeningPoint(
                        "xyz", gui_xyz(screen_axis_two_xyz_vars[3:], "axis-vector point 2")
                    ),
                )
            elif vector_source == "Two overlay atoms":
                axis = ScreeningAxis(
                    source="geometric",
                    point=axis_point,
                    vector_source="two_overlay_atoms",
                    point1=gui_overlay_point(
                        screen_axis_two_atom_vars[0].get(), "axis-vector overlay atom 1"
                    ),
                    point2=gui_overlay_point(
                        screen_axis_two_atom_vars[1].get(), "axis-vector overlay atom 2"
                    ),
                )
            elif vector_source == "Normal to two vectors":
                axis = ScreeningAxis(
                    source="geometric",
                    point=axis_point,
                    vector_source="normal_vectors",
                    normal_vector1=gui_xyz(
                        screen_axis_normal_vars[:3], "normal-defining vector 1"
                    ),
                    normal_vector2=gui_xyz(
                        screen_axis_normal_vars[3:], "normal-defining vector 2"
                    ),
                )
            else:
                raise ValueError("Choose a valid geometric rotation-axis vector source.")
        else:
            raise ValueError("Choose a valid rotation-axis source.")

        return ScreeningRequest(
            mode=mode,
            target=target,
            point1=point1,
            point2=point2,
            axis=axis,
        )

    def overlay_mapping_summary() -> str:
        input_pdb = input_var.get().strip()
        if not input_pdb:
            return "Choose an input PDB to display its origin-overlay chain mapping."
        try:
            _records, residues = read_pdb(input_pdb)
            chains: List[str] = []
            for residue in residues.values():
                if residue.p_atom is not None and residue.chain_id not in chains:
                    chains.append(residue.chain_id)
            if len(chains) != 2:
                return f"The selected PDB has {len(chains)} P-containing chains; Bend Helix requires 2."
            model1, model2 = build_origin_overlay_chain_maps(chains)
            original = ", ".join(
                f"{overlay} ← input {source}" for source, overlay in model1.items()
            )
            transformed = ", ".join(
                f"{overlay} ← input {source}" for source, overlay in model2.items()
            )
            return (
                f"Origin overlay: original model [{original}]; "
                f"transformed model [{transformed}]."
            )
        except Exception as exc:
            return f"Origin-overlay mapping is unavailable: {exc}"

    def refresh_screen_angle_state(*_args) -> None:
        selected = selected_screen_angles()
        for name in ("phi", "beta", "tau"):
            angle_entries[name].configure(
                state="disabled" if name in selected else "normal"
            )
            screen_checkbuttons[name].configure(
                state="disabled"
                if len(selected) >= 2 and name not in selected
                else "normal"
            )
        screening_button.configure(
            state="normal" if len(selected) in (1, 2) else "disabled"
        )
        screen_dialog_state["configured"] = False
        if selected:
            screen_status_var.set(
                "Screening settings need review for " + ", ".join(selected) + "."
            )
        else:
            screen_status_var.set(
                "Check one or two angles, then configure the screening target."
            )

    def open_screening_dialog() -> None:
        selected = selected_screen_angles()
        if len(selected) not in (1, 2):
            messagebox.showerror(
                APP_TITLE,
                "Check Screen for exactly one or two angles.",
                parent=root,
            )
            return

        existing = screen_dialog_state.get("window")
        if existing is not None:
            try:
                existing.deiconify()
                existing.lift()
                existing.focus_force()
                return
            except Exception:
                screen_dialog_state["window"] = None

        dialog = tk.Toplevel(root)
        screen_dialog_state["window"] = dialog
        persisted_screen_vars = [
            screen_mode_var,
            screen_target_var,
            screen_point1_atom_var,
            screen_point2_source_var,
            screen_point2_atom_var,
            *screen_point2_xyz_vars,
            screen_axis_source_var,
            screen_axis_point_source_var,
            *screen_axis_point_xyz_vars,
            screen_axis_point_atom_var,
            screen_axis_vector_source_var,
            *screen_axis_vector_vars,
            *screen_axis_two_xyz_vars,
            *screen_axis_two_atom_vars,
            *screen_axis_normal_vars,
        ]
        for range_values in screen_range_vars.values():
            persisted_screen_vars.extend(range_values.values())
        screen_setting_snapshot = [
            (variable, variable.get()) for variable in persisted_screen_vars
        ]
        dialog.title(f"{APP_TITLE} - Screening to achieve...")
        apply_optional_icon(dialog, __file__)
        dialog.geometry("980x760")
        dialog.minsize(900, 650)
        dialog.transient(root)
        dialog.grab_set()
        dialog.columnconfigure(0, weight=1)
        dialog.rowconfigure(0, weight=1)

        def screening_help_button(parent, key: str):
            return tk.Button(
                parent,
                text="?",
                command=lambda: messagebox.showinfo(
                    f"{APP_TITLE} - Screening help",
                    SCREENING_GUI_HELP[key],
                    parent=dialog,
                ),
                bg="#d9ecff",
                activebackground="#c4e0ff",
                highlightbackground="#d9ecff",
                width=2,
                relief="raised",
            )

        outer = ttk.Frame(dialog, padding=10)
        outer.grid(row=0, column=0, sticky="nsew")
        outer.columnconfigure(0, weight=1)
        outer.rowconfigure(3, weight=1)

        ttk.Label(
            outer,
            text="Screening " + " and ".join(selected),
            font=("TkDefaultFont", 11, "bold"),
        ).grid(row=0, column=0, sticky="w")
        ttk.Label(
            outer,
            text=overlay_mapping_summary(),
            wraplength=920,
        ).grid(row=1, column=0, sticky="ew", pady=(3, 8))

        ranges_box = ttk.LabelFrame(outer, text="Coarse candidate angle grid", padding=8)
        ranges_box.grid(row=2, column=0, sticky="ew", pady=(0, 8))
        ttk.Label(ranges_box, text="Angle").grid(
            row=0, column=0, sticky="w", padx=(0, 8)
        )
        for column, label in (
            (1, "From (deg)"),
            (3, "To (deg)"),
            (5, "Step (deg)"),
        ):
            ttk.Label(ranges_box, text=label).grid(
                row=0, column=column, sticky="w", padx=(0, 8)
            )
        for row_index, name in enumerate(selected, start=1):
            ttk.Label(ranges_box, text=name.capitalize()).grid(
                row=row_index, column=0, sticky="w", padx=(0, 8), pady=2
            )
            values = screen_range_vars[name]
            for column, key, help_key in (
                (1, "start", "grid_from"),
                (3, "stop", "grid_to"),
                (5, "step", "grid_step"),
            ):
                ttk.Entry(
                    ranges_box, textvariable=values[key], width=14
                ).grid(row=row_index, column=column, sticky="w", padx=(0, 8), pady=2)
                screening_help_button(ranges_box, help_key).grid(
                    row=row_index,
                    column=column + 1,
                    sticky="w",
                    padx=(0, 10),
                    pady=2,
                )
        ttk.Label(
            ranges_box,
            text="From, To, and Step are degrees. Both endpoints are tested; Step is a "
            "positive magnitude. The strongest coarse results are adaptively refined "
            "between steps to 0.001°. Descending ranges are allowed. Maximum coarse grid: "
            "250,000 candidates.",
            wraplength=840,
        ).grid(
            row=len(selected) + 1,
            column=0,
            columnspan=7,
            sticky="w",
            pady=(5, 0),
        )
        screen_grid_preview_var = tk.StringVar()
        ttk.Label(
            ranges_box,
            textvariable=screen_grid_preview_var,
            wraplength=840,
            foreground="#335f85",
        ).grid(
            row=len(selected) + 2,
            column=0,
            columnspan=7,
            sticky="w",
            pady=(6, 0),
        )

        def refresh_grid_preview(*_args) -> None:
            try:
                preview = format_screening_grid_preview(screening_ranges_from_gui())
            except Exception as exc:
                preview = f"Calculated grid preview unavailable: {exc}"
            screen_grid_preview_var.set(preview)

        grid_trace_handles = []
        for name in selected:
            for variable in screen_range_vars[name].values():
                trace_id = variable.trace_add("write", refresh_grid_preview)
                grid_trace_handles.append((variable, trace_id))
        refresh_grid_preview()

        target_box = ttk.LabelFrame(outer, text="Target", padding=8)
        target_box.grid(row=3, column=0, sticky="nsew", pady=(0, 8))
        target_box.columnconfigure(1, weight=1)

        ttk.Label(target_box, text="Mode").grid(row=0, column=0, sticky="w", pady=2)
        mode_combo = ttk.Combobox(
            target_box,
            textvariable=screen_mode_var,
            values=("Screening for distance", "Screening for rotation"),
            state="readonly",
            width=28,
        )
        mode_combo.grid(row=0, column=1, sticky="w", padx=6, pady=2)
        screening_help_button(target_box, "mode").grid(
            row=0, column=2, sticky="w", padx=6, pady=2
        )

        target_label = ttk.Label(target_box, text="Target distance (angstroms)")
        target_label.grid(row=1, column=0, sticky="w", pady=2)
        ttk.Entry(target_box, textvariable=screen_target_var, width=18).grid(
            row=1, column=1, sticky="w", padx=6, pady=2
        )
        screening_help_button(target_box, "target").grid(
            row=1, column=2, sticky="w", padx=6, pady=2
        )

        ttk.Label(target_box, text="Endpoint 1 overlay atom").grid(
            row=2, column=0, sticky="w", pady=2
        )
        ttk.Entry(target_box, textvariable=screen_point1_atom_var, width=34).grid(
            row=2, column=1, sticky="w", padx=6, pady=2
        )
        ttk.Label(
            target_box,
            text="Use origin-overlay identity CHAIN:RESIDUE:ATOM, e.g. A:36:P or C:36:O5'.",
        ).grid(row=2, column=2, sticky="w", padx=6, pady=2)
        screening_help_button(target_box, "endpoint1_atom").grid(
            row=2, column=3, sticky="w", padx=6, pady=2
        )

        ttk.Label(target_box, text="Endpoint 2 source").grid(
            row=3, column=0, sticky="w", pady=2
        )
        point2_source_combo = ttk.Combobox(
            target_box,
            textvariable=screen_point2_source_var,
            state="readonly",
            width=28,
        )
        point2_source_combo.grid(row=3, column=1, sticky="w", padx=6, pady=2)
        screening_help_button(target_box, "endpoint2_source").grid(
            row=3, column=2, sticky="w", padx=6, pady=2
        )

        point2_host = ttk.Frame(target_box)
        point2_host.grid(row=4, column=0, columnspan=4, sticky="ew", pady=(2, 6))
        point2_atom_frame = ttk.Frame(point2_host)
        ttk.Label(point2_atom_frame, text="Endpoint 2 overlay atom").pack(side="left")
        ttk.Entry(
            point2_atom_frame, textvariable=screen_point2_atom_var, width=34
        ).pack(side="left", padx=6)
        screening_help_button(point2_atom_frame, "endpoint2_atom").pack(side="left")

        point2_xyz_frame = ttk.Frame(point2_host)
        ttk.Label(point2_xyz_frame, text="Endpoint 2 XYZ").pack(side="left")
        for label, variable in zip(("x", "y", "z"), screen_point2_xyz_vars):
            ttk.Label(point2_xyz_frame, text=label).pack(side="left", padx=(6, 0))
            ttk.Entry(point2_xyz_frame, textvariable=variable, width=12).pack(
                side="left", padx=(2, 0)
            )
        screening_help_button(point2_xyz_frame, "endpoint2_xyz").pack(
            side="left", padx=(8, 0)
        )

        point2_pivot_frame = ttk.Frame(point2_host)
        ttk.Label(
            point2_pivot_frame,
            text="Endpoint 2 is the candidate-dependent phi-corrected pivot P position.",
        ).pack(side="left")
        screening_help_button(point2_pivot_frame, "endpoint2_pivot").pack(
            side="left", padx=(8, 0)
        )

        axis_box = ttk.LabelFrame(target_box, text="Rotation axis", padding=8)
        axis_box.grid(row=5, column=0, columnspan=4, sticky="ew", pady=(5, 0))
        axis_box.columnconfigure(1, weight=1)
        ttk.Label(axis_box, text="Axis source").grid(row=0, column=0, sticky="w", pady=2)
        axis_source_combo = ttk.Combobox(
            axis_box,
            textvariable=screen_axis_source_var,
            values=("Geometric definition", "Local axis range(s)"),
            state="readonly",
            width=28,
        )
        axis_source_combo.grid(row=0, column=1, sticky="w", padx=6, pady=2)
        ttk.Label(
            axis_box,
            text="Local axis range(s) uses the range selected in the main window.",
        ).grid(row=0, column=2, sticky="w", padx=6, pady=2)
        screening_help_button(axis_box, "axis_source").grid(
            row=0, column=3, sticky="w", padx=6, pady=2
        )

        geometric_axis_frame = ttk.Frame(axis_box)
        geometric_axis_frame.grid(row=1, column=0, columnspan=4, sticky="ew", pady=(3, 0))
        geometric_axis_frame.columnconfigure(1, weight=1)

        ttk.Label(geometric_axis_frame, text="Axis point source").grid(
            row=0, column=0, sticky="w", pady=2
        )
        axis_point_source_combo = ttk.Combobox(
            geometric_axis_frame,
            textvariable=screen_axis_point_source_var,
            values=("XYZ point", "Overlay atom"),
            state="readonly",
            width=22,
        )
        axis_point_source_combo.grid(row=0, column=1, sticky="w", padx=6, pady=2)
        screening_help_button(geometric_axis_frame, "axis_point_source").grid(
            row=0, column=2, sticky="w", padx=6, pady=2
        )

        axis_point_host = ttk.Frame(geometric_axis_frame)
        axis_point_host.grid(row=1, column=0, columnspan=4, sticky="ew", pady=2)
        axis_point_xyz_frame = ttk.Frame(axis_point_host)
        ttk.Label(axis_point_xyz_frame, text="Axis point XYZ").pack(side="left")
        for label, variable in zip(("x", "y", "z"), screen_axis_point_xyz_vars):
            ttk.Label(axis_point_xyz_frame, text=label).pack(side="left", padx=(6, 0))
            ttk.Entry(axis_point_xyz_frame, textvariable=variable, width=12).pack(
                side="left", padx=(2, 0)
            )
        screening_help_button(axis_point_xyz_frame, "axis_point_xyz").pack(
            side="left", padx=(8, 0)
        )
        axis_point_atom_frame = ttk.Frame(axis_point_host)
        ttk.Label(axis_point_atom_frame, text="Axis point overlay atom").pack(side="left")
        ttk.Entry(
            axis_point_atom_frame, textvariable=screen_axis_point_atom_var, width=34
        ).pack(side="left", padx=6)
        screening_help_button(axis_point_atom_frame, "axis_point_atom").pack(
            side="left"
        )

        ttk.Label(geometric_axis_frame, text="Axis vector source").grid(
            row=2, column=0, sticky="w", pady=2
        )
        axis_vector_source_combo = ttk.Combobox(
            geometric_axis_frame,
            textvariable=screen_axis_vector_source_var,
            values=(
                "Direct vector",
                "Two XYZ points",
                "Two overlay atoms",
                "Normal to two vectors",
            ),
            state="readonly",
            width=28,
        )
        axis_vector_source_combo.grid(row=2, column=1, sticky="w", padx=6, pady=2)
        screening_help_button(geometric_axis_frame, "axis_vector_source").grid(
            row=2, column=2, sticky="w", padx=6, pady=2
        )

        axis_vector_host = ttk.Frame(geometric_axis_frame)
        axis_vector_host.grid(row=3, column=0, columnspan=4, sticky="ew", pady=2)
        direct_vector_frame = ttk.Frame(axis_vector_host)
        ttk.Label(direct_vector_frame, text="Axis vector").pack(side="left")
        for label, variable in zip(("x", "y", "z"), screen_axis_vector_vars):
            ttk.Label(direct_vector_frame, text=label).pack(side="left", padx=(6, 0))
            ttk.Entry(direct_vector_frame, textvariable=variable, width=12).pack(
                side="left", padx=(2, 0)
            )
        screening_help_button(direct_vector_frame, "direct_vector").pack(
            side="left", padx=(8, 0)
        )

        two_xyz_frame = ttk.Frame(axis_vector_host)
        for point_index in range(2):
            ttk.Label(two_xyz_frame, text=f"Point {point_index + 1}").pack(
                side="left", padx=(0 if point_index == 0 else 12, 0)
            )
            offset = point_index * 3
            for coordinate_index, label in enumerate(("x", "y", "z")):
                ttk.Label(two_xyz_frame, text=label).pack(side="left", padx=(5, 0))
                ttk.Entry(
                    two_xyz_frame,
                    textvariable=screen_axis_two_xyz_vars[offset + coordinate_index],
                    width=10,
                ).pack(side="left", padx=(2, 0))
        screening_help_button(two_xyz_frame, "two_xyz_points").pack(
            side="left", padx=(8, 0)
        )

        two_atoms_frame = ttk.Frame(axis_vector_host)
        for atom_index in range(2):
            ttk.Label(two_atoms_frame, text=f"Overlay atom {atom_index + 1}").pack(
                side="left", padx=(0 if atom_index == 0 else 12, 0)
            )
            ttk.Entry(
                two_atoms_frame,
                textvariable=screen_axis_two_atom_vars[atom_index],
                width=28,
            ).pack(side="left", padx=5)
        screening_help_button(two_atoms_frame, "two_overlay_atoms").pack(
            side="left", padx=(8, 0)
        )

        normal_vectors_frame = ttk.Frame(axis_vector_host)
        for vector_index in range(2):
            ttk.Label(normal_vectors_frame, text=f"Vector {vector_index + 1}").pack(
                side="left", padx=(0 if vector_index == 0 else 12, 0)
            )
            offset = vector_index * 3
            for coordinate_index, label in enumerate(("x", "y", "z")):
                ttk.Label(normal_vectors_frame, text=label).pack(side="left", padx=(5, 0))
                ttk.Entry(
                    normal_vectors_frame,
                    textvariable=screen_axis_normal_vars[offset + coordinate_index],
                    width=9,
                ).pack(side="left", padx=(2, 0))
        ttk.Label(normal_vectors_frame, text="axis = vector 1 × vector 2").pack(
            side="left", padx=(12, 0)
        )
        screening_help_button(normal_vectors_frame, "normal_vectors").pack(
            side="left", padx=(8, 0)
        )

        screening_note_label = ttk.Label(
            outer,
            text=(
                "Rotation is measured from endpoint 1 toward endpoint 2 by the right-hand "
                "rule around +axis and compared circularly across ±180°. Screening uses "
                "the origin-overlay geometry and automatically writes the winning -ori PDB."
            ),
            wraplength=920,
        )
        screening_note_label.grid(row=4, column=0, sticky="ew", pady=(0, 8))

        button_row = ttk.Frame(outer)
        button_row.grid(row=5, column=0, sticky="e")

        def close_dialog(restore_snapshot: bool = False) -> None:
            for variable, trace_id in grid_trace_handles:
                variable.trace_remove("write", trace_id)
            if restore_snapshot:
                for variable, saved_value in screen_setting_snapshot:
                    variable.set(saved_value)
            try:
                dialog.grab_release()
            except Exception:
                pass
            screen_dialog_state["window"] = None
            dialog.destroy()

        def save_settings() -> None:
            try:
                ranges = screening_ranges_from_gui()
                screening_request_from_gui()
                _normalized, _values, candidate_count = validate_screen_angle_ranges(
                    ranges, candidate_cap=250000
                )
            except Exception as exc:
                messagebox.showerror(APP_TITLE, str(exc), parent=dialog)
                return
            screen_dialog_state["configured"] = True
            mode_short = (
                "rotation"
                if screen_mode_var.get().strip() == "Screening for rotation"
                else "distance"
            )
            screen_status_var.set(
                f"Configured {mode_short} screening: {candidate_count} coarse candidate(s) "
                "plus adaptive refinement."
            )
            close_dialog(restore_snapshot=False)

        ttk.Button(button_row, text="Save settings", command=save_settings).pack(
            side="left", padx=(0, 6)
        )
        ttk.Button(
            button_row,
            text="Cancel",
            command=lambda: close_dialog(restore_snapshot=True),
        ).pack(side="left")

        point2_frames = {
            "Overlay atom": point2_atom_frame,
            "XYZ point": point2_xyz_frame,
            "Phi-corrected pivot P": point2_pivot_frame,
        }
        axis_point_frames = {
            "XYZ point": axis_point_xyz_frame,
            "Overlay atom": axis_point_atom_frame,
        }
        axis_vector_frames = {
            "Direct vector": direct_vector_frame,
            "Two XYZ points": two_xyz_frame,
            "Two overlay atoms": two_atoms_frame,
            "Normal to two vectors": normal_vectors_frame,
        }

        def refresh_dialog_fields(_event=None) -> None:
            is_rotation = screen_mode_var.get().strip() == "Screening for rotation"
            target_label.configure(
                text="Target rotation (degrees)" if is_rotation else "Target distance (angstroms)"
            )
            screening_note_label.configure(
                text=(
                    "Rotation is measured from endpoint 1 toward endpoint 2 by the right-hand "
                    "rule around +axis and compared circularly across ±180°. Screening uses "
                    "the origin-overlay geometry and automatically writes the winning -ori PDB."
                    if is_rotation
                    else
                    "Distance is measured directly between endpoint 1 and endpoint 2 in the "
                    "origin-overlay geometry. Screening keeps the closest coarse-to-fine "
                    "candidate and "
                    "automatically writes the winning -ori PDB."
                )
            )
            allowed_point2 = ["Overlay atom", "XYZ point"]
            if is_rotation:
                allowed_point2.append("Phi-corrected pivot P")
            point2_source_combo.configure(values=tuple(allowed_point2))
            if screen_point2_source_var.get() not in allowed_point2:
                screen_point2_source_var.set("Overlay atom")

            for frame in point2_frames.values():
                frame.pack_forget()
            point2_frames[screen_point2_source_var.get()].pack(fill="x")

            if is_rotation:
                axis_box.grid()
            else:
                axis_box.grid_remove()

            if screen_axis_source_var.get().strip() == "Geometric definition":
                geometric_axis_frame.grid()
            else:
                geometric_axis_frame.grid_remove()

            for frame in axis_point_frames.values():
                frame.pack_forget()
            axis_point_frames[screen_axis_point_source_var.get()].pack(fill="x")

            for frame in axis_vector_frames.values():
                frame.pack_forget()
            axis_vector_frames[screen_axis_vector_source_var.get()].pack(fill="x")
            dialog.update_idletasks()

        for combo in (
            mode_combo,
            point2_source_combo,
            axis_source_combo,
            axis_point_source_combo,
            axis_vector_source_combo,
        ):
            combo.bind("<<ComboboxSelected>>", refresh_dialog_fields)
        dialog.protocol(
            "WM_DELETE_WINDOW", lambda: close_dialog(restore_snapshot=True)
        )
        refresh_dialog_fields()
        mode_combo.focus_set()

    screening_button.configure(command=open_screening_dialog)
    for variable in screen_angle_vars.values():
        variable.trace_add("write", refresh_screen_angle_state)
    refresh_screen_angle_state()

    def run_from_gui() -> None:
        cli_cmd = ""
        screening_summary = ""
        try:
            input_pdb = input_var.get().strip()
            output_pdb = normalize_output_path(output_var.get())
            pivot_residue = pivot_var.get().strip()
            if not input_pdb:
                raise ValueError("Please provide an input PDB file.")
            if not pivot_residue:
                raise ValueError("Please provide the pivot P residue.")
            sep_mode = normalize_sep(sep_var.get())
            align_mode = normalize_align(align_var.get())
            origin_mode = normalize_origin(origin_var.get())
            axis_range_specs = get_axis_range_specs_from_gui()

            screened_names = selected_screen_angles()
            if screened_names:
                if not bool(screen_dialog_state.get("configured")):
                    raise ValueError(
                        "Click Screening to achieve..., review the target settings, "
                        "and choose Save settings before running."
                    )
                ranges = screening_ranges_from_gui()
                request = screening_request_from_gui()
                fixed_angles = {
                    name: gui_finite_float(
                        angle_vars[name].get(), f"fixed {name} angle"
                    )
                    for name in ("phi", "beta", "tau")
                    if name not in screened_names
                }
                _normalized, _values, candidate_count = validate_screen_angle_ranges(
                    ranges, candidate_cap=250000
                )
                result_text.delete("1.0", tk.END)
                result_text.insert(
                    "1.0",
                    f"Screening {candidate_count} coarse candidate(s), then refining "
                    f"between steps for {request.mode}...\n",
                )
                result_text.see(tk.END)
                root.update_idletasks()
                print(
                    f"[Bend Helix screening] Evaluating {candidate_count} coarse "
                    f"candidate(s), then refining between steps for {request.mode}.",
                    flush=True,
                )
                context = prepare_screening_context(
                    input_pdb=input_pdb,
                    pivot_residue=pivot_residue,
                    axis_range_specs=axis_range_specs,
                )
                screen_result = screen_bend_angles(
                    context=context,
                    fixed_angles=fixed_angles,
                    ranges=ranges,
                    request=request,
                    align_mode=align_mode,
                    candidate_cap=250000,
                )
                for name, value in screen_result.angles.items():
                    angle_vars[name].set(format_float_for_cli(value))
                phi_deg = screen_result.phi_deg
                beta_deg = screen_result.beta_deg
                tau_deg = screen_result.tau_deg
                origin_mode = "y"
                origin_var.set("y")
                if output_pdb is None:
                    output_pdb = make_output_name(
                        input_pdb,
                        phi_deg,
                        beta_deg,
                        tau_deg,
                        sep_mode=sep_mode,
                        screen_mode=True,
                    )
                unit = "A" if request.mode == "distance" else "deg"
                range_lines = [
                    (
                        f"  {item.name}: {format_float_for_cli(item.start)} to "
                        f"{format_float_for_cli(item.stop)} by "
                        f"{format_float_for_cli(item.step)} deg"
                    )
                    for item in ranges
                ]
                screening_summary = "\n".join(
                    [
                        f"Screening mode: {request.mode}",
                        "Screened coarse angle grids (degrees):",
                        *range_lines,
                        (
                            "Coarse candidates evaluated: "
                            f"{screen_result.candidate_count}"
                        ),
                        (
                            "Refinement candidates evaluated: "
                            f"{screen_result.refinement_candidate_count}"
                        ),
                        (
                            "Total candidates evaluated: "
                            f"{screen_result.evaluated_candidate_count}"
                        ),
                        (
                            "Refinement precision: "
                            f"{SCREEN_REFINEMENT_TOLERANCE_DEG:g} deg"
                        ),
                        f"Target: {request.target:.9g} {unit}",
                        f"Achieved: {screen_result.achieved_value:.9g} {unit}",
                        f"Absolute residual: {screen_result.error:.9g} {unit}",
                        (
                            "Selected angles: "
                            f"phi={format_float_for_cli(phi_deg)}, "
                            f"beta={format_float_for_cli(beta_deg)}, "
                            f"tau={format_float_for_cli(tau_deg)} deg"
                        ),
                        "Origin overlay output was enabled automatically.",
                    ]
                )
                print(screening_summary, flush=True)
            else:
                phi_deg = gui_finite_float(phi_var.get(), "phi angle")
                beta_deg = gui_finite_float(beta_var.get(), "beta angle")
                tau_deg = gui_finite_float(tau_var.get(), "tau angle")

            cli_cmd = build_equivalent_cli_command(
                input_pdb=input_pdb,
                pivot_residue=pivot_residue,
                phi_deg=phi_deg,
                beta_deg=beta_deg,
                tau_deg=tau_deg,
                sep_mode=sep_mode,
                align_mode=align_mode,
                origin_mode=origin_mode,
                output_pdb=output_pdb,
                axis_range_specs=axis_range_specs,
            )

            result_text.delete("1.0", tk.END)
            screening_prefix = f"{screening_summary}\n\n" if screening_summary else ""
            result_text.insert(
                "1.0",
                f"{screening_prefix}CLI: {cli_cmd}\n"
                "Piece #1 is fixed; piece #2 is movable.\nRunning...\n",
            )
            result_text.see(tk.END)
            root.update_idletasks()
            print("Equivalent CLI command:", flush=True)
            print(cli_cmd, flush=True)
            print("Piece #1 is fixed; piece #2 is movable.", flush=True)

            out_path, info = run_bending(
                input_pdb=input_pdb,
                pivot_residue=pivot_residue,
                phi_deg=phi_deg,
                beta_deg=beta_deg,
                tau_deg=tau_deg,
                sep_mode=sep_mode,
                align_mode=align_mode,
                origin_mode=origin_mode,
                output_pdb=output_pdb,
                axis_range_specs=axis_range_specs,
            )
            summary = format_run_summary(out_path, info)
        except Exception as exc:
            if cli_cmd:
                print("Equivalent CLI command:", flush=True)
                print(cli_cmd, flush=True)
            print(f"Error: {exc}", flush=True)
            if cli_cmd:
                result_text.delete("1.0", tk.END)
                result_text.insert("1.0", f"CLI: {cli_cmd}\n\nError: {exc}\n")
            else:
                result_text.delete("1.0", tk.END)
                result_text.insert("1.0", f"Error: {exc}\n")
            messagebox.showerror(APP_TITLE, str(exc), parent=root)
            return

        result_text.delete("1.0", tk.END)
        screening_prefix = f"{screening_summary}\n\n" if screening_summary else ""
        result_text.insert("1.0", f"{screening_prefix}CLI: {cli_cmd}\n\n{summary}")
        result_text.see(tk.END)
        print(summary, flush=True)
        messagebox.showinfo(APP_TITLE, f"Wrote {out_path}", parent=root)

    ttk.Button(button_frame, text="Run", command=run_from_gui).grid(row=0, column=0, padx=(0, 6))
    ttk.Button(button_frame, text="Close", command=root.destroy).grid(row=0, column=1)

    input_entry.focus_set()
    root.mainloop()
    return 0


def main() -> int:
    parser = build_arg_parser()
    args = parser.parse_args()

    if len(sys.argv) == 1 or args.gui:
        defaults = {
            "input_pdb": args.input_pdb_opt or args.input_pdb or "",
            "output_pdb": args.output_pdb or "",
            "pivot_residue": args.pivot_residue_opt or args.pivot_residue or "",
            "phi_deg": "" if (args.phi_deg_opt is None and args.phi_deg is None) else str(args.phi_deg_opt if args.phi_deg_opt is not None else args.phi_deg),
            "beta_deg": "" if (args.beta_deg_opt is None and args.beta_deg is None) else str(args.beta_deg_opt if args.beta_deg_opt is not None else args.beta_deg),
            "tau_deg": "0" if args.tau_deg_opt is None else str(args.tau_deg_opt),
            "sep": args.sep,
            "align": args.align,
            "origin": args.origin,
            "axis_ranges": "\n".join(args.axis_ranges or []),
        }
        return launch_gui(defaults)

    try:
        (
            input_pdb,
            pivot_residue,
            phi_deg,
            beta_deg,
            tau_deg,
            sep_mode,
            align_mode,
            origin_mode,
            output_pdb,
            axis_range_specs,
        ) = resolve_run_parameters(args)
    except Exception as exc:
        parser.print_usage(sys.stderr)
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    try:
        out_path, info = run_bending(
            input_pdb=input_pdb,
            pivot_residue=pivot_residue,
            phi_deg=phi_deg,
            beta_deg=beta_deg,
            tau_deg=tau_deg,
            sep_mode=sep_mode,
            align_mode=align_mode,
            origin_mode=origin_mode,
            output_pdb=output_pdb,
            axis_range_specs=axis_range_specs,
        )
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    print(format_run_summary(out_path, info))
    return 0



if __name__ == "__main__":
    raise SystemExit(main())
