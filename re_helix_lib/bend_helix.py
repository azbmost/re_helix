#!/usr/bin/env python3
"""
bend_helixV2_4.py

Bend a two-chain nucleic-acid helix at a user-selected phosphorus site.

Usage modes
-----------
1) Positional CLI arguments (backward-compatible):
       bend_helixV2_4.py input.pdb A36 0 30
       # tau defaults to 0 degrees in positional mode.
       # align defaults to y unless you set --align n.

2) Named CLI arguments:
       bend_helixV2_4.py --input A60-heli.pdb --pivot A36 --phi 0 --beta 30 --tau 10 --align y --origin y
       bend_helixV2_4.py --input A60-heli.pdb --pivot A36 --phi 0 --beta 30 --axis_range A1-A35,B60-B26 -o bent.pdb

3) GUI mode:
       bend_helixV2_4.py
       bend_helixV2_4.py --gui

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
import math
import os
import shlex
import string
import sys
from collections import OrderedDict
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Tuple

EPS = 1.0e-8
CHAIN_ID_CANDIDATES = string.ascii_uppercase + string.ascii_lowercase + string.digits


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
    script_name = os.path.basename(__file__) if "__file__" in globals() else "bend_helixV2_4.py"
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

def make_output_name(inp: str, phi_deg: float, beta_deg: float, tau_deg: float, sep_mode: str = "n") -> str:
    stem, ext = os.path.splitext(inp)
    if not ext:
        ext = ".pdb"
    suffix = (
        f"_P{format_angle_for_filename(phi_deg)}"
        f"B{format_angle_for_filename(beta_deg)}"
        f"T{format_angle_for_filename(tau_deg)}"
    )
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
            "bend_helixV2_4.py currently expects exactly two chains that contain P atoms; "
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
    duplex = build_duplex_data(residues, selected_key=selected_key, axis_range_defs=axis_range_defs)
    axis_point = duplex["axis_point"]
    axis_dir = duplex["axis_dir"]
    pairs = duplex["pairs"]

    if selected_key not in residues:
        raise ValueError(f"Residue {selected_key[0]}{selected_key[1]} was not found in the input PDB.")
    if residues[selected_key].p_atom is None:
        raise ValueError(f"Residue {selected_key[0]}{selected_key[1]} does not contain a P atom.")

    pair_idx = find_pair_index(pairs, selected_key)
    piece2_keys = set()
    for r1, r2 in pairs[pair_idx:]:
        piece2_keys.add((r1.chain_id, r1.res_seq))
        piece2_keys.add((r2.chain_id, r2.res_seq))

    selected_p = residues[selected_key].p_coord()
    axis_foot = project_point_to_line(selected_p, axis_point, axis_dir)
    radial = v_sub(selected_p, axis_foot)
    radius = v_len(radial)
    if radius < EPS:
        raise ValueError(
            f"Selected P atom at {selected_key[0]}{selected_key[1]} lies too close to the estimated helix axis."
        )

    phi_rad = math.radians(phi_deg)
    beta_rad = math.radians(beta_deg)
    tau_rad = math.radians(tau_deg)
    radial_phi = rotate_vector(radial, axis_dir, phi_rad)
    hinge_point = v_add(axis_foot, radial_phi)
    hinge_dir = v_cross(axis_dir, radial_phi)
    if v_len(hinge_dir) < EPS:
        raise ValueError("Failed to construct a non-zero tangent direction for the hinge.")
    hinge_dir = v_norm(hinge_dir)

    twist_axis_point_pre_align = rotate_point_about_line(axis_foot, hinge_point, hinge_dir, beta_rad)
    twist_axis_dir = v_norm(rotate_vector(axis_dir, hinge_dir, beta_rad))

    def transform_piece2_coord(coord: Tuple[float, float, float]) -> Tuple[float, float, float]:
        new_coord = rotate_point_about_line(coord, hinge_point, hinge_dir, beta_rad)
        if abs(tau_rad) > 0.0:
            new_coord = rotate_point_about_line(new_coord, twist_axis_point_pre_align, twist_axis_dir, tau_rad)
        return new_coord

    align_translation = (0.0, 0.0, 0.0)
    if align_mode == "y":
        transformed_selected_p = transform_piece2_coord(selected_p)
        align_translation = v_sub(selected_p, transformed_selected_p)

    twist_axis_point = v_add(twist_axis_point_pre_align, align_translation)

    for rec in records:
        if isinstance(rec, AtomRecord):
            key = (rec.chain_id, rec.res_seq)
            if key in piece2_keys:
                new_coord = transform_piece2_coord(rec.coord())
                if align_mode == "y":
                    new_coord = v_add(new_coord, align_translation)
                rec.set_coord(new_coord)

    piece2_chain_map = None
    if sep_mode == "y":
        piece2_chain_map = separate_piece2_chains(records, piece2_keys, duplex["chains"])

    out_path = normalize_output_path(output_pdb) or make_output_name(input_pdb, phi_deg, beta_deg, tau_deg, sep_mode=sep_mode)
    info = {
        "pair_idx": pair_idx,
        "n_pairs": len(pairs),
        "piece2_pair_start": pair_idx + 1,
        "piece2_pair_end": len(pairs),
        "radius": radius,
        "axis_point": axis_point,
        "axis_dir": axis_dir,
        "hinge_point": hinge_point,
        "hinge_dir": hinge_dir,
        "beta_deg": beta_deg,
        "tau_deg": tau_deg,
        "twist_axis_point_pre_align": twist_axis_point_pre_align,
        "twist_axis_point": twist_axis_point,
        "twist_axis_dir": twist_axis_dir,
        "align_mode": align_mode,
        "align_translation": align_translation,
        "pair": ((pairs[pair_idx][0].chain_id, pairs[pair_idx][0].res_seq), (pairs[pair_idx][1].chain_id, pairs[pair_idx][1].res_seq)),
        "sep_mode": sep_mode,
        "piece2_chain_map": piece2_chain_map,
        "duplex_chains": tuple(duplex["chains"]),
        "axis_source": duplex.get("axis_source"),
        "axis_range_used": format_axis_range_spec(duplex["axis_range_used"]) if duplex.get("axis_range_used") is not None else None,
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

    overlay_ids = choose_new_chain_ids([], len(chain_order) * 2)
    chain_map_model1 = {old: overlay_ids[i] for i, old in enumerate(chain_order)}
    chain_map_model2 = {old: overlay_ids[len(chain_order) + i] for i, old in enumerate(chain_order)}

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

    root.title("bend_helixV2_4")
    root.geometry("940x760")
    root.resizable(True, True)
    root.columnconfigure(1, weight=1)
    root.rowconfigure(11, weight=1)

    input_var = tk.StringVar(value=defaults.get("input_pdb", ""))
    output_var = tk.StringVar(value=defaults.get("output_pdb", ""))
    pivot_var = tk.StringVar(value=defaults.get("pivot_residue", ""))
    phi_var = tk.StringVar(value=defaults.get("phi_deg", ""))
    beta_var = tk.StringVar(value=defaults.get("beta_deg", ""))
    tau_var = tk.StringVar(value=defaults.get("tau_deg", "0"))
    sep_var = tk.StringVar(value=defaults.get("sep", "n"))
    align_var = tk.StringVar(value=defaults.get("align", "y"))
    origin_var = tk.StringVar(value=defaults.get("origin", "n"))

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
            "main bent model."
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
        messagebox.showinfo("bend_helixV2_4 help", help_text[key], parent=root)

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
    ttk.Entry(root, textvariable=phi_var).grid(row=4, column=1, sticky="ew", padx=8, pady=4)
    help_button(4, "phi")

    ttk.Label(root, text="Beta: bend movable piece #2 (degrees)").grid(row=5, column=0, sticky="w", padx=8, pady=4)
    ttk.Entry(root, textvariable=beta_var).grid(row=5, column=1, sticky="ew", padx=8, pady=4)
    help_button(5, "beta")

    ttk.Label(root, text="Tau: twist movable piece #2 (degrees)").grid(row=6, column=0, sticky="w", padx=8, pady=4)
    ttk.Entry(root, textvariable=tau_var).grid(row=6, column=1, sticky="ew", padx=8, pady=4)
    help_button(6, "tau")

    ttk.Label(root, text="Local axis range(s), one per line").grid(row=7, column=0, sticky="nw", padx=8, pady=4)
    axis_range_text = tk.Text(root, height=3, width=48, wrap="none")
    axis_range_text.grid(row=7, column=1, sticky="ew", padx=8, pady=4)
    axis_default = defaults.get("axis_ranges", "")
    if axis_default:
        axis_range_text.insert("1.0", axis_default)
    ttk.Label(root, text="Example: A1-A35,B60-B26").grid(row=7, column=2, sticky="w", padx=(0, 8), pady=4)
    help_button(7, "axis_range")

    ttk.Label(root, text="Separate piece #2 chain IDs (--sep)").grid(row=8, column=0, sticky="w", padx=8, pady=4)
    sep_box = ttk.Combobox(root, textvariable=sep_var, values=("n", "y"), state="readonly", width=6)
    sep_box.grid(row=8, column=1, sticky="w", padx=8, pady=4)
    sep_box.set(sep_var.get())
    help_button(8, "sep")

    ttk.Label(root, text="Realign pivot P after bending (--align)").grid(row=9, column=0, sticky="w", padx=8, pady=4)
    align_box = ttk.Combobox(root, textvariable=align_var, values=("y", "n"), state="readonly", width=6)
    align_box.grid(row=9, column=1, sticky="w", padx=8, pady=4)
    align_box.set(align_var.get())
    help_button(9, "align")

    ttk.Label(root, text="Write origin overlay PDB (--origin)").grid(row=10, column=0, sticky="w", padx=8, pady=4)
    origin_box = ttk.Combobox(root, textvariable=origin_var, values=("n", "y"), state="readonly", width=6)
    origin_box.grid(row=10, column=1, sticky="w", padx=8, pady=4)
    origin_box.set(origin_var.get())
    help_button(10, "origin")

    result_text = scrolledtext.ScrolledText(root, height=14, width=90, wrap="word")
    result_text.grid(row=11, column=0, columnspan=4, sticky="nsew", padx=8, pady=(8, 4))

    button_frame = ttk.Frame(root)
    button_frame.grid(row=12, column=0, columnspan=4, sticky="e", padx=8, pady=(4, 8))

    def get_axis_range_specs_from_gui() -> List[str]:
        text = axis_range_text.get("1.0", tk.END)
        specs: List[str] = []
        for line in text.replace(";", "\n").splitlines():
            line = line.strip()
            if line:
                specs.append(line)
        return specs

    def run_from_gui() -> None:
        cli_cmd = ""
        try:
            input_pdb = input_var.get().strip()
            output_pdb = normalize_output_path(output_var.get())
            pivot_residue = pivot_var.get().strip()
            if not input_pdb:
                raise ValueError("Please provide an input PDB file.")
            if not pivot_residue:
                raise ValueError("Please provide the pivot P residue.")
            phi_deg = float(phi_var.get().strip())
            beta_deg = float(beta_var.get().strip())
            tau_deg = float(tau_var.get().strip())
            sep_mode = normalize_sep(sep_var.get())
            align_mode = normalize_align(align_var.get())
            origin_mode = normalize_origin(origin_var.get())
            axis_range_specs = get_axis_range_specs_from_gui()
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
            result_text.insert("1.0", f"CLI: {cli_cmd}\nPiece #1 is fixed; piece #2 is movable.\nRunning...\n")
            result_text.see(tk.END)
            root.update_idletasks()

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
                result_text.delete("1.0", tk.END)
                result_text.insert("1.0", f"CLI: {cli_cmd}\n\nError: {exc}\n")
            else:
                result_text.delete("1.0", tk.END)
                result_text.insert("1.0", f"Error: {exc}\n")
            messagebox.showerror("bend_helixV2_4", str(exc), parent=root)
            return

        result_text.delete("1.0", tk.END)
        result_text.insert("1.0", f"CLI: {cli_cmd}\n\n{summary}")
        result_text.see(tk.END)
        messagebox.showinfo("bend_helixV2_4", f"Wrote {out_path}", parent=root)

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
