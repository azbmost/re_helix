#!/usr/bin/env python3
"""Build Heterohelix Junction V1.0 (an RE_Helix tool): join duplex segments of any helix form with fitted junctions.

Give it two or more duplex segments (PDB files of any helix form: B-DNA, A-RNA, Z-DNA, with D or L sugars).  It
stacks them in the order given, fits every junction as one base-pair step, joins the strands across it and writes
the joined model with a JSON report.

How each junction is joined: the chain sense rule
    Each end pair gets a chirality-aware frame: x towards the major groove, y towards the C1' of the strand that
    ends there with its 3' end, z = x cross y.  In B-DNA and A-RNA z points out of the segment, along that strand
    (chain sense +1).  In Z-DNA it points back into the segment (-1).  A mirror image flips the sign: L-DNA and L-RNA
    (the mirror images of B-DNA and A-RNA) are -1, and L-Z-DNA (the mirror image of Z-DNA) is +1.
      same chain sense      matched -> continuous: two ordinary 3'-5' links, both strands run on 5'->3'
                            (B|B, Z|Z, L-DNA|Z-DNA, ...)
      opposite chain sense  opposite -> bowtie: a strand that ran on would meet the base pairs from their other face
                            (the chain sense paradox), so the strands turn back instead: one 3'-3' link through an
                            X33 phosphate and one 5'-5' link through one phosphate (B|Z, D-DNA|L-DNA, ...)
    --link can name the kind, but only the kind the chain senses allow; the other is refused before any fitting.

Checks
    Before fitting: the chain sense of every junction end (printed as a "Chain sense check" block); the chain sense
    all along each segment (the z axes of neighbouring base pairs, y towards strand A, may not turn by more than
    90 deg; a segment with a pair turned over is refused, naming the pairs); at least two base pairs per segment;
    and the atoms each link bonds (O3' and C3' at a 3' end, C4', C5' and O5' at a 5' end).
    After building: every output chain is walked, and every paired nucleotide (X33 skipped) gets its pair frame
    with y towards that chain's residue.  The z axes of neighbouring pairs may not turn by more than 90 deg, inside
    segments or across a junction ("chain sense not matched at ..." otherwise).  The cos between each z axis and
    its segment's helix axis, taken in the chain's file order, shows which way the chain runs through its pairs
    (about +1 in every segment of a B-Z bowtie model, about -1 in an L-Z continuous model); it must keep one sign
    along a chain.

How a junction is fitted
    1. The two end pairs get frames whose z axes both point from the first segment into the second.
    2. The second end is placed on the first by one step: tilt, roll, twist (deg), shift, slide, rise (A).
    3. A step is scored on: heavy atoms across the junction closer than 2.9 A; the gaps the two links must span,
       each with a quick phosphate (bowtie: O3'...O3' 2.40-2.60 A and C5'...C5' 3.5-4.8 A; continuous:
       O3'...O5' 2.40-2.60 A, the O5' free to turn about C4'-C5'); base pairs parallel and 3.5 +- 0.15 A apart;
       base-ring overlap (rewarded).
    4. Powell searches start from a grid: twist -180 to 160 deg in 20-deg steps x shift x slide (-2, 0, 2 A), 162
       starts, and every junction is searched from both of its ends.  A second pass starts again from eight small
       fixed nudges around each of the 12 best different results.
    5. The 16 best different results get real phosphates (force-field-like weights: bonds 300 kcal/mol/A^2, angles
       100 kcal/mol/rad^2, contacts 50), minimized in frames attached to the atoms from many starts, so where the files
       sit in space does not matter; the one with the lowest stacking score + phosphate energies is kept.
    Junctions of one model that join the same two end pairs (in the same order, or turned 180 deg about the duplex
    dyad) share one fit: their searches are pooled and they all get the step with the lowest summed total.
    --separate fits each of them on its own.
    The score changes little with twist, so read a fitted twist as uncertain by several degrees; the report lists
    every candidate.

Segments
    Duplexes as BNP_NA writes them: two chains, the first (strand A) 5'->3' from the bottom end to the top end, the
    second (strand B) 5'->3' back down; strand A's i-th residue pairs with strand B's i-th residue from the other
    end.  DNA (DA DC DG DT DU) or RNA (A C G U), D or L sugars (told apart by the sugar geometry; L-forms keep the D
    residue names, as BNP_NA writes them).  Each segment is a rigid copy of its file; only the junction phosphates
    are built, and the O5' atoms they bond are turned about C4'-C5'.
        SEG.pdb        upright: the junction below it meets its bottom end, the one above it its top end
        SEG.pdb:flip   upside down

Run
    python3 build_heterohelix_junction.py L21.pdb Z24.pdb                   (writes L21_heterohelix.pdb)
    python3 build_heterohelix_junction.py B21.pdb Z24.pdb B21.pdb Z24.pdb -o B21_Z24_B21_Z24.pdb
    python3 build_heterohelix_junction.py D.pdb L.pdb:flip -o D_L.pdb          (the top of D meets the top of L)
    python3 build_heterohelix_junction.py --gui [SEG.pdb ...]                   (the tool window)
    --link auto|bowtie|continuous (default auto; one value, or one per junction separated by commas)
    --cores N (default: all but one), --rise A (default 3.5), --separate, --report FILE
    --reuse             take the fits from this output's own report instead of fitting again, when they were made
                        the same way: this tool version, unchanged segment files, the same search grid and the same
                        junctions fitted together (otherwise the junctions are fitted again and the log says why)
    --reuse-from REPORT the same, from an earlier report (for example the same build under another output name)

Output
    OUT.pdb       default <first segment>_heterohelix.pdb beside the first segment.  Chains A, B numbered from 1;
                  REMARK 900 lines record the segments, the junctions and the L-sugar residues; REMARK 950 RE_SCRIPT
                  lines record the tool, the command and the junction residues (Get Phenix Restraints reads them for
                  its _junctions.params); LINK records for every bond that is not an ordinary 5'->3' step in file
                  order (inverted steps, X33 links, 5'-5' links).  As built, not minimized.
    tmp_file/<OUT stem>_heterohelix_junction.json  next to OUT.pdb (for the default output,
                  <first segment>_heterohelix_junction.json): every fit (all searches, the 16 candidates), the
                  phosphates, the helical axes at each junction, the checks, and what the fits were made from.
Identical junctions (same files, same ends) are fitted once.  About 5 minutes per junction on 9 processes (every
junction is searched from both ends, then the second pass; junctions fitted together pool their searches, so each one
still adds its own).
Needs numpy, scipy and shapely.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing as mp
import os
import shlex
import sys
import time
import traceback
from pathlib import Path

TOOL_NAME = "Build Heterohelix Junction"
VERSION = "1.0"
NEEDED_PACKAGES = "numpy, scipy and shapely"


def _stop_for_missing_packages(exc):  # pragma: no cover - depends on the Python that runs the script
    """Run as a script on a Python without numpy, scipy or shapely: one plain line, a message box when the window
    was asked for (tkinter is in the standard library), and exit status 1."""
    text = f"{TOOL_NAME} needs {NEEDED_PACKAGES} ({exc}). Install them, for example: pip install numpy scipy shapely"
    print(f"Error: {text}", file=sys.stderr, flush=True)
    if not sys.argv[1:] or "--gui" in sys.argv[1:]:
        try:
            import tkinter as tk
            from tkinter import messagebox
            root = tk.Tk()
            root.withdraw()
            messagebox.showerror(TOOL_NAME, text, parent=root)
            root.destroy()
        except Exception:  # noqa: BLE001 - no display: the line above says it
            pass
    sys.exit(1)


try:
    import numpy as np
    from scipy.optimize import minimize
    from scipy.spatial import cKDTree
    from shapely.geometry import Polygon
    from shapely.ops import unary_union
except ImportError as _exc:
    if __name__ != "__main__":
        raise
    _stop_for_missing_packages(_exc)

try:
    from re_helix_lib.gui_icon import apply_optional_icon
except ImportError:  # pragma: no cover - direct script execution fallback
    from gui_icon import apply_optional_icon

# Provenance records, the REMARK 950 RE_SCRIPT convention of the other RE_Helix tools (Get Phenix Restraints reads
# the JUNCTION lines).
REMARK_PREFIX = "REMARK 950 RE_SCRIPT"
SOFTWARE_NAME = "build_heterohelix_junction"
SOFTWARE_DEVELOPER = "DiLiuLab"

# ----------------------------------------------------------------------------- constants
# Phosphate geometry as in the RNA gallery's 13b builder and the B-Z bowtie example (Parkinson et al. 1996)
P_O, P_OP, ANG_OPOP = 1.60, 1.48, 119.6             # P-O bridging and non-bridging (A), OP1-P-OP2 (deg)
ANG_O3P, ANG_O5P, ANG_OPO_BR = 119.7, 120.9, 104.0  # C3'-O3'-P, C5'-O5'-P, O-P-O of the bridging O (deg)
K_BOND, K_ANG, K_CLASH = 300.0, 100.0, 50.0          # kcal/mol per A^2, per rad^2, per A^2 below 2.9 A
GAP = (2.40, 2.60)                                   # O...O spanned by one phosphate (O-P-O 97-109 deg)
RISE = 3.5                                           # stacking distance aimed at (A; +-0.15 A is free)
REBUILT_5 = ("P", "OP1", "OP2", "O5'")               # atoms of a 5'-end residue that a junction link rebuilds
TWISTS = tuple(range(-180, 180, 20))
SHIFTS = (-2.0, 0.0, 2.0)
N_CAND = 16                                          # distinct optima per junction that get real phosphates
N_REFINE = 12                                        # best distinct optima of the grid that are searched again
# perturbed restarts around each of them (tilt roll twist deg, shift slide rise A): a fixed table, so every run of the
# same junction does the same thing
REFINE_KICKS = ((3, 3, 6, 0.7, 0.7, 0.15), (-3, -3, -6, -0.7, -0.7, -0.15), (3, -3, -6, 0.7, -0.7, 0.15),
                (-3, 3, 6, -0.7, 0.7, -0.15), (3, 3, -6, -0.7, 0.7, 0.15), (-3, -3, 6, 0.7, -0.7, -0.15),
                (3, -3, 6, -0.7, -0.7, -0.15), (-3, 3, -6, 0.7, 0.7, 0.15))
NWIN = 6                                             # residues per strand near each end that take part in a fit
SAME_END_RMS = 0.10                                  # A: two ends present the same end pair when their windows
                                                     # (matched residue by residue from the end) superpose this well
BACKBONE = {"P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'"}
ORDER = ["P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'"]
BASE = {"DA": "N9 C8 N7 C5 C6 N6 N1 C2 N3 C4".split(), "DG": "N9 C8 N7 C5 C6 O6 N1 C2 N2 N3 C4".split(),
        "DT": "N1 C2 O2 N3 C4 O4 C5 C7 C6".split(), "DC": "N1 C2 O2 N3 C4 N4 C5 C6".split(),
        "DU": "N1 C2 O2 N3 C4 O4 C5 C6".split()}
BASE.update({"A": BASE["DA"], "G": BASE["DG"], "C": BASE["DC"], "U": BASE["DU"]})
_PUR = [["N9", "C8", "N7", "C5", "C4"], ["C4", "C5", "C6", "N1", "C2", "N3"]]
_PYR = [["N1", "C2", "N3", "C4", "C5", "C6"]]
RINGS = {k: (_PUR if k in ("DA", "DG", "A", "G") else _PYR) for k in BASE}
LINK_WORDS = {"bowtie": "bowtie: one 3'-3' link through an X33 phosphate and one 5'-5' link through one phosphate",
              "continuous": "continuous: two ordinary 3'-5' links (both strands run on 5'->3')"}


class JunctionError(SystemExit):
    pass


# ----------------------------------------------------------------------------- small geometry
def read_atoms(path):
    out = []
    for line in open(path):
        if line.startswith(("ATOM", "HETATM")) and not line[12:16].strip().startswith("H"):
            out.append((line[21], int(line[22:26]), line[17:20].strip(), line[12:16].strip(),
                        np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])))
    return out


def bp_frame(ix, rn, a, b):
    """Chirality-aware base-pair frame (origin, rows x y z): x from the C1'-C1' midpoint towards the base centroid
    (the major groove), y from C1'(b) to C1'(a) (a is strand I), z = x cross y."""
    P = np.array([ix[(a[0], a[1], n)] for n in BASE[rn[a]]] + [ix[(b[0], b[1], n)] for n in BASE[rn[b]]])
    o = P.mean(0)
    nrm = np.linalg.svd(P - o)[2][2]
    c1a, c1b = ix[(a[0], a[1], "C1'")], ix[(b[0], b[1], "C1'")]
    y = c1a - c1b
    y -= (y @ nrm) * nrm
    y /= np.linalg.norm(y)
    x = o - (c1a + c1b) / 2
    x -= (x @ nrm) * nrm
    x -= (x @ y) * y
    x /= np.linalg.norm(x)
    return o, np.array([x, y, np.cross(x, y)])


def Rm(t, r, w):
    t, r, w = np.radians([t, r, w])
    Rx = np.array([[1, 0, 0], [0, np.cos(t), -np.sin(t)], [0, np.sin(t), np.cos(t)]])
    Ry = np.array([[np.cos(r), 0, np.sin(r)], [0, 1, 0], [-np.sin(r), 0, np.cos(r)]])
    Rz = np.array([[np.cos(w), -np.sin(w), 0], [np.sin(w), np.cos(w), 0], [0, 0, 1]])
    return Rz @ Ry @ Rx


def _angles(R):
    roll = float(np.degrees(-np.arcsin(np.clip(R[2, 0], -1, 1))))
    tilt = float(np.degrees(np.arctan2(R[2, 1], R[2, 2])))
    twist = float(np.degrees(np.arctan2(R[1, 0], R[0, 0])))
    return [tilt, roll, twist]


def canonical_step(p):
    """The same step written with |roll| <= 90 deg (Rz(twist) Ry(roll) Rx(tilt) read back from the matrix)."""
    return [*_angles(Rm(*p[:3])), *[float(v) for v in p[3:6]]]


def step_target(o, F, p):
    """Frame reached from (o, F) by the step p = (tilt, roll, twist, shift, slide, rise): origin and rows."""
    Rs = Rm(*p[:3])
    return o + np.asarray(p[3:6]) @ F, Rs.T @ F


def step_between(o1, F1, o2, F2):
    """The step that takes frame (o1, F1) to frame (o2, F2) (the inverse of step_target)."""
    return [*_angles(F1 @ F2.T), *[float(v) for v in (o2 - o1) @ F1.T]]


def move(X, o_cur, F_cur, o_tgt, F_tgt):
    """Rigid motion that takes the frame (o_cur, F_cur) onto (o_tgt, F_tgt), applied to the rows of X."""
    return ((X - o_cur) @ F_cur.T) @ F_tgt + o_tgt


def _perp(ref, axis):
    """Unit vector along ref square to the unit vector axis; when ref runs along the axis, the least aligned
    coordinate axis is used (degenerate geometry only)."""
    u = ref - (ref @ axis) * axis
    n = np.linalg.norm(u)
    if n < 1e-6:
        e = np.eye(3)[int(np.argmin(np.abs(axis)))]
        u = e - (e @ axis) * axis
        n = np.linalg.norm(u)
    return u / n


def _unit(v, fallback):
    n = np.linalg.norm(v)
    return v / n if n > 1e-6 else fallback


def ang3(a, b, c):
    u, v = a - b, c - b
    return float(np.degrees(np.arccos(np.clip(u @ v / max(np.linalg.norm(u) * np.linalg.norm(v), 1e-9), -1, 1))))


def dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v, w = b0 - (b0 @ b1) * b1, b2 - (b2 @ b1) * b1
    return float(np.degrees(np.arctan2(np.cross(b1, v) @ w, v @ w)))


def nonbridging(P, Oa, Ob):
    """OP1/OP2 opposite the two bridging oxygens (OP1-P-OP2 119.6 deg)."""
    u1 = _unit(Oa - P, np.array([1.0, 0, 0]))
    u2 = _unit(Ob - P, np.array([0, 1.0, 0]))
    w = _unit(-(u1 + u2), _unit(np.cross(u1, [0, 0, 1.0]), np.array([0, 0, 1.0])))
    n = _unit(np.cross(u1, u2), _unit(np.cross(w, [1.0, 0, 0]), np.array([0, 1.0, 0])))
    a = np.radians(ANG_OPOP / 2)
    return P + P_OP * (np.cos(a) * w + np.sin(a) * n), P + P_OP * (np.cos(a) * w - np.sin(a) * n)


def rot_about(p, a, b, th):
    """p turned by th (rad) about the axis a -> b."""
    k = (b - a) / np.linalg.norm(b - a)
    v = p - a
    return a + v * np.cos(th) + np.cross(k, v) * np.sin(th) + k * (k @ v) * (1 - np.cos(th))


def sugar_hand(ix, r):
    """'D' or 'L' from the chiral volume at C4' (C3', O4', C5'): positive for D-ribose and D-deoxyribose."""
    try:
        c4 = ix[(r[0], r[1], "C4'")]
        v = (ix[(r[0], r[1], "C3'")] - c4) @ np.cross(ix[(r[0], r[1], "O4'")] - c4, ix[(r[0], r[1], "C5'")] - c4)
    except KeyError:
        return "?"
    return "D" if v > 0 else "L"


FLIP = np.diag([1.0, -1.0, -1.0])                    # the same pair frame with strand I on the other strand
SCREW_SKIP = {"P", "OP1", "OP2", "O5'"}
SUGAR = ("C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'")


def screw_axis(ix, rn, bps, step):
    """Helical axis of a run of base pairs (one residue of each strand per pair, in order along the run), as
    Curve_It's corner mode measures it: the screw that carries every pair onto the pair `step` further on.  Residues
    of the same name share every heavy atom but the phosphate and O5'; different ones only the sugar.  Returns
    (point on the axis, unit direction along the run, rms of the screw fit in A)."""
    names_of = {}
    for (c, r, n) in ix:
        names_of.setdefault((c, r), []).append(n)
    P, Q = [], []
    for bp1, bp2 in zip(bps, bps[step:]):
        for r1, r2 in zip(bp1, bp2):
            names = [n for n in names_of[r1] if n not in SCREW_SKIP]
            if rn[r1] != rn[r2]:
                names = [n for n in names if n in SUGAR]
            for n in names:
                if (r2[0], r2[1], n) in ix:
                    P.append(ix[(r1[0], r1[1], n)])
                    Q.append(ix[(r2[0], r2[1], n)])
    P, Q = np.array(P), np.array(Q)
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = 1.0 if np.linalg.det(Vt.T @ U.T) >= 0 else -1.0
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    t = qc - R @ pc
    rms = float(np.sqrt(np.mean(np.sum((P @ R.T + t - Q) ** 2, 1))))
    w = np.array([R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]])
    u = w / np.linalg.norm(w)
    if u @ t < 0:
        u = -u
    p = np.linalg.lstsq(np.eye(3) - R, t - (u @ t) * u, rcond=None)[0]
    return p - (p @ u) * u, u, rms


def helix_axis(ix, rn, bps):
    """Screw axis with the repeat that fits the run: one pair (B-, A-forms) or two (Z-DNA, C anti / G syn)."""
    p, u, rms = screw_axis(ix, rn, bps, 1)
    if rms > 0.5 and len(bps) > 3:
        p2, u2, rms2 = screw_axis(ix, rn, bps, 2)
        if rms2 < rms:
            return p2, u2, rms2, 2
    return p, u, rms, 1


def axis_relation(c1, d1, c2, d2, at):
    """Bend (deg) between two axis lines (directions made to agree), their sideways offset at the point 'at' and
    their closest approach (A)."""
    if d1 @ d2 < 0:
        d2 = -d2
    bend = float(np.degrees(np.arctan2(np.linalg.norm(np.cross(d1, d2)), d1 @ d2)))
    dm = (d1 + d2) / np.linalg.norm(d1 + d2)
    f1 = c1 + ((at - c1) @ d1) * d1
    f2 = c2 + ((at - c2) @ d2) * d2
    off = f2 - f1
    off -= (off @ dm) * dm
    n = np.cross(d1, d2)
    if np.linalg.norm(n) < 1e-9:
        v = c2 - c1
        skew = float(np.linalg.norm(v - (v @ d1) * d1))
    else:
        skew = float(abs((c2 - c1) @ n) / np.linalg.norm(n))
    return bend, float(np.linalg.norm(off)), skew


# ----------------------------------------------------------------------------- segments and ends
class Segment:
    """A duplex segment as BNP_NA writes it (see the module docstring); residues in file order per chain."""

    def __init__(self, path):
        self.path = Path(path)
        self.atoms = read_atoms(path)
        self.ix = {(a[0], a[1], a[3]): a[4] for a in self.atoms}
        self.rn = {(a[0], a[1]): a[2] for a in self.atoms}
        self.chains = list(dict.fromkeys(a[0] for a in self.atoms))
        if len(self.chains) != 2:
            raise JunctionError(f"{path}: a segment needs exactly two chains (found {self.chains or 'none'})")
        self.res = {c: list(dict.fromkeys((a[0], a[1]) for a in self.atoms if a[0] == c)) for c in self.chains}
        A, B = (self.res[c] for c in self.chains)
        if len(A) != len(B):
            raise JunctionError(f"{path}: the two strands have {len(A)} and {len(B)} residues; give a paired duplex")
        self.n = len(A)
        bad = [r for r in self.rn if self.rn[r] not in BASE]
        if bad:
            raise JunctionError(f"{path}: residue names {sorted({self.rn[r] for r in bad})} are not DNA/RNA bases "
                                f"this script knows ({' '.join(BASE)})")
        for c in self.chains:
            rs = self.res[c]
            for r1, r2 in zip(rs[:3], rs[1:4]):
                fwd = self._d(r1, "O3'", r2, "P")
                back = self._d(r2, "O3'", r1, "P")
                if fwd > 2.2 and back < 2.2:
                    raise JunctionError(f"{path}: chain {c} runs 3'->5' in the file; write each strand 5'->3'")
        for a, b in ((A[0], B[-1]), (A[-1], B[0])):
            d = self._d(a, "C1'", b, "C1'")
            if not 9.0 < d < 12.0:
                raise JunctionError(f"{path}: {a[0]}{a[1]} and {b[0]}{b[1]} should pair at an end (C1'-C1' "
                                    f"{d:.1f} A); strand A's i-th residue must pair with strand B's i-th from the end")
        self.hand = {r: sugar_hand(self.ix, r) for r in self.rn}
        self.names = {}
        for a in self.atoms:
            self.names.setdefault((a[0], a[1]), []).append(a[3])

    def _d(self, r1, n1, r2, n2):
        try:
            return float(np.linalg.norm(self.ix[(*r1, n1)] - self.ix[(*r2, n2)]))
        except KeyError:
            return 99.0

    def end_residues(self, which):
        """(3'-end residue, 5'-end residue) of the end pair: top = strand A's 3' end, bottom = strand B's 3' end."""
        A, B = (self.res[c] for c in self.chains)
        if which == "top":
            return A[-1], B[0]
        if which == "bottom":
            return B[-1], A[0]
        r3, r5 = which
        return tuple(r3), tuple(r5)

    def label(self):
        return self.path.name


class End:
    """One end of a segment: the 3'-end residue r3 and 5'-end residue r5 of its end pair, the residues within
    NWIN of it on both strands (the window a fit uses), the pairs from the end inwards and the chain sense."""

    def __init__(self, seg: Segment, which, nwin=NWIN):
        self.seg, self.which = seg, which
        r3, r5 = seg.end_residues(which)
        s3, s5 = seg.res[r3[0]], seg.res[r5[0]]
        i3, i5 = s3.index(r3), s5.index(r5)
        window = set(s3[max(0, i3 - nwin + 1): i3 + 1]) | set(s5[i5: i5 + nwin])
        at = [a for a in seg.atoms if (a[0], a[1]) in window]
        self.ix, self.rn = seg.ix, seg.rn
        self.r3, self.r5 = r3, r5
        self.key = [(a[0], a[1], a[3]) for a in at]
        self.X = np.array([a[4] for a in at])
        self.k = {kk: i for i, kk in enumerate(self.key)}
        self.pairs = [(s3[i3 - j], s5[i5 + j]) for j in range(min(i3 + 1, len(s5) - i5))]
        # chain sense: with y towards the 3'-end strand, does z point out of the segment (towards the junction)?
        o3, F3 = bp_frame(seg.ix, seg.rn, r3, r5)
        inner = self.pairs[min(3, len(self.pairs) - 1)]
        oi = bp_frame(seg.ix, seg.rn, *inner)[0]
        self.sense = 1 if F3[2] @ (o3 - oi) > 0 else -1
        self.hand = "".join(sorted({seg.hand[r] for r in (r3, r5)}))

    def describe(self):
        w = self.which if isinstance(self.which, str) else f"pair {self.r3[0]}{self.r3[1]}-{self.r5[0]}{self.r5[1]}"
        return (f"{self.seg.label()} {w} (3' end {self.r3[0]}{self.r3[1]} {self.rn[self.r3]}, 5' end "
                f"{self.r5[0]}{self.r5[1]} {self.rn[self.r5]}; {self.hand}-sugars; chain sense {self.sense:+d})")


# ----------------------------------------------------------------------------- the junction
class Junction:
    """The fit of one junction: the step that puts end e2 (of the second segment) next to end e1 (of the first).

    link 'auto' joins the strands as the chain senses allow (module docstring); 'bowtie' or 'continuous' must agree
    with it, since the other way would turn one base pair over."""

    def __init__(self, e1: End, e2: End, link="auto", rise=RISE, gap=GAP):
        natural = "continuous" if e1.sense == e2.sense else "bowtie"
        if link == "auto":
            link = natural
        if link != natural:
            raise JunctionError(
                f"a {link} junction between {e1.describe()} and {e2.describe()} would make one end pair face the "
                f"other way (the chain sense paradox); these ends take a {natural} junction (--link auto)")
        self.e1, self.e2, self.link, self.rise, self.gap = e1, e2, link, rise, gap
        # strand I of the frames: on e1 the strand whose z points into the second segment, on e2 its link partner
        self.y1 = e1.r3 if e1.sense > 0 else e1.r5
        self.x1 = e1.r5 if self.y1 == e1.r3 else e1.r3
        if link == "bowtie":
            self.y2 = e2.r3 if self.y1 == e1.r3 else e2.r5
        else:
            self.y2 = e2.r5 if self.y1 == e1.r3 else e2.r3
        self.x2 = e2.r5 if self.y2 == e2.r3 else e2.r3
        self.o1, self.F1 = bp_frame(e1.ix, e1.rn, self.y1, self.x1)
        self.o2, self.F2 = bp_frame(e2.ix, e2.rn, self.y2, self.x2)
        # atoms out of the contact check: the rebuilt 5' atoms of both end pairs
        self.m1 = np.array([kk[2] not in REBUILT_5 or (kk[0], kk[1]) != e1.r5 for kk in e1.key])
        self.m2 = np.array([kk[2] not in REBUILT_5 or (kk[0], kk[1]) != e2.r5 for kk in e2.key])
        k1, k2 = e1.k, e2.k
        mi1 = lambda i: int(np.sum(self.m1[:i]))                  # index among the kept atoms
        mi2 = lambda i: int(np.sum(self.m2[:i]))
        self.base1 = [k1[(r[0], r[1], n)] for r in (e1.r3, e1.r5) for n in BASE[e1.rn[r]]]
        self.base2 = [k2[(r[0], r[1], n)] for r in (e2.r3, e2.r5) for n in BASE[e2.rn[r]]]
        self.ring1 = [[k1[(r[0], r[1], n)] for n in ring] for r in (e1.r3, e1.r5) for ring in RINGS[e1.rn[r]]]
        self.ring2 = [[k2[(r[0], r[1], n)] for n in ring] for r in (e2.r3, e2.r5) for ring in RINGS[e2.rn[r]]]
        if link == "bowtie":
            self.o3b, self.o3z = k1[(*e1.r3, "O3'")], k2[(*e2.r3, "O3'")]
            self.c3b, self.c3z = k1[(*e1.r3, "C3'")], k2[(*e2.r3, "C3'")]
            self.c5b, self.c5z = k1[(*e1.r5, "C5'")], k2[(*e2.r5, "C5'")]
            self.c4b, self.c4z = k1[(*e1.r5, "C4'")], k2[(*e2.r5, "C4'")]
            self.mb33, self.mz33 = self.m1.copy(), self.m2.copy()
            self.mb33[[self.o3b, self.c3b]] = False
            self.mz33[[self.o3z, self.c3z]] = False
            self.mb55, self.mz55 = self.m1.copy(), self.m2.copy()
            self.mb55[self.c5b] = False
            self.mz55[self.c5z] = False
            self.skip_pairs = [(mi2(self.o3z), mi1(self.o3b))]   # the two O3' bridged by the X33 phosphate
            th = np.radians(np.arange(0, 360, 10))
            self.circ = np.column_stack([np.cos(th), np.sin(th)])
        else:
            # link 0: O3' of e1's 3' end -> P -> O5' of e2's 5' end; link 1: O3' of e2's 3' end -> e1's 5' end
            self.joins = []
            for side3, side5 in ((1, 2), (2, 1)):
                E3, E5 = (e1, e2) if side3 == 1 else (e2, e1)
                m3 = (self.m1 if side3 == 1 else self.m2).copy()
                m5 = (self.m2 if side3 == 1 else self.m1).copy()
                io3, ic3 = E3.k[(*E3.r3, "O3'")], E3.k[(*E3.r3, "C3'")]
                ic4, ic5, io5 = E5.k[(*E5.r5, "C4'")], E5.k[(*E5.r5, "C5'")], E5.k[(*E5.r5, "O5'")]
                m3[[io3, ic3]] = False
                m5[[ic4, ic5]] = False
                self.joins.append(dict(side3=side3, io3=io3, ic3=ic3, ic4=ic4, ic5=ic5, io5=io5,
                                       env1=m3 if side3 == 1 else m5, env2=m5 if side3 == 1 else m3))
            # O3'...C5' across a phosphodiester is 2.6-3.7 A in B- and Z-DNA: not a contact
            j0, j1 = self.joins
            self.skip_pairs = [(mi2(j0["ic5"]), mi1(j0["io3"])), (mi2(j1["io3"]), mi1(j1["ic5"]))]
            self.gam = np.radians(np.arange(0, 360, 15))
            th = np.radians(np.arange(0, 360, 15))
            self.circ = np.column_stack([np.cos(th), np.sin(th)])

    # ------------------------------------------------------------------ placing and scoring
    def place(self, p, X=None):
        o_t, F_t = step_target(self.o1, self.F1, p)
        return move(self.e2.X if X is None else X, self.o2, self.F2, o_t, F_t)

    def _x33(self, XZ):
        """Best X33 position on the circle 1.60 A from both O3' (C3'-O3'-P angles and room), over 36 points; when
        the O3' are too far apart for that, the midpoint (the gap term of the score then dominates)."""
        Oa, Ob, Ca, Cb = self.e1.X[self.o3b], XZ[self.o3z], self.e1.X[self.c3b], XZ[self.c3z]
        mid, dv = (Oa + Ob) / 2, Ob - Oa
        dd = float(np.linalg.norm(dv))
        env = cKDTree(np.vstack([self.e1.X[self.mb33], XZ[self.mz33]]))
        if dd >= 2 * P_O - 1e-3:
            P = mid[None]
        else:
            u = _perp(mid - (Ca + Cb) / 2, dv / dd)          # away from both C3', square to the O3'...O3' line
            v = np.cross(dv / dd, u)
            P = mid + np.sqrt(P_O ** 2 - (dd / 2) ** 2) * (self.circ[:, :1] * u + self.circ[:, 1:] * v)
        unit = lambda V: V / np.maximum(np.linalg.norm(V, axis=-1, keepdims=True), 1e-9)
        u1, u2 = unit(Oa - P), unit(Ob - P)
        w, nn = unit(-(u1 + u2)), unit(np.cross(u1, u2))
        a = np.radians(ANG_OPOP / 2)
        o1 = P + P_OP * (np.cos(a) * w + np.sin(a) * nn)
        o2 = P + P_OP * (np.cos(a) * w - np.sin(a) * nn)
        d, _ = env.query(np.vstack([P, o1, o2]))
        d = d.reshape(3, -1)
        cosang = lambda C, O: np.clip(np.sum(unit(C - O)[None] * unit(P - O), -1), -1, 1)
        cop = np.degrees(np.arccos(np.vstack([cosang(Ca, Oa), cosang(Cb, Ob)])))
        cost = K_ANG * (np.radians(cop - ANG_O3P) ** 2).sum(0) + K_CLASH * (np.maximum(0, 2.9 - d) ** 2).sum(0)
        i = int(np.argmin(cost))
        return dd, (float(cost[i]), P[i], o1[i], o2[i], (float(cop[0, i]), float(cop[1, i])), float(d[:, i].min()))

    def _x35(self, J, X1, X2):
        """Best quick phosphate for one continuous link: the O5' turned about C4'-C5' (24 positions), the P on the
        circle 1.60 A from that O5' and the O3' (24 points), scored on the O...O gap, both C-O-P angles and room
        for P, OP1, OP2 and the O5'."""
        X3, X5 = (X1, X2) if J["side3"] == 1 else (X2, X1)
        Oa, Ca = X3[J["io3"]], X3[J["ic3"]]
        C4b, C5b, O5b = X5[J["ic4"]], X5[J["ic5"]], X5[J["io5"]]
        k = (C5b - C4b) / np.linalg.norm(C5b - C4b)
        v0 = O5b - C4b
        c, s = np.cos(self.gam)[:, None], np.sin(self.gam)[:, None]
        O5 = C4b + v0 * c + np.cross(k, v0) * s + k * (k @ v0) * (1 - c)            # (G, 3)
        dv = O5 - Oa
        dd = np.linalg.norm(dv, axis=1)
        mid = (Oa + O5) / 2
        u = np.array([_perp(r, d) for r, d in zip(mid - (Ca + C5b) / 2, dv / dd[:, None])])   # away from C3', C5'
        v = np.cross(dv / dd[:, None], u)
        r = np.sqrt(np.maximum(P_O ** 2 - (dd / 2) ** 2, 0.0))
        P = mid[:, None, :] + r[:, None, None] * (self.circ[None, :, :1] * u[:, None, :]
                                                    + self.circ[None, :, 1:] * v[:, None, :])          # (G, K, 3)
        unit = lambda V: V / np.maximum(np.linalg.norm(V, axis=-1, keepdims=True), 1e-9)
        u1, u2 = unit(Oa - P), unit(O5[:, None, :] - P)
        w, nn = unit(-(u1 + u2)), unit(np.cross(u1, u2))
        a = np.radians(ANG_OPOP / 2)
        q1 = P + P_OP * (np.cos(a) * w + np.sin(a) * nn)
        q2 = P + P_OP * (np.cos(a) * w - np.sin(a) * nn)
        env = cKDTree(np.vstack([X1[J["env1"]], X2[J["env2"]]]))
        G, K = P.shape[:2]
        d, _ = env.query(np.vstack([P.reshape(-1, 3), q1.reshape(-1, 3), q2.reshape(-1, 3)]))
        d = d.reshape(3, G, K)
        dO5, _ = env.query(O5)
        clash = (np.maximum(0, 2.9 - d) ** 2).sum(0) + (np.maximum(0, 2.9 - dO5) ** 2)[:, None]
        cop3 = np.degrees(np.arccos(np.clip(np.sum(unit(Ca - Oa) * unit(P - Oa), -1), -1, 1)))
        cop5 = np.degrees(np.arccos(np.clip(np.sum(unit(C5b - O5)[:, None, :] * unit(P - O5[:, None, :]), -1), -1, 1)))
        gap = 500 * (np.maximum(0, dd - self.gap[1]) ** 2 + np.maximum(0, self.gap[0] - dd) ** 2)
        cost = K_ANG * (np.radians(cop3 - ANG_O3P) ** 2 + np.radians(cop5 - ANG_O5P) ** 2) + K_CLASH * clash \
            + gap[:, None]
        g, i = np.unravel_index(int(np.argmin(cost)), cost.shape)
        dmin = float(min(d[:, g, i].min(), dO5[g]))
        return float(cost[g, i]), dict(gap=float(dd[g]), cop=(float(cop3[g, i]), float(cop5[g, i])), dmin=dmin,
                                       turn=float(np.degrees(self.gam[g])))

    def metrics(self, p):
        XZ = self.place(p)
        Xb = self.e1.X
        D = np.linalg.norm(XZ[self.m2][:, None, :] - Xb[self.m1][None, :, :], axis=-1)
        for i2, i1 in self.skip_pairs:
            D[i2, i1] = 99.0
        plane = lambda Q: (Q.mean(0), np.linalg.svd(Q - Q.mean(0))[2][2])
        o1, n1 = plane(Xb[self.base1])
        o2, n2 = plane(XZ[self.base2])
        n2 = n2 if n1 @ n2 > 0 else -n2
        nm = (n1 + n2) / np.linalg.norm(n1 + n2)
        e1 = np.cross(nm, self.F1[1])
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(nm, e1)
        om = (o1 + o2) / 2
        poly = lambda Q, rings: unary_union([Polygon([((Q[i] - om) @ e1, (Q[i] - om) @ e2) for i in r]) for r in rings])
        dvec = o2 - o1
        m = dict(dmin=float(D.min()), n26=int((D < 2.6).sum()), clash=float((np.maximum(0, 2.9 - D) ** 2).sum()),
                 ang=float(np.degrees(np.arccos(np.clip(n1 @ n2, -1, 1)))), rise=float(abs(dvec @ nm)),
                 offset=float(np.linalg.norm(dvec - (dvec @ nm) * nm)),
                 ov=float(poly(Xb, self.ring1).intersection(poly(XZ, self.ring2)).area))
        if self.link == "bowtie":
            d33, (c33, P33, q1, q2, cop, x33_dmin) = self._x33(XZ)
            C5a, C5b, C4a, C4b = Xb[self.c5b], XZ[self.c5z], Xb[self.c4b], XZ[self.c4z]
            m55 = (C5a + C5b) / 2
            e5 = m55 - (C4a + C4b) / 2
            c5v = C5b - C5a
            e5 = e5 - (e5 @ c5v) / max(c5v @ c5v, 1e-9) * c5v
            P55 = m55 + 1.0 * e5 / max(np.linalg.norm(e5), 1e-9)
            env55 = cKDTree(np.vstack([Xb[self.mb55], XZ[self.mz55], [P33, q1, q2]]))
            d55c, _ = env55.query(P55)
            m.update(x33_cost=c33, x33_dmin=x33_dmin, x33_cop=list(cop), p55_dmin=float(d55c),
                     p55_clash=float(max(0, 3.0 - d55c) ** 2), d33=d33, d55=float(np.linalg.norm(C5a - C5b)))
        else:
            est = [self._x35(J, Xb, XZ) for J in self.joins]
            m.update(link_cost=float(sum(c for c, _ in est)), links=[e for _, e in est])
        return m

    def score(self, p):
        m = self.metrics(p)
        stack = 0.1 * m["ang"] ** 2 + 10 * max(0, abs(m["rise"] - self.rise) - 0.15) ** 2 - 0.5 * m["ov"]
        if self.link == "bowtie":
            return (20 * m["clash"] + m["x33_cost"] + 20 * m["p55_clash"]
                    + 500 * max(0, m["d33"] - self.gap[1]) ** 2 + 500 * max(0, self.gap[0] - m["d33"]) ** 2
                    + 5 * max(0, m["d55"] - 4.8) ** 2 + 5 * max(0, 3.5 - m["d55"]) ** 2
                    + 0.1 * m["ang"] ** 2 + 10 * max(0, abs(m["rise"] - self.rise) - 0.15) ** 2 - 0.5 * m["ov"])
        return 20 * m["clash"] + m["link_cost"] + stack

    def stack_score(self, m):
        """The fit score without its phosphate estimates: contacts and stacking only."""
        return (20 * m["clash"] + 0.1 * m["ang"] ** 2 + 10 * max(0, abs(m["rise"] - self.rise) - 0.15) ** 2
                - 0.5 * m["ov"])

    def axes(self, p):
        """Helical axes of the two windows with the second end placed by p: bend (deg), sideways offset at the
        junction and closest approach (A), the repeat used on each side (1 or 2 pairs)."""
        XZ = self.place(p)
        ix1 = dict(zip(self.e1.key, self.e1.X))
        ix2 = dict(zip(self.e2.key, XZ))
        inside = lambda ix, pr: all((*r, "C1'") in ix for r in pr)
        run1 = list(reversed([pr for pr in self.e1.pairs if inside(ix1, pr)]))
        run2 = [pr for pr in self.e2.pairs if inside(ix2, pr)]
        p1, u1, _, s1 = helix_axis(ix1, self.e1.rn, run1)
        p2, u2, _, s2 = helix_axis(ix2, self.e2.rn, run2)
        at = (bp_frame(ix1, self.e1.rn, self.e1.r3, self.e1.r5)[0] + bp_frame(ix2, self.e2.rn, self.e2.r3, self.e2.r5)[0]) / 2
        bend, off, skew = axis_relation(p1, u1, p2, u2, at)
        return dict(bend=bend, sideways=off, closest=skew, repeat=[s1, s2])

    # ------------------------------------------------------------------ real phosphates for a candidate step
    def bridges(self, p):
        """Both links built for the step p on the two windows (the second end placed)."""
        XZ = self.place(p)
        B, Z = self.e1, self.e2
        reb_b = {(*B.r5, n) for n in REBUILT_5}
        reb_z = {(*Z.r5, n) for n in REBUILT_5}

        def env(skip_b, skip_z):
            return np.array([B.X[i] for i, kk in enumerate(B.key) if kk not in reb_b | skip_b] +
                            [XZ[i] for i, kk in enumerate(Z.key) if kk not in reb_z | skip_z])
        g = lambda E, X, r, n: X[E.k[(*r, n)]]
        if self.link == "bowtie":
            br = build_bowtie(env({(*B.r3, "O3'"), (*B.r3, "C3'")}, {(*Z.r3, "O3'"), (*Z.r3, "C3'")}),
                              env({(*B.r5, "C5'"), (*B.r5, "C4'")}, {(*Z.r5, "C5'"), (*Z.r5, "C4'")}),
                              g(B, B.X, B.r3, "O3'"), g(B, B.X, B.r3, "C3'"), g(Z, XZ, Z.r3, "O3'"), g(Z, XZ, Z.r3, "C3'"),
                              g(B, B.X, B.r5, "C4'"), g(B, B.X, B.r5, "C5'"), g(B, B.X, B.r5, "O5'"),
                              g(Z, XZ, Z.r5, "C4'"), g(Z, XZ, Z.r5, "C5'"), g(Z, XZ, Z.r5, "O5'"))
            br["geometry"] = [br["g33"], br["g55"]]
            return br
        out = []
        for (E3, X3), (E5, X5) in (((B, B.X), (Z, XZ)), ((Z, XZ), (B, B.X))):
            skip3 = {(*E3.r3, "O3'"), (*E3.r3, "C3'")}
            skip5 = {(*E5.r5, "C5'"), (*E5.r5, "C4'")}
            pts = env(skip3 if E3 is B else skip5, skip5 if E3 is B else skip3)
            tree = cKDTree(pts)
            O5, P, ops, geo = bridge_35_ff(g(E3, X3, E3.r3, "O3'"), g(E3, X3, E3.r3, "C3'"),
                                           g(E5, X5, E5.r5, "C4'"), g(E5, X5, E5.r5, "C5'"), g(E5, X5, E5.r5, "O5'"),
                                           tree)
            out.append(dict(O5=O5, P=P, OP1=ops[0], OP2=ops[1], g=geo))
        return dict(links=out, geometry=[o["g"] for o in out])


# ----------------------------------------------------------------------------- phosphate builders
def _ff_terms(P, Oa, Ob, Ca, Cb, ang_a, ang_b, tree, extra=()):
    o1, o2 = nonbridging(P, Oa, Ob)
    pts = np.array([P, o1, o2, *extra])
    if not np.all(np.isfinite(pts)):
        return 1e9
    d, _ = tree.query(pts)
    rad = np.radians
    return (K_BOND * ((np.linalg.norm(P - Oa) - P_O) ** 2 + (np.linalg.norm(P - Ob) - P_O) ** 2)
            + K_ANG * (rad(ang3(Ca, Oa, P) - ang_a) ** 2 + rad(ang3(Cb, Ob, P) - ang_b) ** 2
                       + rad(ang3(Oa, P, Ob) - ANG_OPO_BR) ** 2)
            + K_CLASH * float((np.maximum(0, 2.9 - d) ** 2).sum()))


def _local_frame(o, a, b):
    """Rows of an orthonormal frame attached to atoms: x along a, y along b square to it, z = x cross y."""
    x = a / np.linalg.norm(a)
    y = _perp(b, x)
    return np.array([x, y, np.cross(x, y)])


def bridge_33_ff(O3a, C3a, O3b, C3b, tree):
    """The X33 phosphate between two O3' (both fixed).  Minimized in a frame attached to the two O3' and C3' (so where
    the molecule sits in space does not matter) from 24 starts around the O3'...O3' midpoint."""
    mid = (O3a + O3b) / 2
    E = _local_frame(mid, O3b - O3a, mid - (C3a + C3b) / 2)
    f = lambda l: _ff_terms(mid + l @ E, O3a, O3b, C3a, C3b, ANG_O3P, ANG_O3P, tree)
    best = None
    dirs = [np.array([0.0, np.cos(t), np.sin(t)]) for t in np.radians(np.arange(0, 360, 45))]
    for s_ in (0.5, 1.0, 1.5):
        for d in dirs:
            r = minimize(f, s_ * d, method="Nelder-Mead", options=dict(xatol=1e-5, fatol=1e-9, maxiter=6000))
            best = r if best is None or r.fun < best.fun else best
    P = mid + best.x @ E
    return P, nonbridging(P, O3a, O3b), dict(bonds=(float(np.linalg.norm(P - O3a)), float(np.linalg.norm(P - O3b))),
                                             OPO=ang3(O3a, P, O3b), COP=(ang3(C3a, O3a, P), ang3(C3b, O3b, P)),
                                             energy=float(best.fun))


def bridge_55_ff(C4a, C5a, O5a, C4b, C5b, O5b, tree):
    """One phosphate between two O5', each turned about its C4'-C5' bond to meet it.  The phosphate is placed in a frame
    attached to the two C5' and C4' atoms; 36 starts (both O5' turned in 60-deg steps)."""
    mid = (C5a + C5b) / 2
    E = _local_frame(mid, C5b - C5a, mid - (C4a + C4b) / 2)

    def geom(v):
        return rot_about(O5a, C4a, C5a, v[0]), rot_about(O5b, C4b, C5b, v[1]), mid + np.asarray(v[2:5]) @ E

    def f(v):
        Oa, Ob, P = geom(v)
        return _ff_terms(P, Oa, Ob, C5a, C5b, ANG_O5P, ANG_O5P, tree, extra=(Oa, Ob))
    best = None
    for ga in np.radians([-120, -60, 0, 60, 120, 180]):
        for gb in np.radians([-120, -60, 0, 60, 120, 180]):
            p0 = (rot_about(O5a, C4a, C5a, ga) + rot_about(O5b, C4b, C5b, gb)) / 2
            r = minimize(f, np.r_[ga, gb, (p0 - mid) @ E.T], method="Powell",
                         options=dict(xtol=1e-5, ftol=1e-9, maxiter=8000))
            best = r if best is None or r.fun < best.fun else best
    Oa, Ob, P = geom(best.x)
    return Oa, Ob, P, nonbridging(P, Oa, Ob), dict(bonds=(float(np.linalg.norm(P - Oa)), float(np.linalg.norm(P - Ob))),
                                                   OPO=ang3(Oa, P, Ob), COP=(ang3(C5a, Oa, P), ang3(C5b, Ob, P)),
                                                   C5C5=float(np.linalg.norm(C5a - C5b)), energy=float(best.fun))


def build_bowtie(env33, env55, O3a, C3a, O3b, C3b, C4a, C5a, O5a, C4b, C5b, O5b):
    """X33 between the two O3', then one phosphate between the two O5' (which must avoid the X33 too)."""
    tree = cKDTree(env33)
    P33, (o1, o2), g33 = bridge_33_ff(O3a, C3a, O3b, C3b, tree)
    tree55 = cKDTree(np.vstack([env55, [P33, o1, o2]]))
    Oa, Ob, P55, (q1, q2), g55 = bridge_55_ff(C4a, C5a, O5a, C4b, C5b, O5b, tree55)
    return dict(P33=P33, OP1_33=o1, OP2_33=o2, g33=g33, O5a=Oa, O5b=Ob, P55=P55, OP1_55=q1, OP2_55=q2, g55=g55)


def bridge_35_ff(O3a, C3a, C4b, C5b, O5b, tree):
    """One ordinary phosphate between the O3' of residue a (fixed) and the O5' of residue b, which turns about its
    C4'-C5' bond (the gamma torsion) to meet it; residue b keeps the phosphate.  The phosphate is placed in a frame
    attached to O3', the O5' and the C3'/C5' atoms (so where the molecule sits in space does not matter); 36 starts
    (the O5' turned in 30-deg steps x three phosphate positions)."""
    E = _local_frame(O3a, O5b - O3a, (O3a + O5b) / 2 - (C3a + C5b) / 2)

    def geom(v):
        return rot_about(O5b, C4b, C5b, v[0]), O3a + np.asarray(v[1:4]) @ E

    def f(v):
        Ob, P = geom(v)
        return _ff_terms(P, O3a, Ob, C3a, C5b, ANG_O3P, ANG_O5P, tree, extra=(Ob,))
    best = None
    for gb in np.radians(np.arange(-150, 181, 30)):
        Ob = rot_about(O5b, C4b, C5b, gb)
        mid = (O3a + Ob) / 2
        away = mid - (C3a + C5b) / 2
        away = away / max(np.linalg.norm(away), 1e-9)
        for s_ in (0.0, 0.8, -0.8):
            r = minimize(f, np.r_[gb, (mid + s_ * away - O3a) @ E.T], method="Powell",
                         options=dict(xtol=1e-5, ftol=1e-9, maxiter=8000))
            best = r if best is None or r.fun < best.fun else best
    Ob, P = geom(best.x)
    turn = float((np.degrees(best.x[0]) + 180) % 360 - 180)
    return Ob, P, nonbridging(P, O3a, Ob), dict(bonds=(float(np.linalg.norm(P - O3a)), float(np.linalg.norm(P - Ob))),
                                                OPO=ang3(O3a, P, Ob), COP=(ang3(C3a, O3a, P), ang3(C5b, Ob, P)),
                                                O3O5=float(np.linalg.norm(O3a - Ob)), O5_turn=turn,
                                                energy=float(best.fun))


# ----------------------------------------------------------------------------- fitting (parallel)
_SEGS: dict = {}
_JX: list = []


def get_segment(path) -> Segment:
    key = str(Path(path).resolve())
    if key not in _SEGS:
        _SEGS[key] = Segment(key)
    return _SEGS[key]


def make_junction(spec) -> Junction:
    """spec: dict(seg1, end1, seg2, end2[, link, nwin, rise, gap]); ends 'top', 'bottom' or [[c, r3], [c, r5]]."""
    e1 = End(get_segment(spec["seg1"]), _end(spec["end1"]), spec.get("nwin", NWIN))
    e2 = End(get_segment(spec["seg2"]), _end(spec["end2"]), spec.get("nwin", NWIN))
    return Junction(e1, e2, spec.get("link", "auto"), spec.get("rise", RISE), tuple(spec.get("gap", GAP)))


def _end(e):
    return e if isinstance(e, str) else tuple(tuple(x) for x in e)


def _init(specs):
    global _JX
    _JX = [make_junction(s) for s in specs]


_INIT_ERROR = None


def _init_worker(specs):
    """Pool initializer: _init, but a failure is kept and raised by the worker's first task as an ordinary error.  A
    worker that dies in its initializer is replaced by one that dies the same way, and the pool would wait forever."""
    global _INIT_ERROR
    try:
        _init(specs)
        _INIT_ERROR = None
    except BaseException as exc:  # noqa: BLE001 - JunctionError is a SystemExit, which would kill the worker
        _INIT_ERROR = f"{type(exc).__name__}: {exc}"


def _check_worker():
    if _INIT_ERROR is not None:
        raise RuntimeError(f"a fitting process could not set up the junctions ({_INIT_ERROR})")


def _search(task):
    _check_worker()
    j, p0 = task
    r = minimize(_JX[j].score, np.array(p0, float), method="Powell", options=dict(maxiter=4000, xtol=1e-3, ftol=1e-4))
    return j, float(r.fun), [float(v) for v in r.x]


def _real(task):
    """Real phosphates for one candidate: stacking/contact score + the force-field energies of both phosphates."""
    _check_worker()
    j, p = task
    J = _JX[j]
    m = J.metrics(p)
    br = J.bridges(p)
    total = J.stack_score(m) + sum(g["energy"] for g in br["geometry"])
    return j, p, m, br["geometry"], float(total)


def distinct(runs, n):
    """The n best runs whose steps differ (canonical angles > 3 deg or translations > 0.5 A apart)."""
    keep = []
    for f, p in runs:
        c = np.array(canonical_step(p))
        if all(np.abs(((c[:3] - q[:3] + 180) % 360) - 180).max() > 3 or np.abs(c[3:] - q[3:]).max() > 0.5
               for _, _, q in keep):
            keep.append((f, p, c))
        if len(keep) == n:
            break
    return keep


def reverse_spec(spec):
    """The same junction read the other way: the second segment's end first."""
    r = dict(spec)
    r["seg1"], r["seg2"], r["end1"], r["end2"] = spec["seg2"], spec["seg1"], spec["end2"], spec["end1"]
    return r


def role(J, which):
    """'r3' or 'r5': which residue of its end pair the y axis of J's frame on end `which` (1 or 2) points to."""
    e, y = (J.e1, J.y1) if which == 1 else (J.e2, J.y2)
    return "r3" if y == e.r3 else "r5"


def forward_step(Jf, Jr, p_r):
    """The step of junction Jf that gives the geometry of step p_r of Jr, where Jr joins the same two end pairs the
    other way round (Jr.e1 is, or matches, Jf.e2 and Jr.e2 matches Jf.e1)."""
    o, F = step_target(Jr.o1, Jr.F1, p_r)                      # Jr's second end (Jf's first end), placed
    F1 = F if role(Jf, 1) == role(Jr, 2) else FLIP @ F          # the same pair; y on Jf's strand I
    F2 = Jr.F1 if role(Jf, 2) == role(Jr, 1) else FLIP @ Jr.F1
    return step_between(o, F1, Jr.o1, F2)


def backward_step(Jf, Jr, p_f):
    """The inverse of forward_step: the step of Jr that gives the geometry of step p_f of Jf."""
    o, F = step_target(Jf.o1, Jf.F1, p_f)                      # Jf's second end (Jr's first end), placed
    F1 = F if role(Jr, 1) == role(Jf, 2) else FLIP @ F
    F2 = Jf.F1 if role(Jr, 2) == role(Jf, 1) else FLIP @ Jf.F1
    return step_between(o, F1, Jf.o1, F2)


def end_match(e, f, tol=SAME_END_RMS):
    """rms (A) when ends e and f present the same end pair (same chain sense, sugars and residue names at the end
    pair, and their windows superpose within tol: residues matched by their place from the end, every atom where the
    residue names agree, backbone and sugar atoms where they do not), else None."""
    if (e.sense, e.hand) != (f.sense, f.hand) or e.rn[e.r3] != f.rn[f.r3] or e.rn[e.r5] != f.rn[f.r5]:
        return None
    P, Q = [], []
    for pe, pf in zip(e.pairs[:NWIN], f.pairs[:NWIN]):
        for ra, rb in zip(pe, pf):
            same = e.rn[ra] == f.rn[rb]
            for nm in e.seg.names[ra]:
                if (rb[0], rb[1], nm) in f.ix and (same or nm in BACKBONE):
                    P.append(e.ix[(ra[0], ra[1], nm)])
                    Q.append(f.ix[(rb[0], rb[1], nm)])
    if len(P) < 3:
        return None
    rms = kabsch_rms(np.array(P), np.array(Q))
    return rms if rms <= tol else None


def kabsch_rms(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = 1.0 if np.linalg.det(Vt.T @ U.T) >= 0 else -1.0
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return float(np.sqrt(np.mean(np.sum(((P - pc) @ R.T - (Q - qc)) ** 2, 1))))


def same_junction(Ja, Jb, tol=SAME_END_RMS):
    """'aligned' when Jb joins the same two end pairs as Ja in the same order, 'reversed' in the other order (a turn
    of 180 deg about the duplex dyad), None when the ends differ or the link kinds do."""
    if Ja.link != Jb.link:
        return None
    if end_match(Ja.e1, Jb.e1, tol) is not None and end_match(Ja.e2, Jb.e2, tol) is not None:
        return "aligned"
    if end_match(Ja.e1, Jb.e2, tol) is not None and end_match(Ja.e2, Jb.e1, tol) is not None:
        return "reversed"
    return None


def fit_junctions(specs, cores=None, twists=TWISTS, shifts=SHIFTS, n_cand=N_CAND, log=print):
    """Fit every junction spec on its own (no sharing between specs); see fit_model_junctions."""
    return fit_model_junctions(specs, cores=cores, twists=twists, shifts=shifts, n_cand=n_cand, share=False, log=log)


def fit_model_junctions(specs, weights=None, cores=None, twists=TWISTS, shifts=SHIFTS, n_cand=N_CAND, share=True,
                        tol=SAME_END_RMS, log=print):
    """Fit the junctions of one model.  specs: distinct junction specs; weights: how often each occurs (default 1).

    Junctions that join the same two end pairs (same_junction: in the same order, or turned 180 deg) share one fit
    when share is True: every search of every such junction, each read from both ends, is pooled in the frames of the
    first of them, the n_cand best distinct steps get real phosphates at every junction of the group, and the step with
    the lowest summed total (weighted by occurrence) is kept for all of them.  A junction alone in its group is fitted
    the same way, from both ends.  Returns one dict per spec: its step, total, metrics, phosphates, candidates and runs
    (all in that junction's own frames), and the group it was fitted with."""
    cores = cores or max(1, (os.cpu_count() or 2) - 1)
    weights = list(weights or [1] * len(specs))
    own = [make_junction(s) for s in specs]
    groups = []                                    # [ref index, [(member index, 'aligned'|'reversed', rms text)]]
    for i, J in enumerate(own):
        for g in groups:
            how = same_junction(own[g[0]], J, tol) if share else None
            if how:
                g[1].append((i, how))
                break
        else:
            groups.append([i, [(i, "aligned")]])
    # every reading: each member forward and reversed, with how it lies relative to its group's reference
    readings = []                                  # (group, member, spec, how relative to the reference)
    for gi, (ref, members) in enumerate(groups):
        for i, how in members:
            readings.append((gi, i, specs[i], how))
            readings.append((gi, i, reverse_spec(specs[i]), "reversed" if how == "aligned" else "aligned"))
    allspecs = [r[2] for r in readings] + list(specs)              # readings, then the members themselves
    base = len(readings)
    for s in allspecs[:base]:                                      # a reading that cannot be set up stops the run
        make_junction(s)                                           # here, before any worker starts
    tasks = [(k, [0, 0, w, sx, sy, 3.6]) for k in range(len(readings)) for w in twists for sx in shifts for sy in shifts]
    t0 = time.time()
    log(f"{len(tasks)} searches ({len(readings)} readings x {len(tasks) // max(1, len(readings))} starts) "
        f"on {cores} process(es)")
    with mp.Pool(cores, initializer=_init_worker, initargs=(allspecs,)) as pool:
        res, every = [], max(1, len(tasks) // 10)
        for r in pool.imap(_search, tasks, chunksize=2):           # the same tasks, results in the same order as map
            res.append(r)
            if len(res) % every == 0 or len(res) == len(tasks):
                log(f"  searches done: {len(res)}/{len(tasks)} ({time.time() - t0:.0f} s)")
        _init(allspecs)
        pooled = {gi: [] for gi in range(len(groups))}             # group -> [score, ref step, start, reading]
        for (k, p0), (kk, f, q) in zip(tasks, res):
            gi, i, _, how = readings[k]
            Jref, R = _JX[base + groups[gi][0]], _JX[k]
            qref = list(q) if how == "aligned" else forward_step(Jref, R, q)
            pooled[gi].append([f, qref, p0, k])
        # second pass: perturbed restarts around the N_REFINE best distinct optima of each group (reference frames)
        kicks = []
        for gi, (ref, members) in enumerate(groups):
            top = distinct(sorted([(f, q) for f, q, _, _ in pooled[gi]], key=lambda t: t[0]), N_REFINE)
            for _, q, _ in top:
                for d in REFINE_KICKS:
                    kicks.append((gi, [a + b for a, b in zip(q, d)]))
        if kicks:
            ref_index = {gi: base + groups[gi][0] for gi in range(len(groups))}
            log(f"  second pass: {len(kicks)} searches around the best distinct optima ({len(REFINE_KICKS)} nudges each)")
            again = pool.map(_search, [(ref_index[gi], p0) for gi, p0 in kicks], chunksize=2)
            for (gi, p0), (_, f, q) in zip(kicks, again):
                pooled[gi].append([f, list(q), p0, -1])
        cand = {gi: distinct(sorted([(f, q) for f, q, _, _ in pooled[gi]], key=lambda t: t[0]), n_cand)
                for gi in pooled}
        jobs = []                                                  # (group, candidate, member, member's step)
        for gi, (ref, members) in enumerate(groups):
            Jref = _JX[base + ref]
            for ci, c in enumerate(cand[gi]):
                for i, how in members:
                    p = list(c[2]) if how == "aligned" else backward_step(Jref, _JX[base + i], list(c[2]))
                    jobs.append((gi, ci, i, p))
        log(f"  real phosphates for {len(jobs)} candidate placement(s)")
        built = pool.map(_real, [(base + i, p) for gi, ci, i, p in jobs], chunksize=1)
    _init(specs)
    real = {}                                                      # (group, candidate, member) -> (p, m, g, total)
    for (gi, ci, i, p), (_, pp, m, g, total) in zip(jobs, built):
        real[(gi, ci, i)] = (pp, m, g, total)
    out = [None] * len(specs)
    for gi, (ref, members) in enumerate(groups):
        sums = []
        for ci in range(len(cand[gi])):
            sums.append(sum(weights[i] * real[(gi, ci, i)][3] for i, _ in members))
        order = sorted(range(len(sums)), key=lambda ci: sums[ci])
        best = order[0]
        Jref = _JX[ref]
        for i, how in members:
            J = _JX[i]
            p, m, g, total = real[(gi, best, i)]
            runs = []
            for f, qref, p0, k in pooled[gi]:
                q = qref if how == "aligned" else backward_step(Jref, J, qref)
                how_read = ("second pass" if k < 0 else f"junction {readings[k][1] + 1} "
                            + ("forward" if readings[k][2] is specs[readings[k][1]] else "reversed"))
                runs.append(dict(start=[float(p0[2]), float(p0[3]), float(p0[4])], start_step=[float(v) for v in p0],
                                 score=f, step=list(q), reading=how_read))
            cands = []
            for ci in order:
                pp, mm, gg, t = real[(gi, ci, i)]
                cands.append(dict(step=list(pp), total=t, stack=J.stack_score(mm), energies=[x["energy"] for x in gg],
                                  metrics=mm, bridges=gg, group_total=sums[ci]))
            out[i] = dict(spec=specs[i], link=J.link, sense=[J.e1.sense, J.e2.sense], step=list(p), total=total,
                          stack=J.stack_score(m), metrics=m, bridges=g, candidates=cands, runs=runs,
                          group=dict(members=[mi for mi, _ in members], reference=ref, how=how,
                                     end_rms=[end_match(J.e1, Jref.e1 if how == "aligned" else Jref.e2, 9),
                                              end_match(J.e2, Jref.e2 if how == "aligned" else Jref.e1, 9)],
                                     group_total=sums[best], searches=len(pooled[gi])))
    log(f"fitted {len(specs)} junction(s) in {len(groups)} group(s) in {time.time() - t0:.0f} s on {cores} cores")
    return out


# ----------------------------------------------------------------------------- the model
class Piece:
    """A rigid copy of a segment in the model (upright or flipped)."""

    def __init__(self, seg: Segment, flip: bool):
        self.seg, self.flip = seg, flip
        self.key = [(a[0], a[1], a[3]) for a in seg.atoms]
        self.X = np.array([a[4] for a in seg.atoms])
        self.k = {kk: i for i, kk in enumerate(self.key)}
        self.entry, self.exit = ("top", "bottom") if flip else ("bottom", "top")

    def ix(self):
        return {kk: self.X[i] for kk, i in self.k.items()}

    def xyz(self, r, name):
        return self.X[self.k[(r[0], r[1], name)]]


def place_pieces(pieces, juncs, fits):
    """Stack the pieces: each next one placed on the one before with its junction's fitted step."""
    for j, (J, fit) in enumerate(zip(juncs, fits)):
        lo, hi = pieces[j], pieces[j + 1]
        o1, F1 = bp_frame(lo.ix(), lo.seg.rn, J.y1, J.x1)
        o2, F2 = bp_frame(hi.ix(), hi.seg.rn, J.y2, J.x2)
        hi.X = move(hi.X, o2, F2, *step_target(o1, F1, fit["step"]))


def build_links(pieces, juncs):
    """The junction phosphates on the placed pieces (environment: both pieces within 16 A of each link)."""
    out = []
    for j, J in enumerate(juncs):
        lo, hi = pieces[j], pieces[j + 1]
        rebuilt = {(0, J.e1.r5, n) for n in REBUILT_5} | {(1, J.e2.r5, n) for n in REBUILT_5}

        def env(skip, centre):
            return np.array([f.X[i] for fi, f in enumerate((lo, hi)) for kk, i in f.k.items()
                             if (fi, (kk[0], kk[1]), kk[2]) not in rebuilt | skip
                             and np.linalg.norm(f.X[i] - centre) < 16])
        if J.link == "bowtie":
            a3, b3, a5, b5 = J.e1.r3, J.e2.r3, J.e1.r5, J.e2.r5
            centre = (lo.xyz(a3, "O3'") + hi.xyz(b3, "O3'")) / 2
            skip33 = {(0, a3, "O3'"), (0, a3, "C3'"), (1, b3, "O3'"), (1, b3, "C3'")}
            skip55 = {(0, a5, "C5'"), (1, b5, "C5'"), (0, a5, "C4'"), (1, b5, "C4'")}
            br = build_bowtie(env(skip33, centre), env(skip55, centre),
                              lo.xyz(a3, "O3'"), lo.xyz(a3, "C3'"), hi.xyz(b3, "O3'"), hi.xyz(b3, "C3'"),
                              lo.xyz(a5, "C4'"), lo.xyz(a5, "C5'"), lo.xyz(a5, "O5'"),
                              hi.xyz(b5, "C4'"), hi.xyz(b5, "C5'"), hi.xyz(b5, "O5'"))
            br["geometry"] = [br["g33"], br["g55"]]
            out.append(br)
            continue
        links = []
        for (i3, P3, r3), (i5, P5, r5) in (((0, lo, J.e1.r3), (1, hi, J.e2.r5)), ((1, hi, J.e2.r3), (0, lo, J.e1.r5))):
            skip = {(i3, r3, "O3'"), (i3, r3, "C3'"), (i5, r5, "C5'"), (i5, r5, "C4'")}
            pts = env(skip, P3.xyz(r3, "O3'"))
            tree = cKDTree(pts)
            O5, P, ops, geo = bridge_35_ff(P3.xyz(r3, "O3'"), P3.xyz(r3, "C3'"),
                                           P5.xyz(r5, "C4'"), P5.xyz(r5, "C5'"), P5.xyz(r5, "O5'"), tree)
            links.append(dict(O5=O5, P=P, OP1=ops[0], OP2=ops[1], g=geo, three=(i3, r3), five=(i5, r5)))
        out.append(dict(links=links, geometry=[l["g"] for l in links]))
    return out


def trace_chains(pieces, juncs):
    """Residue order of the output chains.  Nodes: (piece, chain, resi) and ('X33', j); edges: ('35', u3, v5) an
    O3'(u3)-P(v5) bond, ('X', x, r) an X33 P-O3'(r) bond, ('55', keeper, other) a 5'-5' phosphate kept by 'keeper'."""
    edges = {}
    for pi, f in enumerate(pieces):
        for c in f.seg.chains:
            rs = f.seg.res[c]
            for r1, r2 in zip(rs, rs[1:]):
                edges[frozenset({(pi, *r1), (pi, *r2)})] = ("35", (pi, *r1), (pi, *r2))
    for j, J in enumerate(juncs):
        lo, hi = j, j + 1
        if J.link == "bowtie":
            x = ("X33", j + 1)
            a3, b3 = (lo, *J.e1.r3), (hi, *J.e2.r3)
            edges[frozenset({a3, x})] = ("X", x, a3)
            edges[frozenset({b3, x})] = ("X", x, b3)
            a5, b5 = (lo, *J.e1.r5), (hi, *J.e2.r5)
            edges[frozenset({a5, b5})] = ("55", a5, b5)                    # the lower piece keeps the phosphate
        else:
            edges[frozenset({(lo, *J.e1.r3), (hi, *J.e2.r5)})] = ("35", (lo, *J.e1.r3), (hi, *J.e2.r5))
            edges[frozenset({(hi, *J.e2.r3), (lo, *J.e1.r5)})] = ("35", (hi, *J.e2.r3), (lo, *J.e1.r5))
    adj = {}
    for e in edges:
        u, v = tuple(e)
        adj.setdefault(u, []).append(v)
        adj.setdefault(v, []).append(u)
    for pi, f in enumerate(pieces):                                       # lone residues (none in a duplex)
        for r in f.seg.rn:
            adj.setdefault((pi, *r), [])
    ends = [n for n, v in adj.items() if len(v) <= 1]

    def five_prime(n):
        return n[0] != "X33" and not any(edges[frozenset({n, m})][0] == "35" and edges[frozenset({n, m})][2] == n
                                         for m in adj[n])
    first = (0, *pieces[0].seg.res[pieces[0].seg.chains[0]][0])
    starts = ([first] if first in ends else []) + sorted(
        [e for e in ends if e != first], key=lambda e: (not five_prime(e), -e[0] if e[0] != "X33" else 0))
    seen, paths = set(), []
    for s in starts:
        if s in seen:
            continue
        path, prev, cur = [], None, s
        while cur is not None:
            path.append(cur)
            seen.add(cur)
            nxt = [v for v in adj[cur] if v != prev and v not in seen]
            prev, cur = cur, (nxt[0] if nxt else None)
        paths.append(path)
    return paths, edges


def write_model(pieces, juncs, fits, links, out: Path, remarks, extra_head=()):
    paths, edges = trace_chains(pieces, juncs)
    coords = [dict(zip(f.key, f.X)) for f in pieces]
    xatoms = {}
    for j, (J, br) in enumerate(zip(juncs, links)):
        lo, hi = coords[j], coords[j + 1]
        if J.link == "bowtie":
            lo[(*J.e1.r5, "O5'")], hi[(*J.e2.r5, "O5'")] = br["O5a"], br["O5b"]
            lo[(*J.e1.r5, "P")], lo[(*J.e1.r5, "OP1")], lo[(*J.e1.r5, "OP2")] = br["P55"], br["OP1_55"], br["OP2_55"]
            for n in ("P", "OP1", "OP2"):
                hi.pop((*J.e2.r5, n), None)
            xatoms[("X33", j + 1)] = [("P", br["P33"], "P"), ("OP1", br["OP1_33"], "O"), ("OP2", br["OP2_33"], "O")]
        else:
            for l in br["links"]:
                c = coords[j + l["five"][0]]
                r5 = l["five"][1]
                c[(*r5, "O5'")], c[(*r5, "P")], c[(*r5, "OP1")], c[(*r5, "OP2")] = l["O5"], l["P"], l["OP1"], l["OP2"]
    chain_ids = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    num = {}
    for ch, path in zip(chain_ids, paths):
        for i, item in enumerate(path, 1):
            num[item] = (ch, i)
    lines, serial = [], 1
    hands = {}
    for ch, path in zip(chain_ids, paths):
        for item in path:
            och, ores = num[item]
            if item[0] == "X33":
                atoms, rec, rn = xatoms[item], "HETATM", "X33"
            else:
                pi, c, r = item
                seg = pieces[pi].seg
                rn = seg.rn[(c, r)]
                hands[(och, ores)] = seg.hand[(c, r)]
                names = [n for n in ORDER if (c, r, n) in coords[pi]] + \
                        [kk[2] for kk in pieces[pi].key if (kk[0], kk[1]) == (c, r) and kk[2] not in ORDER
                         and kk in coords[pi]]
                atoms = [(n, coords[pi][(c, r, n)], n[0]) for n in names]
                rec = "ATOM  "
            for n, xyz, el in atoms:
                nm = f" {n:<3s}" if len(n) < 4 else n
                lines.append(f"{rec}{serial:5d} {nm} {rn:>3s} {och}{ores:4d}    {xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}"
                             f"  1.00  0.00          {el:>2s}")
                serial += 1
        lines.append("TER")
    links_out, bonds = [], []
    for path in paths:
        for a, b in zip(path, path[1:]):
            kind, u, v = edges[frozenset({a, b})]
            if kind == "35":
                bonds.append(((num[u], "O3'"), (num[v], "P")))
                if (a, b) != (u, v):                                         # runs 3'->5' in file order
                    links_out.append(((num[v], "P"), (num[u], "O3'")))
            elif kind == "X":
                bonds.append(((num[u], "P"), (num[v], "O3'")))
                links_out.append(((num[u], "P"), (num[v], "O3'")))
            else:
                bonds.append(((num[u], "P"), (num[v], "O5'")))
                links_out.append(((num[v], "O5'"), (num[u], "P")))
    rname, xyz = {}, {}
    for l in lines:
        if l.startswith(("ATOM", "HETATM")):
            k = (l[21], int(l[22:26]))
            rname[k] = l[17:20].strip()
            xyz[(k, l[12:16].strip())] = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    link_lines = []
    for (k1, n1), (k2, n2) in links_out:
        d = np.linalg.norm(xyz[(k1, n1)] - xyz[(k2, n2)])
        a1, a2 = (f" {n:<3s}" if len(n) < 4 else n for n in (n1, n2))
        link_lines.append(f"LINK        {a1} {rname[k1]:>3s} {k1[0]}{k1[1]:4d}                {a2} {rname[k2]:>3s} "
                          f"{k2[0]}{k2[1]:4d}     1555   1555 {d:5.2f}")
    lres = sorted(k for k, h in hands.items() if h == "L")
    if lres:
        remarks = remarks + ["L-sugar residues (mirror-image nucleotides; residue names as for D): "
                             + "; ".join(f"chain {c} {ranges([r for cc, r in lres if cc == c])}"
                                         for c in dict.fromkeys(c for c, _ in lres))]
    head = [f"REMARK 900 {r}" for r in remarks] + list(extra_head)
    if xatoms:
        head.append("HETNAM     X33 3'-3' PHOSPHODIESTER LINKER PHOSPHATE")
    out.write_text("\n".join(head + link_lines + lines) + "\nEND\n")
    return paths, num, bonds, xyz


def ranges(nums):
    out, start, prev = [], None, None
    for n in sorted(nums) + [None]:
        if start is None:
            start = prev = n
        elif n is not None and n == prev + 1:
            prev = n
        else:
            out.append(f"{start}-{prev}" if prev != start else f"{start}")
            start = prev = n
    return ", ".join(out)


# ----------------------------------------------------------------------------- checks of the written model
def check_model(path: Path, bonds, xyz, pieces, juncs, num):
    """Every inter-residue bond, heavy-atom contacts closer than 2.2 A between residues that are not bonded to
    each other, and the helical axes on both sides of every junction."""
    blen = [float(np.linalg.norm(xyz[a] - xyz[b])) for a, b in bonds]
    keys = list(xyz)
    X = np.array([xyz[k] for k in keys])
    bonded = {frozenset({a[0], b[0]}) for a, b in bonds}
    bad = []
    for i, j in cKDTree(X).query_pairs(2.2):
        ri, rj = keys[i][0], keys[j][0]
        if ri != rj and frozenset({ri, rj}) not in bonded:
            bad.append((f"{ri[0]}{ri[1]}:{keys[i][1]}", f"{rj[0]}{rj[1]}:{keys[j][1]}",
                        float(np.linalg.norm(X[i] - X[j]))))
    axes = []
    for j, J in enumerate(juncs):
        lo, hi = pieces[j], pieces[j + 1]
        run1 = list(reversed(J.e1.pairs[:10]))
        run2 = J.e2.pairs[:10]
        p1, u1, rms1, s1 = helix_axis(lo.ix(), lo.seg.rn, run1)
        p2, u2, rms2, s2 = helix_axis(hi.ix(), hi.seg.rn, run2)
        at = (bp_frame(lo.ix(), lo.seg.rn, J.e1.r3, J.e1.r5)[0] + bp_frame(hi.ix(), hi.seg.rn, J.e2.r3, J.e2.r5)[0]) / 2
        bend, off, skew = axis_relation(p1, u1, p2, u2, at)
        axes.append(dict(bend=bend, sideways=off, closest=skew, repeat=[s1, s2], screw_rms=[rms1, rms2]))
    return dict(bonds=dict(n=len(blen), min=min(blen) if blen else None, max=max(blen) if blen else None),
                contacts=bad, axes=axes)


def backbone_torsions(xyz, r3, r5):
    """epsilon, zeta (3' residue r3), alpha, beta, gamma (5' residue r5) across one O3'-P link (model keys)."""
    g = lambda r, n: xyz[(r, n)]
    try:
        return dict(epsilon=dihedral(g(r3, "C4'"), g(r3, "C3'"), g(r3, "O3'"), g(r5, "P")),
                    zeta=dihedral(g(r3, "C3'"), g(r3, "O3'"), g(r5, "P"), g(r5, "O5'")),
                    alpha=dihedral(g(r3, "O3'"), g(r5, "P"), g(r5, "O5'"), g(r5, "C5'")),
                    beta=dihedral(g(r5, "P"), g(r5, "O5'"), g(r5, "C5'"), g(r5, "C4'")),
                    gamma=dihedral(g(r5, "O5'"), g(r5, "C5'"), g(r5, "C4'"), g(r5, "C3'")))
    except KeyError:
        return {}



# ----------------------------------------------------------------------------- command line
def parse_segment(text):
    path, _, opt = text.partition(":flip")
    if opt or (_ and not text.endswith(":flip")):
        raise JunctionError(f"segment '{text}': write FILE.pdb or FILE.pdb:flip")
    return path, bool(_)


OUT_TAG = "_heterohelix"                             # default OUT: <first segment stem>_heterohelix.pdb
REPORT_TAG = "_heterohelix_junction"                 # report: tmp_file/<OUT stem>_heterohelix_junction.json
MAX_TURN = 90.0                                      # deg: neighbouring base pairs turned further face opposite ways
MIN_PAIRS = 2                                        # base pairs a segment needs to have a chain sense at all
LINK_ATOMS = (("r3", ("O3'", "C3'")), ("r5", ("C4'", "C5'", "O5'")))   # what the links bond at each end residue
KIND_WORDS = {"continuous": "continuous (two ordinary 3'-5' links)",
              "bowtie": "bowtie (a 3'-3' link through an X33 phosphate and a 5'-5' link)"}
SENSE_WORDS = {"continuous": "matched", "bowtie": "opposite"}


def default_output_path(first_segment) -> Path:
    """<first segment stem>_heterohelix.pdb in the first segment's folder."""
    p = Path(first_segment).expanduser().resolve()
    return p.with_name(f"{p.stem}{OUT_TAG}.pdb")


def default_report_path(out) -> Path:
    """tmp_file/<OUT stem>_heterohelix_junction.json next to OUT.  When the stem already ends in _heterohelix (the
    default output), tmp_file/<OUT stem>_junction.json, so the word is not written twice."""
    out = Path(out)
    tag = REPORT_TAG[len(OUT_TAG):] if out.stem.endswith(OUT_TAG) else REPORT_TAG
    return out.parent / "tmp_file" / f"{out.stem}{tag}.json"


def _res_label(r):
    return f"{r[0]}{r[1]}"


def _missing_atom(exc: KeyError) -> str:
    key = exc.args[0] if exc.args else "?"
    if isinstance(key, tuple) and len(key) == 3:
        return f"residue {key[0]}{key[1]} has no atom {key[2]}"
    return f"an atom is missing ({key})"


def _sugar_words(hands) -> str:
    hs = sorted(set(hands) - {"?"})
    return {"D": "D sugars", "L": "L sugars", "DL": "D and L sugars"}.get("".join(hs), "sugars not recognised")


def segment_sense(seg: Segment) -> dict:
    """The chain sense along one segment, checked before any fitting.

    Every base pair (strand A's i-th residue with strand B's i-th from the other end) gets its frame with y towards
    strand A's residue.  In a duplex whose pairs all face the same way the z axes of neighbouring pairs turn by a few
    degrees; a turn of more than MAX_TURN (90 deg) means one pair is turned over.  Also the sugars and the chain sense
    at both ends (as End measures it)."""
    A, B = (seg.res[c] for c in seg.chains)
    n = seg.n
    # one base pair has no direction: End would compare the end pair with itself and always find -1
    short = (f"too short: at least {MIN_PAIRS} base pairs are needed to measure the chain sense"
             if n < MIN_PAIRS else None)
    zs, missing = [], []
    for i in range(n):
        try:
            zs.append(bp_frame(seg.ix, seg.rn, A[i], B[n - 1 - i])[1][2])
        except KeyError as exc:
            zs.append(None)
            missing.append(_missing_atom(exc))
    turns, bad = [], []
    for i in range(n - 1):
        if zs[i] is None or zs[i + 1] is None:
            continue
        turn = float(np.degrees(np.arccos(np.clip(zs[i] @ zs[i + 1], -1, 1))))
        turns.append(turn)
        if turn > MAX_TURN:
            bad.append(dict(pairs=[f"{_res_label(A[i])}-{_res_label(B[n - 1 - i])}",
                                   f"{_res_label(A[i + 1])}-{_res_label(B[n - 2 - i])}"], turn=turn))
    ends = {}
    if not missing and not short:
        for which in ("bottom", "top"):
            try:
                ends[which] = End(seg, which).sense
            except KeyError as exc:
                missing.append(_missing_atom(exc))
    return dict(file=str(seg.path), pairs=n, sugars=_sugar_words(seg.hand.values()), sense=ends,
                max_turn=max(turns) if turns else None, uniform=not bad and not missing and not short, bad=bad,
                missing=missing, short=short)


def segment_text(info: dict) -> str:
    """One plain line for a segment: residues per strand, sugars, chain sense at both ends, uniform or not."""
    s = info["sense"]
    text = f"{info['pairs']} residue{'' if info['pairs'] == 1 else 's'} per strand, {info['sugars']}"
    if "bottom" in s and "top" in s:
        text += f", chain sense {s['bottom']:+d} at the bottom and {s['top']:+d} at the top"
    if info.get("short"):
        return text + "; " + info["short"]
    if info["missing"]:
        return text + "; " + "; ".join(dict.fromkeys(info["missing"]))
    if info["bad"]:
        return text + "; chain sense NOT uniform: " + "; ".join(
            f"pairs {b['pairs'][0]} and {b['pairs'][1]} turn {b['turn']:.0f} deg" for b in info["bad"])
    if info["max_turn"] is not None:
        text += f"; neighbouring pairs turn at most {info['max_turn']:.1f} deg (uniform)"
    return text


def junction_ends(lo: Piece, hi: Piece, nwin=NWIN):
    """The two ends a junction joins: the exit end of the lower piece and the entry end of the upper one."""
    return End(lo.seg, lo.exit, nwin), End(hi.seg, hi.entry, nwin)


def natural_link(e1: End, e2: End) -> str:
    return "continuous" if e1.sense == e2.sense else "bowtie"


def end_words(k, piece: Piece, e: End, detail=True) -> str:
    """'segment k (FILE[, upside down]) END' and, with detail, the end pair's residues, sugars and chain sense."""
    head = f"segment {k} ({piece.seg.label()}{', upside down' if piece.flip else ''}) {e.which}"
    if not detail:
        return head
    return (f"{head} (3' end {_res_label(e.r3)} {e.rn[e.r3]}, 5' end {_res_label(e.r5)} {e.rn[e.r5]}; "
            f"{e.hand}-sugars; chain sense {e.sense:+d})")


def missing_link_atoms(*ends: End) -> list:
    """The atoms the links of a junction bond, looked up before any fitting (for both link kinds): O3' and C3' of
    each 3'-end residue, C4', C5' and O5' of each 5'-end residue.  Returns one plain line per missing atom."""
    out = []
    for e in ends:
        for role, names in LINK_ATOMS:
            r = getattr(e, role)
            out += [f"residue {_res_label(r)} of {e.seg.label()} has no atom {n}; the links need it"
                    for n in names if (r[0], r[1], n) not in e.seg.ix]
    return out


def plan_model(entries, separate=False, cache=None):
    """What the tool will do with these segments, worked out from the files alone (no fitting); the GUI shows it.

    entries: [(path, upside_down)] in stacking order.  Returns dict(segments=[dict(text, error)], junctions=[text],
    junction_errors=[bool], sharing=[text]).  cache: an optional dict (path, mtime, size) -> Segment or error text,
    kept by the caller."""
    cache = {} if cache is None else cache
    segs = []
    out = dict(segments=[], junctions=[], junction_errors=[], sharing=[])
    for path, flip in entries:
        path = str(path).strip()
        if not path:
            segs.append(None)
            out["segments"].append(dict(text="", error=False))
            continue
        p = Path(path).expanduser()
        if not p.is_file():
            segs.append(None)
            out["segments"].append(dict(text="file not found", error=True))
            continue
        st = p.stat()
        key = (str(p.resolve()), st.st_mtime_ns, st.st_size)
        if key not in cache:
            try:
                seg = Segment(p.resolve())
                cache[key] = (seg, segment_sense(seg))
            except JunctionError as exc:
                msg = str(exc)
                cache[key] = (None, msg.split(": ", 1)[1] if msg.startswith(str(p.resolve())) else msg)
            except Exception as exc:  # noqa: BLE001 - an unreadable file is shown, not raised
                cache[key] = (None, f"not a usable duplex ({type(exc).__name__}: {exc})")
        seg, info = cache[key]
        if seg is None:
            segs.append(None)
            out["segments"].append(dict(text=info, error=True))
            continue
        segs.append(Piece(seg, bool(flip)) if info["uniform"] else None)
        out["segments"].append(dict(text=segment_text(info), error=not info["uniform"]))
    juncs, keys = [], []
    for j in range(len(entries) - 1):
        lo, hi = segs[j], segs[j + 1]
        if lo is None or hi is None:
            out["junctions"].append(f"Junction {j + 1}: needs two usable segments")
            out["junction_errors"].append(False)
            juncs.append(None)
            keys.append(None)
            continue
        try:
            e1, e2 = junction_ends(lo, hi)
            miss = missing_link_atoms(e1, e2)
            if miss:
                raise JunctionError("; ".join(miss))
            J = Junction(e1, e2)
        except (KeyError, JunctionError) as exc:
            out["junctions"].append(f"Junction {j + 1}: {_missing_atom(exc) if isinstance(exc, KeyError) else exc}")
            out["junction_errors"].append(True)
            juncs.append(None)
            keys.append(None)
            continue
        out["junctions"].append(f"Junction {j + 1}: {end_words(j + 1, lo, e1, False)} -> "
                                f"{end_words(j + 2, hi, e2, False)}; chain senses {e1.sense:+d} and {e2.sense:+d}: "
                                f"{SENSE_WORDS[J.link]} -> {KIND_WORDS[J.link]}")
        out["junction_errors"].append(False)
        juncs.append(J)
        keys.append((str(lo.seg.path.resolve()), e1.which, str(hi.seg.path.resolve()), e2.which))
    # sharing, as the build does it: identical junctions are one; others are grouped with same_junction
    refs = []                                          # (junction number, Junction) of each group's first member
    first_of = {}
    for j, (J, key) in enumerate(zip(juncs, keys), 1):
        if J is None:
            continue
        if key in first_of:
            out["sharing"].append(f"Junction {j} is junction {first_of[key]} again (same files and ends): "
                                  f"fitted once.")
            continue
        first_of[key] = j
        for m, Jm in refs:
            how = same_junction(Jm, J)
            if how:
                turned = ", turned 180 deg" if how == "reversed" else ""
                out["sharing"].append(f"Junction {j} joins the same end pairs as junction {m}{turned}: "
                                      + ("fitted on its own (--separate)." if separate else "fitted together."))
                if separate:
                    refs.append((j, J))
                break
        else:
            refs.append((j, J))
    if sum(J is not None for J in juncs) > 1 and not out["sharing"]:
        out["sharing"].append("Each junction joins different end pairs and is fitted on its own.")
    return out


def output_numbering(paths):
    """(piece, chain, resi) or ('X33', j) -> (output chain, residue number), as write_model numbers them."""
    num = {}
    for ch, path in zip("ABCDEFGHIJKLMNOPQRSTUVWXYZ", paths):
        for i, item in enumerate(path, 1):
            num[item] = (ch, i)
    return num


def _clean(value) -> str:
    return str(value).replace("\n", " ").replace("\r", " ").strip()


def provenance_remark_lines(command=None):
    """REMARK 950 RE_SCRIPT SOFTWARE, COMMAND and OUTPUT_STAGE, as the other RE_Helix tools write them."""
    lines = [f"{REMARK_PREFIX} SOFTWARE name={SOFTWARE_NAME} version={VERSION} developer={SOFTWARE_DEVELOPER}"]
    if command:
        lines.append(f"{REMARK_PREFIX} COMMAND text={_clean(command)}")
    lines.append(f"{REMARK_PREFIX} OUTPUT_STAGE name={SOFTWARE_NAME}")
    return lines


def junction_residues(pieces, juncs, num):
    """Per junction, the residues of each link in output labels (chain:resi:name): continuous, the 3' end and the
    5' end of both 3'-5' links; bowtie, the two 3' ends with the X33 between them and the two 5' ends."""
    def lab(j, r):
        ch, i = num[(j, *r)]
        return f"{ch}:{i}:{pieces[j].seg.rn[r]}"
    out = []
    for j, J in enumerate(juncs):
        if J.link == "bowtie":
            ch, i = num[("X33", j + 1)]
            out.append([("3to3", "end3_1,x33,end3_2", [lab(j, J.e1.r3), f"{ch}:{i}:X33", lab(j + 1, J.e2.r3)]),
                        ("5to5", "end5_1,end5_2", [lab(j, J.e1.r5), lab(j + 1, J.e2.r5)])])
        else:
            out.append([("3to5", "end3,end5", [lab(j, J.e1.r3), lab(j + 1, J.e2.r5)]),
                        ("3to5", "end3,end5", [lab(j + 1, J.e2.r3), lab(j, J.e1.r5)])])
    return out


def junction_remark_lines(juncs, residues):
    """REMARK 950 RE_SCRIPT JUNCTION lines, one per link (two per junction); Get Phenix Restraints frees the sugars
    of these residues in its _junctions.params."""
    lines = []
    for j, (J, links) in enumerate(zip(juncs, residues), 1):
        for link, roles, res in links:
            lines.append(f"{REMARK_PREFIX} JUNCTION op={j} kind={J.link} link={link} roles={roles} "
                         f"residues={','.join(res)} segments={j},{j + 1}")
    return lines


def chain_sense_after(paths, num, xyz, pieces):
    """The chain sense of the written model.  Along every output chain, each paired nucleotide (X33 skipped) gets its
    pair frame with y towards that chain's residue.  Two things are measured:
      the turn between the z axes of neighbouring pairs, inside segments and across junctions (every pair of a
        correct model faces the same way, so every turn stays far below MAX_TURN, 90 deg);
      the cos between each z axis and the helix axis of its segment, oriented in the chain's file order (which way
        the chain runs through its pairs: about +1 in every segment of a B-Z bowtie model, about -1 in every segment
        of an L-Z continuous model).  A chain whose pairs face along it in one segment and against it in another
        is not matched either."""
    ix = {(k[0][0], k[0][1], k[1]): v for k, v in xyz.items()}
    rn = {}
    partner = {}
    for pi, f in enumerate(pieces):
        A, B = (f.seg.res[c] for c in f.seg.chains)
        for i in range(f.seg.n):
            a, b = (pi, *A[i]), (pi, *B[f.seg.n - 1 - i])
            partner[a], partner[b] = b, a
            rn[num[a]], rn[num[b]] = f.seg.rn[A[i]], f.seg.rn[B[f.seg.n - 1 - i]]
    chains, junctions, errors = [], {}, []
    for path in paths:
        nodes = [n for n in path if n[0] != "X33"]
        if not nodes:
            continue
        ch = num[nodes[0]][0]
        z = [bp_frame(ix, rn, num[n], num[partner[n]])[1][2] for n in nodes]
        within = []
        for k in range(len(nodes) - 1):
            u, v = nodes[k], nodes[k + 1]
            c = float(np.clip(z[k] @ z[k + 1], -1, 1))
            ang = float(np.degrees(np.arccos(c)))
            where = f"chain {ch} {ch}{num[u][1]} -> {ch}{num[v][1]}"
            if u[0] == v[0]:
                within.append((ang, c))
                if ang > MAX_TURN:
                    errors.append(f"chain sense not matched at {where} (inside segment {u[0] + 1}): {ang:.1f} deg")
            else:
                j = min(u[0], v[0]) + 1
                junctions.setdefault(j, []).append(dict(chain=ch, step=f"{ch}{num[u][1]} -> {ch}{num[v][1]}",
                                                        angle=ang, cos=c))
                if ang > MAX_TURN:
                    errors.append(f"chain sense not matched at {where} (junction {j}): {ang:.1f} deg")
        # the runs of the chain through each segment, and the cos of every z with that run's helix axis
        runs, k0 = [], 0
        for k in range(1, len(nodes) + 1):
            if k == len(nodes) or nodes[k][0] != nodes[k0][0]:
                runs.append(nodes[k0:k])
                k0 = k
        segs, at = [], 0
        for run in runs:
            zr = z[at: at + len(run)]
            at += len(run)
            row = dict(segment=run[0][0] + 1, residues=f"{ch}{num[run[0]][1]}-{ch}{num[run[-1]][1]}", axis_cos=None,
                       faces=None, repeat=None, screw_rms=None)
            if len(run) >= MIN_PAIRS:
                _, u_ax, rms, rep = helix_axis(ix, rn, [(num[n], num[partner[n]]) for n in run])
                cs = [float(np.clip(zz @ u_ax, -1, 1)) for zz in zr]
                row.update(axis_cos=[min(cs), max(cs)], repeat=rep, screw_rms=rms,
                           faces="along" if min(cs) > 0 else "against" if max(cs) < 0 else "mixed")
            segs.append(row)
        faces = {s["faces"] for s in segs if s["faces"]}
        if len(faces) > 1 or "mixed" in faces:
            where = lambda f: ", ".join(str(s["segment"]) for s in segs if s["faces"] == f)
            words = {"along": "along it", "against": "against it", "mixed": "both ways"}
            parts = [f"{words[f]} in segment(s) {where(f)}" for f in ("along", "against", "mixed") if f in faces]
            errors.append(f"chain sense not matched along chain {ch}: its pairs face " + "; ".join(parts))
        chains.append(dict(chain=ch, nucleotides=len(nodes), within_max=max((a for a, _ in within), default=0.0),
                           within_cos=[min((c for _, c in within), default=1.0),
                                       max((c for _, c in within), default=1.0)],
                           faces=faces.pop() if len(faces) == 1 else ("mixed" if faces else None), segments=segs))
    jrows = [dict(j=j, steps=steps, max_angle=max(s["angle"] for s in steps)) for j, steps in sorted(junctions.items())]
    return dict(chains=chains, junctions=jrows, errors=errors, matched=not errors, max_turn=MAX_TURN)


def chain_sense_after_lines(after) -> list:
    """The run-log lines of the check after building."""
    lines = []
    words = {"along": "z points along the chain (file order) in every segment",
             "against": "z points against the chain (file order) in every segment",
             "mixed": "z does NOT point the same way along the chain in every segment"}
    for c in after["chains"]:
        cos = ", ".join(f"segment {s['segment']} {s['axis_cos'][0]:+.3f} to {s['axis_cos'][1]:+.3f}"
                        for s in c["segments"] if s["axis_cos"])
        lines.append(f"  chain {c['chain']} ({c['nucleotides']} nucleotides): {words.get(c['faces'], 'no segment run')}"
                     f" (cos with the helix axis: {cos}); turn between neighbouring pairs inside segments at most "
                     f"{c['within_max']:.1f} deg")
    for jr in after["junctions"]:
        lines.append(f"  junction {jr['j']}: turn between neighbouring pairs "
                     + ", ".join(f"{s['step']} {s['angle']:.1f} deg" for s in jr["steps"]))
    lines += [f"  Error: {e}" for e in after["errors"]]
    lines.append("  matched: every base pair faces the same way along its chain (no turn above 90 deg)"
                 if after["matched"] else "  chain sense NOT matched")
    return lines


def _grid_values(text, default, integral=False):
    if text in (None, ""):
        return tuple(default)
    vals = []
    for t in str(text).replace(";", ",").split(","):
        if t.strip():
            v = float(t)
            vals.append(int(v) if integral and v.is_integer() else v)
    if not vals:
        raise JunctionError(f"empty search grid '{text}'")
    return tuple(vals)


def build_cli_command(segments, out=None, link="auto", cores=None, report=None, reuse=False, rise=RISE,
                      separate=False, twists=None, shifts=None, candidates=None, executable="python3", script=None,
                      reuse_from=None):
    """The command line that repeats a run (shell-quoted)."""
    parts = [executable, script or Path(__file__).name, *[str(s) for s in segments]]
    if out:
        parts += ["-o", str(out)]
    if link and link != "auto":
        parts += ["--link", link]
    if cores:
        parts += ["--cores", str(cores)]
    if rise is not None and float(rise) != RISE:
        parts += ["--rise", f"{float(rise):g}"]
    if separate:
        parts.append("--separate")
    if reuse:
        parts.append("--reuse")
    if reuse_from:
        parts += ["--reuse-from", str(reuse_from)]
    if report:
        parts += ["--report", str(report)]
    for flag, val in (("--twists", twists), ("--shifts", shifts), ("--candidates", candidates)):
        if val not in (None, ""):
            parts.append(f"{flag}={val}")                  # attached, so a list like -40,-20 is not read as a flag
    return " ".join(shlex.quote(p) for p in parts)


def build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="build_heterohelix_junction.py", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog="More in the RE_Helix README, section Build Heterohelix Junction Tool.")
    ap.add_argument("segments", nargs="*", help="segment PDB files in stacking order, bottom first (FILE.pdb or "
                                                "FILE.pdb:flip); at least two")
    ap.add_argument("-o", "--out", default=None, help="the joined model (PDB; default <first segment>_heterohelix.pdb "
                                                      "beside the first segment)")
    ap.add_argument("--link", default="auto", help="auto, bowtie or continuous: one value, or one per junction "
                                                   "separated by commas (default auto; a kind the chain senses do "
                                                   "not allow is refused)")
    ap.add_argument("--cores", type=int, default=None, help="processes for the fit (default: all but one)")
    ap.add_argument("--report", default=None, help="the JSON report (default tmp_file/<OUT stem>_heterohelix_junction"
                                                   ".json next to OUT; <first segment>_heterohelix_junction.json for "
                                                   "the default OUT)")
    ap.add_argument("--reuse", action="store_true", help="take the fits from this output's own report (--report, or "
                                                         "the default) when they were made the same way: this tool "
                                                         "version, unchanged segment files, the same search grid and "
                                                         "the same junctions fitted together; otherwise fit again")
    ap.add_argument("--reuse-from", default=None, metavar="REPORT", help="the same as --reuse, with the fits taken "
                                                                         "from an earlier report of this tool")
    ap.add_argument("--rise", type=float, default=RISE, help=f"stacking distance aimed at (A, default {RISE})")
    ap.add_argument("--separate", action="store_true", help="fit every junction on its own, also those that join "
                                                            "the same two end pairs (by default they share one fit)")
    # A smaller search grid, for tests and quick trials only (the default grid is the fit the tool is made for).
    ap.add_argument("--twists", default=None, help=argparse.SUPPRESS)
    ap.add_argument("--shifts", default=None, help=argparse.SUPPRESS)
    ap.add_argument("--candidates", type=int, default=None, help=argparse.SUPPRESS)
    ap.add_argument("--gui", action="store_true", help="open the tool window; segments given here fill its rows")
    ap.add_argument("-v", "--version", action="version", version=f"{TOOL_NAME} V{VERSION}")
    return ap


def main(argv=None):
    """Command line entry: with --gui (or no arguments) the tool window, otherwise a build (returns OUT)."""
    argv = sys.argv[1:] if argv is None else list(argv)
    if not argv:
        return run_gui()
    a = build_arg_parser().parse_args(argv)
    if a.gui:
        return run_gui(a.segments, out=a.out, cores=a.cores, rise=a.rise, separate=a.separate, reuse=a.reuse,
                       reuse_from=a.reuse_from)
    return build_model(a)


# ----------------------------------------------------------------------------- reusing fits
def file_md5(path) -> str:
    return hashlib.md5(Path(path).read_bytes()).hexdigest()


def fit_groups(juncs, share=True, tol=SAME_END_RMS) -> list:
    """The groups fit_model_junctions makes of these junctions (indices; the first of each group is its reference):
    junctions that join the same two end pairs together when share is on, every junction alone when it is off."""
    groups = []
    for i, J in enumerate(juncs):
        for g in groups:
            if share and same_junction(juncs[g[0]], J, tol):
                g.append(i)
                break
        else:
            groups.append([i])
    return groups


def fitting_record(uniq, weights, groups, grid, share) -> dict:
    """What the fits of a run are made from (written to the report): the fitter, the search grid, the distinct
    junctions with how often each occurs, the groups fitted together and an md5 of every segment file.  --reuse and
    --reuse-from take fits only from a report with the same record, so a reused fit is the fit this run would make."""
    files = {}
    for s in uniq:
        for f in (s["seg1"], s["seg2"]):
            if f not in files:
                files[f] = file_md5(f)
    return json.loads(json.dumps(dict(fitter=f"{TOOL_NAME} V{VERSION}", grid=grid, same_end_rms=SAME_END_RMS,
                                      share=bool(share), junctions=uniq, weights=weights, groups=groups, files=files)))


def reusable_fits(path, fitting):
    """(None, fits in the order of fitting['junctions']) when the report at path holds fits made exactly as this run
    would make them; otherwise (the reason in plain words, None)."""
    path = Path(path)
    if not path.is_file():
        return f"no earlier report at {path}", None
    try:
        earlier = json.loads(path.read_text())
    except (OSError, ValueError) as exc:
        return f"{path.name} is not a readable report ({exc})", None
    if not isinstance(earlier, dict) or earlier.get("tool") != TOOL_NAME:
        return f"{path.name} was not written by {TOOL_NAME}", None
    was = earlier.get("fitting")
    if not isinstance(was, dict) or was.get("fitter") != fitting["fitter"]:
        return f"{path.name} was written by another version of {TOOL_NAME}", None
    same = lambda k: json.dumps(was.get(k), sort_keys=True) == json.dumps(fitting[k], sort_keys=True)
    if not same("grid"):
        return f"{path.name} was fitted on another search grid", None
    if not (same("junctions") and same("weights") and same("same_end_rms")):
        return f"{path.name} holds the fits of other junctions", None
    if not same("groups"):
        together = lambda groups: any(len(g) > 1 for g in groups or [])
        if together(was.get("groups")) and not together(fitting["groups"]):
            return (f"{path.name} fitted the junctions that join the same end pairs together, and this run fits each "
                    f"one on its own (--separate)"), None
        if together(fitting["groups"]) and not together(was.get("groups")):
            return (f"{path.name} fitted each junction on its own (--separate), and this run fits the junctions that "
                    f"join the same end pairs together"), None
        return f"{path.name} grouped the junctions differently", None
    changed = [Path(f).name for f, md5 in fitting["files"].items() if (was.get("files") or {}).get(f) != md5]
    if changed:
        return f"segment file {', '.join(changed)} changed after {path.name} was written", None
    key = lambda s: json.dumps(s, sort_keys=True)
    fits = {key(f.get("spec")): f for f in earlier.get("fits", []) if isinstance(f, dict)}
    if not all(key(s) in fits for s in fitting["junctions"]):
        return f"{path.name} holds the fits of other junctions", None
    return None, [fits[key(s)] for s in fitting["junctions"]]


# ----------------------------------------------------------------------------- the build
def build_model(a, log=print):
    """The build itself, from parsed arguments: checks, fits, model, report.  Returns the model's path."""
    if len(a.segments) < 2:
        raise JunctionError("give at least two segments (FILE.pdb or FILE.pdb:flip)")
    segs = [parse_segment(s) for s in a.segments]
    for p, _ in segs:
        if not Path(p).expanduser().is_file():
            raise JunctionError(f"segment file not found: {p}")
    out = Path(a.out).expanduser().resolve() if a.out else default_output_path(segs[0][0])
    report = Path(a.report).expanduser().resolve() if a.report else default_report_path(out)
    reuse_from = getattr(a, "reuse_from", None)
    if out in {Path(p).expanduser().resolve() for p, _ in segs}:
        raise JunctionError(f"the output {out.name} would overwrite an input segment; choose another -o")
    linkopt = [s.strip() for s in a.link.split(",")]
    if len(linkopt) == 1:
        linkopt = linkopt * (len(segs) - 1)
    if len(linkopt) != len(segs) - 1 or any(l not in ("auto", "bowtie", "continuous") for l in linkopt):
        raise JunctionError(f"--link: give auto, bowtie or continuous once, or once per junction ({len(segs) - 1})")
    twists = _grid_values(a.twists, TWISTS, integral=True)
    shifts = _grid_values(a.shifts, SHIFTS)
    n_cand = a.candidates or N_CAND
    grid = dict(twists=[float(v) for v in twists], shifts=[float(v) for v in shifts], candidates=n_cand)
    command = build_cli_command(a.segments, out, a.link, a.cores, a.report, a.reuse, a.rise, a.separate,
                                a.twists, a.shifts, a.candidates, reuse_from=reuse_from)
    log("Equivalent CLI command:")
    log(command)
    pieces = [Piece(get_segment(str(Path(p).expanduser())), flip) for p, flip in segs]
    # ---- checks before fitting: each segment uniform and long enough, each junction's two ends, the link atoms
    log("Chain sense check (before fitting)")
    problems, before = [], dict(segments=[], junctions=[])
    infos = {}
    for k, pc in enumerate(pieces, 1):
        if pc.seg.path not in infos:
            infos[pc.seg.path] = segment_sense(pc.seg)
        info = infos[pc.seg.path]
        log(f"  segment {k}: {pc.seg.label()}{' (upside down)' if pc.flip else ''}: {segment_text(info)}")
        before["segments"].append(dict(segment=k, file=str(pc.seg.path), upside_down=pc.flip,
                                       **{kk: info[kk] for kk in ("pairs", "sugars", "sense", "max_turn", "uniform",
                                                                  "bad", "missing", "short")}))
        if info["short"]:
            problems.append(f"segment {k} ({pc.seg.label()}) has {info['pairs']} base pair"
                            f"{'' if info['pairs'] == 1 else 's'}; at least {MIN_PAIRS} are needed to measure its "
                            f"chain sense")
        elif info["missing"]:
            problems.append(f"segment {k} ({pc.seg.label()}): {'; '.join(dict.fromkeys(info['missing']))}")
        elif info["bad"]:
            problems.append(f"segment {k} ({pc.seg.label()}): the chain sense is not uniform: "
                            + "; ".join(f"pairs {b['pairs'][0]} and {b['pairs'][1]} turn {b['turn']:.1f} deg"
                                        for b in info["bad"])
                            + " (a base pair is turned over); every base pair of a segment must face the same way")
    for j in range(len(pieces) - 1):
        lo, hi = pieces[j], pieces[j + 1]
        if not (infos[lo.seg.path]["uniform"] and infos[hi.seg.path]["uniform"]):
            log(f"  junction {j + 1}: not checked (a segment failed)")
            continue
        e1, e2 = junction_ends(lo, hi)
        natural = natural_link(e1, e2)
        ends = [end_words(j + 1, lo, e1), end_words(j + 2, hi, e2)]
        line = f"  junction {j + 1}: {ends[0]} | {ends[1]}: {SENSE_WORDS[natural]} -> {KIND_WORDS[natural]}"
        refused = linkopt[j] not in ("auto", natural)
        if refused:
            line += f"; --link {linkopt[j]} refused"
            try:
                Junction(e1, e2, linkopt[j])
            except JunctionError as exc:
                problems.append(f"junction {j + 1}: {exc}")
        missing = missing_link_atoms(e1, e2)
        problems += [f"junction {j + 1}: {m}" for m in missing]
        log(line)
        before["junctions"].append(dict(j=j + 1, ends=ends, sense=[e1.sense, e2.sense],
                                        result=f"{SENSE_WORDS[natural]} -> {natural}", link_asked=linkopt[j],
                                        refused=refused, missing_atoms=missing))
    if problems:
        raise JunctionError("refused before fitting:\n  " + "\n  ".join(problems))
    specs, juncs = [], []
    for j in range(len(pieces) - 1):
        lo, hi = pieces[j], pieces[j + 1]
        spec = dict(seg1=str(lo.seg.path.resolve()), end1=lo.exit, seg2=str(hi.seg.path.resolve()), end2=hi.entry,
                    link=linkopt[j], nwin=NWIN, rise=a.rise, gap=list(GAP))
        try:
            J = make_junction(spec)
        except KeyError as exc:
            raise JunctionError(f"junction {j + 1}: {_missing_atom(exc)}; the links need it") from None
        spec["link"] = J.link
        specs.append(spec)
        juncs.append(J)
    # distinct junctions (how often each occurs), reused from an earlier report when asked and made the same way
    uniq, weights = [], []
    for s in specs:
        if s in uniq:
            weights[uniq.index(s)] += 1
        else:
            uniq.append(s)
            weights.append(1)
    key = lambda s: json.dumps(s, sort_keys=True)
    groups = fit_groups([juncs[specs.index(s)] for s in uniq], share=not a.separate)
    fitting = fitting_record(uniq, weights, groups, grid, not a.separate)
    source = Path(reuse_from).expanduser().resolve() if reuse_from else (report if a.reuse else None)
    old = {}
    if source is not None:
        reason, earlier = reusable_fits(source, fitting)
        if reason is None:
            old = {key(s): f for s, f in zip(uniq, earlier)}
            log(f"{len(uniq)} distinct junction(s); reusing the fits from {source}")
        else:
            log(f"not reusing earlier fits: {reason}")
    if not old:
        log(f"{len(uniq)} distinct junction(s); fitting them" + ("" if a.separate else
            " (junctions that join the same two end pairs share one fit)") + "; a few minutes per junction")
        for s, f in zip(uniq, fit_model_junctions(uniq, weights, a.cores, twists=twists, shifts=shifts, n_cand=n_cand,
                                                  share=not a.separate, log=log)):
            old[key(s)] = f
    for f in (old[key(s)] for s in uniq):
        gr = f["group"]
        if len(gr["members"]) > 1 and gr["reference"] != uniq.index(f["spec"]):
            at = [j + 1 for j, s in enumerate(specs) if s == f["spec"]]
            ref_at = [j + 1 for j, s in enumerate(specs) if s == uniq[gr["reference"]]]
            who = f"junction {at[0]} joins" if len(at) == 1 else f"junctions {', '.join(map(str, at))} join"
            log(f"{who} the same two end pairs as junction {ref_at[0]}"
                + (", turned 180 deg" if gr["how"] == "reversed" else "")
                + f" (ends match within {max(gr['end_rms']):.3f} A): fitted together, one step")
    fits = [old[key(s)] for s in specs]
    # the model
    place_pieces(pieces, juncs, fits)
    links = build_links(pieces, juncs)
    remarks = [f"Built by RE_Helix {TOOL_NAME} V{VERSION} (as built, not minimized): "
               + " | ".join(f"{p.seg.label()}{' (flipped)' if p.flip else ''}" for p in pieces)]
    for j, (J, f) in enumerate(zip(juncs, fits)):
        st = canonical_step(f["step"])
        remarks.append(f"Junction {j + 1}: {J.e1.seg.label()} {J.e1.which} -> {J.e2.seg.label()} {J.e2.which}, "
                       f"{J.link}; step " + " ".join(f"{v:.2f}" for v in st)
                       + " (tilt roll twist deg, shift slide rise A)")
    num = output_numbering(trace_chains(pieces, juncs)[0])
    jres = junction_residues(pieces, juncs, num)
    head950 = provenance_remark_lines(command) + junction_remark_lines(juncs, jres)
    out.parent.mkdir(parents=True, exist_ok=True)
    paths, num, bonds, xyz = write_model(pieces, juncs, fits, links, out, remarks, extra_head=head950)
    checks = check_model(out, bonds, xyz, pieces, juncs, num)
    rows = []
    for j, (J, f, br, ax) in enumerate(zip(juncs, fits, links, checks["axes"])):
        m = f["metrics"]
        row = dict(j=j + 1, link=J.link, ends=[J.e1.describe(), J.e2.describe()], step=canonical_step(f["step"]),
                   planes=m["ang"], rise=m["rise"], overlap=m["ov"], closest=m["dmin"], axes=ax,
                   phosphates=br["geometry"],
                   residues=[dict(link=link, roles=roles.split(","), residues=res) for link, roles, res in jres[j]])
        if J.link == "continuous":
            row["torsions"] = [backbone_torsions(xyz, num[(j + l["three"][0], *l["three"][1])],
                                                 num[(j + l["five"][0], *l["five"][1])]) for l in br["links"]]
            row["joins"] = [f"{num[(j + l['three'][0], *l['three'][1])][0]}{num[(j + l['three'][0], *l['three'][1])][1]}"
                            f" O3' -> {num[(j + l['five'][0], *l['five'][1])][0]}"
                            f"{num[(j + l['five'][0], *l['five'][1])][1]} P" for l in br["links"]]
        rows.append(row)
        g = br["geometry"]
        what = ("X33 P-O3' {:.2f}/{:.2f}, 5'-5' P-O5' {:.2f}/{:.2f}".format(*g[0]["bonds"], *g[1]["bonds"])
                if J.link == "bowtie" else
                "P-O3'/P-O5' {:.2f}/{:.2f} and {:.2f}/{:.2f}".format(*g[0]["bonds"], *g[1]["bonds"]))
        st = row["step"]
        log(f"junction {j + 1} ({J.link}): tilt {st[0]:.1f} roll {st[1]:.1f} twist {st[2]:.1f} deg, shift "
            f"{st[3]:.2f} slide {st[4]:.2f} rise {st[5]:.2f} A | planes {m['ang']:.1f} deg, rise {m['rise']:.2f} A, "
            f"overlap {m['ov']:.1f} A2, closest {m['dmin']:.2f} A | {what} | axes bend {ax['bend']:.1f} deg, "
            f"{ax['sideways']:.1f} A sideways")
    # ---- chain sense check after building
    after = chain_sense_after(paths, num, xyz, pieces)
    log("Chain sense check (after building): every base pair along each chain, its frame with y towards that "
        "chain's residue")
    for line in chain_sense_after_lines(after):
        log(line)
    seg_of = {}
    for item, (ch, r) in num.items():
        seg_of[f"{ch}{r}"] = "X33" if item[0] == "X33" else f"{item[0] + 1}:{item[1]}{item[2]}"
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(json.dumps(dict(tool=TOOL_NAME, version=VERSION, command=command, model=str(out),
                                      segments=[[str(p.seg.path), p.flip] for p in pieces], grid=grid,
                                      fitting=fitting, junctions=rows,
                                      fits=[old[json.dumps(s, sort_keys=True)] for s in uniq],
                                      checks=dict(bonds=checks["bonds"], contacts=checks["contacts"]),
                                      chain_sense=dict(before=before, after=after),
                                      segment_of=seg_of, chains=[len(p) for p in paths]),
                                 indent=1, default=lambda o: o.tolist() if isinstance(o, np.ndarray) else float(o)))
    c = checks["bonds"]
    log(f"wrote {out}: chains {', '.join(f'{ch} {len(p)}' for ch, p in zip('ABCDEFGH', paths))} residues; "
        f"{c['n']} inter-residue bonds {c['min']:.2f}-{c['max']:.2f} A; {len(checks['contacts'])} contacts < 2.2 A")
    log(f"report {report}")
    if after["errors"]:
        raise JunctionError("chain sense not matched in the written model (kept for inspection): "
                            + "; ".join(after["errors"]))
    return out


def _flush_stdout():
    try:
        sys.stdout.flush()                             # so the error comes after the log it refers to
    except Exception:  # noqa: BLE001 - nobody reads it any more
        pass


def cli(argv=None) -> int:
    """main() as a command: a clear one-line error and exit status 1 when the build stops."""
    try:
        result = main(argv)
    except JunctionError as exc:
        _flush_stdout()
        print(f"Error: {exc}", file=sys.stderr, flush=True)
        return 1
    except KeyboardInterrupt:
        _flush_stdout()
        print("Error: stopped", file=sys.stderr, flush=True)
        return 1
    except Exception as exc:  # noqa: BLE001 - shown as one clear line after the traceback
        _flush_stdout()
        traceback.print_exc()
        print(f"Error: {type(exc).__name__}: {exc}", file=sys.stderr, flush=True)
        return 1
    return result if isinstance(result, int) else 0


# ----------------------------------------------------------------------------- the tool window
HELP_TEXTS = {
    "tool": """Build Heterohelix Junction joins two or more duplex segments end to end and builds the junctions between them. A segment can be any duplex form, such as B-DNA, A-RNA or Z-DNA, with D or L sugars, as BNP_NA writes it. Each segment stays a rigid copy of its file; only the phosphates of the junctions are built.

The chain sense rule
Every end base pair gets a frame: x towards the major groove, y towards the nucleotide that ends there with its 3' end, and z = x cross y. In B-DNA and A-RNA, z points out of the segment (chain sense +1). In Z-DNA, z points back into it (chain sense -1). A mirror image flips the sign, so L-DNA and L-RNA are -1, and L-Z-DNA is +1.
  Matched (same chain sense): a continuous junction, two ordinary 3'-5' links.
  Opposite chain sense: a bowtie junction, one 3'-3' link through an X33 phosphate and one 5'-5' link. A strand that ran straight on would meet the next base pairs from their other face.

Checks
Before any fitting, the tool lists the chain senses of every junction. It also checks that every base pair of a segment faces the same way, that each segment has at least two base pairs, and that the residues each link joins have the atoms it bonds. Anything wrong stops the run before fitting.
After building, it walks each output chain. Neighbouring base pairs may not turn by more than 90 deg, and the pairs must face the same way along the chain in every segment (the cos between each pair's z axis and its segment's helix axis keeps one sign).

How a junction is fitted
1. The two end pairs get frames whose z axes point from the lower segment into the upper one.
2. The upper end is placed on the lower one by one base-pair step: tilt, roll, twist, shift, slide and rise.
3. A step is scored on atoms closer than 2.9 A, the gaps the two links must span (each with a quick phosphate), how parallel and how far apart (3.5 A) the two end pairs are, and how much their bases overlap.
4. Searches start from a grid of twists (-180 to 160 deg in 20-deg steps), shifts and slides: 162 starts, run from both ends of the junction. A second pass starts again from small nudges around the 12 best results.
5. The 16 best different results get real phosphates, and the one with the lowest total is kept. Neither the search nor the phosphates depend on where the files sit in space.
Junctions that join the same two end pairs (also turned 180 deg) are fitted together and get the same step.
The score changes little with twist, so read the fitted twist as uncertain by several degrees; the report lists every candidate.

Inputs: the segments in stacking order, bottom first; tick "upside down" to turn one over.
Outputs: the joined model (chains A, B, ... numbered from 1, with LINK records and REMARK lines that say how it was built) and a JSON report in tmp_file/ beside it. The model is as built, not minimized.

Time: about 5 minutes per junction on 9 processes. A junction that repeats another one (same files, same ends) is fitted once.

On the command line the same build is
  python3 build_heterohelix_junction.py SEG1.pdb SEG2.pdb[:flip] ... -o OUT.pdb""",
    "segments": """The duplex segments, in stacking order: row 1 is the bottom segment, the last row the top one. Use at least two; a segment can appear more than once.

Each file must be a duplex as BNP_NA writes it: two chains, strand A running 5'->3' from the bottom end to the top end and strand B back down, with strand A's i-th residue paired to strand B's i-th from the other end. B-, A- and Z-forms work, DNA or RNA, with D or L sugars (L-forms keep the D residue names). A segment needs at least two base pairs.

Under each file the tool shows what it read: residues per strand, the sugars, the chain sense at the bottom and at the top end, and whether every base pair of the segment faces the same way. A file that cannot be used says why, in red.

Upside down turns a segment over before it is stacked: the junction below it then meets its top end, and the one above it its bottom end. On the command line this is FILE.pdb:flip.

Select a row with its number to remove it or move it up or down.""",
    "junctions": """Worked out from the segments; nothing to set here.

Each row names the two ends a junction joins (segment number, file and the file's end) and their chain senses:
  matched (the same sign) -> continuous: two ordinary 3'-5' links;
  opposite signs -> bowtie: a 3'-3' link through an X33 phosphate and a 5'-5' link.

The chain senses decide the kind. On the command line --link bowtie or --link continuous can name it, but only the kind the chain senses allow; the other is refused before fitting, because it would turn one base pair over. A junction whose residues lack an atom the links bond is shown in red.

Junctions that join the same two end pairs, in the same order or turned 180 deg, are fitted together and get the same step. The rows below the junctions say which.""",
    "rise": """The stacking distance the fit aims at between the two end pairs of a junction, in A. Default 3.5 A; distances within 0.15 A of it cost nothing.""",
    "processes": """How many processes share the searches. The default is all processors but one. More processes finish sooner; the result is the same.""",
    "share": """On (the default): junctions of the model that join the same two end pairs, in the same order or turned 180 deg about the duplex dyad, are fitted together. Their searches are pooled and all of them get the step with the lowest summed total, so equal junctions come out equal.

Off: each junction is fitted on its own (--separate). Junctions that use the same files and the same ends are still fitted once.""",
    "reuse": """Off (the default): every Run fits the junctions.

On: the fits are taken from an earlier report instead of fitting again. This works only when that report holds fits of exactly these junctions, made the same way: by this version of the tool, from the same segment files (unchanged), on the same search grid, with the same junctions fitted together. Otherwise the junctions are fitted again, and the run log says why.

The earlier report starts as this output's own report, so a second Run with the same output name writes the model again without fitting. Browse... picks another report, for example the report of the same build under another output name. The new report is always written beside the new output.""",
    "report": """The JSON report is written automatically beside the output PDB, in tmp_file/<output>_heterohelix_junction.json; for the default output, <first segment>_heterohelix.pdb, it is tmp_file/<first segment>_heterohelix_junction.json. It holds every search, the candidates, the chosen steps, the phosphates, the helical axes at each junction, the chain sense checks, and what the fits were made from.""",
    "output": """The joined model. Default: <first segment>_heterohelix.pdb in the first segment's folder; it follows the first segment until you type or choose another name. A name without a folder goes in the first segment's folder, and an empty field goes back to the default.

Chains are named A, B, ... and numbered from 1. LINK records are written for every bond that is not an ordinary 5'->3' step in file order (X33 links, 5'-5' links, inverted steps). REMARK 900 lines record the segments, the junctions and the L-sugar residues; REMARK 950 RE_SCRIPT lines record the tool, the command and the junction residues, which Get Phenix Restraints reads. The model is as built, not minimized.""",
    "log": """Run starts the build in its own process; this log shows its output as it comes: the exact command, the chain sense check, the search progress, the result for each junction and the check after building. When a run fails, the footer shows the error and its first reason.""",
}


def run_gui(segments=None, out=None, cores=None, rise=RISE, separate=False, reuse=False, reuse_from=None) -> int:
    """The tool window.  segments: initial rows (FILE.pdb or FILE.pdb:flip); a Run starts the command line in its
    own process (sys.executable this_file ...) and streams its output into the run log."""
    try:
        import queue
        import signal
        import subprocess
        import threading
        import tkinter as tk
        from tkinter import filedialog, messagebox, ttk
    except Exception as exc:  # pragma: no cover - depends on the Python build
        print(f"Error: the tool window needs tkinter ({exc}).", file=sys.stderr)
        return 1
    try:
        from re_helix_lib import azbmost_ui as ui
    except ImportError:  # pragma: no cover - direct script execution fallback
        import azbmost_ui as ui
    try:
        root = tk.Tk()
    except tk.TclError as exc:
        print(f"Error: cannot open the tool window ({exc}).", file=sys.stderr)
        return 1
    root.title(f"{TOOL_NAME} V{VERSION}")
    apply_optional_icon(root, __file__)
    header = ui.tool_window(root, TOOL_NAME, tool="heterohelix_junction", version=f"V{VERSION}")
    root.geometry("1300x880")
    root.minsize(1000, 620)

    def show_help(title, key):
        messagebox.showinfo(title, HELP_TEXTS[key], parent=root)

    def help_button(parent, title, key, **kw):
        return ui.HelpButton(parent, command=lambda: show_help(title, key), **kw)

    help_row = ttk.Frame(header.right, style="Header.TFrame")
    help_row.pack(anchor="e")
    ttk.Label(help_row, text="How it works", style="Muted.TLabel").pack(side="left", padx=(0, 8))
    help_button(help_row, TOOL_NAME, "tool", size=22).pack(side="left")

    cpu_default = max(1, (os.cpu_count() or 2) - 1)
    rows = []
    selected = tk.IntVar(value=0)
    rise_var = tk.StringVar(value=f"{float(rise):g}")
    cores_var = tk.StringVar(value=str(cores or cpu_default))
    share_var = tk.BooleanVar(value=not separate)
    reuse_var = tk.BooleanVar(value=bool(reuse or reuse_from))
    out_var = tk.StringVar(value=out or "")
    report_var = tk.StringVar()
    reuse_path_var = tk.StringVar(value=str(reuse_from or ""))
    status_var = tk.StringVar(value="Choose the segments, bottom first, and press Run.")
    q = queue.Queue()
    state = dict(proc=None, running=False, stopped=False, auto_out=None, auto_reuse=None, pending=None, cache={},
                 error_lines=[], error_open=False, target=None, mirror=True, wrap=520)

    footer = ui.Footer(root, textvariable=status_var)
    footer.status.configure(wraplength=620, justify="left")
    footer.pack(side="bottom", fill="x")

    body = ttk.Frame(root, padding=(18, 12, 18, 12))
    body.pack(fill="both", expand=True)
    body.columnconfigure(0, weight=3, minsize=600)
    body.columnconfigure(1, weight=2, minsize=360)
    body.rowconfigure(0, weight=1)
    scroll = ui.ScrollFrame(body)
    scroll.grid(row=0, column=0, sticky="nsew", padx=(0, 16))
    left = scroll.body
    left.columnconfigure(0, weight=1)

    # ---- segments (the inputs)
    seg_card = ui.Card(left, "Segments, bottom to top", help_command=lambda: show_help("Segments", "segments"))
    seg_card.grid(row=0, column=0, sticky="ew", pady=(0, 10), padx=(0, 6))
    seg_card.columnconfigure(0, weight=1)
    rows_frame = ttk.Frame(seg_card)
    rows_frame.grid(row=0, column=0, sticky="ew")
    rows_frame.columnconfigure(1, weight=1)
    seg_bar = ttk.Frame(seg_card)
    seg_bar.grid(row=1, column=0, sticky="w", pady=(10, 0))

    # ---- output, with the report path it gives (greyed: worked out from the output)
    out_card = ui.Card(left, "Output")
    out_card.grid(row=1, column=0, sticky="ew", pady=(0, 10), padx=(0, 6))
    out_card.columnconfigure(1, weight=1)
    ttk.Label(out_card, text="Output PDB", style="Field.TLabel").grid(row=0, column=0, sticky="w", padx=(0, 12), pady=3)
    out_entry = ttk.Entry(out_card, textvariable=out_var)
    out_entry.grid(row=0, column=1, sticky="ew", pady=3)

    def save_as():
        cur = resolved_output()
        path = filedialog.asksaveasfilename(
            parent=root, title="Save the joined model as", defaultextension=".pdb",
            initialdir=str(cur.parent) if cur else None, initialfile=cur.name if cur else None,
            filetypes=[("PDB files", "*.pdb"), ("All files", "*.*")])
        if path:
            out_var.set(path)
            show_end(out_entry)

    ttk.Button(out_card, text="Save as...", command=save_as).grid(row=0, column=2, padx=(8, 0), pady=3)
    help_button(out_card, "Output PDB", "output").grid(row=0, column=3, padx=(8, 0))
    ttk.Label(out_card, text="Report file", style="Field.TLabel").grid(row=1, column=0, sticky="w", padx=(0, 12),
                                                                       pady=3)
    report_entry = ttk.Entry(out_card, textvariable=report_var)
    report_entry.grid(row=1, column=1, columnspan=2, sticky="ew", pady=3)
    report_entry.state(["disabled"])                   # worked out from the output name, not an input
    help_button(out_card, "Report file", "report").grid(row=1, column=3, padx=(8, 0))

    # ---- junctions (worked out, greyed)
    junc_card = ui.Card(left, "Junctions", help_command=lambda: show_help("Junctions", "junctions"),
                        subtitle="worked out from the segments")
    junc_card.grid(row=2, column=0, sticky="ew", pady=(0, 10), padx=(0, 6))
    junc_card.columnconfigure(0, weight=1)
    junc_frame = ttk.Frame(junc_card)
    junc_frame.grid(row=0, column=0, sticky="ew")

    # ---- options
    opt_card = ui.Card(left, "Options")
    opt_card.grid(row=3, column=0, sticky="ew", pady=(0, 10), padx=(0, 6))
    opt_card.columnconfigure(1, weight=1)
    ttk.Label(opt_card, text="Rise (A)", style="Field.TLabel").grid(row=0, column=0, sticky="w", padx=(0, 12), pady=3)
    rise_row = ttk.Frame(opt_card)
    rise_row.grid(row=0, column=1, sticky="w", pady=3)
    ttk.Entry(rise_row, textvariable=rise_var, width=8).pack(side="left")
    help_button(rise_row, "Rise", "rise").pack(side="left", padx=(8, 0))
    ttk.Label(opt_card, text="Processes", style="Field.TLabel").grid(row=1, column=0, sticky="w", padx=(0, 12), pady=3)
    cores_row = ttk.Frame(opt_card)
    cores_row.grid(row=1, column=1, sticky="w", pady=3)
    ttk.Entry(cores_row, textvariable=cores_var, width=8).pack(side="left")
    help_button(cores_row, "Processes", "processes").pack(side="left", padx=(8, 0))
    ttk.Label(cores_row, text=f"default: all but one ({cpu_default})", style="Hint.TLabel").pack(side="left",
                                                                                                 padx=(10, 0))
    share_row = ttk.Frame(opt_card)
    share_row.grid(row=2, column=0, columnspan=2, sticky="w", pady=(6, 0))
    ttk.Checkbutton(share_row, text="Fit junctions that join the same end pairs together",
                    variable=share_var).pack(side="left")
    help_button(share_row, "Fit equal junctions together", "share").pack(side="left", padx=(8, 0))
    reuse_row = ttk.Frame(opt_card)
    reuse_row.grid(row=3, column=0, columnspan=2, sticky="w", pady=(2, 0))
    ttk.Checkbutton(reuse_row, text="Reuse fits from an earlier report", variable=reuse_var).pack(side="left")
    help_button(reuse_row, "Reuse fits", "reuse").pack(side="left", padx=(8, 0))
    ttk.Label(opt_card, text="Earlier report", style="Field.TLabel").grid(row=4, column=0, sticky="w", padx=(0, 12),
                                                                          pady=(4, 3))
    reuse_path_row = ttk.Frame(opt_card)
    reuse_path_row.grid(row=4, column=1, sticky="ew", pady=(4, 3))
    reuse_path_row.columnconfigure(0, weight=1)
    reuse_entry = ttk.Entry(reuse_path_row, textvariable=reuse_path_var)
    reuse_entry.grid(row=0, column=0, sticky="ew")

    def browse_report():
        cur = reuse_path_var.get().strip()
        start = Path(cur).expanduser().parent if cur else None
        path = filedialog.askopenfilename(parent=root, title="Choose the earlier report",
                                          initialdir=str(start) if start and start.is_dir() else None,
                                          filetypes=[("JSON reports", "*.json"), ("All files", "*.*")])
        if path:
            reuse_path_var.set(path)
            show_end(reuse_entry)

    reuse_browse = ttk.Button(reuse_path_row, text="Browse...", command=browse_report)
    reuse_browse.grid(row=0, column=1, padx=(8, 0))

    # ---- run log
    log_card = ui.Card(body, "Run log", help_command=lambda: show_help("Run log", "log"))
    log_card.grid(row=0, column=1, sticky="nsew")
    log_card.columnconfigure(0, weight=1)
    log_card.rowconfigure(0, weight=1)
    log_view = ui.LogView(log_card, wrap="word", height=20)
    log_view.grid(row=0, column=0, sticky="nsew")
    log_text = log_view.text

    def append(text):
        log_text.insert("end", text)
        log_text.see("end")
        if state["mirror"]:                            # the RE_Helix main window mirrors the tool's output
            try:
                sys.stdout.write(text)
                sys.stdout.flush()
            except Exception:  # noqa: BLE001 - nobody reads it any more
                state["mirror"] = False

    # ---- segment rows
    def update_buttons():
        n = len(rows)
        i = selected.get()
        remove_btn.state(["!disabled"] if n > 2 else ["disabled"])
        up_btn.state(["!disabled"] if i > 0 else ["disabled"])
        down_btn.state(["!disabled"] if i < n - 1 else ["disabled"])

    def browse(i):
        cur = rows[i]["path"].get().strip() or (rows[0]["path"].get().strip() if rows else "")
        start = Path(cur).expanduser().parent if cur else None
        path = filedialog.askopenfilename(parent=root, title=f"Choose segment {i + 1}",
                                          initialdir=str(start) if start and start.is_dir() else None,
                                          filetypes=[("PDB files", "*.pdb"), ("All files", "*.*")])
        if path:
            rows[i]["path"].set(path)
            selected.set(i)
            update_buttons()
            if rows[i].get("entry") is not None:
                show_end(rows[i]["entry"])

    def show_end(entry):
        """Scroll a path entry to its end, so the file name shows."""
        try:
            entry.icursor("end")
            entry.xview_moveto(1.0)
        except tk.TclError:
            pass

    def render_rows():
        for w in rows_frame.winfo_children():
            w.destroy()
        for i, row in enumerate(rows):
            r, top = 2 * i, (0 if i == 0 else 8)
            ttk.Radiobutton(rows_frame, text=str(i + 1), value=i, variable=selected, command=update_buttons).grid(
                row=r, column=0, sticky="w", padx=(0, 8), pady=(top, 0))
            entry = ttk.Entry(rows_frame, textvariable=row["path"])
            entry.grid(row=r, column=1, sticky="ew", pady=(top, 0))
            entry.after_idle(lambda e=entry: show_end(e))
            row["entry"] = entry
            entry.bind("<FocusIn>", lambda _e, i=i: (selected.set(i), update_buttons()), add="+")
            ttk.Checkbutton(rows_frame, text="upside down", variable=row["flip"]).grid(
                row=r, column=2, sticky="w", padx=(10, 0), pady=(top, 0))
            ttk.Button(rows_frame, text="Browse...", command=lambda i=i: browse(i)).grid(
                row=r, column=3, sticky="ew", padx=(8, 0), pady=(top, 0))
            row["label"] = ttk.Label(rows_frame, textvariable=row["info"], justify="left", wraplength=state["wrap"],
                                     style="Error.TLabel" if row["error"] else "Hint.TLabel")
            row["label"].grid(row=r + 1, column=1, columnspan=3, sticky="w", pady=(2, 0))
        update_buttons()
        scroll.refresh()

    def new_row(path="", flip=False):
        row = dict(path=tk.StringVar(value=path), flip=tk.BooleanVar(value=flip), info=tk.StringVar(), error=False)
        row["path"].trace_add("write", schedule_refresh)
        row["flip"].trace_add("write", schedule_refresh)
        return row

    def add_row():
        rows.append(new_row())
        selected.set(len(rows) - 1)
        render_rows()
        schedule_refresh()

    def remove_row():
        if len(rows) <= 2:
            return
        i = min(selected.get(), len(rows) - 1)
        del rows[i]
        selected.set(min(i, len(rows) - 1))
        render_rows()
        schedule_refresh()

    def move(delta):
        i = selected.get()
        k = i + delta
        if not (0 <= i < len(rows) and 0 <= k < len(rows)):
            return
        rows[i], rows[k] = rows[k], rows[i]
        selected.set(k)
        render_rows()
        schedule_refresh()

    ttk.Button(seg_bar, text="Add segment", command=add_row).pack(side="left")
    remove_btn = ttk.Button(seg_bar, text="Remove", command=remove_row)
    remove_btn.pack(side="left", padx=(8, 0))
    up_btn = ttk.Button(seg_bar, text="Move up", command=lambda: move(-1))
    up_btn.pack(side="left", padx=(8, 0))
    down_btn = ttk.Button(seg_bar, text="Move down", command=lambda: move(1))
    down_btn.pack(side="left", padx=(8, 0))

    # ---- what the segments give: info lines, junction rows, default output, report path
    def first_segment():
        return rows[0]["path"].get().strip() if rows else ""

    def resolved_output():
        """The model a Run writes: the Output entry, where a name without a folder goes beside the first segment and
        an empty entry means the default; None while there is no first segment either."""
        text, first = out_var.get().strip(), first_segment()
        if not text:
            return default_output_path(first) if first else None
        path = Path(text).expanduser()
        if not path.is_absolute():
            path = (Path(first).expanduser().resolve().parent if first else Path.cwd()) / path
        return path.resolve()

    def update_report(*_):
        target = resolved_output()
        report_var.set(str(default_report_path(target)) if target else
                       "tmp_file/<output>_heterohelix_junction.json beside the output")
        report_entry.after_idle(lambda: show_end(report_entry))
        auto = str(default_report_path(target)) if target else ""
        if reuse_path_var.get().strip() in ("", state["auto_reuse"]):
            reuse_path_var.set(auto)                   # the earlier report follows the output until edited
            reuse_entry.after_idle(lambda: show_end(reuse_entry))
        state["auto_reuse"] = auto

    def update_reuse_state(*_):
        on = reuse_var.get()
        reuse_entry.state(["!disabled"] if on else ["disabled"])
        reuse_browse.state(["!disabled"] if on else ["disabled"])

    def update_output_default():
        first = first_segment()
        auto = str(default_output_path(first)) if first else ""
        current = out_var.get().strip()
        if current in ("", state["auto_out"]):
            out_var.set(auto)
            out_entry.after_idle(lambda: show_end(out_entry))
        state["auto_out"] = auto
        update_report()

    def refill_output(_event=None):
        if not out_var.get().strip() and state["auto_out"]:
            out_var.set(state["auto_out"])             # an emptied Output field goes back to the default
            show_end(out_entry)

    def refresh():
        state["pending"] = None
        entries = [(r["path"].get(), r["flip"].get()) for r in rows]
        plan = plan_model(entries, separate=not share_var.get(), cache=state["cache"])
        for r, s in zip(rows, plan["segments"]):
            r["error"] = s["error"]
            name = Path(r["path"].get().strip()).name
            r["info"].set(f"{name}: {s['text']}" if s["text"] else "Choose a duplex PDB file.")
            if r.get("label") is not None:
                r["label"].configure(style="Error.TLabel" if s["error"] else "Hint.TLabel")
        for w in junc_frame.winfo_children():
            w.destroy()
        for text, bad in zip(plan["junctions"], plan["junction_errors"]):
            ttk.Label(junc_frame, text=text, style="Error.TLabel" if bad else "Muted.TLabel", justify="left",
                      wraplength=state["wrap"]).pack(anchor="w", pady=1)
        for k, text in enumerate(plan["sharing"]):
            ttk.Label(junc_frame, text=text, style="Hint.TLabel", justify="left",
                      wraplength=state["wrap"]).pack(anchor="w", pady=((6 if k == 0 else 0), 0))
        update_output_default()
        scroll.refresh()
        return plan

    def schedule_refresh(*_):
        if state["pending"] is not None:
            try:
                root.after_cancel(state["pending"])
            except tk.TclError:
                pass
        state["pending"] = root.after(250, refresh)

    def on_width(event):
        wrap = max(240, event.width - 70)
        if abs(wrap - state["wrap"]) < 4:
            return                                     # a Configure for a height change (or a few pixels): no re-wrap
        state["wrap"] = wrap
        for r in rows:
            if r.get("label") is not None:
                r["label"].configure(wraplength=state["wrap"])
        for w in junc_frame.winfo_children():
            w.configure(wraplength=state["wrap"])

    seg_card.bind("<Configure>", on_width, add="+")
    out_var.trace_add("write", update_report)
    out_entry.bind("<FocusOut>", refill_output, add="+")
    share_var.trace_add("write", schedule_refresh)
    reuse_var.trace_add("write", update_reuse_state)

    # ---- run
    def build_command():
        entries = [(r["path"].get().strip(), r["flip"].get()) for r in rows]
        if any(not p for p, _ in entries):
            raise ValueError("Every segment row needs a PDB file; fill it in or remove the row.")
        for p, _ in entries:
            if not Path(p).expanduser().is_file():
                raise ValueError(f"Segment file not found: {p}")
        args = [p + (":flip" if f else "") for p, f in entries]
        target = resolved_output()
        if out_var.get().strip() != str(target):
            out_var.set(str(target))                   # show the path the run writes
            show_end(out_entry)
        state["target"] = str(target)
        args += ["-o", str(target)]
        try:
            r = float(rise_var.get().strip())
        except ValueError:
            raise ValueError(f"Rise must be a number in A; got {rise_var.get()!r}") from None
        if not 2.0 <= r <= 6.0:
            raise ValueError("Rise should be between 2 and 6 A.")
        if r != RISE:
            args += ["--rise", f"{r:g}"]
        try:
            n = int(cores_var.get().strip())
        except ValueError:
            raise ValueError(f"Processes must be a whole number; got {cores_var.get()!r}") from None
        if n < 1:
            raise ValueError("Processes must be 1 or more.")
        if n != cpu_default:
            args += ["--cores", str(n)]
        if not share_var.get():
            args.append("--separate")
        if reuse_var.get():
            own = default_report_path(target)
            text = reuse_path_var.get().strip()
            source = Path(text).expanduser() if text else own
            if not source.is_absolute():
                source = target.parent / source
            source = source.resolve()
            args += ["--reuse"] if source == own.resolve() else ["--reuse-from", str(source)]
        return [sys.executable, str(Path(__file__).resolve()), *args]

    def set_running(on):
        state["running"] = on
        run_btn.state(["disabled"] if on else ["!disabled"])
        stop_btn.state(["!disabled"] if on else ["disabled"])

    def start_run():
        if state["running"]:
            return
        try:
            cmd = build_command()
        except ValueError as exc:
            status_var.set(f"Error: {exc}")
            messagebox.showerror(TOOL_NAME, str(exc), parent=root)
            return
        if log_text.get("1.0", "end-1c").strip():
            append("\n")
        append("Command:\n" + shlex.join(cmd) + "\n\n")
        env = dict(os.environ, PYTHONUNBUFFERED="1")
        try:
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1,
                                    env=env, start_new_session=(os.name == "posix"))
        except OSError as exc:
            status_var.set(f"Error: could not start the build ({exc})")
            append(f"Error: could not start the build ({exc})\n")
            return
        state.update(proc=proc, stopped=False, error_lines=[], error_open=False)
        set_running(True)
        status_var.set("Running: a few minutes per junction; the run log shows the progress.")

        def reader():
            for line in proc.stdout:
                q.put(("line", line))
            proc.stdout.close()
            q.put(("done", proc.wait()))

        threading.Thread(target=reader, daemon=True).start()

    def stop_run():
        proc = state["proc"]
        if not state["running"] or proc is None or proc.poll() is not None:
            return
        state["stopped"] = True
        try:
            if os.name == "posix":
                os.killpg(os.getpgid(proc.pid), signal.SIGTERM)   # the build and its worker processes
            else:
                proc.terminate()
        except (ProcessLookupError, PermissionError, OSError):
            pass

    def failure_text():
        """The error of a failed run for the footer: its first line and the first reason under it."""
        lines = state["error_lines"]
        if not lines:
            return "see the run log."
        text = lines[0][len("Error:"):].strip()
        if len(lines) > 1:
            text = f"{text} {lines[1]}"
        if len(text) > 320:
            text = text[:317].rstrip() + "..."
        return text + (" More in the run log." if len(lines) > 2 else "")

    def finish(code):
        set_running(False)
        if state["stopped"]:
            status_var.set("Stopped; nothing more is written.")
            append("[stopped]\n")
        elif code == 0:
            status_var.set(f"Done: wrote {state['target']} and its report in tmp_file/.")
        else:
            status_var.set(f"Failed (exit code {code}): {failure_text()}")

    def take_line(line):
        """One line of the run's output: into the log, and an error with its indented reasons kept for the footer."""
        append(line)
        if line.startswith("Error:"):
            state["error_lines"], state["error_open"] = [line.strip()], True
        elif state["error_open"] and line.startswith((" ", "\t")) and line.strip():
            state["error_lines"].append(line.strip())
        else:
            state["error_open"] = False

    def poll():
        """Take the run's output from the queue; runs every 100 ms (one chain, also when called directly)."""
        if state.get("poll") is not None:
            root.after_cancel(state["poll"])
        try:
            while True:
                kind, value = q.get_nowait()
                if kind == "line":
                    take_line(value)
                else:
                    finish(value)
        except queue.Empty:
            pass
        state["poll"] = root.after(100, poll)

    def close():
        if state["running"]:
            if not messagebox.askyesno(TOOL_NAME, "A build is running. Stop it and close the window?", parent=root):
                return
            stop_run()
        root.destroy()

    footer.add_button("Close", close)
    stop_btn = footer.add_button("Stop", stop_run)
    run_btn = footer.add_button("Run", start_run, primary=True)
    stop_btn.state(["disabled"])
    root.protocol("WM_DELETE_WINDOW", close)

    # ---- initial rows: the segments given on the command line (the main window passes its input PDB first)
    for text in list(segments or []):
        try:
            path, flip = parse_segment(text)
        except JunctionError:
            path, flip = text, False
        rows.append(new_row(path, flip))
    while len(rows) < 2:
        rows.append(new_row())
    render_rows()
    refresh()
    update_report()
    update_reuse_state()
    # handles for tests and scripts that drive the window
    root._heterohelix = dict(rows=rows, selected=selected, refresh=refresh, build_command=build_command,
                             add_row=add_row, remove_row=remove_row, move=move, out_var=out_var,
                             report_var=report_var, share_var=share_var, reuse_var=reuse_var,
                             reuse_path_var=reuse_path_var, reuse_entry=reuse_entry, rise_var=rise_var,
                             cores_var=cores_var, status_var=status_var, junc_frame=junc_frame, run_btn=run_btn,
                             stop_btn=stop_btn, remove_btn=remove_btn, up_btn=up_btn, down_btn=down_btn,
                             log_text=log_text, start_run=start_run, refill_output=refill_output, poll=poll,
                             running=lambda: state["running"])
    state["poll"] = root.after(100, poll)
    root.mainloop()
    return 0


if __name__ == "__main__":
    sys.exit(cli())
