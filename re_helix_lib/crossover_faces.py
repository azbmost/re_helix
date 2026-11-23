#!/usr/bin/env python3
"""crossover_faces.py V1.0 - keep each strand of a parallel crossover on its own face (RE_Helix --xo_faces).

A strand that crosses from one helix to another comes in on one face of the plane through the two helix axes and
leaves either on the same face (a parallel crossover, Han et al. 2017 Fig. 1B) or on the other face (an antiparallel
crossover, which goes through the plane by design and is left alone).  A reciprocal exchange rewires bonds without
moving atoms, so at a parallel crossover one atom of each crossing strand, the leaving O3' or the entering phosphate,
is left 0.1-2 A on the other side of that plane, where the two crossing strands meet.

This tool traces every crossing strand: "in" is its phosphate one step before the crossing, "out" its phosphate one
step after.  When both lie on the same face, the atoms of the crossing are put on that face (at least min_face,
default 0.5 A, from the plane):
  1. If the leaving O3' is on the other side, the leaving base pair is turned a few degrees about the helix axis, or
     else the leaving sugar about its glycosidic bond (chi); the phosphates next to it are turned (gamma, beta) so they
     stay bonded with the model's own O3'-P length and angles.
  2. If the entering C5' is on the other side, the entering sugar is turned the same way and the next nucleotide's
     phosphate keeps the bond to it.
  3. The entering phosphate (O5', P, OP1, OP2) is turned about C4'-C5' (gamma) and C5'-O5' (beta) onto the face.
Only rigid turns are used, so every bond length and angle inside a nucleotide is kept and bases stay paired; the
crossover link is never made longer than RE_Helix left it (+0.3 A) and is not closed; moved atoms keep 2.5 A from
every atom more than three bonds away.  A crossover is placed whole or left exactly as built (when a strand cannot
reach its face, or the junction would clash), and every crossing strand is reported.  Leaving-O3' crossings usually
need a local twist change over several base pairs, which this tool does not make.  Minimize afterwards
(get_phenix_restraints, Phenix) when bond geometry matters.

Standalone use (RE_Helix runs it with --xo_faces):
    python3 crossover_faces.py model_aligned_rex.pdb --reference model_aligned.pdb --helix AB --helix CD
writes model_aligned_rex_faces.pdb (or -o) and prints one line per crossing strand.
"""
from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

TOOL_NAME = "Crossover Faces"
TOOL_VERSION = "V1.0"
PHOSPHATE = {"P", "OP1", "OP2", "OP3", "O1P", "O2P", "O3P"}
REMARK = "REMARK 950 RE_SCRIPT XO_FACES"
GRID = np.radians(np.arange(-180.0, 180.0, 2.0))
CHI_STEPS = sorted((d for d in range(-60, 61, 5) if d), key=lambda d: (abs(d), d))
CLASH = 2.5          # A, heavy atoms more than three bonds apart
LINK_SLACK = 0.3     # A the crossover link may grow


# ----------------------------------------------------------------------------- the model
def _name(line: str) -> str:
    return line[12:16].strip().replace("*", "'")


def _is_h(name: str, element: str) -> bool:
    return element.upper() in ("H", "D") or (not element and name[:1] in ("H", "D"))


@dataclass
class Residue:
    key: Tuple[str, int, str]                 # chain, resSeq, iCode
    resname: str
    atoms: Dict[str, int] = field(default_factory=dict)   # atom name -> index into Model.xyz
    order: int = 0


class Model:
    """ATOM/HETATM coordinates of a PDB file, kept with the text so only coordinates change on writing."""

    def __init__(self, lines: Sequence[str]):
        self.lines = list(lines)
        self.line_of: List[int] = []
        self.names: List[str] = []
        self.heavy: List[bool] = []
        self.res_of_atom: List[int] = []
        xyz = []
        self.residues: List[Residue] = []
        last = None
        for i, l in enumerate(self.lines):
            if not l.startswith(("ATOM", "HETATM")):
                continue
            key = (l[21], int(l[22:26]), l[26])
            if key != last:
                self.residues.append(Residue(key, l[17:20].strip(), order=len(self.residues)))
                last = key
            res = self.residues[-1]
            nm = _name(l)
            res.atoms.setdefault(nm, len(xyz))
            self.line_of.append(i)
            self.names.append(nm)
            self.heavy.append(not _is_h(nm, l[76:78].strip()))
            self.res_of_atom.append(len(self.residues) - 1)
            xyz.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
        self.xyz = np.array(xyz, float)
        self.heavy = np.array(self.heavy)
        self.res_of_atom = np.array(self.res_of_atom)

    def get(self, res: Residue, name: str) -> Optional[np.ndarray]:
        i = res.atoms.get(name)
        return None if i is None else self.xyz[i]

    def text(self) -> List[str]:
        out = list(self.lines)
        for k, i in enumerate(self.line_of):
            l = out[i]
            x, y, z = self.xyz[k]
            out[i] = f"{l[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{l[54:]}"
        return out


def connected_next(model: Model, max_gap: float = 8.0) -> Dict[int, int]:
    """Residue index -> next residue index along a strand: consecutive residues of one chain whose O3'-P is shorter
    than max_gap (exchange links are 4-7 A)."""
    nxt = {}
    for a, b in zip(model.residues, model.residues[1:]):
        if a.key[0] != b.key[0]:
            continue
        o3, p = model.get(a, "O3'"), model.get(b, "P")
        if o3 is not None and p is not None and np.linalg.norm(o3 - p) < max_gap:
            nxt[a.order] = b.order
    return nxt


def bond_graph(model: Model, nxt: Dict[int, int]) -> List[set]:
    """Covalent bonds: atom pairs of one residue closer than 1.9 A (1.25 A with H), plus O3'-P along each strand."""
    adj = [set() for _ in range(len(model.xyz))]
    for res in model.residues:
        idx = list(res.atoms.values())
        P = model.xyz[idx]
        d = np.linalg.norm(P[:, None] - P[None], axis=2)
        for a in range(len(idx)):
            for b in range(a + 1, len(idx)):
                lim = 1.9 if (model.heavy[idx[a]] and model.heavy[idx[b]]) else 1.25
                if d[a, b] < lim:
                    adj[idx[a]].add(idx[b])
                    adj[idx[b]].add(idx[a])
    for a, b in nxt.items():
        i, j = model.residues[a].atoms.get("O3'"), model.residues[b].atoms.get("P")
        if i is not None and j is not None:
            adj[i].add(j)
            adj[j].add(i)
    return adj


def within3(adj: List[set], i: int) -> set:
    s1 = adj[i]
    s2 = set().union(*(adj[j] for j in s1)) if s1 else set()
    s3 = set().union(*(adj[j] for j in s2)) if s2 else set()
    return s1 | s2 | s3 | {i}


# ----------------------------------------------------------------------------- geometry helpers
def _unit(v):
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


def rotate(points: np.ndarray, origin: np.ndarray, axis: np.ndarray, angles: np.ndarray) -> np.ndarray:
    """Rodrigues rotation of points (n, 3) about the line (origin, axis) by each angle -> (len(angles), n, 3).
    origin/axis may also carry a leading batch dimension (G, 3) matching points (G, n, 3): -> (G, len(angles), n, 3)."""
    k = _unit(axis)
    if points.ndim == 2:
        v = points - origin
        c, s = np.cos(angles)[:, None, None], np.sin(angles)[:, None, None]
        kv = v @ k
        return origin + v[None] * c + np.cross(k, v)[None] * s + (kv[:, None] * k)[None] * (1 - c)
    v = points - origin[:, None, :]                                       # (G, n, 3)
    kv = np.einsum("gnc,gc->gn", v, k)
    kxv = np.cross(k[:, None, :], v)
    c, s = np.cos(angles)[None, :, None, None], np.sin(angles)[None, :, None, None]
    return (origin[:, None, None, :] + v[:, None] * c + kxv[:, None] * s
            + (kv[..., None] * k[:, None, :])[:, None] * (1 - c))


def _angle(a, b, c):
    v1, v2 = a - b, c - b
    cos = np.sum(v1 * v2, -1) / np.linalg.norm(v1, axis=-1) / np.linalg.norm(v2, axis=-1)
    return np.degrees(np.arccos(np.clip(cos, -1, 1)))


def template_link(model: Model, nxt: Dict[int, int], helix_of: Dict[int, object]) -> Tuple[float, float, float]:
    """The model's own O3'-P length and C3'-O3'-P, O3'-P-O5' angles (medians over intact steps inside a helix)."""
    L, a1, a2 = [], [], []
    for a, b in nxt.items():
        ra, rb = model.residues[a], model.residues[b]
        if helix_of.get(a) is None or helix_of.get(a) != helix_of.get(b):
            continue
        o3, p, c3, o5 = model.get(ra, "O3'"), model.get(rb, "P"), model.get(ra, "C3'"), model.get(rb, "O5'")
        if o3 is None or p is None or c3 is None or o5 is None or np.linalg.norm(o3 - p) > 2.0:
            continue
        L.append(np.linalg.norm(o3 - p))
        a1.append(float(_angle(c3, o3, p)))
        a2.append(float(_angle(o3, p, o5)))
    if not L:
        return 1.607, 119.7, 104.0                                      # Phenix rna3p/bna values
    return float(np.median(L)), float(np.median(a1)), float(np.median(a2))


def helix_axes(model: Model, helix_of: Dict[int, object]) -> Dict[object, Tuple[np.ndarray, np.ndarray]]:
    """Axis of every helix: direction from its P atoms (first principal axis), point from a circle fit across it."""
    pts: Dict[object, list] = {}
    for res in model.residues:
        h = helix_of.get(res.order)
        p = model.get(res, "P")
        if h is not None and p is not None:
            pts.setdefault(h, []).append(p)
    axes = {}
    for h, P in pts.items():
        P = np.array(P)
        if len(P) < 6:
            continue
        c = P.mean(0)
        u = np.linalg.svd(P - c)[2][0]
        e1 = np.cross(u, [0, 0, 1.0]) if abs(u[2]) < 0.9 else np.cross(u, [1.0, 0, 0])
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(u, e1)
        q = np.c_[(P - c) @ e1, (P - c) @ e2]
        sol = np.linalg.lstsq(np.c_[2 * q, np.ones(len(q))], (q ** 2).sum(1), rcond=None)[0]
        axes[h] = (u, c + sol[0] * e1 + sol[1] * e2)
    return axes


# ----------------------------------------------------------------------------- one crossing strand
@dataclass
class Crossing:
    D: int                      # residue leaving its helix
    A: int                      # residue entering the other
    helices: Tuple[object, object]
    mid: np.ndarray
    n: np.ndarray               # normal to the plane of the two axes
    h_in: float
    h_out: float
    kind: str                   # "same face" (parallel), "through" (antiparallel), "unclear"
    sg: int = 0                 # the face: +1 or -1 along n
    link0: float = 0.0
    before: float = 0.0         # lowest face height of the crossing atoms before / after
    after: float = 0.0
    moves: List[str] = field(default_factory=list)
    placed: bool = True

    def h(self, p: np.ndarray) -> np.ndarray:
        return (p - self.mid) @ self.n


def find_crossings(model: Model, nxt: Dict[int, int], helix_of: Dict[int, object], axes) -> List[Crossing]:
    out = []
    for a, b in nxt.items():
        ha, hb = helix_of.get(a), helix_of.get(b)
        if ha is None or hb is None or ha == hb or ha not in axes or hb not in axes:
            continue
        rD, rA = model.residues[a], model.residues[b]
        (uD, cD), (uA, cA) = axes[ha], axes[hb]
        u = _unit(uD + (1 if uD @ uA >= 0 else -1) * uA)
        xD, xA = model.get(rD, "O3'"), model.get(rA, "P")
        pD = cD + ((xD - cD) @ uD) * uD
        pA = cA + ((xA - cA) @ uA) * uA
        ev = pA - pD
        ev = _unit(ev - (ev @ u) * u)
        n = _unit(np.cross(u, ev))
        mid = (pD + pA) / 2
        before = model.get(rD, "P") if model.get(rD, "P") is not None else model.get(rD, "C4'")
        nA = nxt.get(b)
        after = None
        if nA is not None and helix_of.get(nA) == hb:
            after = model.get(model.residues[nA], "P")
        if after is None:
            after = model.get(rA, "C4'")
        c = Crossing(a, b, (ha, hb), mid, n, float((before - mid) @ n), float((after - mid) @ n), "unclear",
                     link0=float(np.linalg.norm(xD - xA)))
        if min(abs(c.h_in), abs(c.h_out)) >= 1.0:
            c.kind = "same face" if c.h_in * c.h_out > 0 else "through"
            c.sg = 1 if c.h_in > 0 else -1
        out.append(c)
    return out


class Placer:
    def __init__(self, model: Model, nxt: Dict[int, int], helix_of: Dict[int, object]):
        self.m = model
        self.nxt = nxt
        self.prv = {b: a for a, b in nxt.items()}
        self.adj = bond_graph(model, nxt)
        self.tmpl = template_link(model, nxt, helix_of)
        self.helix_of = helix_of
        self.axes = helix_axes(model, helix_of)

    # -- atom groups
    def idx(self, r: int, names: Iterable[str]) -> List[int]:
        res = self.m.residues[r]
        return [res.atoms[n] for n in names if n in res.atoms]

    def gamma_group(self, r: int) -> List[int]:
        """Atoms that turn with gamma (about C4'-C5'): O5', the phosphate and the H atoms on C5'/O5'/P."""
        res = self.m.residues[r]
        out = []
        for nm, i in res.atoms.items():
            if nm in PHOSPHATE or nm in ("O5'", "H5'", "H5''", "H5'1", "H5'2", "HO5'"):
                out.append(i)
        return out

    def beta_group(self, r: int) -> List[int]:
        res = self.m.residues[r]
        return [i for nm, i in res.atoms.items() if nm in PHOSPHATE]

    def sugar_group(self, r: int) -> List[int]:
        """Atoms that turn with chi: everything but the base and C1' (on the axis)."""
        res = self.m.residues[r]
        return [i for nm, i in res.atoms.items() if ("'" in nm or nm in PHOSPHATE) and nm != "C1'"]

    def glyco(self, r: int):
        res = self.m.residues[r]
        n = "N9" if "N9" in res.atoms else "N1"
        return self.m.get(res, "C1'"), self.m.get(res, n)

    # -- clashes of moved heavy atoms
    def clash_free(self, moved: List[int], new: np.ndarray, exclude: Iterable[int] = ()):
        """new: (K, len(moved), 3) candidate positions.  Returns the closest distance per candidate to a heavy atom
        more than three bonds away, ignoring `exclude` (atoms that are about to be moved themselves)."""
        m = self.m
        skip = set(moved) | set(exclude)
        center = new.reshape(-1, 3).mean(0)
        near = np.where(np.linalg.norm(m.xyz - center, axis=1) < 16.0)[0]
        near = [j for j in near if m.heavy[j] and j not in skip]
        best = np.full(new.shape[0], 99.0)
        for col, i in enumerate(moved):
            if not m.heavy[i]:
                continue
            excl = within3(self.adj, i)
            js = [j for j in near if j not in excl]
            if not js:
                continue
            d = np.linalg.norm(new[:, col, None, :] - m.xyz[js][None], axis=2).min(1)
            best = np.minimum(best, d)
        return best

    # -- re-close a nucleotide's 5' phosphate to the O3' before it, by gamma/beta
    def reclose(self, r: int, xyz: np.ndarray):
        """xyz: current coordinates (after a sugar turn).  Best gamma/beta placement of residue r's 5' phosphate group
        bonded to the previous residue's O3' with the model's own geometry.  Returns (score, new coords of the group)
        or (inf, None) when r has no previous residue."""
        p = self.prv.get(r)
        if p is None:
            return 0.0, {}
        m = self.m
        rp, rr = m.residues[p], m.residues[r]
        if "O3'" not in rp.atoms or "C3'" not in rp.atoms or "P" not in rr.atoms:
            return 0.0, {}
        O3, C3 = xyz[rp.atoms["O3'"]], xyz[rp.atoms["C3'"]]
        g_idx, b_idx = self.gamma_group(r), self.beta_group(r)
        C4, C5 = xyz[rr.atoms["C4'"]], xyz[rr.atoms["C5'"]]
        G = rotate(xyz[g_idx], C5, C5 - C4, GRID)                          # (G, ng, 3)
        o5 = g_idx.index(rr.atoms["O5'"])
        bsel = [g_idx.index(i) for i in b_idx]
        B = rotate(G[:, bsel, :], G[:, o5, :], G[:, o5, :] - C5[None], GRID)   # (G, B, nb, 3)
        pc = bsel.index(g_idx.index(rr.atoms["P"]))
        Pp = B[:, :, pc, :]
        L, a1, a2 = self.tmpl
        eL = np.linalg.norm(Pp - O3, axis=-1) - L
        e1 = _angle(np.broadcast_to(C3, Pp.shape), np.broadcast_to(O3, Pp.shape), Pp) - a1
        e2 = _angle(np.broadcast_to(O3, Pp.shape), Pp, np.broadcast_to(G[:, None, o5, :], Pp.shape)) - a2
        score = (eL / 0.02) ** 2 + (e1 / 2.0) ** 2 + (e2 / 2.0) ** 2
        gi, bi = np.unravel_index(np.argmin(score), score.shape)
        out = {i: G[gi, t] for t, i in enumerate(g_idx)}
        for t, i in enumerate(b_idx):
            out[i] = B[gi, bi, t]
        return float(score[gi, bi]), out

    # -- 1/2: turn a sugar about chi until the given atom is on the face
    def turn_sugar(self, c: Crossing, r: int, atom: str, h_min: float, reclose_res: Optional[int],
                   exclude: Iterable[int] = ()):
        m = self.m
        res = m.residues[r]
        grp = self.sugar_group(r)
        a0, a1 = self.glyco(r)
        if a0 is None or a1 is None or atom not in res.atoms:
            return False
        target = res.atoms[atom]
        for d in CHI_STEPS:
            new = rotate(m.xyz[grp], a0, a1 - a0, np.radians([d]))[0]
            xyz = m.xyz.copy()
            xyz[grp] = new
            if c.h(xyz[target]) * c.sg < h_min:
                continue
            moved = dict(zip(grp, new))
            if reclose_res is not None:                     # the next nucleotide's phosphate follows this O3'
                sc, grp2 = self.reclose(reclose_res, xyz)
            else:                                           # this nucleotide's own phosphate follows the O3' before it
                sc, grp2 = self.reclose(r, xyz)
            if sc > 60:
                continue
            for i, p in grp2.items():
                xyz[i] = p
                moved[i] = p
            ids = list(moved)
            if self.clash_free(ids, np.array([[moved[i] for i in ids]]), exclude)[0] < CLASH:
                continue
            m.xyz[ids] = np.array([moved[i] for i in ids])
            c.moves.append(f"chi {res.key[0]}{res.key[1]} {d:+d}")
            return True
        return False

    # -- 1b: turn the leaving nucleotide's base pair about the helix axis (a local twist) until its O3' is on the face
    def partner(self, r: int) -> Optional[int]:
        """The nucleotide paired with r: another chain of the same helix, closest base atoms (< 3.5 A)."""
        m = self.m
        res = m.residues[r]
        base = [i for nm, i in res.atoms.items() if "'" not in nm and nm not in PHOSPHATE and m.heavy[i]]
        if not base:
            return None
        h = self.helix_of.get(r)
        best = None
        for o in m.residues:
            if o.key[0] == res.key[0] or self.helix_of.get(o.order) != h:
                continue
            ob = [i for nm, i in o.atoms.items() if "'" not in nm and nm not in PHOSPHATE and m.heavy[i]]
            if not ob:
                continue
            d = float(np.linalg.norm(m.xyz[base][:, None] - m.xyz[ob][None], axis=2).min())
            if best is None or d < best[0]:
                best = (d, o.order)
        return best[1] if best and best[0] < 3.5 else None

    def turn_pair(self, c: Crossing, h_min: float, axes, exclude: Iterable[int] = ()):
        m = self.m
        r = c.D
        rD = m.residues[r]
        q = self.partner(r)
        h = self.helix_of.get(r)
        if q is None or h not in axes or "O3'" not in rD.atoms:
            return False
        u, pt = axes[h]
        grp = list(rD.atoms.values()) + list(m.residues[q].atoms.values())
        o3 = rD.atoms["O3'"]
        for d in sorted((x for x in range(-20, 21, 2) if x), key=lambda x: (abs(x), x)):
            new = rotate(m.xyz[grp], pt, u, np.radians([d]))[0]
            xyz = m.xyz.copy()
            xyz[grp] = new
            if c.h(xyz[o3]) * c.sg < h_min:
                continue
            moved = dict(zip(grp, new))
            ok = True
            for rr in (r, q, self.nxt.get(q)):           # the three backbone bonds around the turned pair
                if rr is None:
                    continue
                sc, g2 = self.reclose(rr, xyz)
                if sc > 60:
                    ok = False
                    break
                for i, p_ in g2.items():
                    xyz[i] = p_
                    moved[i] = p_
            if not ok:
                continue
            ids = list(moved)
            if self.clash_free(ids, np.array([[moved[i] for i in ids]]), exclude)[0] < CLASH:
                continue
            m.xyz[ids] = np.array([moved[i] for i in ids])
            c.moves.append(f"pair {rD.key[0]}{rD.key[1]}:{m.residues[q].key[0]}{m.residues[q].key[1]} turned {d:+d} deg")
            return True
        return False

    # -- 3: the entering phosphate onto the face
    def place_phosphate(self, c: Crossing, h_min: float, exclude: Iterable[int] = (), h_cap: float = 1.5):
        m = self.m
        rD, rA = m.residues[c.D], m.residues[c.A]
        g_idx, b_idx = self.gamma_group(c.A), self.beta_group(c.A)
        C4, C5 = m.get(rA, "C4'"), m.get(rA, "C5'")
        O3 = m.get(rD, "O3'")
        G = rotate(m.xyz[g_idx], C5, C5 - C4, GRID)
        o5 = g_idx.index(rA.atoms["O5'"])
        bsel = [g_idx.index(i) for i in b_idx]
        B = rotate(G[:, bsel, :], G[:, o5, :], G[:, o5, :] - C5[None], GRID)    # (G, B, nb, 3)
        Gb = np.broadcast_to(G[:, None], (len(GRID), len(GRID)) + G.shape[1:]).copy()
        Gb[:, :, bsel, :] = B                                                   # every gamma-group atom, per (g, b)
        heavy = [t for t, i in enumerate(g_idx) if m.heavy[i]]
        hmin = (((Gb[:, :, heavy, :] - c.mid) @ c.n) * c.sg).min(-1)
        link = np.linalg.norm(Gb[:, :, g_idx.index(rA.atoms["P"]), :] - O3, axis=-1)
        ok = link <= c.link0 + LINK_SLACK
        disp = np.linalg.norm(Gb - m.xyz[g_idx], axis=-1).max(-1)
        value = np.where(ok, np.minimum(hmin, h_cap), -99.0) - 0.02 * link - 0.01 * disp
        order = np.argsort(value, axis=None)[::-1]
        flat = Gb.reshape(-1, len(g_idx), 3)
        as_built = float((((m.xyz[g_idx][heavy] - c.mid) @ c.n) * c.sg).min())
        for start in range(0, min(len(order), 6000), 1500):
            cand = order[start:start + 1500]
            cand = cand[value.reshape(-1)[cand] > -50]
            if len(cand) == 0:
                break
            clear = self.clash_free(g_idx, flat[cand], exclude)
            good = np.where(clear >= CLASH)[0]
            if len(good):
                k = cand[good[0]]
                new = flat[k]
                if float((((new[heavy] - c.mid) @ c.n) * c.sg).min()) <= as_built + 0.05:
                    return False                                                # no better than as built
                m.xyz[g_idx] = new
                gi, bi = np.unravel_index(k, value.shape)
                c.moves.append(f"phosphate {rA.key[0]}{rA.key[1]} gamma {math.degrees(GRID[gi]):+.0f} "
                               f"beta {math.degrees(GRID[bi]):+.0f}")
                return True
        return False

    def crossing_atoms(self, c: Crossing) -> List[int]:
        return self.idx(c.D, ["O3'"]) + self.idx(c.A, ["P", "OP1", "OP2", "O5'"])

    def lowest(self, c: Crossing) -> float:
        return float(min(c.h(self.m.xyz[i]) * c.sg for i in self.crossing_atoms(c) if self.m.heavy[i]))

    def junction(self, c: Crossing) -> set:
        """Atoms a crossing may move at its junction: the leaving O3' and the entering phosphate group."""
        return set(self.idx(c.D, ["O3'"])) | set(self.gamma_group(c.A))

    def fix_crossover(self, cs: Sequence[Crossing], h_min: float):
        """Put the crossing strands of one crossover (usually two) on their faces.  While they are being placed the
        atoms the others will still move are not counted as clashes; the finished junction is checked at the end."""
        todo = [c for c in cs if c.kind == "same face"]
        for c in todo:
            c.before = self.lowest(c)
        J = set().union(*(self.junction(c) for c in cs)) if cs else set()
        start = self.m.xyz.copy()
        sep0 = self.separation(cs)
        for c in todo:
            rD = self.m.residues[c.D]
            if c.h(self.m.get(rD, "O3'")) * c.sg < h_min:
                if not self.turn_pair(c, h_min, self.axes, exclude=J):
                    self.turn_sugar(c, c.D, "O3'", h_min, reclose_res=None, exclude=J)
        for c in todo:
            rA = self.m.residues[c.A]
            if "C5'" in rA.atoms and c.h(self.m.get(rA, "C5'")) * c.sg < h_min:
                self.turn_sugar(c, c.A, "C5'", h_min, reclose_res=self.nxt.get(c.A), exclude=J)
        for k, c in enumerate(todo):
            if self.lowest(c) < h_min:
                later = set().union(*(set(self.gamma_group(o.A)) for o in todo[k + 1:])) if todo[k + 1:] else set()
                self.place_phosphate(c, h_min, exclude=later)
        for c in todo:
            c.after = self.lowest(c)
            c.placed = c.after >= h_min
        moved = np.where(np.linalg.norm(self.m.xyz - start, axis=1) > 1e-6)[0]
        reason = None
        if len(moved) and (self.closest_contact(moved) < CLASH or self.separation(cs) < min(sep0, CLASH)):
            reason = "the junction would clash"
        elif any(not c.placed for c in todo):
            reason = "not every strand of this crossover reached its face"
        if reason and len(moved):
            self.m.xyz[:] = start                       # a crossover is placed whole or left exactly as built
        if reason:
            for c in todo:
                c.after, c.placed = c.before, False
                c.moves = [f"left as built: {reason}"]

    def path_atoms(self, c: Crossing) -> List[int]:
        """The crossing path of one strand: C3', O3' of the leaving nucleotide; P, OP1, OP2, O5', C5' of the
        entering one."""
        return self.idx(c.D, ["C3'", "O3'"]) + self.idx(c.A, ["P", "OP1", "OP2", "O5'", "C5'"])

    def separation(self, cs: Sequence[Crossing]) -> float:
        """Closest atoms of the two crossing paths of one crossover."""
        if len(cs) != 2:
            return 99.0
        ia, ib = self.path_atoms(cs[0]), self.path_atoms(cs[1])
        return float(np.linalg.norm(self.m.xyz[ia][:, None] - self.m.xyz[ib][None], axis=2).min())

    def closest_contact(self, atoms: Iterable[int]) -> float:
        m = self.m
        heavy = np.where(m.heavy)[0]
        best = 99.0
        for i in atoms:
            if not m.heavy[i]:
                continue
            d = np.linalg.norm(m.xyz[heavy] - m.xyz[i], axis=1)
            excl = within3(self.adj, i)
            for t in np.argsort(d)[:40]:
                j = heavy[t]
                if j != i and j not in excl:
                    best = min(best, float(d[t]))
                    break
        return best


# ----------------------------------------------------------------------------- driver
def pair_crossings(model: Model, cs: List[Crossing]) -> List[Tuple[Crossing, Crossing]]:
    """The two strands of each crossover: in each helix one leaves where the other enters."""
    out, used = [], set()
    for i, a in enumerate(cs):
        if i in used:
            continue
        best = None
        for j, b in enumerate(cs):
            if j == i or j in used or {a.helices[0], a.helices[1]} != {b.helices[0], b.helices[1]}:
                continue
            d = np.linalg.norm(model.get(model.residues[a.D], "O3'") - model.get(model.residues[b.A], "P"))
            if best is None or d < best[0]:
                best = (d, j)
        if best is not None and best[0] < 3.0:
            used |= {i, best[1]}
            out.append((a, cs[best[1]]))
    return out


def place_crossovers_on_faces(lines: Sequence[str], residue_helix, min_face: float = 0.5):
    """lines: the PDB text of an exchanged model.  residue_helix(chain, resSeq, iCode, coords dict) -> helix id or
    None.  Returns (new lines, crossings, pairs, separations)."""
    model = Model(lines)
    helix_of = {}
    for res in model.residues:
        helix_of[res.order] = residue_helix(res.key, {nm: model.xyz[i] for nm, i in res.atoms.items()})
    nxt = connected_next(model)
    axes = helix_axes(model, helix_of)
    cs = find_crossings(model, nxt, helix_of, axes)
    pairs = pair_crossings(model, cs)
    placer = Placer(model, nxt, helix_of)
    start = model.xyz.copy()
    paired = {id(c) for ab in pairs for c in ab}
    for ab in pairs:
        placer.fix_crossover(list(ab), min_face)
    for c in cs:
        if id(c) not in paired:
            placer.fix_crossover([c], min_face)
    moved = np.where(np.linalg.norm(model.xyz - start, axis=1) > 1e-6)[0]
    placer.final_contact = placer.closest_contact(moved) if len(moved) else 99.0
    placer.n_moved = int(len(moved))
    sep = [placer.separation([a, b]) for a, b in pairs]
    new = model.text()
    rem = []
    for c in cs:
        rD, rA = model.residues[c.D], model.residues[c.A]
        tag = f"{rD.key[0]}{rD.key[1]}->{rA.key[0]}{rA.key[1]}"
        if c.kind != "same face":
            rem.append(f"{REMARK} crossing={tag} kind={c.kind.replace(' ', '_')} left_as_built\n")
        else:
            rem.append(f"{REMARK} crossing={tag} kind=same_face lowest_before={c.before:+.2f} "
                       f"lowest_after={c.after:+.2f} {'placed' if c.placed else 'not_placed'}"
                       + (" moves=" + ";".join(m_.replace(' ', '_') for m_ in c.moves) if c.moves else "") + "\n")
    first = next((i for i, l in enumerate(new) if l.startswith(("ATOM", "HETATM"))), 0)
    new = new[:first] + rem + new[first:]
    return new, cs, pairs, sep, model, placer


def summary_lines(model: Model, cs: List[Crossing], pairs, sep, placer=None) -> List[str]:
    out = []
    for c in cs:
        rD, rA = model.residues[c.D], model.residues[c.A]
        tag = f"{rD.key[0]}{rD.key[1]}->{rA.key[0]}{rA.key[1]}"
        if c.kind != "same face":
            out.append(f"  {tag}: {c.kind} (in {c.h_in:+.1f}, out {c.h_out:+.1f} A); left as built")
        else:
            out.append(f"  {tag}: same face ({'+' if c.sg > 0 else '-'}); lowest crossing atom {c.before:+.2f} -> "
                       f"{c.after:+.2f} A; {', '.join(c.moves) if c.moves else 'unchanged'}"
                       + ("" if c.placed else "  [not placed]"))
    if sep:
        out.append(f"  the two crossing paths of each crossover come no closer than {min(sep):.2f}-{max(sep):.2f} A")
    if placer is not None and getattr(placer, "n_moved", 0):
        out.append(f"  {placer.n_moved} atoms moved; closest atom to any of them, more than three bonds away: "
                   f"{placer.final_contact:.2f} A")
    return out


# ----------------------------------------------------------------------------- standalone
KEY_ATOMS = ("C1'", "N9", "N1", "P")


def helix_lookup_from_atoms(atoms: Iterable[Tuple[str, str, float, float, float]], helix_of_chain: Dict[str, object]):
    """Helix of each nucleotide of an exchanged model.  atoms: (chain, atom name, x, y, z) of the same model before
    the exchange, whose chains each lie in one helix; an exchange moves no atom, so a nucleotide is found by its C1'
    (or N9, N1, P) coordinates."""
    where = {nm: {} for nm in KEY_ATOMS}
    for ch, nm, x, y, z in atoms:
        nm = nm.replace("*", "'")
        if nm in where and ch in helix_of_chain:
            where[nm][(round(x, 3), round(y, 3), round(z, 3))] = helix_of_chain[ch]

    def lookup(key, coords):
        for nm in KEY_ATOMS:
            if nm in coords:
                x, y, z = (float(v) for v in coords[nm])
                return where[nm].get((round(x, 3), round(y, 3), round(z, 3)))
        return None
    return lookup


def helix_lookup_from_reference(reference_lines: Sequence[str], helix_of_chain: Dict[str, object]):
    """As helix_lookup_from_atoms, from the PDB text of the pre-exchange model."""
    atoms = [(l[21], _name(l), float(l[30:38]), float(l[38:46]), float(l[46:54]))
             for l in reference_lines if l.startswith(("ATOM", "HETATM"))]
    return helix_lookup_from_atoms(atoms, helix_of_chain)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=f"{TOOL_NAME} {TOOL_VERSION}: keep each strand of a parallel crossover "
                                             "on the face of the helix plane it comes in and leaves on.")
    ap.add_argument("pdb_in", help="exchanged model (for example model_aligned_rex.pdb)")
    ap.add_argument("--reference", help="the same model before the exchange (model_aligned.pdb); without it the "
                                        "input's own chains define the helices")
    ap.add_argument("--helix", action="append", required=True,
                    help="chains of one helix in the reference, for example --helix AB --helix CD")
    ap.add_argument("--min_face", type=float, default=0.5, help="clearance of the crossing atoms (A, default 0.5)")
    ap.add_argument("-o", "--output", help="output PDB (default <input>_faces.pdb)")
    args = ap.parse_args(argv)
    helix_of_chain = {ch: tuple(g) for g in args.helix for ch in g}
    lines = open(args.pdb_in).readlines()
    ref_lines = open(args.reference).readlines() if args.reference else lines
    lookup = helix_lookup_from_reference(ref_lines, helix_of_chain)
    new, cs, pairs, sep, model, placer = place_crossovers_on_faces(lines, lookup, args.min_face)
    out = args.output or (args.pdb_in[:-4] if args.pdb_in.lower().endswith(".pdb") else args.pdb_in) + "_faces.pdb"
    with open(out, "w") as f:
        f.writelines(new)
    print(f"{TOOL_NAME} {TOOL_VERSION}: {len(cs)} crossing strands")
    for l in summary_lines(model, cs, pairs, sep, placer):
        print(l)
    print(f"Wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
