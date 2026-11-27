"""Build Heterohelix Junction (re_helix_lib/build_heterohelix_junction.py).

Fixtures in tests/data, made once from the Junction_Fit masters (junction_fit_azbmost/tmp_file/masters):
  Z8_Z-DNA_CG4.pdb         5'-(CG)4-3' Z-DNA, strand A 1-8 and strand B 17-24 of Z24_master.pdb, renumbered 1-8.
  B8_B-DNA_palindrome.pdb  5'-GCTTAAGC-3' B-DNA: each residue is the B21 master residue of the same name, moved along
                           the master's own screw (one base-pair step) to its place.
  L8_L-DNA_palindrome.pdb  the x -> -x mirror of B8 (L-DNA, D residue names), as BNP_NA makes L-forms.
All three are self-complementary, so both ends of each segment present the same end pair (within 0.04 A).
R24_A-RNA_duplex.pdb (a BNP_NA A-RNA duplex, shared with other tests) and mirror images made here cover A-RNA, L-RNA
and L-Z-DNA.

The end-to-end builds use a small search grid (twists -40 and -20 deg, shift and slide 0, two candidates) through the
tool's hidden --twists/--shifts/--candidates options.  Every test that must not fit replaces the fit with a stand-in
that fails if it is called.
"""
from __future__ import annotations

import contextlib
import copy
import hashlib
import io
import json
import multiprocessing as mp
import os
import shutil
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

from re_helix_lib import build_heterohelix_junction as JF
from re_helix_lib import get_phenix_restraints as gpr


ROOT = Path(__file__).resolve().parents[1]
TOOL = ROOT / "re_helix_lib" / "build_heterohelix_junction.py"
DATA = ROOT / "tests" / "data"
B8 = DATA / "B8_B-DNA_palindrome.pdb"
L8 = DATA / "L8_L-DNA_palindrome.pdb"
Z8 = DATA / "Z8_Z-DNA_CG4.pdb"
R24 = DATA / "R24_A-RNA_duplex.pdb"
SMALL = ["--twists=-40,-20", "--shifts=0", "--candidates=2"]
CORES = str(max(1, min(4, (os.cpu_count() or 2) - 1)))
# steps near the fitted ones on the 21/24-bp masters (same end pairs as these fixtures)
BOWTIE_STEP = [-5.45, -0.32, -2.82, 0.95, 3.24, 3.93]
CONTINUOUS_STEP = [1.18, 1.91, -41.84, 0.77, 0.26, 3.90]
# The small-grid builds of BuildTest, made with jf_fixed2.py (the corrected Junction_Fit the tool was ported from) on
# the same fixtures: every fit's step and total, the group total, and a sha256 of the model's ATOM/HETATM/LINK lines.
REFERENCE = {
    "B8_Z8_B8": dict(
        steps=[[-5.465724555546593, -0.27374780463711057, -2.8143927805759623, 0.9570640405921141, 3.2435047241981936, 3.9311173403284827],
               [5.472553120492157, -0.004209567386879206, -2.8276609764306855, -0.8154243390875852, 2.8975689074548425, 4.222458477107419]],
        totals=[8.184854542638144, 8.350642965956277], group_totals=[16.53549750859442, 16.53549750859442],
        body_sha256="0dfed37a02735516bd707f53845ce32a9e5f8237895aa376dde2e78ff7387019", body_lines=1004),
    "L8_Z8": dict(
        steps=[[1.4625415412616702, 3.4572792158461714, -34.541817516215595, 1.4773231715447654, -0.4571293770860852, 4.0150906865532]],
        totals=[2.4102082768609048], group_totals=[2.4102082768609048],
        body_sha256="9cd96d09e15409cb3951a471e926e407f4e211e474445dca39677b710b8593fc", body_lines=656),
}


class Refit(Exception):
    """Raised by the stand-in fit: the build reached the fitting."""


def no_fit():
    return mock.patch.object(JF, "fit_model_junctions", side_effect=Refit("fitting started"))


def spec(a, ea, b, eb, link="auto"):
    return dict(seg1=str(a), end1=ea, seg2=str(b), end2=eb, link=link, nwin=JF.NWIN, rise=JF.RISE, gap=list(JF.GAP))


def run_cli(*args, merged=False):
    if merged:
        return subprocess.run([sys.executable, str(TOOL), *map(str, args)], stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, text=True, cwd=ROOT)
    return subprocess.run([sys.executable, str(TOOL), *map(str, args)], capture_output=True, text=True, cwd=ROOT)


def run_main(argv):
    """JF.cli in this process: (exit status, stdout, stderr)."""
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        code = JF.cli([str(a) for a in argv])
    return code, out.getvalue(), err.getvalue()


def read_pdb(path):
    atoms, resname, order = {}, {}, {}
    for line in Path(path).read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            key = (line[21], int(line[22:26]))
            atoms[(key, line[12:16].strip())] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            resname[key] = line[17:20].strip()
            order.setdefault(key[0], [])
            if key not in order[key[0]]:
                order[key[0]].append(key)
    return atoms, resname, order


def read_links(path):
    """LINK records: [((chain, resi), atom, residue name, (chain, resi), atom, residue name)]."""
    out = []
    for line in Path(path).read_text().splitlines():
        if line.startswith("LINK"):
            out.append(((line[21], int(line[22:26])), line[12:16].strip(), line[17:20].strip(),
                        (line[51], int(line[52:56])), line[42:46].strip(), line[47:50].strip()))
    return out


def junction_remarks(path):
    out = []
    for line in Path(path).read_text().splitlines():
        if line.startswith("REMARK 950 RE_SCRIPT JUNCTION"):
            kv = gpr.parse_remark_key_values(line)
            kv["res"] = [(c, int(r)) for c, r, _ in (x.split(":") for x in kv["residues"].split(","))]
            out.append(kv)
    return out


def body_sha256(path):
    body = [l for l in Path(path).read_text().splitlines() if l.startswith(("ATOM", "HETATM", "LINK"))]
    return hashlib.sha256("\n".join(body).encode()).hexdigest(), len(body)


def turned(seg, R, t):
    """A copy of a Segment turned and moved rigidly in memory (no rounding through a file)."""
    out = copy.copy(seg)
    out.atoms = [(a[0], a[1], a[2], a[3], R @ a[4] + t) for a in seg.atoms]
    out.ix = {(a[0], a[1], a[3]): a[4] for a in out.atoms}
    return out


def random_rotation(seed):
    q = np.random.default_rng(seed).normal(size=4)
    w, x, y, z = q / np.linalg.norm(q)
    return np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
                     [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
                     [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)]])


def write_with_pair_turned_over(src, dst, a, b):
    """src with the base pair (a, b) turned 180 deg about its own x axis (towards the major groove)."""
    seg = JF.Segment(src)
    o, F = JF.bp_frame(seg.ix, seg.rn, a, b)
    R = 2 * np.outer(F[0], F[0]) - np.eye(3)
    lines = []
    for line in Path(src).read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")) and (line[21], int(line[22:26])) in (a, b):
            x = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            x = o + R @ (x - o)
            line = line[:30] + f"{x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}" + line[54:]
        lines.append(line)
    Path(dst).write_text("\n".join(lines) + "\n")


def write_mirror(src, dst):
    """The x -> -x mirror image of src (L sugars for D, D for L; residue names kept, as BNP_NA writes L-forms)."""
    lines = []
    for line in Path(src).read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            line = line[:30] + f"{-float(line[30:38]):8.3f}" + line[38:]
        lines.append(line)
    Path(dst).write_text("\n".join(lines) + "\n")


def write_filtered(src, dst, keep):
    """src with only the atom lines for which keep(chain, resi, atom name) is true."""
    lines = [l for l in Path(src).read_text().splitlines()
             if not l.startswith(("ATOM", "HETATM")) or keep(l[21], int(l[22:26]), l[12:16].strip())]
    Path(dst).write_text("\n".join(lines) + "\n")


class ChainSenseTest(unittest.TestCase):
    """The chain sense rule and the checks before fitting."""

    @classmethod
    def setUpClass(cls):
        cls._tmp = tempfile.TemporaryDirectory()
        cls.d = Path(cls._tmp.name).resolve()
        cls.R24_mirror = cls.d / "R24_L-RNA_mirror.pdb"
        cls.Z8_mirror = cls.d / "Z8_L-Z-DNA_mirror.pdb"
        write_mirror(R24, cls.R24_mirror)
        write_mirror(Z8, cls.Z8_mirror)

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def test_chain_sense_of_every_form_and_its_mirror(self):
        # a mirror image flips the sign: L-DNA and L-RNA are -1 like Z-DNA, and L-Z-DNA is +1
        for path, sense, sugars in ((B8, +1, "D sugars"), (Z8, -1, "D sugars"), (L8, -1, "L sugars"),
                                    (R24, +1, "D sugars"), (self.R24_mirror, -1, "L sugars"),
                                    (self.Z8_mirror, +1, "L sugars")):
            seg = JF.Segment(path)
            for which in ("bottom", "top"):
                self.assertEqual(JF.End(seg, which).sense, sense, f"{path.name} {which}")
            info = JF.segment_sense(seg)
            self.assertTrue(info["uniform"], path.name)
            self.assertEqual(info["sugars"], sugars, path.name)
            self.assertLess(info["max_turn"], 10.0, path.name)
            self.assertEqual(info["sense"], {"bottom": sense, "top": sense}, path.name)

    def test_link_kind_follows_the_chain_senses(self):
        cases = ((B8, Z8, "bowtie"), (L8, Z8, "continuous"), (B8, L8, "bowtie"), (B8, B8, "continuous"),
                 (Z8, Z8, "continuous"), (R24, Z8, "bowtie"), (R24, self.R24_mirror, "bowtie"),
                 (self.R24_mirror, Z8, "continuous"), (B8, self.Z8_mirror, "continuous"),
                 (L8, self.Z8_mirror, "bowtie"))
        for a, b, kind in cases:
            J = JF.make_junction(spec(a, "top", b, "bottom"))
            self.assertEqual(J.link, kind, f"{a.name} | {b.name}")
            self.assertEqual(JF.natural_link(J.e1, J.e2), kind)

    def test_forced_link_against_the_chain_senses_is_refused_before_fitting(self):
        with self.assertRaises(JF.JunctionError) as ctx:
            JF.make_junction(spec(B8, "top", Z8, "bottom", link="continuous"))
        self.assertIn("chain sense paradox", str(ctx.exception))
        with tempfile.TemporaryDirectory() as tmp, no_fit() as fit:
            out = Path(tmp) / "BZ.pdb"
            with self.assertRaises(JF.JunctionError) as ctx:
                JF.main([str(B8), str(Z8), str(B8), "-o", str(out), "--link", "bowtie,continuous"])
            fit.assert_not_called()
            self.assertFalse(out.exists())
        msg = str(ctx.exception)
        self.assertIn("refused before fitting", msg)
        self.assertIn("junction 2", msg)
        self.assertIn("take a bowtie junction", msg)
        with self.assertRaises(JF.JunctionError) as ctx:
            JF.main([str(B8), str(Z8), str(B8), "--link", "bowtie,continuous,bowtie"])
        self.assertIn("once, or once per junction (2)", str(ctx.exception))

    def test_command_line_prints_the_check_and_stops_with_an_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "BZ.pdb"
            result = run_cli(B8, Z8, "--link", "continuous", "-o", out)
            merged = run_cli(B8, Z8, "--link", "continuous", "-o", out, merged=True)
            self.assertFalse(out.exists())
        self.assertEqual(result.returncode, 1)
        self.assertIn("Chain sense check (before fitting)", result.stdout)
        self.assertIn("junction 1: segment 1 (B8_B-DNA_palindrome.pdb) top (3' end A8 DC, 5' end B1 DG; D-sugars; "
                      "chain sense +1) | segment 2 (Z8_Z-DNA_CG4.pdb) bottom", result.stdout)
        self.assertIn("opposite -> bowtie", result.stdout)
        self.assertIn("--link continuous refused", result.stdout)
        self.assertNotIn("searches", result.stdout)
        self.assertTrue(result.stderr.strip().startswith("Error: refused before fitting"), result.stderr)
        # with both streams in one log, the error comes after the check it refers to
        self.assertLess(merged.stdout.index("Chain sense check"), merged.stdout.index("Error: refused before fitting"))

    def test_upside_down_segment_is_named_in_the_check(self):
        with no_fit():
            code, out, _ = run_main([L8, f"{Z8}:flip", "-o", self.d / "LZflip.pdb"])
        self.assertEqual(code, 1)                                   # the stand-in fit stopped it
        self.assertIn("| segment 2 (Z8_Z-DNA_CG4.pdb, upside down) top (3' end A8 DG", out)

    def test_segment_with_a_pair_turned_over_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            bad = Path(tmp) / "B8_pair4_turned.pdb"
            write_with_pair_turned_over(B8, bad, ("A", 4), ("B", 5))
            info = JF.segment_sense(JF.Segment(bad))
            self.assertFalse(info["uniform"])
            named = [tuple(b["pairs"]) for b in info["bad"]]
            self.assertEqual(named, [("A3-B6", "A4-B5"), ("A4-B5", "A5-B4")])
            self.assertTrue(all(b["turn"] > 150 for b in info["bad"]))
            with no_fit() as fit:
                with self.assertRaises(JF.JunctionError) as ctx:
                    JF.main([str(bad), str(Z8), "-o", str(Path(tmp) / "out.pdb")])
                fit.assert_not_called()
        msg = str(ctx.exception)
        self.assertIn("not uniform", msg)
        self.assertIn("pairs A3-B6 and A4-B5", msg)
        self.assertIn("pairs A4-B5 and A5-B4", msg)

    def test_one_base_pair_is_too_short_and_two_are_enough(self):
        one, two = self.d / "B1_cut.pdb", self.d / "B2_cut.pdb"
        write_filtered(B8, one, lambda c, r, n: (c, r) in (("A", 1), ("B", 8)))
        write_filtered(B8, two, lambda c, r, n: (c, r) in (("A", 1), ("A", 2), ("B", 7), ("B", 8)))
        info = JF.segment_sense(JF.Segment(one))
        self.assertFalse(info["uniform"])
        self.assertEqual(info["sense"], {})                        # no chain sense is read from one pair
        self.assertIn("1 residue per strand, D sugars; too short", JF.segment_text(info))
        info2 = JF.segment_sense(JF.Segment(two))
        self.assertTrue(info2["uniform"])
        self.assertEqual(info2["sense"], {"bottom": +1, "top": +1})
        with no_fit() as fit:
            with self.assertRaises(JF.JunctionError) as ctx:
                JF.main([str(one), str(Z8), "-o", str(self.d / "B1_Z8.pdb")])
            fit.assert_not_called()
        self.assertIn("segment 1 (B1_cut.pdb) has 1 base pair; at least 2 are needed to measure its chain sense",
                      str(ctx.exception))
        plan = JF.plan_model([(one, False), (Z8, False)])
        self.assertTrue(plan["segments"][0]["error"])
        self.assertEqual(plan["junctions"], ["Junction 1: needs two usable segments"])

    def test_missing_link_atom_is_refused_before_fitting(self):
        # B8's top end pair is A8-B1: B1 is the 5' end there, and the links rebuild its phosphate from its O5'
        no_o5 = self.d / "B8_noO5_B1.pdb"
        write_filtered(B8, no_o5, lambda c, r, n: not (c == "B" and r == 1 and n == "O5'"))
        self.assertTrue(JF.segment_sense(JF.Segment(no_o5))["uniform"])
        for other, kind in ((Z8, "bowtie"), (B8, "continuous")):
            with no_fit() as fit:
                with self.assertRaises(JF.JunctionError) as ctx:
                    JF.main([str(no_o5), str(other), "-o", str(self.d / f"noO5_{kind}.pdb")])
                fit.assert_not_called()
            self.assertIn("refused before fitting", str(ctx.exception), kind)
            self.assertIn("junction 1: residue B1 of B8_noO5_B1.pdb has no atom O5'; the links need it",
                          str(ctx.exception), kind)
            plan = JF.plan_model([(no_o5, False), (other, False)])
            self.assertEqual(plan["junction_errors"], [True], kind)
            self.assertIn("has no atom O5'", plan["junctions"][0], kind)

    def test_plan_lists_junctions_and_shared_fits(self):
        plan = JF.plan_model([(B8, False), (Z8, False), (B8, False)])
        self.assertIn("8 residues per strand, D sugars, chain sense +1 at the bottom and +1 at the top",
                      plan["segments"][0]["text"])
        self.assertEqual(plan["junctions"][0],
                         "Junction 1: segment 1 (B8_B-DNA_palindrome.pdb) top -> segment 2 (Z8_Z-DNA_CG4.pdb) bottom; "
                         "chain senses +1 and -1: opposite -> bowtie (a 3'-3' link through an X33 phosphate and a "
                         "5'-5' link)")
        self.assertIn("segment 2 (Z8_Z-DNA_CG4.pdb) top -> segment 3 (B8_B-DNA_palindrome.pdb) bottom; chain senses "
                      "-1 and +1: opposite", plan["junctions"][1])
        self.assertEqual(plan["junction_errors"], [False, False])
        self.assertEqual(plan["sharing"],
                         ["Junction 2 joins the same end pairs as junction 1, turned 180 deg: fitted together."])
        alone = JF.plan_model([(B8, False), (Z8, False), (B8, False)], separate=True)
        self.assertIn("fitted on its own (--separate)", alone["sharing"][0])
        flipped = JF.plan_model([(L8, False), (Z8, True)])
        self.assertIn("segment 2 (Z8_Z-DNA_CG4.pdb, upside down) top; chain senses -1 and -1: matched -> continuous",
                      flipped["junctions"][0])
        missing = JF.plan_model([(B8, False), (DATA / "no_such_file.pdb", False)])
        self.assertEqual(missing["segments"][1], dict(text="file not found", error=True))
        self.assertEqual(missing["junctions"], ["Junction 1: needs two usable segments"])

    def test_default_paths(self):
        out = JF.default_output_path(B8)
        self.assertEqual(out, B8.resolve().with_name("B8_B-DNA_palindrome_heterohelix.pdb"))
        # the default output's report does not repeat the word; other outputs get <stem>_heterohelix_junction.json
        self.assertEqual(JF.default_report_path(out),
                         out.parent / "tmp_file" / "B8_B-DNA_palindrome_heterohelix_junction.json")
        self.assertEqual(JF.default_report_path(Path("/x/model.pdb")),
                         Path("/x/tmp_file/model_heterohelix_junction.json"))
        with tempfile.TemporaryDirectory() as tmp, no_fit() as fit:
            z = Path(tmp) / "Z8_copy.pdb"
            shutil.copy(Z8, z)
            before = z.read_bytes()
            with self.assertRaises(JF.JunctionError) as ctx:
                JF.main([str(B8), str(z), "-o", str(z)])
            fit.assert_not_called()
            self.assertEqual(z.read_bytes(), before)
        self.assertIn("would overwrite an input segment", str(ctx.exception))

    def test_missing_packages_give_one_plain_line(self):
        code = ("import runpy, sys; sys.modules['numpy'] = None; sys.argv = [sys.argv[1], '--version']; "
                "runpy.run_path(sys.argv[0], run_name='__main__')")
        result = subprocess.run([sys.executable, "-c", code, str(TOOL)], capture_output=True, text=True, cwd=ROOT)
        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        self.assertIn("Error: Build Heterohelix Junction needs numpy, scipy and shapely", result.stderr)
        self.assertNotIn("Traceback", result.stderr)


class ScoreTest(unittest.TestCase):
    """The fit score sees only the geometry: the same placement scores the same read from either end, and with a
    segment turned anywhere in space."""

    def _check(self, a, b, p):
        J = JF.make_junction(spec(a, "top", b, "bottom"))
        Jr = JF.make_junction(JF.reverse_spec(spec(a, "top", b, "bottom")))
        q = JF.backward_step(J, Jr, p)
        s, sr = J.score(p), Jr.score(q)
        self.assertLess(abs(s - sr), 1e-6 * max(1.0, abs(s)), f"forward {s} vs backward {sr}")
        R, t = random_rotation(7), np.array([13.7, -42.1, 25.3])
        e1 = JF.End(turned(J.e1.seg, R, t), J.e1.which)
        Jt = JF.Junction(e1, J.e2)
        st = Jt.score(p)
        self.assertLess(abs(s - st), 1e-6 * max(1.0, abs(s)), f"as read {s} vs turned in space {st}")
        return J, Jr, q

    def test_bowtie_score_is_the_same_read_backwards_and_turned(self):
        rng = np.random.default_rng(3)
        for k in range(3):
            p = list(np.array(BOWTIE_STEP) + (0 if k == 0 else rng.normal(scale=[4, 4, 8, 0.5, 0.5, 0.3])))
            self._check(B8, Z8, p)

    def test_continuous_score_is_the_same_read_backwards_and_turned(self):
        rng = np.random.default_rng(4)
        for k in range(3):
            p = list(np.array(CONTINUOUS_STEP) + (0 if k == 0 else rng.normal(scale=[4, 4, 8, 0.5, 0.5, 0.3])))
            self._check(L8, Z8, p)

    def test_forward_and_backward_steps_round_trip_to_the_same_geometry(self):
        J, Jr, q = self._check(B8, Z8, BOWTIE_STEP)
        back = JF.forward_step(J, Jr, q)
        np.testing.assert_allclose(JF.canonical_step(back), JF.canonical_step(BOWTIE_STEP), atol=1e-9)
        # the same placement: every distance between the two windows is the same either way
        XZ = J.place(BOWTIE_STEP)                      # J: B fixed, Z placed
        XB = Jr.place(q)                               # Jr: Z fixed, B placed
        o3_b = [J.e1.k[k] for k in J.e1.key if k[2] in ("O3'", "C1'", "P")]
        o3_z = [J.e2.k[k] for k in J.e2.key if k[2] in ("O3'", "C1'", "P")]
        d_forward = np.linalg.norm(J.e1.X[o3_b][:, None] - XZ[o3_z][None], axis=-1)
        d_backward = np.linalg.norm(XB[[Jr.e2.k[J.e1.key[i]] for i in o3_b]][:, None]
                                    - Jr.e1.X[[Jr.e1.k[J.e2.key[i]] for i in o3_z]][None], axis=-1)
        np.testing.assert_allclose(d_forward, d_backward, atol=1e-9)

    def test_same_junction_finds_the_reversed_match(self):
        J1 = JF.make_junction(spec(B8, "top", Z8, "bottom"))
        J2 = JF.make_junction(spec(Z8, "top", B8, "bottom"))
        self.assertEqual(JF.same_junction(J1, J2), "reversed")
        self.assertEqual(JF.same_junction(J1, JF.make_junction(spec(B8, "bottom", Z8, "top", ))), "aligned")
        self.assertIsNone(JF.same_junction(J1, JF.make_junction(spec(L8, "top", Z8, "bottom"))))
        self.assertLess(JF.end_match(J1.e1, J2.e2, 9), 0.05)
        self.assertLess(JF.end_match(J1.e2, J2.e1, 9), 0.01)
        # the groups the fit makes: together by default, alone with --separate
        self.assertEqual(JF.fit_groups([J1, J2]), [[0, 1]])
        self.assertEqual(JF.fit_groups([J1, J2], share=False), [[0], [1]])


class WorkerSetupTest(unittest.TestCase):
    def test_a_worker_that_cannot_set_up_ends_the_run_with_an_error(self):
        bad = spec(B8, "top", Z8, "bottom", link="continuous")      # refused: JunctionError, a SystemExit
        ctx = mp.get_context("spawn")
        with ctx.Pool(1, initializer=JF._init_worker, initargs=([bad],)) as pool:
            with self.assertRaises(RuntimeError) as err:
                pool.map_async(JF._search, [(0, [0, 0, -20, 0, 0, 3.6])]).get(timeout=120)
        self.assertIn("could not set up the junctions", str(err.exception))
        self.assertIn("chain sense paradox", str(err.exception))


class BuildTest(unittest.TestCase):
    """End-to-end builds through the command line, on the small grid."""

    @classmethod
    def setUpClass(cls):
        cls._tmp = tempfile.TemporaryDirectory()
        cls.d = Path(cls._tmp.name).resolve()            # the tool writes resolved paths in its log
        cls.bzb = cls.d / "B8_Z8_B8.pdb"
        cls.lz = cls.d / "L8_Z8.pdb"
        cls.bzb_run = run_cli(B8, Z8, B8, "-o", cls.bzb, "--cores", CORES, *SMALL)
        cls.lz_run = run_cli(L8, Z8, "--cores", CORES, *SMALL, "--report", cls.d / "lz.json", "-o", cls.lz)
        for r in (cls.bzb_run, cls.lz_run):
            if r.returncode:
                raise AssertionError(r.stdout + r.stderr)
        cls.bzb_report_path = cls.d / "tmp_file" / "B8_Z8_B8_heterohelix_junction.json"
        cls.bzb_report = json.loads(cls.bzb_report_path.read_text())
        cls.lz_report = json.loads((cls.d / "lz.json").read_text())

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def _chain_steps_bonded(self, path):
        """Every step along every written chain is bonded (1.5-1.7 A), whatever kind of link it is."""
        atoms, resname, order = read_pdb(path)
        lengths = []
        for chain, keys in order.items():
            for u, v in zip(keys, keys[1:]):
                pairs = [((u, "O3'"), (v, "P")), ((u, "P"), (v, "O3'")), ((u, "P"), (v, "O5'")),
                         ((u, "O5'"), (v, "P"))]
                d = min(np.linalg.norm(atoms[a] - atoms[b]) for a, b in pairs if a in atoms and b in atoms)
                lengths.append(d)
        return lengths, resname

    def test_fits_and_model_match_the_fixed_fitter(self):
        # pinned: the tool must give the numbers of the corrected fitter (scratchpad jf_fixed3.py, 2026-11-26)
        for name, report, model in (("B8_Z8_B8", self.bzb_report, self.bzb), ("L8_Z8", self.lz_report, self.lz)):
            ref = REFERENCE[name]
            fits = report["fits"]
            self.assertEqual(len(fits), len(ref["steps"]), name)
            np.testing.assert_allclose([f["step"] for f in fits], ref["steps"], rtol=0, atol=1e-6, err_msg=name)
            np.testing.assert_allclose([f["total"] for f in fits], ref["totals"], rtol=1e-6, err_msg=name)
            np.testing.assert_allclose([f["group"]["group_total"] for f in fits], ref["group_totals"], rtol=1e-6,
                                       err_msg=name)
            self.assertEqual(body_sha256(model), (ref["body_sha256"], ref["body_lines"]), name)

    def test_bowtie_model_links_bonds_and_remarks(self):
        text = self.bzb.read_text()
        self.assertIn("Chain sense check (before fitting)", self.bzb_run.stdout)
        self.assertIn("opposite -> bowtie", self.bzb_run.stdout)
        self.assertIn("junction 2 joins the same two end pairs as junction 1, turned 180 deg", self.bzb_run.stdout)
        links = [l for l in text.splitlines() if l.startswith("LINK")]
        # per junction two X33 P-O3' and one 5'-5' P-O5'; each chain then runs through the Z segment against its
        # file order, so its 7 steps there are written as LINKs too
        self.assertEqual(sum("X33" in l for l in links), 4)
        self.assertEqual(sum("O5'" in l for l in links), 2)
        self.assertEqual(sum("X33" not in l and "O5'" not in l for l in links), 14)
        for l in links:
            self.assertTrue(1.5 <= float(l[73:78]) <= 1.7, l)
        lengths, resname = self._chain_steps_bonded(self.bzb)
        self.assertTrue(all(1.5 <= d <= 1.7 for d in lengths), (min(lengths), max(lengths)))
        self.assertEqual(sum(1 for n in resname.values() if n == "X33"), 2)
        bonds = self.bzb_report["checks"]["bonds"]
        self.assertTrue(1.5 <= bonds["min"] <= bonds["max"] <= 1.7, bonds)
        head = [l for l in text.splitlines() if l.startswith("REMARK 950 RE_SCRIPT")]
        self.assertTrue(head[0].startswith("REMARK 950 RE_SCRIPT SOFTWARE name=build_heterohelix_junction version=1.0"))
        self.assertTrue(head[1].startswith("REMARK 950 RE_SCRIPT COMMAND text=python3 build_heterohelix_junction.py"))
        self.assertEqual(head[2], "REMARK 950 RE_SCRIPT OUTPUT_STAGE name=build_heterohelix_junction")
        junctions = [l for l in head if " JUNCTION " in l]
        self.assertEqual(len(junctions), 4)
        self.assertEqual([gpr.parse_remark_key_values(l)["link"] for l in junctions], ["3to3", "5to5"] * 2)
        self.assertTrue(all(":X33" in gpr.parse_remark_key_values(l)["residues"] for l in junctions[::2]))
        self.assertTrue(text.index("REMARK 900 Built by RE_Helix Build Heterohelix Junction V1.0")
                        < text.index("REMARK 950") < text.index("LINK"))

    def test_junction_residues_are_the_residues_the_links_join(self):
        # bowtie: the 3to3 residues are exactly the two residues bonded to that X33, the 5to5 pair is the O5'-P LINK
        links = read_links(self.bzb)
        x33 = {}
        for r1, a1, n1, r2, a2, _ in links:
            if (a1, n1, a2) == ("P", "X33", "O3'"):
                x33.setdefault(r1, set()).add(r2)
        five = [{r1, r2} for r1, a1, _, r2, a2, _ in links if (a1, a2) == ("O5'", "P")]
        rem = junction_remarks(self.bzb)
        self.assertEqual(len(rem), 4)
        for kv in rem:
            if kv["link"] == "3to3":
                end3_1, x, end3_2 = kv["res"]
                self.assertEqual(x33[x], {end3_1, end3_2}, kv["residues"])
            else:
                self.assertEqual(kv["link"], "5to5")
                self.assertIn(set(kv["res"]), five, kv["residues"])
        self.assertEqual(len(x33), 2)
        self.assertEqual(len(five), 2)
        # continuous: the 3to5 pairs are exactly the O3'-P bonds that cross from one segment to the next
        atoms, _, order = read_pdb(self.lz)
        seg_of = self.lz_report["segment_of"]
        crossing = set()
        for keys in order.values():
            for u, v in zip(keys, keys[1:]):
                if seg_of[f"{u[0]}{u[1]}"].split(":")[0] != seg_of[f"{v[0]}{v[1]}"].split(":")[0]:
                    crossing.add((u, v))
                    self.assertTrue(1.5 <= np.linalg.norm(atoms[(u, "O3'")] - atoms[(v, "P")]) <= 1.7, (u, v))
        listed = {tuple(kv["res"]) for kv in junction_remarks(self.lz)}
        self.assertEqual(len(crossing), 2)
        self.assertEqual(listed, crossing)

    def test_get_phenix_restraints_reads_the_junction_residues(self):
        listed = []
        for line in self.bzb.read_text().splitlines():
            if line.startswith("REMARK 950 RE_SCRIPT JUNCTION"):
                listed += gpr.parse_remark_key_values(line)["residues"].split(",")
        found = gpr.read_junction_residues_from_remarks(self.bzb)
        self.assertEqual([f"{c}:{r}:{n}" for c, r, n in found], list(dict.fromkeys(listed)))
        self.assertEqual(len(found), 10)                         # 2 x (two 3' ends, X33, two 5' ends)
        out_base = self.d / "gpr" / "B8_Z8_B8"
        out_base.parent.mkdir(exist_ok=True)
        result = subprocess.run([sys.executable, str(ROOT / "re_helix_lib" / "get_phenix_restraints.py"),
                                 str(self.bzb), "--output-base", str(out_base)], capture_output=True, text=True,
                                cwd=ROOT)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        params = Path(str(out_base) + "_junctions.params").read_text()
        nucleotides = [x for x in dict.fromkeys(listed) if not x.endswith(":X33")]
        self.assertIn("# sugar_expanded_residues = " + ",".join(nucleotides), params)

    def test_chain_sense_after_building(self):
        after = self.bzb_report["chain_sense"]["after"]
        self.assertTrue(after["matched"])
        self.assertEqual(after["errors"], [])
        self.assertEqual([j["j"] for j in after["junctions"]], [1, 2])
        for j in after["junctions"]:
            self.assertEqual(len(j["steps"]), 2)                 # one step across the junction per chain
            self.assertLess(j["max_angle"], 30.0)
        # B-Z bowtie: every pair's z points along the chain's file order, in every segment
        for c in after["chains"]:
            self.assertGreater(c["within_cos"][0], 0.98)
            self.assertEqual(c["faces"], "along")
            self.assertEqual([s["segment"] for s in c["segments"]], [1, 2, 3] if c["chain"] == "A" else [3, 2, 1])
            for s in c["segments"]:
                self.assertGreater(s["axis_cos"][0], 0.9, s)
        self.assertIn("z points along the chain (file order) in every segment", self.bzb_run.stdout)
        self.assertIn("matched: every base pair faces the same way along its chain", self.bzb_run.stdout)
        # L-Z continuous: against it in every segment
        lz = self.lz_report["chain_sense"]
        self.assertTrue(lz["after"]["matched"])
        self.assertEqual(lz["before"]["junctions"][0]["result"], "matched -> continuous")
        for c in lz["after"]["chains"]:
            self.assertEqual(c["faces"], "against")
            for s in c["segments"]:
                self.assertLess(s["axis_cos"][1], -0.9, s)
        self.assertIn("z points against the chain (file order) in every segment", self.lz_run.stdout)

    def test_a_pair_turned_over_after_building_is_an_error(self):
        # rebuild B8 Z8 B8 from this class's fits (no fitting); the check reads coordinates with Z8's bottom end pair
        # (A1-B8, at junction 1) turned over
        out = self.d / "bad" / "B8_Z8_B8_bad.pdb"
        original = JF.write_model

        def write_and_turn(pieces, juncs, fits, links, path, remarks, extra_head=()):
            paths, num, bonds, xyz = original(pieces, juncs, fits, links, path, remarks, extra_head)
            ka, kb = num[(1, "A", 1)], num[(1, "B", 8)]
            ix = {(k[0][0], k[0][1], k[1]): v for k, v in xyz.items()}
            o, F = JF.bp_frame(ix, {ka: pieces[1].seg.rn[("A", 1)], kb: pieces[1].seg.rn[("B", 8)]}, ka, kb)
            R = 2 * np.outer(F[0], F[0]) - np.eye(3)
            return paths, num, bonds, {k: (o + R @ (v - o) if k[0] in (ka, kb) else v) for k, v in xyz.items()}

        with no_fit() as fit, mock.patch.object(JF, "write_model", write_and_turn):
            code, stdout, stderr = run_main([B8, Z8, B8, "-o", out, "--reuse-from", self.bzb_report_path, *SMALL])
            fit.assert_not_called()
        self.assertEqual(code, 1, stdout + stderr)
        self.assertIn(f"reusing the fits from {self.bzb_report_path}", stdout)
        self.assertIn("  Error: chain sense not matched at chain", stdout)
        self.assertIn("chain sense NOT matched", stdout)
        self.assertTrue(stderr.startswith("Error: chain sense not matched in the written model (kept for inspection)"))
        self.assertTrue(out.exists())                            # kept for inspection
        report = json.loads((out.parent / "tmp_file" / "B8_Z8_B8_bad_heterohelix_junction.json").read_text())
        after = report["chain_sense"]["after"]
        self.assertFalse(after["matched"])
        self.assertGreaterEqual(len(after["errors"]), 2)
        self.assertTrue(any("(junction 1)" in e for e in after["errors"]), after["errors"])
        # and the turned pair faces the other way along its chain: segment 2 is read both ways
        self.assertTrue(any("its pairs face along it in segment(s) 1, 3; both ways in segment(s) 2" in e
                            for e in after["errors"]), after["errors"])
        self.assertEqual({c["faces"] for c in after["chains"]}, {"mixed"})
        # the written model itself is the reused build: the same atoms as the first build
        self.assertEqual(body_sha256(out), body_sha256(self.bzb))

    def test_reuse_takes_only_fits_made_the_same_way(self):
        rep = self.bzb_report_path
        other_grid = ["--twists=-40", "--shifts=0", "--candidates=2"]
        cases = [  # (arguments after the segments, an edit of the earlier report or None, what the log must say)
            ([*SMALL, "--separate"], None, "fitted the junctions that join the same end pairs together, and this run "
                                           "fits each one on its own (--separate)"),
            (SMALL, lambda r: r["fitting"].update(groups=[[0], [1]], share=False),     # as a --separate run writes it
             "fitted each junction on its own (--separate), and this run fits the junctions that join the same end "
             "pairs together"),
            (other_grid, None, "was fitted on another search grid"),
            ([*SMALL, "--rise", "3.6"], None, "holds the fits of other junctions"),
            (SMALL, lambda r: r.pop("tool"), "was not written by Build Heterohelix Junction"),
            (SMALL, lambda r: r.pop("fitting"), "was written by another version of Build Heterohelix Junction"),
            (SMALL, lambda r: r["fitting"].update(fitter="Build Heterohelix Junction V0.9"), "another version"),
            (SMALL, lambda r: r["fitting"]["files"].update({str(B8.resolve()): "0" * 32}),
             "segment file B8_B-DNA_palindrome.pdb changed after"),
        ]
        for k, (args, edit, words) in enumerate(cases):
            source = rep
            if edit is not None:
                r = json.loads(rep.read_text())
                edit(r)
                source = self.d / f"edited_{k}.json"
                source.write_text(json.dumps(r))
            with no_fit() as fit:
                code, stdout, _ = run_main([B8, Z8, B8, *args, "-o", self.d / f"reuse_{k}.pdb", "--reuse-from", source])
                fit.assert_called_once()                                 # it went on to fit
            self.assertEqual(code, 1, (args, stdout))                     # stopped by the stand-in fit
            self.assertIn("not reusing earlier fits: ", stdout, args)
            self.assertIn(words, stdout, args)
        # --reuse reads the output's own report, which a new output name does not have yet
        with no_fit() as fit:
            code, stdout, _ = run_main([B8, Z8, B8, *SMALL, "-o", self.d / "new_name.pdb", "--reuse"])
            fit.assert_called_once()
        expected = self.d / "tmp_file" / "new_name_heterohelix_junction.json"
        self.assertIn(f"not reusing earlier fits: no earlier report at {expected}", stdout)
        # a report of the old fitter (Junction_Fit V1.0 wrote no 'tool' or 'fitting') is never reused
        old = {k: v for k, v in self.bzb_report.items() if k not in ("tool", "fitting", "grid")}
        old_path = self.d / "old_junction_fit.json"
        old_path.write_text(json.dumps(old))
        with no_fit() as fit:
            code, stdout, _ = run_main([B8, Z8, B8, *SMALL, "-o", self.d / "from_old.pdb", "--reuse-from", old_path])
            fit.assert_called_once()
        self.assertIn("old_junction_fit.json was not written by Build Heterohelix Junction", stdout)

    def test_reuse_of_the_outputs_own_report(self):
        out = self.d / "own" / "BZB_again.pdb"
        own = JF.default_report_path(out)
        own.parent.mkdir(parents=True)
        shutil.copy(self.bzb_report_path, own)
        with no_fit() as fit:
            code, stdout, stderr = run_main([B8, Z8, B8, *SMALL, "-o", out, "--reuse"])
            fit.assert_not_called()
        self.assertEqual(code, 0, stdout + stderr)
        self.assertIn(f"2 distinct junction(s); reusing the fits from {own}", stdout)
        self.assertEqual(body_sha256(out), body_sha256(self.bzb))
        report = json.loads(own.read_text())
        self.assertEqual(report["fitting"], self.bzb_report["fitting"])
        self.assertIn("--reuse", report["command"])

    def test_continuous_model(self):
        text = self.lz.read_text()
        self.assertNotIn("LINK", text)                           # two ordinary 3'-5' links in file order
        lengths, _ = self._chain_steps_bonded(self.lz)
        self.assertTrue(all(1.5 <= d <= 1.7 for d in lengths), (min(lengths), max(lengths)))
        lines = [l for l in text.splitlines() if l.startswith("REMARK 950 RE_SCRIPT JUNCTION")]
        self.assertEqual(len(lines), 2)
        for l in lines:
            kv = gpr.parse_remark_key_values(l)
            self.assertEqual((kv["kind"], kv["link"], len(kv["residues"].split(","))), ("continuous", "3to5", 2))
        self.assertIn("L-sugar residues", text)

    def test_equal_junctions_get_the_same_step(self):
        fits = self.bzb_report["fits"]
        self.assertEqual(len(fits), 2)
        self.assertEqual(fits[0]["group"]["members"], [0, 1])
        self.assertEqual(fits[1]["group"]["how"], "reversed")
        self.assertEqual(self.bzb_report["fitting"]["groups"], [[0, 1]])
        J1, J2 = (JF.make_junction(f["spec"]) for f in fits)
        # junction 2's step, carried into junction 1's frames, is junction 1's step
        np.testing.assert_allclose(JF.canonical_step(JF.forward_step(J1, J2, fits[1]["step"])),
                                   JF.canonical_step(fits[0]["step"]), atol=1e-7)
        # and in the written model the two junctions have the same geometry: the step between the end-pair frames
        # (y on the residues that one link joins) measured across junction 1 and, read from its other end, junction 2
        atoms, resname, _ = read_pdb(self.bzb)
        ix = {(k[0][0], k[0][1], k[1]): v for k, v in atoms.items()}
        seg_of = self.bzb_report["segment_of"]
        label = {v: k for k, v in seg_of.items()}

        def frame(piece, a, b):
            ka, kb = label[f"{piece}:{a}"], label[f"{piece}:{b}"]
            ra, rb = (ka[0], int(ka[1:])), (kb[0], int(kb[1:]))
            return JF.bp_frame(ix, {ra: resname[ra], rb: resname[rb]}, ra, rb)
        # junction 1: B (piece 1) top pair A8-B1 meets Z (piece 2) bottom pair B8-A1; junction 2 read from the
        # top: B (piece 3) bottom pair B8-A1 meets Z (piece 2) top pair A8-B1.  y on the 3'-3' (X33) side.
        s1 = JF.step_between(*frame(1, "A8", "B1"), *frame(2, "B8", "A1"))
        s2 = JF.step_between(*frame(3, "B8", "A1"), *frame(2, "A8", "B1"))
        np.testing.assert_allclose(s1, s2, atol=0.05)


def _tk():
    try:
        import tkinter as tk
        from tkinter import ttk
    except Exception:  # pragma: no cover
        return None, None
    return tk, ttk


class WindowTest(unittest.TestCase):
    def test_window_rows_junctions_and_command(self):
        tk, ttk = _tk()
        if tk is None:
            self.skipTest("no tkinter")
        seen = {}
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        turned_b8 = Path(tmp.name) / "B8_turned.pdb"
        write_with_pair_turned_over(B8, turned_b8, ("A", 4), ("B", 5))

        def walk(w):
            for c in w.winfo_children():
                yield c
                yield from walk(c)

        def drive(root):
            root.withdraw()
            h = root._heterohelix
            rows = h["rows"]
            seen["rows"] = len(rows)
            seen["first"] = rows[0]["path"].get()
            rows[1]["path"].set(str(Z8))
            h["refresh"]()
            seen["info"] = [r["info"].get() for r in rows]
            seen["out"] = h["out_var"].get()
            seen["report"] = h["report_var"].get()
            report_entries = [w for w in walk(root) if isinstance(w, ttk.Entry)
                              and str(w.cget("textvariable")) == str(h["report_var"])]
            seen["report_disabled"] = all(w.instate(["disabled"]) for w in report_entries) and bool(report_entries)
            seen["remove_enabled_2"] = h["remove_btn"].instate(["!disabled"])
            h["add_row"]()
            h["rows"][2]["path"].set(str(B8))
            h["refresh"]()
            seen["junctions"] = [w.cget("text") for w in h["junc_frame"].winfo_children()]
            seen["remove_enabled_3"] = h["remove_btn"].instate(["!disabled"])
            h["share_var"].set(False)
            h["rows"][1]["flip"].set(True)
            h["refresh"]()
            seen["command"] = h["build_command"]()
            # a name without a folder goes beside the first segment; an emptied field goes back to the default
            h["out_var"].set("joined.pdb")
            seen["relative"] = h["build_command"]()
            seen["relative_shown"] = h["out_var"].get()
            h["out_var"].set("")
            seen["empty_report"] = h["report_var"].get()
            h["refill_output"]()
            seen["refilled"] = h["out_var"].get()
            # reuse: the earlier report is an input only when the box is ticked; it starts as the output's own
            seen["reuse_entry_off"] = h["reuse_entry"].instate(["disabled"])
            h["reuse_var"].set(True)
            seen["reuse_entry_on"] = h["reuse_entry"].instate(["!disabled"])
            seen["reuse_path"] = h["reuse_path_var"].get()
            seen["reuse_cmd"] = h["build_command"]()
            h["reuse_path_var"].set(str(Path(tmp.name) / "earlier.json"))
            seen["reuse_from_cmd"] = h["build_command"]()
            h["reuse_var"].set(False)
            # a real run that is refused before fitting: the footer shows the error and its first reason
            while len(h["rows"]) > 2:
                h["remove_row"]()
            h["rows"][0]["path"].set(str(turned_b8))
            h["rows"][1]["path"].set(str(Z8))
            h["rows"][1]["flip"].set(False)
            h["out_var"].set(str(Path(tmp.name) / "refused.pdb"))
            h["refresh"]()
            h["start_run"]()
            t_end = time.time() + 120
            while h["running"]() and time.time() < t_end:
                # what the window's 100-ms timer does.  No root.update(): on macOS, after an earlier Tk root (the
                # Check Clashes tests make one when they are imported), update() on this withdrawn window can go
                # round the scroll area's geometry handling for ever.
                h["poll"]()
                time.sleep(0.05)
            seen["status"] = h["status_var"].get()
            root.destroy()

        with mock.patch.object(tk.Tk, "mainloop", lambda self, n=0: drive(self)), \
                contextlib.redirect_stdout(io.StringIO()):
            if JF.run_gui([str(B8)]) == 1:
                self.skipTest("no display for Tk")
        self.assertEqual(seen["rows"], 2)
        self.assertEqual(seen["first"], str(B8))
        self.assertIn("chain sense +1 at the bottom and +1 at the top", seen["info"][0])
        self.assertIn("chain sense -1 at the bottom and -1 at the top", seen["info"][1])
        self.assertEqual(seen["out"], str(B8.resolve().with_name("B8_B-DNA_palindrome_heterohelix.pdb")))
        self.assertTrue(seen["report"].endswith("tmp_file/B8_B-DNA_palindrome_heterohelix_junction.json"))
        self.assertTrue(seen["report_disabled"])
        self.assertFalse(seen["remove_enabled_2"])
        self.assertTrue(seen["remove_enabled_3"])
        self.assertTrue(seen["junctions"][0].startswith(
            "Junction 1: segment 1 (B8_B-DNA_palindrome.pdb) top -> segment 2 (Z8_Z-DNA_CG4.pdb) bottom; chain senses "
            "+1 and -1: opposite -> bowtie"))
        self.assertIn("Junction 2 joins the same end pairs as junction 1, turned 180 deg: fitted together.",
                      seen["junctions"])
        cmd = seen["command"]
        self.assertEqual(cmd[:2], [sys.executable, str(TOOL)])
        self.assertEqual(cmd[2:5], [str(B8), f"{Z8}:flip", str(B8)])
        self.assertIn("--separate", cmd)
        self.assertEqual(cmd[cmd.index("-o") + 1], seen["out"])
        joined = str(B8.resolve().parent / "joined.pdb")
        self.assertEqual(seen["relative"][seen["relative"].index("-o") + 1], joined)
        self.assertEqual(seen["relative_shown"], joined)
        self.assertTrue(seen["empty_report"].endswith("tmp_file/B8_B-DNA_palindrome_heterohelix_junction.json"))
        self.assertEqual(seen["refilled"], seen["out"])
        self.assertTrue(seen["reuse_entry_off"])
        self.assertTrue(seen["reuse_entry_on"])
        self.assertEqual(seen["reuse_path"], seen["report"])
        self.assertIn("--reuse", seen["reuse_cmd"])
        self.assertNotIn("--reuse-from", seen["reuse_cmd"])
        rf = seen["reuse_from_cmd"]
        self.assertEqual(rf[rf.index("--reuse-from") + 1], str((Path(tmp.name) / "earlier.json").resolve()))
        self.assertNotIn("--reuse", [a for a in rf if a != "--reuse-from"])
        self.assertTrue(seen["status"].startswith("Failed (exit code 1): refused before fitting: segment 1 "
                                                  "(B8_turned.pdb): the chain sense is not uniform: pairs A3-B6 and "
                                                  "A4-B5 turn"), seen["status"])


if __name__ == "__main__":
    unittest.main()
