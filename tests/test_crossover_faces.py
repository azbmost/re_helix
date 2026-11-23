from __future__ import annotations

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

import re_helix as re_helix_app
from re_helix_lib import crossover_faces as xof


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "re_helix.py"
DUPLEX = ROOT / "tests" / "data" / "R24_A-RNA_duplex.pdb"     # BNP_NA A-RNA, 24 bp
PARALLEL = "18B 18C d 7B 7C d"         # a parallel double crossover (DPE type): strands come in and leave on one face
ANTIPARALLEL = "8A 16D d 19A 5D d"     # an antiparallel one (DAE type): strands go through the plane by design
COMMON = ["--replicate", "--axis_dist", "22", "--axis_parallel", "y",
          "--axis_range", "A,B", "--axis_range", "C,D", "--cir_shift", "8"]
HELICES = {"A": "AB", "B": "AB", "C": "CD", "D": "CD"}


def _run(ops: str, out: Path, *extra: str) -> subprocess.CompletedProcess:
    cmd = [sys.executable, str(SCRIPT), str(DUPLEX), *ops.split(), *COMMON, *extra, "-o", str(out)]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT)
    if result.returncode:
        raise AssertionError(result.stdout + result.stderr)
    return result


def _atoms(path: Path) -> dict:
    out = {}
    for line in open(path):
        if line.startswith(("ATOM", "HETATM")):
            out[(line[21], int(line[22:26]), line[12:16].strip())] = np.array(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return out


def _xo_remarks(path: Path) -> list:
    return [line for line in open(path) if line.startswith(xof.REMARK)]


class CrossoverFacesTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls._tmp = tempfile.TemporaryDirectory()
        cls.d = Path(cls._tmp.name)
        _run(PARALLEL, cls.d / "par")
        cls.par_log = _run(PARALLEL, cls.d / "par_xo", "--xo_faces").stdout
        _run(ANTIPARALLEL, cls.d / "anti")
        _run(ANTIPARALLEL, cls.d / "anti_xo", "--xo_faces")

    @classmethod
    def tearDownClass(cls) -> None:
        cls._tmp.cleanup()

    def _frames(self, plain: Path, aligned: Path):
        """Crossings and their faces measured on the model as built (independent of the option's own run)."""
        model = xof.Model(open(plain).readlines())
        lookup = xof.helix_lookup_from_reference(open(aligned).readlines(), HELICES)
        helix_of = {r.order: lookup(r.key, {n: model.xyz[i] for n, i in r.atoms.items()}) for r in model.residues}
        nxt = xof.connected_next(model)
        return model, xof.find_crossings(model, nxt, helix_of, xof.helix_axes(model, helix_of))

    def test_option_is_off_by_default(self) -> None:
        self.assertEqual(_xo_remarks(self.d / "par_aligned_rex.pdb"), [])
        self.assertTrue(_xo_remarks(self.d / "par_xo_aligned_rex.pdb"))

    def test_parallel_crossings_end_on_their_faces(self) -> None:
        plain, xo = self.d / "par_aligned_rex.pdb", self.d / "par_xo_aligned_rex.pdb"
        model, crossings = self._frames(plain, self.d / "par_aligned.pdb")
        self.assertEqual(len(crossings), 4)
        self.assertTrue(all(c.kind == "same face" for c in crossings))
        after = _atoms(xo)
        for c in crossings:
            d_res, a_res = model.residues[c.D], model.residues[c.A]
            keys = [(d_res.key[0], d_res.key[1], "O3'")] + [(a_res.key[0], a_res.key[1], n)
                                                             for n in ("P", "OP1", "OP2", "O5'")]
            built = min(float(c.h(model.xyz[model.residues[r.order].atoms[k[2]]])) * c.sg
                        for r, k in zip([d_res] + [a_res] * 4, keys))
            placed = min(float(c.h(after[k])) * c.sg for k in keys)
            self.assertLess(built, 0.0, "as built, an atom of the crossing sits on the other face")
            self.assertGreaterEqual(placed, 0.5, f"crossing {d_res.key}->{a_res.key} is not on its face")
        remarks = _xo_remarks(xo)
        self.assertEqual(len(remarks), 4)
        self.assertTrue(all("kind=same_face" in r and " placed" in r for r in remarks))
        self.assertIn("--xo_faces", self.par_log)

    def test_only_crossing_nucleotides_move_and_keep_their_bonds(self) -> None:
        plain, xo = _atoms(self.d / "par_aligned_rex.pdb"), _atoms(self.d / "par_xo_aligned_rex.pdb")
        self.assertEqual(set(plain), set(xo))
        model, crossings = self._frames(self.d / "par_aligned_rex.pdb", self.d / "par_aligned.pdb")
        near = set()
        for c in crossings:
            for r in (c.D, c.A):
                ch, resi = model.residues[r].key[:2]
                near |= {(ch, resi - 1), (ch, resi), (ch, resi + 1)}
        moved = {k for k in plain if np.linalg.norm(plain[k] - xo[k]) > 1e-3}
        self.assertTrue(moved)
        self.assertTrue({k[:2] for k in moved} <= near, sorted({k[:2] for k in moved} - near))
        by_res = {}
        for k in plain:
            by_res.setdefault(k[:2], []).append(k)
        for res, keys in by_res.items():
            for i, a in enumerate(keys):
                for b in keys[i + 1:]:
                    d0 = np.linalg.norm(plain[a] - plain[b])
                    if d0 < 1.9:
                        self.assertAlmostEqual(np.linalg.norm(xo[a] - xo[b]), d0, delta=0.002, msg=f"{a}-{b}")
        moved_xyz = np.array([xo[k] for k in moved])
        for k, p in xo.items():
            if k in moved or k[2].startswith("H"):
                continue
            close = np.linalg.norm(moved_xyz - p, axis=1).min()
            if close < 2.5:
                # allowed only for atoms bonded to, or within three bonds of, a moved atom
                self.assertLess(np.linalg.norm(plain[k] - p), 1e-3)
                self.assertTrue(any(np.linalg.norm(plain[m_] - plain[k]) < 4.0 for m_ in moved), k)

    def test_antiparallel_crossings_are_left_as_built(self) -> None:
        self.assertEqual(_atoms(self.d / "anti_aligned_rex.pdb").keys(), _atoms(self.d / "anti_xo_aligned_rex.pdb").keys())
        a, b = _atoms(self.d / "anti_aligned_rex.pdb"), _atoms(self.d / "anti_xo_aligned_rex.pdb")
        self.assertTrue(all(np.array_equal(a[k], b[k]) for k in a))
        remarks = _xo_remarks(self.d / "anti_xo_aligned_rex.pdb")
        self.assertEqual(len(remarks), 4)
        self.assertTrue(all("kind=through" in r and "left_as_built" in r for r in remarks))

    def test_gui_command_carries_the_option(self) -> None:
        common = dict(
            script_path=str(SCRIPT), pdb_in="input.pdb", helix_def_text="(AB)", pair_rows=[],
            pair_args_text="A2 B2 d", axis_range_rows=[], axis_move_rows=[], use_user_axis=False,
            user_axis_dir=["0", "0", "1"], user_axis_point=["0", "0", "0"], output_base="output",
            axis_dist="22", axis_parallel="y", fix_chain="", replicate=True, re_only=False, cir_shift="8",
            linker_phosphate_resname="X33",
        )
        self.assertNotIn("--xo_faces", re_helix_app._build_equivalent_cli_command(**common))
        self.assertIn("--xo_faces", re_helix_app._build_equivalent_cli_command(**common, xo_faces=True))


if __name__ == "__main__":
    unittest.main()
