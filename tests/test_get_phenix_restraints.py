from __future__ import annotations

import contextlib
import io
import re
import tempfile
import unittest
from pathlib import Path

from re_helix_lib import get_phenix_restraints as gpr


SUGAR_PHOSPHATE = ("P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'")


def _atom_line(serial: int, name: str, res_name: str, chain_id: str, res_seq: int, record: str = "ATOM") -> str:
    x, y, z = float(serial), float(res_seq), 0.0  # the tool never reads coordinates
    return (
        f"{record:<6s}{serial:5d} {name:<4s} {res_name:>3s} {chain_id}{res_seq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}          {name[0]:>2s}"
    )


def _link_line(atom_1: str, res_1: str, chain_1: str, seq_1: int, atom_2: str, res_2: str, chain_2: str, seq_2: int) -> str:
    return (
        f"LINK        {atom_1:<4s} {res_1:>3s} {chain_1}{seq_1:4d}                "
        f"{atom_2:<4s} {res_2:>3s} {chain_2}{seq_2:4d}     1555   1555  1.61"
    )


def _parse_blocks(text: str, kind: str):
    """Return the key = value fields of every bond/angle/dihedral edit block."""
    blocks = []
    for body in re.findall(rf"^\s*{kind} {{\n(.*?)\n\s*}}", text, flags=re.MULTILINE | re.DOTALL):
        fields = {}
        for line in body.splitlines():
            key, _, value = line.strip().partition(" = ")
            fields[key] = value
        blocks.append(fields)
    return blocks


def _selection(chain_id: str, res_seq: int, name: str) -> str:
    return f"chain {chain_id} and resid {res_seq} and name {name}"


class BackboneLinkRestraintTests(unittest.TestCase):
    def _write_pdb(self, residues, links, hetatm_residues=()) -> Path:
        """residues: (chain, seq, resname[, atom names]); hetatm_residues: (chain, seq, resname, names)."""
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        lines = [_link_line(*link) for link in links]
        serial = 1
        for chain_id, res_seq, res_name, *names in residues:
            for name in (names[0] if names else SUGAR_PHOSPHATE):
                lines.append(_atom_line(serial, name, res_name, chain_id, res_seq))
                serial += 1
        for chain_id, res_seq, res_name, names in hetatm_residues:
            for name in names:
                lines.append(_atom_line(serial, name, res_name, chain_id, res_seq, record="HETATM"))
                serial += 1
        path = Path(tmp.name) / "model.pdb"
        path.write_text("\n".join(lines) + "\nEND\n", encoding="utf-8")
        return path

    def _params(self, pdb: Path, **kwargs):
        return gpr.build_params_text(pdb, gpr.read_link_records(pdb), return_reports=True, **kwargs)

    def _closure_pdb(self, **kwargs) -> Path:
        residues = [("A", seq, "DA") for seq in (1, 2, 3)]
        return self._write_pdb(residues, [("P", "DA", "A", 1, "O3'", "DA", "A", 3)], **kwargs)

    def test_closure_link_gets_phenix_rna3p_restraint_set(self):
        text, reports = self._params(self._closure_pdb())

        bonds = _parse_blocks(text, "bond")
        self.assertEqual(len(bonds), 1)
        self.assertEqual((bonds[0]["distance_ideal"], bonds[0]["sigma"]), ("1.607", "0.015"))

        angles = {
            (a["atom_selection_1"], a["atom_selection_2"], a["atom_selection_3"]): (a["angle_ideal"], a["sigma"])
            for a in _parse_blocks(text, "angle")
        }
        o3, p = _selection("A", 3, "O3'"), _selection("A", 1, "P")
        self.assertEqual(
            angles,
            {
                (_selection("A", 1, "OP1"), p, o3): ("108.000", "3.000"),
                (_selection("A", 1, "OP2"), p, o3): ("108.000", "3.000"),
                (_selection("A", 1, "O5'"), p, o3): ("104.000", "1.500"),
                (_selection("A", 3, "C3'"), o3, p): ("120.200", "1.500"),
            },
        )

        dihedrals = [
            (
                tuple(d[f"atom_selection_{i}"] for i in range(1, 5)),
                (d["angle_ideal"], d["sigma"], d["periodicity"]),
            )
            for d in _parse_blocks(text, "dihedral")
        ]
        self.assertEqual(
            dihedrals,
            [
                ((_selection("A", 3, "C4'"), _selection("A", 3, "C3'"), o3, p), ("-140.000", "35.000", "1")),
                ((_selection("A", 3, "C3'"), o3, p, _selection("A", 1, "O5'")), ("172.000", "30.000", "3")),
                ((o3, p, _selection("A", 1, "O5'"), _selection("A", 1, "C5'")), ("-60.000", "20.000", "3")),
            ],
        )
        self.assertTrue(any("rna3p backbone restraints for 1 standard" in msg for msg in reports))

    def test_legacy_mode_keeps_v1_0_bond_and_p_centered_angles(self):
        text, reports = self._params(self._closure_pdb(), full_backbone_restraints=False)

        bonds = _parse_blocks(text, "bond")
        self.assertEqual([(b["distance_ideal"], b["sigma"]) for b in bonds], [("1.600", "0.020")])
        angles = sorted((a["atom_selection_1"], a["angle_ideal"], a["sigma"]) for a in _parse_blocks(text, "angle"))
        self.assertEqual(
            angles,
            [
                (_selection("A", 1, "O5'"), "103.000", "3.000"),
                (_selection("A", 1, "OP1"), "108.000", "3.000"),
                (_selection("A", 1, "OP2"), "108.000", "3.000"),
            ],
        )
        self.assertEqual(_parse_blocks(text, "dihedral"), [])
        self.assertNotIn("rna3p", text)
        self.assertFalse(any("rna3p" in msg for msg in reports))

    def test_include_builtin_angles_mode_still_adds_backbone_terms(self):
        text, _ = self._params(self._closure_pdb(), skip_phenix_builtin_angles=False)

        angles = _parse_blocks(text, "angle")
        self.assertEqual(len(angles), 7)  # six pairs around P plus C3'-O3'-P
        op1_op2 = [a for a in angles if a["atom_selection_1"].endswith("OP1") and a["atom_selection_3"].endswith("OP2")]
        self.assertEqual([(a["angle_ideal"], a["sigma"]) for a in op1_op2], [("119.000", "3.000")])
        self.assertEqual(len(_parse_blocks(text, "dihedral")), 3)

    def test_inter_chain_link_gets_rna3p_restraints(self):
        residues = [("A", 1, "DA"), ("A", 2, "DA"), ("B", 1, "DT"), ("B", 2, "DT")]
        pdb = self._write_pdb(residues, [("P", "DT", "B", 1, "O3'", "DA", "A", 2)])
        text, _ = self._params(pdb)

        self.assertEqual(len(_parse_blocks(text, "dihedral")), 3)
        c3_o3_p = [a for a in _parse_blocks(text, "angle") if a["angle_ideal"] == "120.200"]
        self.assertEqual(
            [(a["atom_selection_1"], a["atom_selection_3"]) for a in c3_o3_p],
            [(_selection("A", 2, "C3'"), _selection("B", 1, "P"))],
        )

    def test_consecutive_link_is_left_to_phenix(self):
        residues = [("A", seq, "DA") for seq in (1, 2, 3)]
        pdb = self._write_pdb(residues, [("O3'", "DA", "A", 1, "P", "DA", "A", 2)])
        text, reports = self._params(pdb)

        self.assertEqual([(b["distance_ideal"], b["sigma"]) for b in _parse_blocks(text, "bond")], [("1.600", "0.020")])
        self.assertEqual(_parse_blocks(text, "angle"), [])
        self.assertEqual(_parse_blocks(text, "dihedral"), [])
        self.assertTrue(any("1 P-O3' LINK(s) between consecutive residues" in msg for msg in reports))

    def test_nonstandard_residue_keeps_v1_0_handling(self):
        residues = [("A", 1, "DA"), ("A", 2, "DA"), ("A", 3, "MOD")]
        pdb = self._write_pdb(residues, [("P", "DA", "A", 1, "O3'", "MOD", "A", 3)])
        text, reports = self._params(pdb)

        self.assertEqual(text, self._params(pdb, full_backbone_restraints=False)[0])
        self.assertTrue(any("nonstandard residue MOD" in msg for msg in reports))

    def test_missing_sugar_atom_keeps_v1_0_handling(self):
        no_c4 = tuple(name for name in SUGAR_PHOSPHATE if name != "C4'")
        residues = [("A", 1, "DA"), ("A", 2, "DA"), ("A", 3, "DA", no_c4)]
        pdb = self._write_pdb(residues, [("P", "DA", "A", 1, "O3'", "DA", "A", 3)])
        text, reports = self._params(pdb)

        self.assertEqual(text, self._params(pdb, full_backbone_restraints=False)[0])
        self.assertTrue(any("missing C4' of A3" in msg for msg in reports))

    def test_3to3_linker_output_is_unchanged(self):
        residues = [("A", 1, "DA"), ("A", 2, "DA"), ("B", 1, "DT"), ("B", 2, "DT")]
        links = [("P", "X33", "C", 1, "O3'", "DA", "A", 2), ("P", "X33", "C", 1, "O3'", "DT", "B", 2)]
        pdb = self._write_pdb(residues, links, hetatm_residues=[("C", 1, "X33", ("P", "OP1", "OP2"))])
        text, reports = self._params(pdb)

        self.assertEqual(text, self._params(pdb, full_backbone_restraints=False)[0])
        self.assertIn("Internal restraints for standalone X33", text)
        self.assertFalse(any(msg.startswith("NOTE") for msg in reports))

    def test_5to5_link_output_is_unchanged(self):
        residues = [("A", 1, "DA"), ("A", 2, "DA"), ("B", 1, "DT"), ("B", 2, "DT")]
        pdb = self._write_pdb(residues, [("P", "DA", "A", 1, "O5'", "DT", "B", 1)])
        text, _ = self._params(pdb)

        self.assertEqual(text, self._params(pdb, full_backbone_restraints=False)[0])
        self.assertEqual(_parse_blocks(text, "dihedral"), [])

    def test_legacy_atom_names_are_written_as_found(self):
        legacy = tuple({"OP1": "O1P", "OP2": "O2P"}.get(n, n.replace("'", "*")) for n in SUGAR_PHOSPHATE)
        residues = [("A", seq, "DA", legacy) for seq in (1, 2, 3)]
        pdb = self._write_pdb(residues, [("P", "DA", "A", 1, "O3*", "DA", "A", 3)])
        text, _ = self._params(pdb)

        self.assertIn(f"atom_selection_1 = {_selection('A', 1, 'O1P')}", text)
        self.assertIn(f"atom_selection_1 = {_selection('A', 3, 'C4*')}", text)
        self.assertEqual(len(_parse_blocks(text, "dihedral")), 3)

    def test_cli_flag_writes_legacy_file_and_equivalent_command(self):
        pdb = self._closure_pdb()
        base = pdb.with_name("out")
        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout):
            self.assertEqual(gpr.main([str(pdb), "--output-base", str(base), "--legacy-link-restraints"]), 0)

        command = stdout.getvalue().splitlines()[1]
        self.assertIn("--legacy-link-restraints", command)
        self.assertEqual(_parse_blocks(Path(f"{base}_links.params").read_text(encoding="utf-8"), "dihedral"), [])

        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout):
            self.assertEqual(gpr.main([str(pdb), "--output-base", str(base)]), 0)
        self.assertNotIn("--legacy-link-restraints", stdout.getvalue().splitlines()[1])
        self.assertEqual(len(_parse_blocks(Path(f"{base}_links.params").read_text(encoding="utf-8"), "dihedral")), 3)


if __name__ == "__main__":
    unittest.main()
