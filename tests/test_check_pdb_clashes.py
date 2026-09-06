from __future__ import annotations

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from re_helix_lib import check_pdb_clashes as clash_tool


ROOT = Path(__file__).resolve().parents[1]
TOOL = ROOT / "re_helix_lib" / "check_pdb_clashes.py"


def _atom_line(
    serial: int,
    atom_name: str,
    chain_id: str,
    res_seq: int,
    x: float,
    y: float,
    z: float,
    altloc: str = " ",
    res_name: str = "DA",
    element: str = "",
    icode: str = " ",
) -> str:
    el = element or (atom_name[0] if atom_name else "C")
    return (
        f"ATOM  {serial:5d} {atom_name:<4s}{altloc}{res_name:>3s} "
        f"{chain_id}{res_seq:4d}{icode}   "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}          {el:>2s}"
    )


def _link_line(chain1: str, res1: int, atom1: str, chain2: str, res2: int, atom2: str) -> str:
    return (
        f"LINK        {atom1:>4s}  DA {chain1}{res1:4d}"
        f"                {atom2:>4s}  DA {chain2}{res2:4d}     1555   1555  1.20"
    )


class ClashToolTestCase(unittest.TestCase):
    def _write_pdb(self, lines) -> Path:
        handle = tempfile.NamedTemporaryFile("w", suffix=".pdb", delete=False)
        handle.write("\n".join(lines) + "\nEND\n")
        handle.close()
        path = Path(handle.name)
        self.addCleanup(path.unlink)
        return path

    # -- core geometry ----------------------------------------------------

    def test_reports_nonadjacent_clash(self):
        pdb = self._write_pdb(
            [
                _atom_line(1, "C1'", "A", 5, 0.0, 0.0, 0.0),
                _atom_line(2, "C1'", "A", 20, 1.2, 0.0, 0.0),
            ]
        )
        report = clash_tool.check_clashes(pdb)
        self.assertEqual(len(report.clashes), 1)
        self.assertAlmostEqual(report.clashes[0][0], 1.2, places=6)

    def test_pair_at_or_beyond_cutoff_is_not_a_clash(self):
        pdb = self._write_pdb(
            [
                _atom_line(1, "C1'", "A", 5, 0.0, 0.0, 0.0),
                _atom_line(2, "C1'", "A", 20, 1.6, 0.0, 0.0),
            ]
        )
        self.assertEqual(len(clash_tool.check_clashes(pdb).clashes), 0)
        self.assertEqual(len(clash_tool.check_clashes(pdb, cutoff=2.0).clashes), 1)

    def test_sequence_adjacent_pair_is_excluded(self):
        pdb = self._write_pdb(
            [
                _atom_line(1, "P", "A", 1, 0.0, 0.0, 0.0),
                _atom_line(2, "P", "A", 2, 1.0, 0.0, 0.0),
            ]
        )
        self.assertEqual(len(clash_tool.check_clashes(pdb).clashes), 0)
        # A zero window compares every non-identical residue pair.
        self.assertEqual(len(clash_tool.check_clashes(pdb, adjacent_window=0).clashes), 1)

    def test_atoms_in_the_same_residue_are_excluded(self):
        pdb = self._write_pdb(
            [
                _atom_line(1, "OP1", "A", 7, 0.0, 0.0, 0.0),
                _atom_line(2, "OP2", "A", 7, 0.9, 0.0, 0.0),
            ]
        )
        self.assertEqual(len(clash_tool.check_clashes(pdb).clashes), 0)

    def test_hydrogens_are_ignored(self):
        pdb = self._write_pdb(
            [
                _atom_line(1, "H1", "A", 1, 0.0, 0.0, 0.0, element="H"),
                _atom_line(2, "H2", "A", 9, 0.3, 0.0, 0.0, element="H"),
            ]
        )
        report = clash_tool.check_clashes(pdb)
        self.assertEqual(report.heavy_atoms, 0)
        self.assertEqual(len(report.clashes), 0)

    def test_different_chains_are_compared(self):
        pdb = self._write_pdb(
            [
                _atom_line(1, "P", "A", 1, 0.0, 0.0, 0.0),
                _atom_line(2, "P", "B", 1, 1.0, 0.0, 0.0),
            ]
        )
        # Same residue number, different chain: not adjacent, so it is a clash.
        self.assertEqual(len(clash_tool.check_clashes(pdb).clashes), 1)

    # -- LINK-aware exclusion ---------------------------------------------

    def test_link_bonded_residues_are_excluded(self):
        lines = [
            _link_line("A", 5, "C1'", "A", 20, "C1'"),
            _atom_line(1, "C1'", "A", 5, 0.0, 0.0, 0.0),
            _atom_line(2, "C1'", "A", 20, 1.2, 0.0, 0.0),
        ]
        pdb = self._write_pdb(lines)
        report = clash_tool.check_clashes(pdb)
        self.assertEqual(report.link_records, 1)
        self.assertEqual(report.link_excluded, 1)
        self.assertEqual(len(report.clashes), 0)

    def test_link_exclusion_can_be_disabled(self):
        lines = [
            _link_line("A", 5, "C1'", "A", 20, "C1'"),
            _atom_line(1, "C1'", "A", 5, 0.0, 0.0, 0.0),
            _atom_line(2, "C1'", "A", 20, 1.2, 0.0, 0.0),
        ]
        pdb = self._write_pdb(lines)
        report = clash_tool.check_clashes(pdb, use_link_exclusion=False)
        self.assertEqual(report.link_excluded, 0)
        self.assertEqual(len(report.clashes), 1)

    def test_cross_chain_link_is_excluded(self):
        """A reciprocal-exchange junction joins two chains."""
        lines = [
            _link_line("A", 8, "O3'", "B", 3, "P"),
            _atom_line(1, "O3'", "A", 8, 0.0, 0.0, 0.0),
            _atom_line(2, "P", "B", 3, 1.3, 0.0, 0.0),
        ]
        pdb = self._write_pdb(lines)
        self.assertEqual(len(clash_tool.check_clashes(pdb).clashes), 0)

    # -- alternate conformations ------------------------------------------

    def test_alternate_conformations_are_excluded(self):
        pdb = self._write_pdb(
            [
                _atom_line(1, "C1'", "A", 5, 0.0, 0.0, 0.0, altloc="A"),
                _atom_line(2, "C1'", "A", 9, 0.5, 0.0, 0.0, altloc="B"),
            ]
        )
        report = clash_tool.check_clashes(pdb)
        self.assertEqual(report.altloc_excluded, 1)
        self.assertEqual(len(report.clashes), 0)

    def test_same_altloc_still_clashes(self):
        pdb = self._write_pdb(
            [
                _atom_line(1, "C1'", "A", 5, 0.0, 0.0, 0.0, altloc="A"),
                _atom_line(2, "C1'", "A", 9, 0.5, 0.0, 0.0, altloc="A"),
            ]
        )
        self.assertEqual(len(clash_tool.check_clashes(pdb).clashes), 1)

    # -- cyclic closure ----------------------------------------------------

    def test_circular_auto_uses_each_chain_own_bounds(self):
        lines = [
            _atom_line(1, "P", "A", 1, 0.0, 0.0, 0.0),
            _atom_line(2, "P", "A", 30, 1.1, 0.0, 0.0),
            _atom_line(3, "P", "B", 1, 40.0, 0.0, 0.0),
            _atom_line(4, "P", "B", 12, 41.1, 0.0, 0.0),
        ]
        pdb = self._write_pdb(lines)
        self.assertEqual(len(clash_tool.check_clashes(pdb).clashes), 2)
        report = clash_tool.check_clashes(pdb, circular="auto")
        self.assertEqual(len(report.clashes), 0)
        self.assertEqual(report.circular_chains, {"A": (1, 30), "B": (1, 12)})

    def test_circular_fixed_number_applies_to_every_chain(self):
        lines = [
            _atom_line(1, "P", "A", 1, 0.0, 0.0, 0.0),
            _atom_line(2, "P", "A", 63, 1.1, 0.0, 0.0),
        ]
        pdb = self._write_pdb(lines)
        self.assertEqual(len(clash_tool.check_clashes(pdb, circular="63").clashes), 0)
        self.assertEqual(len(clash_tool.check_clashes(pdb, circular="40").clashes), 1)

    def test_circular_defaults_to_off(self):
        report = clash_tool.check_clashes(
            self._write_pdb([_atom_line(1, "P", "A", 1, 0.0, 0.0, 0.0)])
        )
        self.assertEqual(report.circular_mode, "off")
        self.assertEqual(report.circular_chains, {})

    # -- insertion codes ---------------------------------------------------

    def test_insertion_codes_distinguish_residues(self):
        pdb = self._write_pdb(
            [
                _atom_line(1, "C1'", "A", 5, 0.0, 0.0, 0.0, icode="A"),
                _atom_line(2, "C1'", "A", 5, 0.9, 0.0, 0.0, icode="B"),
            ]
        )
        # Same resSeq but different iCode: still within the adjacency window.
        self.assertEqual(len(clash_tool.check_clashes(pdb).clashes), 0)
        report = clash_tool.check_clashes(pdb, adjacent_window=0)
        self.assertEqual(len(report.clashes), 1)

    # -- multi-model files -------------------------------------------------

    def test_model_selection(self):
        lines = [
            "MODEL        1",
            _atom_line(1, "P", "A", 1, 0.0, 0.0, 0.0),
            _atom_line(2, "P", "A", 10, 50.0, 0.0, 0.0),
            "ENDMDL",
            "MODEL        2",
            _atom_line(1, "P", "A", 1, 0.0, 0.0, 0.0),
            _atom_line(2, "P", "A", 10, 1.1, 0.0, 0.0),
            "ENDMDL",
        ]
        pdb = self._write_pdb(lines)
        first = clash_tool.check_clashes(pdb)
        self.assertEqual(first.model_used, 1)
        self.assertEqual(len(first.clashes), 0)
        second = clash_tool.check_clashes(pdb, model=2)
        self.assertEqual(second.model_used, 2)
        self.assertEqual(len(second.clashes), 1)

    def test_skipped_line_count_is_scoped_to_the_checked_model(self):
        """A malformed line in another MODEL must not be reported here."""
        bad = "ATOM      1  P    DA A   1       bad_x   0.000   0.000  1.00  0.00           P"
        lines = [
            "MODEL        1",
            _atom_line(1, "P", "A", 1, 0.0, 0.0, 0.0),
            "ENDMDL",
            "MODEL        2",
            _atom_line(1, "P", "A", 1, 0.0, 0.0, 0.0),
            bad,
            "ENDMDL",
        ]
        pdb = self._write_pdb(lines)
        self.assertEqual(clash_tool.check_clashes(pdb).skipped_atom_lines, 0)
        self.assertEqual(clash_tool.check_clashes(pdb, model=2).skipped_atom_lines, 1)

    def test_unknown_model_raises(self):
        lines = [
            "MODEL        1",
            _atom_line(1, "P", "A", 1, 0.0, 0.0, 0.0),
            "ENDMDL",
        ]
        pdb = self._write_pdb(lines)
        with self.assertRaises(ValueError):
            clash_tool.check_clashes(pdb, model=7)

    def test_model_on_single_model_file_raises(self):
        pdb = self._write_pdb([_atom_line(1, "P", "A", 1, 0.0, 0.0, 0.0)])
        with self.assertRaises(ValueError):
            clash_tool.check_clashes(pdb, model=1)

    # -- option parsing ----------------------------------------------------

    def test_parse_residue_set_accepts_plain_and_chain_qualified(self):
        self.assertEqual(clash_tool.parse_residue_set("none"), set())
        self.assertEqual(clash_tool.parse_residue_set("11,32"), {(None, 11), (None, 32)})
        self.assertEqual(clash_tool.parse_residue_set("A:11, B:32"), {("A", 11), ("B", 32)})
        with self.assertRaises(ValueError):
            clash_tool.parse_residue_set("A:notanumber")

    def test_parse_circular_option(self):
        self.assertEqual(clash_tool.parse_circular_option("off"), ("off", None))
        self.assertEqual(clash_tool.parse_circular_option("0"), ("off", None))
        self.assertEqual(clash_tool.parse_circular_option("auto"), ("auto", None))
        self.assertEqual(clash_tool.parse_circular_option("63"), ("fixed", 63))
        with self.assertRaises(ValueError):
            clash_tool.parse_circular_option("sometimes")

    def test_invalid_cutoff_and_window_raise(self):
        pdb = self._write_pdb([_atom_line(1, "P", "A", 1, 0.0, 0.0, 0.0)])
        with self.assertRaises(ValueError):
            clash_tool.check_clashes(pdb, cutoff=0.0)
        with self.assertRaises(ValueError):
            clash_tool.check_clashes(pdb, adjacent_window=-1)

    def test_composition_labels_selected_residues(self):
        lines = [
            _atom_line(1, "C1'", "A", 11, 0.0, 0.0, 0.0),
            _atom_line(2, "C1'", "A", 30, 1.2, 0.0, 0.0),
        ]
        pdb = self._write_pdb(lines)
        plain = clash_tool.check_clashes(pdb)
        self.assertEqual(dict(plain.composition), {"heavy-heavy": 1})
        labeled = clash_tool.check_clashes(pdb, labeled_residues_text="11")
        self.assertEqual(dict(labeled.composition), {"labeled-other": 1})
        chain_scoped = clash_tool.check_clashes(pdb, labeled_residues_text="B:11")
        self.assertEqual(dict(chain_scoped.composition), {"other-other": 1})

    # -- Chimera / ChimeraX select commands --------------------------------

    def _clash_pdb(self) -> Path:
        return self._write_pdb(
            [
                _atom_line(1, "C1'", "A", 5, 0.0, 0.0, 0.0),
                _atom_line(2, "C1'", "A", 20, 1.2, 0.0, 0.0),
            ]
        )

    def test_chimera_atom_spec_syntax(self):
        atom = clash_tool.PdbAtom(
            serial=3, atom_name="C1'", altloc="", resname="DA", chain="A",
            resseq=20, icode="", coord=(0.0, 0.0, 0.0), element="C",
        )
        # Chimera: #model:residue.chain@atom
        self.assertEqual(clash_tool.chimera_atom_spec(atom), "#0:20.A@C1'")
        self.assertEqual(clash_tool.chimera_atom_spec(atom, base_model=3), "#3:20.A@C1'")
        self.assertEqual(
            clash_tool.chimera_atom_spec(atom, model_used=2), "#0.2:20.A@C1'"
        )

    def test_chimerax_atom_spec_syntax(self):
        atom = clash_tool.PdbAtom(
            serial=3, atom_name="C1'", altloc="", resname="DA", chain="A",
            resseq=20, icode="", coord=(0.0, 0.0, 0.0), element="C",
        )
        # ChimeraX: #model/chain:residue@atom
        self.assertEqual(clash_tool.chimerax_atom_spec(atom), "#1/A:20@C1'")
        self.assertEqual(clash_tool.chimerax_atom_spec(atom, base_model=7), "#7/A:20@C1'")
        self.assertEqual(
            clash_tool.chimerax_atom_spec(atom, model_used=2), "#1.2/A:20@C1'"
        )

    def test_specs_handle_blank_chain(self):
        atom = clash_tool.PdbAtom(
            serial=1, atom_name="P", altloc="", resname="DA", chain=".",
            resseq=12, icode="", coord=(0.0, 0.0, 0.0), element="P",
        )
        # Chimera writes a bare period for a no-ID chain; ChimeraX omits the chain.
        self.assertEqual(clash_tool.chimera_atom_spec(atom), "#0:12.@P")
        self.assertEqual(clash_tool.chimerax_atom_spec(atom), "#1:12@P")

    def test_specs_append_insertion_code(self):
        atom = clash_tool.PdbAtom(
            serial=1, atom_name="P", altloc="", resname="DA", chain="A",
            resseq=52, icode="B", coord=(0.0, 0.0, 0.0), element="P",
        )
        self.assertEqual(clash_tool.chimera_atom_spec(atom), "#0:52B.A@P")
        self.assertEqual(clash_tool.chimerax_atom_spec(atom), "#1/A:52B@P")

    def test_select_all_deduplicates_shared_atoms(self):
        """One atom clashing with two partners must be listed once."""
        pdb = self._write_pdb(
            [
                _atom_line(1, "C1'", "A", 5, 0.0, 0.0, 0.0),
                _atom_line(2, "C1'", "A", 20, 1.2, 0.0, 0.0),
                _atom_line(3, "C1'", "A", 40, -1.2, 0.0, 0.0),
            ]
        )
        report = clash_tool.check_clashes(pdb)
        self.assertEqual(len(report.clashes), 2)
        _chimera, chimerax = clash_tool.select_all_commands(report)
        self.assertEqual(chimerax.count("|"), 2)  # three unique atoms
        self.assertEqual(chimerax.count(":5@"), 1)

    def test_all_select_commands_sit_alone_on_their_line(self):
        """A line selection must pick up the command and nothing else."""
        text = clash_tool.format_report(clash_tool.check_clashes(self._clash_pdb()))
        lines = text.splitlines()
        for key, expected in (
            ("chimera_select_all_clashes", "select #0:5.A@C1'|#0:20.A@C1'"),
            ("chimerax_select_all_clashes", "select #1/A:5@C1'|#1/A:20@C1'"),
        ):
            self.assertIn(key, lines)
            command_line = lines[lines.index(key) + 1]
            self.assertEqual(command_line, expected)
            self.assertNotIn("\t", command_line)

    def test_report_includes_select_commands_by_default(self):
        text = clash_tool.format_report(clash_tool.check_clashes(self._clash_pdb()))
        self.assertIn("chimera_select_all_clashes", text)
        self.assertIn("chimerax_select_all_clashes", text)
        header = next(
            line for line in text.splitlines() if line.startswith("closest_clashes")
        )
        self.assertEqual(
            header, "closest_clashes\tdistance_A\tatom_1\tatom_2\tchimera\tchimerax"
        )
        row = text.splitlines()[-1].split("\t")
        self.assertEqual(len(row), 5)
        self.assertTrue(row[3].startswith("select #0:"))
        self.assertTrue(row[4].startswith("select #1/"))

    def test_select_commands_can_be_suppressed(self):
        text = clash_tool.format_report(
            clash_tool.check_clashes(self._clash_pdb()), include_select_commands=False
        )
        self.assertNotIn("select", text)
        self.assertIn("closest_clashes\tdistance_A\tatom_1\tatom_2", text)

    def test_custom_model_numbers_flow_into_report(self):
        text = clash_tool.format_report(
            clash_tool.check_clashes(self._clash_pdb()), chimera_model=3, chimerax_model=7
        )
        self.assertIn("select #3:5.A@C1'", text)
        self.assertIn("select #7/A:5@C1'", text)

    def test_selected_model_becomes_a_submodel_in_the_commands(self):
        lines = [
            "MODEL        1",
            _atom_line(1, "P", "A", 1, 0.0, 0.0, 0.0),
            _atom_line(2, "P", "A", 10, 50.0, 0.0, 0.0),
            "ENDMDL",
            "MODEL        2",
            _atom_line(1, "P", "A", 1, 0.0, 0.0, 0.0),
            _atom_line(2, "P", "A", 10, 1.1, 0.0, 0.0),
            "ENDMDL",
        ]
        pdb = self._write_pdb(lines)
        text = clash_tool.format_report(clash_tool.check_clashes(pdb, model=2))
        self.assertIn("select #0.2:1.A@P", text)
        self.assertIn("select #1.2/A:1@P", text)

    def test_no_select_lines_when_there_are_no_clashes(self):
        pdb = self._write_pdb(
            [
                _atom_line(1, "C1'", "A", 5, 0.0, 0.0, 0.0),
                _atom_line(2, "C1'", "A", 20, 30.0, 0.0, 0.0),
            ]
        )
        text = clash_tool.format_report(clash_tool.check_clashes(pdb))
        self.assertNotIn("select_all_clashes", text)

    def test_cli_emits_select_commands(self):
        pdb = self._clash_pdb()
        result = subprocess.run(
            [sys.executable, str(TOOL), str(pdb)],
            capture_output=True, text=True, cwd=str(ROOT),
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("chimera_select_all_clashes", result.stdout)
        self.assertIn("chimerax_select_all_clashes", result.stdout)

        suppressed = subprocess.run(
            [sys.executable, str(TOOL), str(pdb), "--no-select-commands"],
            capture_output=True, text=True, cwd=str(ROOT),
        )
        self.assertEqual(suppressed.returncode, 0, suppressed.stderr)
        self.assertNotIn("select_all_clashes", suppressed.stdout)

    # -- backend parity ----------------------------------------------------

    def test_pure_python_matches_kdtree(self):
        lines = []
        serial = 1
        for chain in ("A", "B"):
            for res in range(1, 21):
                for index, name in enumerate(("P", "O5'", "C1'")):
                    lines.append(
                        _atom_line(
                            serial,
                            name,
                            chain,
                            res,
                            0.7 * index + (0.0 if chain == "A" else 1.4),
                            0.9 * res,
                            0.3 * index,
                        )
                    )
                    serial += 1
        pdb = self._write_pdb(lines)

        fast = clash_tool.check_clashes(pdb)
        saved_np, saved_tree = clash_tool.np, clash_tool.cKDTree
        clash_tool.np, clash_tool.cKDTree = None, None
        try:
            slow = clash_tool.check_clashes(pdb)
        finally:
            clash_tool.np, clash_tool.cKDTree = saved_np, saved_tree

        self.assertEqual(slow.method, "pure Python pair scan")
        self.assertEqual(len(fast.clashes), len(slow.clashes))
        self.assertEqual(
            [(round(d, 9), a.serial, b.serial) for d, a, b in fast.clashes],
            [(round(d, 9), a.serial, b.serial) for d, a, b in slow.clashes],
        )

    # -- report and CLI ----------------------------------------------------

    def test_report_contains_expected_keys(self):
        pdb = self._write_pdb(
            [
                _atom_line(1, "C1'", "A", 5, 0.0, 0.0, 0.0),
                _atom_line(2, "C1'", "A", 20, 1.2, 0.0, 0.0),
            ]
        )
        text = clash_tool.format_report(clash_tool.check_clashes(pdb))
        for key in (
            "heavy_atoms",
            "cutoff_A",
            "link_excluded_pairs",
            "altloc_excluded_pairs",
            "heavy_atom_clashes_lt_cutoff",
            "min_clash_distance_A",
            "closest_clashes",
        ):
            self.assertIn(key, text)

    def test_report_says_na_when_clean(self):
        pdb = self._write_pdb(
            [
                _atom_line(1, "C1'", "A", 5, 0.0, 0.0, 0.0),
                _atom_line(2, "C1'", "A", 20, 30.0, 0.0, 0.0),
            ]
        )
        text = clash_tool.format_report(clash_tool.check_clashes(pdb))
        self.assertIn("min_clash_distance_A\tNA", text)
        self.assertNotIn("closest_clashes", text)

    def test_cli_runs_and_writes_report_file(self):
        pdb = self._write_pdb(
            [
                _atom_line(1, "C1'", "A", 5, 0.0, 0.0, 0.0),
                _atom_line(2, "C1'", "A", 20, 1.2, 0.0, 0.0),
            ]
        )
        out_path = pdb.with_name(pdb.stem + "_report.txt")
        self.addCleanup(lambda: out_path.exists() and out_path.unlink())

        result = subprocess.run(
            [sys.executable, str(TOOL), str(pdb), "-o", str(out_path)],
            capture_output=True,
            text=True,
            cwd=str(ROOT),
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("heavy_atom_clashes_lt_cutoff\t1", result.stdout)
        self.assertTrue(out_path.exists())
        self.assertIn("min_clash_distance_A", out_path.read_text())

    def test_cli_reports_error_for_missing_file(self):
        result = subprocess.run(
            [sys.executable, str(TOOL), str(ROOT / "does_not_exist.pdb")],
            capture_output=True,
            text=True,
            cwd=str(ROOT),
        )
        self.assertEqual(result.returncode, 1)
        self.assertIn("Error:", result.stderr)

    def test_build_cli_command_omits_defaults(self):
        command = clash_tool.build_cli_command(
            "check_pdb_clashes.py",
            "model.pdb",
            clash_tool.DEFAULT_CUTOFF,
            clash_tool.DEFAULT_ADJACENT_WINDOW,
            "off",
            "none",
            True,
            clash_tool.DEFAULT_TOP,
            None,
            "",
        )
        self.assertNotIn("--cutoff", command)
        self.assertNotIn("--circular", command)
        self.assertNotIn("--no-link-exclusion", command)

        verbose = clash_tool.build_cli_command(
            "check_pdb_clashes.py",
            "model.pdb",
            2.2,
            0,
            "auto",
            "A:11",
            False,
            5,
            2,
            "out.txt",
        )
        for fragment in (
            "--cutoff 2.2",
            "--adjacent-window 0",
            "--circular auto",
            "--label-residues A:11",
            "--no-link-exclusion",
            "--top 5",
            "--model 2",
            "-o out.txt",
        ):
            self.assertIn(fragment, verbose)

    def test_build_cli_command_round_trips_select_options(self):
        default = clash_tool.build_cli_command(
            "check_pdb_clashes.py", "model.pdb",
            clash_tool.DEFAULT_CUTOFF, clash_tool.DEFAULT_ADJACENT_WINDOW,
            "off", "none", True, clash_tool.DEFAULT_TOP, None, "",
        )
        self.assertNotIn("--no-select-commands", default)
        self.assertNotIn("--chimera-model", default)

        custom = clash_tool.build_cli_command(
            "check_pdb_clashes.py", "model.pdb",
            clash_tool.DEFAULT_CUTOFF, clash_tool.DEFAULT_ADJACENT_WINDOW,
            "off", "none", True, clash_tool.DEFAULT_TOP, None, "",
            include_select_commands=False, chimera_model=3, chimerax_model=7,
        )
        self.assertIn("--no-select-commands", custom)
        self.assertIn("--chimera-model 3", custom)
        self.assertIn("--chimerax-model 7", custom)


def _tk_available() -> bool:
    try:
        import tkinter as tk

        tk.Tk().destroy()
        return True
    except Exception:
        return False


@unittest.skipUnless(_tk_available(), "requires a usable Tk display")
class DoubleClickCopyTestCase(unittest.TestCase):
    """The GUI must select a whole select command on a double-click.

    Tk's own double-click picks a single word, and these commands are full of
    word-breaking punctuation, so the tool binds its own handler.
    """

    def _run_gui_and_double_click(self, line_predicate):
        import tkinter as tk

        captured = {}
        real_tag_bind = tk.Text.tag_bind

        def spy(self, tagName=None, sequence=None, func=None, add=None):
            if tagName == "select_command" and sequence == "<Double-Button-1>":
                captured["func"] = func
                captured["text"] = self
            return real_tag_bind(self, tagName, sequence, func, add)

        class Click:
            def __init__(self, x, y):
                self.x, self.y = x, y

        pdb = ClashToolTestCase._write_pdb(
            self,
            [
                _atom_line(1, "C1'", "A", 5, 0.0, 0.0, 0.0),
                _atom_line(2, "C1'", "A", 20, 1.2, 0.0, 0.0),
            ],
        )

        result = {}
        real_tk = tk.Tk

        class FakeRoot(real_tk):
            def mainloop(inner, n=0):
                inner.update_idletasks()
                widgets = []

                def walk(w):
                    for child in w.winfo_children():
                        widgets.append(child)
                        walk(child)

                walk(inner)
                run = [
                    w
                    for w in widgets
                    if w.winfo_class() == "TButton" and w.cget("text") == "Run"
                ][0]
                run.invoke()
                inner.update_idletasks()

                text_widget = captured["text"]
                content = text_widget.get("1.0", "end-1c")
                target = None
                for number, line in enumerate(content.splitlines(), 1):
                    hit = line_predicate(line)
                    if hit is not None:
                        target = (number, hit[0], hit[1])
                        break
                self_number, column, expected = target
                middle = column + len(expected) // 2
                text_widget.see(f"{self_number}.{middle}")
                inner.update_idletasks()
                box = text_widget.bbox(f"{self_number}.{middle}")
                result["returned"] = captured["func"](Click(box[0] + 1, box[1] + 1))
                inner.update_idletasks()
                result["expected"] = expected
                result["selected"] = text_widget.get("sel.first", "sel.last")
                result["clipboard"] = inner.clipboard_get()
                result["tag_ranges"] = len(text_widget.tag_ranges("select_command")) // 2
                inner.destroy()

        tk.Text.tag_bind = spy
        tk.Tk = FakeRoot
        try:
            from re_helix_lib import check_pdb_clashes as tool

            tool.run_gui(initial_path=str(pdb))
        finally:
            tk.Tk = real_tk
            tk.Text.tag_bind = real_tag_bind
        return result

    def test_double_click_on_tab_separated_command_selects_whole_command(self):
        """The hard case: a command wedged between tabs in the clash table."""

        def predicate(line):
            if line.startswith("1.200000") and "\tselect " in line:
                column = line.index("\tselect ") + 1
                return column, line[column:].split("\t")[0]
            return None

        result = self._run_gui_and_double_click(predicate)
        self.assertEqual(result["selected"], result["expected"])
        self.assertEqual(result["clipboard"], result["expected"])
        # "break" stops Tk from replacing the selection with a single word.
        self.assertEqual(result["returned"], "break")

    def test_double_click_on_standalone_command_line(self):
        def predicate(line):
            if line.startswith("select #1/"):
                return 0, line
            return None

        result = self._run_gui_and_double_click(predicate)
        self.assertEqual(result["selected"], result["expected"])
        self.assertEqual(result["clipboard"], result["expected"])

    def test_every_command_in_the_report_is_tagged(self):
        def predicate(line):
            if line.startswith("select #0:"):
                return 0, line
            return None

        result = self._run_gui_and_double_click(predicate)
        # Two standalone all-select lines plus two in the clash table.
        self.assertEqual(result["tag_ranges"], 4)


if __name__ == "__main__":
    unittest.main()
