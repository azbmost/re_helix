from __future__ import annotations

import subprocess
import sys
import tempfile
import unittest
from collections import Counter
from io import StringIO
from pathlib import Path

from re_helix_lib import edit_pdb_link
from re_helix_lib import reciprocal_exchange_pdbV3_3 as rex
from re_helix_lib import reverse_strand_direction as reverse_tool


ROOT = Path(__file__).resolve().parents[1]
RE_HELIX = ROOT / "re_helix.py"


def _atom_line(
    serial: int,
    atom_name: str,
    chain_id: str,
    res_seq: int,
    x: float,
    y: float,
    z: float,
) -> str:
    element = "P" if atom_name == "P" else "O"
    return (
        f"ATOM  {serial:5d} {atom_name:>4s}  DA {chain_id:1s}{res_seq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}          {element:>2s}\n"
    )


def _synthetic_duplex_text(residue_count: int = 6) -> str:
    lines = []
    serial = 1
    for chain_index, chain_id in enumerate(("A", "B")):
        x = float(chain_index * 10)
        for res_seq in range(1, residue_count + 1):
            p_z = float(res_seq * 4)
            for atom_name, atom_x, atom_y, atom_z in (
                ("P", x, 0.0, p_z),
                ("OP1", x + 0.8, 0.0, p_z),
                ("OP2", x - 0.8, 0.0, p_z),
                ("O5'", x, 1.6, p_z),
                ("O3'", x, 0.0, p_z + 2.4),
            ):
                lines.append(
                    _atom_line(serial, atom_name, chain_id, res_seq, atom_x, atom_y, atom_z)
                )
                serial += 1
        lines.append(f"TER   {serial:5d}       DA {chain_id:1s}{residue_count:4d}\n")
        serial += 1
    return "".join(lines)


def _parse_text(text: str):
    records = []
    links = []
    edit_pdb_link.file2rec_link(StringIO(text), records, links)
    atoms = [
        record
        for record in records
        if getattr(record, "recordName", "") in {"ATOM", "HETATM"}
    ]
    return atoms, links


def _assert_link_endpoints_resolve(testcase: unittest.TestCase, atoms, links) -> None:
    atom_index = edit_pdb_link.build_atom_index(atoms)
    for link in links:
        testcase.assertIn(
            (link.chainID1, link.resSeq1, link.name1.strip()),
            atom_index,
            link.string,
        )
        testcase.assertIn(
            (link.chainID2, link.resSeq2, link.name2.strip()),
            atom_index,
            link.string,
        )


def _coordinate_lines(text: str):
    return [line for line in text.splitlines() if line[:6] in {"ATOM  ", "HETATM"}]


def _chain_source_order(text: str, chain_id: str):
    order = []
    seen = set()
    for line in _coordinate_lines(text):
        if line[21] != chain_id:
            continue
        key = int(line[22:26])
        if key in seen:
            continue
        seen.add(key)
        # The synthetic z coordinate identifies the original residue number.
        original_resseq = int(round(float(line[46:54]) / 4.0))
        order.append((key, original_resseq))
    return order


def _endpoint_coordinate(node, endpoint: str):
    candidates = {
        "P": {"P"},
        "O3": {"O3'", "O3*", "O3"},
        "O5": {"O5'", "O5*", "O5"},
    }[endpoint]
    for atom in node.atoms:
        if atom.name.strip() in candidates:
            return (endpoint, round(atom.x, 3), round(atom.y, 3), round(atom.z, 3))
    raise AssertionError(f"Missing {endpoint} in {node.orig_label()}")


def _typed_topology_by_coordinates(text: str):
    atoms, links = _parse_text(text)
    nodes, labels = rex.build_residue_nodes(atoms)
    topology = rex.build_input_backbone_topology(nodes, labels, links)
    result = set()
    for edge in topology.graph.edges.values():
        endpoint_a = _endpoint_coordinate(nodes[edge.a], edge.end_a)
        endpoint_b = _endpoint_coordinate(nodes[edge.b], edge.end_b)
        result.add((edge.kind, frozenset((endpoint_a, endpoint_b))))
    return result


class ReverseStrandDirectionTests(unittest.TestCase):
    def test_natural_chain_reversal_adds_inverted_links_and_round_trips(self) -> None:
        source = _synthetic_duplex_text(3)
        reversed_lines, stats, topology_count, preserved_count, total_count = (
            reverse_tool.reverse_strand_lines(source.splitlines(True), ["A"])
        )
        reversed_text = "".join(reversed_lines)
        self.assertEqual(_chain_source_order(reversed_text, "A"), [(1, 3), (2, 2), (3, 1)])
        self.assertEqual(_chain_source_order(reversed_text, "B"), [(1, 1), (2, 2), (3, 3)])
        self.assertEqual(len(stats), 1)
        self.assertFalse(stats[0].is_cycle)
        self.assertEqual((topology_count, preserved_count, total_count), (2, 0, 2))
        atoms, links = _parse_text(reversed_text)
        _assert_link_endpoints_resolve(self, atoms, links)
        self.assertCountEqual(
            [
                (link.chainID1, link.resSeq1, link.name1.strip(), link.chainID2, link.resSeq2, link.name2.strip())
                for link in links
            ],
            [
                ("A", 1, "P", "A", 2, "O3'"),
                ("A", 2, "P", "A", 3, "O3'"),
            ],
        )

        roundtrip_lines, _stats, _topology, _preserved, roundtrip_links = (
            reverse_tool.reverse_strand_lines(reversed_lines, ["A"])
        )
        roundtrip_text = "".join(roundtrip_lines)
        self.assertEqual(_coordinate_lines(roundtrip_text), _coordinate_lines(source))
        self.assertEqual(roundtrip_links, 0)
        self.assertEqual(
            _typed_topology_by_coordinates(roundtrip_text),
            _typed_topology_by_coordinates(source),
        )

    def test_multiple_selection_is_case_sensitive_and_order_independent(self) -> None:
        source = _synthetic_duplex_text(3)
        first = reverse_tool.reverse_strand_lines(source.splitlines(True), ["B", "A"])[0]
        second = reverse_tool.reverse_strand_lines(source.splitlines(True), ["A", "B"])[0]
        self.assertEqual(first, second)
        text = "".join(first)
        self.assertEqual(_chain_source_order(text, "A"), [(1, 3), (2, 2), (3, 1)])
        self.assertEqual(_chain_source_order(text, "B"), [(1, 3), (2, 2), (3, 1)])
        with self.assertRaisesRegex(ValueError, "not found"):
            reverse_tool.reverse_strand_lines(source.splitlines(True), ["a"])

    def test_blank_chain_selector_and_cli_alias(self) -> None:
        named_source = _synthetic_duplex_text(3)
        atoms, _links = _parse_text(named_source)
        atom_index = edit_pdb_link.build_atom_index(atoms)
        blank_passthrough = rex._format_link_line(
            atom_index[("A", 1, "OP1")][0],
            atom_index[("B", 1, "OP2")][0],
        )
        blank_passthrough = blank_passthrough[:21] + " " + blank_passthrough[22:]
        blank_passthrough_endpoint2 = rex._format_link_line(
            atom_index[("B", 2, "OP2")][0],
            atom_index[("A", 1, "OP1")][0],
        )
        blank_passthrough_endpoint2 = (
            blank_passthrough_endpoint2[:51]
            + " "
            + blank_passthrough_endpoint2[52:]
        )

        source_lines = named_source.splitlines(True)
        for index, line in enumerate(source_lines):
            if line[:6].strip() in {"ATOM", "HETATM", "TER"} and line[21] == "A":
                source_lines[index] = line[:21] + " " + line[22:]
        source_lines.insert(0, blank_passthrough)
        source_lines.insert(1, blank_passthrough_endpoint2)

        reversed_lines, stats, topology_count, preserved_count, total_count = (
            reverse_tool.reverse_strand_lines(source_lines, ["blank"])
        )
        reversed_text = "".join(reversed_lines)
        self.assertEqual(_chain_source_order(reversed_text, " "), [(1, 3), (2, 2), (3, 1)])
        self.assertEqual(stats[0].chain_id, " ")
        self.assertEqual((topology_count, preserved_count, total_count), (2, 2, 4))
        self.assertIn("chains=blank topology=preserved", reversed_text)
        parsed_output_links = [
            reverse_tool._parse_link_line(line + "\n")
            for line in reversed_text.splitlines()
            if line.startswith("LINK")
        ]
        remapped_passthrough = next(
            link
            for link in parsed_output_links
            if link.name1.strip() == "OP1" and link.chainID1 == " "
        )
        self.assertEqual((remapped_passthrough.chainID1, remapped_passthrough.resSeq1), (" ", 3))
        remapped_endpoint2 = next(
            link
            for link in parsed_output_links
            if link.name2.strip() == "OP1" and link.chainID2 == " "
        )
        self.assertEqual((remapped_endpoint2.chainID2, remapped_endpoint2.resSeq2), (" ", 3))

        roundtrip_lines = reverse_tool.reverse_strand_lines(reversed_lines, ["_"])[0]
        roundtrip_text = "".join(roundtrip_lines)
        self.assertEqual(_coordinate_lines(roundtrip_text), _coordinate_lines("".join(source_lines)))
        command = reverse_tool.build_cli_command(
            "reverse_strand_direction.py",
            "input.pdb",
            [" "],
            "output.pdb",
        )
        self.assertIn("--strand blank", command)

    def test_unselected_topology_is_not_validated_and_links_stay_raw(self) -> None:
        source = _synthetic_duplex_text(5)
        atoms, _links = _parse_text(source)
        atom_index = edit_pdb_link.build_atom_index(atoms)

        # This non-adjacent supported link makes chain B branched (degree 3),
        # but B is outside the requested edit and must remain untouched.
        unselected_branch = rex._format_link_line(
            atom_index[("B", 1, "P")][0],
            atom_index[("B", 3, "O3'")][0],
        ).rstrip("\n") + "   \n"
        unselected_other = rex._format_link_line(
            atom_index[("B", 2, "OP1")][0],
            atom_index[("B", 4, "OP2")][0],
        ).rstrip("\n") + "  \n"
        changed_cross_link = rex._format_link_line(
            atom_index[("A", 1, "OP1")][0],
            atom_index[("B", 1, "OP2")][0],
        )
        unchanged_cross_link = rex._format_link_line(
            atom_index[("A", 3, "OP1")][0],
            atom_index[("B", 3, "OP2")][0],
        ).rstrip("\n") + "    \n"

        coordinate_lines = source.splitlines(True)
        first_b = next(
            index
            for index, line in enumerate(coordinate_lines)
            if line[:6] == "ATOM  " and line[21] == "B"
        )
        # An insertion code on unselected B is another unsupported condition
        # that must not prevent the scoped reversal of A.
        coordinate_lines[first_b] = (
            coordinate_lines[first_b][:26] + "Q" + coordinate_lines[first_b][27:]
        )
        source_with_links = "".join(
            [
                unselected_branch,
                unselected_other,
                changed_cross_link,
                unchanged_cross_link,
                *coordinate_lines,
            ]
        )

        reversed_lines, _stats, topology_count, preserved_count, total_count = (
            reverse_tool.reverse_strand_lines(source_with_links.splitlines(True), ["A"])
        )
        reversed_text = "".join(reversed_lines)
        self.assertEqual(_chain_source_order(reversed_text, "A"), [(1, 5), (2, 4), (3, 3), (4, 2), (5, 1)])
        self.assertEqual((topology_count, preserved_count, total_count), (4, 4, 8))

        output_link_lines = [
            line for line in reversed_text.splitlines() if line.startswith("LINK")
        ]
        self.assertIn(unselected_branch.rstrip("\n"), output_link_lines)
        self.assertIn(unselected_other.rstrip("\n"), output_link_lines)
        self.assertIn(unchanged_cross_link.rstrip("\n"), output_link_lines)

        _output_atoms, output_links = _parse_text(reversed_text)
        remapped_cross = next(
            link
            for link in output_links
            if link.name1.strip() == "OP1"
            and link.chainID1 == "A"
            and link.chainID2 == "B"
            and link.resSeq2 == 1
        )
        self.assertEqual(remapped_cross.resSeq1, 5)

    def test_selected_backbone_link_to_unselected_chain_is_rejected(self) -> None:
        source = _synthetic_duplex_text(3)
        atoms, _links = _parse_text(source)
        atom_index = edit_pdb_link.build_atom_index(atoms)
        cross_chain_backbone = rex._format_link_line(
            atom_index[("A", 1, "P")][0],
            atom_index[("B", 1, "O3'")][0],
        )
        with self.assertRaisesRegex(ValueError, "backbone connection to another chain"):
            reverse_tool.reverse_strand_lines(
                (cross_chain_backbone + source).splitlines(True),
                ["A"],
            )

    def test_re_helix_bowtie_topology_survives_reversal(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            input_path = tmp / "input.pdb"
            input_path.write_text(_synthetic_duplex_text(12))
            first_base = tmp / "first"
            run = subprocess.run(
                [
                    sys.executable,
                    str(RE_HELIX),
                    str(input_path),
                    "A7",
                    "B7",
                    "b",
                    "--re_only",
                    "--cir_shift",
                    "0",
                    "-o",
                    str(first_base),
                ],
                cwd=ROOT,
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(run.returncode, 0, run.stderr)
            generated = (tmp / "first_rex.pdb").read_text()
            reversed_lines, _stats, _topology, _preserved, _total = (
                reverse_tool.reverse_strand_lines(generated.splitlines(True), ["A"])
            )
            reversed_text = "".join(reversed_lines)
            self.assertEqual(
                _typed_topology_by_coordinates(reversed_text),
                _typed_topology_by_coordinates(generated),
            )
            generated_signatures = Counter(
                line[:17] + line[27:] for line in _coordinate_lines(generated)
            )
            reversed_signatures = Counter(
                line[:17] + line[27:] for line in _coordinate_lines(reversed_text)
            )
            self.assertEqual(reversed_signatures, generated_signatures)
            atoms, links = _parse_text(reversed_text)
            _assert_link_endpoints_resolve(self, atoms, links)
            self.assertGreater(len(links), 0)

    def test_closed_cycle_keeps_anchor_and_passthrough_symmetry(self) -> None:
        source = _synthetic_duplex_text(3)
        atoms, _links = _parse_text(source)
        atom_index = edit_pdb_link.build_atom_index(atoms)
        closure = rex._format_link_line(
            atom_index[("A", 1, "P")][0],
            atom_index[("A", 3, "O3'")][0],
        )
        crystal_link = edit_pdb_link.pdb_link_record(
            "LINK OP1 DA A 1 OP2 DA B 1 1555 2555 9.99\n"
        ).string
        source_with_links = closure + crystal_link + source

        reversed_lines, stats, _topology, preserved, total = (
            reverse_tool.reverse_strand_lines(
                source_with_links.splitlines(True),
                ["A"],
            )
        )
        reversed_text = "".join(reversed_lines)
        self.assertTrue(stats[0].is_cycle)
        self.assertEqual(_chain_source_order(reversed_text, "A"), [(1, 1), (2, 3), (3, 2)])
        self.assertEqual(preserved, 1)
        self.assertEqual(total, 4)
        reversed_atoms, reversed_links = _parse_text(reversed_text)
        _assert_link_endpoints_resolve(self, reversed_atoms, reversed_links)
        remapped_crystal = next(link for link in reversed_links if link.sym2 == "2555")
        self.assertEqual((remapped_crystal.chainID1, remapped_crystal.resSeq1), ("A", 1))
        self.assertAlmostEqual(remapped_crystal.distance or 0.0, 9.99, places=2)
        self.assertEqual(
            _typed_topology_by_coordinates(reversed_text),
            _typed_topology_by_coordinates(source_with_links),
        )

        roundtrip_lines, roundtrip_stats, _topology2, preserved2, total2 = (
            reverse_tool.reverse_strand_lines(reversed_lines, ["A"])
        )
        roundtrip_text = "".join(roundtrip_lines)
        self.assertTrue(roundtrip_stats[0].is_cycle)
        self.assertEqual(_coordinate_lines(roundtrip_text), _coordinate_lines(source))
        self.assertEqual((preserved2, total2), (1, 2))
        self.assertEqual(
            _typed_topology_by_coordinates(roundtrip_text),
            _typed_topology_by_coordinates(source_with_links),
        )

    def test_rejects_multimodel_and_insertion_codes(self) -> None:
        source = _synthetic_duplex_text(2)
        with self.assertRaisesRegex(ValueError, "Multi-model"):
            reverse_tool.reverse_strand_lines(
                ["MODEL        1\n", *source.splitlines(True), "ENDMDL\n"],
                ["A"],
            )
        insertion_lines = source.splitlines(True)
        first_a = next(index for index, line in enumerate(insertion_lines) if line[:6] == "ATOM  ")
        line = insertion_lines[first_a]
        insertion_lines[first_a] = line[:26] + "A" + line[27:]
        insertion = "".join(insertion_lines)
        with self.assertRaisesRegex(ValueError, "Insertion codes"):
            reverse_tool.reverse_strand_lines(insertion.splitlines(True), ["A"])

    def test_rejects_ambiguous_two_residue_cycle_before_writing(self) -> None:
        source = _synthetic_duplex_text(2)
        atoms, _links = _parse_text(source)
        atom_index = edit_pdb_link.build_atom_index(atoms)
        ambiguous_link = rex._format_link_line(
            atom_index[("A", 1, "P")][0],
            atom_index[("A", 2, "O3'")][0],
        )
        ambiguous_source = ambiguous_link + source

        with self.assertRaisesRegex(ValueError, "ambiguous two-residue topology"):
            reverse_tool.reverse_strand_lines(
                ambiguous_source.splitlines(True),
                ["A"],
            )

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            input_path = tmp / "ambiguous.pdb"
            output_path = tmp / "must_not_exist.pdb"
            input_path.write_text(ambiguous_source)
            with self.assertRaisesRegex(ValueError, "ambiguous two-residue topology"):
                reverse_tool.reverse_strands(
                    input_path,
                    ["A"],
                    output_pdb=output_path,
                    verbose=False,
                )
            self.assertFalse(output_path.exists())

if __name__ == "__main__":
    unittest.main()
