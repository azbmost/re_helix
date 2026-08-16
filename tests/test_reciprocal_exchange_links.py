from __future__ import annotations

import re
import subprocess
import sys
import tempfile
import unittest
from io import StringIO
from pathlib import Path

import re_helix as re_helix_app
from re_helix_lib import edit_pdb_link
from re_helix_lib import reciprocal_exchange_pdbV3_3 as rex


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "re_helix.py"


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


def _synthetic_duplex_text(residue_count: int = 12) -> str:
    lines = []
    serial = 1
    for chain_index, chain_id in enumerate(("A", "B")):
        x = float(chain_index * 10)
        for res_seq in range(1, residue_count + 1):
            p_z = float(res_seq * 4)
            atoms = (
                ("P", x, 0.0, p_z),
                ("OP1", x + 0.8, 0.0, p_z),
                ("OP2", x - 0.8, 0.0, p_z),
                ("O5'", x, 1.6, p_z),
                ("O3'", x, 0.0, p_z + 2.4),
            )
            for atom_name, atom_x, atom_y, atom_z in atoms:
                lines.append(
                    _atom_line(serial, atom_name, chain_id, res_seq, atom_x, atom_y, atom_z)
                )
                serial += 1
        lines.append(f"TER   {serial:5d}       DA {chain_id:1s}{residue_count:4d}\n")
        serial += 1
    return "".join(lines)


def _read_pdb(path: Path):
    atoms = []
    links = []
    with path.open() as handle:
        edit_pdb_link.file2rec_link(handle, atoms, links)
    return atoms, links


def _remark_link_counts(path: Path):
    pattern = re.compile(r"(\w+)=(\d+)")
    for line in path.read_text().splitlines():
        if line.startswith(f"{rex.REMARK_PREFIX} SPECIAL event=link_records "):
            return {key: int(value) for key, value in pattern.findall(line)}
    raise AssertionError("No link_records REMARK found")


def _assert_link_endpoints_resolve(testcase: unittest.TestCase, atoms, links) -> None:
    atom_index = edit_pdb_link.build_atom_index(atoms)
    for link in links:
        endpoint1 = (link.chainID1, link.resSeq1, link.name1.strip())
        endpoint2 = (link.chainID2, link.resSeq2, link.name2.strip())
        testcase.assertIn(endpoint1, atom_index, link.string)
        testcase.assertIn(endpoint2, atom_index, link.string)


class ExistingLinkTopologyTests(unittest.TestCase):
    def test_cir_shift_accepts_optional_closure_suffix(self) -> None:
        self.assertEqual(rex.parse_cir_shift("10"), rex.CirShiftSpec(10, False))
        self.assertEqual(rex.parse_cir_shift("10c"), rex.CirShiftSpec(10, True))
        self.assertEqual(rex.parse_cir_shift("-4C"), rex.CirShiftSpec(-4, True))
        self.assertEqual(rex.parse_cir_shift(0), rex.CirShiftSpec(0, False))
        with self.assertRaises(ValueError):
            rex.parse_cir_shift("c")

    def test_gui_cli_builder_emits_min_link_policy_only_when_selected(self) -> None:
        common = dict(
            script_path=str(SCRIPT),
            pdb_in="input.pdb",
            helix_def_text="(AB)",
            pair_rows=[],
            pair_args_text="A2 B2 d",
            axis_range_rows=[],
            axis_move_rows=[],
            use_user_axis=False,
            user_axis_dir=["0", "0", "1"],
            user_axis_point=["0", "0", "0"],
            output_base="output",
            axis_dist="0",
            axis_parallel="y",
            fix_chain="",
            replicate=False,
            re_only=True,
            cir_shift="0",
            linker_phosphate_resname="X33",
        )
        default_command = re_helix_app._build_equivalent_cli_command(**common)
        minimized_command = re_helix_app._build_equivalent_cli_command(
            **common,
            min_link_records=True,
        )
        self.assertNotIn("--min_link_records", default_command)
        self.assertIn("--min_link_records", minimized_command)

    def test_link_free_input_keeps_numeric_predecessors(self) -> None:
        atoms = []
        links = []
        edit_pdb_link.file2rec_link(StringIO(_synthetic_duplex_text(4)), atoms, links)
        atom_records = [record for record in atoms if record.recordName in {"ATOM", "HETATM"}]
        nodes, label_to_idx = rex.build_residue_nodes(atom_records)

        topology = rex.build_input_backbone_topology(nodes, label_to_idx, links)

        self.assertEqual(topology.recognized_link_count, 0)
        self.assertEqual(len(topology.graph.edges), 6)
        self.assertEqual(
            topology.orig_prev[label_to_idx[("A", 3)]],
            label_to_idx[("A", 2)],
        )

    def test_component_orientation_preserves_input_direction_unless_minimized(self) -> None:
        positional_node = rex.ResidueNode("Z", 9, [], True, ("Z", 8), True, "Q", 7)
        self.assertTrue(positional_node.is_phos_bridge)
        self.assertEqual(positional_node.phos_source, ("Z", 8))
        self.assertEqual(positional_node.new_chain_id, "Q")
        self.assertEqual(positional_node.new_res_seq, 7)
        self.assertEqual(positional_node.input_chain_rank, 0)

        nodes = [
            rex.ResidueNode(
                orig_chain_id="A",
                orig_res_seq=index + 1,
                atoms=[],
                input_chain_rank=index,
            )
            for index in range(5)
        ]

        # In forward input order this cycle contains two inverted standard
        # edges and two special edges. Reversing it saves one LINK, which must
        # not outweigh preserving four directed input adjacencies.
        graph = rex.BackboneGraph(len(nodes))
        graph.add_edge(0, 1, kind="std", end_a="O3", end_b="P")
        graph.add_edge(
            1,
            2,
            kind="3to3",
            end_a="O3",
            end_b="O3",
            phos_key=("B", 1),
        )
        graph.add_edge(2, 3, kind="std", end_a="P", end_b="O3")
        graph.add_edge(3, 4, kind="std", end_a="P", end_b="O3")
        graph.add_edge(4, 0, kind="5to5", end_a="P", end_b="O5")

        for circularized in (False, True):
            components = rex.build_ordered_components(
                graph,
                nodes,
                junction_nodes=set(),
                cir_shift=0,
                circularize_cycles=circularized,
            )
            self.assertEqual(len(components), 1)
            component = components[0]
            self.assertTrue(component["is_cycle"])
            self.assertEqual(component["order"], [0, 1, 2, 3, 4])
            self.assertEqual(component["rotation"], 0)
            self.assertFalse(component["was_reversed"])
            self.assertGreater(
                rex._inverted_cost(component["order"], graph, is_cycle=circularized),
                rex._inverted_cost(
                    list(reversed(component["order"])),
                    graph,
                    is_cycle=circularized,
                ),
            )
            minimized_cycle = rex.build_ordered_components(
                graph,
                nodes,
                junction_nodes=set(),
                cir_shift=0,
                circularize_cycles=circularized,
                min_link_records=True,
            )[0]
            self.assertEqual(minimized_cycle["order"], [0, 4, 3, 2, 1])
            self.assertTrue(minimized_cycle["was_reversed"])
            for shift in (1, -1, 6):
                shifted_minimized = rex.build_ordered_components(
                    graph,
                    nodes,
                    junction_nodes=set(),
                    cir_shift=shift,
                    circularize_cycles=circularized,
                    min_link_records=True,
                )[0]
                shifted_order = shifted_minimized["order"]
                anchor = shifted_order.index(0)
                normalized = shifted_order[anchor:] + shifted_order[:anchor]
                self.assertEqual(normalized, [0, 4, 3, 2, 1])
                self.assertEqual(shifted_minimized["rotation"], shift % 5)

        path_graph = rex.BackboneGraph(4)
        path_nodes = nodes[:4]
        path_graph.add_edge(0, 1, kind="std", end_a="P", end_b="O3")
        path_graph.add_edge(1, 2, kind="std", end_a="P", end_b="O3")
        path_graph.add_edge(2, 3, kind="std", end_a="P", end_b="O3")
        path_components = rex.build_ordered_components(
            path_graph,
            path_nodes,
            junction_nodes=set(),
            cir_shift=0,
        )
        self.assertEqual(path_components[0]["order"], [0, 1, 2, 3])
        minimized_path = rex.build_ordered_components(
            path_graph,
            path_nodes,
            junction_nodes=set(),
            cir_shift=0,
            min_link_records=True,
        )[0]
        self.assertEqual(minimized_path["order"], [3, 2, 1, 0])
        self.assertLess(
            rex._inverted_cost(minimized_path["order"], path_graph),
            rex._inverted_cost(path_components[0]["order"], path_graph),
        )

        # A phosphate-only bridge adjacency always emits a LINK, regardless of
        # traversal direction. The legacy inverted-edge cost incorrectly made
        # the reverse order look cheaper for this path; writer-aware scoring
        # recognizes the tie and retains the forward terminal provenance.
        bridge_records = []
        edit_pdb_link.file2rec_link(
            StringIO(_synthetic_duplex_text(4)),
            bridge_records,
            [],
        )
        bridge_atoms = [
            record
            for record in bridge_records
            if record.recordName in {"ATOM", "HETATM"} and record.chainID == "A"
        ]
        bridge_nodes, _bridge_labels = rex.build_residue_nodes(bridge_atoms)
        bridge_nodes[0].is_phos_bridge = True
        bridge_graph = rex.BackboneGraph(4)
        bridge_graph.add_edge(0, 1, kind="std", end_a="P", end_b="O3")
        bridge_graph.add_edge(1, 2, kind="std", end_a="O3", end_b="P")
        bridge_graph.add_edge(2, 3, kind="std", end_a="P", end_b="O3")
        bridge_minimized = rex.build_ordered_components(
            bridge_graph,
            bridge_nodes,
            junction_nodes=set(),
            cir_shift=0,
            min_link_records=True,
        )[0]
        self.assertEqual(bridge_minimized["order"], [0, 1, 2, 3])
        self.assertLess(
            rex._inverted_cost([3, 2, 1, 0], bridge_graph),
            rex._inverted_cost([0, 1, 2, 3], bridge_graph),
        )
        for candidate_order in ([0, 1, 2, 3], [3, 2, 1, 0]):
            candidate = {
                "order": candidate_order,
                "is_cycle": False,
                "circularized": False,
            }
            output_atoms, _phosphate_labels = rex.assign_new_labels_and_collect_atoms(
                [candidate],
                bridge_nodes,
            )
            emitted_links, _emitted_counts = rex.build_link_records(
                [candidate],
                bridge_graph,
                bridge_nodes,
                output_atoms,
            )
            self.assertEqual(
                rex._output_link_cost(
                    candidate_order,
                    bridge_graph,
                    bridge_nodes,
                    is_cycle=False,
                ),
                2,
            )
            self.assertEqual(len(emitted_links), 2)

        # Model the D10 conflict: both terminal source fragments run forward
        # in the raw direction, while a longer interior fragment makes the
        # reverse direction win both continuity and LINK-count scoring. The
        # default preserves the visible terminal direction; minimization is an
        # explicit opt-in.
        terminal_labels = [
            ("A", 1, 0),
            ("A", 2, 1),
            ("C", 5, 4),
            ("C", 4, 3),
            ("C", 3, 2),
            ("C", 2, 1),
            ("C", 1, 0),
            ("E", 1, 0),
            ("E", 2, 1),
        ]
        terminal_nodes = [
            rex.ResidueNode(chain, res_seq, [], input_chain_rank=rank)
            for chain, res_seq, rank in terminal_labels
        ]
        terminal_graph = rex.BackboneGraph(len(terminal_nodes))
        terminal_graph.add_edge(0, 1, kind="std", end_a="O3", end_b="P")
        terminal_graph.add_edge(1, 2, kind="5to5", end_a="P", end_b="O5")
        for index in range(2, 6):
            terminal_graph.add_edge(index, index + 1, kind="std", end_a="P", end_b="O3")
        terminal_graph.add_edge(6, 7, kind="5to5", end_a="P", end_b="O5")
        terminal_graph.add_edge(7, 8, kind="std", end_a="O3", end_b="P")

        preserved_terminal = rex.build_ordered_components(
            terminal_graph,
            terminal_nodes,
            junction_nodes=set(),
            cir_shift=0,
        )[0]
        minimized_terminal = rex.build_ordered_components(
            terminal_graph,
            terminal_nodes,
            junction_nodes=set(),
            cir_shift=0,
            min_link_records=True,
        )[0]
        self.assertEqual(preserved_terminal["order"], list(range(9)))
        self.assertEqual(minimized_terminal["order"], list(reversed(range(9))))
        self.assertEqual(
            rex._terminal_direction_score(preserved_terminal["order"], terminal_nodes),
            2,
        )
        self.assertGreater(
            rex._input_forward_continuity_score(
                minimized_terminal["order"], terminal_nodes
            ),
            rex._input_forward_continuity_score(
                preserved_terminal["order"], terminal_nodes
            ),
        )
        self.assertLess(
            rex._inverted_cost(minimized_terminal["order"], terminal_graph),
            rex._inverted_cost(preserved_terminal["order"], terminal_graph),
        )

        # A terminal fragment must be established by the pair directly at the
        # exposed end. Do not scan inward and count the same interior fragment
        # at both ends when the true terminal fragments are single residues.
        singleton_terminal_nodes = [
            rex.ResidueNode("A", 1, [], input_chain_rank=0),
            rex.ResidueNode("B", 1, [], input_chain_rank=0),
            rex.ResidueNode("B", 2, [], input_chain_rank=1),
            rex.ResidueNode("C", 1, [], input_chain_rank=0),
        ]
        self.assertEqual(
            rex._terminal_direction_score(
                list(range(4)),
                singleton_terminal_nodes,
            ),
            0,
        )
        self.assertEqual(
            rex._terminal_direction_score(
                list(reversed(range(4))),
                singleton_terminal_nodes,
            ),
            0,
        )

        # A same-chain exchange jump is not a contiguous input fragment merely
        # because its provenance ranks have the same sign.
        jumped_terminal_nodes = [
            rex.ResidueNode("A", 1, [], input_chain_rank=0),
            rex.ResidueNode("A", 9, [], input_chain_rank=8),
            rex.ResidueNode("B", 1, [], input_chain_rank=0),
            rex.ResidueNode("C", 1, [], input_chain_rank=0),
            rex.ResidueNode("C", 9, [], input_chain_rank=8),
        ]
        self.assertEqual(
            rex._terminal_direction_score(
                list(range(5)),
                jumped_terminal_nodes,
            ),
            0,
        )

        # Provenance comes from input record order, not numeric residue order.
        descending_nodes = [
            rex.ResidueNode(
                orig_chain_id="A",
                orig_res_seq=index + 1,
                atoms=[],
                input_chain_rank=2 - index,
            )
            for index in range(3)
        ]
        descending_graph = rex.BackboneGraph(3)
        descending_graph.add_edge(0, 1, kind="std", end_a="O3", end_b="P")
        descending_graph.add_edge(1, 2, kind="std", end_a="O3", end_b="P")
        descending_component = rex.build_ordered_components(
            descending_graph,
            descending_nodes,
            junction_nodes=set(),
            cir_shift=0,
        )[0]
        self.assertEqual(descending_component["order"], [2, 1, 0])

        # When total provenance continuity ties, retain the legacy preference
        # against paths whose fragments at both ends run backward.
        tie_labels = [
            ("A", 2, 1),
            ("A", 1, 0),
            ("B", 1, 0),
            ("B", 2, 1),
            ("B", 3, 2),
            ("C", 2, 1),
            ("C", 1, 0),
        ]
        tie_nodes = [
            rex.ResidueNode(chain, res_seq, [], input_chain_rank=rank)
            for chain, res_seq, rank in tie_labels
        ]
        tie_graph = rex.BackboneGraph(len(tie_nodes))
        for index in range(len(tie_nodes) - 1):
            tie_graph.add_edge(index, index + 1, kind="std", end_a="O3", end_b="P")
        tie_component = rex.build_ordered_components(
            tie_graph,
            tie_nodes,
            junction_nodes=set(),
            cir_shift=0,
        )[0]
        self.assertEqual(tie_component["order"], list(reversed(range(7))))

        # Changing cir_shift may rotate a provenance/link-cost tie, but must
        # never flip the selected cycle direction.
        balanced_nodes = [
            rex.ResidueNode(chain, 1, [], input_chain_rank=0)
            for chain in ("A", "B", "C", "D")
        ]
        balanced_graph = rex.BackboneGraph(4)
        balanced_graph.add_edge(0, 1, kind="std", end_a="O3", end_b="P")
        balanced_graph.add_edge(1, 2, kind="std", end_a="P", end_b="O3")
        balanced_graph.add_edge(2, 3, kind="std", end_a="O3", end_b="P")
        balanced_graph.add_edge(3, 0, kind="std", end_a="P", end_b="O3")
        expected_shift_orders = {
            0: [0, 1, 2, 3],
            1: [1, 2, 3, 0],
            2: [2, 3, 0, 1],
            3: [3, 0, 1, 2],
            4: [0, 1, 2, 3],
            -1: [3, 0, 1, 2],
            -4: [0, 1, 2, 3],
            5: [1, 2, 3, 0],
            -5: [3, 0, 1, 2],
        }
        for shift, expected_order in expected_shift_orders.items():
            balanced_component = rex.build_ordered_components(
                balanced_graph,
                balanced_nodes,
                junction_nodes=set(),
                cir_shift=shift,
                circularize_cycles=True,
            )[0]
            balanced_order = balanced_component["order"]
            self.assertEqual(balanced_order, expected_order)
            self.assertEqual(balanced_component["rotation"], shift % 4)
            anchor = balanced_order.index(0)
            normalized = balanced_order[anchor:] + balanced_order[:anchor]
            self.assertEqual(normalized, [0, 1, 2, 3])

        # Min-LINK mode scores the exact positioned output, including whether
        # the cycle is open or closed. The canonical direction below has the
        # same rotation-independent inversion count as its reverse, but writes
        # one extra LINK at the shift-0 break.
        actual_link_records = []
        edit_pdb_link.file2rec_link(
            StringIO(_synthetic_duplex_text(4)),
            actual_link_records,
            [],
        )
        actual_link_atoms = [
            record
            for record in actual_link_records
            if record.recordName in {"ATOM", "HETATM"} and record.chainID == "A"
        ]
        actual_link_nodes, _actual_labels = rex.build_residue_nodes(actual_link_atoms)
        actual_link_graph = rex.BackboneGraph(4)
        for index, inverted in enumerate((False, True, True, False)):
            neighbor = (index + 1) % 4
            endpoints = ("P", "O3") if inverted else ("O3", "P")
            actual_link_graph.add_edge(
                index,
                neighbor,
                kind="std",
                end_a=endpoints[0],
                end_b=endpoints[1],
            )
        for circularized, expected_cost, alternative_cost in (
            (False, 1, 2),
            (True, 2, 3),
        ):
            actual_minimized = rex.build_ordered_components(
                actual_link_graph,
                actual_link_nodes,
                junction_nodes=set(),
                cir_shift=0,
                circularize_cycles=circularized,
                min_link_records=True,
            )[0]
            self.assertEqual(actual_minimized["order"], [0, 3, 2, 1])
            self.assertEqual(
                rex._output_link_cost(
                    actual_minimized["order"],
                    actual_link_graph,
                    actual_link_nodes,
                    is_cycle=circularized,
                ),
                expected_cost,
            )
            output_atoms, _phosphate_labels = rex.assign_new_labels_and_collect_atoms(
                [actual_minimized],
                actual_link_nodes,
            )
            emitted_links, _emitted_counts = rex.build_link_records(
                [actual_minimized],
                actual_link_graph,
                actual_link_nodes,
                output_atoms,
            )
            self.assertEqual(len(emitted_links), expected_cost)
            self.assertEqual(
                rex._output_link_cost(
                    [0, 1, 2, 3],
                    actual_link_graph,
                    actual_link_nodes,
                    is_cycle=circularized,
                ),
                alternative_cost,
            )

        # The exact-shift contract must not silently move a zero shift merely
        # because the requested break touches reciprocal-exchange junctions.
        junction_component = rex.build_ordered_components(
            balanced_graph,
            balanced_nodes,
            junction_nodes={0, 3},
            cir_shift=0,
            circularize_cycles=False,
        )[0]
        self.assertEqual(junction_component["order"], [0, 1, 2, 3])
        self.assertEqual(junction_component["rotation"], 0)

        # The historical helper remains callable for external users.
        self.assertIsInstance(
            rex.choose_cycle_orientation([0, 1, 2, 3], balanced_graph),
            list,
        )

    def test_chain_ids_and_symmetry_operators_are_distinct_link_identities(self) -> None:
        atoms = []
        links = []
        mixed_case_text = _synthetic_duplex_text(2).replace("B   1", "a   1").replace(
            "B   2", "a   2"
        )
        edit_pdb_link.file2rec_link(StringIO(mixed_case_text), atoms, links)
        atom_records = [record for record in atoms if record.recordName in {"ATOM", "HETATM"}]
        nodes, label_to_idx = rex.build_residue_nodes(atom_records)
        self.assertEqual(len(nodes), 4)
        self.assertIn(("A", 1), label_to_idx)
        self.assertIn(("a", 1), label_to_idx)

        def make_link(chain_id: str, sym2: str = "1555", distance: float = 1.60):
            line = (
                f"LINK          O3'  DA {chain_id:1s}   1                 "
                f"P    DA {chain_id:1s}   2     1555   {sym2:4s} {distance:5.2f}"
            )
            return edit_pdb_link.pdb_link_record(line)

        case_sensitive = rex.merge_link_records(
            [make_link("A")],
            [make_link("a")],
        )
        self.assertEqual(len(case_sensitive), 2)

        symmetry_sensitive = rex.merge_link_records(
            [make_link("A", "1555")],
            [make_link("A", "2555")],
        )
        self.assertEqual(len(symmetry_sensitive), 2)

        for node in nodes:
            node.new_chain_id = node.orig_chain_id
            node.new_res_seq = node.orig_res_seq
        identity, crystal_mate = rex.remap_passthrough_link_records(
            [make_link("A", "1555", 9.99), make_link("A", "2555", 9.99)],
            nodes,
            label_to_idx,
            atom_records,
        )
        self.assertAlmostEqual(identity.distance or 0.0, 1.60, places=2)
        self.assertAlmostEqual(crystal_mate.distance or 0.0, 9.99, places=2)

    def test_generated_rex_output_is_valid_input_for_second_rex(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            input_path = tmp / "input.pdb"
            input_path.write_text(_synthetic_duplex_text())

            first_base = tmp / "first"
            first_run = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
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
            self.assertEqual(first_run.returncode, 0, first_run.stderr)
            first_output = tmp / "first_rex.pdb"
            first_atoms, first_links = _read_pdb(first_output)
            self.assertGreater(len(first_links), 0)
            _assert_link_endpoints_resolve(self, first_atoms, first_links)
            self.assertIn(
                "orientation_policy mode=preserve_terminal_direction min_link_records=0",
                first_output.read_text(),
            )

            first_min_base = tmp / "first_min"
            first_min_run = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    str(input_path),
                    "A7",
                    "B7",
                    "b",
                    "--re_only",
                    "--cir_shift",
                    "0",
                    "--min-link-records",
                    "-o",
                    str(first_min_base),
                ],
                cwd=ROOT,
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(first_min_run.returncode, 0, first_min_run.stderr)
            self.assertIn(
                "orientation_policy mode=min_link_records min_link_records=1",
                (tmp / "first_min_rex.pdb").read_text(),
            )

            # Prove an inherited inverted link changes the chemical predecessor
            # relative to simple residue-number order.
            first_atom_records = [
                record for record in first_atoms if record.recordName in {"ATOM", "HETATM"}
            ]
            first_nodes, first_label_to_idx = rex.build_residue_nodes(first_atom_records)
            first_topology = rex.build_input_backbone_topology(
                first_nodes, first_label_to_idx, first_links
            )
            inherited_inverted = next(
                link
                for link in first_links
                if rex._normalize_backbone_atom_name(link.name1) == "P"
                and rex._normalize_backbone_atom_name(link.name2) == "O3"
                and link.resName1 != rex.X33_HETID
            )
            p_label = (inherited_inverted.chainID1, inherited_inverted.resSeq1)
            o3_label = (inherited_inverted.chainID2, inherited_inverted.resSeq2)
            self.assertEqual(
                first_topology.orig_prev[first_label_to_idx[p_label]],
                first_label_to_idx[o3_label],
            )

            second_open_base = tmp / "second_open"
            second_open_run = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    str(first_output),
                    "A2",
                    "A6",
                    "d",
                    "--re_only",
                    "--cir_shift",
                    "0",
                    "-o",
                    str(second_open_base),
                ],
                cwd=ROOT,
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(second_open_run.returncode, 0, second_open_run.stderr)
            second_open_output = tmp / "second_open_rex.pdb"
            second_open_atoms, second_open_links = _read_pdb(second_open_output)
            second_open_counts = _remark_link_counts(second_open_output)

            self.assertGreater(len(second_open_links), 0)
            self.assertEqual(second_open_counts["total"], len(second_open_links))
            self.assertEqual(second_open_counts["bowtie_5to5"], 1)
            self.assertEqual(second_open_counts["bowtie_3to3"], 2)
            self.assertEqual(second_open_counts["circular_closure"], 0)
            self.assertIn("closure_link=omitted", second_open_output.read_text())
            self.assertEqual(
                sum(
                    1
                    for record in second_open_atoms
                    if record.recordName == "HETATM" and record.resName == rex.X33_HETID
                ),
                3,
            )
            _assert_link_endpoints_resolve(self, second_open_atoms, second_open_links)

            second_closed_base = tmp / "second_closed"
            second_closed_run = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    str(first_output),
                    "A2",
                    "A6",
                    "d",
                    "--re_only",
                    "--cir_shift",
                    "0C",
                    "-o",
                    str(second_closed_base),
                ],
                cwd=ROOT,
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(second_closed_run.returncode, 0, second_closed_run.stderr)
            second_closed_output = tmp / "second_closed_rex.pdb"
            second_closed_atoms, second_closed_links = _read_pdb(second_closed_output)
            second_closed_counts = _remark_link_counts(second_closed_output)

            self.assertEqual(second_closed_counts["total"], len(second_closed_links))
            self.assertEqual(second_closed_counts["bowtie_5to5"], 1)
            self.assertEqual(second_closed_counts["bowtie_3to3"], 2)
            self.assertGreaterEqual(second_closed_counts["circular_closure"], 1)
            self.assertGreater(len(second_closed_links), len(second_open_links))
            self.assertIn("closure_link=written", second_closed_output.read_text())
            _assert_link_endpoints_resolve(self, second_closed_atoms, second_closed_links)


if __name__ == "__main__":
    unittest.main()
