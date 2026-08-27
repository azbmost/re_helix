import unittest

import re_helix
from re_helix_lib.edit_pdb_atom import pdb_atom_record


def _atom_line(serial, atom_name, chain_id, res_seq, x, y, z):
    element = atom_name.strip()[0]
    return (
        f"ATOM  {serial:5d} {atom_name:>4s}  DA {chain_id:1s}{res_seq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}\n"
    )


def _atom(serial, atom_name, chain_id, res_seq, x, y, z):
    return pdb_atom_record(
        _atom_line(serial, atom_name, chain_id, res_seq, x, y, z)
    )


class RestrainedTranslationTests(unittest.TestCase):
    def tearDown(self):
        re_helix.reset_helix_axis_overrides()

    def test_direction_can_come_from_vector_points_or_atoms(self):
        self.assertEqual(
            re_helix.normalize_restrained_translation_direction((0, 0, 4)),
            (0.0, 0.0, 1.0),
        )
        self.assertEqual(
            re_helix.restrained_translation_direction_from_points(
                (1, 2, 3), (1, 5, 3)
            ),
            (0.0, 1.0, 0.0),
        )
        self.assertEqual(
            re_helix.direction_normal_to_vectors((1, 0, 0), (0, 1, 0)),
            (0.0, 0.0, 1.0),
        )

        atoms = [
            _atom(10, "P", "A", 30, 1, 2, 3),
            _atom(11, "C4'", "A", 30, 1, 2, 8),
        ]
        self.assertIs(
            re_helix.resolve_pdb_atom_selector(atoms, "A:30:P"), atoms[0]
        )
        self.assertIs(
            re_helix.resolve_pdb_atom_selector(atoms, "#11"), atoms[1]
        )
        self.assertEqual(
            re_helix.restrained_translation_direction_from_atoms(
                atoms, "A30:P", "30A:C4*"
            ),
            (0.0, 0.0, 1.0),
        )

    def test_zero_length_directions_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "cannot be zero"):
            re_helix.normalize_restrained_translation_direction((0, 0, 0))
        with self.assertRaisesRegex(ValueError, "cannot be zero"):
            re_helix.restrained_translation_direction_from_points(
                (1, 2, 3), (1, 2, 3)
            )
        with self.assertRaisesRegex(ValueError, "zero or parallel"):
            re_helix.direction_normal_to_vectors((1, 0, 0), (2, 0, 0))

    def test_least_squares_translation_uses_all_pairs(self):
        atoms = [
            _atom(1, "P", "A", 1, 0, 0, 0),
            _atom(2, "P", "A", 2, 2, 0, 0),
            _atom(3, "P", "C", 1, 5, 3, 0),
            _atom(4, "P", "C", 2, 9, -1, 0),
        ]
        _chain_atoms, residue_to_p, _chain_p = re_helix.build_nucleic_acid_maps(
            atoms
        )
        translation, signed_distance, sum_sq = (
            re_helix.calculate_restrained_translation(
                [(('A', 1), ('C', 1)), (('A', 2), ('C', 2))],
                residue_to_p,
                (1, 0, 0),
            )
        )
        self.assertEqual(translation, (-6.0, -0.0, -0.0))
        self.assertEqual(signed_distance, -6.0)
        self.assertAlmostEqual(sum_sq, 12.0)

    def test_alignment_honors_move_with_axis_selection_without_rotation(self):
        atoms = [
            _atom(1, "P", "A", 1, 0, 0, 0),
            _atom(2, "P", "B", 1, 0, 1, 0),
            _atom(3, "P", "C", 1, 5, 2, 0),
            _atom(4, "C4'", "C", 1, 6, 4, 1),
            _atom(5, "P", "D", 1, 7, 3, 0),
            _atom(6, "P", "D", 2, 8, 5, 2),
        ]
        original = {
            atom.serial: (atom.x, atom.y, atom.z)
            for atom in atoms
        }
        re_helix.set_helix_axis_range_definitions([{"C": (1, 1)}])
        re_helix.set_helix_axis_move_definitions([{"D": (1, 1)}])
        re_helix.align_helices_for_exchanges(
            atoms,
            [(('A', 1), ('C', 1), "double", None)],
            axis_dist=22.0,
            axis_parallel_flag=True,
            explicit_helices=[('A', 'B'), ('C', 'D')],
            fix_chain="A",
            restrained_translation_direction=(1, 0, 0),
        )

        self.assertEqual((atoms[0].x, atoms[0].y, atoms[0].z), original[1])
        self.assertEqual((atoms[1].x, atoms[1].y, atoms[1].z), original[2])
        for atom in atoms[2:5]:
            old_x, old_y, old_z = original[atom.serial]
            self.assertAlmostEqual(atom.x, old_x - 5.0)
            self.assertAlmostEqual(atom.y, old_y)
            self.assertAlmostEqual(atom.z, old_z)
        self.assertEqual((atoms[5].x, atoms[5].y, atoms[5].z), original[6])

    def test_gui_command_builder_emits_selected_direction_only(self):
        command = re_helix._build_equivalent_cli_command(
            "re_helix.py",
            "input.pdb",
            "(AB) (CD)",
            [{"pos1": "1A", "pos2": "1C", "kind": "d"}],
            "",
            ["A,B"],
            ["D1-D1"],
            False,
            ["0", "0", "1"],
            ["0", "0", "0"],
            "output",
            "22",
            "y",
            "A",
            False,
            False,
            "8",
            "X33",
            alignment_mode="Restrained translation",
            restrained_translation_source="Two PDB atoms",
            restrained_translation_atoms=["#10", "A:30:P"],
        )
        self.assertIn("--alignment_mode", command)
        self.assertIn("--translation_atoms", command)
        self.assertIn("--axis_range", command)
        self.assertIn("--axis_move", command)
        self.assertNotIn("--axis_dist", command)
        self.assertNotIn("--axis_parallel", command)

    def test_restrained_rotation_rotates_without_translation(self):
        atoms = [
            _atom(1, "P", "A", 1, 0, 1, 0),
            _atom(2, "P", "B", 1, 0, 2, 0),
            _atom(3, "P", "C", 1, 1, 0, 0),
            _atom(4, "C4'", "C", 1, 0, 0, 2),
            _atom(5, "P", "D", 1, 2, 0, 0),
            _atom(6, "P", "D", 2, 3, 0, 0),
        ]
        fixed_before = [(atom.x, atom.y, atom.z) for atom in atoms[:2]]
        excluded_payload_before = (atoms[5].x, atoms[5].y, atoms[5].z)
        re_helix.set_helix_axis_range_definitions([{"C": (1, 1)}])
        re_helix.set_helix_axis_move_definitions([{"D": (1, 1)}])
        re_helix.align_helices_for_exchanges(
            atoms,
            [(('A', 1), ('C', 1), "double", None)],
            axis_dist=22.0,
            axis_parallel_flag=True,
            explicit_helices=[('A', 'B'), ('C', 'D')],
            fix_chain="A",
            restrained_rotation_axis=((0, 0, 1), (0, 0, 0)),
        )

        self.assertEqual(
            [(atom.x, atom.y, atom.z) for atom in atoms[:2]],
            fixed_before,
        )
        self.assertAlmostEqual(atoms[2].x, 0.0, places=6)
        self.assertAlmostEqual(atoms[2].y, 1.0, places=6)
        self.assertAlmostEqual(atoms[2].z, 0.0, places=6)
        self.assertAlmostEqual(atoms[3].x, 0.0, places=6)
        self.assertAlmostEqual(atoms[3].y, 0.0, places=6)
        self.assertAlmostEqual(atoms[3].z, 2.0, places=6)
        self.assertAlmostEqual(atoms[4].x, 0.0, places=6)
        self.assertAlmostEqual(atoms[4].y, 2.0, places=6)
        self.assertAlmostEqual(atoms[4].z, 0.0, places=6)
        self.assertEqual(
            (atoms[5].x, atoms[5].y, atoms[5].z),
            excluded_payload_before,
        )

    def test_gui_command_builder_emits_restrained_rotation_axis(self):
        command = re_helix._build_equivalent_cli_command(
            "re_helix.py",
            "input.pdb",
            "(AB) (CD)",
            [{"pos1": "1A", "pos2": "1C", "kind": "d"}],
            "",
            ["C"],
            ["D"],
            False,
            ["0", "0", "1"],
            ["0", "0", "0"],
            "output",
            "22",
            "y",
            "A",
            False,
            False,
            "8",
            "X33",
            alignment_mode="Restrained rotation",
            restrained_rotation_point_source="PDB atom",
            restrained_rotation_point_atom="A:1:P",
            restrained_rotation_vector_source="Two XYZ points",
            restrained_rotation_vector_points=[
                "0", "0", "0", "0", "0", "1"
            ],
        )
        self.assertIn("restrained_rotation", command)
        self.assertIn("--rotation_axis_point_atom", command)
        self.assertIn("--rotation_axis_vector_points", command)
        self.assertIn("--axis_move", command)
        self.assertNotIn("--axis_dist", command)
        self.assertNotIn("--axis_parallel", command)

    def test_gui_command_builder_emits_normal_vector_sources(self):
        common_args = (
            "re_helix.py",
            "input.pdb",
            "(AB) (CD)",
            [{"pos1": "1A", "pos2": "1C", "kind": "d"}],
            "",
            [],
            [],
            False,
            ["0", "0", "1"],
            ["0", "0", "0"],
            "output",
            "22",
            "y",
            "A",
            False,
            False,
            "8",
            "X33",
        )
        translation_command = re_helix._build_equivalent_cli_command(
            *common_args,
            alignment_mode="Restrained translation",
            restrained_translation_source="Normal to two vectors",
            restrained_translation_normal_vectors=[
                "1", "0", "0", "0", "1", "0"
            ],
        )
        self.assertIn("--translation_normal_vectors", translation_command)

        rotation_command = re_helix._build_equivalent_cli_command(
            *common_args,
            alignment_mode="Restrained rotation",
            restrained_rotation_point_source="XYZ point",
            restrained_rotation_point=["0", "0", "0"],
            restrained_rotation_vector_source="Normal to two vectors",
            restrained_rotation_normal_vectors=[
                "1", "0", "0", "0", "1", "0"
            ],
        )
        self.assertIn("--rotation_axis_normal_vectors", rotation_command)


if __name__ == "__main__":
    unittest.main()
