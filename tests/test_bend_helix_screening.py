import inspect
import math
import re
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from re_helix_lib import bend_helix as bend


def _atom_line(serial, atom_name, chain_id, res_seq, x, y, z):
    element = "P" if atom_name.strip().upper() == "P" else "C"
    return (
        f"ATOM  {serial:5d} {atom_name:>4s}  DA {chain_id:1s}{res_seq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}\n"
    )


def _screening_fixture_pdb():
    """Return an ideal X/Y duplex whose overlay namespace must become A/B/C/D."""
    atoms = [
        ("P", "X", 1, 1.0, 0.0, 0.0),
        ("C4'", "X", 1, 0.0, 0.0, 0.0),
        ("P", "X", 2, 1.0, 0.0, 2.0),
        ("C4'", "X", 2, 0.0, 0.0, 2.0),
        ("P", "X", 3, 1.0, 0.0, 4.0),
        ("C4'", "X", 3, 0.0, 0.0, 4.0),
        ("P", "Y", 3, -1.0, 0.0, 0.0),
        ("C4'", "Y", 3, 0.0, 0.0, 0.0),
        ("P", "Y", 2, -1.0, 0.0, 2.0),
        ("C4'", "Y", 2, 0.0, 0.0, 2.0),
        ("P", "Y", 1, -1.0, 0.0, 4.0),
        ("C4'", "Y", 1, 0.0, 0.0, 4.0),
    ]
    lines = [
        _atom_line(index, *atom)
        for index, atom in enumerate(atoms, start=1)
    ]
    lines.append("END\n")
    return "".join(lines)


class BendHelixScreeningTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.input_path = Path(self.temp_dir.name) / "xy_duplex.pdb"
        self.input_path.write_text(_screening_fixture_pdb(), encoding="utf-8")

    def tearDown(self):
        self.temp_dir.cleanup()

    def context(self, *, local_axis=False):
        ranges = ["X1-X3,Y3-Y1"] if local_axis else None
        return bend.prepare_screening_context(
            str(self.input_path),
            "X2",
            axis_range_specs=ranges,
        )

    @staticmethod
    def atom(selector):
        return bend.ScreeningPoint("overlay_atom", selector)

    @staticmethod
    def xyz(x, y, z):
        return bend.ScreeningPoint("xyz", (x, y, z))

    def test_overlay_namespace_maps_nonstandard_input_chains(self):
        context = self.context()
        self.assertEqual(context.origin_chain_map_model1, {"X": "A", "Y": "B"})
        self.assertEqual(context.origin_chain_map_model2, {"X": "C", "Y": "D"})
        self.assertEqual(
            context.overlay_chain_map,
            {
                "A": ("X", False),
                "B": ("Y", False),
                "C": ("X", True),
                "D": ("Y", True),
            },
        )

        transform = bend.build_bend_transform(
            context.preparation, 0.0, 0.0, 90.0, align_mode="n"
        )
        original = bend.resolve_overlay_atom_coordinate(context, "A:1:P", transform)
        transformed = bend.resolve_overlay_atom_coordinate(context, "C:1:P", transform)
        self.assertEqual(original, (1.0, 0.0, 0.0))
        self.assertAlmostEqual(transformed[0], 0.0, places=7)
        self.assertAlmostEqual(transformed[1], 1.0, places=7)
        self.assertAlmostEqual(transformed[2], 0.0, places=7)

    def test_distance_screening_supports_atom_atom_and_atom_xyz(self):
        context = self.context()
        beta_range = [bend.ScreenAngleRange("beta", 0.0, 90.0, 45.0)]
        fixed = {"phi": 0.0, "beta": 0.0, "tau": 0.0}
        target = math.sqrt(8.0)

        atom_atom = bend.ScreeningRequest(
            mode="distance",
            target=target,
            point1=self.atom("A:1:P"),
            point2=self.atom("C:1:P"),
        )
        result = bend.screen_bend_angles(
            context, fixed, beta_range, atom_atom, align_mode="y"
        )
        self.assertEqual(result.candidate_count, 3)
        self.assertAlmostEqual(result.beta_deg, 90.0)
        self.assertAlmostEqual(result.achieved_value, target, places=7)
        self.assertAlmostEqual(result.error, 0.0, places=7)

        atom_xyz = bend.ScreeningRequest(
            mode="distance",
            target=target,
            point1=self.atom("C:1:P"),
            point2=self.xyz(1.0, 0.0, 0.0),
        )
        result = bend.screen_bend_angles(
            context, fixed, beta_range, atom_xyz, align_mode="y"
        )
        self.assertAlmostEqual(result.beta_deg, 90.0)
        self.assertAlmostEqual(result.achieved_value, target, places=7)

    def test_distance_screening_rejects_negative_target(self):
        request = bend.ScreeningRequest(
            mode="distance",
            target=-1.0,
            point1=self.atom("A:1:P"),
            point2=self.atom("C:1:P"),
        )
        with self.assertRaisesRegex(ValueError, "cannot be negative"):
            bend.screen_bend_angles(
                self.context(),
                {"phi": 0.0, "tau": 0.0},
                [bend.ScreenAngleRange("beta", 0.0, 90.0, 45.0)],
                request,
                align_mode="y",
            )

    def test_rotation_screening_uses_geometric_overlay_atom_axis(self):
        context = self.context()
        axis = bend.ScreeningAxis(
            source="geometric",
            point=self.atom("A:1:C4'"),
            vector_source="two_overlay_atoms",
            point1=self.atom("A:1:C4'"),
            point2=self.atom("A:3:C4'"),
        )
        request = bend.ScreeningRequest(
            mode="rotation",
            target=80.0,
            point1=self.atom("A:1:P"),
            point2=self.atom("C:1:P"),
            axis=axis,
        )
        result = bend.screen_bend_angles(
            context,
            {"phi": 0.0, "beta": 0.0, "tau": 0.0},
            [bend.ScreenAngleRange("tau", 0.0, 90.0, 45.0)],
            request,
            align_mode="n",
        )
        self.assertAlmostEqual(result.tau_deg, 80.0, delta=0.001)
        self.assertAlmostEqual(result.achieved_value, 80.0, delta=0.001)
        self.assertLess(result.error, 0.001)
        self.assertGreater(result.refinement_candidate_count, 0)

    def test_rotation_screening_supports_atom_xyz(self):
        context = self.context()
        request = bend.ScreeningRequest(
            mode="rotation",
            target=90.0,
            point1=self.atom("A:1:P"),
            point2=self.xyz(0.0, 1.0, 0.0),
            axis=bend.ScreeningAxis(
                source="geometric",
                point=self.xyz(0.0, 0.0, 0.0),
                vector_source="direct_vector",
                vector=(0.0, 0.0, 1.0),
            ),
        )
        result = bend.screen_bend_angles(
            context,
            {"phi": 0.0, "beta": 0.0, "tau": 0.0},
            [bend.ScreenAngleRange("tau", 0.0, 0.0, 1.0)],
            request,
            align_mode="n",
        )
        self.assertEqual(result.candidate_count, 1)
        self.assertAlmostEqual(result.achieved_value, 90.0, places=7)
        self.assertAlmostEqual(result.error, 0.0, places=7)

    def test_invalid_rotation_candidate_is_skipped_and_later_candidate_wins(self):
        context = self.context()
        request = bend.ScreeningRequest(
            mode="rotation",
            target=180.0,
            point1=self.atom("A:1:P"),
            point2=self.atom("C:1:P"),
            axis=bend.ScreeningAxis(
                source="geometric",
                point=self.xyz(0.0, 0.0, 0.0),
                vector_source="direct_vector",
                vector=(0.0, 0.0, 1.0),
            ),
        )
        result = bend.screen_bend_angles(
            context,
            {"phi": 0.0, "tau": 0.0},
            [bend.ScreenAngleRange("beta", 30.0, 90.0, 60.0)],
            request,
            align_mode="y",
        )

        self.assertEqual(result.candidate_count, 2)
        self.assertGreater(result.beta_deg, 30.0)
        self.assertLess(result.beta_deg, 31.0)
        self.assertAlmostEqual(result.achieved_value, 180.0, places=7)
        self.assertAlmostEqual(result.error, 0.0, places=7)
        self.assertGreater(result.refinement_candidate_count, 0)

    def test_all_invalid_rotation_candidates_raise_clear_error(self):
        request = bend.ScreeningRequest(
            mode="rotation",
            target=90.0,
            point1=self.atom("A:1:C4'"),
            point2=self.atom("C:1:P"),
            axis=bend.ScreeningAxis(
                source="geometric",
                point=self.xyz(0.0, 0.0, 0.0),
                vector_source="direct_vector",
                vector=(0.0, 0.0, 1.0),
            ),
        )
        with self.assertRaisesRegex(
            ValueError, "did not contain any geometrically valid candidates"
        ):
            bend.screen_bend_angles(
                self.context(),
                {"phi": 0.0, "tau": 0.0},
                [bend.ScreenAngleRange("beta", 0.0, 90.0, 90.0)],
                request,
                align_mode="y",
            )

    def test_rotation_screening_uses_local_axis_range(self):
        context = self.context(local_axis=True)
        request = bend.ScreeningRequest(
            mode="rotation",
            target=-90.0,
            point1=self.atom("A:1:P"),
            point2=self.atom("C:1:P"),
            axis=bend.ScreeningAxis(source="local_axis"),
        )
        result = bend.screen_bend_angles(
            context,
            {"phi": 0.0, "beta": 0.0, "tau": 0.0},
            [bend.ScreenAngleRange("tau", -90.0, 0.0, 45.0)],
            request,
            align_mode="n",
        )
        self.assertEqual(context.preparation.axis_range_used, "X1-X3,Y3-Y1")
        self.assertAlmostEqual(result.tau_deg, -90.0)
        self.assertAlmostEqual(result.achieved_value, -90.0, places=7)

    def test_popup_local_axis_ranges_override_main_ranges_for_screening(self):
        popup_specs = bend.split_axis_range_spec_text(
            "X1-X3,Y3-Y1; X2-X3,Y2-Y1"
        )
        local_request = bend.ScreeningRequest(
            mode="rotation",
            target=0.0,
            point1=self.atom("A:1:P"),
            point2=self.atom("C:1:P"),
            axis=bend.ScreeningAxis(source="local_axis"),
        )

        self.assertEqual(
            popup_specs,
            ["X1-X3,Y3-Y1", "X2-X3,Y2-Y1"],
        )
        self.assertEqual(
            bend.select_screening_axis_range_specs(
                local_request,
                main_axis_range_specs=["X1-X2,Y3-Y2"],
                popup_local_axis_range_specs=popup_specs,
            ),
            popup_specs,
        )
        geometric_request = bend.ScreeningRequest(
            mode="rotation",
            target=0.0,
            point1=self.atom("A:1:P"),
            point2=self.atom("C:1:P"),
            axis=bend.ScreeningAxis(source="geometric"),
        )
        self.assertEqual(
            bend.select_screening_axis_range_specs(
                geometric_request,
                main_axis_range_specs=["X1-X2,Y3-Y2"],
                popup_local_axis_range_specs=popup_specs,
            ),
            ["X1-X2,Y3-Y2"],
        )
        with self.assertRaisesRegex(ValueError, "Screening to achieve window"):
            bend.select_screening_axis_range_specs(
                local_request,
                main_axis_range_specs=["X1-X2,Y3-Y2"],
                popup_local_axis_range_specs=[],
            )

    def test_phi_corrected_pivot_is_recomputed_for_every_candidate(self):
        context = self.context(local_axis=True)
        request = bend.ScreeningRequest(
            mode="rotation",
            target=90.0,
            point1=self.atom("A:2:P"),
            point2=bend.ScreeningPoint("phi_corrected_pivot"),
            axis=bend.ScreeningAxis(source="local_axis"),
        )
        result = bend.screen_bend_angles(
            context,
            {"phi": 0.0, "beta": 0.0, "tau": 0.0},
            [bend.ScreenAngleRange("phi", 0.0, 90.0, 45.0)],
            request,
            align_mode="n",
        )
        self.assertAlmostEqual(result.phi_deg, 90.0)
        self.assertAlmostEqual(result.achieved_value, 90.0, places=7)

    def test_grid_is_inclusive_and_requires_exactly_one_or_two_angles(self):
        angle_range = bend.ScreenAngleRange("phi", 0.0, 1.0, 0.5)
        self.assertEqual(bend.inclusive_angle_values(angle_range), (0.0, 0.5, 1.0))
        self.assertEqual(
            bend.inclusive_angle_values(
                bend.ScreenAngleRange("phi", 0.0, 1.0, 0.4)
            ),
            (0.0, 0.4, 0.8, 1.0),
        )

        with self.assertRaisesRegex(ValueError, "one or two"):
            bend.validate_screen_angle_ranges([], candidate_cap=100)
        with self.assertRaisesRegex(ValueError, "one or two"):
            bend.validate_screen_angle_ranges(
                [
                    bend.ScreenAngleRange("phi", 0.0, 1.0, 1.0),
                    bend.ScreenAngleRange("beta", 0.0, 1.0, 1.0),
                    bend.ScreenAngleRange("tau", 0.0, 1.0, 1.0),
                ],
                candidate_cap=100,
            )

        _ranges, values, total = bend.validate_screen_angle_ranges(
            [
                bend.ScreenAngleRange("phi", 0.0, 1.0, 1.0),
                bend.ScreenAngleRange("beta", 0.0, 2.0, 1.0),
            ],
            candidate_cap=100,
        )
        self.assertEqual(values, ((0.0, 1.0), (0.0, 1.0, 2.0)))
        self.assertEqual(total, 6)

        with self.assertRaisesRegex(ValueError, "exceeding"):
            bend.validate_screen_angle_ranges(
                [
                    bend.ScreenAngleRange("phi", 0.0, 1.0, 1.0),
                    bend.ScreenAngleRange("beta", 0.0, 1.0, 1.0),
                ],
                candidate_cap=3,
            )

        two_angle_result = bend.screen_bend_angles(
            self.context(),
            {"phi": 12.0, "beta": 34.0, "tau": 0.0},
            [
                bend.ScreenAngleRange("phi", 0.0, 90.0, 90.0),
                bend.ScreenAngleRange("beta", 0.0, 90.0, 90.0),
            ],
            bend.ScreeningRequest(
                mode="distance",
                target=0.0,
                point1=self.atom("A:1:P"),
                point2=self.atom("C:1:P"),
            ),
            align_mode="y",
        )
        self.assertEqual(two_angle_result.candidate_count, 4)
        self.assertEqual(two_angle_result.angles, {"phi": 0.0, "beta": 0.0, "tau": 0.0})

    def test_adaptive_refinement_finds_value_between_coarse_steps(self):
        target_tau = 30.0
        target_distance = 2.0 * math.sin(math.radians(target_tau / 2.0))
        result = bend.screen_bend_angles(
            self.context(),
            {"phi": 0.0, "beta": 0.0},
            [bend.ScreenAngleRange("tau", 0.0, 90.0, 45.0)],
            bend.ScreeningRequest(
                mode="distance",
                target=target_distance,
                point1=self.atom("A:1:P"),
                point2=self.atom("C:1:P"),
            ),
            align_mode="n",
        )

        self.assertEqual(result.candidate_count, 3)
        self.assertGreater(result.refinement_candidate_count, 0)
        self.assertEqual(
            result.evaluated_candidate_count,
            result.candidate_count + result.refinement_candidate_count,
        )
        self.assertNotIn(result.tau_deg, (0.0, 45.0, 90.0))
        self.assertAlmostEqual(result.tau_deg, target_tau, delta=0.001)
        self.assertLess(result.error, 1.0e-5)

    def test_multiple_solution_branches_are_refined_reported_and_sorted(self):
        target_distance = 2.0 * math.sin(math.radians(15.0))
        result = bend.screen_bend_angles(
            self.context(),
            {"phi": 0.0, "beta": 0.0},
            [bend.ScreenAngleRange("tau", -90.0, 90.0, 45.0)],
            bend.ScreeningRequest(
                mode="distance",
                target=target_distance,
                point1=self.atom("A:1:P"),
                point2=self.atom("C:1:P"),
            ),
            align_mode="n",
            solution_tolerance=0.001,
        )

        self.assertTrue(result.target_tolerance_met)
        self.assertEqual(result.refinement_region_count, 2)
        self.assertEqual(result.solution_count, 2)
        self.assertEqual(result.solutions[0].angles, result.angles)
        self.assertAlmostEqual(result.solutions[0].tau_deg, -30.0, delta=0.001)
        self.assertAlmostEqual(result.solutions[1].tau_deg, 30.0, delta=0.001)
        self.assertTrue(all(solution.error <= 0.001 for solution in result.solutions))
        table = bend.format_screening_solution_table(result, "A")
        self.assertIn("2 distinct solution(s) within tolerance", table)
        self.assertIn("phi (deg)", table)
        self.assertIn("residual", table)

    def test_solution_tolerance_falls_back_to_best_when_none_qualify(self):
        request = bend.ScreeningRequest(
            mode="distance",
            target=2.0 * math.sin(math.radians(15.0)),
            point1=self.atom("A:1:P"),
            point2=self.atom("C:1:P"),
        )
        result = bend.screen_bend_angles(
            self.context(),
            {"phi": 0.0, "beta": 0.0},
            [bend.ScreenAngleRange("tau", -90.0, 90.0, 45.0)],
            request,
            align_mode="n",
            solution_tolerance=1.0e-12,
        )

        self.assertFalse(result.target_tolerance_met)
        self.assertEqual(result.solution_count, 1)
        self.assertEqual(result.solutions[0].angles, result.angles)
        self.assertIn(
            "no solution met tolerance; closest fallback shown",
            bend.format_screening_solution_table(result, "A"),
        )

        with self.assertRaisesRegex(ValueError, "tolerance.*nonnegative"):
            bend.screen_bend_angles(
                self.context(),
                {"phi": 0.0, "beta": 0.0},
                [bend.ScreenAngleRange("tau", -90.0, 90.0, 45.0)],
                request,
                align_mode="n",
                solution_tolerance=-0.1,
            )

    def test_adaptive_refinement_adjusts_two_angles_between_coarse_steps(self):
        result = bend.screen_bend_angles(
            self.context(),
            {"tau": 0.0},
            [
                bend.ScreenAngleRange("phi", 0.0, 90.0, 45.0),
                bend.ScreenAngleRange("beta", 0.0, 90.0, 45.0),
            ],
            bend.ScreeningRequest(
                mode="distance",
                target=1.5,
                point1=self.atom("A:1:P"),
                point2=self.atom("C:1:P"),
            ),
            align_mode="y",
        )

        self.assertEqual(result.candidate_count, 9)
        self.assertGreater(result.refinement_candidate_count, 0)
        self.assertNotIn(result.phi_deg, (0.0, 45.0, 90.0))
        self.assertNotIn(result.beta_deg, (0.0, 45.0, 90.0))
        self.assertAlmostEqual(result.achieved_value, 1.5, delta=2.0e-5)

    def test_grid_preview_shows_values_counts_degrees_and_default_step(self):
        self.assertEqual(bend.DEFAULT_SCREEN_STEP_DEG, 6.0)
        self.assertEqual(bend.DEFAULT_SCREEN_STEP_A, 0.5)
        self.assertEqual(bend.DEFAULT_SCREEN_SOLUTION_TOLERANCE, 0.001)
        self.assertEqual(
            bend.DEFAULT_SCREEN_RANGES,
            {
                "phi": (-90.0, 90.0),
                "beta": (-180.0, 180.0),
                "tau": (-180.0, 180.0),
                "shift_axial": (-10.0, 10.0),
                "shift_radial": (-10.0, 10.0),
            },
        )
        self.assertEqual(
            bend.SCREEN_VARIABLE_NAMES,
            ("phi", "beta", "tau", "shift_axial", "shift_radial"),
        )
        for name in ("phi", "beta", "tau"):
            self.assertEqual(bend.screen_variable_unit(name), "deg")
            self.assertEqual(bend.screen_variable_default_step(name), 6.0)
            self.assertEqual(bend.screen_variable_tolerance(name), 0.001)
        for name in ("shift_axial", "shift_radial"):
            self.assertEqual(bend.screen_variable_unit(name), "A")
            self.assertEqual(bend.screen_variable_default_step(name), 0.5)
            self.assertEqual(bend.screen_variable_tolerance(name), 0.001)

        preview = bend.format_screening_grid_preview(
            [
                bend.ScreenAngleRange("phi", 0.0, 10.0, 4.0),
                bend.ScreenAngleRange("beta", -6.0, 6.0, 6.0),
            ]
        )

        self.assertIn("Phi grid (4): 0, 4, 8, 10 deg", preview)
        self.assertIn("Beta grid (3): -6, 0, 6 deg", preview)
        self.assertIn("Total coarse candidates: 12", preview)

        shift_preview = bend.format_screening_grid_preview(
            [
                bend.ScreenAngleRange("shift_axial", 0.0, 2.0, 1.0),
                bend.ScreenAngleRange("shift_radial", -1.0, 1.0, 1.0),
            ]
        )
        self.assertIn("Shift axial grid (3): 0, 1, 2 A", shift_preview)
        self.assertIn("Shift radial grid (3): -1, 0, 1 A", shift_preview)
        self.assertIn("Total coarse candidates: 9", shift_preview)

    def test_screening_popup_help_covers_every_argument(self):
        expected_help_keys = {
            "grid_from",
            "grid_to",
            "grid_step",
            "mode",
            "target",
            "solution_tolerance",
            "write_all_solutions",
            "endpoint1_atom",
            "endpoint2_source",
            "endpoint2_atom",
            "endpoint2_xyz",
            "endpoint2_pivot",
            "axis_source",
            "local_axis_ranges",
            "axis_point_source",
            "axis_point_xyz",
            "axis_point_atom",
            "axis_vector_source",
            "direct_vector",
            "two_xyz_points",
            "two_overlay_atoms",
            "normal_vectors",
        }
        launch_gui_source = inspect.getsource(bend.launch_gui)

        self.assertEqual(set(bend.SCREENING_GUI_HELP), expected_help_keys)
        for key in expected_help_keys:
            self.assertIn(f'"{key}"', launch_gui_source)
        self.assertIn('bg="#d9ecff"', launch_gui_source)
        # Column headers carry no unit; each row label supplies its own, because
        # the two pivot shifts are angstroms while the three angles are degrees.
        self.assertIn('(1, "From")', launch_gui_source)
        self.assertIn('(3, "To")', launch_gui_source)
        self.assertIn('(5, "Step")', launch_gui_source)
        self.assertIn("screen_variable_gui_label(name)", launch_gui_source)
        self.assertIn("screen_solution_tolerance_var", launch_gui_source)
        self.assertIn("screen_write_all_solutions_var", launch_gui_source)
        self.assertIn("screen_local_axis_ranges_var", launch_gui_source)
        step_help = bend.SCREENING_GUI_HELP["grid_step"]
        self.assertIn("searches between nearby grid values", step_help)
        self.assertIn("0.001 degree", step_help)
        self.assertIn("0.001 angstrom", step_help)

    def test_main_window_help_covers_the_pivot_shift_fields(self):
        launch_gui_source = inspect.getsource(bend.launch_gui)

        # Assert on the wiring, never on a literal grid row: pinning the row
        # number would fail whenever a row is inserted above these fields, yet
        # stay silent if the two "?" buttons were swapped onto each other's rows.
        # Comparing each field's own row against its help button's row catches
        # the misalignment without caring what the row number happens to be.
        for key in ("shift_axial", "shift_radial"):
            self.assertIn(f'"{key}": (', launch_gui_source)
            self.assertIn(f"{key}_var", launch_gui_source)
            self.assertIn(f'screen_select_vars["{key}"]', launch_gui_source)

            entry_row = re.search(
                rf"{key}_entry = ttk\.Entry\(root, textvariable={key}_var\)\s*\n"
                rf"\s*{key}_entry\.grid\(row=(\d+),",
                launch_gui_source,
            )
            self.assertIsNotNone(entry_row, f"{key} entry row not found")
            help_row = re.search(
                rf'help_button\((\d+), "{key}"\)', launch_gui_source
            )
            self.assertIsNotNone(help_row, f"{key} help button not found")
            self.assertEqual(
                entry_row.group(1),
                help_row.group(1),
                f"the {key} help button is not on the {key} row",
            )
        self.assertIn("variable_entries", launch_gui_source)

    def test_inclusive_grid_keeps_endpoints_when_step_exceeds_span(self):
        self.assertEqual(
            bend.inclusive_angle_values(
                bend.ScreenAngleRange("phi", 0.0, 1.0, 1.0e12)
            ),
            (0.0, 1.0),
        )
        self.assertEqual(
            bend.inclusive_angle_values(
                bend.ScreenAngleRange("phi", 1.0, 0.0, 1.0e12)
            ),
            (1.0, 0.0),
        )

    def test_tiny_inclusive_grid_has_unique_values_and_both_endpoints(self):
        values = bend.inclusive_angle_values(
            bend.ScreenAngleRange("phi", 0.0, 1.0e-11, 1.0e-12)
        )

        self.assertEqual(len(values), 11)
        self.assertEqual(len(set(values)), 11)
        self.assertEqual(values[0], 0.0)
        self.assertEqual(values[-1], 1.0e-11)
        self.assertTrue(all(first < second for first, second in zip(values, values[1:])))

    def test_signed_rotation_and_wrapped_error(self):
        axis_point = (0.0, 0.0, 0.0)
        axis_dir = (0.0, 0.0, 1.0)
        self.assertAlmostEqual(
            bend.signed_projected_angle_deg(
                (1.0, 0.0, 0.0),
                (0.0, 1.0, 0.0),
                axis_point,
                axis_dir,
            ),
            90.0,
        )
        self.assertAlmostEqual(
            bend.signed_projected_angle_deg(
                (0.0, 1.0, 0.0),
                (1.0, 0.0, 0.0),
                axis_point,
                axis_dir,
            ),
            -90.0,
        )
        self.assertAlmostEqual(bend.wrapped_angle_error_deg(179.0, -179.0), 2.0)
        self.assertAlmostEqual(bend.wrapped_angle_error_deg(-179.0, 179.0), 2.0)

    def test_ties_are_deterministic_and_screening_is_nonmutating(self):
        context = self.context()
        atoms_before = context.atoms
        file_before = self.input_path.read_bytes()
        request = bend.ScreeningRequest(
            mode="distance",
            target=math.sqrt(8.0),
            point1=self.atom("A:1:P"),
            point2=self.atom("C:1:P"),
        )
        args = (
            context,
            {"phi": 0.0, "beta": 0.0, "tau": 0.0},
            [bend.ScreenAngleRange("beta", -90.0, 90.0, 180.0)],
            request,
        )
        first = bend.screen_bend_angles(*args, align_mode="y")
        second = bend.screen_bend_angles(*args, align_mode="y")
        reverse_order = bend.screen_bend_angles(
            context,
            {"phi": 0.0, "beta": 0.0, "tau": 0.0},
            [bend.ScreenAngleRange("beta", 90.0, -90.0, 180.0)],
            request,
            align_mode="y",
        )

        self.assertAlmostEqual(first.beta_deg, -90.0)
        self.assertEqual(first, second)
        self.assertEqual(first, reverse_order)
        self.assertEqual(context.atoms, atoms_before)
        self.assertEqual(self.input_path.read_bytes(), file_before)

    def test_version_option_uses_centralized_version(self):
        completed = subprocess.run(
            [sys.executable, str(Path(bend.__file__).resolve()), "--version"],
            check=False,
            capture_output=True,
            text=True,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(
            completed.stdout.strip(),
            f"{bend.TOOL_NAME} {bend.VERSION}",
        )
        self.assertEqual(bend.VERSION, "V2.7")

    def test_screening_automatic_names_append_scr_before_optional_sep(self):
        self.assertEqual(
            bend.make_output_name(
                "model.pdb",
                phi_deg=0.0,
                beta_deg=30.0,
                tau_deg=0.0,
                sep_mode="n",
                screen_mode=True,
            ),
            "model_P0B30T0_scr.pdb",
        )
        self.assertEqual(
            bend.make_output_name(
                "model.pdb",
                phi_deg=0.0,
                beta_deg=30.0,
                tau_deg=0.0,
                sep_mode="y",
                screen_mode=True,
            ),
            "model_P0B30T0_scr_sep.pdb",
        )

        solution = bend.ScreeningSolution(10.0, -20.0, 30.0, 5.0, 0.0)
        self.assertEqual(
            bend.make_screen_solution_output_name(
                "model.pdb", solution, solution_index=2
            ),
            "model_P10Bm20T30_scr_sol002.pdb",
        )
        self.assertEqual(
            bend.make_screen_solution_output_name(
                "model.pdb", solution, solution_index=2, sep_mode="y"
            ),
            "model_P10Bm20T30_scr_sol002_sep.pdb",
        )
        self.assertEqual(
            bend.make_screen_solution_output_name(
                "model.pdb",
                solution,
                solution_index=2,
                explicit_primary_output="custom.pdb",
            ),
            "custom_sol002.pdb",
        )

    def test_screening_winner_writes_scr_output_and_origin_overlay(self):
        result = bend.screen_bend_angles(
            self.context(),
            {"phi": 0.0, "tau": 0.0},
            [bend.ScreenAngleRange("beta", 0.0, 90.0, 90.0)],
            bend.ScreeningRequest(
                mode="distance",
                target=math.sqrt(8.0),
                point1=self.atom("A:1:P"),
                point2=self.atom("C:1:P"),
            ),
            align_mode="y",
        )
        output_path = bend.make_output_name(
            str(self.input_path),
            result.phi_deg,
            result.beta_deg,
            result.tau_deg,
            screen_mode=True,
        )

        actual_path, info = bend.run_bending(
            input_pdb=str(self.input_path),
            pivot_residue="X2",
            phi_deg=result.phi_deg,
            beta_deg=result.beta_deg,
            tau_deg=result.tau_deg,
            origin_mode="y",
            output_pdb=output_path,
        )

        self.assertEqual(actual_path, output_path)
        self.assertTrue(actual_path.endswith("_P0B90T0_scr.pdb"))
        self.assertTrue(Path(actual_path).is_file())
        origin_path = Path(bend.make_origin_output_name(actual_path))
        self.assertEqual(info["origin_out_path"], str(origin_path))
        self.assertTrue(origin_path.is_file())
        overlay_chains = {
            line[21]
            for line in origin_path.read_text(encoding="utf-8").splitlines()
            if line.startswith(("ATOM  ", "HETATM"))
        }
        self.assertEqual(overlay_chains, {"A", "B", "C", "D"})

    def test_write_all_reported_solutions_writes_numbered_models_and_overlays(self):
        result = bend.screen_bend_angles(
            self.context(),
            {"phi": 0.0, "beta": 0.0},
            [bend.ScreenAngleRange("tau", -90.0, 90.0, 45.0)],
            bend.ScreeningRequest(
                mode="distance",
                target=2.0 * math.sin(math.radians(15.0)),
                point1=self.atom("A:1:P"),
                point2=self.atom("C:1:P"),
            ),
            align_mode="n",
            solution_tolerance=0.001,
        )
        primary_path = bend.make_output_name(
            str(self.input_path),
            result.phi_deg,
            result.beta_deg,
            result.tau_deg,
            screen_mode=True,
        )
        bend.run_bending(
            input_pdb=str(self.input_path),
            pivot_residue="X2",
            phi_deg=result.phi_deg,
            beta_deg=result.beta_deg,
            tau_deg=result.tau_deg,
            align_mode="n",
            origin_mode="y",
            output_pdb=primary_path,
        )

        additional = bend.write_additional_screening_solution_outputs(
            result=result,
            input_pdb=str(self.input_path),
            pivot_residue="X2",
            sep_mode="n",
            align_mode="n",
        )

        self.assertEqual(len(additional), 1)
        additional_path, origin_path = additional[0]
        self.assertTrue(additional_path.endswith("_scr_sol002.pdb"))
        self.assertEqual(origin_path, bend.make_origin_output_name(additional_path))
        self.assertTrue(Path(primary_path).is_file())
        self.assertTrue(Path(bend.make_origin_output_name(primary_path)).is_file())
        self.assertTrue(Path(additional_path).is_file())
        self.assertTrue(Path(origin_path).is_file())


class BendHelixPivotShiftTests(unittest.TestCase):
    """Geometry of the Sa/Sr pivot shifts on the ideal X/Y duplex.

    Fixture reference frame for pivot X2: axis_point (0,0,2), axis_dir (0,0,1),
    pivot P (1,0,2), axis_foot (0,0,2), radial_dir (1,0,0), radius 1.0, and
    movable piece #2 = {X2, X3, Y1, Y2}, so shift_vector == (Sr, 0, Sa).
    """

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.input_path = Path(self.temp_dir.name) / "xy_duplex.pdb"
        self.input_path.write_text(_screening_fixture_pdb(), encoding="utf-8")
        self.prep = bend.prepare_screening_context(
            str(self.input_path), "X2"
        ).preparation

    def tearDown(self):
        self.temp_dir.cleanup()

    def transform(self, phi, beta, tau, sa, sr, align="y"):
        return bend.build_bend_transform(
            self.prep, phi, beta, tau, sa, sr, align_mode=align
        )

    def assertPointAlmostEqual(self, actual, expected, places=7):
        for axis, (got, want) in enumerate(zip(actual, expected)):
            self.assertAlmostEqual(got, want, places=places, msg=f"axis {axis}")

    def test_fixture_reference_frame_is_what_the_shift_math_assumes(self):
        self.assertPointAlmostEqual(self.prep.axis_dir, (0.0, 0.0, 1.0))
        self.assertPointAlmostEqual(self.prep.selected_p, (1.0, 0.0, 2.0))
        self.assertPointAlmostEqual(self.prep.axis_foot, (0.0, 0.0, 2.0))
        self.assertPointAlmostEqual(self.prep.radial_dir, (1.0, 0.0, 0.0))
        self.assertAlmostEqual(self.prep.radius, 1.0)
        self.assertEqual(
            sorted(self.prep.piece2_keys),
            [("X", 2), ("X", 3), ("Y", 1), ("Y", 2)],
        )

    def test_zero_shift_reproduces_the_pre_shift_geometry_exactly(self):
        """A zero shift must be bit-identical, not merely close, to V2.6."""
        for phi, beta, tau in (
            (0.0, 0.0, 0.0),
            (37.5, -22.0, 13.0),
            (-180.0, 90.0, -45.0),
            (91.25, 179.0, 61.5),
        ):
            for align in ("y", "n"):
                transform = self.transform(phi, beta, tau, 0.0, 0.0, align)
                self.assertEqual(transform.shift_vector, (0.0, 0.0, 0.0))
                self.assertEqual(transform.pivot_shifted, self.prep.selected_p)
                self.assertEqual(transform.shifted_radius, self.prep.radius)

                # Inlined pre-shift reference implementation.
                beta_rad = math.radians(beta)
                tau_rad = math.radians(tau)
                radial_phi = bend.rotate_vector(
                    self.prep.radial, self.prep.axis_dir, math.radians(phi)
                )
                hinge_point = bend.v_add(self.prep.axis_foot, radial_phi)
                hinge_dir = bend.v_norm(
                    bend.v_cross(self.prep.axis_dir, radial_phi)
                )
                twist_point = bend.rotate_point_about_line(
                    self.prep.axis_foot, hinge_point, hinge_dir, beta_rad
                )
                twist_dir = bend.v_norm(
                    bend.rotate_vector(self.prep.axis_dir, hinge_dir, beta_rad)
                )

                def rotate_only(coord):
                    moved = bend.rotate_point_about_line(
                        coord, hinge_point, hinge_dir, beta_rad
                    )
                    if abs(tau_rad) > 0.0:
                        moved = bend.rotate_point_about_line(
                            moved, twist_point, twist_dir, tau_rad
                        )
                    return moved

                self.assertEqual(transform.hinge_point, hinge_point)
                self.assertEqual(transform.hinge_dir, hinge_dir)
                self.assertEqual(transform.twist_axis_point_pre_align, twist_point)
                self.assertEqual(transform.twist_axis_dir, twist_dir)

                for coord in ((1.0, 0.0, 4.0), (-7.5, 3.25, 12.0)):
                    expected = rotate_only(coord)
                    if align == "y":
                        expected = bend.v_add(
                            expected,
                            bend.v_sub(
                                self.prep.selected_p,
                                rotate_only(self.prep.selected_p),
                            ),
                        )
                    self.assertEqual(
                        transform.transform_coord(coord),
                        expected,
                        msg=f"phi={phi} beta={beta} tau={tau} align={align}",
                    )

    def test_shift_vector_is_axial_plus_radial_and_moves_only_piece_two(self):
        transform = self.transform(0.0, 0.0, 0.0, 3.0, 2.0)
        self.assertPointAlmostEqual(transform.shift_vector, (2.0, 0.0, 3.0))
        self.assertPointAlmostEqual(transform.pivot_shifted, (3.0, 0.0, 5.0))
        self.assertAlmostEqual(transform.shifted_radius, 3.0)
        # No rotation, so align y needs no correction at all.
        self.assertEqual(transform.align_translation, (0.0, 0.0, 0.0))
        self.assertPointAlmostEqual(
            transform.transform_coord((1.0, 0.0, 4.0)), (3.0, 0.0, 7.0)
        )
        self.assertPointAlmostEqual(
            transform.transform_coord((-1.0, 0.0, 2.0)), (1.0, 0.0, 5.0)
        )

        out_path, info = bend.run_bending(
            input_pdb=str(self.input_path),
            pivot_residue="X2",
            phi_deg=0.0,
            beta_deg=0.0,
            shift_axial=3.0,
            shift_radial=2.0,
        )
        coords = {}
        for line in Path(out_path).read_text(encoding="utf-8").splitlines():
            if line.startswith("ATOM  ") and line[12:16].strip() == "P":
                coords[(line[21], int(line[22:26]))] = (
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                )
        # Fixed piece #1 is untouched.
        self.assertPointAlmostEqual(coords[("X", 1)], (1.0, 0.0, 0.0), places=3)
        self.assertPointAlmostEqual(coords[("Y", 3)], (-1.0, 0.0, 0.0), places=3)
        # Movable piece #2 is displaced by exactly the shift vector.
        self.assertPointAlmostEqual(coords[("X", 3)], (3.0, 0.0, 7.0), places=3)
        self.assertPointAlmostEqual(coords[("Y", 2)], (1.0, 0.0, 5.0), places=3)
        self.assertPointAlmostEqual(info["shift_vector"], (2.0, 0.0, 3.0))

    def test_align_y_makes_the_shift_a_rigid_translation_of_the_bent_piece(self):
        """With --align y the shift never changes the shape of the bend.

        Aligning back to the *unshifted* pivot instead would cancel the shift
        outright, which is why the shifted pivot is the alignment target.
        """
        for phi, beta, tau in (
            (0.0, 30.0, 0.0),
            (90.0, -45.0, 20.0),
            (180.0, 90.0, 0.0),
            (-37.0, 12.5, -63.0),
        ):
            unshifted = self.transform(phi, beta, tau, 0.0, 0.0, "y")
            shifted = self.transform(phi, beta, tau, 3.0, 2.0, "y")
            for coord in ((1.0, 0.0, 4.0), (-1.0, 0.0, 2.0), (0.5, -2.0, 9.0)):
                self.assertPointAlmostEqual(
                    shifted.transform_coord(coord),
                    bend.v_add(
                        unshifted.transform_coord(coord), shifted.shift_vector
                    ),
                )

    def test_align_n_rebuilds_the_hinge_from_the_shifted_pivot(self):
        """With --align n the relocated hinge is observable, not a translation."""
        transform = self.transform(180.0, 90.0, 0.0, 0.0, 2.0, "n")
        self.assertAlmostEqual(transform.shifted_radius, 3.0)
        self.assertPointAlmostEqual(transform.hinge_point, (-3.0, 0.0, 2.0))
        self.assertPointAlmostEqual(
            transform.transform_coord((1.0, 0.0, 4.0)), (-5.0, 0.0, 8.0)
        )
        # A naive "translate the old hinge too" implementation gives (-1,0,4)
        # and a "leave the hinge alone" one gives (-3,0,6).
        unshifted = self.transform(180.0, 90.0, 0.0, 0.0, 0.0, "n")
        self.assertNotAlmostEqual(
            unshifted.transform_coord((1.0, 0.0, 4.0))[0] + 2.0, -5.0
        )

    def test_radial_shift_does_not_rotate_with_phi(self):
        """Sr follows the pivot's own radial direction, fixed at the phi=0 one."""
        transform = self.transform(180.0, 0.0, 0.0, 0.0, 2.0, "y")
        self.assertPointAlmostEqual(transform.shift_vector, (2.0, 0.0, 0.0))
        # Rotating the radial reference with phi would instead give (-1,0,4).
        self.assertPointAlmostEqual(
            transform.transform_coord((1.0, 0.0, 4.0)), (3.0, 0.0, 4.0)
        )

    def test_axial_shift_uses_the_prebend_axis_not_piece_twos_bent_axis(self):
        transform = self.transform(0.0, 90.0, 0.0, 3.0, 0.0, "y")
        # Sliding along piece #2's own bent axis would give (6,0,2) instead.
        self.assertPointAlmostEqual(
            transform.transform_coord((1.0, 0.0, 4.0)), (3.0, 0.0, 5.0)
        )

    def test_positive_axial_shift_opens_the_junction_for_either_pivot_partner(self):
        partner = bend.prepare_screening_context(
            str(self.input_path), "Y2"
        ).preparation
        self.assertEqual(partner.pair_idx, self.prep.pair_idx)
        self.assertEqual(
            sorted(partner.piece2_keys), sorted(self.prep.piece2_keys)
        )
        for preparation in (self.prep, partner):
            transform = bend.build_bend_transform(
                preparation, 0.0, 0.0, 0.0, 3.0, 0.0, align_mode="y"
            )
            self.assertPointAlmostEqual(transform.shift_vector, (0.0, 0.0, 3.0))

    def test_naming_the_partner_residue_reverses_the_radial_shift(self):
        partner = bend.prepare_screening_context(
            str(self.input_path), "Y2"
        ).preparation
        self.assertPointAlmostEqual(partner.radial_dir, (-1.0, 0.0, 0.0))
        transform = bend.build_bend_transform(
            partner, 0.0, 0.0, 0.0, 0.0, 2.0, align_mode="y"
        )
        self.assertPointAlmostEqual(
            transform.transform_coord((1.0, 0.0, 4.0)), (-1.0, 0.0, 4.0)
        )

    def test_tau_twists_about_piece_twos_own_translated_axis(self):
        transform = self.transform(0.0, 0.0, 90.0, 0.0, 2.0, "n")
        self.assertPointAlmostEqual(
            transform.twist_axis_point_pre_align, (2.0, 0.0, 2.0)
        )
        # Twisting about the unshifted axis foot would give (0,3,4).
        self.assertPointAlmostEqual(
            transform.transform_coord((1.0, 0.0, 4.0)), (2.0, 1.0, 4.0)
        )

    def test_hinge_direction_stays_continuous_through_and_past_the_axis(self):
        """Sr <= -radius is allowed and must not reverse the sense of beta."""
        previous = None
        for shift_radial in (0.5, 0.0, -0.5, -0.999, -1.0, -1.001, -2.0, -3.0):
            transform = self.transform(0.0, 90.0, 0.0, 0.0, shift_radial, "n")
            self.assertPointAlmostEqual(transform.hinge_dir, (0.0, 1.0, 0.0))
            self.assertAlmostEqual(transform.shifted_radius, 1.0 + shift_radial)
            self.assertPointAlmostEqual(
                transform.hinge_point, (1.0 + shift_radial, 0.0, 2.0)
            )
            moved = transform.transform_coord((0.0, 0.0, 12.0))
            self.assertPointAlmostEqual(moved, (11.0 + shift_radial, 0.0, 3.0))
            if previous is not None:
                self.assertLess(moved[0], previous)
            previous = moved[0]

    def test_hinge_on_the_axis_still_lets_phi_orient_the_bend(self):
        transform = self.transform(90.0, 90.0, 0.0, 0.0, -1.0, "n")
        self.assertAlmostEqual(transform.shifted_radius, 0.0)
        self.assertPointAlmostEqual(transform.hinge_point, (0.0, 0.0, 2.0))
        self.assertPointAlmostEqual(transform.hinge_dir, (-1.0, 0.0, 0.0))
        self.assertPointAlmostEqual(
            transform.transform_coord((1.0, 0.0, 4.0)), (0.0, 2.0, 2.0)
        )

    def test_non_finite_shifts_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "Axial pivot shift must be finite"):
            self.transform(0.0, 0.0, 0.0, float("inf"), 0.0)
        with self.assertRaisesRegex(ValueError, "Radial pivot shift must be finite"):
            self.transform(0.0, 0.0, 0.0, 0.0, float("nan"))

    def test_origin_overlay_applies_the_same_shift_as_the_main_output(self):
        out_path, info = bend.run_bending(
            input_pdb=str(self.input_path),
            pivot_residue="X2",
            phi_deg=25.0,
            beta_deg=40.0,
            tau_deg=15.0,
            shift_axial=2.0,
            shift_radial=-1.5,
            align_mode="n",
            origin_mode="y",
        )
        main_coords = _p_coords(Path(out_path))
        overlay_coords = _p_coords(Path(info["origin_out_path"]))
        model1 = info["origin_chain_map_model1"]
        model2 = info["origin_chain_map_model2"]
        source_coords = _p_coords(self.input_path)

        for (chain, res_seq), original in source_coords.items():
            # Overlay model 1 is the untouched input.
            self.assertPointAlmostEqual(
                overlay_coords[(model1[chain], res_seq)], original, places=3
            )
            # Overlay model 2 is the same rigid transform piece #2 received, so
            # for piece-2 residues it must equal the main bent output exactly.
            transformed = overlay_coords[(model2[chain], res_seq)]
            if (chain, res_seq) in self.prep.piece2_keys:
                self.assertPointAlmostEqual(
                    transformed, main_coords[(chain, res_seq)], places=3
                )
            else:
                # Piece-1 residues stay put in the main output but are still
                # transformed in the overlay, so they must differ by the shift.
                self.assertNotAlmostEqual(
                    transformed[2], main_coords[(chain, res_seq)][2], places=3
                )

    def test_output_name_appends_the_shift_block_only_when_nonzero(self):
        self.assertEqual(
            bend.make_output_name("model.pdb", 0.0, 30.0, 0.0),
            "model_P0B30T0.pdb",
        )
        self.assertEqual(
            bend.make_output_name("model.pdb", 0.0, 30.0, 0.0, 0.0, 0.0),
            "model_P0B30T0.pdb",
        )
        self.assertEqual(
            bend.make_output_name("model.pdb", 0.0, 30.0, 0.0, 2.0, -1.5),
            "model_P0B30T0Sa2Srm1p5.pdb",
        )
        self.assertEqual(
            bend.make_output_name(
                "model.pdb", 0.0, 30.0, 0.0, 0.0, -1.5, sep_mode="y", screen_mode=True
            ),
            "model_P0B30T0Sa0Srm1p5_scr_sep.pdb",
        )
        solution = bend.ScreeningSolution(
            10.0, -20.0, 30.0, 5.0, 0.0, shift_axial=1.25, shift_radial=0.0
        )
        self.assertEqual(
            bend.make_screen_solution_output_name(
                "model.pdb", solution, solution_index=2
            ),
            "model_P10Bm20T30Sa1p25Sr0_scr_sol002.pdb",
        )

    def test_equivalent_cli_command_reports_both_shifts(self):
        command = bend.build_equivalent_cli_command(
            input_pdb="in.pdb",
            pivot_residue="A36",
            phi_deg=0.0,
            beta_deg=30.0,
            tau_deg=0.0,
            sep_mode="n",
            align_mode="y",
            origin_mode="n",
            shift_axial=2.0,
            shift_radial=-1.5,
        )
        self.assertIn("--shift_axial 2", command)
        self.assertIn("--shift_radial -1.5", command)

    def test_cli_accepts_and_applies_both_shifts(self):
        completed = subprocess.run(
            [
                sys.executable,
                str(Path(bend.__file__).resolve()),
                "--input",
                str(self.input_path),
                "--pivot",
                "X2",
                "--phi",
                "0",
                "--beta",
                "0",
                "--shift_axial",
                "3",
                "--shift_radial",
                "2",
            ],
            check=False,
            capture_output=True,
            text=True,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("axial 3.000000 A, radial 2.000000 A", completed.stdout)
        written = self.input_path.with_name("xy_duplex_P0B0T0Sa3Sr2.pdb")
        self.assertTrue(written.is_file(), completed.stdout)
        self.assertPointAlmostEqual(
            _p_coords(written)[("X", 3)], (3.0, 0.0, 7.0), places=3
        )


class BendHelixShiftScreeningTests(unittest.TestCase):
    """The two pivot shifts as screened variables alongside phi/beta/tau."""

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.input_path = Path(self.temp_dir.name) / "xy_duplex.pdb"
        self.input_path.write_text(_screening_fixture_pdb(), encoding="utf-8")

    def tearDown(self):
        self.temp_dir.cleanup()

    def context(self):
        return bend.prepare_screening_context(str(self.input_path), "X2")

    @staticmethod
    def atom(selector):
        return bend.ScreeningPoint("overlay_atom", selector)

    def test_variable_names_and_aliases_normalize(self):
        for raw, expected in (
            ("phi", "phi"),
            ("Beta_deg", "beta"),
            ("shift_axial", "shift_axial"),
            ("shift-radial", "shift_radial"),
            ("Sa", "shift_axial"),
            ("sr", "shift_radial"),
            ("axial_shift", "shift_axial"),
            ("shift_axial_a", "shift_axial"),
        ):
            self.assertEqual(bend._normalize_screen_variable_name(raw), expected)
        with self.assertRaisesRegex(ValueError, "Unknown screened variable"):
            bend._normalize_screen_variable_name("gamma")

    def test_unscreened_shifts_default_to_zero_for_angle_only_callers(self):
        """Pre-V2.7 callers pass only phi/beta/tau and must keep working."""
        result = bend.screen_bend_angles(
            self.context(),
            {"phi": 0.0, "beta": 0.0},
            [bend.ScreenAngleRange("tau", 0.0, 90.0, 45.0)],
            bend.ScreeningRequest(
                mode="distance",
                target=2.0 * math.sin(math.radians(15.0)),
                point1=self.atom("A:1:P"),
                point2=self.atom("C:1:P"),
            ),
            align_mode="n",
        )
        self.assertEqual(result.shift_axial, 0.0)
        self.assertEqual(result.shift_radial, 0.0)
        self.assertEqual(result.angles, {"phi": 0.0, "beta": 0.0, "tau": result.tau_deg})
        self.assertEqual(
            result.variables,
            {
                "phi": 0.0,
                "beta": 0.0,
                "tau": result.tau_deg,
                "shift_axial": 0.0,
                "shift_radial": 0.0,
            },
        )

    def test_a_nonzero_unscreened_shift_is_held_not_discarded(self):
        """Screening one variable must carry a fixed shift through untouched.

        The GUI reaches this whenever a shift is typed in but not checked for
        screening, so a fixed shift silently reset to 0 would bend the wrong
        structure while still reporting a matched target.
        """
        result = bend.screen_bend_angles(
            self.context(),
            {"phi": 0.0, "beta": 0.0, "shift_axial": 4.0, "shift_radial": 0.0},
            [bend.ScreenAngleRange("tau", 0.0, 90.0, 45.0)],
            bend.ScreeningRequest(
                mode="distance",
                target=4.0,
                point1=self.atom("A:2:P"),
                point2=self.atom("C:2:P"),
            ),
            align_mode="y",
        )
        # A:2:P is the untouched original pivot and C:2:P the shifted one, so
        # the separation is the axial shift itself and tau cannot change it.
        self.assertEqual(result.shift_axial, 4.0)
        self.assertEqual(result.shift_radial, 0.0)
        self.assertAlmostEqual(result.achieved_value, 4.0, places=7)
        self.assertEqual(result.variables["shift_axial"], 4.0)
        for solution in result.solutions:
            self.assertEqual(solution.shift_axial, 4.0)

    def test_screening_an_axial_shift_finds_a_distance_between_coarse_steps(self):
        """Sa is the only free variable, so the answer is exact and known."""
        result = bend.screen_bend_angles(
            self.context(),
            {"phi": 0.0, "beta": 0.0, "tau": 0.0, "shift_radial": 0.0},
            [bend.ScreenAngleRange("shift_axial", 0.0, 6.0, 2.0)],
            bend.ScreeningRequest(
                mode="distance",
                target=3.7,
                point1=self.atom("A:2:P"),
                point2=self.atom("C:2:P"),
            ),
            align_mode="y",
        )
        # A:2:P is the fixed original pivot; C:2:P is the shifted one, so the
        # separation is exactly |Sa| along the axis.
        self.assertAlmostEqual(result.shift_axial, 3.7, delta=1.0e-3)
        self.assertAlmostEqual(result.achieved_value, 3.7, delta=1.0e-3)
        self.assertGreater(result.refinement_candidate_count, 0)
        self.assertTrue(result.target_tolerance_met)

    def test_screening_a_radial_shift_reports_both_signed_branches(self):
        result = bend.screen_bend_angles(
            self.context(),
            {"phi": 0.0, "beta": 0.0, "tau": 0.0, "shift_axial": 0.0},
            [bend.ScreenAngleRange("shift_radial", -4.0, 4.0, 2.0)],
            bend.ScreeningRequest(
                mode="distance",
                target=2.5,
                point1=self.atom("A:2:P"),
                point2=self.atom("C:2:P"),
            ),
            align_mode="y",
        )
        self.assertTrue(result.target_tolerance_met)
        self.assertEqual(result.solution_count, 2)
        found = sorted(round(s.shift_radial, 3) for s in result.solutions)
        self.assertAlmostEqual(found[0], -2.5, delta=1.0e-3)
        self.assertAlmostEqual(found[1], 2.5, delta=1.0e-3)
        table = bend.format_screening_solution_table(result, "A")
        self.assertIn("shift axial (A)", table)
        self.assertIn("shift radial (A)", table)

    def test_shift_and_angle_screen_together_and_cache_on_all_variables(self):
        """A stale 3-angle cache key would collapse every Sa candidate into one."""
        result = bend.screen_bend_angles(
            self.context(),
            {"phi": 0.0, "tau": 0.0, "shift_radial": 0.0},
            [
                bend.ScreenAngleRange("beta", 0.0, 90.0, 45.0),
                bend.ScreenAngleRange("shift_axial", 0.0, 4.0, 2.0),
            ],
            bend.ScreeningRequest(
                mode="distance",
                target=5.0,
                point1=self.atom("A:1:P"),
                point2=self.atom("C:3:P"),
            ),
            align_mode="y",
        )
        self.assertEqual(result.candidate_count, 9)
        self.assertGreater(result.refinement_candidate_count, 0)
        self.assertLess(result.error, 1.0e-3)
        # Rebuilding the winning candidate must reproduce the reported metric.
        transform = bend.build_bend_transform(
            self.context().preparation,
            result.phi_deg,
            result.beta_deg,
            result.tau_deg,
            result.shift_axial,
            result.shift_radial,
            align_mode="y",
        )
        self.assertAlmostEqual(
            bend.distance_between_points(
                (1.0, 0.0, 0.0), transform.transform_coord((1.0, 0.0, 4.0))
            ),
            result.achieved_value,
            places=9,
        )

    def test_screened_shift_solutions_write_models_with_shifted_names(self):
        result = bend.screen_bend_angles(
            self.context(),
            {"phi": 0.0, "beta": 0.0, "tau": 0.0, "shift_axial": 0.0},
            [bend.ScreenAngleRange("shift_radial", -4.0, 4.0, 2.0)],
            bend.ScreeningRequest(
                mode="distance",
                target=2.5,
                point1=self.atom("A:2:P"),
                point2=self.atom("C:2:P"),
            ),
            align_mode="y",
        )
        self.assertEqual(result.solution_count, 2)
        primary = bend.make_output_name(
            str(self.input_path),
            result.phi_deg,
            result.beta_deg,
            result.tau_deg,
            result.shift_axial,
            result.shift_radial,
            screen_mode=True,
        )
        self.assertIn("Sa0Sr", Path(primary).name)
        bend.run_bending(
            input_pdb=str(self.input_path),
            pivot_residue="X2",
            phi_deg=result.phi_deg,
            beta_deg=result.beta_deg,
            tau_deg=result.tau_deg,
            shift_axial=result.shift_axial,
            shift_radial=result.shift_radial,
            origin_mode="y",
            output_pdb=primary,
        )
        additional = bend.write_additional_screening_solution_outputs(
            result=result,
            input_pdb=str(self.input_path),
            pivot_residue="X2",
            sep_mode="n",
            align_mode="y",
        )
        self.assertEqual(len(additional), 1)
        additional_path, origin_path = additional[0]
        self.assertIn("Sa0Sr", Path(additional_path).name)
        self.assertTrue(additional_path.endswith("_scr_sol002.pdb"))
        for path in (primary, additional_path, origin_path):
            self.assertTrue(Path(path).is_file(), path)
        # The extra model must carry the other branch's shift, not the winner's.
        other = [
            solution
            for solution in result.solutions[1:]
        ][0]
        self.assertAlmostEqual(
            _p_coords(Path(additional_path))[("X", 2)][0],
            1.0 + other.shift_radial,
            places=3,
        )

    def test_still_refuses_more_than_two_screened_variables(self):
        with self.assertRaisesRegex(ValueError, "one or two"):
            bend.validate_screen_angle_ranges(
                [
                    bend.ScreenAngleRange("phi", 0.0, 1.0, 1.0),
                    bend.ScreenAngleRange("shift_axial", 0.0, 1.0, 1.0),
                    bend.ScreenAngleRange("shift_radial", 0.0, 1.0, 1.0),
                ],
                candidate_cap=100,
            )
        with self.assertRaisesRegex(ValueError, "more than once"):
            bend.validate_screen_angle_ranges(
                [
                    bend.ScreenAngleRange("shift_axial", 0.0, 1.0, 1.0),
                    bend.ScreenAngleRange("sa", 0.0, 1.0, 1.0),
                ],
                candidate_cap=100,
            )


def _p_coords(path: Path):
    """Return {(chain, res_seq): xyz} for every P atom in a written PDB."""
    coords = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(("ATOM  ", "HETATM")) and line[12:16].strip() == "P":
            coords[(line[21], int(line[22:26]))] = (
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54]),
            )
    return coords


if __name__ == "__main__":
    unittest.main()
