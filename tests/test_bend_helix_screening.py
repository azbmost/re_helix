import inspect
import math
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
        self.assertEqual(bend.DEFAULT_SCREEN_SOLUTION_TOLERANCE, 0.001)
        self.assertEqual(
            bend.DEFAULT_SCREEN_RANGES,
            {
                "phi": (-90.0, 90.0),
                "beta": (-180.0, 180.0),
                "tau": (-180.0, 180.0),
            },
        )
        preview = bend.format_screening_grid_preview(
            [
                bend.ScreenAngleRange("phi", 0.0, 10.0, 4.0),
                bend.ScreenAngleRange("beta", -6.0, 6.0, 6.0),
            ]
        )

        self.assertIn("Phi grid (4): 0, 4, 8, 10 deg", preview)
        self.assertIn("Beta grid (3): -6, 0, 6 deg", preview)
        self.assertIn("Total coarse candidates: 12", preview)

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
        self.assertIn('"From (deg)"', launch_gui_source)
        self.assertIn('"To (deg)"', launch_gui_source)
        self.assertIn('"Step (deg)"', launch_gui_source)
        self.assertIn("screen_solution_tolerance_var", launch_gui_source)
        self.assertIn("screen_write_all_solutions_var", launch_gui_source)
        self.assertIn("screen_local_axis_ranges_var", launch_gui_source)
        step_help = bend.SCREENING_GUI_HELP["grid_step"]
        self.assertIn("searches between nearby grid values", step_help)
        self.assertIn("0.001 degree", step_help)

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
        self.assertEqual(bend.VERSION, "V2.6")

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


if __name__ == "__main__":
    unittest.main()
