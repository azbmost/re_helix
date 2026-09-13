from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from re_helix_lib import insert_virtual_resi as ivr


ROOT = Path(__file__).resolve().parents[1]

HEADERED_PDB = """\
HEADER    DNA NANOSTRUCTURE
REMARK   1 UPSTREAM STEP
LINK         O3'  DA A  55                 P    DT A  56     1555   1555  1.60
ATOM      1  P    DA A  55      10.000  10.000  10.000  1.00  0.00           P
ATOM      2  P    DT A  56      12.000  10.000  10.000  1.00  0.00           P
ATOM      3  P    DG A  70      13.000  10.000  10.000  1.00  0.00           P
ATOM      4  P    DC A  71      14.000  10.000  10.000  1.00  0.00           P
TER       5       DC A  71
END
"""


def _specs(*pairs):
    return [ivr.parse_insertion_spec(token, str(count)) for token, count in pairs]


def _run(text, specs, **kwargs):
    with tempfile.TemporaryDirectory() as tmp:
        src = Path(tmp) / "in.pdb"
        src.write_text(text)
        out, stats = ivr.insert_virtual_residues(
            src, specs, output_pdb=Path(tmp) / "out.pdb", verbose=False, **kwargs
        )
        return out.read_text(), stats


def _remark_lines(text):
    return [line for line in text.splitlines() if line.startswith(ivr.REMARK_PREFIX)]


class VirtualRangeTests(unittest.TestCase):
    def test_single_spec_range_follows_anchor(self):
        (item,) = ivr.compute_virtual_ranges(_specs(("A55", 3)))
        self.assertEqual((item.after_new, item.start, item.end), (55, 56, 58))

    def test_later_spec_in_chain_is_pushed_by_earlier_gap(self):
        first, second = ivr.compute_virtual_ranges(_specs(("A55", 3), ("A70", 2)))
        self.assertEqual((first.start, first.end), (56, 58))
        self.assertEqual((second.after_new, second.start, second.end), (73, 74, 75))

    def test_specs_sharing_an_anchor_get_adjacent_spans(self):
        first, second = ivr.compute_virtual_ranges(_specs(("A55", 3), ("A55", 2)))
        self.assertEqual((first.start, first.end), (56, 58))
        self.assertEqual((second.start, second.end), (59, 60))

    def test_range_abuts_the_shifted_next_residue(self):
        specs = _specs(("A55", 3), ("A70", 2), ("B10", 1))
        shift_map = ivr.build_shift_map(specs)
        for item in ivr.compute_virtual_ranges(specs):
            chain = item.spec.chain_id
            anchor = item.spec.after_resseq
            self.assertEqual(item.after_new, ivr.shifted_resseq(chain, anchor, shift_map))
            self.assertEqual(item.end + 1, ivr.shifted_resseq(chain, anchor + 1, shift_map))
            self.assertEqual(item.end - item.start + 1, item.spec.count)


class RemarkWritingTests(unittest.TestCase):
    def test_virtual_insert_remark_reports_chain_anchor_count_and_range(self):
        text, stats = _run(HEADERED_PDB, _specs(("A55", 3), ("A70", 2)))
        inserts = [line for line in _remark_lines(text) if "VIRTUAL_INSERT" in line]
        self.assertEqual(len(inserts), 2)
        self.assertIn(
            "VIRTUAL_INSERT op=1 chain=A after_orig=A:55 after_new=A:55 "
            "count=3 start=A:56 end=A:58",
            inserts[0],
        )
        self.assertIn(
            "VIRTUAL_INSERT op=2 chain=A after_orig=A:70 after_new=A:73 "
            "count=2 start=A:74 end=A:75",
            inserts[1],
        )
        self.assertEqual(stats.remark_lines_written, len(_remark_lines(text)))

    def test_header_records_software_and_command(self):
        text, _ = _run(HEADERED_PDB, _specs(("A55", 3)), command="demo command")
        remarks = _remark_lines(text)
        self.assertIn(f"{ivr.REMARK_PREFIX} SOFTWARE name=insert_virtual_resi", remarks[0])
        self.assertIn(f"{ivr.REMARK_PREFIX} COMMAND text=demo command", remarks[1])
        self.assertIn(f"{ivr.REMARK_PREFIX} OUTPUT_STAGE name=insert_virtual_resi", remarks[2])

    def test_remarks_land_after_existing_remarks_and_before_link(self):
        text, _ = _run(HEADERED_PDB, _specs(("A55", 3)))
        records = [line[:6].strip() for line in text.splitlines()]
        first_new = next(
            i for i, line in enumerate(text.splitlines()) if line.startswith(ivr.REMARK_PREFIX)
        )
        self.assertEqual(records[first_new - 1], "REMARK")
        self.assertLess(first_new, records.index("LINK"))
        self.assertNotIn("ATOM", records[:first_new])

    def test_remarks_precede_coordinates_when_no_header_remark_exists(self):
        text, _ = _run(
            "ATOM      1  P    DA A  56      1.000  1.000  1.000  1.00  0.00           P\n",
            _specs(("A55", 3)),
        )
        self.assertTrue(text.splitlines()[0].startswith(ivr.REMARK_PREFIX))

    def test_no_remark_option_leaves_the_file_untouched_apart_from_numbering(self):
        with_remarks, _ = _run(HEADERED_PDB, _specs(("A55", 3)))
        without, stats = _run(HEADERED_PDB, _specs(("A55", 3)), write_remarks=False)
        self.assertEqual(_remark_lines(without), [])
        self.assertEqual(stats.remark_lines_written, 0)
        self.assertEqual(
            [line for line in with_remarks.splitlines() if not line.startswith(ivr.REMARK_PREFIX)],
            without.splitlines(),
        )

    def test_renumbering_is_unchanged_by_the_remark_block(self):
        text, _ = _run(HEADERED_PDB, _specs(("A55", 3), ("A70", 2)))
        body = [line for line in text.splitlines() if line[:6].strip() in {"ATOM", "TER", "LINK"}]
        self.assertIn("P    DT A  59", body[2])
        self.assertIn("P    DG A  73", body[3])
        self.assertIn("P    DC A  76", body[4])
        self.assertIn("P    DT A  59", body[0])

    def test_summary_reports_the_virtual_ranges(self):
        _, stats = _run(HEADERED_PDB, _specs(("A55", 3), ("A70", 2)))
        summary = ivr.format_summary(Path("out.pdb"), stats)
        self.assertIn("Virtual residue ranges (output numbering): A56-A58, A74-A75", summary)

    def test_cli_command_round_trips_the_no_remark_flag(self):
        specs = _specs(("A55", 3))
        self.assertNotIn("--no-remark", ivr.build_cli_command("t.py", "in.pdb", specs, "out.pdb"))
        self.assertIn(
            "--no-remark", ivr.build_cli_command("t.py", "in.pdb", specs, "out.pdb", False)
        )


if __name__ == "__main__":
    unittest.main()
