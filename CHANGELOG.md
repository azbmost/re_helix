# Changelog

## V3.7 - 2026-06-16

- Bumped `re_helix` to V3.7.
- Bundled Bend Helix V2.4 as `re_helix_lib/bend_helix.py`.
- Added a `Bend Helix` button to the `re_helix` GUI. The button opens the Bend Helix GUI in a separate process and pre-fills the current input PDB when available.
- Moved tool launchers into a dedicated bold-titled `Other tools` GUI area for future expansion.
- Used bold styling for the main GUI title and section titles.
- Updated the README with Bend Helix GUI and CLI usage.

## V3.6 - 2026-06-16

- Renamed the main entry script from `align_re_helicesV3_5.py` to `re_helix.py`.
- Added `--re_only` / `--re-only` mode to apply reciprocal exchanges directly to the input PDB without running helix alignment.
- Reused the same reciprocal-exchange writer for RE-only and aligned+RE outputs, including REMARK 950 metadata, LINK records, and X33 linker records.
- Moved helper scripts into `re_helix_lib/` and made the helper folder importable as a Python package.
- Moved the optional GUI icon to `assets/icon.png`; the GUI loads it when present and runs normally when it is absent.
- Updated the Tk GUI title and controls for `re_helix V3.6`, including an RE-only checkbox.
- Added AZBMOST Package Module #2 wording to the GUI title/header and README.
- Expanded the README explanation of reciprocal exchange as a virtual DNA-nanostructure topology operation in the Seeman design tradition.
- Added `-v` / `--version` to print the app version from the command line.
- Added GitHub-ready project files: README, MIT license, changelog, and `.gitignore`.

## V3.5 - 2026-06-14

- Added optional fixed rho-angle syntax for single-site inter-helix reciprocal exchanges.
- Supported `<pos1> <pos2> <rho_deg> <kind>` for a single reciprocal-exchange site when `--axis_parallel n` is selected.
- Ignored fixed rho definitions with warnings for multi-site helix pairs or `--axis_parallel y`.

## V3.4 - 2026-05-23

- Loaded the V3.3 reciprocal-exchange helper in memory.
- Wrote parse-friendly REMARK 950 RE_SCRIPT metadata into output PDB files.
- Included X33 HET/HETNAM records for standalone 3'-3' linker phosphates.

## V3.3 - 2026-03-29

- Improved GUI row-count controls for pair rows and axis-range rows.
- Added a vertically scrollable GUI container for large inputs.

## V3.2 - 2026-03-29

- Ignored blank or incomplete GUI pair rows and axis-range rows.
- Improved GUI completion messages.

## V3.1 - 2026-03-29

- Added GUI help buttons for core fields and option groups.
- Added unit labels in the GUI.
- Changed replication handling for `--axis_range` definitions so final post-replication chain IDs can be targeted explicitly.

## V3.0 - 2026-03-28

- Added manual helical-axis residue-range definitions with repeatable `--axis_range` options.
- Added optional Tk GUI mode.
