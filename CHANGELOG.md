# Changelog

## V3.16 - 2026-06-26

- Bumped `re_helix` to V3.16.
- Bundled Get Phenix Restraints V1.0 as `re_helix_lib/get_phenix_restraints.py`, based on `link_to_geometry_restraintsV4.py`.
- Added a `Get Phenix Restraints` button to the `re_helix` GUI `Other tools` area. The button opens the bundled GUI in a separate process and pre-fills the current input PDB when available.
- Added GUI/CLI controls for `link_distance_cutoff`, standalone 3'-to-3' linker residue name, optional `*_junctions.params` generation, and linker CIF/safe-interpretation support file generation.
- Updated the README with the Phenix minimization workflow and the warning to use exactly one movement-selection params file.

## V3.15 - 2026-06-25

- Bumped `re_helix` to V3.15.
- Updated the GUI so the default `Output base` follows input PDB changes unless the user has entered a custom output base.
- Added a GUI `CLI pair args` field below the individual Exchange pairs rows. When filled, it replaces the row widgets and lets large exchange specs be pasted as one command-line-style token string.
- Fixed the bundled Add PDB LINK Record tool so topology rebuilds preserve existing input `LINK` records, remapping their endpoints to the rebuilt chain/residue labels when possible.

## V3.14 - 2026-06-25

- Bumped `re_helix` to V3.14.
- Added configurable phosphate-only residue output for bowtie 3'-3' linkages: default `HETATM X33`, custom `HETATM <name>`, or regular `ATOM DA` for Phenix-friendly relaxation without a custom residue definition.
- Added CLI and GUI controls for the 3'-3' linker phosphate residue name.

## V3.13 - 2026-06-24

- Bumped `re_helix` to V3.13.
- Internal package maintenance for bundled library utilities.

## V3.12 - 2026-06-23

- Bumped `re_helix` to V3.12.
- Added the bundled Generate Lattice V3.1 tool as `re_helix_lib/generate_lattice.py`.
- Added a `Generate Lattice` button to the `re_helix` GUI `Other tools` area. The button opens the bundled GUI in a separate process and pre-fills the current input PDB when available.
- Increased the Generate Lattice GUI window height and minimum size so all controls fit more comfortably.
- Updated the bundled tool title, `--gui` prefill handling, `-v` / `--version` output, and stdout run summaries for the main `re_helix` log.
- Updated the README with Generate Lattice GUI and CLI usage.

## V3.11 - 2026-06-20

- Bumped `re_helix` to V3.11.
- Added the bundled Insert Virtual Resi V1.0 tool as `re_helix_lib/insert_virtual_resi.py`.
- Added an `Insert Virtual Resi` button to the `re_helix` GUI `Other tools` area. The button opens the bundled GUI in a separate process and pre-fills the current input PDB when available.
- Supported insertion specs such as `A55 3`, `A.55 3`, `55A 3`, and `55.A 3`, with repeated specs for multiple residue-numbering gaps.
- Updated coordinate-like records, `TER` records, and both residue endpoints of fixed-column `LINK` records using the same virtual-residue renumbering map.
- Updated the README with Insert Virtual Resi GUI and CLI usage.

## V3.10 - 2026-06-18

- Bumped `re_helix` to V3.10.
- Mirrored stdout/stderr from bundled `Other tools` processes into the main `re_helix` run log after launch.
- Updated Bend Helix and Do Symmetry GUI runs to print equivalent CLI commands and completion/error summaries for the main log.
- Updated Add PDB LINK Record GUI runs to print related CLI information when applicable, selected LINK records, output paths, and errors for the main log.
- Updated the README with the main-log mirroring behavior.

## V3.9 - 2026-06-18

- Bumped `re_helix` to V3.9.
- Bundled Add PDB LINK Record V1.0 as `re_helix_lib/add_pdb_link_record.py`.
- Added an `Add PDB LINK Record` button to the `re_helix` GUI `Other tools` area. The button opens the bundled GUI in a separate process and pre-fills the current input PDB when available.
- Updated the bundled tool title, CLI program name, and `-v` / `--version` output for package-friendly presentation.
- Updated the README with Add PDB LINK Record GUI and CLI usage.

## V3.8 - 2026-06-18

- Bumped `re_helix` to V3.8.
- Bundled Do Symmetry V1.0 as `re_helix_lib/do_symmetry.py`.
- Added a `Do Symmetry` button to the `re_helix` GUI `Other tools` area. The button opens the Do Symmetry GUI in a separate process and pre-fills the current input PDB when available.
- Improved the bundled Do Symmetry GUI with a bold title and bold option-section headings.
- Updated the README with Do Symmetry GUI/CLI usage and its symmetry-averaging working principle.

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
