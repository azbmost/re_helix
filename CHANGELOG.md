# Changelog

## V4.5 - 2026-08-27

- Bumped `re_helix` to V4.5.
- Added **Normal to two vectors** as a fourth vector source for both Restrained translation and Restrained rotation.
- Defined the normal using the normalized right-hand cross product `vector 1 × vector 2`, with validation that rejects non-finite, zero, or parallel vector pairs.
- Added conditional GUI fields, CLI options, help text, documentation, and regression coverage for both restrained modes.

## V4.4 - 2026-08-27

- Bumped `re_helix` to V4.4.
- Added the **Restrained rotation** alignment mode. It leaves the fixed helix unchanged, applies no translation, and optimizes only one rotation angle about a supplied point-and-vector axis using all usable alignment pairs.
- Allowed the rotation-axis point to come from XYZ coordinates or one input-PDB atom.
- Allowed the rotation-axis vector to come from a direct vector, two XYZ points, or two input-PDB atoms.
- Added conditional GUI fields for the rotation point and vector sources while retaining Axis definition grouping and `move with axis` payload behavior.
- Added CLI validation, run-log reporting, documentation, and regression coverage confirming rotation without translation.

## V4.3 - 2026-08-27

- Bumped `re_helix` to V4.3.
- Added the **Restrained translation** alignment mode. It leaves the fixed helix unchanged, applies no rotation to an unfixed helix, and chooses the least-squares translation constrained to the supplied direction using all usable alignment pairs.
- Added direction input by a vector, two XYZ points, or two atoms in the input PDB. Atom selectors accept serial syntax such as `#123` and chain/residue/name syntax such as `A:30:P`, `A30:P`, or `30A:P`.
- Added the alignment-mode and direction-source fields between **Axis definitions** and **Other tools** in the GUI, with conditional display of only the fields applicable to the selected mode and direction source, plus equivalent CLI generation.
- Kept Axis definition rows active in Restrained translation mode so `move with axis` chains or residue windows receive the same constrained translation; the fitted axis itself remains unused in this mode.
- Added CLI validation, help text, run-log reporting, and regression coverage for the new transform and all direction sources.

## V4.2 - 2026-08-15

- Bumped `re_helix` to V4.2 and the bundled reciprocal-exchange engine to V3.6.
- Made the user-visible directions of a path's terminal input fragments the primary default orientation rule. Only the consecutive same-chain provenance pair directly exposed at each end establishes that direction; singleton ends and exchange jumps are neutral. This keeps closely related reciprocal-exchange products consistently directed when changing a cut position changes the relative lengths of their internal forward/reverse fragments.
- Added the default-off `--min_link_records` / `--min-link-records` CLI option and a **Min LINK records** GUI checkbox. When enabled, predicted topology-LINK count becomes the primary orientation criterion and may reverse a whole output strand. Cycles are compared after their exact shift with the requested open or closed output mode.
- Added parse-friendly orientation-policy metadata and regression coverage for default terminal preservation, explicit LINK minimization, modulo circular shifts, and junction-adjacent zero shifts.
- Added the topology-aware **Reverse Strand Direction V1.0** helper under GUI **Other tools**. It reverses selected path serializations (or reverses closed cycles around their existing first residue), continuously renumbers complete residue blocks, regenerates supported backbone LINK topology, remaps other LINKs, supports blank chain IDs, and preserves all coordinates and atom identities. Validation is scoped to selected strands so unrelated chains and their LINK records remain untouched.

## V4.1 - 2026-08-15

- Bumped `re_helix` to V4.1 and the bundled reciprocal-exchange engine to V3.5.
- Fixed RE-only processing of prior re_helix outputs by reading existing `LINK` records and reconstructing inverted P--O3', bowtie O5'--P, and standalone phosphate-bridge topology before applying new exchanges.
- Detected phosphate-only bridges structurally, so default `HETATM X33`, custom linker residue names, and phosphate-only `ATOM DA` styles remain composable across repeated reciprocal-exchange runs.
- Regenerated inherited backbone links with their new chain/residue labels, remapped non-backbone links through residue provenance, recomputed identity-symmetry distances while preserving crystallographic-mate distances, and suppressed duplicate atom-pair/symmetry links.
- Extended `--cir_shift` to accept an optional `c`/`C` suffix: plain integers retain an open serialized break, while values such as `10c` add an explicit closing `LINK` across the `TER` boundary and keep resulting cyclic strands covalently circularized.
- Preserved directed input residue serialization when choosing path and cycle orientation, using LINK count only as a secondary tie-breaker; this prevents repeated RE from reversing long strands merely to save one LINK record.
- Used a rotation-independent secondary cycle LINK cost so changing `cir_shift` moves only the break and never flips strand direction; recorded whether each cyclic component's closure LINK was written or omitted.
- Made `cir_shift` an exact signed rotation modulo cycle length: `0` now performs no strand permutation, rather than being advanced automatically to avoid a junction-adjacent break.
- Added regression coverage for link-free compatibility, LINK-aware chemical predecessors, repeated reciprocal exchange of generated output, phosphate-bridge preservation, cycle closure, and output endpoint integrity.

## V3.25 - 2026-08-13

- Bumped `re_helix` to V3.25.
- Added a compact GUI `Axis definitions line` below the individual rows. When filled, it replaces those rows; semicolons separate rows and `|` separates an axis definition from its optional `move with axis` value.
- Expanded the Axis definitions and new single-line field `?` help with the compact syntax, override behavior, and pasteable examples.

## V3.24 - 2026-07-19

- Bumped `re_helix` to V3.24.
- Made explicit Helix defs directional: the first chain in each group orients the positive fitted axis from smaller to larger residue numbers.
- Made `(AB)` and `(BA)` retain identical rigid-group membership while producing opposite axis directions for antiparallel chains, giving fixed single-site beta angles a deterministic sign.
- Made Axis definition rows directional in their entered chain order, with the first chain taking precedence over a conflicting Helix defs direction.
- Preserved ordered Helix defs and their first-chain direction references through replication and axis-coupled group construction.
- Fixed replication with explicit future-copy Helix defs such as `(AB) (CD)`: only input-present groups are treated as templates, while future groups are matched to generated copies and retain their entered direction order. Coordinate copying remains alphabetical and independent of that order, so `(DC)` still copies A to C and B to D.
- Fixed directed-axis alignment so reversing a moving Helix def, such as `(CD)` to `(DC)`, physically turns the moving helix end-for-end before beta fitting instead of silently discarding the direction difference.
- Expanded the Helix defs and Axis definitions `?` help plus CLI/README documentation with the direction convention, precedence rule, single-site beta scope, and right-/left-handed beta sign.

## V3.23 - 2026-07-16

- Bumped `re_helix` to V3.23 and Add PDB LINK Record to V1.2.
- Added automatic peptide cyclization by linking the N-terminal `N` atom to the C-terminal carbonyl `C` atom of each selected chain.
- Added nucleic-acid/peptide molecule-type controls to the Add PDB LINK Record GUI with dynamic terminal labels, availability checks, and N-C or P-O3' distance reporting.
- Extended Manual LINK creation with an independent nucleic-acid/peptide selector, dynamic P/O3' or N/C endpoint rows, and atom-specific validation for peptide links.
- Added `--peptide-chains` for peptide terminal LINK insertion from the command line and a distinct `_peptide_circ.pdb` default output name.
- Generalized LINK-driven topology rebuilding so peptide C-to-N connectivity follows the same segment ordering, residue remapping, and existing-LINK preservation rules as nucleic-acid O3'-to-P connectivity.
- Added a light-blue `?` help button with expanded positive- and negative-shift explanations and examples to the Permute Chain GUI.

## V3.22 - 2026-07-15

- Bumped `re_helix` to V3.22.
- Added the bundled Permute Chain V1.0 tool as `re_helix_lib/permute_chain.py` with standalone GUI and CLI workflows.
- Added a `Permute Chain` button to the `re_helix` GUI `Other tools` area and pre-filled the selected input PDB when available.
- Added signed cyclic residue-block shifts: positive values move residues from the chain start to its end, while negative values move residues from the end to the start.
- Added a permutation-site count with dynamic GUI rows and repeatable CLI `--permute CHAIN SHIFT` arguments so multiple chains can be rearranged in one output.
- Continuously renumbered the permuted chain from its original minimum residue number and applied the same mapping to coordinate-like records, `TER`, `HET`, `LINK`, and recognizable `REMARK` residue references.
- Rebuilt re_helix `REMARK 950 RE_SCRIPT CHAIN_RANGE` and `CHAIN_RESIDUES` inventories in the new chain order.

## V3.21 - 2026-07-15

- Bumped `re_helix` to V3.21.
- Added chain-associated virtual exchange-pair endpoints using either `A(x,y,z)` or `(x,y,z)A` syntax; virtual and real endpoints may be mixed.
- Applied virtual coordinates in the same alignment objectives as real P atoms and moved virtual points rigidly with their assigned helices.
- Made every run containing a virtual endpoint alignment-only: exchange kinds are ignored and reciprocal-exchange output is skipped.

## V3.20 - 2026-07-13

- Bumped `re_helix` to V3.20.
- Corrected user-defined-axis alignment to optimize a full XYZ translation after rotation, rather than restricting translation to the axis direction.

## V3.19 - 2026-07-13

- Bumped `re_helix` to V3.19.
- Fixed user-defined-axis alignment so movable helices can translate along the supplied axis while rotating around it.

## V3.18 - 2026-07-09

- Bumped `re_helix` to V3.18.
- Added user-defined alignment axes with `--user_axis_dir X Y Z` and `--user_axis_point X Y Z`.
- Added GUI controls in Axis definitions for direction-plus-point axis mode, including dynamic greying of the residue-range axis rows when the checkbox is enabled.
- In user-defined-axis mode, skipped fixed/moving helix-axis estimation and optimized movable helices by rotation around the supplied axis.

## V3.17 - 2026-06-28

- Bumped `re_helix` to V3.17.
- Bumped `re_helix_ccg` and `re_helix_cck` to V3.1.
- Renamed the cyclic alignment scripts to `re_helix_lib/re_helix_ccgV3_1.py` and `re_helix_lib/re_helix_cckV3_1.py`.
- Renamed alignment angle terminology from `theta` to `tau` and from `rho` to `beta`; legacy rho-angle input remains accepted as a beta alias.
- Extended `--axis_range` so whole-chain letters such as `A,B` can define an axis from all P atoms on those chains, alongside residue-window terms such as `A1-A35,B60-B26`.
- Added paired `--axis_move` definitions for moving additional whole chains or residue windows with a defined axis, for example `--axis_range A,B --axis_move C,D` or `--axis_range A,B --axis_move C1-C50,D`.
- Added GUI fields for axis-coupled move definitions, allowing triplex-like groups to be defined without stdin prompts.
- Added optional `assets/icon.png` window icons to the bundled helper GUIs.
- Updated Get Phenix Restraints so X33/internal linker support files and suggested-command arguments are generated only when a standalone linker phosphate has two explicit P--O3' LINK records.

## V3.16 - 2026-06-26

- Bumped `re_helix` to V3.16.
- Bundled Get Phenix Restraints V1.0 as `re_helix_lib/get_phenix_restraints.py`, based on `link_to_geometry_restraintsV4.py`.
- Added a `Get Phenix Restraints` button to the `re_helix` GUI `Other tools` area. The button opens the bundled GUI in a separate process and pre-fills the current input PDB when available.
- Added GUI/CLI controls for `link_distance_cutoff`, standalone 3'-to-3' linker residue name, optional `*_junctions.params` generation, and linker CIF/safe-interpretation support file generation.
- Changed the Get Phenix Restraints default `link_distance_cutoff` to `6.5` and collapsed its advanced restraint controls behind a toggle button.
- Updated Get Phenix Restraints phosphate defaults to better match Phenix/CCP4 DNA/RNA link restraints: linker P-OP distance `1.495`, mixed O-P-O angle `108.0`, and bridging O-P-O angle `103.0`.
- Updated the bundled Do Symmetry GUI to dynamically grey out fields that do not apply to the selected symmetry/alignment mode.
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

- Added optional fixed beta-angle syntax for single-site inter-helix reciprocal exchanges.
- Supported `<pos1> <pos2> <beta_deg> <kind>` for a single reciprocal-exchange site when `--axis_parallel n` is selected. This field was formerly documented as legacy `rho_deg`.
- Ignored fixed beta definitions with warnings for multi-site helix pairs or `--axis_parallel y`.

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
