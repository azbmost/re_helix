# Changelog

## V4.10 - 2026-11-26

- Bumped `re_helix` to V4.10 for Build Heterohelix Junction **V1.1**, which adds `--pick J=K[,J=K]`: build junction J, and the junctions fitted together with it, from candidate K of its fit (the report's candidate list, best first) instead of the best. With `--reuse` or `--reuse-from` nothing is fitted again, so a pick takes seconds.
- Why: the score cannot tell some candidates apart. A D|L junction (a D segment joined to its own L mirror image) is its own mirror image, so its best step comes as two mirror-image solutions with the same score, their helix axes set off to opposite sides. In the D/L minicircles of the AZBMOST junction_fit_azbmost example (Example 3) the tool kept, for both junction kinds, the solutions whose offsets point the same way round the ring, 2.4-2.6 A each, so the ring cannot lie flat, and Curve_It's corner mode bought flatness with up to 3.5 deg per step of extra twist on the B-DNA ring. With the runner-up of the L->D junctions 2 and 4 (`--pick 2=2`: +0.007 in total for B-DNA, +0.657 for A-RNA) the offsets alternate, as in the RNA gallery's own D/L fit, both junctions lie within 0.7-1.2 A (RMS) of that fit, and the rings lie flat (within 0.1 deg) with at most 1.4 deg (B-DNA) and 0.14 deg (A-RNA) of extra twist per step.
- The fits are made exactly as before: the reports record fitter V1.0 (`FIT_VERSION`), so every report written by V1.0 stays reusable. A report written with `--pick` records the picks (`picks`) and keeps the fits as fitted, so a later `--reuse` without `--pick` builds the best candidates again. The model's REMARK 900 junction lines and the log name the candidate used, and the recorded command includes `--pick`.
- The tool window builds with the best candidates; take a runner-up on the command line.
- Tests: 5 new (`PickTest`: parsing, one pick for every junction fitted together, other groups left alone, refusals, the recorded command); the two version checks now read the version from the tool. 32 tests pass.

## V4.9 - 2026-11-26

- Bumped `re_helix` to V4.9 for a new bundled tool, **Build Heterohelix Junction** (V1.0, `re_helix_lib/build_heterohelix_junction.py`). It is the last entry of the `Other tools` list. Every other function, option and output of `re_helix` is unchanged.
- The tool joins two or more duplex segments end to end: B-DNA, A-RNA, Z-DNA, with D or L sugars. Each junction is fitted as one base-pair step, and only its phosphates are built. The segments stay rigid copies of their files.
- The chain senses of the two ends decide each junction. Matched gives a continuous junction (two 3'-5' links). Opposite gives a bowtie (a 3'-3' link through an X33 phosphate and a 5'-5' link). A mirror image flips the chain sense, so L-DNA is -1 like Z-DNA, and L-Z-DNA is +1. A `--link` that names the other kind is refused before any fitting, since it would turn one base pair over.
- Checks before fitting. A `Chain sense check` block lists both ends of every junction (segment number and file, end, 3' and 5' residues, sugars, chain sense) with the result. Every base pair of a segment must face the same way: neighbouring pair frames may not turn by more than 90 degrees, and the message names the pairs. A segment needs at least two base pairs. The residues each link joins need the atoms it bonds. A run that fails any of these stops with a plain message before anything is fitted.
- Check after building: every output chain is walked. The turn between neighbouring base pairs is measured inside the segments and across each junction, and a turn above 90 degrees is an error (`chain sense not matched at ...`). The cos between each pair's z axis and its segment's helix axis, in the chain's file order, shows which way the chain runs through its pairs: about +1 in a B-Z bowtie model, about -1 in an L-Z continuous model. A chain with both signs is an error too. The results go into the run log and the report.
- The fitter is Junction_Fit V1.0 with the fixes made the same day. A junction is now fitted the same way whichever end is read first and wherever its files sit in space: the contact term counts every close pair once, the quick phosphates are sampled in frames attached to the atoms, and the real phosphates (X33, 5'-5' and 3'-5') are minimized in such frames from more starts. The two links of a continuous junction are built independently. Every junction is searched from both ends, and a second pass starts again from small nudges around the 12 best results; the 16 best (was 10) get real phosphates. Junctions that join the same two end pairs share one fit.
- The tool gives the same numbers as the fixed fitter it was taken from. On the B21, L21 and Z24 master segments of the AZBMOST junction_fit_azbmost example, the complete fit results are identical for a group of two L21|Z24 continuous junctions and for a B21|Z24 bowtie junction on a small grid. A new test pins the fitted steps, totals and model coordinates of two small builds to the fixed fitter's values. With the default grid, an L21|Z24 continuous junction came out the same fitted from either end and with the L21 file turned in space (to 0.004 degrees).
- A fitting process that cannot set up its junctions now ends the run with a clear error. Before, the process pool could wait forever. Every reading of every junction is set up in the main process before the processes start, and the log reports the search progress.
- `--reuse` takes the fits from the output's own report, and `--reuse-from REPORT` from an earlier one. Fits are reused only when they were made the same way: the same tool version, the same segment files (checked by md5), the same search grid, and the same junctions fitted together. Otherwise the tool fits again and the log says why. Reports of other tools, such as Junction_Fit V1.0, are never reused.
- Output: `<first segment>_heterohelix.pdb` beside the first segment by default. It has the REMARK 900 lines of the fitter, `REMARK 950 RE_SCRIPT` `SOFTWARE`, `COMMAND` and `OUTPUT_STAGE` lines, and one `REMARK 950 RE_SCRIPT JUNCTION` line per link. Continuous links list their 3' and 5' residues; a bowtie lists the two 3' ends with the X33 between them, and the two 5' ends. Get Phenix Restraints reads these residues into its `_junctions.params`.
- The report is `tmp_file/<output stem>_heterohelix_junction.json`; for the default output it is `tmp_file/<first segment>_heterohelix_junction.json`. It also records what the fits were made from, which `--reuse` compares.
- The tool window (`--gui`, or from the main window with the current input PDB as the first segment) has four cards. Segments: one row per segment with an `upside down` box and a grey line read from the file (red when the file cannot be used). Output: the output PDB, with the report path shown grey. Junctions: grey rows worked out from the segments, and which junctions are fitted together. Options: rise, processes, fitting equal junctions together, and reusing an earlier report (its path becomes an input when ticked). `Run` streams the command-line build into the run log, and `Stop` ends it. A failed run shows its error and first reason in the footer.
- Without numpy, scipy or shapely the tool stops with a one-line message, and with a message box when it was opened as a window.
- New tool icon in the shared AZBMOST look: two duplex segments on one axis, a right-handed ladder below and a left-handed zig-zag above, joined by two links marked with dots. `azbmost_icons.py` (2.1.0) and `azbmost_ui.py` (1.1.0, the tool's one-line purpose) are copied unchanged to the AZBMOST launcher, BNP_NA and Curve_It.
- Added `tests/test_build_heterohelix_junction.py` (27 tests, about a minute) with three 8-bp fixtures in `tests/data`: (CG)4 Z-DNA cut from the Z24 master, a self-complementary B-DNA (GCTTAAGC) built from the B21 master's residues along its own screw, and its x -> -x L-DNA mirror. The A-RNA duplex `R24_A-RNA_duplex.pdb` and mirror images made in the tests cover A-RNA, L-RNA and L-Z-DNA.
- The tests cover the chain senses and link kinds, and every refusal before fitting (the fit is replaced, so nothing is fitted). They check the search score read from either end and with a segment moved in space, and the forward/backward step round trip. Small command-line builds are checked for their LINK records, bonds, `RE_SCRIPT` lines and the residues each link joins, and Get Phenix Restraints reads those lines. Other tests cover the check after building (also with one pair turned over), `--reuse` and `--reuse-from`, the fitted numbers against the fixed fitter's, the missing-package message, and the window.

## V4.8 - 2026-11-23

- Bumped `re_helix` to V4.8 for a new option, `--xo_faces` (GUI: **Crossovers on faces**, with its `?` help), off by default, so every existing command and output is unchanged.
- At a parallel crossover a strand comes in and leaves on the same face of the plane through the two helix axes, so it should stay on that face at the crossover (Han et al. 2017, Fig. 1B). A plain exchange moves no atoms and leaves one atom of each crossing strand, the leaving O3' or the entering phosphate, 0.1-2 A on the other side, where the two crossing strands meet. With `--xo_faces`, every crossing strand is traced after the exchange (in = its phosphate one step before the crossing, out = one step after) and, when both lie on one face, the crossing is put on that face: the entering phosphate is turned about its gamma and beta bonds, with a small sugar turn about the glycosidic bond or a small turn of the leaving base pair about the helix axis when those suffice, the neighbouring phosphates keeping their bonds at the model's own O3'-P geometry. Antiparallel crossings (in on one face, out on the other) are left as built.
- Only rigid turns are used: bond lengths and angles inside each nucleotide are kept, bases stay paired, links are not lengthened (more than 0.3 A) or closed, and moved atoms keep 2.5 A from every atom more than three bonds away. A crossover is placed whole or left exactly as built, and each crossing strand is reported in the run log and in `REMARK 950 RE_SCRIPT XO_FACES` header lines. `--xo_min_face` (default 0.5 A) sets how far the crossing atoms must sit on their face.
- Crossovers whose leaving O3' is on the wrong face need a local twist change over several base pairs, which this option does not make; they are left as built and reported. On the lab's models: A-RNA Han-type sheet 34 of 36 crossing strands placed, DPE 4 of 4, DPOW, DPON and the 7:4 PX none; B-DNA DPON 4 of 4, DPE 2 of 4, DPOW none; DAE and DAO (antiparallel) untouched, atoms identical.
- Added `re_helix_lib/crossover_faces.py` (Crossover Faces V1.0), which the option calls and which also runs on its own on any exchanged model with its pre-exchange reference.
- Added `tests/test_crossover_faces.py` with a 24-bp A-RNA duplex fixture (`tests/data/R24_A-RNA_duplex.pdb`): the option is off by default; on a parallel double crossover every crossing ends at least 0.5 A on its face, only crossing nucleotides move, their bonds keep their lengths and moved atoms make no new contacts; an antiparallel double crossover keeps every atom; the GUI command carries the flag. Each test was mutation-checked against a broken implementation (phosphates never placed, antiparallel crossings treated as parallel, GUI flag dropped).

## V4.7 - 2026-11-19

- Bumped `re_helix` to V4.7 for a redesign of every window in the shared AZBMOST look, the same one the launcher, BNP_NA and Curve_It now use. The look lives in `re_helix_lib/azbmost_ui.py` (theme and window parts) and `re_helix_lib/azbmost_icons.py` (icons drawn as vector art), identical copies of the files shipped with the other modules.
- The main window has a header with the RE_Helix icon, version and purpose; Inputs, Options, Exchange pairs, Axis definitions and Alignment mode as cards on the left; the `Other tools` as a list with an icon, a one-line purpose and a round `?` per tool on the right, with the `Run log` beneath it, so the log is always in view instead of below the fold; and `Close` and `Run` in a footer.
- Restyled all nine bundled tool windows to match, keeping each tool's own version: Bend Helix (V3.1, including its screening window, which now scrolls instead of clipping its lowest rows on a short screen), Do Symmetry, Add PDB LINK Record, Insert Virtual Resi, Permute Chain, Generate Lattice, Get Phenix Restraints, Reverse Strand Direction and Check Clashes. Each has a tool header, its settings in cards, round `?` buttons that open the same help as before, its log in a quiet monospace box, and its buttons in a footer with the main action on the right.
- Redrew the RE_Helix icon in the new AZBMOST icon family.
- Only the look changed: every field, option, default, help text, command, output file and log line is as before. Two checks in `tests/test_bend_helix_screening.py` that pinned the old look (the light-blue `?` colour and the parent frame of the shift entries) now pin the new one.

## Bend Helix V3.1 - 2026-10-17

- Bumped the bundled Bend Helix tool to V3.1 while keeping the main `re_helix` application at V4.6.
- Named screening atoms in the final bent model by default. A new **Atom names** choice in the screening window offers **Final bent model**, where `A:150:P` is residue A150 where it ends up after every bend, by the written model's own chain IDs, and **Origin overlay helices**, the previous naming, where A/B is the unmoved input helix and each bend adds the two chains of the full helix it generates. With several pivots the overlay naming made a residue on a moved piece reachable only through the helix that supplies that piece, such as `E:150:P` for A150 behind pivots A118 and A149, while `A:150:P` silently named the unmoved input copy.
- Followed `--sep y` in final-model names: a moved piece is named by the chain IDs the written model gives it, computed by the same routine that writes them, and a name on the wrong chain is refused with the right one, such as `use E:12:P`.
- Refused a screen whose target cannot change: when both endpoints, and for rotation the axis, all stay on the unmoved input or all move rigidly with one helix, every candidate gives the same value and the screen could only report its first grid value. The error names the selector that follows the bent model, for example `use E:150:P, or name it in the final model as A:150:P`. A one-candidate grid only measures and is still allowed.
- Kept refining a promising region below the 0.001-degree or 0.001-angstrom precision, down to 1e-6, while its residual is above the solution tolerance. One value shared by two pivots moves the target about twice as fast as one bend, so a tau screen of A150 against A91 behind pivots A118 and A149 stopped 0.0009 degree from the root, with a 0.0019-degree residual, and reported no solution within the 0.001-degree tolerance. It now reports tau = -3.71973 degrees within tolerance, and the written model puts A150 in register with A91 to 0.0001 degree. Regions that already meet the tolerance stop where they did before.
- Listed, under the **Atom names** choice, which residues each piece holds and which chain names reach them, for example `#3 A149-A235, B1-B87`, for the current input, pivots, and `--sep` setting.
- Renamed **Overlay atom** to **Atom** throughout the screening window, since atoms can now be named either way, and made selector parse errors name-neutral.
- Added regression coverage for final-model names with and without `--sep y` against the written files, refused and allowed unchangeable targets, the corrected selector in the error, refinement to the tolerance for a value shared by two pivots, the piece descriptions in both namings, and a GUI run that screens with final-model names. Every new test was mutation-checked against a deliberately broken implementation. Stopped the GUI tests from creating a throwaway Tk root to probe for a display, after which a later root's `update()` can stall on macOS.

## Bend Helix V3.0 - 2026-10-16

- Bumped the bundled Bend Helix tool to V3.0 while keeping the main `re_helix` application at V4.6.
- Allowed several pivots in `--pivot` and the GUI **Pivot P residue(s)** field, separated by spaces, commas, or semicolons, for example `A20, A40`. The bends are applied in the listed order.
- Made each bend act on the full helix the previous bend generated. Bend 1 acts on the input helix as before; bend k acts on the rigidly bent copy of the whole helix that bend k - 1 produced. With n pivots the duplex splits into n + 1 pieces: piece #1 comes from the input helix and piece #k + 1 from the helix bend k generated, and the output links them and drops the unused part of every generated helix.
- Required every pivot after the first to lie in the movable piece of the bend before it, i.e. after it along the helix axis. A pivot listed out of order, or two pivots on one base-pair border (including the two residues of one pair), is rejected with an error naming both pivots rather than producing a torn model.
- Built every bend from its own pivot's frame on the helix it acts on, carried rigidly from the input. A generated helix is a rigid copy of the input, so this is the frame Bend Helix would fit to that helix itself; re-fitting a whole-duplex axis to the partly bent model instead would tilt every later bend by an amount that depends on the shape so far. Phi, beta, tau, Sa, and Sr therefore keep their single-bend meaning at every pivot, and the angle between consecutive pieces' axes is that bend's beta.
- Let `--phi`, `--beta`, `--tau`, `--shift_axial`, and `--shift_radial`, and the matching GUI fields, take one value per pivot, in pivot order, or a single value used at every pivot, with the same separators. Any other count is an error, and omitted tau, Sa, and Sr stay 0 at every pivot.
- Made the bend options repeatable, so `--pivot A20 --pivot A40` equals `--pivot A20,A40` instead of silently keeping the last value, and attached a list starting with a minus sign to its option before parsing, so `--tau -10,5` works as typed. argparse would otherwise read `-10,5` as an unknown option.
- Recorded several pivots and per-pivot values in the `COMMAND` remark as comma lists, writing a negative list as `--tau=-10,5` so the recorded command replays with any argparse. A one-pivot run writes the same files as V2.9, with the same names, coordinates, and `COMMAND` record; only the version in the `SOFTWARE` record changes.
- Named multi-pivot outputs in bending order: `model_PvA20PvA40_P0B30T0.pdb` when every value is shared, or one block per pivot such as `model_PvA20_P0B30Tm10_PvA40_P90B30T5.pdb`. A name that would pass 240 bytes is shortened to the first and last pivot, the pivot count, and a digest, such as `model_PvA10-A140x14_5e3c1d9a0b.pdb`, keeping room for the `-ori` and `_solNNN` suffixes under the common 255-byte filename limit.
- Applied `--align y` at every bend, so each pivot P returns to its shifted position on the helix that bend acts on and every junction stays connected. Each pivot uses the `--axis_range` containing it, and ranges that order the base pairs differently, such as two ranges running in opposite directions, are rejected.
- Gave every movable piece its own new chain IDs under `--sep y`, in piece order, and wrote the input helix plus every generated full helix to the `-ori` PDB under `--origin y`: A/B input, C/D from bend 1, E/F from bend 2, and so on. Wrapped overlay atom serials past 99999, which a long helix with many bends can reach, so the five-column field never spills into the next one.
- Kept screening for several pivots while phi, beta, tau, Sa, and Sr each hold a single value: every fixed value and every screened candidate applies at each pivot, the overlay namespace grows by two chains per bend, and **Phi-corrected pivot P** and the **Local axis range(s)** rotation axis refer to the first pivot. Entering more than one value in any of these fields turns screening off, clears the **Screen** checkboxes, and says why in the screening status line until every field holds one value again.
- Reported every piece with its base-pair range and source helix, and every bend with its own axis, hinge, twist axis, and alignment translation, in the run summary for a multi-pivot run. The single-pivot summary is unchanged.
- Updated the pivot, value, screening, axis-range, `--sep`, `--origin`, and output-name help in the GUI and CLI, the screening popup help for the larger overlay namespace, and the README.
- Added regression coverage for the chained geometry against a chain built step by step with the single-pivot code, the per-bend kink angles fitted from the overlay helices, junction connectivity, piece sources, shared and per-pivot values, pivot-order and axis-range validation, value-count errors, byte identity between one-pivot list and scalar inputs, the overlay and `--sep` chain IDs, naming and the long-name digest, the recorded-command round trip with negative lists, repeated and negative CLI options, multi-pivot screening, and the GUI switch that turns screening off. Every new test was mutation-checked against a deliberately broken implementation.

## Bend Helix V2.9 - 2026-10-14

- Allowed multiple atoms in both screening endpoint overlay fields, separated by commas or semicolons. Each endpoint uses the arithmetic mean of its atoms' XYZ coordinates, with equal weight per atom, for both distance and signed-rotation screening.
- Resolved each atom in every candidate's origin-overlay geometry before averaging, including groups that mix original and transformed overlay chains. Preserved single-atom entries and the XYZ/pivot endpoint choices.
- Updated the endpoint labels, help, and README with the list syntax and averaging behavior, and added regression coverage for centroids, candidate-dependent screening, invalid lists, and existing selector formats.

## Get Phenix Restraints V1.1 - 2026-10-11

- Bumped the bundled Get Phenix Restraints tool to V1.1 while keeping the main `re_helix` application at V4.6.
- Gave every standard 3'-5' P--O3' LINK that Phenix does not link by itself, such as a strand closure from `add_pdb_link_record.py` or a link between non-consecutive residues, the restraints Phenix applies to each ordinary backbone link (its `rna3p` chain link): the O3'-P bond (1.607 A, sigma 0.015), O3'-P-O5' (104.0, sigma 1.5), O3'-P-OP1/OP2 (108.0, sigma 3.0), the O3'-side angle C3'-O3'-P (120.2, sigma 1.5), and the epsilon, zeta and alpha torsions (-140/35/period 1, 172/30/period 3, -60/20/period 3).
- Fixed the closure geometry this left behind. V1.0 wrote only the bond and the P-centered angles, so the unrestrained C3'-O3'-P angle absorbed the strain of the junction: in three minimized 210-bp closed duplexes it ended at 105-106 degrees, with O3'-P-O5' at 94-95, while the other 418 links sat at 118 and 101. With V1.1 both closures of each model match the other links.
- Wrote the restraints as explicit `geometry_restraints.edits` bonds, angles, and periodic dihedrals rather than a `pdb_interpretation.apply_cif_link`. Phenix drops an `apply_cif_link` bond whenever the input O3'-P distance exceeds `link_distance_cutoff` (3.0 A unless a selection file raises it), which would silently unbond a widely gapped closure. The explicit edits minimize to the same coordinates as `apply_cif_link rna3p` (largest difference 0.001 A on the three test models).
- Applied the new set only where it cannot duplicate Phenix's own restraints. LINKs between consecutive residues of one chain, which Phenix links itself, keep the V1.0 bond restraint. 3'-to-3' linker phosphates and 5'-to-5' linkages are unchanged. LINKs involving a nonstandard residue name or a missing backbone atom fall back to V1.0 output with a NOTE in the run summary.
- Took the bond, angle, and torsion values for these LINKs from `rna3p`, so `--sigma`, `--angle-sigma`, and the `--angle-*` ideals now apply to the remaining LINK-derived restraints only.
- Added `--legacy-link-restraints` and a matching checkbox under **Advanced restraint values** to write the V1.0 restraints. Its output is byte-identical to V1.0's on the three test models.
- Added regression coverage for the restraint set, legacy output, the consecutive-residue, nonstandard-residue, missing-atom, 3'-to-3', and 5'-to-5' cases, legacy atom names, and the CLI flag. Each test was mutation-checked against a deliberately broken implementation.

## Bend Helix V2.8 - 2026-09-12

- Bumped the bundled Bend Helix tool to V2.8 while keeping the main `re_helix` application at V4.6.
- Led automatic output names with a `Pv` block naming the pivot residue: `model_PvA36_P0B30T0.pdb`, and `model_PvA36_P0B30T0Sa2Srm1p5.pdb` when a pivot shift is nonzero. Bends pivoted at different residues previously collided on one filename whenever the angles matched.
- Took the pivot label from the residue as resolved in the input rather than the typed token, so `A36`, `36A`, and `A.36` all name the same file. Blank chains and negative residue numbers stay filename-safe as `_36` and `Am12`.
- Carried the pivot block through the `-ori` overlay and every numbered `_solNNN` screening solution, so a screening run's files all agree on the pivot they came from.
- Recorded each run in every output PDB as `REMARK 950 RE_SCRIPT` SOFTWARE, COMMAND, and OUTPUT_STAGE lines, the same convention `reciprocal_exchange_pdb` and `insert_virtual_resi` already use and that `get_phenix_restraints` and `reverse_strand_direction` parse. The COMMAND record holds the full equivalent CLI command, so a bent model carries the recipe that produced it.
- Made the GUI pass the command it prints into the writer, so the command in the result pane and the command in the file are the same string rather than two independent reconstructions.
- Placed the records after any existing `REMARK` lines and before the first structural record, so they never land among `LINK` or coordinate records even when a file appends `REMARK` lines after the coordinates. Re-running the tool on its own output appends a fresh block and keeps the earlier one.
- Wrote the same block into the `-ori` origin-overlay PDB, which is built from scratch rather than from the input records and would otherwise have carried no provenance at all.
- Flattened embedded newlines in the COMMAND value so one run can never produce more than one record, and documented that `read_pdb` reads universal newlines, so a CRLF input has always been written back as LF.
- Added regression coverage for the naming block, the resolved-pivot normalization, record content and placement, overlay parity, the GUI-supplied command, re-run accumulation, and a round trip that replays the recorded command and compares the result byte for byte against the original output.

## Bend Helix V2.7 - 2026-09-12

- Bumped the bundled Bend Helix tool to V2.7 while keeping the main `re_helix` application at V4.6.
- Added two pivot-shift variables that translate movable piece #2 before any rotation is built: `--shift_axial`/`--sa` (**Sa**) along the local helix axis and `--shift_radial`/`--sr` (**Sr**) along the pivot P atom's own radial direction, both in angstroms and both defaulting to 0. Positive Sa moves piece #2 forward, away from fixed piece #1; positive Sr moves it outward, away from the axis.
- Rebuilt the hinge, the phi circle, and the tau twist axis from the shifted pivot. Sa slides the hinge circle along the axis, Sr changes its radius to `radius + Sr`, and the twist axis follows piece #2's own translated helical axis so tau still spins the piece about its own centerline.
- Retargeted `--align y` at the shifted pivot rather than the original one. Aligning back to the unshifted pivot would cancel Sa and Sr for every phi, beta, and tau, leaving the two settings with no effect at all.
- Documented the consequence that follows from `--align y`: with the default alignment, the net effect of Sa and Sr is exactly a rigid translation of the bent piece #2 by the shift vector, because aligning the pivot cancels every hinge *location* from the result and leaves only the alignment target. The rebuilt hinge position is observable only with `--align n`. This is a mathematical identity, not a policy, so a future change that retargets the alignment would silently disable the feature.
- Kept a zero shift bit-identical to V2.6. The shift offsets are exactly zero-preserving, the hinge tangent still comes from the same expression whenever the shifted radius is physically meaningful, and a regression test compares against an inlined pre-shift implementation with exact equality rather than a tolerance.
- Allowed `Sr = -radius` and beyond. The hinge then sits on the helix axis, which turns beta into a pure axial kink, or passes through to the far side. The hinge direction is derived from the unit radial in that regime, so the sense of beta stays continuous across the crossing instead of silently reversing, and a near-zero shifted radius no longer trips the degenerate-tangent guard.
- Made both shifts screenable alongside phi, beta, and tau. The screening pool is now five variables, still one or two at a time, with per-variable units: degrees for the angles and angstroms for the shifts. Coarse ranges default to -10 to 10 A with a 0.5 A step, refinement runs to 0.001 A, and solutions are de-duplicated at each screened variable's own precision.
- Extended the candidate cache key, the tie-break sort key, and the reported solution record to all five variables. A cache still keyed on the three angles would have collapsed every shift candidate in a grid onto one evaluation and reported the first grid point as the optimum, silently.
- Defaulted an unscreened shift to 0 instead of requiring it, so existing callers that supply only phi, beta, and tau keep working. `ScreeningResult.angles` and `ScreeningSolution.angles` still return the three angles; a new `.variables` property returns all five.
- Added **Sa** and **Sr** rows to the main GUI window, each with a **Screen** checkbox and a light-blue contextual **?** button, and moved the screening-target row and everything below it down two rows. The screening popup's grid headers dropped their `(deg)` suffix and each row now carries its own unit, since a mixed pair such as beta and Sr spans two units.
- Appended an `Sa...Sr...` block to automatic output names whenever either shift is nonzero, for example `model_P0B30T0Sa2Srm1p5.pdb`, and left the historical `model_P0B30T0.pdb` name untouched at zero shift.
- Reported the shift values, shift vector, shifted pivot position, and shifted radius in the run summary, and added the two shifts to the equivalent CLI command the GUI prints.
- Documented the conventions that are easy to get wrong: Sr does not rotate with phi, so at phi = 180 a positive Sr moves piece #2 toward the side opposite the hinge; Sa follows the pre-bend axis rather than piece #2's bent axis; and the radial direction comes from the P atom of the residue actually named, so naming the partner residue of the same pivot base pair gives the same two pieces but reverses Sr.
- Added regression coverage for the shift geometry, the alignment identity, the axis-crossing continuity, the origin-overlay agreement, output naming, the CLI, and shift screening. Every test was mutation-checked against a deliberately broken implementation.

## Insert Virtual Resi V1.1 - 2026-09-12

- Recorded every run in the output header as parse-friendly `REMARK 950 RE_SCRIPT` lines, matching the convention `reciprocal_exchange_pdb` already uses. Previously the output carried no trace of which gaps were opened, and a virtual gap was indistinguishable from a gap already present in the input.
- Wrote one `VIRTUAL_INSERT` record per insertion spec, giving the chain, the anchor residue in both input and output numbering, the count, and the output-numbering range the virtual residues occupy: `REMARK 950 RE_SCRIPT VIRTUAL_INSERT op=1 chain=A after_orig=A:55 after_new=A:55 count=3 start=A:56 end=A:58`.
- Reported `after_new` alongside `after_orig` because an earlier gap in the same chain shifts the anchor of every later spec, so an anchor quoted in input numbering alone is ambiguous. For `A55 3` followed by `A70 2`, the second anchor is written out as `A:73` and its virtual range as `A:74`-`A:75`.
- Gave two specs anchored at the same residue adjacent, non-overlapping ranges rather than two overlapping claims on the same span.
- Added the `SOFTWARE`, `COMMAND`, and `OUTPUT_STAGE` header records used by the other bundled tools, so a re-run appends a new block and the file keeps a readable history of successive insertions.
- Placed new records after any existing `REMARK` lines and before the first structural record, so they never land among `LINK` or coordinate records even when a file appends `REMARK` lines after the coordinates.
- Added `--no-remark` and a matching GUI checkbox to write the plain renumbered file. The equivalent CLI command shown by the GUI carries the flag.
- Reported the virtual residue ranges in the run summary, which previously listed only the specs and the per-record change counts.
- Added regression coverage for the range arithmetic, the record format, header placement, the `--no-remark` path, and the guarantee that the REMARK block leaves the renumbering itself unchanged.

## V4.6 - 2026-09-05

- Bumped `re_helix` to V4.6 and added the bundled Check Clashes V1.0 tool at `re_helix_lib/check_pdb_clashes.py`.
- Added **Check Clashes** to the GUI **Other tools** area. The button launches the bundled tool with the currently selected input PDB pre-filled, and mirrors its stdout/stderr into the main run log like the other helpers.
- Counted close heavy-atom contacts using a purely geometric, structure-agnostic check that assumes nothing about chain length, residue numbering, or composition, so nucleic acids, proteins, ligands, and mixed assemblies are all supported. Hydrogens and deuteriums are ignored.
- Excluded same-residue pairs, alternate conformations carrying different non-blank `altLoc` labels, same-chain sequence neighbors within `--adjacent-window`, and residue pairs joined by a `LINK` record. Insertion-code siblings such as `5A` and `5B` count as one residue apart, so the default window excludes them while `--adjacent-window 0` still compares every non-identical residue pair.
- Made the LINK exclusion the default because cyclization, reciprocal exchange, and permuted chains create junctions where residue numbering no longer tracks real connectivity; without it every such junction is reported as a false clash. `--no-link-exclusion` restores the raw geometric contacts.
- Added `--circular off|auto|N` for cyclic chains. `auto` uses each chain's own lowest and highest residue number so chains of different lengths are handled correctly; the default is `off` because LINK records already cover most cyclizations.
- Added `--label-residues` for the composition breakdown, accepting plain residue numbers or chain-qualified `A:11` entries, and `--model` for selecting one MODEL from a multi-model file so ensemble members are never compared against each other.
- Reported per-exclusion pair counts, the clash count, the minimum clash distance, the composition breakdown, and the closest clash pairs as a tab-separated key/value block, with an optional `-o` report file.
- Added ready-to-paste UCSF Chimera and UCSF ChimeraX selection commands to the report. `chimera_select_all_clashes` and `chimerax_select_all_clashes` select every clashing atom at once, and each listed pair carries its own commands as the last two columns of the `closest_clashes` table.
- Wrote each program's own atom-specifier grammar: Chimera `#model:residue.chain@atom` with a bare period for a blank chain, and ChimeraX `#model/chain:residue@atom` with the chain part dropped for a blank chain. Insertion codes are appended to the residue number in both, and atoms shared by several clashes are listed once.
- Defaulted the model numbers to Chimera `#0` and ChimeraX `#1`, matching how each program numbers the first opened structure, with `--chimera-model` and `--chimerax-model` to override. A `--model` selection from a multi-model file becomes the submodel automatically, giving `#0.2` and `#1.2`.
- Added a **Chimera / ChimeraX select commands** GUI section with its own **?** help and model-number fields, and `--no-select-commands` to suppress the commands entirely.
- Made the commands easy to copy. Selection commands contain punctuation that a text widget treats as word boundaries, so an ordinary double-click would select only a fragment such as `C1`. Double-clicking any command in the GUI report now highlights the entire command and copies it to the clipboard, and commands are drawn on a pale blue background so they are easy to find.
- Added **Copy Chimera select** and **Copy ChimeraX select** buttons beside **Run**. They copy the select-everything command without touching the report and stay disabled until a run produces clashes.
- Placed each select-everything command alone on its own line, under its `chimera_select_all_clashes` or `chimerax_select_all_clashes` key, so a triple-click or any line selection in the GUI or the saved report file picks up the whole command and nothing else.
- Used `scipy.spatial.cKDTree` for neighbor search when SciPy is present and fell back to a pure-Python pair scan with identical results when it is not, keeping the command-line workflow free of third-party requirements.
- Added light-blue contextual **?** buttons to every GUI field and regression coverage for the geometry, each exclusion rule, option parsing, backend parity, the report format, and the CLI.
- Added a light-blue contextual **?** button beside every button in the main **Other tools** area, covering all nine bundled tools. Each explains what the tool does, the inputs it expects, and the files it writes, so the tools are documented without leaving the main window.
- Rebalanced **Other tools** to five buttons on the first row and four on the second. With the added **?** buttons the measured rows are 987 px and 878 px against the 1320 px window, where the previous seven-button first row would have overflowed.

## Bend Helix V2.6 - 2026-08-30

- Bumped the bundled Bend Helix tool to V2.6 while keeping the main `re_helix` application at V4.5.
- Replaced the single-winner screening report with a deterministic table of every distinct coarse or refined solution within a configurable target tolerance; when no candidate qualifies, the closest result is reported as an explicit fallback.
- Refined every promising local coarse region independently so separate solution branches can be discovered, while retaining the adaptive 0.001-degree search precision and the selected `--align y` or `--align n` geometry.
- Added **Solution tolerance** and **Write all reported solutions** controls with light-blue contextual help. Additional models use numbered `_solNNN` names and each receives its own origin-overlay PDB.
- Added a dedicated **Local axis range(s)** field inside the screening popup. When local-axis rotation screening is selected, these popup ranges override the main-window ranges throughout candidate evaluation, equivalent CLI reporting, and every screened output.

## Bend Helix V2.5 - 2026-08-30

- Bumped the bundled Bend Helix tool to V2.5 while keeping the main `re_helix` application at V4.5.
- Added a **Screen** checkbox beside each phi, beta, and tau field. A screening run may vary exactly one or two of these angles over inclusive **From**/**To** coarse grids with a positive **Step**, all in degrees, while keeping every unchecked angle fixed at its entered value.
- Added the **Screening to achieve...** window with distance and signed-rotation objectives. Distance screening supports two overlay atoms or one overlay atom and one XYZ point; rotation screening also supports an overlay atom and the phi-corrected pivot P position.
- Added rotation-axis definitions using the same geometric point/vector sources as restrained rotation in `re_helix`, or the Bend Helix **Local axis range(s)** definition.
- Made screening atom selectors use the origin-overlay namespace and `CHAIN:RESIDUE:ATOM` syntax. In the standard two-chain overlay, chains A/B are the original model and C/D are the fully transformed model.
- Evaluated every one- or two-angle coarse-grid candidate, then efficiently searched between steps around up to eight of the strongest results with adaptive refinement down to 0.001-degree precision. The closest coarse or refined candidate is retained even when the exact target is unavailable, and screening automatically writes its origin-overlay PDB.
- Added light-blue contextual **?** buttons with explanations and examples for every argument and source choice in the screening popup.
- Changed the default coarse screening step from 5° to 6° and added a live calculated-grid preview with per-angle values/counts and the total Cartesian-product candidate count.
- Set the default screening range to −90° through 90° for phi; beta and tau use −180° through 180°.
- Added `_scr` to automatic screening output names after the selected P/B/T angle values, for example `model_P0B30T0_scr.pdb` or `model_P0B30T0_scr_sep.pdb`.

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
