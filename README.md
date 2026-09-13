# re_helix

`re_helix` is AZBMOST Package Module #2: Align Helices and Performing Reciprocal Exchanges.

It aligns nucleic-acid helices from real P-atom or chain-associated virtual-atom pairs and then applies reciprocal exchanges when every endpoint is real. It can also run reciprocal exchange only, without alignment.

Current version: V4.5

## Contents

- `re_helix.py`: main CLI/GUI script.
- `re_helix_lib/`: helper modules for PDB parsing, LINK records, and reciprocal-exchange graph handling.
- `re_helix_lib/bend_helix.py`: bundled Bend Helix tool for bending a straight two-chain helix.
- `re_helix_lib/do_symmetry.py`: bundled Do Symmetry tool for averaging a pseudosymmetric assembly into an idealized symmetric PDB.
- `re_helix_lib/add_pdb_link_record.py`: bundled Add PDB LINK Record tool for nucleic-acid P/O3' and peptide N/C terminal links, manual P/O3' links, and chain-topology rebuilding.
- `re_helix_lib/insert_virtual_resi.py`: bundled Insert Virtual Resi tool for inserting residue-numbering gaps, updating LINK endpoints, and recording each gap's range as `REMARK 950 RE_SCRIPT VIRTUAL_INSERT` header lines.
- `re_helix_lib/permute_chain.py`: bundled Permute Chain tool for cyclically rearranging and continuously renumbering one or more chains.
- `re_helix_lib/reverse_strand_direction.py`: bundled topology-aware Reverse Strand Direction tool for reversing selected nucleic-acid chain serializations without moving atoms.
- `re_helix_lib/generate_lattice.py`: bundled Generate Lattice tool for writing a P1 CRYST1 lattice record from three lattice vectors.
- `re_helix_lib/get_phenix_restraints.py`: bundled Get Phenix Restraints tool for converting LINK records into Phenix geometry restraints, optional junction movement-selection params, and linker support files.
- `re_helix_lib/check_pdb_clashes.py`: bundled Check Clashes tool for counting close heavy-atom contacts in any PDB structure, with LINK-aware, alternate-conformation, and sequence-adjacency exclusions.
- `assets/icon.png`: optional GUI/task-menu icon. The main GUI and bundled helper GUIs use it when present and fall back to the default Tk icon when it is missing.

## Requirements

- Python 3.9 or newer is recommended.
- No third-party Python packages are required for the command-line workflow.
- Tkinter is required only for GUI mode. Most python.org and system Python installs include it.

## Quick Start

Launch the GUI:

```bash
python3 re_helix.py
```

In the GUI, use the `Other tools` area to open bundled helper tools. `Bend Helix` opens the helix-bending GUI, `Do Symmetry` opens the symmetry-averaging GUI, `Add PDB LINK Record` opens the LINK-record/topology helper, `Insert Virtual Resi` opens the residue-numbering-gap helper, `Permute Chain` opens the cyclic chain-rearrangement helper, `Reverse Strand Direction` opens the topology-preserving strand-serialization helper, `Generate Lattice` opens the P1 lattice/CRYST1 helper, `Get Phenix Restraints` opens the Phenix restraint-generation helper, and `Check Clashes` opens the heavy-atom clash checker. Every tool button has a light-blue `?` button beside it that explains what that tool does, the inputs it expects, and the files it writes. If an input PDB is already selected in `re_helix`, the helper window is opened with that input pre-filled.

When the input PDB changes, the GUI updates the default `Output base` automatically unless that field has been changed to a custom value. For large exchange specifications, the `CLI pair args` field below the pair rows can be filled with the same concatenated pair tokens used on the command line; when it is filled, the individual pair rows are ignored. Likewise, `Axis definitions line` accepts a compact value such as `A,B | C,D; E,F`: semicolons separate rows, and `|` separates an axis definition from its optional `move with axis` value. When filled, it replaces the individual Axis definition rows. The **Alignment mode** area between **Axis definitions** and **Other tools** selects **Standard**, **Restrained translation**, or **Restrained rotation**. Standard mode hides the restrained controls. Each restrained mode displays only its applicable point, direction, and vector-source fields.

The main `re_helix` run log also mirrors stdout/stderr from bundled tools launched through `Other tools`, so equivalent CLI commands, selected LINK summaries, completion messages, and errors remain visible in the main window even when the helper window has no log box.

Run alignment plus reciprocal exchange from the command line:

```bash
python3 re_helix.py input.pdb '(AB)' '(CD)' 30A 8D d 13B 24C s -o model
```

This writes:

- `model_aligned.pdb`: aligned structure before reciprocal exchange.
- `model_aligned_rex.pdb`: aligned structure after reciprocal exchange.

Run restrained-translation alignment along a vector:

```bash
python3 re_helix.py input.pdb '(AB)' '(CD)' 30A 8D d 13B 24C s \
  --fix A --alignment_mode restrained_translation \
  --translation_vector 1 0 0 -o translated_model
```

This leaves the fixed helix unchanged and translates each unfixed helix without rotation. The signed distance along the direction is the least-squares value for all usable alignment pairs. The direction can instead be defined by two XYZ points or two input-PDB atoms:

```bash
--translation_points 0 0 0 1 1 0
--translation_atoms A:10:P "C:20:C4'"
```

PDB atom selectors also accept compact forms such as `A10:P` and `10A:P`, or a quoted serial such as `'#123'`. Use `_` as the chain ID for a blank-chain atom.

The direction can also be the normal to two XYZ vectors:

```bash
--translation_normal_vectors 1 0 0 0 1 0
```

This uses the normalized right-hand cross product `vector 1 × vector 2`. The vectors cannot be zero or parallel.

Run restrained-rotation alignment about the Z-axis through an XYZ point:

```bash
python3 re_helix.py input.pdb '(AB)' '(CD)' 30A 8D d 13B 24C s \
  --fix A --alignment_mode restrained_rotation \
  --rotation_axis_point 0 0 0 \
  --rotation_axis_vector 0 0 1 -o rotated_model
```

This leaves the fixed helix unchanged and optimizes only the rotation angle of each unfixed helix around the supplied line. It never translates the helix. The point may instead be selected from the input PDB, and the vector may use two XYZ points or two PDB atoms:

```bash
--rotation_axis_point_atom A:10:P
--rotation_axis_vector_points 0 0 0 0 0 1
--rotation_axis_vector_atoms A:10:P "A:20:C4'"
--rotation_axis_normal_vectors 1 0 0 0 1 0
```

Choose exactly one point option and exactly one vector option.

Run reciprocal exchange only, with no alignment:

```bash
python3 re_helix.py input.pdb 9C 23A d 23C 23F b --re_only -o model
```

This writes:

- `model_rex.pdb`: reciprocal-exchanged structure generated directly from `input.pdb`.

## Bend Helix Tool

The bundled Bend Helix V2.6 tool bends a straight two-chain nucleic-acid helix at a selected phosphorus residue. It treats the helix as two rigid pieces: piece #1 stays fixed, while piece #2 is moved by a beta bend and optional tau twist.

Open its GUI directly:

```bash
python3 re_helix_lib/bend_helix.py --gui
```

Run it from the command line:

```bash
python3 re_helix_lib/bend_helix.py --input straight_helix.pdb --pivot A36 --phi 0 --beta 30 --tau 0
```

Useful Bend Helix options:

- `--pivot A36`: P-bearing residue that marks the border between fixed piece #1 and movable piece #2.
- `--phi 0`: hinge direction around the helix axis, in degrees.
- `--beta 30`: bend angle for movable piece #2, in degrees.
- `--tau 10`: optional twist of movable piece #2 around its bent axis, in degrees.
- `--axis_range A1-A35,B60-B26`: optional local-axis range for already-bent inputs.
- `--sep y`: give movable piece #2 new chain IDs in the output.
- `--origin y`: also write an origin-overlay PDB for comparing the original and transformed helix.

### Angle screening in the GUI

Each phi, beta, and tau field has a **Screen** checkbox. Check exactly one or two angles, then click **Screening to achieve...** to define the target and the candidate grid. An unchecked angle is held at the value in its main-window field. For each checked angle, provide **From**, **To**, and **Step** values, all in degrees: the From and To values are included in the coarse search, and Step must be positive. The default coarse ranges are phi −90° to 90°, beta −180° to 180°, and tau −180° to 180°; Step defaults to 6°. With two checked angles, Bend Helix evaluates their Cartesian-product coarse grid. A live preview beneath the fields shows each calculated coarse-grid value (compacted for long grids), the value count for each angle, and the total number of coarse candidates.

After the coarse pass, Bend Helix identifies every local coarse minimum as a promising region and refines each region independently. It adjusts one or both screened angles with an adaptive pattern search, halves the refinement spacing when no improvement is found, and stops at 0.001-degree precision. Thus, From = 0, To = 10, Step = 4 starts with 0, 4, 8, and 10 degrees but can select an in-between result such as 6.35 degrees. Refining the local regions independently allows separate solution branches to be found without constructing an extremely dense full grid.

**Solution tolerance** determines which distinct results are reported. Its unit is angstroms in distance mode and degrees in rotation mode, and its default is 0.001. Qualifying solutions are de-duplicated at the 0.001-degree angle refinement precision and sorted by residual, then phi, beta, and tau. If none reaches the requested tolerance, the closest candidate is reported explicitly as a fallback. The result pane and run log show a table containing every reported solution.

By default, only the best solution and its origin-overlay PDB are written. Enable **Write all reported solutions** to write every additional solution under numbered names such as `model_P10Bm20T30_scr_sol002.pdb`; each additional model also receives a corresponding `-ori.pdb` overlay. All reported solutions use the selected `--align y` or `--align n` behavior during screening and output generation.

Every argument or source choice in the screening popup has a light-blue **?** button that opens a focused explanation and example.

The screening window provides two target modes:

- **Screening for distance** minimizes the absolute difference from a requested distance in angstroms. Define the two endpoints as either two atoms or one atom plus one XYZ point.
- **Screening for rotation** minimizes the circular difference from a requested signed angle in degrees around a defined axis. Define the endpoints as two atoms, one atom plus one XYZ point, or one atom plus the phi-corrected pivot P position. The latter is the pivot P position after applying the candidate phi correction, so it is recalculated for every candidate.

Atom selectors refer to the origin-overlay PDB, not to the chain IDs in the input PDB. Use `CHAIN:RESIDUE:ATOM` syntax, for example `A:36:P` or `C:36:O5'`. For a standard two-chain origin overlay, A and B identify the original model, while C and D identify the fully transformed model. Residue numbers and atom names remain those of the corresponding overlay atoms.

For rotation screening, the axis can be supplied geometrically with the same source choices used by the main tool's restrained-rotation mode: an XYZ point or overlay atom for the axis point, together with a direct vector, two XYZ points, two overlay atoms, or the right-hand normal to two vectors for its direction. Alternatively, choose **Local axis range(s)** and enter the range definitions directly in the screening popup; multiple definitions can be separated by semicolons. These popup-owned ranges override the main-window Local axis range(s) for candidate generation and for the best/additional screened outputs, so the evaluated and written geometries remain identical without copying values from the main window.

The measured rotation is signed from endpoint 1 toward endpoint 2 by the right-hand rule about the positive axis direction. Angular differences wrap across -180/180 degrees, so equivalent directions near the wrap boundary compare correctly. Bend Helix reports every distinct coarse or refined solution within tolerance; an exact target is not required because the closest fallback is retained when necessary. Automatic screening output names retain the selected P/B/T angle values and add `_scr`, for example `model_P0B30T0_scr.pdb` or `model_P0B30T0_scr_sep.pdb`; an explicit **Save as** path is honored without automatically adding `_scr`. A screening run automatically writes the best model's origin-overlay PDB even when origin-overlay output was not otherwise selected, so its name inherits `_scr` as well (for example, `model_P0B30T0_scr-ori.pdb`).

## Do Symmetry Tool

The bundled Do Symmetry tool generates an averaged, idealized symmetric PDB model from a pseudosymmetric homomeric assembly. It is useful when a model is expected to have rotational symmetry, such as C3 symmetry, but the coordinates are only approximately symmetric after building, editing, minimization, or conversion.

Open its GUI directly:

```bash
python3 re_helix_lib/do_symmetry.py --gui
```

Run it with explicit symmetry-related chain groups:

```bash
python3 re_helix_lib/do_symmetry.py model.pdb --groups ABCDMNOP EFGHQRST IJKLUVWX -o model_C3
```

Run it with a symmetry fold and continuous chain range:

```bash
python3 re_helix_lib/do_symmetry.py model.pdb --fold 3 --chains A-X -o model_C3
```

This writes `model_C3_symmetric.pdb`. Without `-o`, it writes `<input>_symmetric.pdb`.

Useful Do Symmetry options:

- `--groups ABCDMNOP EFGHQRST IJKLUVWX`: safest mode for noncontinuous or custom chain organization.
- `--fold 3 --chains A-X`: convenience mode when chains are continuous and evenly divisible by the fold number.
- `--fit-atoms all|p|phosphorus|backbone|ca|calpha`: choose atoms used for rigid-body alignment before averaging.
- `--keep-intermediate`: write reordered and aligned intermediate PDB files for visual checking.
- `--no-align`: average symmetry-permuted structures without fitting; use only when copies are already in the same coordinate frame.
- `--ignore-resname`: match atoms without requiring residue names to be identical.
- `--allow-missing`: average available matching atoms instead of stopping on missing atoms.

Working principle: the script builds cyclic chain permutations from the symmetry definition, reorders each symmetry-equivalent copy into the same chain organization, rigidly aligns each copy to a reference using a pure-Python quaternion/Kabsch-style least-squares fit, and averages matching atom coordinates. The result is a consensus structure that is closer to the intended symmetry than the original pseudosymmetric input.

## Add PDB LINK Record Tool

The bundled Add PDB LINK Record tool creates terminal cyclization links for nucleic acids and peptides. Its automatic GUI mode can link a nucleic-acid 5'-terminal `P` to its 3'-terminal `O3'`/`O3*`/`O3`, or a peptide N-terminal `N` to its C-terminal carbonyl `C`. Select the molecule type, then check the chains to process; the corresponding terminal-atom distance is displayed for each chain.

The GUI can also manually stage internal or inter-chain links. In the Manual LINK creation section, choose `Nucleic acid (P to O3')` or `Peptide (N to C)`, then provide the residue number and chain for both endpoints. Manual peptide mode validates the selected `N` and carbonyl `C` atoms and stages them in C-to-N topology order. GUI runs rebuild chain topology, chain IDs, TER records, residue numbering, and LINK records in one pass. Existing input `LINK` records are preserved and remapped to the rebuilt chain/residue labels when possible.

Open its GUI directly:

```bash
python3 re_helix_lib/add_pdb_link_record.py --gui
```

Run automatic chain circularization from the command line:

```bash
python3 re_helix_lib/add_pdb_link_record.py input.pdb --chains A B -o input_linked.pdb
```

Run automatic peptide cyclization by linking N-terminal `N` to C-terminal `C`:

```bash
python3 re_helix_lib/add_pdb_link_record.py peptide.pdb --peptide-chains A -o peptide_circ.pdb
```

If `--chains` is omitted, command-line mode attempts nucleic-acid circularization for every chain with usable terminal `P` and `O3'` atoms. Supplying `--peptide-chains` without chain IDs processes every chain with usable terminal `N` and `C` atoms. The default peptide output inserts `_peptide_circ` before the extension; the nucleic-acid default remains `_circ`.

Useful Add PDB LINK Record options:

- `--gui`: open the GUI for automatic and manual LINK staging.
- `--chains A B` or `--chains A,B`: choose chains for automatic terminal circularization.
- `--peptide-chains A B` or `--peptide-chains A,B`: choose peptide chains for automatic N-to-C cyclization.
- `-o output.pdb`: choose the output file.
- `-q`: suppress console messages.
- `-v` or `--version`: show the bundled tool version.

## Insert Virtual Resi Tool

The bundled Insert Virtual Resi tool inserts virtual residue-numbering gaps after selected residues. It does not create new atom records; it shifts residue numbers after each specified point. For example, inserting `3` virtual residues after `A55` changes original `A56` to `A59`, then `A57` to `A60`, and so on.

Open its GUI directly:

```bash
python3 re_helix_lib/insert_virtual_resi.py --gui
```

Run it from the command line:

```bash
python3 re_helix_lib/insert_virtual_resi.py input.pdb --insert A55 3 -o input_vresi.pdb
```

Multiple insertion points can be repeated:

```bash
python3 re_helix_lib/insert_virtual_resi.py input.pdb --insert A55 3 --insert B.20 2 -o input_vresi.pdb
```

Accepted residue token formats are `A55`, `A.55`, `55A`, and `55.A`. Multiple insertions are interpreted against the original input residue numbering. The tool updates coordinate-like records, `TER` records, and both residue endpoints of fixed-column `LINK` records using the same renumbering map.

Because no atoms are written, a virtual gap would otherwise be indistinguishable from a gap already present in the input. Each run therefore records what it did in the output header as parse-friendly `REMARK 950 RE_SCRIPT` lines, the same convention `reciprocal_exchange_pdb` uses:

```
REMARK 950 RE_SCRIPT SOFTWARE name=insert_virtual_resi version=1.1 developer=DiLiuLab
REMARK 950 RE_SCRIPT COMMAND text=...
REMARK 950 RE_SCRIPT OUTPUT_STAGE name=insert_virtual_resi
REMARK 950 RE_SCRIPT VIRTUAL_INSERT op=1 chain=A after_orig=A:55 after_new=A:55 count=3 start=A:56 end=A:58
REMARK 950 RE_SCRIPT VIRTUAL_INSERT op=2 chain=A after_orig=A:70 after_new=A:73 count=2 start=A:74 end=A:75
```

One `VIRTUAL_INSERT` record is written per insertion spec. `start` and `end` give the output-numbering range the virtual residues occupy, so the second record above reserves `A74` and `A75`. `after_orig` names the anchor residue in the input numbering and `after_new` in the output numbering; the two differ whenever an earlier gap in the same chain has already shifted that anchor. Two specs anchored at the same residue receive adjacent, non-overlapping ranges.

New records are placed after any existing `REMARK` lines and before the first structural record. A re-run appends a fresh block, so a file that passed through the tool twice keeps a readable history of both insertions.

Useful Insert Virtual Resi options:

- `--insert A55 3`: add a gap of 3 residue numbers after residue 55 in chain A.
- Repeat `--insert` for more chains or residue positions.
- `-o output.pdb`: choose the output file. Without `-o`, the default output inserts `_vresi` before the input extension.
- `--no-remark`: write the renumbered file without the `REMARK 950 RE_SCRIPT` records. The GUI has a matching checkbox, enabled by default.
- `-v` or `--version`: show the bundled tool version.

## Permute Chain Tool

The bundled Permute Chain tool cyclically rearranges complete residue blocks in one or more specified chains. A positive shift moves that many residues from the start to the end, while a negative shift moves that many residues from the end to the start. For example, shift `5` makes the original sixth residue the new first residue; shift `-5` makes the original last five residues the new first five.

The output residue numbers remain continuous and increase from top to bottom. Numbering begins at the smallest residue number used by the input chain, so a chain originally numbered `10-50` remains numbered `10-50` after permutation. The same old-to-new mapping is applied to coordinate-like records, `TER`, `HET`, both endpoints of `LINK` records, and recognizable residue references in `REMARK` lines. re_helix `REMARK 950 RE_SCRIPT CHAIN_RANGE` and `CHAIN_RESIDUES` inventories are rebuilt in the new order.

In the GUI, set `Number of permutation sites` to create that many dynamic rows, then enter one chain ID and signed shift per row. Each chain may occur once. All listed sites are applied to the same output PDB. The light-blue `?` button beside `Signed shift (+n / -n)` opens detailed movement and renumbering examples for both signs.

Open its GUI directly:

```bash
python3 re_helix_lib/permute_chain.py --gui
```

Run it from the command line:

```bash
python3 re_helix_lib/permute_chain.py input.pdb \
  --permute A 5 --permute B -5 \
  -o input_permuted.pdb
```

For compatibility, one site can also be written as positional arguments: `python3 re_helix_lib/permute_chain.py input.pdb B 5`.

Useful Permute Chain arguments:

- `--permute B 5`: one chain ID and signed shift; repeat the option for additional chains.
- Positive shifts move residues from start to end; negative shifts move residues from end to start.
- `-o output.pdb`: choose the output file. Without `-o`, the default output inserts `_permuted` before the input extension.
- `-v` or `--version`: show the bundled tool version.

## Reverse Strand Direction Tool

The bundled Reverse Strand Direction tool reverses the serialized residue order of one or more selected nucleic-acid chains while leaving every atom and coordinate unchanged. Complete residue blocks move together and are relabeled continuously in ascending order from the chain's original minimum residue number. For an open path `A10,A11,A12,A13`, the source blocks become `A13,A12,A11,A10` and are relabeled as `A10,A11,A12,A13`. For a covalently closed cycle, the current first residue remains anchored while the remaining blocks reverse, so direction changes without an unintended circular shift.

This is a PDB serialization edit, not a geometric flip, reverse complement, or chemical rebuilding of each nucleotide. Reversing residue blocks changes which backbone bonds can be represented implicitly, so the tool reconstructs the input nucleic-acid topology and regenerates required P--O3' LINK records. Existing inverted, 5'-5', 3'-3'/standalone-phosphate, cycle-closure, and other LINK topology is preserved and remapped. TER, HET, structured RE_SCRIPT residue references, and chain-inventory remarks are updated; atom serials, atom names, coordinates, ANISOU-like records, and CONECT records remain unchanged.

Launch the GUI directly:

```bash
python3 re_helix_lib/reverse_strand_direction.py --gui
```

Reverse chains A and D from the command line:

```bash
python3 re_helix_lib/reverse_strand_direction.py input.pdb \
  --strand A --strand D \
  -o input_strand_reversed.pdb
```

Useful options:

- `--strand A` or `--reverse A`: select one case-sensitive chain ID; repeat for multiple strands. Use `blank`, `<blank>`, or `_` for a blank PDB chain ID.
- `-o output.pdb`: choose the output path. Without `-o`, `_strand_reversed` is inserted before the input extension.
- `--gui`: open the graphical chain selector, which lists detected chains and residue counts.
- `-v` or `--version`: show Reverse Strand Direction V1.0.

For safety, V1.0 rejects multi-model inputs and, on selected strands, insertion codes, interleaved or TER-split serialization, cross-chain backbone components, and non-nucleic-acid chains whose required P/O3' endpoints cannot be resolved. Unselected chains are not reordered or subjected to strand-topology validation, and their LINK lines remain unchanged unless a selected endpoint must be renumbered. The tool also rejects an opposite P--O3' LINK on a selected two-residue chain, because that PDB representation is ambiguous between an inverted open path and a covalently closed two-residue cycle. The input file is never overwritten in place.

## Generate Lattice Tool

The bundled Generate Lattice tool writes or replaces the PDB `CRYST1` record for a P1 lattice from three user-provided lattice directions and distances. By default, it also rotates `ATOM`/`HETATM` coordinates and `ANISOU` tensors into the standard PDB crystallographic Cartesian frame, where `a` is along +X, `b` is in the XY plane, and `c` has positive Z. It preserves non-coordinate records, including `REMARK`, `LINK`, `TITLE`, and `SEQRES`, unless an option explicitly changes that behavior.

Open its GUI directly:

```bash
python3 re_helix_lib/generate_lattice.py --gui
```

Run it from the command line:

```bash
python3 re_helix_lib/generate_lattice.py input.pdb \
  --u1 1 0 0 --d1 80 \
  --u2 0 1 0 --d2 80 \
  --u3 0 0 1 --d3 80 \
  -o input_cryst.pdb
```

Useful Generate Lattice options:

- `--u1 X Y Z`, `--u2 X Y Z`, `--u3 X Y Z`: lattice direction vectors; each is normalized before use.
- `--d1`, `--d2`, `--d3`: distances along the three lattice directions.
- `--no-rotate-to-cryst-frame`: write/update `CRYST1` without rotating coordinates.
- `--allow-reflection`: allow an improper transform if the supplied lattice-vector order is left-handed.
- `--no-cryst1-update`: preserve existing `CRYST1` records.
- `--drop-conect`: remove `CONECT` records from the output.
- `-v` or `--version`: show the bundled tool version.

## Get Phenix Restraints Tool

The bundled Get Phenix Restraints tool is based on `link_to_geometry_restraintsV4.py`. It converts PDB `LINK` records into a `*_links.params` file containing Phenix `geometry_restraints.edits` bond restraints and phosphate-centered angle restraints. For standalone 3'-to-3' linker phosphate residues detected from two explicit P--O3' LINK records, it also writes internal P-OP1/P-OP2 bond and OP1-P-OP2 angle restraints.

Open its GUI directly:

```bash
python3 re_helix_lib/get_phenix_restraints.py --gui
```

Run it from the command line:

```bash
python3 re_helix_lib/get_phenix_restraints.py model_rex.pdb --output-base model_rex
```

By default, this writes `model_rex_links.params` when LINK records are present, `model_rex_junctions.params` when `REMARK 950 RE_SCRIPT JUNCTION` lines are present, and nonstandard-linker support files such as `X33_phenix_atomtypes.cif` and `X33_safe_interpretation.params` only when a true standalone 3'-to-3' linker phosphate is detected.

Recommended Phenix command:

```bash
phenix.geometry_minimization \
  model_rex.pdb \
  model_rex_links.params \
  model_rex_junctions.params
```

If X33 support files are generated, the suggested command printed by the tool appends `X33_phenix_atomtypes.cif` and then `X33_safe_interpretation.params`.

Useful Get Phenix Restraints options:

- `--linker-resname X33`: choose the standalone 3'-to-3' linker phosphate residue name. The default is `X33`; custom nonstandard names get matching CIF/safe params files.
- `--link-distance-cutoff 6.5`: choose the generated `pdb_interpretation.link_distance_cutoff` value. The default is `6.5`.
- `--no-junctions-params`: skip `*_junctions.params`. If you do this, use exactly one movement-selection file such as `min_P_C5.params` or a carefully prepared `min.params`.
- `--no-linker-support-files`: skip writing `<resname>_phenix_atomtypes.cif` and `<resname>_safe_interpretation.params`.
- `--include-phenix-builtin-angles`: diagnostic mode that can reproduce older duplicate-prone angle output.

Use exactly one movement-selection file for `phenix.geometry_minimization`. Usually this should be `*_junctions.params`; do not combine it with `min_P_C5.params` or `min.params` unless you deliberately want to test which top-level `selection = ...` Phenix uses.

## Check Clashes Tool

The bundled Check Clashes tool counts close heavy-atom contacts in a PDB file. The check is purely geometric and structure agnostic, so it works on nucleic acids, proteins, ligands, and mixed assemblies, and it assumes nothing about chain length, residue numbering, or composition. Hydrogens and deuteriums are ignored.

A pair is excluded from the clash count when it is not an independent contact:

- both atoms belong to the same residue;
- the atoms are alternate conformations of one site, meaning they carry different non-blank `altLoc` labels;
- the atoms lie in the same chain within `--adjacent-window` residues of each other, so ordinary covalent neighbors are not counted. Insertion-code siblings such as `5A` and `5B` count as one residue apart;
- the two residues are joined by a `LINK` record;
- the two residues close a chain that `--circular` marks as cyclic.

The LINK exclusion matters for `re_helix` output. Cyclization, reciprocal exchange, and permuted chains all create junctions where residue numbering no longer tracks real connectivity, so a bonded pair can look far apart in sequence. Without the exclusion every such junction is reported as a false clash. Use `--no-link-exclusion` to see the raw geometric contacts instead.

Open its GUI directly:

```bash
python3 re_helix_lib/check_pdb_clashes.py --gui
```

Run it from the command line:

```bash
python3 re_helix_lib/check_pdb_clashes.py model_aligned_rex.pdb
```

Check a cyclic model with a looser cutoff and save the report:

```bash
python3 re_helix_lib/check_pdb_clashes.py model.pdb \
  --cutoff 2.2 --circular auto -o model_clashes.txt
```

Useful Check Clashes options:

- `--cutoff A`: heavy-atom distance cutoff in Angstrom. Default: `1.60`, a tight geometric check for atoms driven essentially on top of each other. Raise it to flag softer nonbonded contacts.
- `--adjacent-window N`: exclude same-chain residue pairs separated by at most `N` residues. Default: `1`. Use `0` to compare every non-identical residue pair.
- `--circular off|auto|N`: close cyclic chains so the first-to-last junction is not reported. `auto` uses each chain's own lowest and highest residue number, so chains of different lengths are handled correctly. Default: `off`, because `LINK` records already cover most cyclizations.
- `--no-link-exclusion`: do not exclude residue pairs joined by a `LINK` record.
- `--label-residues LIST`: label residues for the composition breakdown, as `11,32,53` or chain-qualified `A:11,B:32`. Labeled residues report as `labeled` and the rest as `other`. Default: `none`, giving a single `heavy` class.
- `--model N`: MODEL to check in a multi-model file such as an NMR ensemble. Default: the first MODEL. Only one model is checked at a time, so ensemble members are never compared against each other.
- `--top N`: how many of the closest clash pairs to list. Default: `10`. Use `0` for counts only.
- `--no-select-commands`: omit the UCSF Chimera and UCSF ChimeraX select commands.
- `--chimera-model N`: Chimera model number used in the select commands. Default: `0`.
- `--chimerax-model N`: ChimeraX model number used in the select commands. Default: `1`.
- `-o` or `--output`: also write the report to a text file.
- `-v` or `--version`: show the bundled tool version.

The report is a tab-separated key/value block covering the atom counts, the active criteria, how many pairs each exclusion removed, the clash count, the minimum clash distance, the composition breakdown, and the closest clash pairs. `scipy.spatial.cKDTree` is used for neighbor search when SciPy is installed, and a pure-Python pair scan produces identical results when it is not.

### Selecting clashing atoms in UCSF Chimera and ChimeraX

By default the report also carries ready-to-paste selection commands. Each select-everything command sits alone on its own line, so a triple-click or any line selection picks up the whole command and nothing else:

```
chimera_select_all_clashes
select #0:5.A@C1'|#0:20.A@C1'
chimerax_select_all_clashes
select #1/A:5@C1'|#1/A:20@C1'
```

Each listed clash pair additionally gets its own pair of commands, as the last two columns of the `closest_clashes` table, so a single contact can be isolated:

```
closest_clashes	distance_A	atom_1	atom_2	chimera	chimerax
1.200000	A:DA5:C1'#3	A:DA20:C1'#4	select #0:5.A@C1'|#0:20.A@C1'	select #1/A:5@C1'|#1/A:20@C1'
```

The two programs use different atom-specifier grammars, so both forms are written out:

- Chimera uses `#model:residue.chain@atom`. A blank chain ID is a bare period, as in `#0:12.@P`.
- ChimeraX uses `#model/chain:residue@atom`. A blank chain ID has no specifier, so the chain part is dropped, as in `#1:12@P`.
- Residue insertion codes are appended to the residue number in both, as in `#0:52B.A@P` and `#1/A:52B@P`.

Chimera numbers the first opened structure `#0` and ChimeraX numbers it `#1`, which is why the defaults differ. Use `--chimera-model` or `--chimerax-model` when the structure is not the first one open. When `--model` selects a MODEL from a multi-model file, the submodel is included automatically, giving `#0.2` in Chimera and `#1.2` in ChimeraX.

#### Copying a command

These commands are full of punctuation that a text widget treats as word boundaries, so an ordinary double-click would select only a fragment such as `C1`. The GUI therefore offers three ways to copy a whole command:

- **Double-click any command in the report.** The tool highlights the entire command, from `select` to the end of that command, and copies it to the clipboard. Commands are shown on a pale blue background so they are easy to spot.
- **Use the `Copy Chimera select` and `Copy ChimeraX select` buttons** next to `Run`. These copy the select-everything command directly and stay disabled until a run produces clashes.
- **Triple-click a select-everything line**, in the GUI or in the saved report file. Those commands sit alone on their own line for exactly this reason.

The per-pair commands in the `closest_clashes` table are separated by tabs, so a line selection there would also pick up the distance and atom labels. Double-click those instead.

## What Reciprocal Exchange Means

In this package, a reciprocal exchange is a virtual topology edit on a PDB model. The operation identifies residues on two DNA strands or helices, cuts the original backbone graph at the specified sites, and reconnects the graph so the strand continuities are exchanged. This is a design operation, not an enzymatic simulation: it is meant to help build the intended crossover, junction, or bowtie connectivity before later structural refinement, sequence design, synthesis, or visualization.

The idea sits in the design tradition introduced by Nadrian C. Seeman, who founded structural DNA nanotechnology by treating DNA as a programmable construction material rather than only as genetic information. Seeman's key move was to use designed sequence asymmetry to make immobile branched junctions, avoiding the branch migration of natural Holliday junctions, so junctions could serve as predictable vertices for DNA objects, arrays, and lattices. A reciprocal-exchange operation is useful for that style of design because it gives a compact way to say: "these helices meet here, and their backbone routes trade partners here." In practice, that lets a designer specify the intended strand routing and junction topology without manually rebuilding every atom record, residue number, chain break, TER record, and LINK record.

`re_helix` supports three exchange kinds:

- `double`: exchange both local backbone continuities between two specified residues.
- `single`: exchange one strand-continuity relationship while leaving the complementary local relationship unchanged.
- `bowtie`: create paired 3'-3' and 5'-5' junction behavior, including LINK records and phosphate-only linker residues where needed. By default these linker phosphates are written as `HETATM X33`; they can instead use a custom residue name or regular `ATOM DA` records.

Background reading:

- Nadrian C. Seeman, "Nucleic acid junctions and lattices," Journal of Theoretical Biology, 1982.
- [DNA nanotechnology history overview](https://en.wikipedia.org/wiki/DNA_nanotechnology#History), including Seeman's motivation for designed immobile junctions and lattices.

## Exchange Syntax

Residue tokens can be written as `30A`, `A30`, `A.30`, or `30.A`.

Explicit Helix defs are ordered. The first chain in each parenthesized group defines the positive helical-axis direction by progressing from its smaller to larger residue numbers. Thus `(AB)` follows chain A, whereas `(BA)` follows chain B; for an antiparallel duplex these definitions point in opposite directions. Before beta fitting, the moving helix is physically oriented so its directed axis aligns with the fixed directed axis. Consequently, `(AB) (CD)` and `(AB) (DC)` can differ by an end-for-end turn and produce different structures for the same signed beta angle. The convention applies to an optional beta angle for exactly one exchange site under `--axis_parallel n`; fixed beta definitions are ignored for helix pairs with multiple exchange sites. A matching Axis definition row takes precedence over Helix defs.

Virtual atoms can be written with the chain before or after their coordinates: `A(x,y,z)` or `(x,y,z)A`. The chain ID assigns the virtual point to that helix, and coordinates are in angstroms. Virtual endpoints can be paired with real or other virtual endpoints.

```bash
python3 re_helix.py input.pdb '(AB)' '(CD)' 'A(1,2,3)' 8D d '(4,5,6)B' 24C s -o virtual_model
```

If any endpoint is virtual, the run writes only `virtual_model_aligned.pdb`. Reciprocal exchange is skipped because a virtual atom has no PDB residue; consequently, `d`, `s`, and `b` are accepted for syntax compatibility but ignored for that run.

Each exchange is:

```text
<pos1> <pos2> <kind>
```

Accepted kinds:

- `d` or `double`
- `s` or `single`
- `b` or `bowtie`

For alignment mode, a single-site inter-helix pair can also include a fixed beta angle when `--axis_parallel n` is used:

```text
<pos1> <pos2> <beta_deg> <kind>
```

Legacy examples may call this optional field `rho_deg`; it is accepted as an alias for `beta_deg`.

Example:

```bash
python3 re_helix.py input.pdb '(AB)' '(CD)' 26A 9C 90 d --axis_parallel n -o angled_model
```

## Common Options

- `-o, --output`: output base path. A `.pdb` suffix is stripped before output suffixes are added.
- `--gui`: launch the Tk GUI explicitly.
- `-v, --version`: show the app version and exit.
- `--re_only` or `--re-only`: apply reciprocal exchange only and write `<base>_rex.pdb`. When the input is a prior re_helix reciprocal-exchange output, its P--O3', O5'--P, and standalone-phosphate `LINK` records are used to reconstruct the current backbone topology before new cuts are applied. Non-backbone `LINK` records are remapped through output renumbering. By default, output paths preserve the user-visible directions of the consecutive same-chain fragments directly exposed at their ends; singleton ends and exchange jumps are neutral. Cycles preserve input serialization continuity.
- `--axis_dist 22.0`: target helix-axis distance in angstroms during alignment.
- `--axis_parallel y|n`: keep axes parallel (`y`) or allow a beta interhelical tilt (`n`). Angle terminology is `tau` = axial twist/spin of the moving helix, `phi` = orbital azimuth around the fixed helix, `beta` = interhelical tilt/bend, and `d` = axial slide.
- Beta handedness: positive beta follows the right-hand rule around `L_beta`, directed from the fixed-axis anchor toward the nearest point on the moving axis. Negative beta is the corresponding left-handed rotation.
- Helix defs `(AB)` and `(BA)`: define the same rigid chain membership but opposite directional references when A and B are antiparallel; the first chain runs low-to-high along the positive axis.
- `--axis_range B26-B60,A1-A35` or `--axis_range A,B`: define residue windows or whole chains for helical-axis estimation. Bare chain letters include the whole chain's P atoms in the axis fit. The first listed chain defines the positive axis direction from low-to-high residue number and overrides a conflicting Helix defs direction. Repeat as needed. In either restrained mode, the row remains available for grouping with `--axis_move`, but its fitted direction is not used.
- `--axis_move C,D` or `--axis_move C1-C50,D`: move additional whole chains or residue windows with the corresponding `--axis_range` row. For example, `--axis_range A,B --axis_move C` fits the axis from A/B and moves C with that axis, avoiding triplex stdin prompts. Move selections are also honored in both restrained modes, where the selected payload receives the same constrained transform.
- `--user_axis_dir X Y Z --user_axis_point X Y Z`: define a single alignment axis from a direction vector and point. When used, helix-axis estimation from P atoms is skipped and each movable helix is optimized by rotation around that line plus a full XYZ translation.
- `--alignment_mode restrained_translation`: translate each unfixed helix without rotation, constrained to exactly one direction source. Axis definitions can still group a `--axis_move` payload with the helix, but they do not determine the translation direction. The direction-plus-point axis, `axis_dist`, `axis_parallel`, and beta angles are not used in this mode.
- `--translation_vector X Y Z`: provide the restrained-translation direction directly.
- `--translation_points X1 Y1 Z1 X2 Y2 Z2`: define the direction from point 1 toward point 2.
- `--translation_atoms ATOM1 ATOM2`: define the direction from two input-PDB atoms. Selectors accept `A:30:P`, `A30:P`, `30A:P`, or a quoted serial such as `'#123'`.
- `--translation_normal_vectors U1 U2 U3 V1 V2 V3`: define the translation direction as the normalized right-hand normal `U × V`.
- `--alignment_mode restrained_rotation`: rotate each unfixed helix without translation, optimizing only its angle about the supplied point-and-vector axis.
- `--rotation_axis_point X Y Z` or `--rotation_axis_point_atom ATOM`: define the point on the restrained-rotation axis using XYZ coordinates or one input-PDB atom.
- `--rotation_axis_vector X Y Z`: provide the restrained-rotation axis vector directly.
- `--rotation_axis_vector_points X1 Y1 Z1 X2 Y2 Z2`: define the rotation-axis vector from point 1 toward point 2.
- `--rotation_axis_vector_atoms ATOM1 ATOM2`: define the rotation-axis vector from atom 1 toward atom 2 in the input PDB.
- `--rotation_axis_normal_vectors U1 U2 U3 V1 V2 V3`: define the rotation-axis vector as the normalized right-hand normal `U × V`.
- `--fix A`: keep the helix containing chain `A` fixed during alignment.
- `--replicate`: replicate the full input chain set before alignment or RE-only processing. Helix defs may name future copies: for a two-chain A/B input, `(AB) (CD)` uses A/B as the base template and assigns C/D to the generated copy; `(DC)` reverses only the copied helix's axis-direction reference. Coordinate copying always follows alphabetical base-chain order, so C copies A and D copies B regardless of Helix-def order.
- `--cir_shift 8`: apply an exact signed residue rotation, modulo cycle length, when serializing cyclic reciprocal-exchange strands; the shifted break remains open. `--cir_shift 0` preserves the canonical input-provenance start with no circular permutation.
- `--cir_shift 8c` or `--cir_shift 8C`: apply the same shift and add an explicit closing `LINK` across the output `TER` boundary, keeping each resulting cyclic strand covalently circularized. The suffix also works with other signed integers, such as `0c` or `-4C`.
- `--min_link_records` or `--min-link-records`: opt in to choosing each strand direction primarily by the topology `LINK` records that output would contain. This can reverse an entire output strand. Cycles are compared after the exact `cir_shift`, counting the omitted edge of a plain open cycle or the explicit closure of `c`/`C` mode correctly; generated inverted, bowtie, and phosphate-bridge links participate, while preserved non-topology input links do not. The default is off so closely related structures retain consistent terminal strand directions even when different cut positions change their interior fragment lengths.
- `--linker_phosphate_resname X33|NAME|DA`: choose the residue name for phosphate-only 3'-3' bowtie linker residues. `X33` is the default `HETATM` custom residue; any other 1-3 character name is written as `HETATM` by default; `DA`/`dA` writes regular `ATOM DA` while keeping only `P`, `OP1`, and `OP2`.
- `--linker_phosphate_record ATOM|HETATM`: advanced override for the inserted linker phosphate record type.

## Clone And Update

Clone creates a local copy of the GitHub repository:

```bash
git clone https://github.com/azbmost/re_helix.git
cd re_helix
```

Pull updates an existing local copy with the latest commits from GitHub:

```bash
git pull origin main
```

Run `git status` before pulling if you have local edits. If you changed files locally, commit or stash them before pulling so Git can update cleanly.

## Make The Script Executable

Make the script directly executable:

```bash
chmod +x re_helix.py
./re_helix.py input.pdb 9C 23A d --re_only -o model
```

Optionally add a short command on your PATH:

```bash
mkdir -p ~/bin
ln -s "$(pwd)/re_helix.py" ~/bin/re_helix
export PATH="$HOME/bin:$PATH"
re_helix input.pdb 9C 23A d --re_only -o model
```

For a standalone executable, PyInstaller is a practical option:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install pyinstaller
pyinstaller --onefile --name re_helix --add-data "re_helix_lib:re_helix_lib" --add-data "assets:assets" re_helix.py
```

The built executable will be under `dist/`. Platform-native app icons may require converting `assets/icon.png` to `.icns` on macOS or `.ico` on Windows.

## License

This project is released under the MIT License. See `LICENSE`.
