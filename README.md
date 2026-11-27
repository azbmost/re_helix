# re_helix

`re_helix` is AZBMOST Package Module #2: Align Helices and Performing Reciprocal Exchanges.

It aligns nucleic-acid helices from real P-atom or chain-associated virtual-atom pairs and then applies reciprocal exchanges when every endpoint is real. It can also run reciprocal exchange only, without alignment.

Current version: V4.9

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
- `re_helix_lib/crossover_faces.py`: the `--xo_faces` step (also usable on its own), which keeps each strand of a parallel crossover on its own face of the plane through the two helix axes.
- `re_helix_lib/check_pdb_clashes.py`: bundled Check Clashes tool for counting close heavy-atom contacts in any PDB structure, with LINK-aware, alternate-conformation, and sequence-adjacency exclusions.
- `re_helix_lib/build_heterohelix_junction.py`: bundled Build Heterohelix Junction tool for joining duplex segments of any helix form (B-, A-, Z-forms, D or L sugars) with junctions fitted as one base-pair step, continuous or bowtie as the chain senses of the two ends decide.
- `assets/icon.png`: optional GUI/task-menu icon. The main GUI and bundled helper GUIs use it when present and fall back to the default Tk icon when it is missing.

## Requirements

- Python 3.9 or newer is recommended.
- No third-party Python packages are required for the command-line workflow.
- The bundled Build Heterohelix Junction tool needs numpy, scipy and shapely.
- Tkinter is required only for GUI mode. Most python.org and system Python installs include it.

## Quick Start

Launch the GUI:

```bash
python3 re_helix.py
```

In the GUI, use the `Other tools` area to open bundled helper tools. `Bend Helix` opens the helix-bending GUI, `Do Symmetry` opens the symmetry-averaging GUI, `Add PDB LINK Record` opens the LINK-record/topology helper, `Insert Virtual Resi` opens the residue-numbering-gap helper, `Permute Chain` opens the cyclic chain-rearrangement helper, `Reverse Strand Direction` opens the topology-preserving strand-serialization helper, `Generate Lattice` opens the P1 lattice/CRYST1 helper, `Get Phenix Restraints` opens the Phenix restraint-generation helper, `Check Clashes` opens the heavy-atom clash checker, and `Build Heterohelix Junction` opens the junction builder that joins duplex segments of any helix form. Every tool button has a light-blue `?` button beside it that explains what that tool does, the inputs it expects, and the files it writes. If an input PDB is already selected in `re_helix`, the helper window is opened with that input pre-filled.

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

The bundled Bend Helix V3.1 tool bends a straight two-chain nucleic-acid helix at one or more selected phosphorus residues. With one pivot it treats the helix as two rigid pieces: piece #1 stays fixed, while piece #2 is moved by an optional pivot shift, a beta bend, and an optional tau twist. With several pivots it bends in pivot order, each bend acting on the helix the previous bend generated; see [Bending at several pivots](#bending-at-several-pivots).

Open its GUI directly:

```bash
python3 re_helix_lib/bend_helix.py --gui
```

Run it from the command line:

```bash
python3 re_helix_lib/bend_helix.py --input straight_helix.pdb --pivot A36 --phi 0 --beta 30 --tau 0
```

Shift the movable piece before bending it:

```bash
python3 re_helix_lib/bend_helix.py --input straight_helix.pdb --pivot A36 --phi 0 --beta 30 --shift_axial 2 --shift_radial -1.5
```

Bend at two pivots, with per-pivot phi and tau and one beta for both:

```bash
python3 re_helix_lib/bend_helix.py --input straight_helix.pdb --pivot A20,A40 --phi 0,90 --beta 30 --tau -10,5
```

Useful Bend Helix options:

- `--pivot A36`: P-bearing residue that marks the border between fixed piece #1 and movable piece #2. Several pivots, such as `A20,A40`, are bent in the listed order.
- `--phi 0`: hinge direction around the helix axis, in degrees.
- `--beta 30`: bend angle for movable piece #2, in degrees.
- `--tau 10`: optional twist of movable piece #2 around its bent axis, in degrees.
- `--shift_axial 2` (`--sa`): shift of movable piece #2 along the helix axis before bending, in angstroms. Positive is forward, away from fixed piece #1.
- `--shift_radial -1.5` (`--sr`): shift of movable piece #2 along the pivot radius before bending, in angstroms. Positive is outward, away from the helix axis.
- `--axis_range A1-A35,B60-B26`: optional local-axis range for already-bent inputs.
- `--sep y`: give movable piece #2 new chain IDs in the output; with several pivots, every movable piece gets its own.
- `--origin y`: also write an origin-overlay PDB for comparing the original and transformed helix; with several pivots, it holds every helix the bends generate.

With several pivots, `--phi`, `--beta`, `--tau`, `--shift_axial`, and `--shift_radial` each take one value per pivot, in pivot order, or a single value used at every pivot.

### Output file names and provenance

Automatic output names lead with a `Pv` block naming the pivot residue, then the angle block: `model_PvA36_P0B30T0.pdb`, or `model_PvA36_P0B30T0_sep.pdb` with `--sep y`. The pivot label comes from the residue as resolved in the input, so `A36`, `36A`, and `A.36` all produce the same file, and bends pivoted at different residues no longer overwrite one another. An explicit `-o`/**Save as** path is used as given.

Several pivots are named in bending order. When every value is shared they lead one value block, as in `model_PvA20PvA40_P0B30T0.pdb`; when any value differs between pivots, each pivot carries its own block, as in `model_PvA20_P0B30Tm10_PvA40_P90B30T5.pdb`. A name that would pass 240 bytes, which many pivots with detailed values can reach, is shortened to the first and last pivot, the pivot count, and a 10-character digest of the full name, as in `model_PvA10-A140x14_5e3c1d9a0b.pdb`; the `COMMAND` record below still carries every value.

Every output PDB — the bent model, the `-ori` overlay, and each screened solution — records how it was made, using the same `REMARK 950 RE_SCRIPT` convention as the other bundled tools:

```
REMARK 950 RE_SCRIPT SOFTWARE name=bend_helix version=V3.1 developer=DiLiuLab
REMARK 950 RE_SCRIPT COMMAND text=python bend_helix.py --input A60-heli.pdb --pivot A36 --phi 0 --beta 30 --tau 0 --shift_axial 0 --shift_radial 0 --sep n --align y --origin n
REMARK 950 RE_SCRIPT OUTPUT_STAGE name=bend_helix
```

The `COMMAND` record holds the full equivalent CLI command, so a model carries the recipe that produced it and can be regenerated from the file alone. Several pivots and per-pivot values are recorded as comma lists, such as `--pivot A20,A40 --phi 0,90`, and a list starting with a minus sign is attached to its option, as in `--tau=-10,5`. When the run came from the GUI, this is the same command string the GUI prints in its result pane. Records are placed after any existing `REMARK` lines and before the first structural record, so they never land among `LINK` or coordinate records, and re-running the tool on its own output appends a fresh block rather than replacing the old one.

### Bending at several pivots

`--pivot` and the GUI **Pivot P residue(s)** field accept several residues separated by spaces, commas, or semicolons, for example `A20, A40`. The bends are applied in the order the pivots are listed:

1. Bend 1 acts on the input helix at the first pivot, exactly like a single-pivot run. Conceptually it produces two full helices: the input helix, which supplies fixed piece #1, and a rigidly bent copy of the whole helix.
2. Bend 2 acts on that generated full helix at the second pivot and produces the next full helix, and so on.
3. Finally the pieces are linked: with n pivots the duplex splits into n + 1 pieces, piece #1 is taken from the input helix and piece #k + 1 from the helix bend k generated. The unused part of every generated helix is dropped, and the pieces keep their chain IDs and residue numbering, so the output is one continuous duplex.

Because each bend acts on the helix the previous bend generated, every pivot after the first must lie in the movable piece of the bend before it, that is, after it along the helix axis. List the pivots from the fixed end toward the free end; a pivot that lies before the previous one, or two pivots on the same base pair (including the two residues of one pair), is rejected with an error that names both pivots.

Each bend is built from its own pivot's frame on the helix it acts on: the local axis, the pivot P position, and the phi = 0 radial direction all travel with that helix. Phi, beta, tau, Sa, and Sr therefore mean exactly what they mean for a single bend, measured at that pivot, and the angle between two consecutive pieces' axes is that bend's beta. A generated helix is a rigid copy of the input, so carrying the input frame along gives the frame Bend Helix would fit to the generated helix itself; the regression suite checks this against a chain built with the single-pivot code alone, step by step, to within 1e-6 A.

`--phi`, `--beta`, `--tau`, `--shift_axial`, and `--shift_radial` (and the matching GUI fields) each take either one value per pivot, in pivot order, or a single value used at every pivot. Any other count is an error. Omitted tau, Sa, and Sr are 0 at every pivot. On the command line, an option may also be repeated (`--pivot A20 --pivot A40`), and a list starting with a minus sign works as typed (`--tau -10,5`); positional arguments take lists too, but a negative list there needs the named option.

`--align y` applies to every bend, so each pivot P returns to its (shifted) position on the helix that bend acts on, and the backbone stays connected at every junction. With `--axis_range`, each pivot uses the range containing it; the ranges must all run in the same start-to-end direction along the helix so that every pivot sees the same base-pair order. `--sep y` gives every movable piece its own new chain IDs, in piece order. `--origin y` writes the input helix and every generated full helix to the `-ori` PDB, in bend order: A/B input, C/D from bend 1, E/F from bend 2, and so on. The run summary lists every piece with its base-pair range and source helix, then every bend with its own axis, hinge, twist axis, and alignment translation.

Screening still works with several pivots as long as phi, beta, tau, Sa, and Sr each hold a single value: every fixed value and every screened candidate is then applied at each pivot. Entering more than one value in any of these five fields turns screening off and clears the **Screen** checkboxes until every field holds one value again. Atom names default to the final bent model, so `A:150:P` is residue A150 where it ends up after every bend (see [Atom names](#atom-names) below); the overlay naming adds two chains per bend (A/B input, C/D bend 1, E/F bend 2, ...). **Phi-corrected pivot P** and the **Local axis range(s)** rotation axis refer to the first pivot, whose bend acts on the input helix.

### Shifting the pivot

`--shift_axial` (**Sa**) and `--shift_radial` (**Sr**) translate movable piece #2 in angstroms *before* any rotation is built, and every rotation is then rebuilt from the shifted pivot. Both default to 0, so a run without them behaves exactly as it did in earlier versions.

- **Sa**, along the axis. Positive Sa moves piece #2 forward, away from fixed piece #1, opening the junction; negative Sa moves it backward, toward piece #1. Piece #2 is by construction the axis-downstream side of the pivot, so this sign does not depend on which chain appears first in the PDB the way phi does, and it is not the 5'-to-3' direction of the pivot strand. The slide follows the local axis as it stands *before* the bend, not piece #2's own bent axis, so at a large beta the displacement will not look parallel to piece #2. With `--axis_range`, it follows the range you supplied, which for an already-bent input may describe piece #1.
- **Sr**, along the radius. Positive Sr moves piece #2 out, away from the helix axis; negative Sr moves it in. The direction is the pivot P atom's own radial direction, fixed at the phi = 0 direction: it does **not** rotate with phi, because the shift happens before the rotation. At phi = 180 a positive Sr therefore moves piece #2 toward the side opposite the hinge. The direction is taken from the P atom of the residue you named, so naming the other residue of the same pivot base pair gives the same two pieces but reverses Sr.

Sr also changes the radius of the circle that phi sweeps: the hinge sits at distance `radius + Sr` from the axis. `Sr = -radius` places the hinge exactly on the helix axis, turning beta into a pure axial kink, and values below that push the hinge through to the far side. Both are allowed, and neither reverses the sense of beta. The run log reports the original radius, the shifted radius, the shift vector, and the shifted pivot position.

One consequence is worth knowing before you use the shifts:

> With `--align y` (the default), the net effect of Sa and Sr is exactly to translate the bent piece #2 by the shift vector. The bend and twist keep the same shape and direction they would have without a shift. The rebuilt hinge position changes the result only with `--align n`.

This follows from what `--align y` does: it puts the pivot P atom back on a fixed target, which cancels every hinge *location* from the result and leaves only the alignment target itself. The target is the *shifted* pivot; aligning back to the unshifted pivot instead would cancel Sa and Sr entirely and make the two settings do nothing.

Automatic output names gain an `Sa...Sr...` block after the angle block whenever either shift is nonzero, for example `model_PvA36_P0B30T0Sa2Srm1p5.pdb`. A run with both shifts at 0 omits that block.

### Screening in the GUI

Each phi, beta, tau, **Sa**, and **Sr** field has a **Screen** checkbox. Check exactly one or two of the five, then click **Screening to achieve...** to define the target and the candidate grid. An unchecked variable is held at the value in its main-window field. For each checked variable, provide **From**, **To**, and **Step** values: the From and To values are included in the coarse search, and Step must be positive. Each row uses its own unit, shown in its row label — degrees for phi, beta, and tau; angstroms for Sa and Sr — so a mixed pair such as beta and Sr is screened in degrees and angstroms respectively. The default coarse ranges are phi −90° to 90°, beta −180° to 180°, tau −180° to 180°, and −10 Å to 10 Å for both shifts; Step defaults to 6° for an angle and 0.5 Å for a shift. With two checked variables, Bend Helix evaluates their Cartesian-product coarse grid. A live preview beneath the fields shows each calculated coarse-grid value (compacted for long grids), the value count for each variable, and the total number of coarse candidates.

After the coarse pass, Bend Helix identifies every local coarse minimum as a promising region and refines each region independently. It adjusts one or both screened variables with an adaptive pattern search, halves the refinement spacing when no improvement is found, and stops at 0.001-degree or 0.001-angstrom precision. A region whose residual is still above the solution tolerance at that precision keeps halving, down to 1e-6 degree or angstrom, until it meets the tolerance; this matters when one screened value is shared by several pivots, because the target then moves several times as fast as with one bend, and a 0.001-degree step can leave a residual just above a 0.001-degree tolerance. Thus, From = 0, To = 10, Step = 4 starts with 0, 4, 8, and 10 degrees but can select an in-between result such as 6.35 degrees. Refining the local regions independently allows separate solution branches to be found without constructing an extremely dense full grid.

**Solution tolerance** determines which distinct results are reported. Its unit is angstroms in distance mode and degrees in rotation mode, and its default is 0.001. Qualifying solutions are de-duplicated at each screened variable's own refinement precision and sorted by residual, then phi, beta, tau, Sa, and Sr. If none reaches the requested tolerance, the closest candidate is reported explicitly as a fallback. The result pane and run log show a table containing every reported solution; the table gains **shift axial** and **shift radial** columns when a shift is in play.

By default, only the best solution and its origin-overlay PDB are written. Enable **Write all reported solutions** to write every additional solution under numbered names such as `model_P10Bm20T30_scr_sol002.pdb`; each additional model also receives a corresponding `-ori.pdb` overlay. All reported solutions use the selected `--align y` or `--align n` behavior during screening and output generation.

Every argument or source choice in the screening popup has a round **?** button that opens a focused explanation and example.

The screening window provides two target modes:

- **Screening for distance** minimizes the absolute difference from a requested distance in angstroms. Define the two endpoints as either two atoms or one atom plus one XYZ point.
- **Screening for rotation** minimizes the circular difference from a requested signed angle in degrees around a defined axis. Define the endpoints as two atoms, one atom plus one XYZ point, or one atom plus the phi-corrected pivot P position. The latter is the pivot P position after applying the candidate's pivot shifts and phi correction, so it is recalculated for every candidate; note that it degenerates to a single phi-independent point when `radius + Sr` reaches 0.

#### Atom names

Atoms use `CHAIN:RESIDUE:ATOM` syntax, for example `A:36:P` or `B:24:O5'`. The **Atom names** choice at the top of the screening window sets what the chain IDs mean:

- **Final bent model** (default): an atom is named where it ends up in the written model, by that model's chain IDs and residue numbers. With pivots A118 and A149, `A:150:P` is residue A150 after both bends have moved it, and `A:91:P` is A91, which stays in fixed piece #1. With `--sep y`, the moved pieces are written under new chain IDs, and atoms use those IDs: there, A150 is `E:150:P`, and `A:150:P` is refused with the right name.
- **Origin overlay helices**: names follow the `-ori` overlay, which holds full helices. A/B is the unmoved input helix, C/D the full helix bend 1 generates, E/F the full helix bend 2 generates, and so on, so a residue follows the bent model only on the helix that supplies its piece: with pivots A118 and A149, A150 is `E:150:P`, while `A:150:P` stays at its input position. Use these names to measure against the unmoved input, for example how far a residue moved (`A:36:P` to `C:36:P` with one pivot).

The text under the choice lists, for the current input and pivots, which residues each piece holds and which chain names reach them, for example `#3 A149-A235, B1-B87`.

A screen whose target cannot change is refused before it starts: when both endpoints, and for rotation the axis, all stay on the unmoved input or all move rigidly with one helix, every candidate would give the same value and the screen could only report its first grid value. The error names the selector that would follow the bent model, for example `use E:150:P, or name it in the final model as A:150:P`. A one-candidate grid only measures and is not refused.

Both **Endpoint 1 atom(s)** and **Endpoint 2 atom(s)** accept one or more atoms, separated by commas or semicolons. For example, enter `A:36:P, A:37:P, A:38:P` for endpoint 1 and `A:40:P, A:41:P, A:42:P` for endpoint 2. Each endpoint uses the arithmetic mean of its atoms' XYZ coordinates (equal weight per atom), recalculated after resolving every atom in each candidate's geometry. This works for both distance and rotation screening; the two endpoints may contain different numbers of atoms, and a group may span chains and pieces. Single-atom entries and endpoint 2's XYZ/pivot choices continue to work as before.

For rotation screening, the axis can be supplied geometrically with the same source choices used by the main tool's restrained-rotation mode: an XYZ point or atom for the axis point, together with a direct vector, two XYZ points, two atoms, or the right-hand normal to two vectors for its direction. Alternatively, choose **Local axis range(s)** and enter the range definitions directly in the screening popup; multiple definitions can be separated by semicolons. These popup-owned ranges override the main-window Local axis range(s) for candidate generation and for the best/additional screened outputs, so the evaluated and written geometries remain identical without copying values from the main window.

The measured rotation is signed from endpoint 1 toward endpoint 2 by the right-hand rule about the positive axis direction. Angular differences wrap across -180/180 degrees, so equivalent directions near the wrap boundary compare correctly. Bend Helix reports every distinct coarse or refined solution within tolerance; an exact target is not required because the closest fallback is retained when necessary. Automatic screening output names retain the pivot block and the selected P/B/T angle values, add the `Sa...Sr...` block when a shift is nonzero, and add `_scr`, for example `model_PvA36_P0B30T0_scr.pdb`, `model_PvA36_P0B30T0_scr_sep.pdb`, or `model_PvA36_P0B30T0Sa2Sr0_scr.pdb`; an explicit **Save as** path is honored without automatically adding `_scr`. A screening run automatically writes the best model's origin-overlay PDB even when origin-overlay output was not otherwise selected, so its name inherits `_scr` as well (for example, `model_PvA36_P0B30T0_scr-ori.pdb`).

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

A standard 3'-5' P--O3' LINK that Phenix does not link by itself, such as a strand closure or a link between non-consecutive residues, gets the restraints Phenix applies to every ordinary backbone link (its `rna3p` chain link): the O3'-P bond, the O3'-P-O5'/OP1/OP2 angles, the O3'-side C3'-O3'-P angle, and the epsilon, zeta and alpha torsions, all with Phenix's ideals and sigmas. A closure therefore relaxes like the rest of the backbone. LINKs between consecutive residues are left to Phenix, and LINKs involving a nonstandard residue or a missing backbone atom keep the bond and P-centered angles only.

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
- `--legacy-link-restraints`: write the V1.0 restraints for standard 3'-5' backbone LINKs: the bond and P-centered angles only, without the C3'-O3'-P angle and torsions. The GUI has the same switch under **Advanced restraint values**.

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

## Build Heterohelix Junction Tool

The bundled Build Heterohelix Junction tool (V1.0) joins two or more duplex segments end to end and builds the junctions between them. A segment can be any duplex form: B-DNA, A-RNA, Z-DNA or another, with D or L sugars, written as BNP_NA writes duplexes (two chains, strand A running 5'->3' from the bottom end to the top end, strand B back down, strand A's i-th residue paired with strand B's i-th from the other end). Each segment stays a rigid copy of its file. Only the phosphates of the junctions are built, and the O5' atoms they bond are turned about their C4'-C5' bond.

The fitter started as Junction_Fit V1.0, written for the AZBMOST B-Z bowtie and L-Z junction examples. It was fixed on 2026-11-26: a junction is now fitted the same way whichever end is read first and wherever its files sit in space (the search score counts contacts from both sides and samples its quick phosphates in frames attached to the atoms; the real phosphates are minimized in such frames from many starts), every junction is searched from both of its ends with a second pass, and equal junctions share one fit.

### The chain sense rule

Every base pair gets a frame of three axes, much like the standard base-pair frame (Olson et al. 2001): x runs from the midpoint of the two C1' atoms towards the bases (the major-groove edge of the pair), y runs towards the C1' of one chosen strand, called strand I, and z = x × y stands perpendicular to the pair. At the end of a segment, strand I is the strand whose 3' end is there, so that strand I runs 5'->3' out of the segment; the end's chain sense is +1 when z also points out of the segment, and -1 when it points back in. In B-DNA and A-RNA, z points out of the segment: chain sense +1. In Z-DNA, z points back into the segment: chain sense -1. A mirror image flips the sign. So L-DNA and L-RNA, the mirror images of B-DNA and A-RNA, are -1, and L-Z-DNA, the mirror image of Z-DNA, is +1. The chain senses of the two ends decide each junction:

- matched (B|B, Z|Z, L-DNA|Z-DNA): a **continuous** junction, two ordinary 3'-5' links. Both strands run on 5'->3'.
- opposite (B|Z, D-DNA|L-DNA): a **bowtie** junction, one 3'-3' link through an X33 phosphate and one 5'-5' link through one phosphate. A strand that ran straight on would meet the next base pairs from their other face (the chain sense paradox), so the strands turn back instead.

`--link bowtie` or `--link continuous` can name the kind, but only the kind the chain senses allow. The other kind would turn one base pair over, so it is refused before any fitting.

The chain sense is checked before fitting and after building. Every result goes into the run log and the JSON report.

- Before fitting, a `Chain sense check` block lists every junction with both of its ends (segment number and file, end, 3' and 5' residues, sugars, chain sense) and the result, such as `matched -> continuous`.
- Before fitting, every segment must have the same chain sense all along. Each base pair gets its frame with y towards strand A, and the z axes of neighbouring pairs may not turn by more than 90 degrees. A segment with a pair turned over is refused, and the message names the pairs. A segment also needs at least two base pairs, since one pair alone has no direction.
- Before fitting, the residues each link joins must have the atoms the link bonds: O3' and C3' at a 3' end, C4', C5' and O5' at a 5' end. A missing atom is named.
- After building, each output chain is walked. Every paired nucleotide (X33 skipped) gets its base-pair frame with y towards that chain's residue. Two things are measured. First, the turn between the z axes of neighbouring pairs, inside the segments and across each junction. In a correct model every pair faces the same way along a chain, so every turn stays small; a turn above 90 degrees is an error (`chain sense not matched at ...`). Second, the cos between each z axis and the helix axis of its segment, taken in the chain's file order. It shows which way each chain runs through its pairs: about +1 in every segment of a B-Z bowtie model, and about -1 in every segment of an L-Z continuous model. A chain whose pairs face along it in one segment and against it in another is an error too.

### How a junction is fitted

1. The two end pairs get frames whose z axes both point from the lower segment into the upper one.
2. The upper end is placed on the lower one by one base-pair step: tilt, roll, twist, shift, slide and rise.
3. A step is scored on heavy atoms closer than 2.9 A across the junction, on the gaps the two links must span (each with a quick phosphate), on how parallel and how far apart (3.5 A) the two end pairs are, and on how much their bases overlap.
4. Powell searches start from a grid of twists (-180 to 160 degrees in 20-degree steps), shifts and slides (-2, 0 and 2 A): 162 starts, and every junction is searched from both of its ends. A second pass starts again from eight small fixed nudges around each of the 12 best different results.
5. The 16 best different results get real phosphates, minimized in frames attached to the atoms from many starts, and the one with the lowest stacking score plus phosphate energies is kept.

Junctions of one model that join the same two end pairs, in the same order or turned 180 degrees about the duplex dyad (for example junctions 1 and 2 of L21-Z24-L21-Z24), are fitted together: their searches are pooled and all of them get the step with the lowest summed total, so equal junctions come out equal. `--separate` fits each of them on its own. A junction that repeats another one with the same files and the same ends is fitted once.

### Using it

Open it from the `Other tools` list of the main window, which passes the current input PDB as the first segment, or directly:

```bash
python3 re_helix_lib/build_heterohelix_junction.py --gui
python3 re_helix_lib/build_heterohelix_junction.py --gui L21.pdb Z24.pdb
```

The window asks only for what the tool cannot work out itself:

- **Segments, bottom to top**: one row per segment, with the PDB file, an `upside down` box (the junction below it then meets its top end) and, once a file is chosen, a grey line saying what was read: residues per strand, D or L sugars, the chain sense at the bottom and at the top, and whether every base pair faces the same way. A file that cannot be used says why, in red. `Add segment`, `Remove`, `Move up` and `Move down` act on the selected row.
- **Output**: the output PDB, by default `<first segment>_heterohelix.pdb` beside the first segment. A name without a folder goes beside the first segment too, and an empty field goes back to the default. The report file is shown grey under it, because it follows the output name.
- **Junctions**: grey rows worked out from the segments, such as `Junction 1: segment 1 (L21.pdb) top -> segment 2 (Z24.pdb) bottom; chain senses -1 and -1: matched -> continuous (two ordinary 3'-5' links)`, and which junctions are fitted together. The kind of each junction is not an input.
- **Options**: the rise aimed at (3.5 A), the number of processes (all but one), `Fit junctions that join the same end pairs together` (on), and `Reuse fits from an earlier report` (off). When reuse is on, the earlier report becomes an input: it starts as the output's own report, and `Browse...` picks another one.

`Run` starts the command line in its own process, shows the exact command first and streams the output into the run log; `Stop` ends it. When a run fails, the footer shows the error with its first reason. The `?` buttons explain each part, and `How it works` explains the whole tool.

From the command line:

```bash
python3 re_helix_lib/build_heterohelix_junction.py L21.pdb Z24.pdb
python3 re_helix_lib/build_heterohelix_junction.py B21.pdb Z24.pdb B21.pdb Z24.pdb -o B21_Z24_B21_Z24.pdb
python3 re_helix_lib/build_heterohelix_junction.py D.pdb L.pdb:flip -o D_L.pdb
```

The first writes `L21_heterohelix.pdb` beside `L21.pdb`. In the last, `:flip` turns `L.pdb` upside down, so the top of D meets the top of L.

Useful Build Heterohelix Junction options:

- `-o` or `--out`: the joined model. Default: `<first segment stem>_heterohelix.pdb` in the first segment's folder.
- `--link auto|bowtie|continuous`: one value, or one per junction separated by commas. Default `auto`; a kind the chain senses do not allow is refused before fitting.
- `--cores N`: processes for the searches. Default: all but one.
- `--rise A`: the stacking distance aimed at. Default `3.5`; within 0.15 A of it costs nothing.
- `--separate`: fit every junction on its own, also those that join the same two end pairs.
- `--reuse`: take the fits from this output's own report (the `--report` file, or the default one) instead of fitting again. This works only when that report holds fits of exactly the same junctions, made the same way: by this version of the tool, from the same segment files (unchanged), on the same search grid, with the same junctions fitted together. Use it to write a model again under the same name. Otherwise the junctions are fitted again, and the log says why.
- `--reuse-from REPORT`: the same, with the fits taken from an earlier report, for example the report of the same build under another output name. The new report is still written to the output's own report path.
- `--report FILE`: the JSON report. Default: `tmp_file/<output stem>_heterohelix_junction.json` beside the output. For the default output, `<first segment>_heterohelix.pdb`, it is `tmp_file/<first segment>_heterohelix_junction.json`.
- `--gui`: open the window; segments given on the command line fill its rows.
- `-v` or `--version`: show the tool version.

### Outputs

- The joined model, with chains A, B, ... numbered from 1. `LINK` records are written for every bond that is not an ordinary 5'->3' step in file order: X33 links, 5'-5' links, and the steps of a strand that runs through a segment against its file order after a bowtie. `REMARK 900` lines record the segments, each junction's step and the L-sugar residues. `REMARK 950 RE_SCRIPT` lines record the tool (`SOFTWARE`, `COMMAND`, `OUTPUT_STAGE`) and one `JUNCTION` line per link: the 3' and 5' residues of each 3'-5' link, the two 3' ends with the X33 between them, and the two 5' ends. Get Phenix Restraints reads these lines and frees the sugars of those residues in its `_junctions.params`.
- The JSON report: every search, the candidates with real phosphates, the chosen steps, the phosphate geometry (and the backbone torsions of continuous links), the helical axes at each junction, the bond and contact checks, the chain sense checks before and after building, and what the fits were made from (tool version, search grid, junctions fitted together, and a checksum of each segment file), which `--reuse` compares.

The run log gives one line per junction (step, stacking, phosphates, axis bend), the chain sense check after building, and the bond and contact summary of the written model.

### Run time and limits

- About 5 minutes per junction on 9 processes: every junction is searched from both of its ends (324 searches), then the second pass (up to 96 more); junctions fitted together pool their searches, so each one still adds its own. A junction that repeats another one (same files, same ends) is fitted once.
- The score changes slowly with twist, so read a fitted twist with care; the report lists every candidate. For the L21|Z24 continuous junction the best candidate (twist -30.8 degrees) is more than 1.0 below every other one; for the B21|Z24 bowtie a second candidate lies within 1.0 of the best, 4 degrees of twist away (-5.9 against -10.1).
- Fitted from either end, an L21|Z24 continuous junction gave the same step to 0.001 degrees, and with the L21 file turned and moved in space to 0.004 degrees and 0.0004 A. A B21|Z24 bowtie junction gave the same step from either end; with the B21 file turned in space it moved by 1.0 degree and 0.22 A, because the turned file's coordinates are rounded to 0.001 A and the bowtie score is very flat near its best step.
- The search can miss a better step. Fitted on its own, the B21 bottom | Z24 top bowtie junction reached a total of 6.98 (twist -6.3 degrees), although the step it shares with the B21 top | Z24 bottom junction in the B21-Z24-B21-Z24 model gives it 6.28 (twist -10.1 degrees): only one of its 420 searches ended near that step. The quick 5'-5' phosphate of the search score is a rough guess (one phosphate 1 A out from the two C5'), so the step with the best real energy does not always sit at a minimum of the search score. Junctions fitted together pool their searches, which is what found the better step here; when a junction matters, fit its equal partners with it or compare runs.
- Models are as built, not minimized.
- L residues keep the D residue names, as BNP_NA writes them; the `REMARK 900` lines list them. Phenix restrains those names with D chirality, so L residues would need L restraints before minimization.
- Before minimizing a bowtie model, check `link_distance_cutoff` in the params files. In a B21-Z24-B21-Z24 bowtie model, 6.5 A (Get Phenix Restraints) or 7.0 A (`min_P_C5.params`) let Phenix bond an O3' to the P of the next residue in file order across each 5'-5' junction, 5.5 A apart; 3.0 A avoided it.
- Needs numpy, scipy and shapely. Without them the tool stops with a one-line message (and a message box when it was opened as a window).

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

## Crossovers On Faces (--xo_faces)

A strand that crosses from one helix to another comes in on one face of the plane through the two helix axes. At a
parallel crossover it also leaves on that face, so it should stay on that face at the crossover: the two crossing
strands pass one above the other (Han et al. 2017, Science 358, eaao2648, Fig. 1B). At an antiparallel crossover it
leaves on the other face and goes through the plane by design. Because a reciprocal exchange rewires bonds without
moving atoms, at a parallel crossover one atom of each crossing strand, the leaving O3' or the entering phosphate, is
left 0.1-2 A on the other side of the plane, where the two crossing strands meet.

`--xo_faces` (GUI: **Crossovers on faces**, off by default) traces every crossing strand after the exchange. *In* is
its phosphate one step before the crossing and *out* its phosphate one step after; when both lie on the same face, the
atoms of the crossing are put on that face (at least 0.5 A from the plane; `--xo_min_face` changes this):

- the entering phosphate (O5', P, OP1, OP2) is turned about C4'-C5' (gamma) and C5'-O5' (beta);
- when that is not enough, the leaving sugar is turned a little about its glycosidic bond, or the leaving base pair a
  few degrees about the helix axis, and the phosphates next to it are turned to keep their bonds at the model's own
  O3'-P length and angles.

Every move is a rigid turn about a bond or the helix axis, so bond lengths and angles inside each nucleotide are kept
and bases stay paired. Crossover links are never lengthened (by more than 0.3 A) and are not closed; moved atoms keep
2.5 A from every atom more than three bonds away. Antiparallel crossings are left as built. A crossover is placed
whole or left exactly as built, and every crossing strand is reported in the run log and in
`REMARK 950 RE_SCRIPT XO_FACES` lines. Each input chain must lie in one helix, as in an aligned model.

Crossovers whose leaving O3' sits on the wrong face would need a local change of helix twist over several base
pairs; the option does not make that change and leaves them as built (on the lab's A-RNA models: the Han-type sheet
34 of 36 crossing strands placed, DPE 4 of 4, DPOW, DPON and the 7:4 PX none; B-DNA DPON 4 of 4, DPE 2 of 4, DPOW
none). Minimize afterwards (Get Phenix Restraints, Phenix) when bond geometry matters.

On its own:

```bash
python3 re_helix_lib/crossover_faces.py model_aligned_rex.pdb --reference model_aligned.pdb --helix AB --helix CD
```

writes `model_aligned_rex_faces.pdb` and prints one line per crossing strand.

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
- `--xo_faces` or `--xo-faces`: after the exchange, keep each strand of a parallel crossover on its own face (see Crossovers On Faces). Default off.
- `--xo_min_face 0.5`: with `--xo_faces`, how far every crossing atom must sit on its face (A).
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
