# re_helix

`re_helix` is AZBMOST Package Module #2: Align Helices and Performing Reciprocal Exchanges.

It aligns nucleic-acid helices from reciprocal-exchange-style P-atom pairs and then applies reciprocal exchanges to the aligned structure. It can also run reciprocal exchange only, without alignment.

Current version: V3.16

## Contents

- `re_helix.py`: main CLI/GUI script.
- `re_helix_lib/`: helper modules for PDB parsing, LINK records, and reciprocal-exchange graph handling.
- `re_helix_lib/bend_helix.py`: bundled Bend Helix tool for bending a straight two-chain helix.
- `re_helix_lib/do_symmetry.py`: bundled Do Symmetry tool for averaging a pseudosymmetric assembly into an idealized symmetric PDB.
- `re_helix_lib/add_pdb_link_record.py`: bundled Add PDB LINK Record tool for staging P/O3' LINK records and rebuilding chain topology.
- `re_helix_lib/insert_virtual_resi.py`: bundled Insert Virtual Resi tool for inserting residue-numbering gaps and updating LINK endpoints.
- `re_helix_lib/generate_lattice.py`: bundled Generate Lattice tool for writing a P1 CRYST1 lattice record from three lattice vectors.
- `re_helix_lib/get_phenix_restraints.py`: bundled Get Phenix Restraints tool for converting LINK records into Phenix geometry restraints, optional junction movement-selection params, and linker support files.
- `assets/icon.png`: optional GUI/task-menu icon. The script uses it when present and falls back to the default Tk icon when it is missing.

## Requirements

- Python 3.9 or newer is recommended.
- No third-party Python packages are required for the command-line workflow.
- Tkinter is required only for GUI mode. Most python.org and system Python installs include it.

## Quick Start

Launch the GUI:

```bash
python3 re_helix.py
```

In the GUI, use the `Other tools` area to open bundled helper tools. `Bend Helix` opens the helix-bending GUI, `Do Symmetry` opens the symmetry-averaging GUI, `Add PDB LINK Record` opens the LINK-record/topology helper, `Insert Virtual Resi` opens the residue-renumbering helper, `Generate Lattice` opens the P1 lattice/CRYST1 helper, and `Get Phenix Restraints` opens the Phenix restraint-generation helper. If an input PDB is already selected in `re_helix`, the helper window is opened with that input pre-filled.

When the input PDB changes, the GUI updates the default `Output base` automatically unless that field has been changed to a custom value. For large exchange specifications, the `CLI pair args` field below the pair rows can be filled with the same concatenated pair tokens used on the command line; when it is filled, the individual pair rows are ignored.

The main `re_helix` run log also mirrors stdout/stderr from bundled tools launched through `Other tools`, so equivalent CLI commands, selected LINK summaries, completion messages, and errors remain visible in the main window even when the helper window has no log box.

Run alignment plus reciprocal exchange from the command line:

```bash
python3 re_helix.py input.pdb '(AB)' '(CD)' 30A 8D d 13B 24C s -o model
```

This writes:

- `model_aligned.pdb`: aligned structure before reciprocal exchange.
- `model_aligned_rex.pdb`: aligned structure after reciprocal exchange.

Run reciprocal exchange only, with no alignment:

```bash
python3 re_helix.py input.pdb 9C 23A d 23C 23F b --re_only -o model
```

This writes:

- `model_rex.pdb`: reciprocal-exchanged structure generated directly from `input.pdb`.

## Bend Helix Tool

The bundled Bend Helix tool bends a straight two-chain nucleic-acid helix at a selected phosphorus residue. It treats the helix as two rigid pieces: piece #1 stays fixed, while piece #2 is moved by a beta bend and optional tau twist.

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

The bundled Add PDB LINK Record tool helps create PDB `LINK` records between phosphate `P` atoms and `O3'`/`O3*`/`O3` atoms. In GUI mode, it can automatically stage terminal-chain circularization links and manually stage internal or inter-chain P/O3' links, then rebuild the chain topology, chain IDs, TER records, residue numbering, and LINK records in one pass. Existing input `LINK` records are preserved and remapped to the rebuilt chain/residue labels when possible.

Open its GUI directly:

```bash
python3 re_helix_lib/add_pdb_link_record.py --gui
```

Run automatic chain circularization from the command line:

```bash
python3 re_helix_lib/add_pdb_link_record.py input.pdb --chains A B -o input_linked.pdb
```

If `--chains` is omitted, the command-line mode attempts automatic circularization for every chain with usable terminal `P` and `O3'` atoms. Without `-o`, the default output inserts `_circ` before the input extension, for example `input.pdb` becomes `input_circ.pdb`.

Useful Add PDB LINK Record options:

- `--gui`: open the GUI for automatic and manual LINK staging.
- `--chains A B` or `--chains A,B`: choose chains for automatic terminal circularization.
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

Useful Insert Virtual Resi options:

- `--insert A55 3`: add a gap of 3 residue numbers after residue 55 in chain A.
- Repeat `--insert` for more chains or residue positions.
- `-o output.pdb`: choose the output file. Without `-o`, the default output inserts `_vresi` before the input extension.
- `-v` or `--version`: show the bundled tool version.

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

The bundled Get Phenix Restraints tool is based on `link_to_geometry_restraintsV4.py`. It converts PDB `LINK` records into a `*_links.params` file containing Phenix `geometry_restraints.edits` bond restraints and phosphate-centered angle restraints. For standalone 3'-to-3' linker phosphate residues, it also writes internal P-OP1/P-OP2 bond and OP1-P-OP2 angle restraints.

Open its GUI directly:

```bash
python3 re_helix_lib/get_phenix_restraints.py --gui
```

Run it from the command line:

```bash
python3 re_helix_lib/get_phenix_restraints.py model_rex.pdb --output-base model_rex
```

By default, this writes `model_rex_links.params`, `model_rex_junctions.params` when `REMARK 950 RE_SCRIPT JUNCTION` lines are present, and nonstandard-linker support files such as `X33_phenix_atomtypes.cif` and `X33_safe_interpretation.params`.

Recommended Phenix command:

```bash
phenix.geometry_minimization \
  model_rex.pdb \
  model_rex_links.params \
  model_rex_junctions.params \
  X33_phenix_atomtypes.cif \
  X33_safe_interpretation.params
```

Useful Get Phenix Restraints options:

- `--linker-resname X33`: choose the standalone 3'-to-3' linker phosphate residue name. The default is `X33`; custom nonstandard names get matching CIF/safe params files.
- `--link-distance-cutoff 6.5`: choose the generated `pdb_interpretation.link_distance_cutoff` value. The default is `6.5`.
- `--no-junctions-params`: skip `*_junctions.params`. If you do this, use exactly one movement-selection file such as `min_P_C5.params` or a carefully prepared `min.params`.
- `--no-linker-support-files`: skip writing `<resname>_phenix_atomtypes.cif` and `<resname>_safe_interpretation.params`.
- `--include-phenix-builtin-angles`: diagnostic mode that can reproduce older duplicate-prone angle output.

Use exactly one movement-selection file for `phenix.geometry_minimization`. Usually this should be `*_junctions.params`; do not combine it with `min_P_C5.params` or `min.params` unless you deliberately want to test which top-level `selection = ...` Phenix uses.

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

Each exchange is:

```text
<pos1> <pos2> <kind>
```

Accepted kinds:

- `d` or `double`
- `s` or `single`
- `b` or `bowtie`

For alignment mode, a single-site inter-helix pair can also include a fixed rho angle when `--axis_parallel n` is used:

```text
<pos1> <pos2> <rho_deg> <kind>
```

Example:

```bash
python3 re_helix.py input.pdb '(AB)' '(CD)' 26A 9C 90 d --axis_parallel n -o angled_model
```

## Common Options

- `-o, --output`: output base path. A `.pdb` suffix is stripped before output suffixes are added.
- `--gui`: launch the Tk GUI explicitly.
- `-v, --version`: show the app version and exit.
- `--re_only` or `--re-only`: apply reciprocal exchange only and write `<base>_rex.pdb`.
- `--axis_dist 22.0`: target helix-axis distance in angstroms during alignment.
- `--axis_parallel y|n`: keep axes parallel (`y`) or allow a rho tilt (`n`).
- `--axis_range B26-B60,A1-A35`: define residue windows for helical-axis estimation. Repeat as needed.
- `--fix A`: keep the helix containing chain `A` fixed during alignment.
- `--replicate`: replicate the full input chain set before alignment or RE-only processing.
- `--cir_shift 8`: choose the residue shift used when writing circular reciprocal-exchange strands.
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
