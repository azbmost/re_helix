# re_helix

`re_helix` is AZBMOST Package Module #2: Align Helices and Performing Reciprocal Exchanges.

It aligns nucleic-acid helices from reciprocal-exchange-style P-atom pairs and then applies reciprocal exchanges to the aligned structure. It can also run reciprocal exchange only, without alignment.

Current version: V3.8

## Contents

- `re_helix.py`: main CLI/GUI script.
- `re_helix_lib/`: helper modules for PDB parsing, LINK records, and reciprocal-exchange graph handling.
- `re_helix_lib/bend_helix.py`: bundled Bend Helix tool for bending a straight two-chain helix.
- `re_helix_lib/do_symmetry.py`: bundled Do Symmetry tool for averaging a pseudosymmetric assembly into an idealized symmetric PDB.
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

In the GUI, use the `Other tools` area to open bundled helper tools. `Bend Helix` opens the helix-bending GUI, and `Do Symmetry` opens the symmetry-averaging GUI. If an input PDB is already selected in `re_helix`, the helper window is opened with that input pre-filled.

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

## What Reciprocal Exchange Means

In this package, a reciprocal exchange is a virtual topology edit on a PDB model. The operation identifies residues on two DNA strands or helices, cuts the original backbone graph at the specified sites, and reconnects the graph so the strand continuities are exchanged. This is a design operation, not an enzymatic simulation: it is meant to help build the intended crossover, junction, or bowtie connectivity before later structural refinement, sequence design, synthesis, or visualization.

The idea sits in the design tradition introduced by Nadrian C. Seeman, who founded structural DNA nanotechnology by treating DNA as a programmable construction material rather than only as genetic information. Seeman's key move was to use designed sequence asymmetry to make immobile branched junctions, avoiding the branch migration of natural Holliday junctions, so junctions could serve as predictable vertices for DNA objects, arrays, and lattices. A reciprocal-exchange operation is useful for that style of design because it gives a compact way to say: "these helices meet here, and their backbone routes trade partners here." In practice, that lets a designer specify the intended strand routing and junction topology without manually rebuilding every atom record, residue number, chain break, TER record, and LINK record.

`re_helix` supports three exchange kinds:

- `double`: exchange both local backbone continuities between two specified residues.
- `single`: exchange one strand-continuity relationship while leaving the complementary local relationship unchanged.
- `bowtie`: create paired 3'-3' and 5'-5' junction behavior, including LINK records and X33 linker-phosphate records where needed.

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
