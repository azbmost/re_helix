#!/usr/bin/env python3
"""check_pdb_clashes.py

Count close heavy-atom contacts (clashes) in a PDB file.

The check is intentionally geometric and structure agnostic. It works on
nucleic acids, proteins, ligands, and mixed assemblies, and it makes no
assumption about chain length, residue numbering, or composition.

Pairs are excluded from the clash count when they are not really independent
contacts:

- the two atoms belong to the same residue;
- the two atoms are alternate conformations of the same site (different,
  non-blank altLoc labels);
- the two atoms sit in the same chain within ``--adjacent-window`` residues of
  each other, so ordinary covalent neighbors are not counted;
- the two residues are joined by a ``LINK`` record. re_helix builds topology
  with LINK records for cyclization, reciprocal exchange, and permuted chains,
  where residue numbering no longer tracks real connectivity, so a junction
  would otherwise be reported as a false clash;
- the two residues close a chain that ``--circular`` marks as cyclic.

Everything that narrows the check is opt-in or derived from the file itself.
With default options the tool reports every heavy-atom pair closer than the
cutoff that is not a same-residue, alternate-conformation, sequence-adjacent,
or LINK-bonded contact.
"""

from __future__ import annotations

import argparse
import math
import shlex
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

try:
    from re_helix_lib.gui_icon import apply_optional_icon
except ImportError:  # pragma: no cover - direct script execution fallback
    from gui_icon import apply_optional_icon

try:
    from re_helix_lib.get_phenix_restraints import parse_link_line
except ImportError:  # pragma: no cover - direct script execution fallback
    from get_phenix_restraints import parse_link_line

try:  # pragma: no cover - exercised only where SciPy is installed.
    import numpy as np
    from scipy.spatial import cKDTree
except Exception:  # pragma: no cover - fallback is for non-SciPy environments.
    np = None  # type: ignore[assignment]
    cKDTree = None  # type: ignore[assignment]

TOOL_NAME = "Check Clashes"
VERSION = "1.0"

DEFAULT_CUTOFF = 1.60
DEFAULT_ADJACENT_WINDOW = 1
DEFAULT_TOP = 10

# Chimera numbers the first opened structure #0; ChimeraX numbers it #1.
DEFAULT_CHIMERA_MODEL = 0
DEFAULT_CHIMERAX_MODEL = 1

# (chain, resSeq, iCode) identifies a residue everywhere in this module.
ResidueKey = Tuple[str, int, str]


@dataclass(frozen=True)
class PdbAtom:
    serial: int
    atom_name: str
    altloc: str
    resname: str
    chain: str
    resseq: int
    icode: str
    coord: Tuple[float, float, float]
    element: str

    @property
    def residue_key(self) -> ResidueKey:
        return (self.chain, self.resseq, self.icode)


@dataclass
class ParsedPdb:
    atoms: List[PdbAtom]
    link_pairs: Set[frozenset]
    model_numbers: List[int]
    model_used: Optional[int]
    skipped_atom_lines: int


@dataclass
class ClashReport:
    pdb_path: Path
    method: str
    atoms_total: int
    heavy_atoms: int
    cutoff: float
    adjacent_window: int
    circular_mode: str
    circular_chains: Dict[str, Tuple[int, int]]
    link_records: int
    link_excluded: int
    altloc_excluded: int
    adjacent_excluded: int
    clashes: List[Tuple[float, PdbAtom, PdbAtom]]
    composition: "Counter[str]"
    model_numbers: List[int]
    model_used: Optional[int]
    skipped_atom_lines: int


def read_pdb(path: Path, model: Optional[int] = None) -> ParsedPdb:
    """Read atoms and LINK-bonded residue pairs from a PDB file.

    ``model`` selects one MODEL from a multi-model file. When it is None the
    first MODEL is used. Files without MODEL records are read whole.
    """
    text = Path(path).read_text(errors="replace")

    model_numbers: List[int] = []
    current_model: Optional[int] = None
    atoms: List[PdbAtom] = []
    link_pairs: Set[frozenset] = set()
    skipped = 0

    for line in text.splitlines():
        record = line[:6]

        if record.startswith("MODEL"):
            try:
                current_model = int(line[10:14])
            except ValueError:
                current_model = len(model_numbers) + 1
            model_numbers.append(current_model)
            continue

        if record.startswith("ENDMDL"):
            current_model = None
            continue

        if record.startswith("LINK"):
            parsed = parse_link_line(line)
            if parsed is not None:
                first = (parsed.chain_id_1, parsed.res_seq_1, "")
                second = (parsed.chain_id_2, parsed.res_seq_2, "")
                if first != second:
                    link_pairs.add(frozenset((first, second)))
            continue

        if not record.startswith(("ATOM  ", "HETATM")):
            continue

        atom = _parse_atom_line(line, fallback_serial=len(atoms) + 1)
        if atom is None:
            skipped += 1
            continue
        atoms.append(atom)

    model_used: Optional[int] = None
    if model_numbers:
        model_used = model if model is not None else model_numbers[0]
        if model_used not in model_numbers:
            raise ValueError(
                "MODEL %d is not in %s. Available models: %s"
                % (model_used, Path(path).name, ", ".join(str(n) for n in model_numbers))
            )
        atoms, skipped = _atoms_for_model(text, model_used)
    elif model is not None:
        raise ValueError("%s has no MODEL records, so --model cannot be used." % Path(path).name)

    return ParsedPdb(
        atoms=atoms,
        link_pairs=link_pairs,
        model_numbers=model_numbers,
        model_used=model_used,
        skipped_atom_lines=skipped,
    )


def _atoms_for_model(text: str, model_number: int) -> Tuple[List[PdbAtom], int]:
    """Collect the atoms of one MODEL block, with that block's skipped count."""
    atoms: List[PdbAtom] = []
    skipped = 0
    in_target = False
    for line in text.splitlines():
        if line.startswith("MODEL"):
            try:
                current = int(line[10:14])
            except ValueError:
                current = -1
            in_target = current == model_number
            continue
        if line.startswith("ENDMDL"):
            in_target = False
            continue
        if not in_target:
            continue
        if not line[:6].startswith(("ATOM  ", "HETATM")):
            continue
        atom = _parse_atom_line(line, fallback_serial=len(atoms) + 1)
        if atom is None:
            skipped += 1
            continue
        atoms.append(atom)
    return atoms, skipped


def _parse_atom_line(line: str, fallback_serial: int) -> Optional[PdbAtom]:
    """Parse one ATOM/HETATM line, or return None when it is unusable."""
    padded = line.ljust(80)
    try:
        x = float(padded[30:38])
        y = float(padded[38:46])
        z = float(padded[46:54])
    except ValueError:
        return None
    try:
        resseq = int(padded[22:26])
    except ValueError:
        return None
    try:
        serial = int(padded[6:11])
    except ValueError:
        serial = fallback_serial
    return PdbAtom(
        serial=serial,
        atom_name=padded[12:16].strip(),
        altloc=padded[16:17].strip(),
        resname=padded[17:20].strip(),
        chain=padded[21].strip() or ".",
        resseq=resseq,
        icode=padded[26:27].strip(),
        coord=(x, y, z),
        element=padded[76:78].strip(),
    )


def is_hydrogen(atom: PdbAtom) -> bool:
    element = atom.element.upper()
    if element:
        return element in {"H", "D"}
    name = atom.atom_name.upper()
    while name and name[0].isdigit():
        name = name[1:]
    return name.startswith(("H", "D"))


def parse_residue_set(text: str) -> Set[Tuple[Optional[str], int]]:
    """Parse a residue selection such as ``11,32,53`` or ``A:11,B:32``.

    A bare number matches that residue in every chain. ``CHAIN:NUMBER``
    matches only that chain.
    """
    if text is None or text.strip().lower() in {"", "none", "off", "false"}:
        return set()
    residues: Set[Tuple[Optional[str], int]] = set()
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        if ":" in item:
            chain_text, number_text = item.split(":", 1)
            chain = chain_text.strip() or None
        else:
            chain = None
            number_text = item
        try:
            residues.add((chain, int(number_text.strip())))
        except ValueError:
            raise ValueError(
                "Residue selection entries must be NUMBER or CHAIN:NUMBER; got %r" % item
            )
    return residues


def parse_circular_option(text: str) -> Tuple[str, Optional[int]]:
    """Parse ``--circular`` into (mode, fixed_last_residue).

    Modes are ``off``, ``auto`` (detect each chain's first/last residue), and
    ``fixed`` (an explicit last residue number, applied to every chain).
    """
    value = (text or "").strip().lower()
    if value in {"", "off", "none", "false", "0"}:
        return "off", None
    if value == "auto":
        return "auto", None
    try:
        number = int(value)
    except ValueError:
        raise ValueError("--circular must be 'off', 'auto', or a residue number; got %r" % text)
    if number <= 0:
        return "off", None
    return "fixed", number


def circular_bounds(
    atoms: Sequence[PdbAtom],
    mode: str,
    fixed_last: Optional[int],
) -> Dict[str, Tuple[int, int]]:
    """Return the (first, last) residue pair that closes each cyclic chain."""
    if mode == "off":
        return {}
    bounds: Dict[str, Tuple[int, int]] = {}
    if mode == "fixed" and fixed_last is not None:
        for atom in atoms:
            bounds[atom.chain] = (1, fixed_last)
        return bounds
    for atom in atoms:
        low, high = bounds.get(atom.chain, (atom.resseq, atom.resseq))
        bounds[atom.chain] = (min(low, atom.resseq), max(high, atom.resseq))
    return {chain: pair for chain, pair in bounds.items() if pair[0] != pair[1]}


def distance(a: PdbAtom, b: PdbAtom) -> float:
    return math.sqrt(sum((a.coord[idx] - b.coord[idx]) ** 2 for idx in range(3)))


def is_alternate_conformation(a: PdbAtom, b: PdbAtom) -> bool:
    """True when the two atoms are different altLoc states of one site."""
    if not a.altloc or not b.altloc:
        return False
    return a.altloc != b.altloc


def is_linked_pair(a: PdbAtom, b: PdbAtom, link_pairs: Set[frozenset]) -> bool:
    if not link_pairs:
        return False
    return frozenset((a.residue_key, b.residue_key)) in link_pairs


def sequence_separation(a: PdbAtom, b: PdbAtom) -> int:
    """Residue separation within one chain.

    Insertion-code siblings such as 5A and 5B share a residue number but are
    consecutive residues, so they count as separation 1 rather than 0. That
    keeps them excluded at the default window while still letting
    ``--adjacent-window 0`` compare every non-identical residue pair.
    """
    separation = abs(a.resseq - b.resseq)
    if separation == 0 and a.icode != b.icode:
        return 1
    return separation


def is_sequence_adjacent(
    a: PdbAtom,
    b: PdbAtom,
    adjacent_window: int,
    circular: Dict[str, Tuple[int, int]],
) -> bool:
    if a.chain != b.chain:
        return False
    if sequence_separation(a, b) <= adjacent_window:
        return True
    bounds = circular.get(a.chain)
    if bounds is not None and {a.resseq, b.resseq} == {bounds[0], bounds[1]}:
        return True
    return False


def close_pairs_with_method(
    atoms: Sequence[PdbAtom], cutoff: float
) -> Tuple[List[Tuple[int, int]], str]:
    if cKDTree is not None and np is not None and atoms:
        coords = np.array([atom.coord for atom in atoms], dtype=float)
        pairs = [(int(i), int(j)) for i, j in cKDTree(coords).query_pairs(cutoff)]
        return pairs, "scipy.spatial.cKDTree"

    pairs: List[Tuple[int, int]] = []
    cutoff_sq = cutoff * cutoff
    for i, atom_i in enumerate(atoms):
        xi, yi, zi = atom_i.coord
        for j in range(i + 1, len(atoms)):
            xj, yj, zj = atoms[j].coord
            dist_sq = (xi - xj) ** 2 + (yi - yj) ** 2 + (zi - zj) ** 2
            if dist_sq < cutoff_sq:
                pairs.append((i, j))
    return pairs, "pure Python pair scan"


def residue_class(atom: PdbAtom, labeled_residues: Set[Tuple[Optional[str], int]]) -> str:
    if not labeled_residues:
        return "heavy"
    if (None, atom.resseq) in labeled_residues or (atom.chain, atom.resseq) in labeled_residues:
        return "labeled"
    return "other"


def atom_label(atom: PdbAtom) -> str:
    altloc = f".{atom.altloc}" if atom.altloc else ""
    icode = atom.icode if atom.icode else ""
    return f"{atom.chain}:{atom.resname}{atom.resseq}{icode}:{atom.atom_name}{altloc}#{atom.serial}"


def model_spec(base_model: int, model_used: Optional[int]) -> str:
    """Build the leading model specifier shared by both programs.

    Chimera and ChimeraX open a multi-MODEL PDB as submodels, so MODEL 2 of the
    first opened structure is ``#0.2`` in Chimera and ``#1.2`` in ChimeraX.
    """
    if model_used is not None:
        return f"#{base_model}.{model_used}"
    return f"#{base_model}"


def chimera_atom_spec(atom: PdbAtom, base_model: int = 0, model_used: Optional[int] = None) -> str:
    """UCSF Chimera atom specifier: ``#model:residue.chain@atom``.

    The insertion code is appended to the residue number, and a blank chain ID
    is written as a bare period, so ``:12.`` means residue 12 in a no-ID chain.
    """
    chain = "" if atom.chain == "." else atom.chain
    return (
        f"{model_spec(base_model, model_used)}"
        f":{atom.resseq}{atom.icode}.{chain}@{atom.atom_name}"
    )


def chimerax_atom_spec(atom: PdbAtom, base_model: int = 1, model_used: Optional[int] = None) -> str:
    """UCSF ChimeraX atom specifier: ``#model/chain:residue@atom``.

    A blank chain ID has no ChimeraX specifier, so the chain part is omitted
    and the residue and atom still pin the selection down.
    """
    prefix = model_spec(base_model, model_used)
    chain = "" if atom.chain == "." else f"/{atom.chain}"
    return f"{prefix}{chain}:{atom.resseq}{atom.icode}@{atom.atom_name}"


def _unique_atoms(atoms: Iterable[PdbAtom]) -> List[PdbAtom]:
    """De-duplicate atoms by serial while preserving first-seen order."""
    seen: Set[int] = set()
    unique: List[PdbAtom] = []
    for atom in atoms:
        if atom.serial in seen:
            continue
        seen.add(atom.serial)
        unique.append(atom)
    return unique


def chimera_select_command(
    atoms: Iterable[PdbAtom], base_model: int = 0, model_used: Optional[int] = None
) -> str:
    """A Chimera ``select`` command covering every supplied atom."""
    specs = [chimera_atom_spec(atom, base_model, model_used) for atom in _unique_atoms(atoms)]
    return "select " + "|".join(specs) if specs else ""


def chimerax_select_command(
    atoms: Iterable[PdbAtom], base_model: int = 1, model_used: Optional[int] = None
) -> str:
    """A ChimeraX ``select`` command covering every supplied atom."""
    specs = [chimerax_atom_spec(atom, base_model, model_used) for atom in _unique_atoms(atoms)]
    return "select " + "|".join(specs) if specs else ""


def summarize_clashes(
    atoms: Sequence[PdbAtom],
    candidate_pairs: Iterable[Tuple[int, int]],
    cutoff: float,
    adjacent_window: int,
    circular: Dict[str, Tuple[int, int]],
    link_pairs: Set[frozenset],
    labeled_residues: Set[Tuple[Optional[str], int]],
) -> Tuple[List[Tuple[float, PdbAtom, PdbAtom]], "Counter[str]", Dict[str, int]]:
    """Classify every candidate pair into a clash or a named exclusion."""
    clashes: List[Tuple[float, PdbAtom, PdbAtom]] = []
    composition: Counter[str] = Counter()
    excluded = {"same_residue": 0, "altloc": 0, "adjacent": 0, "link": 0}

    for i, j in candidate_pairs:
        atom_i = atoms[i]
        atom_j = atoms[j]
        dist = distance(atom_i, atom_j)
        if dist >= cutoff:
            continue

        if atom_i.residue_key == atom_j.residue_key:
            excluded["same_residue"] += 1
            continue
        if is_alternate_conformation(atom_i, atom_j):
            excluded["altloc"] += 1
            continue
        if is_sequence_adjacent(atom_i, atom_j, adjacent_window, circular):
            excluded["adjacent"] += 1
            continue
        if is_linked_pair(atom_i, atom_j, link_pairs):
            excluded["link"] += 1
            continue

        clashes.append((dist, atom_i, atom_j))
        key = "-".join(
            sorted((residue_class(atom_i, labeled_residues), residue_class(atom_j, labeled_residues)))
        )
        composition[key] += 1

    clashes.sort(key=lambda item: (item[0], item[1].serial, item[2].serial))
    return clashes, composition, excluded


def check_clashes(
    pdb_path: Path,
    cutoff: float = DEFAULT_CUTOFF,
    adjacent_window: int = DEFAULT_ADJACENT_WINDOW,
    circular: str = "off",
    labeled_residues_text: str = "none",
    use_link_exclusion: bool = True,
    model: Optional[int] = None,
) -> ClashReport:
    """Run the full clash check and return a structured report."""
    pdb_path = Path(pdb_path)
    if cutoff <= 0:
        raise ValueError("--cutoff must be positive; got %g" % cutoff)
    if adjacent_window < 0:
        raise ValueError("--adjacent-window must be zero or positive; got %d" % adjacent_window)

    parsed = read_pdb(pdb_path, model=model)
    heavy_atoms = [atom for atom in parsed.atoms if not is_hydrogen(atom)]

    circular_mode, fixed_last = parse_circular_option(circular)
    circular_map = circular_bounds(heavy_atoms, circular_mode, fixed_last)
    labeled_residues = parse_residue_set(labeled_residues_text)
    link_pairs = parsed.link_pairs if use_link_exclusion else set()

    candidate_pairs, method = close_pairs_with_method(heavy_atoms, cutoff)
    clashes, composition, excluded = summarize_clashes(
        heavy_atoms,
        candidate_pairs,
        cutoff,
        adjacent_window,
        circular_map,
        link_pairs,
        labeled_residues,
    )

    return ClashReport(
        pdb_path=pdb_path,
        method=method,
        atoms_total=len(parsed.atoms),
        heavy_atoms=len(heavy_atoms),
        cutoff=cutoff,
        adjacent_window=adjacent_window,
        circular_mode=circular_mode if circular_mode != "fixed" else f"1-{fixed_last}",
        circular_chains=circular_map,
        link_records=len(parsed.link_pairs),
        link_excluded=excluded["link"],
        altloc_excluded=excluded["altloc"],
        adjacent_excluded=excluded["adjacent"] + excluded["same_residue"],
        clashes=clashes,
        composition=composition,
        model_numbers=parsed.model_numbers,
        model_used=parsed.model_used,
        skipped_atom_lines=parsed.skipped_atom_lines,
    )


def select_all_commands(
    report: ClashReport,
    chimera_model: int = DEFAULT_CHIMERA_MODEL,
    chimerax_model: int = DEFAULT_CHIMERAX_MODEL,
) -> Tuple[str, str]:
    """Return the (Chimera, ChimeraX) commands selecting every clashing atom."""
    every_atom = [atom for _dist, a, b in report.clashes for atom in (a, b)]
    return (
        chimera_select_command(every_atom, chimera_model, report.model_used),
        chimerax_select_command(every_atom, chimerax_model, report.model_used),
    )


def format_report(
    report: ClashReport,
    top: int = DEFAULT_TOP,
    include_select_commands: bool = True,
    chimera_model: int = DEFAULT_CHIMERA_MODEL,
    chimerax_model: int = DEFAULT_CHIMERAX_MODEL,
) -> str:
    """Render the report as the tab-separated key/value text block."""
    lines = [
        f"pdb\t{report.pdb_path}",
        f"method\t{report.method}",
    ]
    if report.model_numbers:
        lines.append(f"models_in_file\t{len(report.model_numbers)}")
        lines.append(f"model_used\t{report.model_used}")
    lines.extend(
        [
            f"atoms_total\t{report.atoms_total}",
            f"heavy_atoms\t{report.heavy_atoms}",
            f"cutoff_A\t{report.cutoff:.3f}",
            f"same_chain_adjacent_window\t{report.adjacent_window}",
            f"circular_closure\t{report.circular_mode}",
        ]
    )
    if report.circular_chains:
        closures = ", ".join(
            f"{chain}:{low}-{high}" for chain, (low, high) in sorted(report.circular_chains.items())
        )
        lines.append(f"circular_chains\t{closures}")
    lines.append(f"link_records_parsed\t{report.link_records}")
    lines.append(f"link_excluded_pairs\t{report.link_excluded}")
    lines.append(f"altloc_excluded_pairs\t{report.altloc_excluded}")
    lines.append(f"adjacent_excluded_pairs\t{report.adjacent_excluded}")
    if report.skipped_atom_lines:
        lines.append(f"skipped_unparsable_atom_lines\t{report.skipped_atom_lines}")
    lines.append(f"heavy_atom_clashes_lt_cutoff\t{len(report.clashes)}")
    if report.clashes:
        lines.append(f"min_clash_distance_A\t{report.clashes[0][0]:.6f}")
    else:
        lines.append("min_clash_distance_A\tNA")

    if report.composition:
        lines.append("composition")
        for key in sorted(report.composition):
            lines.append(f"{key}\t{report.composition[key]}")

    if report.clashes and include_select_commands:
        chimera_all, chimerax_all = select_all_commands(report, chimera_model, chimerax_model)
        # Each command sits alone on its line so a line selection, such as a
        # triple-click, picks up the whole command and nothing else.
        lines.append("chimera_select_all_clashes")
        lines.append(chimera_all)
        lines.append("chimerax_select_all_clashes")
        lines.append(chimerax_all)

    if top > 0 and report.clashes:
        header = "closest_clashes\tdistance_A\tatom_1\tatom_2"
        if include_select_commands:
            header += "\tchimera\tchimerax"
        lines.append(header)
        for dist, atom_i, atom_j in report.clashes[:top]:
            row = f"{dist:.6f}\t{atom_label(atom_i)}\t{atom_label(atom_j)}"
            if include_select_commands:
                pair = (atom_i, atom_j)
                row += "\t" + chimera_select_command(pair, chimera_model, report.model_used)
                row += "\t" + chimerax_select_command(pair, chimerax_model, report.model_used)
            lines.append(row)

    return "\n".join(lines)


def build_cli_command(
    script_name: str,
    pdb_path: str,
    cutoff: float,
    adjacent_window: int,
    circular: str,
    labeled_residues_text: str,
    use_link_exclusion: bool,
    top: int,
    model: Optional[int],
    output: str,
    include_select_commands: bool = True,
    chimera_model: int = DEFAULT_CHIMERA_MODEL,
    chimerax_model: int = DEFAULT_CHIMERAX_MODEL,
) -> str:
    parts: List[str] = [sys.executable, script_name, pdb_path]
    if cutoff != DEFAULT_CUTOFF:
        parts.extend(["--cutoff", f"{cutoff:g}"])
    if adjacent_window != DEFAULT_ADJACENT_WINDOW:
        parts.extend(["--adjacent-window", str(adjacent_window)])
    if (circular or "off").strip().lower() not in {"", "off"}:
        parts.extend(["--circular", circular.strip()])
    if (labeled_residues_text or "none").strip().lower() not in {"", "none"}:
        parts.extend(["--label-residues", labeled_residues_text.strip()])
    if not use_link_exclusion:
        parts.append("--no-link-exclusion")
    if top != DEFAULT_TOP:
        parts.extend(["--top", str(top)])
    if model is not None:
        parts.extend(["--model", str(model)])
    if not include_select_commands:
        parts.append("--no-select-commands")
    if chimera_model != DEFAULT_CHIMERA_MODEL:
        parts.extend(["--chimera-model", str(chimera_model)])
    if chimerax_model != DEFAULT_CHIMERAX_MODEL:
        parts.extend(["--chimerax-model", str(chimerax_model)])
    if output:
        parts.extend(["-o", output])
    return " ".join(shlex.quote(str(part)) for part in parts)


def default_output_path(input_path: Path) -> Path:
    return input_path.with_name(f"{input_path.stem}_clashes.txt")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="check_pdb_clashes.py",
        description="Count close heavy-atom contacts (clashes) in a PDB file.",
    )
    parser.add_argument("pdb", nargs="?", help="Input PDB file.")
    parser.add_argument(
        "--cutoff",
        type=float,
        default=DEFAULT_CUTOFF,
        help="Heavy-atom distance cutoff in Angstrom for a clash. Default: %(default)s.",
    )
    parser.add_argument(
        "--adjacent-window",
        type=int,
        default=DEFAULT_ADJACENT_WINDOW,
        help=(
            "Exclude same-chain residue pairs with sequence separation <= this value. "
            "Default: %(default)s."
        ),
    )
    parser.add_argument(
        "--circular",
        default="off",
        help=(
            "Close cyclic chains: 'off', 'auto' to use each chain's own first and last "
            "residue, or a residue number applied to every chain. Default: %(default)s."
        ),
    )
    parser.add_argument(
        "--label-residues",
        default="none",
        help=(
            "Residues to label for composition counts, as NUMBER or CHAIN:NUMBER, "
            "comma separated. Use 'none' to disable. Default: %(default)s."
        ),
    )
    parser.add_argument(
        "--no-link-exclusion",
        action="store_true",
        help="Do not exclude residue pairs joined by a LINK record.",
    )
    parser.add_argument(
        "--model",
        type=int,
        help="MODEL number to check in a multi-model file. Default: the first MODEL.",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=DEFAULT_TOP,
        help="Print this many closest clash pairs. Default: %(default)s.",
    )
    parser.add_argument(
        "--no-select-commands",
        action="store_true",
        help="Omit the UCSF Chimera and ChimeraX select commands from the report.",
    )
    parser.add_argument(
        "--chimera-model",
        type=int,
        default=DEFAULT_CHIMERA_MODEL,
        help=(
            "Chimera model number for the select commands. Default: %(default)s, the "
            "first opened structure."
        ),
    )
    parser.add_argument(
        "--chimerax-model",
        type=int,
        default=DEFAULT_CHIMERAX_MODEL,
        help=(
            "ChimeraX model number for the select commands. Default: %(default)s, the "
            "first opened structure."
        ),
    )
    parser.add_argument("-o", "--output", help="Write the report to this file as well as stdout.")
    parser.add_argument("--gui", action="store_true", help="Open the Tk GUI.")
    parser.add_argument("-v", "--version", action="version", version=f"{TOOL_NAME} V{VERSION}")
    return parser


HELP_TEXTS = {
    "input": (
        "PDB file to check. Any structure works: nucleic acid, protein, ligand, or a "
        "mixed assembly. Hydrogens and deuteriums are ignored, so only heavy atoms are "
        "compared."
    ),
    "cutoff": (
        "Two heavy atoms closer than this distance, in Angstrom, count as a clash.\n\n"
        "1.60 A is a tight geometric check that finds atoms driven essentially on top of "
        "each other by a build or an exchange. Raise it, for example to 2.2 A, to flag "
        "softer nonbonded contacts as well."
    ),
    "adjacent": (
        "Exclude same-chain residue pairs whose residue numbers differ by at most this "
        "value, so ordinary covalent neighbors are not reported.\n\n"
        "1 excludes i and i+1. Use 0 to compare every non-identical residue pair, "
        "including sequence neighbors."
    ),
    "circular": (
        "Close cyclic chains so the first-to-last junction is not reported as a clash.\n\n"
        "off   : no closure (default).\n"
        "auto  : each chain uses its own lowest and highest residue number, so chains of "
        "different lengths are handled correctly.\n"
        "N     : treat residue 1 and residue N as neighbors in every chain.\n\n"
        "LINK records already cover most cyclizations, so 'auto' is only needed when a "
        "cyclic model has no LINK record for its closing bond."
    ),
    "link": (
        "Exclude residue pairs joined by a LINK record, in addition to sequence "
        "neighbors.\n\n"
        "re_helix writes LINK records for cyclization, reciprocal exchange, and permuted "
        "chains. At those junctions residue numbers no longer track real connectivity, so "
        "a bonded pair can look far apart in sequence and would otherwise be counted as a "
        "clash. Turn this off to see the raw geometric contacts."
    ),
    "label": (
        "Optional residue labeling for the composition breakdown. Enter residue numbers "
        "such as 11,32,53, or chain-qualified entries such as A:11,B:32.\n\n"
        "Labeled residues are reported as 'labeled' and the rest as 'other', so the "
        "composition lines show whether clashes fall inside a region of interest such as "
        "a linker. Leave it as 'none' for a single 'heavy' class."
    ),
    "model": (
        "MODEL number to check in a multi-model file such as an NMR ensemble.\n\n"
        "Leave blank to use the first MODEL. Only one model is checked at a time so "
        "separate ensemble members are never compared against each other."
    ),
    "top": (
        "How many of the closest clash pairs to list, shortest distance first. Use 0 to "
        "report only the counts and the minimum distance."
    ),
    "output": (
        "Optional text file for the report. The same report is always printed to the run "
        "log. Default when left blank and 'Save report' is used: <input>_clashes.txt."
    ),
    "select": (
        "Add ready-to-paste selection commands for UCSF Chimera and UCSF ChimeraX.\n\n"
        "Two lines select every clashing atom at once:\n"
        "  chimera_select_all_clashes\n"
        "  chimerax_select_all_clashes\n\n"
        "Each listed clash pair also gets its own pair of commands, so a single contact "
        "can be isolated and inspected.\n\n"
        "Syntax used:\n"
        "  Chimera   select #0:20.A@C1'\n"
        "  ChimeraX  select #1/A:20@C1'\n\n"
        "Chimera numbers the first opened structure #0 and ChimeraX numbers it #1. "
        "Change the model numbers below if the structure is not the first one open. When "
        "a MODEL is selected from a multi-model file, the submodel is included "
        "automatically, for example #1.2 in ChimeraX.\n\n"
        "Copying a command:\n"
        "  Double-click any command in the report to highlight the whole command and "
        "copy it to the clipboard. A plain double-click would otherwise stop at the "
        "punctuation and grab only one word.\n"
        "  Copy Chimera select and Copy ChimeraX select, next to Run, copy the "
        "select-everything command without touching the report.\n"
        "  Each select-everything command also sits alone on its own line, so a "
        "triple-click or a line selection in the saved report file picks up the whole "
        "command and nothing else."
    ),
    "chimera_model": (
        "Model number used in the Chimera select commands. Chimera numbers the first "
        "opened structure #0, so leave this at 0 unless the structure is opened later in "
        "a session that already holds other models."
    ),
    "chimerax_model": (
        "Model number used in the ChimeraX select commands. ChimeraX numbers the first "
        "opened structure #1, so leave this at 1 unless the structure is opened later in "
        "a session that already holds other models."
    ),
}


def run_gui(initial_path: Optional[str] = None) -> int:
    try:
        import tkinter as tk
        from tkinter import filedialog, messagebox, scrolledtext
        from tkinter import ttk
    except Exception as exc:
        print("Error: GUI mode requires tkinter (%s)." % exc, file=sys.stderr)
        return 1

    root = tk.Tk()
    root.title(f"{TOOL_NAME} V{VERSION}")
    apply_optional_icon(root, __file__)
    root.geometry("940x760")

    style = ttk.Style(root)
    style.configure("ToolTitle.TLabel", font=("TkDefaultFont", 12, "bold"))
    style.configure("Tool.TLabelframe.Label", font=("TkDefaultFont", 10, "bold"))

    def make_help_button(parent, title: str, message: str):
        return tk.Button(
            parent,
            text="?",
            width=2,
            padx=0,
            pady=0,
            bg="#d9ecff",
            activebackground="#c4e0ff",
            highlightbackground="#d9ecff",
            relief="raised",
            bd=1,
            command=lambda: messagebox.showinfo(title, message, parent=root),
        )

    input_var = tk.StringVar(value=initial_path or "")
    cutoff_var = tk.StringVar(value=f"{DEFAULT_CUTOFF:g}")
    adjacent_var = tk.StringVar(value=str(DEFAULT_ADJACENT_WINDOW))
    circular_var = tk.StringVar(value="off")
    label_var = tk.StringVar(value="none")
    link_var = tk.BooleanVar(value=True)
    select_var = tk.BooleanVar(value=True)
    chimera_model_var = tk.StringVar(value=str(DEFAULT_CHIMERA_MODEL))
    chimerax_model_var = tk.StringVar(value=str(DEFAULT_CHIMERAX_MODEL))
    model_var = tk.StringVar(value="")
    top_var = tk.StringVar(value=str(DEFAULT_TOP))
    output_var = tk.StringVar()
    status_var = tk.StringVar(value="Choose a PDB file and press Run.")

    outer = ttk.Frame(root, padding=10)
    outer.pack(fill="both", expand=True)
    outer.columnconfigure(1, weight=1)

    ttk.Label(outer, text=f"{TOOL_NAME} V{VERSION}", style="ToolTitle.TLabel").grid(
        row=0, column=0, columnspan=4, sticky="w", pady=(0, 8)
    )

    ttk.Label(outer, text="Input PDB").grid(row=1, column=0, sticky="e", padx=(0, 6), pady=4)
    ttk.Entry(outer, textvariable=input_var, width=70).grid(row=1, column=1, sticky="ew", pady=4)

    def browse_input() -> None:
        path = filedialog.askopenfilename(
            title="Choose input PDB",
            filetypes=[("PDB files", "*.pdb"), ("All files", "*.*")],
        )
        if path:
            input_var.set(path)

    ttk.Button(outer, text="Browse...", command=browse_input).grid(
        row=1, column=2, sticky="ew", padx=(6, 0), pady=4
    )
    make_help_button(outer, "Input PDB", HELP_TEXTS["input"]).grid(row=1, column=3, padx=(6, 0))

    options_box = ttk.LabelFrame(outer, text="Clash criteria", padding=8, style="Tool.TLabelframe")
    options_box.grid(row=2, column=0, columnspan=4, sticky="ew", pady=(8, 4))
    options_box.columnconfigure(1, weight=1)
    options_box.columnconfigure(4, weight=1)

    def add_option_row(row: int, column: int, label: str, variable, help_key: str, width: int = 16):
        ttk.Label(options_box, text=label).grid(
            row=row, column=column, sticky="e", padx=(0, 6), pady=4
        )
        ttk.Entry(options_box, textvariable=variable, width=width).grid(
            row=row, column=column + 1, sticky="w", pady=4
        )
        make_help_button(options_box, label, HELP_TEXTS[help_key]).grid(
            row=row, column=column + 2, padx=(6, 12)
        )

    add_option_row(0, 0, "Cutoff (A)", cutoff_var, "cutoff")
    add_option_row(0, 3, "Adjacent window", adjacent_var, "adjacent")
    add_option_row(1, 0, "Circular closure", circular_var, "circular")
    add_option_row(1, 3, "Label residues", label_var, "label")
    add_option_row(2, 0, "MODEL number", model_var, "model")
    add_option_row(2, 3, "Top clashes", top_var, "top")

    link_row = ttk.Frame(options_box)
    link_row.grid(row=3, column=0, columnspan=6, sticky="w", pady=(6, 0))
    ttk.Checkbutton(
        link_row,
        text="Exclude residue pairs joined by a LINK record",
        variable=link_var,
    ).pack(side="left")
    make_help_button(link_row, "LINK exclusion", HELP_TEXTS["link"]).pack(side="left", padx=(6, 0))

    select_box = ttk.LabelFrame(
        outer, text="Chimera / ChimeraX select commands", padding=8, style="Tool.TLabelframe"
    )
    select_box.grid(row=3, column=0, columnspan=4, sticky="ew", pady=(4, 4))

    select_row = ttk.Frame(select_box)
    select_row.pack(fill="x")
    ttk.Checkbutton(
        select_row,
        text="Include select commands in the report",
        variable=select_var,
    ).pack(side="left")
    make_help_button(select_row, "Select commands", HELP_TEXTS["select"]).pack(
        side="left", padx=(6, 16)
    )
    ttk.Label(select_row, text="Chimera model #").pack(side="left")
    ttk.Entry(select_row, textvariable=chimera_model_var, width=5).pack(side="left", padx=(4, 0))
    make_help_button(select_row, "Chimera model number", HELP_TEXTS["chimera_model"]).pack(
        side="left", padx=(4, 16)
    )
    ttk.Label(select_row, text="ChimeraX model #").pack(side="left")
    ttk.Entry(select_row, textvariable=chimerax_model_var, width=5).pack(side="left", padx=(4, 0))
    make_help_button(select_row, "ChimeraX model number", HELP_TEXTS["chimerax_model"]).pack(
        side="left", padx=(4, 0)
    )

    ttk.Label(outer, text="Report file").grid(row=4, column=0, sticky="e", padx=(0, 6), pady=4)
    ttk.Entry(outer, textvariable=output_var, width=70).grid(row=4, column=1, sticky="ew", pady=4)

    def browse_output() -> None:
        path = filedialog.asksaveasfilename(
            title="Save clash report",
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")],
        )
        if path:
            output_var.set(path)

    ttk.Button(outer, text="Save as...", command=browse_output).grid(
        row=4, column=2, sticky="ew", padx=(6, 0), pady=4
    )
    make_help_button(outer, "Report file", HELP_TEXTS["output"]).grid(row=4, column=3, padx=(6, 0))

    results_box = ttk.LabelFrame(outer, text="Report", padding=8, style="Tool.TLabelframe")
    results_box.grid(row=5, column=0, columnspan=4, sticky="nsew", pady=(8, 4))
    results_box.columnconfigure(0, weight=1)
    results_box.rowconfigure(0, weight=1)
    outer.rowconfigure(5, weight=1)
    results_text = scrolledtext.ScrolledText(results_box, height=16, wrap="none")
    results_text.grid(row=0, column=0, sticky="nsew")
    ttk.Label(
        results_box,
        text="Double-click any select command to highlight the whole command and copy it.",
    ).grid(row=1, column=0, sticky="w", pady=(4, 0))

    # Select commands are tagged so a double-click can take the whole command
    # instead of the single word Tk would pick out of the punctuation.
    results_text.tag_configure("select_command", background="#e8f2ff")

    def tag_select_commands() -> None:
        results_text.tag_remove("select_command", "1.0", tk.END)
        last_line = int(results_text.index("end-1c").split(".")[0])
        for line_number in range(1, last_line + 1):
            line = results_text.get(f"{line_number}.0", f"{line_number}.end")
            search_from = 0
            while True:
                start = line.find("select ", search_from)
                if start == -1:
                    break
                end = line.find("\t", start)
                if end == -1:
                    end = len(line)
                results_text.tag_add(
                    "select_command", f"{line_number}.{start}", f"{line_number}.{end}"
                )
                search_from = end
                if search_from >= len(line):
                    break

    def copy_to_clipboard(text: str, description: str) -> None:
        root.clipboard_clear()
        root.clipboard_append(text)
        root.update_idletasks()
        shown = text if len(text) <= 70 else text[:67] + "..."
        status_var.set(f"Copied {description}: {shown}")

    def on_command_double_click(event) -> str:
        """Select the whole command under the pointer and copy it."""
        clicked = results_text.index(f"@{event.x},{event.y}")
        ranges = results_text.tag_ranges("select_command")
        for index in range(0, len(ranges), 2):
            start, end = ranges[index], ranges[index + 1]
            if results_text.compare(start, "<=", clicked) and results_text.compare(
                clicked, "<", end
            ):
                results_text.tag_remove("sel", "1.0", tk.END)
                results_text.tag_add("sel", start, end)
                results_text.mark_set("insert", start)
                results_text.focus_set()
                copy_to_clipboard(results_text.get(start, end), "command")
                return "break"
        return ""

    results_text.tag_bind("select_command", "<Double-Button-1>", on_command_double_click)

    ttk.Label(outer, textvariable=status_var, wraplength=840).grid(
        row=6, column=0, columnspan=4, sticky="ew", pady=4
    )

    buttons = ttk.Frame(outer)
    buttons.grid(row=7, column=0, columnspan=4, sticky="ew", pady=(8, 0))

    # Filled in by each run so the copy buttons always hold the latest commands.
    last_commands: Dict[str, str] = {"chimera": "", "chimerax": ""}
    copy_buttons: Dict[str, object] = {}

    def refresh_copy_buttons() -> None:
        for key, button in copy_buttons.items():
            button.configure(state="normal" if last_commands[key] else "disabled")

    def make_copy_command(key: str, description: str):
        def handler() -> None:
            command = last_commands[key]
            if command:
                copy_to_clipboard(command, description)

        return handler

    def run_from_gui() -> None:
        try:
            input_pdb = input_var.get().strip()
            if not input_pdb:
                raise ValueError("Please choose an input PDB file.")
            if not Path(input_pdb).is_file():
                raise ValueError("Input PDB file does not exist: %s" % input_pdb)

            try:
                cutoff = float(cutoff_var.get().strip())
            except ValueError:
                raise ValueError("Cutoff must be a number; got %r" % cutoff_var.get())
            try:
                adjacent_window = int(adjacent_var.get().strip())
            except ValueError:
                raise ValueError("Adjacent window must be an integer; got %r" % adjacent_var.get())
            try:
                top = int(top_var.get().strip())
            except ValueError:
                raise ValueError("Top clashes must be an integer; got %r" % top_var.get())

            model_text = model_var.get().strip()
            if model_text:
                try:
                    model = int(model_text)
                except ValueError:
                    raise ValueError("MODEL number must be an integer; got %r" % model_text)
            else:
                model = None

            try:
                chimera_model = int(chimera_model_var.get().strip())
            except ValueError:
                raise ValueError(
                    "Chimera model number must be an integer; got %r" % chimera_model_var.get()
                )
            try:
                chimerax_model = int(chimerax_model_var.get().strip())
            except ValueError:
                raise ValueError(
                    "ChimeraX model number must be an integer; got %r" % chimerax_model_var.get()
                )

            output_path = output_var.get().strip()

            cli_cmd = build_cli_command(
                Path(__file__).name,
                input_pdb,
                cutoff,
                adjacent_window,
                circular_var.get(),
                label_var.get(),
                link_var.get(),
                top,
                model,
                output_path,
                include_select_commands=select_var.get(),
                chimera_model=chimera_model,
                chimerax_model=chimerax_model,
            )
            print("Equivalent CLI command:", flush=True)
            print(cli_cmd, flush=True)

            report = check_clashes(
                Path(input_pdb),
                cutoff=cutoff,
                adjacent_window=adjacent_window,
                circular=circular_var.get(),
                labeled_residues_text=label_var.get(),
                use_link_exclusion=link_var.get(),
                model=model,
            )
            text = format_report(
                report,
                top=top,
                include_select_commands=select_var.get(),
                chimera_model=chimera_model,
                chimerax_model=chimerax_model,
            )
            print(text, flush=True)

            results_text.delete("1.0", tk.END)
            results_text.insert("1.0", text)
            tag_select_commands()

            if report.clashes and select_var.get():
                last_commands["chimera"], last_commands["chimerax"] = select_all_commands(
                    report, chimera_model, chimerax_model
                )
            else:
                last_commands["chimera"] = last_commands["chimerax"] = ""
            refresh_copy_buttons()

            if output_path:
                Path(output_path).write_text(text + "\n")
                print("Wrote: %s" % output_path, flush=True)

            summary = "%d clash%s below %.3f A" % (
                len(report.clashes),
                "" if len(report.clashes) == 1 else "es",
                report.cutoff,
            )
            if report.clashes:
                summary += " (closest %.3f A)" % report.clashes[0][0]
            if output_path:
                summary += ". Wrote %s" % output_path
            status_var.set(summary)
        except Exception as exc:
            print("Error: %s" % exc, flush=True)
            status_var.set("Error: %s" % exc)
            messagebox.showerror(TOOL_NAME, str(exc), parent=root)

    ttk.Button(buttons, text="Run", command=run_from_gui).pack(side="left")
    ttk.Button(buttons, text="Close", command=root.destroy).pack(side="left", padx=(6, 0))
    copy_buttons["chimera"] = ttk.Button(
        buttons,
        text="Copy Chimera select",
        command=make_copy_command("chimera", "Chimera command"),
        state="disabled",
    )
    copy_buttons["chimera"].pack(side="left", padx=(18, 0))
    copy_buttons["chimerax"] = ttk.Button(
        buttons,
        text="Copy ChimeraX select",
        command=make_copy_command("chimerax", "ChimeraX command"),
        state="disabled",
    )
    copy_buttons["chimerax"].pack(side="left", padx=(6, 0))

    root.mainloop()
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    if argv is None:
        argv = sys.argv[1:]
    if len(argv) == 0:
        return run_gui()

    parser = build_arg_parser()
    args = parser.parse_args(list(argv))
    if args.gui:
        return run_gui(initial_path=args.pdb)
    if not args.pdb:
        parser.print_help(sys.stderr)
        return 2

    try:
        cli_cmd = build_cli_command(
            Path(__file__).name,
            args.pdb,
            args.cutoff,
            args.adjacent_window,
            args.circular,
            args.label_residues,
            not args.no_link_exclusion,
            args.top,
            args.model,
            args.output or "",
            include_select_commands=not args.no_select_commands,
            chimera_model=args.chimera_model,
            chimerax_model=args.chimerax_model,
        )
        print("Equivalent CLI command:", flush=True)
        print(cli_cmd, flush=True)

        report = check_clashes(
            Path(args.pdb),
            cutoff=args.cutoff,
            adjacent_window=args.adjacent_window,
            circular=args.circular,
            labeled_residues_text=args.label_residues,
            use_link_exclusion=not args.no_link_exclusion,
            model=args.model,
        )
        text = format_report(
            report,
            top=args.top,
            include_select_commands=not args.no_select_commands,
            chimera_model=args.chimera_model,
            chimerax_model=args.chimerax_model,
        )
        print(text, flush=True)
        if args.output:
            Path(args.output).write_text(text + "\n")
            print("Wrote: %s" % args.output, flush=True)
        return 0
    except Exception as exc:
        print("Error: %s" % exc, file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
