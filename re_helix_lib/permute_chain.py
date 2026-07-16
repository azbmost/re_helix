#!/usr/bin/env python3
"""Cyclically permute and continuously renumber one or more PDB chains."""

from __future__ import annotations

import argparse
import re
import shlex
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

try:
    from re_helix_lib.gui_icon import apply_optional_icon
except ImportError:  # pragma: no cover - direct script execution fallback
    from gui_icon import apply_optional_icon


TOOL_NAME = "Permute Chain"
VERSION = "1.0"

COORD_RECORDS = {"ATOM", "HETATM", "ANISOU", "SIGATM", "SIGUIJ"}


@dataclass
class ResidueBlock:
    old_resseq: int
    old_icode: str
    resname: str
    lines: List[str] = field(default_factory=list)
    new_resseq: Optional[int] = None


@dataclass
class PermuteStats:
    chain_id: str
    requested_shift: int
    effective_shift: int
    residue_count: int
    numbering_start: int
    numbering_end: int
    old_first_resseq: int
    new_first_source_resseq: int
    coordinate_lines_changed: int = 0
    ter_lines_changed: int = 0
    het_lines_changed: int = 0
    link_endpoints_changed: int = 0
    remark_references_changed: int = 0
    remark_inventory_rebuilt: bool = False
    warnings: List[str] = field(default_factory=list)


def split_line_ending(line: str) -> Tuple[str, str]:
    if line.endswith("\r\n"):
        return line[:-2], "\r\n"
    if line.endswith("\n"):
        return line[:-1], "\n"
    return line, ""


def format_resseq(resseq: int) -> str:
    text = f"{resseq:4d}"
    if len(text) > 4:
        raise ValueError(
            f"Residue number {resseq} is outside the four-column PDB field."
        )
    return text


def replace_field(line: str, start: int, end: int, value: str) -> str:
    body, ending = split_line_ending(line)
    body = body.ljust(end)
    width = end - start
    if len(value) != width:
        raise ValueError(f"Replacement field must be exactly {width} characters: {value!r}")
    return body[:start] + value + body[end:] + ending


def parse_resseq_field(line: str, start: int, end: int) -> Optional[int]:
    body, _ending = split_line_ending(line)
    if len(body) < end:
        return None
    try:
        return int(body[start:end].strip())
    except ValueError:
        return None


def coordinate_identity(line: str) -> Optional[Tuple[str, int, str, str]]:
    if line[:6].strip() not in COORD_RECORDS:
        return None
    body, _ending = split_line_ending(line)
    if len(body) < 27:
        return None
    resseq = parse_resseq_field(line, 22, 26)
    if resseq is None:
        return None
    return body[21], resseq, body[26], body[17:20].strip()


def update_coordinate_residue(line: str, new_resseq: int) -> str:
    updated = replace_field(line, 22, 26, format_resseq(new_resseq))
    return replace_field(updated, 26, 27, " ")


def collect_residue_blocks(
    lines: Sequence[str], chain_id: str
) -> Tuple[List[ResidueBlock], int, int]:
    blocks: List[ResidueBlock] = []
    seen_keys = set()
    first_index: Optional[int] = None
    last_index: Optional[int] = None
    current: Optional[ResidueBlock] = None

    for index, line in enumerate(lines):
        identity = coordinate_identity(line)
        if identity is None or identity[0] != chain_id:
            continue
        _chain, resseq, icode, resname = identity
        key = (resseq, icode)
        if first_index is None:
            first_index = index
        last_index = index

        if current is None or (current.old_resseq, current.old_icode) != key:
            if key in seen_keys:
                raise ValueError(
                    f"Residue {chain_id}{resseq}{icode.strip()} appears in multiple "
                    "non-contiguous blocks."
                )
            seen_keys.add(key)
            current = ResidueBlock(resseq, icode, resname)
            blocks.append(current)
        current.lines.append(line)

    if not blocks or first_index is None or last_index is None:
        raise ValueError(f"No coordinate records found for chain '{chain_id}'.")

    old_numbers = [block.old_resseq for block in blocks]
    if len(set(old_numbers)) != len(old_numbers):
        raise ValueError(
            f"Chain '{chain_id}' uses insertion codes with duplicate residue numbers; "
            "continuous permutation would make REMARK/LINK references ambiguous."
        )

    for line in lines[first_index : last_index + 1]:
        identity = coordinate_identity(line)
        if identity is not None and identity[0] != chain_id:
            raise ValueError(
                f"Chain '{chain_id}' is interleaved with chain '{identity[0]}'; "
                "the chain must occupy one coordinate section."
            )

    return blocks, first_index, last_index


def build_permutation(
    blocks: Sequence[ResidueBlock], shift: int
) -> Tuple[List[ResidueBlock], Dict[int, int], int, int, int]:
    count = len(blocks)
    if count == 0:
        raise ValueError("Cannot permute an empty chain.")
    effective_shift = shift % count
    rotated = list(blocks[effective_shift:]) + list(blocks[:effective_shift])
    numbering_start = min(block.old_resseq for block in blocks)
    numbering_end = numbering_start + count - 1
    format_resseq(numbering_start)
    format_resseq(numbering_end)

    mapping: Dict[int, int] = {}
    for offset, block in enumerate(rotated):
        new_resseq = numbering_start + offset
        block.new_resseq = new_resseq
        mapping[block.old_resseq] = new_resseq
    return rotated, mapping, effective_shift, numbering_start, numbering_end


def update_residue_endpoint(
    line: str,
    chain_index: int,
    resseq_start: int,
    resseq_end: int,
    icode_index: Optional[int],
    chain_id: str,
    mapping: Dict[int, int],
) -> Tuple[str, bool, bool]:
    body, _ending = split_line_ending(line)
    if len(body) <= chain_index or body[chain_index] != chain_id:
        return line, False, False
    old_resseq = parse_resseq_field(line, resseq_start, resseq_end)
    if old_resseq is None or old_resseq not in mapping:
        return line, False, False
    new_resseq = mapping[old_resseq]
    updated = replace_field(line, resseq_start, resseq_end, format_resseq(new_resseq))
    if icode_index is not None:
        updated = replace_field(updated, icode_index, icode_index + 1, " ")
    return updated, True, new_resseq != old_resseq


def update_link_line(
    line: str, chain_id: str, mapping: Dict[int, int]
) -> Tuple[str, int]:
    updated = line
    changed_count = 0
    for chain_index, start, end, icode_index in (
        (21, 22, 26, 26),
        (51, 52, 56, 56),
    ):
        updated, parsed, changed = update_residue_endpoint(
            updated,
            chain_index,
            start,
            end,
            icode_index,
            chain_id,
            mapping,
        )
        if parsed and changed:
            changed_count += 1
    return updated, changed_count


def update_het_line(
    line: str, chain_id: str, mapping: Dict[int, int]
) -> Tuple[str, bool]:
    updated, parsed, changed = update_residue_endpoint(
        line, 12, 13, 17, 17, chain_id, mapping
    )
    return updated, parsed and changed


def update_ter_line(
    line: str,
    chain_id: str,
    last_block: ResidueBlock,
    numbering_end: int,
) -> Tuple[str, bool]:
    body, _ending = split_line_ending(line)
    if len(body) <= 21 or body[21] != chain_id:
        return line, False
    updated = replace_field(line, 17, 20, f"{last_block.resname:>3s}"[-3:])
    updated = replace_field(updated, 22, 26, format_resseq(numbering_end))
    updated = replace_field(updated, 26, 27, " ")
    return updated, updated != line


def replace_remark_residue_references(
    line: str, chain_id: str, mapping: Dict[int, int]
) -> Tuple[str, int]:
    if not line.startswith("REMARK"):
        return line, 0

    body, ending = split_line_ending(line)
    changed = 0
    escaped_chain = re.escape(chain_id)

    colon_pattern = re.compile(
        rf"(?<![A-Za-z0-9]){escaped_chain}:(-?\d+)(?=[:\s,]|$)"
    )

    def replace_colon(match: re.Match[str]) -> str:
        nonlocal changed
        old_resseq = int(match.group(1))
        if old_resseq not in mapping:
            return match.group(0)
        new_resseq = mapping[old_resseq]
        if new_resseq != old_resseq:
            changed += 1
        return f"{chain_id}:{new_resseq}"

    body = colon_pattern.sub(replace_colon, body)

    token_prefix = r"(^|[\s'\"=,])"
    token_suffix = r"(?=$|[\s'\",])"
    suffix_pattern = re.compile(
        token_prefix + rf"(-?\d+)\.?{escaped_chain}" + token_suffix
    )
    prefix_pattern = re.compile(
        token_prefix + rf"{escaped_chain}\.?(-?\d+)" + token_suffix
    )

    def replace_suffix(match: re.Match[str]) -> str:
        nonlocal changed
        old_resseq = int(match.group(2))
        if old_resseq not in mapping:
            return match.group(0)
        new_resseq = mapping[old_resseq]
        if new_resseq != old_resseq:
            changed += 1
        return f"{match.group(1)}{new_resseq}{chain_id}"

    def replace_prefix(match: re.Match[str]) -> str:
        nonlocal changed
        old_resseq = int(match.group(2))
        if old_resseq not in mapping:
            return match.group(0)
        new_resseq = mapping[old_resseq]
        if new_resseq != old_resseq:
            changed += 1
        return f"{match.group(1)}{chain_id}{new_resseq}"

    body = suffix_pattern.sub(replace_suffix, body)
    body = prefix_pattern.sub(replace_prefix, body)
    return body + ending, changed


def is_chain_inventory_remark(line: str, chain_id: str) -> bool:
    pattern = re.compile(
        rf"^REMARK 950 RE_SCRIPT CHAIN_(?:RANGE|RESIDUES) chain={re.escape(chain_id)}(?:\s|$)"
    )
    return bool(pattern.match(line))


def build_chain_inventory_remarks(
    chain_id: str, rotated: Sequence[ResidueBlock], line_ending: str
) -> List[str]:
    residues = [
        (block.new_resseq, block.resname)
        for block in rotated
        if block.new_resseq is not None
    ]
    if not residues:
        return []
    start_resseq, start_name = residues[0]
    end_resseq, end_name = residues[-1]
    result = [
        "REMARK 950 RE_SCRIPT CHAIN_RANGE "
        f"chain={chain_id} start={chain_id}:{start_resseq}:{start_name} "
        f"end={chain_id}:{end_resseq}:{end_name} count={len(residues)}{line_ending}"
    ]
    chunks = [residues[index : index + 24] for index in range(0, len(residues), 24)]
    part_count = len(chunks)
    for part_index, chunk in enumerate(chunks, start=1):
        residue_text = ",".join(
            f"{chain_id}:{resseq}:{resname}" for resseq, resname in chunk
        )
        result.append(
            "REMARK 950 RE_SCRIPT CHAIN_RESIDUES "
            f"chain={chain_id} part={part_index}/{part_count} "
            f"residues={residue_text}{line_ending}"
        )
    return result


def permute_pdb_lines(
    lines: Sequence[str], chain_id: str, shift: int
) -> Tuple[List[str], PermuteStats]:
    chain_id = str(chain_id).strip()
    if len(chain_id) != 1:
        raise ValueError("Chain ID must be exactly one character.")

    blocks, first_coord_index, _last_coord_index = collect_residue_blocks(lines, chain_id)
    rotated, mapping, effective_shift, numbering_start, numbering_end = build_permutation(
        blocks, shift
    )
    stats = PermuteStats(
        chain_id=chain_id,
        requested_shift=shift,
        effective_shift=effective_shift,
        residue_count=len(blocks),
        numbering_start=numbering_start,
        numbering_end=numbering_end,
        old_first_resseq=blocks[0].old_resseq,
        new_first_source_resseq=rotated[0].old_resseq,
    )

    inventory_indices = [
        index
        for index, line in enumerate(lines)
        if is_chain_inventory_remark(line, chain_id)
    ]
    inventory_first = inventory_indices[0] if inventory_indices else None
    inventory_set = set(inventory_indices)
    inventory_lines: List[str] = []
    if inventory_first is not None:
        _body, inventory_ending = split_line_ending(lines[inventory_first])
        inventory_lines = build_chain_inventory_remarks(
            chain_id, rotated, inventory_ending or "\n"
        )
        stats.remark_inventory_rebuilt = True
    else:
        stats.warnings.append(
            f"No RE_SCRIPT CHAIN_RANGE/CHAIN_RESIDUES inventory found for chain {chain_id}."
        )

    transformed_blocks: List[str] = []
    for block in rotated:
        assert block.new_resseq is not None
        for line in block.lines:
            transformed = update_coordinate_residue(line, block.new_resseq)
            if transformed != line:
                stats.coordinate_lines_changed += 1
            transformed_blocks.append(transformed)

    output: List[str] = []
    emitted_coordinates = False
    emitted_inventory = False
    last_block = rotated[-1]

    for index, line in enumerate(lines):
        identity = coordinate_identity(line)
        if identity is not None and identity[0] == chain_id:
            if not emitted_coordinates and index == first_coord_index:
                output.extend(transformed_blocks)
                emitted_coordinates = True
            continue

        if index in inventory_set:
            if not emitted_inventory and index == inventory_first:
                output.extend(inventory_lines)
                emitted_inventory = True
            continue

        record = line[:6].strip()
        updated = line
        if record == "LINK":
            updated, changed_count = update_link_line(updated, chain_id, mapping)
            stats.link_endpoints_changed += changed_count
        elif record == "HET":
            updated, changed = update_het_line(updated, chain_id, mapping)
            if changed:
                stats.het_lines_changed += 1
        elif record == "TER":
            updated, changed = update_ter_line(
                updated, chain_id, last_block, numbering_end
            )
            if changed:
                stats.ter_lines_changed += 1
        elif record == "REMARK":
            updated, changed_count = replace_remark_residue_references(
                updated, chain_id, mapping
            )
            stats.remark_references_changed += changed_count
        output.append(updated)

    return output, stats


def default_output_path(input_path: Path) -> Path:
    return input_path.with_name(f"{input_path.stem}_permuted{input_path.suffix}")


def permute_chain(
    input_pdb: Path,
    chain_id: str,
    shift: int,
    output_pdb: Optional[Path] = None,
    verbose: bool = True,
) -> Tuple[Path, PermuteStats]:
    output_path, all_stats = permute_chains(
        input_pdb,
        [(chain_id, int(shift))],
        output_pdb=output_pdb,
        verbose=verbose,
    )
    return output_path, all_stats[0]


def validate_permutation_specs(
    specs: Sequence[Tuple[str, int]],
) -> List[Tuple[str, int]]:
    if not specs:
        raise ValueError("At least one permutation site is required.")
    normalized: List[Tuple[str, int]] = []
    seen_chains = set()
    for site_index, (chain_id, shift) in enumerate(specs, start=1):
        chain_id = str(chain_id).strip()
        if len(chain_id) != 1:
            raise ValueError(
                f"Permutation site {site_index}: chain ID must be exactly one character."
            )
        if chain_id in seen_chains:
            raise ValueError(
                f"Chain '{chain_id}' occurs more than once; combine its shifts into one site."
            )
        seen_chains.add(chain_id)
        normalized.append((chain_id, int(shift)))
    return normalized


def permute_chains(
    input_pdb: Path,
    specs: Sequence[Tuple[str, int]],
    output_pdb: Optional[Path] = None,
    verbose: bool = True,
) -> Tuple[Path, List[PermuteStats]]:
    input_pdb = Path(input_pdb)
    if not input_pdb.is_file():
        raise ValueError(f"Input PDB does not exist: {input_pdb}")
    normalized_specs = validate_permutation_specs(specs)
    output_path = Path(output_pdb) if output_pdb is not None else default_output_path(input_pdb)
    output_lines = input_pdb.read_text(errors="replace").splitlines(True)
    all_stats: List[PermuteStats] = []
    for chain_id, shift in normalized_specs:
        output_lines, stats = permute_pdb_lines(output_lines, chain_id, shift)
        all_stats.append(stats)
    output_path.write_text("".join(output_lines))
    if verbose:
        print(format_batch_summary(output_path, all_stats), flush=True)
    return output_path, all_stats


def build_cli_command(
    script_name: str,
    input_pdb: str,
    specs: Sequence[Tuple[str, int]],
    output_pdb: str,
) -> str:
    parts = [sys.executable, script_name, input_pdb]
    for chain_id, shift in specs:
        parts.extend(["--permute", chain_id, str(shift)])
    if output_pdb:
        parts.extend(["-o", output_pdb])
    return " ".join(shlex.quote(str(part)) for part in parts)


def format_summary(output_pdb: Path, stats: PermuteStats) -> str:
    signed_shift = f"{stats.requested_shift:+d}"
    lines = [
        f"Chain: {stats.chain_id}",
        f"Requested shift: {signed_shift}",
        f"Effective forward shift: {stats.effective_shift}",
        f"Residues: {stats.residue_count}",
        f"New first residue source: old {stats.chain_id}{stats.new_first_source_resseq}",
        f"Continuous numbering: {stats.numbering_start}-{stats.numbering_end}",
        f"Coordinate-like records changed: {stats.coordinate_lines_changed}",
        f"TER records changed: {stats.ter_lines_changed}",
        f"HET records changed: {stats.het_lines_changed}",
        f"LINK endpoints changed: {stats.link_endpoints_changed}",
        f"REMARK residue references changed: {stats.remark_references_changed}",
        f"REMARK chain inventory rebuilt: {'yes' if stats.remark_inventory_rebuilt else 'no'}",
        f"Wrote: {output_pdb}",
    ]
    if stats.warnings:
        lines.append("Warnings:")
        lines.extend(f"  {warning}" for warning in stats.warnings)
    return "\n".join(lines)


def format_batch_summary(
    output_pdb: Path, all_stats: Sequence[PermuteStats]
) -> str:
    if len(all_stats) == 1:
        return format_summary(output_pdb, all_stats[0])
    lines = [f"Permutation sites: {len(all_stats)}"]
    for site_index, stats in enumerate(all_stats, start=1):
        lines.extend(
            [
                "",
                f"Site {site_index}: chain {stats.chain_id}, shift {stats.requested_shift:+d}",
                f"  Effective forward shift: {stats.effective_shift}",
                f"  Residues: {stats.residue_count}",
                f"  New first residue source: old {stats.chain_id}{stats.new_first_source_resseq}",
                f"  Continuous numbering: {stats.numbering_start}-{stats.numbering_end}",
                f"  Coordinate-like records changed: {stats.coordinate_lines_changed}",
                f"  TER/HET/LINK changes: {stats.ter_lines_changed}/{stats.het_lines_changed}/{stats.link_endpoints_changed}",
                f"  REMARK references changed: {stats.remark_references_changed}",
            ]
        )
        lines.extend(f"  Warning: {warning}" for warning in stats.warnings)
    lines.extend(["", f"Wrote: {output_pdb}"])
    return "\n".join(lines)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="permute_chain.py",
        description=(
            "Cyclically rotate one or more PDB chains by residue position, "
            "renumber them continuously, and update REMARK/LINK residue references."
        ),
    )
    parser.add_argument("pdb", nargs="?", help="Input PDB file.")
    parser.add_argument("chain", nargs="?", help="One-character chain ID.")
    parser.add_argument(
        "shift",
        nargs="?",
        type=int,
        help="Signed residue shift: +5 moves the first five to the end; -5 moves the last five to the front.",
    )
    parser.add_argument(
        "--permute",
        nargs=2,
        action="append",
        metavar=("CHAIN", "SHIFT"),
        help=(
            "Permutation site. Repeat for multiple chains, for example "
            "--permute A 5 --permute B -5."
        ),
    )
    parser.add_argument("-o", "--output", help="Output PDB. Default: <input>_permuted.pdb")
    parser.add_argument("--gui", action="store_true", help="Open the Tk GUI.")
    parser.add_argument("-v", "--version", action="version", version=f"{TOOL_NAME} V{VERSION}")
    return parser


def run_gui(initial_path: Optional[str] = None) -> int:
    try:
        import tkinter as tk
        from tkinter import filedialog, messagebox, scrolledtext, ttk
    except Exception as exc:
        print(f"Error: GUI mode requires tkinter ({exc}).", file=sys.stderr)
        return 1

    root = tk.Tk()
    root.title(f"{TOOL_NAME} V{VERSION}")
    apply_optional_icon(root, __file__)
    root.geometry("780x500")

    style = ttk.Style(root)
    style.configure("ToolTitle.TLabel", font=("TkDefaultFont", 12, "bold"))
    style.configure("Tool.TLabelframe.Label", font=("TkDefaultFont", 10, "bold"))

    signed_shift_help = (
        "+n: POSITIVE SHIFT\n\n"
        "A positive shift +n moves the first n complete residue blocks of the "
        "selected chain to its end. The original residue at position n+1 becomes "
        "the first residue in the new chain order.\n\n"
        "Example: suppose chain A contains eight residues numbered A10 through A17. "
        "With +3, their source order becomes:\n\n"
        "A13, A14, A15, A16, A17, A10, A11, A12\n\n"
        "The output is then relabeled continuously as A10 through A17. Thus original "
        "A13 becomes new A10, and original A10 becomes new A15. Every ATOM/HETATM "
        "record belonging to a residue moves together.\n\n"
        "Shifts wrap around the chain length. For an eight-residue chain, +10 has "
        "the same ordering effect as +2; +8 and +0 leave the order unchanged. "
        "Related REMARK, HET, LINK, and TER residue references are updated using the "
        "same old-to-new numbering map.\n\n"
        "-n: NEGATIVE SHIFT\n\n"
        "A negative shift -n moves the last n complete residue blocks of the selected "
        "chain to its beginning, preserving the order within that moved block.\n\n"
        "Example: suppose chain A contains eight residues numbered A10 through A17. "
        "With -3, their source order becomes:\n\n"
        "A15, A16, A17, A10, A11, A12, A13, A14\n\n"
        "The output is then relabeled continuously as A10 through A17. Thus original "
        "A15 becomes new A10, and original A10 becomes new A13. Every ATOM/HETATM "
        "record belonging to a residue moves together.\n\n"
        "Shifts wrap around the chain length. For an eight-residue chain, -10 has "
        "the same ordering effect as -2; -8 and 0 leave the order unchanged. Related "
        "REMARK, HET, LINK, and TER residue references are updated using the same "
        "old-to-new numbering map."
    )

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
    output_var = tk.StringVar()
    site_count_var = tk.StringVar(value="1")
    site_vars: List[Tuple[tk.StringVar, tk.StringVar]] = []
    status_var = tk.StringVar(value="Ready")

    if input_var.get().strip():
        output_var.set(str(default_output_path(Path(input_var.get().strip()))))

    outer = ttk.Frame(root, padding=10)
    outer.pack(fill="both", expand=True)
    outer.columnconfigure(1, weight=1)
    outer.rowconfigure(5, weight=1)

    ttk.Label(outer, text=f"{TOOL_NAME} V{VERSION}", style="ToolTitle.TLabel").grid(
        row=0, column=0, columnspan=3, sticky="w", pady=(0, 8)
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
            output_var.set(str(default_output_path(Path(path))))

    ttk.Button(outer, text="Browse...", command=browse_input).grid(
        row=1, column=2, sticky="ew", padx=(6, 0), pady=4
    )

    ttk.Label(outer, text="Output PDB").grid(row=2, column=0, sticky="e", padx=(0, 6), pady=4)
    ttk.Entry(outer, textvariable=output_var, width=70).grid(row=2, column=1, sticky="ew", pady=4)

    def browse_output() -> None:
        path = filedialog.asksaveasfilename(
            title="Choose output PDB",
            defaultextension=".pdb",
            filetypes=[("PDB files", "*.pdb"), ("All files", "*.*")],
        )
        if path:
            output_var.set(path)

    ttk.Button(outer, text="Save as...", command=browse_output).grid(
        row=2, column=2, sticky="ew", padx=(6, 0), pady=4
    )

    settings = ttk.LabelFrame(outer, text="Permutation", padding=8, style="Tool.TLabelframe")
    settings.grid(row=3, column=0, columnspan=3, sticky="ew", pady=(8, 4))
    settings.columnconfigure(0, weight=1)

    site_header = ttk.Frame(settings)
    site_header.grid(row=0, column=0, sticky="ew", pady=(0, 6))
    ttk.Label(site_header, text="Number of permutation sites").pack(side="left")
    site_count_spinbox = ttk.Spinbox(
        site_header,
        from_=1,
        to=50,
        width=6,
        textvariable=site_count_var,
    )
    site_count_spinbox.pack(side="left", padx=(6, 12))
    ttk.Label(site_header, text="Signed shift (+n / -n)").pack(side="left")
    make_help_button(site_header, "Signed shift (+n / -n)", signed_shift_help).pack(
        side="left", padx=(4, 0)
    )

    site_area = ttk.Frame(settings)
    site_area.grid(row=1, column=0, sticky="ew")
    site_area.columnconfigure(0, weight=1)
    site_canvas = tk.Canvas(
        site_area,
        height=132,
        highlightthickness=0,
        background=style.lookup("TFrame", "background"),
    )
    site_scrollbar = ttk.Scrollbar(site_area, orient="vertical", command=site_canvas.yview)
    site_canvas.configure(yscrollcommand=site_scrollbar.set)
    site_canvas.grid(row=0, column=0, sticky="ew")
    site_scrollbar.grid(row=0, column=1, sticky="ns")
    site_rows_frame = ttk.Frame(site_canvas)
    site_window = site_canvas.create_window((0, 0), window=site_rows_frame, anchor="nw")

    def refresh_site_scrollregion(_event=None) -> None:
        site_canvas.configure(scrollregion=site_canvas.bbox("all"))

    def resize_site_rows(event) -> None:
        site_canvas.itemconfigure(site_window, width=event.width)

    site_canvas.bind("<Configure>", resize_site_rows)
    site_rows_frame.bind("<Configure>", refresh_site_scrollregion)

    def render_site_rows() -> None:
        nonlocal site_vars
        try:
            target_count = int(site_count_var.get())
        except ValueError:
            return
        if target_count < 1 or target_count > 50:
            return
        previous_values = [
            (chain_var.get(), shift_var.get())
            for chain_var, shift_var in site_vars
        ]
        for child in site_rows_frame.winfo_children():
            child.destroy()
        ttk.Label(site_rows_frame, text="#", width=4).grid(row=0, column=0, sticky="w")
        ttk.Label(site_rows_frame, text="Chain ID", width=12).grid(
            row=0, column=1, sticky="w"
        )
        ttk.Label(site_rows_frame, text="Signed shift").grid(
            row=0, column=2, sticky="w"
        )
        site_vars = []
        for site_index in range(target_count):
            old_chain, old_shift = (
                previous_values[site_index]
                if site_index < len(previous_values)
                else ("", "0")
            )
            chain_var = tk.StringVar(value=old_chain)
            shift_var = tk.StringVar(value=old_shift)
            ttk.Label(site_rows_frame, text=str(site_index + 1), width=4).grid(
                row=site_index + 1, column=0, sticky="w", pady=2
            )
            ttk.Entry(site_rows_frame, textvariable=chain_var, width=10).grid(
                row=site_index + 1, column=1, sticky="w", pady=2
            )
            ttk.Entry(site_rows_frame, textvariable=shift_var, width=14).grid(
                row=site_index + 1, column=2, sticky="w", pady=2
            )
            site_vars.append((chain_var, shift_var))
        root.after_idle(refresh_site_scrollregion)

    site_count_var.trace_add("write", lambda *_args: root.after_idle(render_site_rows))
    site_count_spinbox.configure(command=render_site_rows)
    render_site_rows()

    ttk.Label(outer, textvariable=status_var, wraplength=740).grid(
        row=4, column=0, columnspan=3, sticky="ew", pady=4
    )
    log_widget = scrolledtext.ScrolledText(outer, wrap="word", height=14)
    log_widget.grid(row=5, column=0, columnspan=3, sticky="nsew", pady=(4, 8))

    buttons = ttk.Frame(outer)
    buttons.grid(row=6, column=0, columnspan=3, sticky="ew")

    def run_from_gui() -> None:
        try:
            input_text = input_var.get().strip()
            if not input_text:
                raise ValueError("Please choose an input PDB file.")
            specs: List[Tuple[str, int]] = []
            for site_index, (chain_var, shift_var) in enumerate(site_vars, start=1):
                chain_id = chain_var.get().strip()
                if len(chain_id) != 1:
                    raise ValueError(
                        f"Permutation site {site_index}: chain ID must be exactly one character."
                    )
                try:
                    shift = int(shift_var.get().strip())
                except ValueError as exc:
                    raise ValueError(
                        f"Permutation site {site_index}: signed shift must be an integer."
                    ) from exc
                specs.append((chain_id, shift))
            specs = validate_permutation_specs(specs)
            output_text = output_var.get().strip() or str(default_output_path(Path(input_text)))
            command = build_cli_command(
                Path(__file__).name, input_text, specs, output_text
            )
            log_widget.delete("1.0", tk.END)
            log_widget.insert(tk.END, "Equivalent CLI command:\n" + command + "\n\n")
            output_path, all_stats = permute_chains(
                Path(input_text),
                specs,
                output_pdb=Path(output_text),
                verbose=False,
            )
            summary = format_batch_summary(output_path, all_stats)
            log_widget.insert(tk.END, summary + "\n")
            log_widget.see(tk.END)
            status_var.set(f"Wrote {output_path}")
            print("Equivalent CLI command:", flush=True)
            print(command, flush=True)
            print(summary, flush=True)
            messagebox.showinfo(TOOL_NAME, summary, parent=root)
        except Exception as exc:
            status_var.set(f"Error: {exc}")
            log_widget.insert(tk.END, f"Error: {exc}\n")
            print(f"Error: {exc}", flush=True)
            messagebox.showerror(TOOL_NAME, str(exc), parent=root)

    ttk.Button(buttons, text="Run", command=run_from_gui).pack(side="left")
    ttk.Button(buttons, text="Close", command=root.destroy).pack(side="left", padx=(6, 0))

    root.mainloop()
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    if argv is None:
        argv = sys.argv[1:]
    if not argv:
        return run_gui()

    parser = build_arg_parser()
    args = parser.parse_args(list(argv))
    if args.gui:
        return run_gui(initial_path=args.pdb)
    if not args.pdb:
        parser.print_help(sys.stderr)
        return 2
    try:
        if args.permute:
            if args.chain is not None or args.shift is not None:
                raise ValueError(
                    "Use either legacy positional CHAIN SHIFT arguments or repeated "
                    "--permute CHAIN SHIFT arguments, not both."
                )
            specs: List[Tuple[str, int]] = []
            for site_index, (chain_id, shift_text) in enumerate(args.permute, start=1):
                try:
                    shift = int(shift_text)
                except ValueError as exc:
                    raise ValueError(
                        f"Permutation site {site_index}: signed shift must be an integer."
                    ) from exc
                specs.append((chain_id, shift))
        else:
            if args.chain is None or args.shift is None:
                raise ValueError(
                    "Provide CHAIN SHIFT or at least one --permute CHAIN SHIFT site."
                )
            specs = [(args.chain, args.shift)]
        specs = validate_permutation_specs(specs)
        command = build_cli_command(
            Path(__file__).name,
            args.pdb,
            specs,
            args.output or "",
        )
        print("Equivalent CLI command:", flush=True)
        print(command, flush=True)
        permute_chains(
            Path(args.pdb),
            specs,
            output_pdb=Path(args.output) if args.output else None,
            verbose=True,
        )
        return 0
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
