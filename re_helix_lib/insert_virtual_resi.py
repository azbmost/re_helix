#!/usr/bin/env python3
"""insert_virtual_resi.py

Insert virtual residue-numbering gaps after selected residues in a PDB file.

This tool does not add ATOM records. Instead, each insertion spec shifts the
residue numbers of records after the designated residue in the same chain. For
example, inserting 3 virtual residues after A55 changes original A56 to A59,
A57 to A60, and so on. LINK record residue numbers are shifted with the same
mapping so topology records continue to point at the correct residues.
"""

from __future__ import annotations

import argparse
import re
import shlex
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

TOOL_NAME = "Insert Virtual Resi"
VERSION = "1.0"

COORD_RECORDS = {"ATOM", "HETATM", "ANISOU", "SIGATM", "SIGUIJ"}


@dataclass(frozen=True)
class InsertionSpec:
    chain_id: str
    after_resseq: int
    count: int
    token: str

    @property
    def label(self) -> str:
        return f"{self.chain_id}{self.after_resseq}+{self.count}"


@dataclass
class RenumberStats:
    specs: List[InsertionSpec]
    coordinate_lines_changed: int = 0
    ter_lines_changed: int = 0
    link_endpoints_changed: int = 0
    warnings: List[str] = field(default_factory=list)


def parse_residue_token(token: str) -> Tuple[str, int]:
    """Parse A55, A.55, 55A, or 55.A into (chain_id, resSeq)."""
    text = token.strip()
    patterns = (
        r"^([A-Za-z0-9])\.(-?\d+)$",
        r"^(-?\d+)\.([A-Za-z0-9])$",
        r"^([A-Za-z])(-?\d+)$",
        r"^(-?\d+)([A-Za-z])$",
    )
    for pattern in patterns:
        match = re.fullmatch(pattern, text)
        if not match:
            continue
        first, second = match.groups()
        if first.lstrip("-").isdigit():
            return second, int(first)
        return first, int(second)
    raise ValueError(
        "Residue token must look like A55, A.55, 55A, or 55.A; got %r" % token
    )


def parse_insertion_spec(token: str, count_text: str) -> InsertionSpec:
    chain_id, after_resseq = parse_residue_token(token)
    try:
        count = int(str(count_text).strip())
    except ValueError:
        raise ValueError("Insertion count must be an integer for %s: %r" % (token, count_text))
    if count <= 0:
        raise ValueError("Insertion count must be positive for %s: %d" % (token, count))
    return InsertionSpec(chain_id=chain_id, after_resseq=after_resseq, count=count, token=token)


def parse_insert_pairs(pairs: Optional[Sequence[Sequence[str]]]) -> List[InsertionSpec]:
    specs: List[InsertionSpec] = []
    for pair in pairs or []:
        if len(pair) != 2:
            raise ValueError("Each --insert option needs a residue token and a count")
        specs.append(parse_insertion_spec(pair[0], pair[1]))
    if not specs:
        raise ValueError("Please provide at least one insertion spec")
    return specs


def parse_specs_text(text: str) -> List[InsertionSpec]:
    """Parse GUI text with one insertion per line: A55 3 or A55:3."""
    specs: List[InsertionSpec] = []
    for line_number, raw_line in enumerate(text.replace(";", "\n").splitlines(), start=1):
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        line = line.replace(",", " ")
        if ":" in line and len(line.split()) == 1:
            token, count_text = line.split(":", 1)
        else:
            parts = line.split()
            if len(parts) != 2:
                raise ValueError(
                    "Line %d must be '<residue> <count>' or '<residue>:<count>': %s"
                    % (line_number, raw_line)
                )
            token, count_text = parts
        specs.append(parse_insertion_spec(token, count_text))
    if not specs:
        raise ValueError("Please enter at least one insertion spec")
    return specs


def default_output_path(input_path: Path) -> Path:
    return input_path.with_name(f"{input_path.stem}_vresi{input_path.suffix}")


def build_shift_map(specs: Sequence[InsertionSpec]) -> Dict[str, List[InsertionSpec]]:
    by_chain: Dict[str, List[InsertionSpec]] = {}
    for spec in specs:
        by_chain.setdefault(spec.chain_id, []).append(spec)
    for chain_id in by_chain:
        by_chain[chain_id].sort(key=lambda item: item.after_resseq)
    return by_chain


def shifted_resseq(chain_id: str, resseq: int, shift_map: Dict[str, List[InsertionSpec]]) -> int:
    shift = 0
    for spec in shift_map.get(chain_id, []):
        if resseq > spec.after_resseq:
            shift += spec.count
    return resseq + shift


def split_line_ending(line: str) -> Tuple[str, str]:
    if line.endswith("\r\n"):
        return line[:-2], "\r\n"
    if line.endswith("\n"):
        return line[:-1], "\n"
    return line, ""


def format_resseq(resseq: int) -> str:
    text = f"{resseq:4d}"
    if len(text) > 4:
        raise ValueError("PDB residue number is outside the 4-column field: %d" % resseq)
    return text


def replace_resseq_field(line: str, start: int, end: int, new_resseq: int) -> str:
    body, ending = split_line_ending(line)
    body = body.ljust(end)
    return body[:start] + format_resseq(new_resseq) + body[end:] + ending


def update_resseq_at_columns(
    line: str,
    chain_index: int,
    resseq_start: int,
    resseq_end: int,
    shift_map: Dict[str, List[InsertionSpec]],
) -> Tuple[str, bool, bool]:
    """Return (updated_line, parsed_endpoint, changed)."""
    body, _ending = split_line_ending(line)
    body = body.ljust(resseq_end)
    chain_id = body[chain_index] if chain_index < len(body) else ""
    field = body[resseq_start:resseq_end]
    try:
        old_resseq = int(field.strip())
    except ValueError:
        return line, False, False
    new_resseq = shifted_resseq(chain_id, old_resseq, shift_map)
    if new_resseq == old_resseq:
        return line, True, False
    return replace_resseq_field(line, resseq_start, resseq_end, new_resseq), True, True


def renumber_pdb_lines(lines: Sequence[str], specs: Sequence[InsertionSpec]) -> Tuple[List[str], RenumberStats]:
    shift_map = build_shift_map(specs)
    stats = RenumberStats(specs=list(specs))
    output: List[str] = []

    for line_number, line in enumerate(lines, start=1):
        record = line[:6].strip()
        updated = line

        if record in COORD_RECORDS:
            updated, parsed, changed = update_resseq_at_columns(updated, 21, 22, 26, shift_map)
            if parsed and changed:
                stats.coordinate_lines_changed += 1
        elif record == "TER":
            updated, parsed, changed = update_resseq_at_columns(updated, 21, 22, 26, shift_map)
            if parsed and changed:
                stats.ter_lines_changed += 1
        elif record == "LINK":
            parsed_any = False
            for chain_index, start, end in ((21, 22, 26), (51, 52, 56)):
                updated, parsed, changed = update_resseq_at_columns(updated, chain_index, start, end, shift_map)
                parsed_any = parsed_any or parsed
                if parsed and changed:
                    stats.link_endpoints_changed += 1
            if not parsed_any:
                stats.warnings.append("Line %d: could not parse LINK residue fields" % line_number)

        output.append(updated)

    return output, stats


def insert_virtual_residues(
    input_pdb: Path,
    specs: Sequence[InsertionSpec],
    output_pdb: Optional[Path] = None,
    verbose: bool = True,
) -> Tuple[Path, RenumberStats]:
    input_pdb = Path(input_pdb)
    if output_pdb is None:
        output_pdb = default_output_path(input_pdb)
    else:
        output_pdb = Path(output_pdb)

    lines = input_pdb.read_text(errors="replace").splitlines(True)
    updated_lines, stats = renumber_pdb_lines(lines, specs)
    output_pdb.write_text("".join(updated_lines))

    if verbose:
        print(format_summary(output_pdb, stats), flush=True)
    return output_pdb, stats


def build_cli_command(script_name: str, input_pdb: str, specs: Sequence[InsertionSpec], output_pdb: str) -> str:
    parts: List[str] = [sys.executable, script_name, input_pdb]
    for spec in specs:
        parts.extend(["--insert", spec.token, str(spec.count)])
    if output_pdb:
        parts.extend(["-o", output_pdb])
    return " ".join(shlex.quote(str(part)) for part in parts)


def format_summary(output_pdb: Path, stats: RenumberStats) -> str:
    lines = [
        "Insertion specs: " + ", ".join(spec.label for spec in stats.specs),
        "Coordinate-like records changed: %d" % stats.coordinate_lines_changed,
        "TER records changed: %d" % stats.ter_lines_changed,
        "LINK endpoints changed: %d" % stats.link_endpoints_changed,
        "Wrote: %s" % output_pdb,
    ]
    if stats.warnings:
        lines.append("Warnings:")
        lines.extend("  " + warning for warning in stats.warnings)
    return "\n".join(lines)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="insert_virtual_resi.py",
        description="Insert virtual residue-numbering gaps after selected residues in a PDB file.",
    )
    parser.add_argument("pdb", nargs="?", help="Input PDB file.")
    parser.add_argument(
        "-i",
        "--insert",
        nargs=2,
        action="append",
        metavar=("RESIDUE", "COUNT"),
        help="Insertion spec such as --insert A55 3. Repeat for multiple sites.",
    )
    parser.add_argument("-o", "--output", help="Output PDB file. Default: <input>_vresi.pdb")
    parser.add_argument("--gui", action="store_true", help="Open the Tk GUI.")
    parser.add_argument("-v", "--version", action="version", version=f"{TOOL_NAME} V{VERSION}")
    return parser


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
    root.geometry("780x500")

    style = ttk.Style(root)
    style.configure("ToolTitle.TLabel", font=("TkDefaultFont", 12, "bold"))
    style.configure("Tool.TLabelframe.Label", font=("TkDefaultFont", 10, "bold"))

    input_var = tk.StringVar(value=initial_path or "")
    output_var = tk.StringVar()
    status_var = tk.StringVar(value="Enter one insertion per line, for example: A55 3")

    def refresh_default_output() -> None:
        text = input_var.get().strip()
        if text and not output_var.get().strip():
            output_var.set(str(default_output_path(Path(text))))

    refresh_default_output()

    outer = ttk.Frame(root, padding=10)
    outer.pack(fill="both", expand=True)
    outer.columnconfigure(1, weight=1)

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

    ttk.Button(outer, text="Browse...", command=browse_input).grid(row=1, column=2, sticky="ew", padx=(6, 0), pady=4)

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

    ttk.Button(outer, text="Save as...", command=browse_output).grid(row=2, column=2, sticky="ew", padx=(6, 0), pady=4)

    specs_box = ttk.LabelFrame(outer, text="Virtual residue insertions", padding=8, style="Tool.TLabelframe")
    specs_box.grid(row=3, column=0, columnspan=3, sticky="nsew", pady=(8, 4))
    specs_box.columnconfigure(0, weight=1)
    specs_box.rowconfigure(1, weight=1)
    outer.rowconfigure(3, weight=1)

    ttk.Label(
        specs_box,
        text="One insertion per line: A55 3, A.55 3, 55A 3, 55.A 3, or A55:3",
    ).grid(row=0, column=0, sticky="w", pady=(0, 4))
    specs_text = scrolledtext.ScrolledText(specs_box, height=8, wrap="word")
    specs_text.grid(row=1, column=0, sticky="nsew")

    ttk.Label(outer, textvariable=status_var, wraplength=720).grid(row=4, column=0, columnspan=3, sticky="ew", pady=4)

    buttons = ttk.Frame(outer)
    buttons.grid(row=5, column=0, columnspan=3, sticky="ew", pady=(8, 0))

    def run_from_gui() -> None:
        try:
            input_pdb = input_var.get().strip()
            if not input_pdb:
                raise ValueError("Please choose an input PDB file.")
            if not Path(input_pdb).is_file():
                raise ValueError("Input PDB file does not exist: %s" % input_pdb)
            specs = parse_specs_text(specs_text.get("1.0", tk.END))
            output_pdb = output_var.get().strip() or str(default_output_path(Path(input_pdb)))
            cli_cmd = build_cli_command(Path(__file__).name, input_pdb, specs, output_pdb)
            print("Equivalent CLI command:", flush=True)
            print(cli_cmd, flush=True)
            out_path, stats = insert_virtual_residues(
                Path(input_pdb),
                specs,
                output_pdb=Path(output_pdb),
                verbose=False,
            )
            summary = format_summary(out_path, stats)
            print(summary, flush=True)
            status_var.set(summary)
            messagebox.showinfo(TOOL_NAME, summary, parent=root)
        except Exception as exc:
            print("Error: %s" % exc, flush=True)
            status_var.set("Error: %s" % exc)
            messagebox.showerror(TOOL_NAME, str(exc), parent=root)

    ttk.Button(buttons, text="Run", command=run_from_gui).pack(side="left")
    ttk.Button(buttons, text="Close", command=root.destroy).pack(side="left", padx=(6, 0))

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
        specs = parse_insert_pairs(args.insert)
        cli_cmd = build_cli_command(Path(__file__).name, args.pdb, specs, args.output or "")
        print("Equivalent CLI command:", flush=True)
        print(cli_cmd, flush=True)
        insert_virtual_residues(
            Path(args.pdb),
            specs,
            output_pdb=Path(args.output) if args.output else None,
            verbose=True,
        )
        return 0
    except Exception as exc:
        print("Error: %s" % exc, file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
