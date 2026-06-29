#!/usr/bin/env python3
"""Hidden Tk launcher for experimental re_helix cyclic alignment methods."""

from __future__ import annotations

import os
import queue
import shlex
import subprocess
import sys
import threading
from pathlib import Path
import tkinter as tk
from tkinter import filedialog, messagebox, ttk


SCRIPT_DIR = Path(__file__).resolve().parent
HELP_BUTTON_BG = "#cfeeff"
HELP_BUTTON_ACTIVE_BG = "#b8e3ff"


METHODS = {
    "ccg": {
        "label": "CCG",
        "script": "re_helix_ccgV3_1.py",
        "help": (
            "CCG uses a weighted global geometric least-squares residual as the "
            "primary cyclic objective. It starts from the regular re_helix tree "
            "alignment, then refines cyclic components with 6 degrees of freedom "
            "per non-root helix. It is usually the simpler geometry-polish path."
        ),
    },
    "cck": {
        "label": "CCK",
        "script": "re_helix_cckV3_1.py",
        "help": (
            "CCK uses an explicit closure residual on cycle edges as the primary "
            "objective. It builds a BFS tree, solves per-tree-edge phi_offset/beta "
            "variables with scipy.root when square and least_squares otherwise, "
            "then ranks by closure first and geometry second. It can also run the "
            "CCG polish and write the twist diagnostic TSV."
        ),
    },
}


def _shlex_join(cmd: list[str]) -> str:
    if hasattr(shlex, "join"):
        return shlex.join(cmd)
    return " ".join(shlex.quote(part) for part in cmd)


class CyclicAlignmentLauncher:
    def __init__(self, root: tk.Tk) -> None:
        self.root = root
        self.root.title("re_helix cyclic alignment V3.1")
        self.root.geometry("1040x650")
        self.root.minsize(880, 520)

        self.output_queue: queue.Queue[tuple[str, object]] = queue.Queue()
        self.proc: subprocess.Popen[str] | None = None

        self.method_var = tk.StringVar(value="cck")
        self.pdb_var = tk.StringVar()
        self.output_var = tk.StringVar()
        self.axis_dist_var = tk.StringVar(value="22.0")
        self.axis_parallel_var = tk.StringVar(value="n")
        self.fix_var = tk.StringVar()
        self.replicate_var = tk.BooleanVar(value=False)
        self.cir_shift_var = tk.StringVar(value="8")
        self.axis_count_var = tk.StringVar(value="1")
        self.axis_row_target = 1
        self.axis_render_pending = False
        self.axis_widgets: list[dict[str, object]] = []
        self.w_pp_var = tk.StringVar(value="1.0")
        self.w_line_var = tk.StringVar(value="1.0")
        self.w_axis_var = tk.StringVar(value="10000")
        self.extra_args_var = tk.StringVar()

        self.geom_attempts_var = tk.StringVar(value="4")
        self.geom_max_nfev_var = tk.StringVar(value="1200")

        self.max_root_attempts_var = tk.StringVar(value="20")
        self.root_maxfev_var = tk.StringVar(value="800")
        self.closure_rel_tol_var = tk.StringVar(value="0.10")
        self.closure_abs_tol_var = tk.StringVar(value="0.001")
        self.geom_polish_var = tk.BooleanVar(value=True)
        self.twist_report_var = tk.BooleanVar(value=True)
        self.twist_report_file_var = tk.StringVar()
        self.twist_helical_repeat_var = tk.StringVar(value="10.5")
        self.twist_pairing_var = tk.StringVar(value="consecutive")

        self._build_ui()
        self.method_var.trace_add("write", lambda *_: self._refresh_method_options())
        self.axis_count_var.trace_add("write", lambda *_: self._schedule_axis_rows())
        self._refresh_method_options()
        self.root.after(100, self._drain_output_queue)

    def _build_ui(self) -> None:
        outer = ttk.Frame(self.root, padding=10)
        outer.grid(row=0, column=0, sticky="nsew")
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        outer.columnconfigure(0, weight=1)
        outer.rowconfigure(4, weight=1)

        top = ttk.LabelFrame(outer, text="Method")
        top.grid(row=0, column=0, sticky="ew", pady=(0, 8))
        top.columnconfigure(4, weight=1)

        ttk.Radiobutton(top, text="CCK", variable=self.method_var, value="cck").grid(row=0, column=0, padx=(8, 2), pady=8)
        self._help_button(top, METHODS["cck"]["help"]).grid(row=0, column=1, padx=(0, 16), pady=8)
        ttk.Radiobutton(top, text="CCG", variable=self.method_var, value="ccg").grid(row=0, column=2, padx=(0, 2), pady=8)
        self._help_button(top, METHODS["ccg"]["help"]).grid(row=0, column=3, padx=(0, 16), pady=8)
        self.script_label = ttk.Label(top, text="")
        self.script_label.grid(row=0, column=4, sticky="e", padx=8)

        form = ttk.LabelFrame(outer, text="Inputs")
        form.grid(row=1, column=0, sticky="ew", pady=(0, 8))
        form.columnconfigure(1, weight=1)

        self._entry_row(
            form,
            0,
            "Input PDB",
            self.pdb_var,
            browse=lambda: self._browse_file(self.pdb_var, [("PDB files", "*.pdb *.pdb.txt"), ("All files", "*")]),
            help_text="Input PDB file to align.",
            entry_colspan=8,
            help_col=11,
        )
        self._entry_row(
            form,
            1,
            "Output base",
            self.output_var,
            browse=lambda: self._browse_save_base(self.output_var),
            help_text="Optional output base path. The method-specific suffixes are added automatically.",
            entry_colspan=8,
            help_col=11,
        )

        ttk.Label(form, text="Exchange ops").grid(row=2, column=0, sticky="nw", padx=(8, 6), pady=4)
        self.ops_text = tk.Text(form, height=3, wrap="word", undo=True)
        self.ops_text.grid(row=2, column=1, columnspan=9, sticky="ew", pady=4)
        self._help_button(
            form,
            "Tokens may include helix definitions and reciprocal-exchange specs, for example: (AB) (CD) 30A 8D d 13B 24C s",
        ).grid(row=2, column=10, sticky="n", padx=(6, 8), pady=4)

        run_options = ttk.LabelFrame(outer, text="Options")
        run_options.grid(row=2, column=0, sticky="ew", pady=(0, 8))
        run_options.columnconfigure(1, weight=1)
        run_options.columnconfigure(7, weight=1)

        self._inline_entry(run_options, 0, 0, "axis_dist", self.axis_dist_var, help_text="Target helix-axis distance in Angstroms.", width=8)
        self._inline_combo(
            run_options,
            0,
            3,
            "axis_parallel",
            self.axis_parallel_var,
            ["n", "y"],
            "Use n for cyclic tilt/refinement; y keeps the original tree alignment behavior.",
            width=5,
        )
        self._inline_entry(run_options, 0, 6, "fix", self.fix_var, help_text="Optional chain ID whose helix should remain fixed.", width=5)
        ttk.Checkbutton(run_options, text="replicate", variable=self.replicate_var).grid(row=0, column=9, sticky="w", padx=(10, 2), pady=4)
        self._help_button(run_options, "Replicate all chains using the same semantics as re_helix.py.").grid(row=0, column=10, sticky="w", padx=(0, 8), pady=4)

        self._inline_entry(run_options, 1, 0, "cir_shift", self.cir_shift_var, help_text="Residue shift for circular strands during reciprocal exchange.", width=8)
        self._inline_entry(run_options, 1, 3, "w_pp", self.w_pp_var, help_text="Weight for phosphate-pair distance residuals or scores.", width=8)
        self._inline_entry(run_options, 1, 6, "w_line", self.w_line_var, help_text="Weight for line-topology residuals or scores.", width=8)
        self._inline_entry(run_options, 1, 9, "w_axis", self.w_axis_var, help_text="Weight for target axis-distance mismatch.", width=9)

        axis_box = ttk.LabelFrame(outer, text="Axis definitions")
        axis_box.grid(row=3, column=0, sticky="ew", pady=(0, 8))
        axis_box.columnconfigure(0, weight=1)

        axis_header = ttk.Frame(axis_box)
        axis_header.grid(row=0, column=0, sticky="ew", padx=8, pady=(8, 0))
        ttk.Label(axis_header, text="Rows").pack(side="left")
        validate_count_cmd = (self.root.register(self._validate_count_text), "%P")
        tk.Spinbox(
            axis_header,
            from_=0,
            to=30,
            width=5,
            textvariable=self.axis_count_var,
            validate="key",
            validatecommand=validate_count_cmd,
            command=self._render_axis_rows,
        ).pack(side="left", padx=6)
        self._help_button(
            axis_header,
            "Each row pairs an axis definition with optional additional chains or residues that move with that axis. Examples: axis A,B with move C,D; or axis B26-B60,A1-A35 with move C1-C50,D.",
        ).pack(side="left")
        ttk.Label(axis_header, text="Axis: A,B or B26-B60,A1-A35; move: C,D or C1-C50,D").pack(side="left", padx=8)

        self.axis_rows_frame = ttk.Frame(axis_box)
        self.axis_rows_frame.grid(row=1, column=0, sticky="ew", padx=8, pady=8)
        self.axis_rows_frame.columnconfigure(1, weight=1)
        self.axis_rows_frame.columnconfigure(2, weight=1)
        ttk.Label(self.axis_rows_frame, text="").grid(row=0, column=0, sticky="w")
        ttk.Label(self.axis_rows_frame, text="axis definition").grid(row=0, column=1, sticky="w", padx=4)
        ttk.Label(self.axis_rows_frame, text="move with axis").grid(row=0, column=2, sticky="w", padx=4)
        self._render_axis_rows()

        self._entry_row(
            run_options,
            2,
            "extra args",
            self.extra_args_var,
            help_text="Optional raw command-line arguments appended after the launcher-generated options.",
            entry_colspan=9,
            help_col=10,
        )

        options = ttk.Frame(outer)
        options.grid(row=4, column=0, sticky="nsew", pady=(0, 8))
        options.columnconfigure(0, weight=1)
        options.rowconfigure(2, weight=1)

        self.ccg_frame = ttk.LabelFrame(options, text="CCG options")
        self.ccg_frame.grid(row=0, column=0, sticky="ew", pady=(0, 8))
        self.ccg_frame.columnconfigure(1, weight=1)
        self._inline_entry(self.ccg_frame, 0, 0, "geom_attempts", self.geom_attempts_var, help_text="Number of geometry least-squares starts per cyclic component.", width=8)
        self._inline_entry(self.ccg_frame, 0, 3, "geom_max_nfev", self.geom_max_nfev_var, help_text="Maximum function evaluations in the coarse geometry stage.", width=8)

        self.cck_frame = ttk.LabelFrame(options, text="CCK options")
        self.cck_frame.grid(row=1, column=0, sticky="ew", pady=(0, 8))
        self.cck_frame.columnconfigure(1, weight=1)
        self._inline_entry(self.cck_frame, 0, 0, "max_root_attempts", self.max_root_attempts_var, help_text="Number of deterministic closure-root starts per cyclic component.", width=8)
        self._inline_entry(self.cck_frame, 0, 3, "root_maxfev", self.root_maxfev_var, help_text="Maximum function evaluations per closure solve attempt.", width=8)
        self._inline_entry(self.cck_frame, 0, 6, "closure_rel_tol", self.closure_rel_tol_var, help_text="Relative tolerance for the near-best closure candidate set.", width=8)
        self._inline_entry(self.cck_frame, 0, 9, "closure_abs_tol", self.closure_abs_tol_var, help_text="Absolute tolerance for the near-best closure candidate set.", width=8)
        ttk.Checkbutton(self.cck_frame, text="geom_polish", variable=self.geom_polish_var).grid(row=1, column=0, sticky="w", padx=(8, 2), pady=4)
        self._help_button(self.cck_frame, "Run final CCG geometry polish after the CCK closure solve.").grid(row=1, column=1, sticky="w", pady=4)
        ttk.Checkbutton(self.cck_frame, text="twist_report", variable=self.twist_report_var).grid(row=1, column=2, sticky="w", padx=(12, 2), pady=4)
        self._help_button(
            self.cck_frame,
            "Write the post-alignment twist TSV. twist_rod is measured from axis-to-axis co-perpendicular connector vectors; twist_helix is measured from RE-site phosphate radial vectors.",
        ).grid(row=1, column=3, sticky="w", pady=4)
        self._inline_combo(self.cck_frame, 1, 4, "twist pairing", self.twist_pairing_var, ["consecutive", "all"], "Report consecutive exchange-site pairs or all pairwise combinations on each helix.", width=12)
        self._inline_entry(self.cck_frame, 1, 8, "helical repeat", self.twist_helical_repeat_var, help_text="Fallback bp/turn used only when observed phosphate-step information is unavailable.", width=8)
        self._entry_row(
            self.cck_frame,
            2,
            "twist file",
            self.twist_report_file_var,
            browse=lambda: self._browse_save_base(self.twist_report_file_var, default_ext=".tsv"),
            help_text="Optional TSV output path for the twist report.",
        )

        log_frame = ttk.LabelFrame(options, text="Run log")
        log_frame.grid(row=2, column=0, sticky="nsew")
        log_frame.columnconfigure(0, weight=1)
        log_frame.rowconfigure(0, weight=1)
        self.log_text = tk.Text(log_frame, height=10, wrap="word", state="disabled")
        self.log_text.grid(row=0, column=0, sticky="nsew")
        scrollbar = ttk.Scrollbar(log_frame, orient="vertical", command=self.log_text.yview)
        scrollbar.grid(row=0, column=1, sticky="ns")
        self.log_text.configure(yscrollcommand=scrollbar.set)

        buttons = ttk.Frame(outer)
        buttons.grid(row=5, column=0, sticky="ew")
        buttons.columnconfigure(0, weight=1)
        self.run_button = ttk.Button(buttons, text="Run", command=self._run)
        self.run_button.grid(row=0, column=1, padx=(8, 0))
        self.stop_button = ttk.Button(buttons, text="Stop", command=self._stop, state="disabled")
        self.stop_button.grid(row=0, column=2, padx=(8, 0))

    def _help_button(self, parent: tk.Widget, text: str) -> tk.Button:
        return tk.Button(
            parent,
            text="?",
            width=2,
            bg=HELP_BUTTON_BG,
            activebackground=HELP_BUTTON_ACTIVE_BG,
            relief="raised",
            command=lambda: messagebox.showinfo("Help", text, parent=self.root),
        )

    def _inline_entry(
        self,
        parent: ttk.Frame,
        row: int,
        col: int,
        label: str,
        variable: tk.StringVar,
        browse=None,
        help_text: str | None = None,
        width: int = 10,
    ) -> None:
        ttk.Label(parent, text=label).grid(row=row, column=col, sticky="w", padx=(8, 4), pady=3)
        entry = ttk.Entry(parent, textvariable=variable, width=width)
        entry.grid(row=row, column=col + 1, sticky="ew", pady=3)
        next_col = col + 2
        if browse is not None:
            ttk.Button(parent, text="Browse", command=browse).grid(row=row, column=next_col, sticky="ew", padx=(4, 0), pady=3)
            next_col += 1
        if help_text is not None:
            self._help_button(parent, help_text).grid(row=row, column=next_col, sticky="w", padx=(4, 8), pady=3)

    def _inline_combo(
        self,
        parent: ttk.Frame,
        row: int,
        col: int,
        label: str,
        variable: tk.StringVar,
        values: list[str],
        help_text: str,
        width: int = 8,
    ) -> None:
        ttk.Label(parent, text=label).grid(row=row, column=col, sticky="w", padx=(8, 4), pady=3)
        combo = ttk.Combobox(parent, textvariable=variable, values=values, state="readonly", width=width)
        combo.grid(row=row, column=col + 1, sticky="w", pady=3)
        self._help_button(parent, help_text).grid(row=row, column=col + 2, sticky="w", padx=(4, 8), pady=3)

    def _entry_row(
        self,
        parent: ttk.Frame,
        row: int,
        label: str,
        variable: tk.StringVar,
        browse=None,
        help_text: str | None = None,
        entry_colspan: int | None = None,
        help_col: int | None = None,
    ) -> None:
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", padx=(0, 6), pady=4)
        span = entry_colspan if entry_colspan is not None else (3 if browse is None else 2)
        entry = ttk.Entry(parent, textvariable=variable)
        entry.grid(row=row, column=1, columnspan=span, sticky="ew", pady=4)
        if browse is not None:
            ttk.Button(parent, text="Browse", command=browse).grid(row=row, column=1 + span, sticky="ew", padx=(6, 0), pady=4)
        if help_text is not None:
            if help_col is None:
                help_col = 5 if entry_colspan is None else 1 + span
                if browse is not None:
                    help_col += 2
            self._help_button(parent, help_text).grid(row=row, column=help_col, sticky="w", padx=(6, 0), pady=4)

    def _validate_count_text(self, proposed: str) -> bool:
        return proposed == "" or proposed.isdigit()

    def _coerce_count(self, raw_value: str, minimum: int, maximum: int) -> int | None:
        value_text = str(raw_value).strip()
        if not value_text or not value_text.isdigit():
            return None
        return max(minimum, min(maximum, int(value_text)))

    def _current_axis_target(self) -> int:
        return self.axis_row_target

    def _schedule_axis_rows(self) -> None:
        if self.axis_render_pending:
            return
        self.axis_render_pending = True

        def _run() -> None:
            self.axis_render_pending = False
            self._render_axis_rows()

        self.root.after_idle(_run)

    def _render_axis_rows(self) -> None:
        if not hasattr(self, "axis_rows_frame"):
            return
        target_opt = self._coerce_count(self.axis_count_var.get(), 0, 30)
        if target_opt is None:
            return
        target = target_opt
        self.axis_row_target = target
        while len(self.axis_widgets) < target:
            row_index = len(self.axis_widgets)
            axis_var = tk.StringVar()
            move_var = tk.StringVar()
            label = ttk.Label(self.axis_rows_frame, text=f"Axis {row_index + 1}")
            axis_entry = ttk.Entry(self.axis_rows_frame, textvariable=axis_var, width=32)
            move_entry = ttk.Entry(self.axis_rows_frame, textvariable=move_var, width=32)
            self.axis_widgets.append(
                {
                    "label_widget": label,
                    "axis_var": axis_var,
                    "move_var": move_var,
                    "widgets": [axis_entry, move_entry],
                }
            )

        for idx, item in enumerate(self.axis_widgets):
            visible = idx < target
            label_widget = item["label_widget"]
            widgets = item["widgets"]
            grid_row = idx + 1
            if visible:
                label_widget.configure(text=f"Axis {idx + 1}")
                label_widget.grid(row=grid_row, column=0, sticky="w", pady=2)
                for col, widget in enumerate(widgets, start=1):
                    widget.grid(row=grid_row, column=col, sticky="ew", padx=4, pady=2)
            else:
                label_widget.grid_remove()
                for widget in widgets:
                    widget.grid_remove()

    def _combo_row(
        self,
        parent: ttk.Frame,
        row: int,
        label: str,
        variable: tk.StringVar,
        values: list[str],
        help_text: str,
    ) -> None:
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", padx=(0, 6), pady=4)
        combo = ttk.Combobox(parent, textvariable=variable, values=values, state="readonly", width=16)
        combo.grid(row=row, column=1, sticky="w", pady=4)
        self._help_button(parent, help_text).grid(row=row, column=5, sticky="w", padx=(6, 0), pady=4)

    def _browse_file(self, variable: tk.StringVar, filetypes: list[tuple[str, str]]) -> None:
        path = filedialog.askopenfilename(parent=self.root, filetypes=filetypes)
        if path:
            variable.set(path)

    def _browse_save_base(self, variable: tk.StringVar, default_ext: str = ".pdb") -> None:
        path = filedialog.asksaveasfilename(parent=self.root, defaultextension=default_ext)
        if path:
            variable.set(path)

    def _refresh_method_options(self) -> None:
        method = self.method_var.get()
        info = METHODS.get(method, METHODS["cck"])
        self.script_label.configure(text=info["script"])
        if method == "ccg":
            self.ccg_frame.grid()
            self.cck_frame.grid_remove()
        else:
            self.ccg_frame.grid()
            self.cck_frame.grid()

    def _append_log(self, text: str) -> None:
        self.log_text.configure(state="normal")
        self.log_text.insert("end", text)
        self.log_text.see("end")
        self.log_text.configure(state="disabled")

    def _clear_log(self) -> None:
        self.log_text.configure(state="normal")
        self.log_text.delete("1.0", "end")
        self.log_text.configure(state="disabled")

    def _parse_ops(self) -> list[str]:
        ops_raw = self.ops_text.get("1.0", "end").strip()
        if not ops_raw:
            raise ValueError("Exchange ops are required.")
        return shlex.split(ops_raw)

    def _add_option(self, cmd: list[str], name: str, value: str) -> None:
        value = value.strip()
        if value:
            cmd.extend([name, value])

    def _build_command(self) -> list[str]:
        method = self.method_var.get()
        info = METHODS.get(method)
        if info is None:
            raise ValueError("Choose CCG or CCK.")

        script = SCRIPT_DIR / info["script"]
        if not script.exists():
            raise FileNotFoundError("Could not find %s" % script)

        pdb_path = self.pdb_var.get().strip()
        if not pdb_path:
            raise ValueError("Input PDB is required.")
        if not Path(pdb_path).exists():
            raise FileNotFoundError("Input PDB does not exist: %s" % pdb_path)

        ops = self._parse_ops()
        cmd = [sys.executable, str(script)]

        self._add_option(cmd, "-o", self.output_var.get())
        self._add_option(cmd, "--axis_dist", self.axis_dist_var.get())
        self._add_option(cmd, "--axis_parallel", self.axis_parallel_var.get())
        self._add_option(cmd, "--fix", self.fix_var.get())
        if self.replicate_var.get():
            cmd.append("--replicate")
        self._add_option(cmd, "--cir_shift", self.cir_shift_var.get())
        for idx, item in enumerate(self.axis_widgets[: self._current_axis_target()]):
            axis_range = item["axis_var"].get().strip()
            axis_move = item["move_var"].get().strip()
            if axis_move and not axis_range:
                raise ValueError(f"Axis {idx + 1}: move with axis requires an axis definition.")
            if axis_range:
                cmd.extend(["--axis_range", axis_range])
                if axis_move:
                    cmd.extend(["--axis_move", axis_move])
        self._add_option(cmd, "--w_pp", self.w_pp_var.get())
        self._add_option(cmd, "--w_line", self.w_line_var.get())
        self._add_option(cmd, "--w_axis", self.w_axis_var.get())

        self._add_option(cmd, "--geom_attempts", self.geom_attempts_var.get())
        self._add_option(cmd, "--geom_max_nfev", self.geom_max_nfev_var.get())

        if method == "cck":
            self._add_option(cmd, "--max_root_attempts", self.max_root_attempts_var.get())
            self._add_option(cmd, "--root_maxfev", self.root_maxfev_var.get())
            self._add_option(cmd, "--closure_rel_tol", self.closure_rel_tol_var.get())
            self._add_option(cmd, "--closure_abs_tol", self.closure_abs_tol_var.get())
            cmd.extend(["--geom_polish", "y" if self.geom_polish_var.get() else "n"])
            cmd.extend(["--twist_report", "y" if self.twist_report_var.get() else "n"])
            self._add_option(cmd, "--twist_report_file", self.twist_report_file_var.get())
            self._add_option(cmd, "--twist_helical_repeat", self.twist_helical_repeat_var.get())
            self._add_option(cmd, "--twist_pairing", self.twist_pairing_var.get())

        extra = self.extra_args_var.get().strip()
        if extra:
            cmd.extend(shlex.split(extra))

        cmd.append(pdb_path)
        cmd.extend(ops)
        return cmd

    def _run(self) -> None:
        try:
            cmd = self._build_command()
        except Exception as exc:
            messagebox.showerror("Cannot run", str(exc), parent=self.root)
            return

        self._clear_log()
        self._append_log("$ %s\n\n" % _shlex_join(cmd))
        self.run_button.configure(state="disabled")
        self.stop_button.configure(state="normal")

        thread = threading.Thread(target=self._worker, args=(cmd,), daemon=True)
        thread.start()

    def _worker(self, cmd: list[str]) -> None:
        cwd = os.getcwd()
        try:
            pdb_path = Path(self.pdb_var.get().strip())
            if pdb_path.is_absolute() and pdb_path.parent.exists():
                cwd = str(pdb_path.parent)
        except Exception:
            pass

        try:
            self.proc = subprocess.Popen(
                cmd,
                cwd=cwd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
            )
            assert self.proc.stdout is not None
            for line in self.proc.stdout:
                self.output_queue.put(("log", line))
            code = self.proc.wait()
            self.output_queue.put(("done", code))
        except Exception as exc:
            self.output_queue.put(("error", exc))
        finally:
            self.proc = None

    def _stop(self) -> None:
        if self.proc is not None and self.proc.poll() is None:
            self.proc.terminate()
            self._append_log("\nRequested process termination.\n")

    def _drain_output_queue(self) -> None:
        try:
            while True:
                kind, payload = self.output_queue.get_nowait()
                if kind == "log":
                    self._append_log(str(payload))
                elif kind == "done":
                    code = int(payload)
                    self.run_button.configure(state="normal")
                    self.stop_button.configure(state="disabled")
                    self._append_log("\nProcess exited with code %d.\n" % code)
                    if code == 0:
                        messagebox.showinfo("Run complete", "Cyclic alignment finished successfully.", parent=self.root)
                    else:
                        messagebox.showerror("Run failed", "Process exited with code %d." % code, parent=self.root)
                elif kind == "error":
                    self.run_button.configure(state="normal")
                    self.stop_button.configure(state="disabled")
                    self._append_log("\nError: %s\n" % payload)
                    messagebox.showerror("Run failed", str(payload), parent=self.root)
        except queue.Empty:
            pass
        self.root.after(100, self._drain_output_queue)


def main() -> None:
    root = tk.Tk()
    CyclicAlignmentLauncher(root)
    root.mainloop()


if __name__ == "__main__":
    main()
