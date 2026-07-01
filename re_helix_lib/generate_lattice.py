#!/usr/bin/env python3
"""
generate_lattice.py

Prepare a PDB file for a P1 crystal/lattice description from three user-provided
lattice vectors. The script writes or replaces the CRYST1 record and, by default,
rotates ATOM/HETATM coordinates into the standard PDB crystallographic Cartesian
frame, where a is along +X, b is in the XY plane, and c has positive Z.

Important V3_1 simplification:
  - No repeat ranges and no explicit coordinate expansion are performed.
  - The output remains one asymmetric/unit-cell coordinate set.
  - Non-coordinate PDB records, including REMARK, LINK, TITLE, SEQRES, etc., are
    preserved exactly as text, except CRYST1 is inserted/replaced unless disabled.

Inputs:
  - One PDB file.
  - Three direction vectors u1/u2/u3.
  - Three distances d1/d2/d3 along those directions.

Output:
  - A PDB file with one coordinate set, updated CRYST1, and rotated coordinates
    by default.

Example:
  python generate_lattice.py input.pdb \
      --u1 0.00000735294 1 0.00000735294 --d1 136 \
      --u2 0.483332 -0.339442 -0.806951 --d2 136 \
      --u3 -0.940643 -0.339398 0 --d3 136

GUI mode:
  python generate_lattice.py
  python generate_lattice.py --gui
"""

from __future__ import print_function

import argparse
import math
import os
import sys

try:
    from re_helix_lib.gui_icon import apply_optional_icon
except ImportError:  # pragma: no cover - direct script execution fallback
    from gui_icon import apply_optional_icon

try:
    import tkinter as tk
    from tkinter import filedialog, messagebox
except Exception:
    tk = None
    filedialog = None
    messagebox = None


EPS = 1.0e-10
COORD_RECORDS = ("ATOM  ", "HETATM")
TOOL_NAME = "Generate Lattice"
VERSION = "V3.1"


def norm(v):
    return math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def normalize(v, name):
    n = norm(v)
    if n == 0.0:
        raise ValueError("{} cannot be the zero vector.".format(name))
    return (v[0] / n, v[1] / n, v[2] / n)


def scale(v, s):
    return (v[0] * s, v[1] * s, v[2] * s)


def angle_deg(a, b):
    na = norm(a)
    nb = norm(b)
    if na == 0.0 or nb == 0.0:
        raise ValueError("Cannot calculate an angle involving a zero vector.")
    c = dot(a, b) / (na * nb)
    c = max(-1.0, min(1.0, c))
    return math.degrees(math.acos(c))


def determinant3(m):
    return (
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    )


def inverse3(m):
    det = determinant3(m)
    if abs(det) < EPS:
        raise ValueError(
            "The three lattice vectors are nearly coplanar or linearly dependent; "
            "a valid 3D unit cell cannot be built."
        )

    return [
        [
            (m[1][1] * m[2][2] - m[1][2] * m[2][1]) / det,
            (m[0][2] * m[2][1] - m[0][1] * m[2][2]) / det,
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) / det,
        ],
        [
            (m[1][2] * m[2][0] - m[1][0] * m[2][2]) / det,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) / det,
            (m[0][2] * m[1][0] - m[0][0] * m[1][2]) / det,
        ],
        [
            (m[1][0] * m[2][1] - m[1][1] * m[2][0]) / det,
            (m[0][1] * m[2][0] - m[0][0] * m[2][1]) / det,
            (m[0][0] * m[1][1] - m[0][1] * m[1][0]) / det,
        ],
    ]


def matmul3(a, b):
    out = [[0.0, 0.0, 0.0] for _ in range(3)]
    for i in range(3):
        for j in range(3):
            out[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j]
    return out


def transpose3(m):
    return [
        [m[0][0], m[1][0], m[2][0]],
        [m[0][1], m[1][1], m[2][1]],
        [m[0][2], m[1][2], m[2][2]],
    ]


def matvec3(m, v):
    return (
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    )


def columns_to_matrix(a_vec, b_vec, c_vec):
    return [
        [a_vec[0], b_vec[0], c_vec[0]],
        [a_vec[1], b_vec[1], c_vec[1]],
        [a_vec[2], b_vec[2], c_vec[2]],
    ]


def max_abs_orthogonality_error(m):
    mtm = matmul3(transpose3(m), m)
    max_err = 0.0
    for i in range(3):
        for j in range(3):
            expected = 1.0 if i == j else 0.0
            max_err = max(max_err, abs(mtm[i][j] - expected))
    return max_err


def is_atom_record(line):
    return line.startswith(COORD_RECORDS)


def is_anisou_record(line):
    return line.startswith("ANISOU")


def is_coordinate_related_record(line):
    return is_atom_record(line) or is_anisou_record(line) or line.startswith("TER")


def parse_xyz(line):
    try:
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
    except ValueError:
        raise ValueError("Could not parse XYZ coordinates from line:\n{}".format(line.rstrip("\n")))
    return (x, y, z)


def format_coord(value):
    text = "{:8.3f}".format(value)
    if len(text) > 8:
        raise ValueError(
            "Coordinate {:.3f} does not fit in the PDB 8.3 coordinate field. "
            "Consider translating/recentering the input coordinates or using mmCIF."
            .format(value)
        )
    return text


def update_coord_line(line, xyz):
    base = line.rstrip("\n\r")
    if len(base) < 80:
        base = base.ljust(80)

    return (
        base[:30]
        + format_coord(xyz[0])
        + format_coord(xyz[1])
        + format_coord(xyz[2])
        + base[54:]
        + "\n"
    )


def parse_anisou_tensor(line):
    base = line.rstrip("\n\r")
    if len(base) < 70:
        raise ValueError("ANISOU line is too short to parse tensor fields.")

    u11 = int(base[28:35])
    u22 = int(base[35:42])
    u33 = int(base[42:49])
    u12 = int(base[49:56])
    u13 = int(base[56:63])
    u23 = int(base[63:70])

    return [
        [float(u11), float(u12), float(u13)],
        [float(u12), float(u22), float(u23)],
        [float(u13), float(u23), float(u33)],
    ]


def format_anisou_int(value):
    rounded = int(round(value))
    text = "{:7d}".format(rounded)
    if len(text) > 7:
        raise ValueError("Rotated ANISOU value {} does not fit in a 7-column field.".format(rounded))
    return text


def update_anisou_line(line, coord_transform):
    if coord_transform is None:
        return line

    base = line.rstrip("\n\r")
    if len(base) < 80:
        base = base.ljust(80)

    try:
        tensor = parse_anisou_tensor(base)
    except Exception:
        # Preserve the ANISOU line if it cannot be parsed safely.
        return line

    # Under x' = M x, the anisotropic tensor transforms as U' = M U M^T.
    transformed = matmul3(matmul3(coord_transform, tensor), transpose3(coord_transform))
    u11 = transformed[0][0]
    u22 = transformed[1][1]
    u33 = transformed[2][2]
    u12 = transformed[0][1]
    u13 = transformed[0][2]
    u23 = transformed[1][2]

    return (
        base[:28]
        + format_anisou_int(u11)
        + format_anisou_int(u22)
        + format_anisou_int(u33)
        + format_anisou_int(u12)
        + format_anisou_int(u13)
        + format_anisou_int(u23)
        + base[70:]
        + "\n"
    )


def make_cryst1_line(a_vec, b_vec, c_vec):
    a = norm(a_vec)
    b = norm(b_vec)
    c = norm(c_vec)
    alpha = angle_deg(b_vec, c_vec)
    beta = angle_deg(a_vec, c_vec)
    gamma = angle_deg(a_vec, b_vec)

    return "CRYST1{:9.3f}{:9.3f}{:9.3f}{:7.2f}{:7.2f}{:7.2f} P 1           1\n".format(
        a, b, c, alpha, beta, gamma
    )


def build_cell_vectors(args):
    if args.d1 <= 0.0 or args.d2 <= 0.0 or args.d3 <= 0.0:
        raise ValueError("Distances d1, d2, and d3 must all be positive.")

    u1 = normalize(tuple(args.u1), "u1")
    u2 = normalize(tuple(args.u2), "u2")
    u3 = normalize(tuple(args.u3), "u3")

    a_vec = scale(u1, args.d1)
    b_vec = scale(u2, args.d2)
    c_vec = scale(u3, args.d3)
    return a_vec, b_vec, c_vec


def build_standard_cell_vectors(a_vec, b_vec, c_vec):
    """
    Return the standard PDB crystallographic Cartesian representation of the
    same unit cell: a along +X, b in the XY plane, and c with positive Z.
    """
    a = norm(a_vec)
    b = norm(b_vec)
    c = norm(c_vec)

    alpha = math.radians(angle_deg(b_vec, c_vec))
    beta = math.radians(angle_deg(a_vec, c_vec))
    gamma = math.radians(angle_deg(a_vec, b_vec))

    sin_gamma = math.sin(gamma)
    if abs(sin_gamma) < EPS:
        raise ValueError("gamma is too close to 0 or 180 degrees; cannot build a standard cell frame.")

    va = (a, 0.0, 0.0)
    vb = (b * math.cos(gamma), b * sin_gamma, 0.0)

    cx = c * math.cos(beta)
    cy = c * (math.cos(alpha) - math.cos(beta) * math.cos(gamma)) / sin_gamma
    cz_sq = c * c - cx * cx - cy * cy
    if cz_sq < -1.0e-6:
        raise ValueError("The supplied lattice vectors do not form a physically valid unit cell.")
    vc = (cx, cy, math.sqrt(max(0.0, cz_sq)))

    return va, vb, vc


def build_rotation_to_cryst_frame(source_a, source_b, source_c, allow_reflection):
    target_a, target_b, target_c = build_standard_cell_vectors(source_a, source_b, source_c)

    source = columns_to_matrix(source_a, source_b, source_c)
    target = columns_to_matrix(target_a, target_b, target_c)
    transform = matmul3(target, inverse3(source))

    det_transform = determinant3(transform)
    if det_transform < 0.0 and not allow_reflection:
        raise ValueError(
            "The supplied lattice-vector order is left-handed relative to the standard CRYST1 frame. "
            "A pure rotation cannot fix this without reflection. Try swapping two lattice directions, "
            "or rerun with --allow-reflection if reflection is acceptable."
        )

    ortho_error = max_abs_orthogonality_error(transform)
    if ortho_error > 1.0e-5:
        raise ValueError(
            "The calculated coordinate transform is not close to a pure rotation/reflection "
            "(orthogonality error {:.3g}). Please check the three lattice vectors."
            .format(ortho_error)
        )

    return transform, (target_a, target_b, target_c), det_transform, ortho_error


def derive_output_name(input_path):
    root, ext = os.path.splitext(input_path)
    if ext:
        return root + "_cryst" + ext
    return input_path + "_cryst.pdb"


def format_matrix_for_text(m):
    return (
        "[{:.10f} {:.10f} {:.10f}; {:.10f} {:.10f} {:.10f}; {:.10f} {:.10f} {:.10f}]"
        .format(
            m[0][0], m[0][1], m[0][2],
            m[1][0], m[1][1], m[1][2],
            m[2][0], m[2][1], m[2][2],
        )
    )


def build_cli_command(script_name, args):
    parts = [
        sys.executable,
        script_name,
        args.input_pdb,
    ]
    if args.output:
        parts.extend(["-o", args.output])
    parts.extend(["--u1"] + [str(value) for value in args.u1])
    parts.extend(["--d1", str(args.d1)])
    parts.extend(["--u2"] + [str(value) for value in args.u2])
    parts.extend(["--d2", str(args.d2)])
    parts.extend(["--u3"] + [str(value) for value in args.u3])
    parts.extend(["--d3", str(args.d3)])
    if not args.rotate_to_cryst_frame:
        parts.append("--no-rotate-to-cryst-frame")
    if args.allow_reflection:
        parts.append("--allow-reflection")
    if not args.update_cryst1:
        parts.append("--no-cryst1-update")
    if args.drop_conect:
        parts.append("--drop-conect")
    return " ".join(shlex_quote(part) for part in parts)


def shlex_quote(value):
    text = str(value)
    if not text:
        return "''"
    safe = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_@%+=:,./-"
    if all(ch in safe for ch in text):
        return text
    return "'" + text.replace("'", "'\"'\"'") + "'"


def format_info_summary(info):
    lines = [
        "Wrote PDB: {}".format(info["output_path"]),
        "ATOM/HETATM records processed: {}".format(info["atom_count"]),
        "ANISOU records processed: {}".format(info["anisou_count"]),
        info["cryst1"],
    ]
    if info["rotated"]:
        lines.append("Coordinates were rotated to the standard CRYST1 frame.")
        lines.append("Rotation/reflection determinant: {:.8f}".format(info["det_transform"]))
        lines.append("Orthogonality error: {:.3g}".format(info["orthogonality_error"]))
        lines.append("Rotation matrix: {}".format(info["rotation_matrix"]))
    else:
        lines.append("Coordinates were not rotated.")
    return "\n".join(lines)


def process_pdb(args):
    with open(args.input_pdb, "r", encoding="utf-8", errors="replace") as f:
        lines = f.readlines()

    if not any(is_atom_record(line) for line in lines):
        raise ValueError("No ATOM/HETATM coordinate records were found in the input PDB.")

    source_a, source_b, source_c = build_cell_vectors(args)
    cryst1_line = make_cryst1_line(source_a, source_b, source_c)

    coord_transform = None
    det_transform = None
    ortho_error = None
    target_vectors = None

    if args.rotate_to_cryst_frame:
        coord_transform, target_vectors, det_transform, ortho_error = build_rotation_to_cryst_frame(
            source_a, source_b, source_c, args.allow_reflection
        )

    out_lines = []
    cryst1_written = False
    atom_count = 0
    anisou_count = 0

    for line in lines:
        # If no CRYST1 has appeared before the coordinate section, insert it just
        # before the first coordinate-related record. Existing lines before this
        # point are preserved exactly.
        if args.update_cryst1 and (not cryst1_written) and is_coordinate_related_record(line):
            out_lines.append(cryst1_line)
            cryst1_written = True

        if line.startswith("CRYST1"):
            if args.update_cryst1:
                if not cryst1_written:
                    out_lines.append(cryst1_line)
                    cryst1_written = True
                # Extra/late CRYST1 records are skipped after replacement.
            else:
                out_lines.append(line)
            continue

        if line.startswith("CONECT") and args.drop_conect:
            continue

        if is_atom_record(line):
            xyz = parse_xyz(line)
            if coord_transform is not None:
                xyz = matvec3(coord_transform, xyz)
            out_lines.append(update_coord_line(line, xyz))
            atom_count += 1
        elif is_anisou_record(line):
            out_lines.append(update_anisou_line(line, coord_transform))
            anisou_count += 1
        else:
            out_lines.append(line)

    if args.update_cryst1 and not cryst1_written:
        out_lines.insert(0, cryst1_line)
        cryst1_written = True

    if not out_lines or not out_lines[-1].startswith("END"):
        out_lines.append("END\n")

    output_path = args.output if args.output else derive_output_name(args.input_pdb)
    with open(output_path, "w", encoding="utf-8") as f:
        f.writelines(out_lines)

    info = {
        "output_path": output_path,
        "atom_count": atom_count,
        "anisou_count": anisou_count,
        "cryst1": cryst1_line.strip(),
        "rotated": bool(args.rotate_to_cryst_frame),
        "rotation_matrix": None,
        "det_transform": det_transform,
        "orthogonality_error": ortho_error,
        "target_vectors": target_vectors,
    }

    if coord_transform is not None:
        info["rotation_matrix"] = format_matrix_for_text(coord_transform)

    return info


def parse_vector_text(text, name):
    parts = text.replace(",", " ").split()
    if len(parts) != 3:
        raise ValueError("{} must contain three numbers, for example: 1 0 0".format(name))
    return [float(parts[0]), float(parts[1]), float(parts[2])]


def run_gui(initial_path=None):
    if tk is None:
        sys.stderr.write("ERROR: Tkinter is not available in this Python installation.\n")
        return 1

    root = tk.Tk()
    root.title("{} {}".format(TOOL_NAME, VERSION))
    apply_optional_icon(root, __file__)
    root.geometry("860x720")
    root.minsize(820, 680)

    input_var = tk.StringVar(value=initial_path or "")
    output_var = tk.StringVar()
    if initial_path:
        output_var.set(derive_output_name(initial_path))

    u1_var = tk.StringVar(value="1 0 0")
    u2_var = tk.StringVar(value="0 1 0")
    u3_var = tk.StringVar(value="0 0 1")

    d1_var = tk.StringVar(value="80")
    d2_var = tk.StringVar(value="80")
    d3_var = tk.StringVar(value="80")

    rotate_to_cryst_frame_var = tk.BooleanVar(value=True)
    update_cryst1_var = tk.BooleanVar(value=True)
    allow_reflection_var = tk.BooleanVar(value=False)
    drop_conect_var = tk.BooleanVar(value=False)

    def choose_input():
        path = filedialog.askopenfilename(
            title="Choose input PDB",
            filetypes=[("PDB files", "*.pdb *.ent"), ("All files", "*")]
        )
        if path:
            input_var.set(path)
            if not output_var.get().strip():
                output_var.set(derive_output_name(path))

    def choose_output():
        path = filedialog.asksaveasfilename(
            title="Choose output PDB",
            defaultextension=".pdb",
            filetypes=[("PDB files", "*.pdb"), ("All files", "*")]
        )
        if path:
            output_var.set(path)

    def add_label_entry(row, label, var, width=48):
        tk.Label(root, text=label, anchor="e").grid(row=row, column=0, sticky="e", padx=8, pady=4)
        entry = tk.Entry(root, textvariable=var, width=width)
        entry.grid(row=row, column=1, sticky="we", padx=8, pady=4)
        return entry

    root.columnconfigure(1, weight=1)

    tk.Label(root, text="Input PDB:", anchor="e").grid(row=0, column=0, sticky="e", padx=8, pady=4)
    tk.Entry(root, textvariable=input_var, width=48).grid(row=0, column=1, sticky="we", padx=8, pady=4)
    tk.Button(root, text="Browse", command=choose_input).grid(row=0, column=2, padx=8, pady=4)

    tk.Label(root, text="Output PDB:", anchor="e").grid(row=1, column=0, sticky="e", padx=8, pady=4)
    tk.Entry(root, textvariable=output_var, width=48).grid(row=1, column=1, sticky="we", padx=8, pady=4)
    tk.Button(root, text="Browse", command=choose_output).grid(row=1, column=2, padx=8, pady=4)

    tk.Label(root, text="Direction vectors are normalized before multiplying by distances.", anchor="w").grid(
        row=2, column=1, sticky="w", padx=8, pady=(10, 4)
    )

    add_label_entry(3, "u1 / a direction, x y z:", u1_var)
    add_label_entry(4, "u2 / b direction, x y z:", u2_var)
    add_label_entry(5, "u3 / c direction, x y z:", u3_var)

    add_label_entry(6, "d1 / a distance:", d1_var)
    add_label_entry(7, "d2 / b distance:", d2_var)
    add_label_entry(8, "d3 / c distance:", d3_var)

    tk.Checkbutton(root, text="Rotate coordinates to standard CRYST1 frame", variable=rotate_to_cryst_frame_var).grid(
        row=9, column=1, sticky="w", padx=8, pady=(12, 2)
    )
    tk.Checkbutton(root, text="Write/update CRYST1 as P 1", variable=update_cryst1_var).grid(
        row=10, column=1, sticky="w", padx=8, pady=2
    )
    tk.Checkbutton(root, text="Allow reflection if lattice-vector handedness is reversed", variable=allow_reflection_var).grid(
        row=11, column=1, sticky="w", padx=8, pady=2
    )
    tk.Checkbutton(root, text="Drop CONECT records", variable=drop_conect_var).grid(
        row=12, column=1, sticky="w", padx=8, pady=2
    )

    note = (
        "The script writes one coordinate set with CRYST1. "
        "Non-coordinate records are preserved exactly, including REMARK, LINK, TITLE, and SEQRES. "
        "Rotation is performed about the global origin; choosing a different rotation center would only add "
        "a constant translation and would not change lattice-axis alignment or continuity."
    )
    tk.Label(root, text=note, wraplength=610, justify="left", anchor="w").grid(
        row=13, column=1, sticky="we", padx=8, pady=(12, 4)
    )

    def run_processing():
        try:
            input_path = input_var.get().strip()
            if not input_path:
                raise ValueError("Please choose an input PDB file.")

            output_path = output_var.get().strip()
            if not output_path:
                output_path = derive_output_name(input_path)
                output_var.set(output_path)

            args = argparse.Namespace(
                input_pdb=input_path,
                output=output_path,
                u1=parse_vector_text(u1_var.get(), "u1"),
                u2=parse_vector_text(u2_var.get(), "u2"),
                u3=parse_vector_text(u3_var.get(), "u3"),
                d1=float(d1_var.get()),
                d2=float(d2_var.get()),
                d3=float(d3_var.get()),
                rotate_to_cryst_frame=rotate_to_cryst_frame_var.get(),
                update_cryst1=update_cryst1_var.get(),
                allow_reflection=allow_reflection_var.get(),
                drop_conect=drop_conect_var.get(),
            )

            cli_cmd = build_cli_command(os.path.basename(__file__), args)
            print("Equivalent CLI command:", flush=True)
            print(cli_cmd, flush=True)
            info = process_pdb(args)
            msg = format_info_summary(info)
            print(msg, flush=True)
            messagebox.showinfo("Done", msg)
        except Exception as exc:
            print("Error: {}".format(exc), flush=True)
            messagebox.showerror("Error", str(exc))

    tk.Button(root, text="Write P1 crystallographic PDB", command=run_processing, height=2).grid(
        row=14, column=1, sticky="we", padx=8, pady=16
    )

    root.mainloop()
    return 0


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Write/replace CRYST1 for a P1 PDB and, by default, rotate coordinates "
            "so the supplied lattice vectors match the standard PDB crystallographic frame."
        )
    )

    parser.add_argument("input_pdb", help="Input PDB file.")
    parser.add_argument("-o", "--output", help="Output PDB file. Default: input_cryst.pdb")
    parser.add_argument("--gui", action="store_true", help="Open the Tk GUI.")
    parser.add_argument("-v", "--version", action="version", version="{} {}".format(TOOL_NAME, VERSION))

    parser.add_argument("--u1", nargs=3, type=float, required=True, metavar=("X", "Y", "Z"),
                        help="Direction vector for lattice direction a. It is normalized before use.")
    parser.add_argument("--u2", nargs=3, type=float, required=True, metavar=("X", "Y", "Z"),
                        help="Direction vector for lattice direction b. It is normalized before use.")
    parser.add_argument("--u3", nargs=3, type=float, required=True, metavar=("X", "Y", "Z"),
                        help="Direction vector for lattice direction c. It is normalized before use.")

    parser.add_argument("--d1", type=float, required=True, help="Distance along direction u1/a.")
    parser.add_argument("--d2", type=float, required=True, help="Distance along direction u2/b.")
    parser.add_argument("--d3", type=float, required=True, help="Distance along direction u3/c.")

    parser.add_argument("--no-rotate-to-cryst-frame", dest="rotate_to_cryst_frame", action="store_false", default=True,
                        help="Do not rotate coordinates. Default behavior is to rotate to the standard CRYST1 frame.")
    parser.add_argument("--allow-reflection", action="store_true",
                        help="Allow an improper transform if the supplied lattice-vector order is left-handed. "
                             "By default, the script refuses this because it is not a pure rotation.")
    parser.add_argument("--no-cryst1-update", dest="update_cryst1", action="store_false", default=True,
                        help="Do not insert or replace the CRYST1 record.")
    parser.add_argument("--drop-conect", action="store_true",
                        help="Drop CONECT records. By default they are preserved exactly as text.")

    return parser


def main():
    argv = sys.argv[1:]
    if len(argv) == 0:
        return run_gui()
    if "--gui" in argv:
        remaining = [item for item in argv if item != "--gui"]
        initial_path = remaining[0] if remaining else None
        return run_gui(initial_path=initial_path)

    parser = build_arg_parser()
    args = parser.parse_args(argv)

    try:
        info = process_pdb(args)
    except Exception as exc:
        sys.stderr.write("ERROR: {}\n".format(exc))
        return 1

    print("Equivalent CLI command:")
    print(build_cli_command(os.path.basename(__file__), args))
    print(format_info_summary(info))
    return 0


if __name__ == "__main__":
    sys.exit(main())
