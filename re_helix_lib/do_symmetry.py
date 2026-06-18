#!/usr/bin/env python3
"""
do_symmetry.py

Generate an averaged, idealized symmetric PDB model from a pseudosymmetric
homomeric assembly. This single script combines the old workflow of:
  1) symmetry_2_reorder_args.py
  2) reorder_pdb_atom_by_chain.py
  3) Chimera match/alignment
  4) average_structure_pdb_atom.py

Inputs:
  - A PDB file containing ATOM/HETATM records.
  - A symmetry definition, either as explicit symmetry-related chain groups
    or as an n-fold symmetry plus a continuous chain range.

Outputs:
  - <input>_symmetric.pdb by default, or <output_base>_symmetric.pdb.
  - Optional intermediate reordered/aligned PDB files with --keep-intermediate.

Examples:
  python do_symmetry.py model.pdb --groups ABCDMNOP EFGHQRST IJKLUVWX
  python do_symmetry.py model.pdb --fold 3 --chains A-X -o model_C3
  python do_symmetry.py model.pdb ABCDMNOP EFGHQRST IJKLUVWX --keep-intermediate
  python do_symmetry.py --gui

Notes:
  - With no command-line arguments, the Tk GUI opens automatically.
  - No third-party packages are required. The rigid-body alignment uses a
    pure-Python quaternion/Kabsch-style least-squares fit.
  - Designed to be easy to place later into azbmost/re_helix as a bundled tool.
"""

from __future__ import print_function

import argparse
import copy
import math
import os
import sys

try:
    import tkinter as tk
    from tkinter import filedialog, messagebox
except Exception:  # Tkinter is optional for CLI mode.
    tk = None
    filedialog = None
    messagebox = None

VERSION = "1.0"


class PDBRecord(object):
    """Mutable representation of one ATOM, HETATM, or TER PDB record."""

    def __init__(self, line):
        self.line = line if line.endswith("\n") else line + "\n"
        self.record_name = self.line[:6].strip()
        self.serial = _safe_int(self.line[6:11].strip(), None)
        self.name = self.line[12:16]
        self.name_stripped = self.name.strip()
        self.altloc = self.line[16:17]
        self.resname = self.line[17:20]
        self.resname_stripped = self.resname.strip()
        self.chain_id = self.line[21:22]
        self.resseq = _safe_int(self.line[22:26].strip(), None)
        self.icode = self.line[26:27]
        self.x = None
        self.y = None
        self.z = None
        if self.is_atom():
            self.x = float(self.line[30:38])
            self.y = float(self.line[38:46])
            self.z = float(self.line[46:54])

    def is_atom(self):
        return self.record_name in ("ATOM", "HETATM")

    def is_ter(self):
        return self.record_name == "TER"

    def copy(self):
        return copy.deepcopy(self)

    def identity(self, ignore_resname=False):
        if ignore_resname:
            resname = ""
        else:
            resname = self.resname_stripped
        return (
            self.record_name,
            self.name_stripped,
            self.altloc,
            resname,
            self.chain_id,
            self.resseq,
            self.icode,
        )

    def update_serial(self, serial):
        self.serial = serial
        self.line = self.line[:6] + ("%5d" % serial) + self.line[11:]

    def update_chain_id(self, chain_id):
        if len(chain_id) != 1:
            raise ValueError("PDB chain IDs must be one character: %r" % chain_id)
        self.chain_id = chain_id
        self.line = self.line[:21] + chain_id + self.line[22:]

    def update_resseq(self, resseq):
        self.resseq = int(resseq)
        self.line = self.line[:22] + ("%4d" % self.resseq) + self.line[26:]

    def update_resname(self, resname):
        self.resname_stripped = resname.strip()
        self.resname = "%3s" % self.resname_stripped
        self.line = self.line[:17] + self.resname + self.line[20:]

    def update_xyz(self, x, y, z):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        self.line = self.line[:30] + "%8.3f%8.3f%8.3f" % (self.x, self.y, self.z) + self.line[54:]


def _safe_int(text, default):
    try:
        return int(text)
    except Exception:
        return default


def read_pdb(path):
    """Read ATOM/HETATM/TER records from a PDB file."""
    records = []
    with open(path, "r") as handle:
        for line in handle:
            tag = line[:6].strip()
            if tag in ("ATOM", "HETATM", "TER"):
                try:
                    rec = PDBRecord(line)
                except Exception as exc:
                    raise ValueError("Could not parse PDB line:\n%s\n%s" % (line.rstrip(), exc))
                records.append(rec)
    if not records:
        raise ValueError("No ATOM/HETATM/TER records were found in %s" % path)
    return records


def write_pdb(path, records, reorder_serial=True, include_end=True):
    """Write PDB records to a file."""
    output_records = copy.deepcopy(records)
    if reorder_serial:
        serial = 1
        for rec in output_records:
            if rec.is_atom() or rec.is_ter():
                rec.update_serial(serial)
                serial += 1
    with open(path, "w") as handle:
        for rec in output_records:
            handle.write(rec.line)
        if include_end:
            handle.write("END\n")


def make_ter_from_atom(atom_record):
    """Create a simple TER record after the supplied atom record."""
    if not atom_record.is_atom():
        raise ValueError("TER record can only be created from an atom record")
    serial = atom_record.serial if atom_record.serial is not None else 0
    resname = atom_record.resname_stripped
    chain = atom_record.chain_id
    resseq = atom_record.resseq if atom_record.resseq is not None else 0
    icode = atom_record.icode if atom_record.icode else " "
    line = "TER   %5d      %3s %1s%4d%1s\n" % (serial + 1, resname, chain, resseq, icode)
    return PDBRecord(line)


def unique_chain_order(records):
    """Return the chain IDs seen in atom records, preserving first-seen order."""
    chains = []
    seen = set()
    for rec in records:
        if rec.is_atom() and rec.chain_id not in seen:
            chains.append(rec.chain_id)
            seen.add(rec.chain_id)
    return chains


def reorder_records_by_chain(records, chain_map=None, resseq_shift_map=None):
    """
    Rename chains/residue numbers, then sort records alphabetically by chain.

    chain_map maps old_chain -> new_chain.
    resseq_shift_map maps old_chain -> integer shift.
    """
    chain_map = chain_map or {}
    resseq_shift_map = resseq_shift_map or {}
    atoms = []
    for rec in records:
        if not rec.is_atom():
            continue
        new_rec = rec.copy()
        old_chain = new_rec.chain_id
        if old_chain in chain_map:
            new_rec.update_chain_id(chain_map[old_chain])
            if new_rec.resseq is not None:
                new_rec.update_resseq(new_rec.resseq + resseq_shift_map.get(old_chain, 0))
        atoms.append(new_rec)

    grouped = {}
    chain_first_order = []
    for rec in atoms:
        if rec.chain_id not in grouped:
            grouped[rec.chain_id] = []
            chain_first_order.append(rec.chain_id)
        grouped[rec.chain_id].append(rec)

    # Alphabetical chain order matches reorder_pdb_atom_by_chain.py behavior.
    sorted_chains = sorted(chain_first_order)
    out = []
    for chain in sorted_chains:
        out.extend(grouped[chain])
        out.append(make_ter_from_atom(grouped[chain][-1]))
    return out


def parse_chain_range(text):
    """Parse A-X, A:X, AX, or ABCD into a list of one-character chain IDs."""
    if text is None:
        return None
    text = text.strip()
    if not text:
        return []
    for sep in ("-", ":", ".."):
        if sep in text:
            parts = text.split(sep)
            if len(parts) != 2 or len(parts[0]) != 1 or len(parts[1]) != 1:
                raise ValueError("Invalid chain range: %s" % text)
            start, end = parts[0], parts[1]
            if ord(start) > ord(end):
                raise ValueError("Chain range must go forward alphabetically: %s" % text)
            return [chr(i) for i in range(ord(start), ord(end) + 1)]
    if len(text) == 2 and ord(text[0]) < ord(text[1]):
        # Historical style: AF means A, B, C, D, E, F.
        return [chr(i) for i in range(ord(text[0]), ord(text[1]) + 1)]
    return list(text)


def parse_groups(group_args):
    """Parse group strings supplied as ['AB', 'CD'] or ['AB,CD']."""
    if group_args is None:
        return None
    if isinstance(group_args, str):
        raw = [group_args]
    else:
        raw = list(group_args)
    if len(raw) == 1 and "," in raw[0]:
        raw = [x.strip() for x in raw[0].split(",")]
    groups = []
    for item in raw:
        item = item.strip().strip(",")
        if item:
            groups.append(list(item))
    if not groups:
        return None
    size = len(groups[0])
    if size == 0:
        raise ValueError("Empty symmetry group was provided")
    for group in groups:
        if len(group) != size:
            raise ValueError("All symmetry-related chain groups must have the same length")
    all_chains = [c for group in groups for c in group]
    if len(all_chains) != len(set(all_chains)):
        raise ValueError("The same chain ID appears more than once in the symmetry groups")
    return groups


def groups_from_fold_and_chains(fold, chain_ids):
    if fold is None or chain_ids is None:
        raise ValueError("Both --fold and --chains are required for continuous-chain mode")
    fold = int(fold)
    if fold < 2:
        raise ValueError("Symmetry fold must be at least 2")
    if len(chain_ids) % fold != 0:
        raise ValueError(
            "Number of chains (%d) is not divisible by symmetry fold (%d)" % (len(chain_ids), fold)
        )
    group_size = len(chain_ids) // fold
    return [chain_ids[i * group_size : (i + 1) * group_size] for i in range(fold)]


def make_chain_maps(groups):
    """
    Generate old_chain -> new_chain maps for all cyclic symmetry shifts.

    Returns one identity map followed by shifts 1..n-1.
    """
    n = len(groups)
    if n < 2:
        raise ValueError("At least two symmetry-related groups are needed")
    group_size = len(groups[0])
    chain_maps = []
    identity = {}
    for group in groups:
        for chain in group:
            identity[chain] = chain
    chain_maps.append(identity)
    for shift in range(1, n):
        mapping = {}
        for col in range(group_size):
            column = [groups[row][col] for row in range(n)]
            for row in range(n):
                mapping[column[row]] = column[(row + shift) % n]
        chain_maps.append(mapping)
    return chain_maps


def make_resseq_shift_maps(groups, chain_maps, default_old_resseq=1, default_new_resseq=1):
    """
    Historical scripts allowed tokens such as 1A -> 1E, where residue shifts
    could in principle be nonzero. For chain-only symmetry groups, shifts are 0.
    This function exists to keep the internal data model extensible.
    """
    shift_maps = []
    for mapping in chain_maps:
        shifts = {}
        for old_chain in mapping:
            shifts[old_chain] = int(default_new_resseq) - int(default_old_resseq)
        shift_maps.append(shifts)
    return shift_maps


def select_fit_records(reference_records, moving_records, fit_atoms="all", ignore_resname=False):
    """Return paired atom coordinates for least-squares alignment."""
    ref_atoms = []
    mov_atoms = []
    moving_by_id = {}
    for rec in moving_records:
        if rec.is_atom() and atom_matches_fit_filter(rec, fit_atoms):
            moving_by_id[rec.identity(ignore_resname=ignore_resname)] = rec
    for rec in reference_records:
        if rec.is_atom() and atom_matches_fit_filter(rec, fit_atoms):
            key = rec.identity(ignore_resname=ignore_resname)
            if key in moving_by_id:
                ref_atoms.append(rec)
                mov_atoms.append(moving_by_id[key])
    if len(ref_atoms) < 3:
        raise ValueError(
            "Fewer than 3 matching atoms were found for alignment. "
            "Try --fit-atoms all or --ignore-resname if appropriate."
        )
    fixed = [(rec.x, rec.y, rec.z) for rec in ref_atoms]
    moving = [(rec.x, rec.y, rec.z) for rec in mov_atoms]
    return fixed, moving


def atom_matches_fit_filter(rec, fit_atoms):
    name = rec.name_stripped.upper()
    mode = (fit_atoms or "all").lower()
    if mode == "all":
        return True
    if mode in ("p", "phosphorus"):
        return name == "P"
    if mode in ("ca", "calpha"):
        return name == "CA"
    if mode == "backbone":
        return name in set(["P", "OP1", "OP2", "O1P", "O2P", "O3'", "O5'", "C3'", "C4'", "C5'"])
    raise ValueError("Unknown fit atom mode: %s" % fit_atoms)


def centroid(points):
    n = float(len(points))
    return (
        sum(p[0] for p in points) / n,
        sum(p[1] for p in points) / n,
        sum(p[2] for p in points) / n,
    )


def subtract_point(p, c):
    return (p[0] - c[0], p[1] - c[1], p[2] - c[2])


def add_point(p, c):
    return (p[0] + c[0], p[1] + c[1], p[2] + c[2])


def mat_vec_mul(matrix, vector):
    return (
        matrix[0][0] * vector[0] + matrix[0][1] * vector[1] + matrix[0][2] * vector[2],
        matrix[1][0] * vector[0] + matrix[1][1] * vector[1] + matrix[1][2] * vector[2],
        matrix[2][0] * vector[0] + matrix[2][1] * vector[1] + matrix[2][2] * vector[2],
    )


def best_fit_transform(moving_points, fixed_points):
    """
    Return rotation matrix R and translation t for R*moving + t = fixed.

    Uses Horn's quaternion method with a pure-Python Jacobi eigen solver.
    """
    if len(moving_points) != len(fixed_points):
        raise ValueError("Point lists for alignment must have equal length")
    if len(moving_points) < 3:
        raise ValueError("At least 3 points are required for a stable alignment")

    cm = centroid(moving_points)
    cf = centroid(fixed_points)
    mov = [subtract_point(p, cm) for p in moving_points]
    fix = [subtract_point(p, cf) for p in fixed_points]

    # Cross-covariance S = sum(moving * fixed^T).
    sxx = sum(m[0] * f[0] for m, f in zip(mov, fix))
    sxy = sum(m[0] * f[1] for m, f in zip(mov, fix))
    sxz = sum(m[0] * f[2] for m, f in zip(mov, fix))
    syx = sum(m[1] * f[0] for m, f in zip(mov, fix))
    syy = sum(m[1] * f[1] for m, f in zip(mov, fix))
    syz = sum(m[1] * f[2] for m, f in zip(mov, fix))
    szx = sum(m[2] * f[0] for m, f in zip(mov, fix))
    szy = sum(m[2] * f[1] for m, f in zip(mov, fix))
    szz = sum(m[2] * f[2] for m, f in zip(mov, fix))

    # Quaternion characteristic matrix. Quaternion order is [w, x, y, z].
    k = [
        [sxx + syy + szz, syz - szy, szx - sxz, sxy - syx],
        [syz - szy, sxx - syy - szz, sxy + syx, szx + sxz],
        [szx - sxz, sxy + syx, -sxx + syy - szz, syz + szy],
        [sxy - syx, szx + sxz, syz + szy, -sxx - syy + szz],
    ]
    eigvals, eigvecs = jacobi_eigen_symmetric(k)
    max_index = max(range(len(eigvals)), key=lambda idx: eigvals[idx])
    q = [eigvecs[row][max_index] for row in range(4)]
    q_norm = math.sqrt(sum(v * v for v in q))
    if q_norm == 0:
        raise ValueError("Alignment failed: zero quaternion norm")
    q = [v / q_norm for v in q]
    rotation = quaternion_to_rotation(q)

    rotated_cm = mat_vec_mul(rotation, cm)
    translation = (cf[0] - rotated_cm[0], cf[1] - rotated_cm[1], cf[2] - rotated_cm[2])
    return rotation, translation


def jacobi_eigen_symmetric(matrix, max_iter=100, tolerance=1.0e-12):
    """Eigenvalues/eigenvectors of a small real symmetric matrix."""
    n = len(matrix)
    a = [[float(matrix[i][j]) for j in range(n)] for i in range(n)]
    v = [[0.0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        v[i][i] = 1.0

    for _ in range(max_iter):
        p = 0
        q = 1
        max_off = 0.0
        for i in range(n):
            for j in range(i + 1, n):
                val = abs(a[i][j])
                if val > max_off:
                    max_off = val
                    p = i
                    q = j
        if max_off < tolerance:
            break

        if abs(a[p][p] - a[q][q]) < tolerance:
            angle = math.pi / 4.0
        else:
            angle = 0.5 * math.atan2(2.0 * a[p][q], a[q][q] - a[p][p])
        c = math.cos(angle)
        s = math.sin(angle)

        app = c * c * a[p][p] - 2.0 * s * c * a[p][q] + s * s * a[q][q]
        aqq = s * s * a[p][p] + 2.0 * s * c * a[p][q] + c * c * a[q][q]
        a[p][q] = 0.0
        a[q][p] = 0.0
        a[p][p] = app
        a[q][q] = aqq

        for r in range(n):
            if r != p and r != q:
                arp = c * a[r][p] - s * a[r][q]
                arq = s * a[r][p] + c * a[r][q]
                a[r][p] = arp
                a[p][r] = arp
                a[r][q] = arq
                a[q][r] = arq

        for r in range(n):
            vrp = c * v[r][p] - s * v[r][q]
            vrq = s * v[r][p] + c * v[r][q]
            v[r][p] = vrp
            v[r][q] = vrq

    eigenvalues = [a[i][i] for i in range(n)]
    return eigenvalues, v


def quaternion_to_rotation(q):
    """Convert [w, x, y, z] quaternion to a 3x3 rotation matrix."""
    w, x, y, z = q
    return [
        [1.0 - 2.0 * (y * y + z * z), 2.0 * (x * y - z * w), 2.0 * (x * z + y * w)],
        [2.0 * (x * y + z * w), 1.0 - 2.0 * (x * x + z * z), 2.0 * (y * z - x * w)],
        [2.0 * (x * z - y * w), 2.0 * (y * z + x * w), 1.0 - 2.0 * (x * x + y * y)],
    ]


def transform_records(records, rotation, translation):
    out = []
    for rec in records:
        new_rec = rec.copy()
        if new_rec.is_atom():
            p = mat_vec_mul(rotation, (new_rec.x, new_rec.y, new_rec.z))
            p = add_point(p, translation)
            new_rec.update_xyz(p[0], p[1], p[2])
        out.append(new_rec)
    return out


def rmsd_after_fit(moving_points, fixed_points, rotation, translation):
    n = len(moving_points)
    if n == 0:
        return 0.0
    total = 0.0
    for mov, fix in zip(moving_points, fixed_points):
        p = mat_vec_mul(rotation, mov)
        p = add_point(p, translation)
        dx = p[0] - fix[0]
        dy = p[1] - fix[1]
        dz = p[2] - fix[2]
        total += dx * dx + dy * dy + dz * dz
    return math.sqrt(total / float(n))


def average_aligned_structures(aligned_structures, ignore_resname=False, allow_missing=False):
    """Average coordinates using the first aligned structure as the template."""
    template = aligned_structures[0]
    atom_maps = []
    for records in aligned_structures:
        mapping = {}
        for rec in records:
            if rec.is_atom():
                key = rec.identity(ignore_resname=ignore_resname)
                if key in mapping:
                    raise ValueError("Duplicate atom identity encountered: %r" % (key,))
                mapping[key] = rec
        atom_maps.append(mapping)

    averaged = []
    missing_count = 0
    for rec in template:
        new_rec = rec.copy()
        if rec.is_atom():
            key = rec.identity(ignore_resname=ignore_resname)
            xs = []
            ys = []
            zs = []
            for mapping in atom_maps:
                other = mapping.get(key)
                if other is None:
                    missing_count += 1
                    if not allow_missing:
                        raise ValueError(
                            "Missing atom in at least one symmetry copy: %r. "
                            "Use --allow-missing only if this is expected." % (key,)
                        )
                    continue
                xs.append(other.x)
                ys.append(other.y)
                zs.append(other.z)
            if xs:
                new_rec.update_xyz(sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs))
        averaged.append(new_rec)
    return averaged, missing_count


def derive_output_base(input_path, output_arg):
    if output_arg:
        base = output_arg
        if base.lower().endswith(".pdb"):
            base = base[:-4]
        return base
    root, _ = os.path.splitext(input_path)
    return root


def symmetrize_pdb(
    input_pdb,
    groups,
    output_base=None,
    fit_atoms="all",
    align=True,
    keep_intermediate=False,
    ignore_resname=False,
    allow_missing=False,
    verbose=True,
):
    """Main programmatic API for generating an averaged symmetric PDB."""
    input_records = read_pdb(input_pdb)
    output_base = derive_output_base(input_pdb, output_base)
    final_path = output_base + "_symmetric.pdb"

    chain_maps = make_chain_maps(groups)
    shift_maps = make_resseq_shift_maps(groups, chain_maps)

    reordered_structures = []
    for idx, chain_map in enumerate(chain_maps):
        records = reorder_records_by_chain(input_records, chain_map, shift_maps[idx])
        reordered_structures.append(records)
        if keep_intermediate:
            path = "%s_sym%02d_reordered.pdb" % (output_base, idx + 1)
            write_pdb(path, records, reorder_serial=True)
            if verbose:
                print("wrote", path)

    reference = reordered_structures[0]
    aligned_structures = [reference]
    rmsds = [0.0]

    for idx in range(1, len(reordered_structures)):
        moving_records = reordered_structures[idx]
        if align:
            fixed_points, moving_points = select_fit_records(
                reference, moving_records, fit_atoms=fit_atoms, ignore_resname=ignore_resname
            )
            rotation, translation = best_fit_transform(moving_points, fixed_points)
            aligned = transform_records(moving_records, rotation, translation)
            rmsd = rmsd_after_fit(moving_points, fixed_points, rotation, translation)
        else:
            aligned = moving_records
            rmsd = None
        aligned_structures.append(aligned)
        rmsds.append(rmsd)
        if keep_intermediate:
            path = "%s_sym%02d_aligned.pdb" % (output_base, idx + 1)
            write_pdb(path, aligned, reorder_serial=True)
            if verbose:
                if rmsd is None:
                    print("wrote", path)
                else:
                    print("wrote", path, "RMSD=%.4f Å" % rmsd)

    averaged, missing_count = average_aligned_structures(
        aligned_structures, ignore_resname=ignore_resname, allow_missing=allow_missing
    )
    write_pdb(final_path, averaged, reorder_serial=True)

    if verbose:
        print("symmetry groups:", ["".join(g) for g in groups])
        print("fit atoms:", fit_atoms)
        if align:
            for idx, value in enumerate(rmsds):
                if idx == 0:
                    continue
                print("copy %d alignment RMSD: %.4f Å" % (idx + 1, value))
        if missing_count:
            print("warning: %d missing atom placements were skipped" % missing_count)
        print("wrote", final_path)
    return final_path, rmsds


def build_groups_from_args(args):
    """Resolve explicit argparse options or legacy positional symmetry args."""
    groups = parse_groups(args.groups)
    if groups is not None:
        return groups

    if args.fold is not None or args.chains is not None:
        return groups_from_fold_and_chains(args.fold, parse_chain_range(args.chains))

    legacy = args.symmetry_args or []
    if legacy:
        if legacy[0].isdigit():
            if len(legacy) < 2:
                raise ValueError("Legacy fold mode needs: <fold> <chain_range>, for example: 3 AF")
            return groups_from_fold_and_chains(int(legacy[0]), parse_chain_range(legacy[1]))
        return parse_groups(legacy)

    raise ValueError("Please provide --groups, or --fold with --chains")


def build_gui_defaults_from_args(args):
    """Convert parsed CLI options into initial GUI field values."""
    defaults = {
        "input_pdb": args.input_pdb or "",
        "output_base": args.output or "",
        "mode": "groups",
        "groups": "AB CD EF",
        "fold": str(args.fold) if args.fold is not None else "3",
        "chains": args.chains or "AF",
        "fit_atoms": args.fit_atoms,
        "keep_intermediate": args.keep_intermediate,
        "no_align": args.no_align,
        "ignore_resname": args.ignore_resname,
        "allow_missing": args.allow_missing,
    }

    if args.groups:
        defaults["mode"] = "groups"
        defaults["groups"] = " ".join(args.groups)
    elif args.fold is not None or args.chains is not None:
        defaults["mode"] = "fold"
    elif args.symmetry_args:
        if args.symmetry_args[0].isdigit():
            defaults["mode"] = "fold"
            defaults["fold"] = args.symmetry_args[0]
            if len(args.symmetry_args) > 1:
                defaults["chains"] = args.symmetry_args[1]
        else:
            defaults["mode"] = "groups"
            defaults["groups"] = " ".join(args.symmetry_args)

    if defaults["input_pdb"] and not defaults["output_base"]:
        defaults["output_base"] = derive_output_base(defaults["input_pdb"], None)
    return defaults


def make_arg_parser():
    parser = argparse.ArgumentParser(
        description="Generate an averaged symmetric PDB from a pseudosymmetric assembly."
    )
    parser.add_argument("input_pdb", nargs="?", help="input PDB file")
    parser.add_argument(
        "symmetry_args",
        nargs="*",
        help="optional legacy symmetry syntax: '3 AF' or 'AB CD EF'",
    )
    parser.add_argument(
        "--groups",
        nargs="+",
        help="symmetry-related chain groups, e.g. --groups ABCDMNOP EFGHQRST IJKLUVWX",
    )
    parser.add_argument("--fold", type=int, help="symmetry fold for continuous-chain mode")
    parser.add_argument(
        "--chains",
        help="chain range/list for continuous-chain mode, e.g. A-X, AX, or ABCDEF",
    )
    parser.add_argument("-o", "--output", help="output base path; _symmetric.pdb will be appended")
    parser.add_argument(
        "--fit-atoms",
        default="all",
        choices=["all", "p", "phosphorus", "ca", "calpha", "backbone"],
        help="atoms used for rigid-body alignment; default: all",
    )
    parser.add_argument(
        "--no-align",
        action="store_true",
        help="skip rigid-body alignment before averaging",
    )
    parser.add_argument(
        "--keep-intermediate",
        action="store_true",
        help="write reordered and aligned intermediate PDB files",
    )
    parser.add_argument(
        "--ignore-resname",
        action="store_true",
        help="match atoms without requiring residue names to be identical",
    )
    parser.add_argument(
        "--allow-missing",
        action="store_true",
        help="average available matching atoms instead of failing on missing atoms",
    )
    parser.add_argument("--gui", action="store_true", help="open the Tk GUI")
    parser.add_argument("-v", "--version", action="version", version="do_symmetry.py V%s" % VERSION)
    return parser


def run_cli(argv):
    parser = make_arg_parser()
    args = parser.parse_args(argv)
    if args.gui or args.input_pdb is None:
        launch_gui(build_gui_defaults_from_args(args))
        return 0
    if not os.path.isfile(args.input_pdb):
        raise ValueError("Input PDB file does not exist: %s" % args.input_pdb)
    groups = build_groups_from_args(args)
    symmetrize_pdb(
        input_pdb=args.input_pdb,
        groups=groups,
        output_base=args.output,
        fit_atoms=args.fit_atoms,
        align=not args.no_align,
        keep_intermediate=args.keep_intermediate,
        ignore_resname=args.ignore_resname,
        allow_missing=args.allow_missing,
        verbose=True,
    )
    return 0


class SymmetryGUI(object):
    def __init__(self, root, defaults=None):
        self.root = root
        defaults = defaults or {}
        root.title("Do Symmetry V%s" % VERSION)
        self.input_var = tk.StringVar(value=defaults.get("input_pdb", ""))
        self.output_var = tk.StringVar(value=defaults.get("output_base", ""))
        self.mode_var = tk.StringVar(value=defaults.get("mode", "groups"))
        self.groups_var = tk.StringVar(value=defaults.get("groups", "AB CD EF"))
        self.fold_var = tk.StringVar(value=defaults.get("fold", "3"))
        self.chains_var = tk.StringVar(value=defaults.get("chains", "AF"))
        self.fit_atoms_var = tk.StringVar(value=defaults.get("fit_atoms", "all"))
        self.keep_var = tk.BooleanVar(value=bool(defaults.get("keep_intermediate", False)))
        self.no_align_var = tk.BooleanVar(value=bool(defaults.get("no_align", False)))
        self.ignore_resname_var = tk.BooleanVar(value=bool(defaults.get("ignore_resname", False)))
        self.allow_missing_var = tk.BooleanVar(value=bool(defaults.get("allow_missing", False)))
        self.status_var = tk.StringVar(value="Choose an input PDB and symmetry definition.")
        self._build()

    def _build(self):
        pad = {"padx": 6, "pady": 4}
        row = 0

        tk.Label(
            self.root,
            text="Do Symmetry V%s" % VERSION,
            font=("TkDefaultFont", 11, "bold"),
        ).grid(row=row, column=0, columnspan=3, sticky="w", padx=6, pady=(8, 4))
        row += 1

        tk.Label(self.root, text="Input PDB").grid(row=row, column=0, sticky="e", **pad)
        tk.Entry(self.root, textvariable=self.input_var, width=55).grid(row=row, column=1, sticky="we", **pad)
        tk.Button(self.root, text="Browse", command=self.browse_input).grid(row=row, column=2, **pad)
        row += 1

        tk.Label(self.root, text="Output base").grid(row=row, column=0, sticky="e", **pad)
        tk.Entry(self.root, textvariable=self.output_var, width=55).grid(row=row, column=1, sticky="we", **pad)
        tk.Button(self.root, text="Browse", command=self.browse_output).grid(row=row, column=2, **pad)
        row += 1

        mode_frame = tk.LabelFrame(self.root, text="Symmetry definition", font=("TkDefaultFont", 10, "bold"))
        mode_frame.grid(row=row, column=0, columnspan=3, sticky="we", padx=6, pady=6)
        tk.Radiobutton(mode_frame, text="Explicit chain groups", variable=self.mode_var, value="groups").grid(
            row=0, column=0, sticky="w", **pad
        )
        tk.Entry(mode_frame, textvariable=self.groups_var, width=48).grid(row=0, column=1, columnspan=3, sticky="we", **pad)
        tk.Label(mode_frame, text="Example: ABCDMNOP EFGHQRST IJKLUVWX").grid(
            row=1, column=1, columnspan=3, sticky="w", **pad
        )

        tk.Radiobutton(mode_frame, text="Fold + continuous chain range", variable=self.mode_var, value="fold").grid(
            row=2, column=0, sticky="w", **pad
        )
        tk.Label(mode_frame, text="Fold").grid(row=2, column=1, sticky="e", **pad)
        tk.Entry(mode_frame, textvariable=self.fold_var, width=8).grid(row=2, column=2, sticky="w", **pad)
        tk.Label(mode_frame, text="Chains").grid(row=2, column=3, sticky="e", **pad)
        tk.Entry(mode_frame, textvariable=self.chains_var, width=12).grid(row=2, column=4, sticky="w", **pad)
        row += 1

        options = tk.LabelFrame(self.root, text="Options", font=("TkDefaultFont", 10, "bold"))
        options.grid(row=row, column=0, columnspan=3, sticky="we", padx=6, pady=6)
        tk.Label(options, text="Fit atoms").grid(row=0, column=0, sticky="e", **pad)
        tk.OptionMenu(options, self.fit_atoms_var, "all", "p", "phosphorus", "ca", "calpha", "backbone").grid(
            row=0, column=1, sticky="w", **pad
        )
        tk.Checkbutton(options, text="Skip alignment", variable=self.no_align_var).grid(row=0, column=2, sticky="w", **pad)
        tk.Checkbutton(options, text="Keep intermediate files", variable=self.keep_var).grid(
            row=1, column=0, columnspan=2, sticky="w", **pad
        )
        tk.Checkbutton(options, text="Ignore residue names when matching", variable=self.ignore_resname_var).grid(
            row=1, column=2, columnspan=2, sticky="w", **pad
        )
        tk.Checkbutton(options, text="Allow missing atoms", variable=self.allow_missing_var).grid(
            row=2, column=0, columnspan=2, sticky="w", **pad
        )
        row += 1

        tk.Button(self.root, text="Run symmetry averaging", command=self.run).grid(
            row=row, column=0, columnspan=3, pady=8
        )
        row += 1
        tk.Label(self.root, textvariable=self.status_var, anchor="w", justify="left", wraplength=600).grid(
            row=row, column=0, columnspan=3, sticky="we", padx=6, pady=6
        )
        self.root.columnconfigure(1, weight=1)

    def browse_input(self):
        path = filedialog.askopenfilename(
            title="Choose input PDB",
            filetypes=[("PDB files", "*.pdb"), ("All files", "*")],
        )
        if path:
            self.input_var.set(path)
            if not self.output_var.get().strip():
                root, _ = os.path.splitext(path)
                self.output_var.set(root)

    def browse_output(self):
        path = filedialog.asksaveasfilename(
            title="Choose output base",
            defaultextension=".pdb",
            filetypes=[("PDB files", "*.pdb"), ("All files", "*")],
        )
        if path:
            if path.lower().endswith(".pdb"):
                path = path[:-4]
            self.output_var.set(path)

    def run(self):
        try:
            input_pdb = self.input_var.get().strip()
            if not input_pdb:
                raise ValueError("Please choose an input PDB file")
            if not os.path.isfile(input_pdb):
                raise ValueError("Input PDB file does not exist: %s" % input_pdb)

            if self.mode_var.get() == "groups":
                text = self.groups_var.get().strip()
                groups = parse_groups(text.replace(",", " ").split())
            else:
                fold = int(self.fold_var.get().strip())
                chains = parse_chain_range(self.chains_var.get().strip())
                groups = groups_from_fold_and_chains(fold, chains)

            output_base = self.output_var.get().strip() or None
            final_path, rmsds = symmetrize_pdb(
                input_pdb=input_pdb,
                groups=groups,
                output_base=output_base,
                fit_atoms=self.fit_atoms_var.get(),
                align=not self.no_align_var.get(),
                keep_intermediate=self.keep_var.get(),
                ignore_resname=self.ignore_resname_var.get(),
                allow_missing=self.allow_missing_var.get(),
                verbose=False,
            )
            rmsd_text = ""
            if not self.no_align_var.get():
                values = [r for r in rmsds[1:] if r is not None]
                if values:
                    rmsd_text = "\nAlignment RMSDs: " + ", ".join("%.4f Å" % r for r in values)
            self.status_var.set("Done. Wrote:\n%s%s" % (final_path, rmsd_text))
            messagebox.showinfo("do_symmetry.py", "Done.\n\nWrote:\n%s" % final_path)
        except Exception as exc:
            self.status_var.set("Error: %s" % exc)
            messagebox.showerror("do_symmetry.py", str(exc))


def launch_gui(defaults=None):
    if tk is None:
        raise RuntimeError("Tkinter is not available in this Python installation")
    root = tk.Tk()
    SymmetryGUI(root, defaults=defaults)
    root.mainloop()


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    try:
        return run_cli(argv)
    except Exception as exc:
        print("error: %s" % exc, file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
