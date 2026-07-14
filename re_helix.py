#!/usr/bin/env python3
"""
re_helix.py

V3.20 update (2026-07-13):
- Allow a full XYZ translation in user-defined-axis alignment mode. For every
  candidate rotation angle, use the exact centroid-matching translation of the
  selected P-atom pairs.
- Bump the re_helix app version to V3.20.

V3.19 update (2026-07-13):
- In user-defined-axis alignment mode, optimize translation along the supplied
  axis as well as rotation around it.
- Bump the re_helix app version to V3.19.

V3.18 update (2026-07-09):
- Add a user-defined alignment axis option. When a direction vector and point
  are supplied, re_helix skips fixed/moving helix-axis estimation for alignment
  and optimizes a single rotation of the movable helix around that line.
- Add GUI controls with dynamic greying for the user-defined axis mode.
- Bump the re_helix app version to V3.18.

V3.17 update (2026-06-28):
- Allow --axis_range terms to be whole-chain letters such as A,B in addition
  to residue windows such as A1-A35,B60-B26.
- Add paired --axis_move definitions so additional whole chains or residue
  windows can move with a user-defined axis, avoiding triplex stdin prompts
  and supporting broader axis-coupled payloads.
- Rename alignment-angle terminology: theta is now tau, and rho is now beta.
  Legacy theta/rho wording is still accepted where it affects existing input
  syntax or imported helper compatibility.
- Bump the re_helix app version to V3.17.

V3.16 update (2026-06-26):
- Add Get Phenix Restraints to the GUI Other tools area. The bundled
  re_helix_lib/get_phenix_restraints.py tool generates Phenix LINK geometry
  params, optional junction movement-selection params, and linker support files.
- Bump the re_helix app version to V3.16.

V3.15 update (2026-06-25):
- In GUI mode, keep the default Output base synchronized with the selected or
  typed input PDB path unless the user has entered a custom output base.
- In GUI mode, add a CLI-style Exchange pairs line. When filled, it is used in
  place of the individual pair rows for large reciprocal-exchange specs.
- Bump the re_helix app version to V3.15.

V3.14 update (2026-06-25):
- Add a configurable 3'-3' linker phosphate residue style for bowtie
  exchanges. The default remains HETATM X33; users can supply a custom
  residue name or use ATOM DA for Phenix-friendly relaxation without a custom
  residue definition.
- Bump the re_helix app version to V3.14.

V3.13 update (2026-06-24):
- Internal package maintenance for bundled library utilities.
- Bump the re_helix app version to V3.13.

V3.12 update (2026-06-23):
- Add Generate Lattice to the GUI Other tools area. The bundled
  re_helix_lib/generate_lattice.py tool prepares a P1 CRYST1 lattice record
  from three user-provided lattice directions and distances.
- Bump the re_helix app version to V3.12.

V3.11 update (2026-06-20):
- Add Insert Virtual Resi to the GUI Other tools area. The bundled
  re_helix_lib/insert_virtual_resi.py tool inserts residue-numbering gaps after
  selected residues and updates LINK record endpoints with the same mapping.
- Bump the re_helix app version to V3.11.

V3.10 update (2026-06-18):
- Mirror stdout/stderr from bundled Other tools into the main re_helix run log
  after launch, so tool commands and run summaries remain visible even when
  the helper window has no log box.
- Bump the re_helix app version to V3.10.

V3.9 update (2026-06-18):
- Add the Add PDB LINK Record tool to the GUI Other tools area. The button
  launches the bundled re_helix_lib/add_pdb_link_record.py GUI for staging PDB
  LINK records and rebuilding chain topology.
- Bump the re_helix app version to V3.9.

V3.8 update (2026-06-18):
- Add Do Symmetry to the GUI Other tools area. The button launches the bundled
  re_helix_lib/do_symmetry.py GUI for generating averaged symmetric PDB models.
- Bump the re_helix app version to V3.8.

V3.7 update (2026-06-16):
- Add a Bend Helix tool button to the GUI. The button launches the bundled
  re_helix_lib/bend_helix.py GUI for bending a straight two-chain helix.
- Bump the re_helix app version to V3.7.

V3.6 update (2026-06-16):
- Rename the main script to re_helix.py and move helper modules/resources into
  re_helix_lib/.
- Add --re_only / --re-only mode, which applies reciprocal exchanges directly
  to the input PDB without the helix-alignment stage.
- Load the optional GUI/task-menu icon from assets/icon.png when present;
  the script still runs normally if the icon asset is absent.
- Add -v / --version to print the app version from the command line.

V3.5 update (2026-06-14):
- Add optional fixed beta-angle syntax for single-site inter-helix reciprocal
  exchanges: <pos1> <pos2> <beta_deg> <kind>, e.g. 26A 9C 90 d.
  This was originally documented as the legacy rho-angle syntax.
  When exactly one reciprocal-exchange site connects a helix pair and
  --axis_parallel n is used, beta is held at the requested angle while
  d/tau/phi are optimised.  If the same helix pair has more than one
  exchange site, or if --axis_parallel y is used, the fixed-beta definition
  is ignored with a warning.
- This update helps with generating junctions with controlled interhelical
  angles.

V3.4 update (2026-05-23):
- Load reciprocal_exchange_pdbV3_3.py and write parse-friendly REMARK 950
  RE_SCRIPT header metadata into output PDB files, including software/version,
  developer, command line, chain residue ranges, junction residue IDs for the
  reciprocal-exchange output, and special topology events.
- The aligned+RE output now carries X33 HET/HETNAM records generated by the
  reciprocal-exchange layer for standalone 3'-3' linker phosphates.

V3.3 update (2026-03-29):
- In GUI mode, typing directly into the pair-row and axis-range row-count
  controls now refreshes the displayed rows immediately; arrow-button changes
  still work as before.
- Wrap the GUI contents in a vertically scrollable container and show a
  scrollbar automatically when the content height exceeds the current window
  height, so lower controls remain reachable for large numbers of rows.

V3.2 update (2026-03-29):
- In GUI mode, blank or incomplete alignment-pair rows are ignored rather
  than raising an immediate error; blank axis-range rows are also ignored.
  At least one complete alignment pair is still required to run.
- Give GUI help buttons (?) a light blue shade so they are easier to
  distinguish from the normal input controls.
- Improve the GUI completion message so exit code 0 is reported explicitly
  as a successful run.

V3.1 update (2026-03-29):
- Add GUI help buttons (?) for the main input fields and option groups. Each
  help button opens a pop-up with a short explanation and examples.
- Add explicit unit labels in the GUI for fields such as axis_dist (Å),
  cir_shift (nt), alignment positions (nt+chain), and axis residue ranges
  (nt).
- Change replication handling for --axis_range definitions. These definitions
  are no longer automatically copied to every replicated helix. Instead, in
  replication mode they are interpreted against the final post-replication
  chain IDs so the user can target replicated chains explicitly.

V3 update (2026-03-28):
- Add manual helical-axis residue-range definitions via repeatable --axis_range
  arguments such as "B26-B60,A1-A35". When a helix-pair alignment uses
  P positions that all fall inside one of these user-defined ranges, the axis
  is determined from all P atoms in the specified residue windows rather than
  from the whole helix model.
- Add an optional Tk GUI mode. If no command-line arguments are supplied, or
  if --gui is given, the script opens a GUI for entering the same parameters
  interactively. The GUI supports dynamic numbers of alignment-pair rows and
  axis-range rows, and its run log prints the equivalent CLI command before
  launching the job.

V2.2 update (2026-03-18):
- Detect potential triplex inputs. When exactly three chains with P atoms are
  present, the script treats them as a potential triplex helix, prompts for
  the third-strand chain ID, and uses the other two strands for helical-axis
  determination. Triplex axis-strand assignments are propagated through
  replication so copied helices keep the same duplex-based axis definition.

V2.1 update (2026-03-11):
- Fix for --axis_parallel y: treat helix axes as parallel OR anti-parallel.
  The script now explores both possibilities by optionally flipping the
  moving helix 180 degrees about a perpendicular axis through its helix
  center, runs the fine optimisation for both, and automatically selects
  the better (lower sum(dist^2)) solution.

Align pairs of nucleic-acid helices (or multi-chain helix groups) in a PDB
based on user-specified cross-helix P-atom pairs, and then apply reciprocal-
exchange operations using reciprocal_exchange_pdbV3_3.py.

Highlights
----------
- Input:
    * PDB file with nucleic-acid strands.
    * Reciprocal-exchange-style specification (same as reciprocal_exchange_pdbV3_3.py),
      optionally preceded by explicit helix definitions.  Each operation can
      also include an optional beta angle before the kind, e.g. 26A 9C 90 d.
      The angle is used only for a single-site helix pair when
      --axis_parallel n is selected; it is otherwise ignored with a warning.

          (AB) (CD) 30A 8D d 13B 24C s ...
          (ABMN) 30A 8D d ...
          26A 9C 90 d      # optional fixed beta angle for a single-site pair

      (AB) means chains A and B form one helix group; (ABMN) means chains A, B,
      M, and N form a single rigid helix group. If any such tokens are present,
      helix pairing is NOT auto-detected (unless replication overrides it).

- Helix pairing:
    * If explicit ( ... ) style tokens are present (e.g. (AB), (ABMN)): use them
      as helix groups (each group moves as a rigid block). If a helix group
      contains three chains, the script asks which chain is the third strand
      and uses the other two chains for axis estimation.
    * Otherwise: auto-detect helices as chain pairs whose P atoms are closest
      on average (each auto-detected helix has two chains). If exactly three
      chains with P atoms are present, a potential triplex is detected, the
      user is asked for the third-strand chain ID, and the three chains are
      treated as one rigid triplex helix group.

- Replication mode (replicate all chains):
    * Replication can be triggered in two ways:
        - If the input PDB appears to contain exactly one helix component
          (based on P-atom analysis), replication is enabled automatically; or
        - You explicitly request it with --replicate (in which case it is
          applied regardless of how many helices are present).
    * In replicate mode:
        - ALL chains present in the input are first renamed to consecutive
          letters A, B, C, ... according to the alphabetical order of their
          original chain IDs.
        - Based on the chain IDs appearing in the exchange specification, the
          script determines how many full copies of this base set of chains
          are needed so that all requested chain IDs exist.
        - It then replicates the entire base set of chains that many times
          (copy 0: A.., copy 1: next block of letters, etc.).
        - Helix group definitions (from explicit ( ... ) definitions or from
          auto-detected helices) are propagated across all copies, so that
          each copy has matching helix groups (e.g. (ABMN), (EFQR), ...).
        - Exchange specs are rewritten with uppercase chain IDs consistent
          with the new replicated structure and passed to the alignment stage.

- Alignment per helix pair:
    * Uses P-atom pairs connecting those two helix groups.
    * If there are 1 or 2 P-atom pairs for the helix pair, we use all of them.
    * If there are >2 P-atom pairs, we first select the two P pairs whose
      helix-1 P atoms are most separated along the helix-1 axis (i.e. the
      extreme P positions in projection onto the axis), and we perform the
      alignment objective using ONLY these two P pairs; the others are ignored.
    * **For axis estimation we use only the chains inside each helix group
      that actually appear in these P-pairs.** For example, with helix groups
      (ABCDEFGH) and (IJKLMNOP) and P-pairs 35G–27J, 33H–29I, the axes are
      based on P atoms from chains G,H and I,J respectively, but the resulting
      rigid transform is applied to all chains in each group.

- Degrees of freedom:
    First we run align_axes_for_pair() so helix 1 and 2 axes are parallel
    (direction u) and separated by --axis_dist.

    Then for each helix pair:

      If --axis_parallel y (default):
        - Helices remain parallel.
        - Optimise parameters [d, tau, phi]:
             d     : slide helix 2 along u
             tau   : rotate helix 2 around its own axis
             phi   : rotate helix 2 around helix 1's axis

      If --axis_parallel n:
        - Additionally allow helices to tilt (become non-parallel).
        - Optimise parameters [d, tau, phi, beta], where:
             beta  : rotate helix 2 around a line L_beta that:
                       - is perpendicular to the (initial) common axis direction,
                       - is constructed from the selected P pair(s) on helix 1
                         (see below),
                       - intersects helix 2's axis at the nearest point.
        - For a helix pair with exactly one reciprocal-exchange site, the
          optional syntax <pos1> <pos2> <beta_deg> <kind> fixes beta to that
          user-defined angle while d/tau/phi are optimised.  This helps
          generate junctions with controlled interhelical angles.

        Construction of L_beta:
          * Let u be the common axis direction.
          * For the selected P pair(s), take the helix-1 P coordinates.
            With two selected pairs these are the two extremes along u; with
            one selected pair the single P coordinate is used directly.
          * Their midpoint/centroid is P1_mid.
          * Project P1_mid onto helix-1 axis through center_fixed:
                anchor1 = center_fixed + ((P1_mid - center_fixed)·u) u
          * After applying (phi, d), let C2_phi_d be a point on helix-2 axis.
          * The point on helix 2’s axis nearest to anchor1 is:
                C2_base = C2_phi_d + ((anchor1 - C2_phi_d)·u) u
          * r_perp = C2_base - anchor1 is perpendicular to u.
          * L_beta is the line through anchor1 with direction r_perp.

- Helix graph:
    * Build a graph where nodes are helix groups and edges indicate at least
      one P-atom pair between them.
    * Process each connected component as a tree:
        - pick a root helix (fixed),
        - BFS to align each new helix only once relative to an already-fixed
          neighbour,
        - skip edges between already-fixed helices (cycle edges) so we never
          re-move a helix that’s been aligned.
        - With --fix A, the helix containing chain A is used as the fixed
          root for its component (and is never moved).

- Output:
    * Always writes:  <base>_aligned.pdb     (after alignment, pre-RE)
    * Always writes:  <base>_aligned_rex.pdb (after applying reciprocal
      exchanges to the aligned structure).
    * With --re_only / --re-only, writes: <base>_rex.pdb (reciprocal exchange
      directly on the input structure, with no alignment output).
      <base> is from -o/--output (extension stripped) or from input filename
      (without .pdb).
      The legacy --re flag is kept for backward compatibility but is no longer
      required (reciprocal exchanges are always applied).
"""

import argparse
import importlib
import math
import queue
import shlex
import subprocess
import sys
import threading
from collections import defaultdict
from typing import Dict, List, Tuple, Iterable, Set, Optional

from re_helix_lib.edit_pdb_atom import (
    file2rec,
    pdb_atom_record,
)

import contextlib
import importlib.util
from pathlib import Path

SOFTWARE_NAME = "re_helix"
SOFTWARE_VERSION = "V3.20"
SOFTWARE_DEVELOPER = "DiLiuLab"
APP_TITLE = (
    "re_helix V3.20: AZBMOST Package Module #2 - "
    "Align Helices and Performing Reciprocal Exchanges"
)


def _load_rex_module():
    """Load the reciprocal_exchange_pdbV3_3 helper as a Python module.

    V3.6 requires the V3.3 reciprocal-exchange helper because that module
    provides X33 linker and header-remark utilities used by this script.
    """
    # 1) If there is an already-importable module name, use it.
    for modname in (
        "re_helix_lib.reciprocal_exchange_pdbV3_3",
        "reciprocal_exchange_pdbV3_3",
    ):
        try:
            return importlib.import_module(modname)
        except Exception:
            pass

    # 2) Otherwise, load from a file path in the helper folder or next to this script.
    here = Path(__file__).resolve().parent
    for path in (
        here / "re_helix_lib" / "reciprocal_exchange_pdbV3_3.py",
        here / "re_helix_lib" / "reciprocal_exchange_pdbV3.3.py",
        here / "reciprocal_exchange_pdbV3_3.py",
        here / "reciprocal_exchange_pdbV3.3.py",
    ):
        if not path.exists():
            continue
        spec = importlib.util.spec_from_file_location("reciprocal_exchange_pdbV3_3", str(path))
        if spec is None or spec.loader is None:
            raise ImportError(f"Could not create an import spec for {path}")
        module = importlib.util.module_from_spec(spec)
        # Make the module visible during execution (required by dataclasses, etc.)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)  # type: ignore[attr-defined]
        return module

    raise ImportError(
        "Could not locate reciprocal_exchange_pdbV3_3.py in re_helix_lib/."
    )


rex = _load_rex_module()


# ---------------------------------------------------------------------------
# Exchange-spec parsing (V2-compatible, but with bowtie support)
# ---------------------------------------------------------------------------

def parse_position_token(token: str) -> Tuple[str, int]:
    """Parse a residue position token into (chainID, resSeq).

    Allows formats:
        "30A", "A30", "A.30", "30.A" -> ("A", 30)

    This matches the legacy reciprocal_exchange_pdbV2 parser behavior so the
    rest of re_helix logic remains unchanged.
    """
    t = token.strip()
    t = t.replace(".", "")
    if len(t) < 2:
        raise ValueError(f"Invalid position token '{token}' (too short).")

    # Pattern 1: A30, B8, etc. (letter followed by digits)
    if t[0].isalpha() and t[1:].isdigit():
        chain_id = t[0]
        res_seq = int(t[1:])
        return chain_id, res_seq

    # Pattern 2: 30A, 8D, etc. (digits followed by letter)
    if t[-1].isalpha() and t[:-1].isdigit():
        chain_id = t[-1]
        res_seq = int(t[:-1])
        return chain_id, res_seq

    raise ValueError(
        f"Invalid position token '{token}': expected forms like '30A', 'A30', 'A.30', '30.A'."
    )


def normalize_kind(kind_token: str) -> str:
    """Normalize exchange-kind token to one of: double / single / bowtie."""
    k = kind_token.strip().lower()
    if k in ("double", "d"):
        return "double"
    if k in ("single", "s"):
        return "single"
    if k in ("bowtie", "b"):
        return "bowtie"
    raise ValueError(
        f"Invalid exchange kind '{kind_token}': expected 'double'/'d', 'single'/'s', or 'bowtie'/'b'."
    )


ExchangeSpec = Tuple[Tuple[str, int], Tuple[str, int], str, Optional[float]]

ANGLE_DEFINITIONS_MESSAGE = (
    "Angle definitions: tau = axial twist/spin of moving helix; "
    "phi = orbital azimuth around fixed helix; "
    "beta = interhelical tilt/bend; d = axial slide."
)
_LEGACY_RHO_ALIAS_NOTE_PRINTED = False


def write_angle_definitions_once(prefix: str = "[re_helix]") -> None:
    sys.stderr.write(f"{prefix} {ANGLE_DEFINITIONS_MESSAGE}\n")


def note_legacy_rho_alias_once(prefix: str = "[re_helix]") -> None:
    global _LEGACY_RHO_ALIAS_NOTE_PRINTED
    if _LEGACY_RHO_ALIAS_NOTE_PRINTED:
        return
    _LEGACY_RHO_ALIAS_NOTE_PRINTED = True
    sys.stderr.write(
        f"{prefix} Note: the optional numeric angle field formerly documented as rho "
        "is now beta; interpreting it as beta.\n"
    )


def parse_beta_angle_token(token: str) -> float:
    """Parse an optional fixed beta angle token, in degrees."""
    try:
        value = float(token.strip())
    except ValueError as exc:
        raise ValueError(
            f"Invalid optional beta angle '{token}': expected a numeric degree value "
            f"before the exchange kind, e.g. '26A 9C 90 d'."
        ) from exc
    if not math.isfinite(value):
        raise ValueError(
            f"Invalid optional beta angle '{token}': expected a finite numeric degree value."
        )
    return value


def parse_rho_angle_token(token: str) -> float:
    """Legacy alias for parse_beta_angle_token()."""
    note_legacy_rho_alias_once()
    return parse_beta_angle_token(token)


def unpack_exchange_spec(spec) -> ExchangeSpec:
    """Return (pos1, pos2, kind, beta_deg_or_None) for V3.4/V3.5-style specs."""
    if len(spec) == 3:
        pos1, pos2, kind = spec
        return pos1, pos2, kind, None
    if len(spec) == 4:
        pos1, pos2, kind, beta_deg = spec
        return pos1, pos2, kind, beta_deg
    raise ValueError(f"Invalid exchange spec record with {len(spec)} fields: {spec!r}")


def exchange_spec_without_beta(spec) -> Tuple[Tuple[str, int], Tuple[str, int], str]:
    """Return the reciprocal-exchange part of a spec, dropping alignment-only beta metadata."""
    pos1, pos2, kind, _beta_deg = unpack_exchange_spec(spec)
    return pos1, pos2, kind


def exchange_spec_without_rho(spec) -> Tuple[Tuple[str, int], Tuple[str, int], str]:
    """Legacy alias for exchange_spec_without_beta()."""
    note_legacy_rho_alias_once()
    return exchange_spec_without_beta(spec)


def parse_exchange_specs(
    tokens: List[str],
) -> List[ExchangeSpec]:
    """Parse exchange specs into ((chain1,res1),(chain2,res2),kind,beta_deg).

    Accepted per-operation formats:
        <pos1> <pos2> <kind>
        <pos1> <pos2> <beta_deg> <kind>

    The optional beta angle is alignment-only metadata.  It is used only for
    single-site helix pairs under --axis_parallel n, and is ignored by the
    reciprocal-exchange topology layer.
    """
    if not tokens:
        raise ValueError("No reciprocal-exchange specification tokens provided.")

    ops: List[ExchangeSpec] = []
    i = 0
    while i < len(tokens):
        if i + 2 >= len(tokens):
            raise ValueError(
                "Incomplete exchange specification; expected either "
                "<pos1> <pos2> <kind> or <pos1> <pos2> <beta_deg> <kind>."
            )

        p1 = parse_position_token(tokens[i])
        p2 = parse_position_token(tokens[i + 1])
        third = tokens[i + 2]

        try:
            kind = normalize_kind(third)
            beta_deg: Optional[float] = None
            i += 3
        except ValueError:
            if i + 3 >= len(tokens):
                raise ValueError(
                    "Incomplete exchange specification after optional beta angle; "
                    "expected <pos1> <pos2> <beta_deg> <kind>."
                )
            beta_deg = parse_beta_angle_token(third)
            note_legacy_rho_alias_once()
            kind = normalize_kind(tokens[i + 3])
            i += 4

        ops.append((p1, p2, kind, beta_deg))

    return ops


Point3D = Tuple[float, float, float]
UserAxisDefinition = Tuple[Point3D, Point3D]  # (axis_dir, axis_point)
HelixID = Tuple[str, ...]  # can be 2 or more chains, e.g. ('A', 'B', 'M', 'N')
AxisBounds = Optional[Tuple[int, int]]
AxisSelection = Dict[str, AxisBounds]
ResolvedAxisRange = Dict[str, Tuple[int, int]]
MoveSelection = Dict[str, AxisBounds]


_HELIX_AXIS_OVERRIDES: Dict[HelixID, HelixID] = {}
_HELIX_AXIS_RANGE_DEFINITIONS: List[ResolvedAxisRange] = []
_HELIX_AXIS_MOVE_DEFINITIONS: List[MoveSelection] = []
_HELIX_AXIS_RANGE_OVERRIDES: Dict[HelixID, ResolvedAxisRange] = {}
_HELIX_MOVE_SELECTIONS: Dict[HelixID, MoveSelection] = {}


def _canonical_helix_id(chains: Iterable[str]) -> HelixID:
    return tuple(sorted(chains))


def _format_chain_group(chains: Iterable[str]) -> str:
    return "(" + "".join(sorted(chains)) + ")"


def reset_helix_axis_overrides() -> None:
    _HELIX_AXIS_OVERRIDES.clear()
    _HELIX_AXIS_RANGE_DEFINITIONS.clear()
    _HELIX_AXIS_MOVE_DEFINITIONS.clear()
    _HELIX_AXIS_RANGE_OVERRIDES.clear()
    _HELIX_MOVE_SELECTIONS.clear()


def get_helix_axis_overrides() -> Dict[HelixID, HelixID]:
    return dict(_HELIX_AXIS_OVERRIDES)


def set_helix_axis_overrides(overrides: Optional[Dict[HelixID, HelixID]]) -> None:
    _HELIX_AXIS_OVERRIDES.clear()
    if not overrides:
        return

    for helix_id, axis_chains in overrides.items():
        key = _canonical_helix_id(helix_id)
        value = tuple(sorted(ch for ch in axis_chains if ch in key))
        if key and value:
            _HELIX_AXIS_OVERRIDES[key] = value


def _copy_axis_range_definitions(
    defs: Optional[List[ResolvedAxisRange]],
) -> List[ResolvedAxisRange]:
    copied: List[ResolvedAxisRange] = []
    if not defs:
        return copied
    for spec in defs:
        copied.append(
            {
                str(ch): (int(bounds[0]), int(bounds[1]))
                for ch, bounds in spec.items()
            }
        )
    return copied


def _copy_move_selection(selection: Optional[MoveSelection]) -> MoveSelection:
    copied: MoveSelection = {}
    if not selection:
        return copied
    for ch, bounds in selection.items():
        copied[str(ch)] = None if bounds is None else (int(bounds[0]), int(bounds[1]))
    return copied


def get_helix_axis_range_definitions() -> List[ResolvedAxisRange]:
    return _copy_axis_range_definitions(_HELIX_AXIS_RANGE_DEFINITIONS)


def set_helix_axis_range_definitions(
    defs: Optional[List[ResolvedAxisRange]],
) -> None:
    _HELIX_AXIS_RANGE_DEFINITIONS.clear()
    if not defs:
        return

    for spec in defs:
        normalized: Dict[str, Tuple[int, int]] = {}
        for ch, bounds in spec.items():
            if len(str(ch)) != 1:
                raise ValueError(f"Axis-range chain ID '{ch}' must be a single character.")
            lo = int(bounds[0])
            hi = int(bounds[1])
            if lo > hi:
                lo, hi = hi, lo
            normalized[str(ch)] = (lo, hi)
        if normalized:
            _HELIX_AXIS_RANGE_DEFINITIONS.append(dict(sorted(normalized.items())))


def get_helix_axis_move_definitions() -> List[MoveSelection]:
    return [_copy_move_selection(selection) for selection in _HELIX_AXIS_MOVE_DEFINITIONS]


def set_helix_axis_move_definitions(defs: Optional[List[MoveSelection]]) -> None:
    _HELIX_AXIS_MOVE_DEFINITIONS.clear()
    if not defs:
        return

    for spec in defs:
        normalized: MoveSelection = {}
        for ch, bounds in spec.items():
            if len(str(ch)) != 1:
                raise ValueError(f"Axis-move chain ID '{ch}' must be a single character.")
            if bounds is None:
                normalized[str(ch)] = None
            else:
                lo = int(bounds[0])
                hi = int(bounds[1])
                if lo > hi:
                    lo, hi = hi, lo
                normalized[str(ch)] = (lo, hi)
        _HELIX_AXIS_MOVE_DEFINITIONS.append(dict(sorted(normalized.items())))


def get_helix_axis_range_overrides() -> Dict[HelixID, ResolvedAxisRange]:
    return {h: dict(spec) for h, spec in _HELIX_AXIS_RANGE_OVERRIDES.items()}


def set_helix_axis_range_overrides(
    overrides: Optional[Dict[HelixID, ResolvedAxisRange]],
) -> None:
    _HELIX_AXIS_RANGE_OVERRIDES.clear()
    if not overrides:
        return

    for helix_id, spec in overrides.items():
        key = _canonical_helix_id(helix_id)
        normalized: ResolvedAxisRange = {}
        for ch, bounds in spec.items():
            if str(ch) not in key:
                continue
            lo = int(bounds[0])
            hi = int(bounds[1])
            if lo > hi:
                lo, hi = hi, lo
            normalized[str(ch)] = (lo, hi)
        if key and normalized:
            _HELIX_AXIS_RANGE_OVERRIDES[key] = dict(sorted(normalized.items()))


def get_helix_move_selections() -> Dict[HelixID, MoveSelection]:
    return {h: _copy_move_selection(selection) for h, selection in _HELIX_MOVE_SELECTIONS.items()}


def set_helix_move_selections(
    selections: Optional[Dict[HelixID, MoveSelection]],
) -> None:
    _HELIX_MOVE_SELECTIONS.clear()
    if not selections:
        return

    for helix_id, selection in selections.items():
        key = _canonical_helix_id(helix_id)
        normalized: MoveSelection = {}
        for ch, bounds in selection.items():
            if str(ch) not in key:
                continue
            if bounds is None:
                normalized[str(ch)] = None
            else:
                lo = int(bounds[0])
                hi = int(bounds[1])
                if lo > hi:
                    lo, hi = hi, lo
                normalized[str(ch)] = (lo, hi)
        if key and normalized:
            _HELIX_MOVE_SELECTIONS[key] = dict(sorted(normalized.items()))


def get_axis_chains_for_helix(helix_id: HelixID) -> HelixID:
    key = _canonical_helix_id(helix_id)
    return _HELIX_AXIS_OVERRIDES.get(key, key)


def atom_moves_with_helix(atom: pdb_atom_record, helix_id: HelixID) -> bool:
    """Return whether atom should receive transforms assigned to helix_id."""
    key = _canonical_helix_id(helix_id)
    if atom.chainID not in key:
        return False
    selection = _HELIX_MOVE_SELECTIONS.get(key)
    if not selection:
        return True
    if atom.chainID not in selection:
        return False
    bounds = selection[atom.chainID]
    if bounds is None:
        return True
    lo, hi = bounds
    return lo <= atom.resSeq <= hi


def _resolve_chain_id_for_helix(chain_token: str, helix_id: HelixID) -> str:
    lookup = {ch.upper(): ch for ch in helix_id}
    token = chain_token.strip()
    if not token:
        raise ValueError(
            f"Please provide one of the chain IDs in {_format_chain_group(helix_id)}."
        )

    key = token[0].upper()
    if key not in lookup:
        raise ValueError(
            f"Chain '{token}' is not part of {_format_chain_group(helix_id)}."
        )
    return lookup[key]


def _parse_axis_range_term(term: str) -> Tuple[str, AxisBounds]:
    token = term.strip().replace(".", "")
    if not token:
        raise ValueError("Empty axis-range term.")
    if "-" not in token and len(token) == 1 and token.isalpha():
        return token, None
    if "-" not in token:
        raise ValueError(
            f"Invalid axis-range term '{term}': expected a chain ID like 'B' or a residue range like 'B26-B60'."
        )

    left, right = token.split("-", 1)
    chain1, res1 = parse_position_token(left)

    right = right.strip()
    if not right:
        raise ValueError(
            f"Invalid axis-range term '{term}': missing end residue after '-'."
        )

    if any(ch.isalpha() for ch in right):
        chain2, res2 = parse_position_token(right)
    else:
        chain2, res2 = chain1, int(right)

    if chain1.upper() != chain2.upper():
        raise ValueError(
            f"Axis-range term '{term}' must stay on one chain; got '{chain1}' and '{chain2}'."
        )

    lo = min(res1, res2)
    hi = max(res1, res2)
    return chain1, (lo, hi)


def parse_axis_range_spec(spec: str) -> AxisSelection:
    parts = [part.strip() for part in spec.split(",") if part.strip()]
    if not parts:
        raise ValueError("Axis-range specification cannot be empty.")

    result: AxisSelection = {}
    for part in parts:
        chain, bounds = _parse_axis_range_term(part)
        if chain in result:
            raise ValueError(
                f"Axis-range specification '{spec}' defines chain '{chain}' more than once."
            )
        result[chain] = bounds

    return dict(sorted(result.items()))


def parse_axis_range_specs(specs: Iterable[str]) -> List[AxisSelection]:
    parsed: List[AxisSelection] = []
    for spec in specs:
        parsed.append(parse_axis_range_spec(spec))
    return parsed


def parse_axis_move_specs(specs: Iterable[str]) -> List[MoveSelection]:
    parsed: List[MoveSelection] = []
    for spec in specs:
        parsed.append(parse_axis_range_spec(spec))
    return parsed


def pair_axis_move_definitions(
    axis_defs: List[AxisSelection],
    move_defs: List[MoveSelection],
) -> List[Tuple[AxisSelection, MoveSelection]]:
    if len(move_defs) > len(axis_defs):
        raise ValueError(
            "More --axis_move definitions were provided than --axis_range definitions. "
            "Each --axis_move row is paired with the --axis_range row at the same index."
        )

    paired: List[Tuple[AxisSelection, MoveSelection]] = []
    for idx, axis_def in enumerate(axis_defs):
        move_def = move_defs[idx] if idx < len(move_defs) else {}
        paired.append((axis_def, move_def))
    return paired


def resolve_axis_range_definitions(
    defs: List[AxisSelection],
    available_chains: Iterable[str],
    chain_to_P_atoms: Optional[Dict[str, List[pdb_atom_record]]] = None,
) -> List[ResolvedAxisRange]:
    lookup = {str(ch).upper(): str(ch) for ch in available_chains}
    resolved: List[ResolvedAxisRange] = []

    for spec in defs:
        resolved_spec: ResolvedAxisRange = {}
        for ch, bounds in spec.items():
            key = str(ch).upper()
            if key not in lookup:
                raise ValueError(
                    f"Axis-range definition {_format_axis_range_spec(spec)} refers to chain '{ch}', "
                    "but that chain is not present in the input PDB."
                )
            actual_chain = lookup[key]
            if actual_chain in resolved_spec:
                raise ValueError(
                    f"Axis-range definition {_format_axis_range_spec(spec)} resolves chain '{ch}' "
                    f"to '{actual_chain}' more than once."
                )
            if bounds is None:
                p_atoms = list(chain_to_P_atoms.get(actual_chain, [])) if chain_to_P_atoms else []
                if not p_atoms:
                    raise ValueError(
                        f"Axis-range definition {_format_axis_range_spec(spec)} uses whole chain '{ch}', "
                        "but no P atoms are available on that chain to define an axis."
                    )
                lo = min(atom.resSeq for atom in p_atoms)
                hi = max(atom.resSeq for atom in p_atoms)
            else:
                lo = int(bounds[0])
                hi = int(bounds[1])
                if lo > hi:
                    lo, hi = hi, lo
            resolved_spec[actual_chain] = (lo, hi)
        resolved.append(dict(sorted(resolved_spec.items())))

    return resolved


def resolve_axis_move_definitions(
    defs: List[MoveSelection],
    available_chains: Iterable[str],
) -> List[MoveSelection]:
    lookup = {str(ch).upper(): str(ch) for ch in available_chains}
    resolved: List[MoveSelection] = []

    for spec in defs:
        resolved_spec: MoveSelection = {}
        for ch, bounds in spec.items():
            key = str(ch).upper()
            if key not in lookup:
                raise ValueError(
                    f"Axis-move definition {_format_axis_range_spec(spec)} refers to chain '{ch}', "
                    "but that chain is not present in the input PDB."
                )
            actual_chain = lookup[key]
            if actual_chain in resolved_spec:
                raise ValueError(
                    f"Axis-move definition {_format_axis_range_spec(spec)} resolves chain '{ch}' "
                    f"to '{actual_chain}' more than once."
                )
            if bounds is None:
                resolved_spec[actual_chain] = None
            else:
                lo = int(bounds[0])
                hi = int(bounds[1])
                if lo > hi:
                    lo, hi = hi, lo
                resolved_spec[actual_chain] = (lo, hi)
        resolved.append(dict(sorted(resolved_spec.items())))

    return resolved


def translate_axis_range_definitions_for_replication(
    defs: List[AxisSelection],
    original_chains: Iterable[str],
    final_chains: Iterable[str],
) -> List[AxisSelection]:
    """Translate --axis_range chain IDs into the final post-replication chain space.

    Replication renames the original base-copy chains to consecutive letters
    A, B, C, ... and may create later copies with subsequent letters. For
    axis-range definitions we do *not* duplicate user input across copies.
    Instead we interpret each chain token as follows:

      1) If it already matches a chain ID present after replication, keep it.
      2) Otherwise, if it matches one of the original input chain IDs, remap it
         once to that chain's renamed base-copy chain ID.
      3) Otherwise leave it unchanged so the normal validation path can raise a
         clear error.

    This makes it possible to target replicated chains explicitly while still
    accepting pre-replication base-chain IDs when they are not ambiguous in the
    final chain-ID space.
    """
    original_sorted = sorted(str(ch) for ch in original_chains)
    letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    if len(original_sorted) > len(letters):
        raise ValueError(
            f"Cannot remap {len(original_sorted)} original chains into the final A-Z chain-ID space."
        )

    mapping_upper_to_new = {
        ch.upper(): letters[i]
        for i, ch in enumerate(original_sorted)
    }
    final_lookup = {str(ch).upper(): str(ch) for ch in final_chains}

    translated: List[AxisSelection] = []
    for spec in defs:
        translated_spec: AxisSelection = {}
        for ch, bounds in spec.items():
            key = str(ch).upper()
            if key in final_lookup:
                actual_chain = final_lookup[key]
            elif key in mapping_upper_to_new:
                actual_chain = mapping_upper_to_new[key]
            else:
                actual_chain = str(ch)

            if actual_chain in translated_spec:
                raise ValueError(
                    f"Axis-range definition {_format_axis_range_spec(spec)} maps chain '{ch}' "
                    f"to '{actual_chain}' more than once during replication handling."
                )

            if bounds is None:
                translated_spec[actual_chain] = None
            else:
                lo = int(bounds[0])
                hi = int(bounds[1])
                if lo > hi:
                    lo, hi = hi, lo
                translated_spec[actual_chain] = (lo, hi)
        translated.append(dict(sorted(translated_spec.items())))

    return translated


def _format_axis_range_spec(spec: Dict[str, AxisBounds]) -> str:
    parts = []
    for ch in sorted(spec):
        bounds = spec[ch]
        if bounds is None:
            parts.append(f"{ch}")
        else:
            lo, hi = bounds
            parts.append(f"{ch}{lo}-{ch}{hi}")
    return ",".join(parts)


def _axis_range_total_span(spec: ResolvedAxisRange) -> int:
    return sum((hi - lo) for lo, hi in spec.values())


def validate_axis_range_definitions(
    defs: List[ResolvedAxisRange],
    chain_to_P_atoms: Dict[str, List[pdb_atom_record]],
    chain_to_helix: Dict[str, HelixID],
) -> None:
    for spec in defs:
        actual_helix: Optional[HelixID] = None
        for ch, (lo, hi) in spec.items():
            if ch not in chain_to_P_atoms:
                raise ValueError(
                    f"Axis-range definition {_format_axis_range_spec(spec)} refers to chain '{ch}', "
                    "but that chain has no P atoms."
                )
            if ch not in chain_to_helix:
                raise ValueError(
                    f"Axis-range definition {_format_axis_range_spec(spec)} refers to chain '{ch}', "
                    "but that chain is not assigned to any helix."
                )
            if actual_helix is None:
                actual_helix = chain_to_helix[ch]
            elif chain_to_helix[ch] != actual_helix:
                raise ValueError(
                    f"Axis-range definition {_format_axis_range_spec(spec)} spans multiple helices; "
                    "all chains in one axis-range definition must belong to the same helix group."
                )

            in_range = False
            for atom in chain_to_P_atoms.get(ch, []):
                if lo <= atom.resSeq <= hi:
                    in_range = True
                    break
            if not in_range:
                raise ValueError(
                    f"Axis-range definition {_format_axis_range_spec(spec)} selects no P atoms on chain '{ch}'."
                )


def find_axis_range_definition_for_positions(
    helix_id: HelixID,
    positions: List[Tuple[str, int]],
) -> Optional[ResolvedAxisRange]:
    if not positions or not _HELIX_AXIS_RANGE_DEFINITIONS:
        return None

    helix_set = set(helix_id)
    positions_in_helix = [(ch, res) for ch, res in positions if ch in helix_set]
    if not positions_in_helix:
        return None

    candidates: List[ResolvedAxisRange] = []
    for spec in _HELIX_AXIS_RANGE_DEFINITIONS:
        spec_chains = set(spec.keys())
        if not spec_chains.issubset(helix_set):
            continue

        ok = True
        for ch, res in positions_in_helix:
            if ch not in spec:
                ok = False
                break
            lo, hi = spec[ch]
            if res < lo or res > hi:
                ok = False
                break
        if ok:
            candidates.append(spec)

    if not candidates:
        return None

    ranked = sorted(
        candidates,
        key=lambda spec: (-len(spec), _axis_range_total_span(spec), _format_axis_range_spec(spec)),
    )
    best = ranked[0]
    if len(ranked) > 1:
        best_key = (-len(best), _axis_range_total_span(best), _format_axis_range_spec(best))
        second = ranked[1]
        second_key = (-len(second), _axis_range_total_span(second), _format_axis_range_spec(second))
        if best_key[:2] == second_key[:2] and best != second:
            raise ValueError(
                "Ambiguous axis-range definitions match helix "
                f"{helix_id_str(tuple(sorted(helix_id)))} for positions "
                + ", ".join(f"{res}{ch}" for ch, res in positions_in_helix)
                + ": "
                + _format_axis_range_spec(best)
                + " vs "
                + _format_axis_range_spec(second)
                + ". Please make the axis ranges non-overlapping or more specific."
            )

    return dict(best)


def _selection_chains(selection: Dict[str, AxisBounds]) -> Set[str]:
    return {str(ch) for ch in selection.keys()}


def _merge_move_selection_value(existing: AxisBounds, new: AxisBounds) -> AxisBounds:
    if existing is None or new is None:
        return None
    lo = min(int(existing[0]), int(new[0]))
    hi = max(int(existing[1]), int(new[1]))
    return (lo, hi)


def _find_explicit_helix_for_axis(
    axis_chains: Set[str],
    explicit_helix_defs: Optional[List[HelixID]],
) -> Optional[HelixID]:
    if not axis_chains or not explicit_helix_defs:
        return None
    for raw_h in explicit_helix_defs:
        helix_id = tuple(sorted(raw_h))
        if axis_chains.issubset(set(helix_id)):
            return helix_id
    return None


def build_axis_coupled_helix_defs(
    axis_defs: List[AxisSelection],
    move_defs: List[MoveSelection],
    explicit_helix_defs: Optional[List[HelixID]] = None,
) -> List[HelixID]:
    """Return helix groups implied by axis/move rows.

    These definitions are used before replication to avoid legacy triplex
    stdin prompts and to give replication a template for axis-coupled payloads.
    """
    groups: List[HelixID] = []
    seen: Set[HelixID] = set()
    for axis_def, move_def in pair_axis_move_definitions(axis_defs, move_defs):
        axis_chains = _selection_chains(axis_def)
        if not axis_chains:
            continue
        group_chains = set(axis_chains)
        explicit_group = _find_explicit_helix_for_axis(axis_chains, explicit_helix_defs)
        if explicit_group is not None:
            group_chains.update(explicit_group)
        group_chains.update(_selection_chains(move_def))
        if len(group_chains) < 1:
            continue
        group = tuple(sorted(group_chains))
        if group not in seen:
            groups.append(group)
            seen.add(group)
    return groups


def register_axis_coupling_overrides(
    axis_defs: List[ResolvedAxisRange],
    move_defs: List[MoveSelection],
    explicit_helix_defs: Optional[List[HelixID]] = None,
) -> None:
    """Register axis and movement overrides implied by axis/move rows."""
    axis_overrides = get_helix_axis_overrides()
    axis_range_overrides = get_helix_axis_range_overrides()
    move_selections = get_helix_move_selections()

    for axis_def, move_def in pair_axis_move_definitions(axis_defs, move_defs):
        axis_chains = _selection_chains(axis_def)
        if not axis_chains:
            continue
        group_chains = set(axis_chains)
        full_move_chains = set(axis_chains)
        move_chains = _selection_chains(move_def)
        explicit_group = _find_explicit_helix_for_axis(axis_chains, explicit_helix_defs)
        if explicit_group is not None:
            group_chains.update(explicit_group)
            full_move_chains.update(ch for ch in explicit_group if ch not in move_chains)
        group_chains.update(move_chains)
        group = tuple(sorted(group_chains))
        if not group:
            continue

        axis_overrides[group] = tuple(sorted(axis_chains))
        axis_range_overrides[group] = dict(sorted(axis_def.items()))

        selection: MoveSelection = {ch: None for ch in full_move_chains}
        for ch, bounds in move_def.items():
            if ch in selection:
                selection[ch] = _merge_move_selection_value(selection[ch], bounds)
            else:
                selection[ch] = bounds
        for ch in group:
            if ch not in selection:
                selection[ch] = None
        move_selections[group] = dict(sorted(selection.items()))

    set_helix_axis_overrides(axis_overrides)
    set_helix_axis_range_overrides(axis_range_overrides)
    set_helix_move_selections(move_selections)


def apply_axis_coupling_to_chain_map(
    chain_to_helix: Dict[str, HelixID],
    axis_defs: List[ResolvedAxisRange],
    move_defs: List[MoveSelection],
    chain_to_P_atoms: Dict[str, List[pdb_atom_record]],
) -> Dict[str, HelixID]:
    """Merge axis/move rows into a chain->helix map and register overrides."""
    if not axis_defs:
        return chain_to_helix

    updated = dict(chain_to_helix)
    axis_overrides = get_helix_axis_overrides()
    axis_range_overrides = get_helix_axis_range_overrides()
    move_selections = get_helix_move_selections()

    for axis_def, move_def in pair_axis_move_definitions(axis_defs, move_defs):
        axis_chains = _selection_chains(axis_def)
        if not axis_chains:
            continue

        group_chains: Set[str] = set(axis_chains)
        full_move_chains: Set[str] = set(axis_chains)
        move_chains = _selection_chains(move_def)
        for ch in axis_chains:
            if ch in updated:
                group_chains.update(updated[ch])
                full_move_chains.update(member for member in updated[ch] if member not in move_chains)
        group_chains.update(move_chains)
        group = tuple(sorted(group_chains))
        if not group:
            continue

        for ch, current_group in list(updated.items()):
            current_set = set(current_group)
            if not current_set.intersection(group_chains):
                continue
            if ch in group_chains:
                updated[ch] = group
            else:
                residual = tuple(sorted(current_set - group_chains))
                if residual:
                    updated[ch] = residual

        for ch in group:
            if ch in chain_to_P_atoms:
                updated[ch] = group

        axis_overrides[group] = tuple(sorted(axis_chains))
        axis_range_overrides[group] = dict(sorted(axis_def.items()))

        selection: MoveSelection = {ch: None for ch in full_move_chains}
        for ch, bounds in move_def.items():
            if ch in selection:
                selection[ch] = _merge_move_selection_value(selection[ch], bounds)
            else:
                selection[ch] = bounds
        for ch in group:
            if ch not in selection:
                selection[ch] = None
        move_selections[group] = dict(sorted(selection.items()))

    set_helix_axis_overrides(axis_overrides)
    set_helix_axis_range_overrides(axis_range_overrides)
    set_helix_move_selections(move_selections)
    return updated


def _register_triplex_axis_override(
    helix_id: HelixID,
    third_chain: str,
    context: str,
) -> HelixID:
    helix_key = _canonical_helix_id(helix_id)
    actual_third = _resolve_chain_id_for_helix(third_chain, helix_key)
    axis_chains = tuple(sorted(ch for ch in helix_key if ch != actual_third))
    if len(axis_chains) != 2:
        raise ValueError(
            f"Triplex axis definition for {_format_chain_group(helix_key)} must leave exactly two duplex chains."
        )

    _HELIX_AXIS_OVERRIDES[helix_key] = axis_chains
    sys.stderr.write(
        f"[re_helix] Triplex helix {_format_chain_group(helix_key)} ({context}): "
        f"third strand = '{actual_third}', axis strands = {_format_chain_group(axis_chains)}.\n"
    )
    return helix_key


def _prompt_for_triplex_axis_override(helix_id: HelixID, context: str) -> HelixID:
    helix_key = _canonical_helix_id(helix_id)
    if len(helix_key) != 3:
        return helix_key
    if helix_key in _HELIX_AXIS_OVERRIDES:
        return helix_key

    sys.stderr.write(
        f"[re_helix] Potential triplex detected for {_format_chain_group(helix_key)} ({context}).\n"
    )
    sys.stderr.write(
        "[re_helix] Please provide the chainID for the third strand. "
        "The helical axis will be determined from the other two strands.\n"
    )

    choices = ", ".join(helix_key)
    prompt = f"Third-strand chainID for {_format_chain_group(helix_key)} [{choices}]: "

    while True:
        try:
            answer = input(prompt)
        except EOFError as exc:
            raise ValueError(
                f"Potential triplex detected for {_format_chain_group(helix_key)}, but no interactive input was available to identify the third strand."
            ) from exc

        try:
            return _register_triplex_axis_override(helix_key, answer, context)
        except ValueError as exc:
            sys.stderr.write(f"[re_helix] {exc}\n")
            sys.stderr.write(
                f"[re_helix] Please enter one of the chain IDs: {choices}.\n"
            )


def _replicate_axis_overrides(
    axis_overrides_old: Dict[HelixID, HelixID],
    mapping_upper_to_new: Dict[str, str],
    base_index: Dict[str, int],
    num_copies: int,
    letters: str,
) -> Dict[HelixID, HelixID]:
    if not axis_overrides_old:
        return {}

    remapped_base: Dict[HelixID, HelixID] = {}
    for helix_id, axis_chains in axis_overrides_old.items():
        try:
            base_helix = tuple(sorted(mapping_upper_to_new[ch.upper()] for ch in helix_id))
            base_axis = tuple(sorted(mapping_upper_to_new[ch.upper()] for ch in axis_chains))
        except KeyError:
            continue
        remapped_base[base_helix] = base_axis

    axis_overrides_repl: Dict[HelixID, HelixID] = {}
    n_base = len(base_index)
    for copy_idx in range(num_copies):
        offset = copy_idx * n_base
        for base_helix, base_axis in remapped_base.items():
            new_helix = tuple(sorted(letters[offset + base_index[ch]] for ch in base_helix))
            new_axis = tuple(sorted(letters[offset + base_index[ch]] for ch in base_axis))
            axis_overrides_repl[new_helix] = new_axis

    return axis_overrides_repl


def _replicate_axis_range_overrides(
    axis_range_overrides_old: Dict[HelixID, ResolvedAxisRange],
    mapping_upper_to_new: Dict[str, str],
    base_index: Dict[str, int],
    num_copies: int,
    letters: str,
) -> Dict[HelixID, ResolvedAxisRange]:
    if not axis_range_overrides_old:
        return {}

    remapped_base: Dict[HelixID, ResolvedAxisRange] = {}
    for helix_id, spec in axis_range_overrides_old.items():
        try:
            base_helix = tuple(sorted(mapping_upper_to_new[ch.upper()] for ch in helix_id))
            base_spec = {
                mapping_upper_to_new[ch.upper()]: (int(bounds[0]), int(bounds[1]))
                for ch, bounds in spec.items()
            }
        except KeyError:
            continue
        remapped_base[base_helix] = dict(sorted(base_spec.items()))

    overrides_repl: Dict[HelixID, ResolvedAxisRange] = {}
    n_base = len(base_index)
    for copy_idx in range(num_copies):
        offset = copy_idx * n_base
        for base_helix, base_spec in remapped_base.items():
            new_helix = tuple(sorted(letters[offset + base_index[ch]] for ch in base_helix))
            new_spec = {
                letters[offset + base_index[ch]]: bounds
                for ch, bounds in base_spec.items()
            }
            overrides_repl[new_helix] = dict(sorted(new_spec.items()))

    return overrides_repl


def _replicate_move_selections(
    move_selections_old: Dict[HelixID, MoveSelection],
    mapping_upper_to_new: Dict[str, str],
    base_index: Dict[str, int],
    num_copies: int,
    letters: str,
) -> Dict[HelixID, MoveSelection]:
    if not move_selections_old:
        return {}

    remapped_base: Dict[HelixID, MoveSelection] = {}
    for helix_id, selection in move_selections_old.items():
        try:
            base_helix = tuple(sorted(mapping_upper_to_new[ch.upper()] for ch in helix_id))
            base_selection: MoveSelection = {}
            for ch, bounds in selection.items():
                mapped_ch = mapping_upper_to_new[ch.upper()]
                base_selection[mapped_ch] = None if bounds is None else (int(bounds[0]), int(bounds[1]))
        except KeyError:
            continue
        remapped_base[base_helix] = dict(sorted(base_selection.items()))

    selections_repl: Dict[HelixID, MoveSelection] = {}
    n_base = len(base_index)
    for copy_idx in range(num_copies):
        offset = copy_idx * n_base
        for base_helix, base_selection in remapped_base.items():
            new_helix = tuple(sorted(letters[offset + base_index[ch]] for ch in base_helix))
            new_selection: MoveSelection = {}
            for ch, bounds in base_selection.items():
                new_ch = letters[offset + base_index[ch]]
                new_selection[new_ch] = None if bounds is None else bounds
            selections_repl[new_helix] = dict(sorted(new_selection.items()))

    return selections_repl


# ---------------------------------------------------------------------------
# Small vector helpers
# ---------------------------------------------------------------------------

def v_add(a: Tuple[float, float, float],
          b: Tuple[float, float, float]) -> Tuple[float, float, float]:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def v_sub(a: Tuple[float, float, float],
          b: Tuple[float, float, float]) -> Tuple[float, float, float]:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def v_scale(a: Tuple[float, float, float],
            s: float) -> Tuple[float, float, float]:
    return (a[0] * s, a[1] * s, a[2] * s)


def v_dot(a: Tuple[float, float, float],
          b: Tuple[float, float, float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def v_cross(a: Tuple[float, float, float],
            b: Tuple[float, float, float]) -> Tuple[float, float, float]:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def v_length(a: Tuple[float, float, float]) -> float:
    return math.sqrt(v_dot(a, a))


def v_norm(a: Tuple[float, float, float]) -> Tuple[float, float, float]:
    l = v_length(a)
    if l == 0.0:
        return (0.0, 0.0, 0.0)
    return (a[0] / l, a[1] / l, a[2] / l)


def normalize_user_axis_definition(
    axis_dir: Iterable[float],
    axis_point: Iterable[float],
) -> UserAxisDefinition:
    dir_tuple = tuple(float(value) for value in axis_dir)
    point_tuple = tuple(float(value) for value in axis_point)
    if len(dir_tuple) != 3 or len(point_tuple) != 3:
        raise ValueError("User-defined axis direction and point must each have exactly three numeric values.")
    if not all(math.isfinite(value) for value in dir_tuple + point_tuple):
        raise ValueError("User-defined axis direction and point values must be finite numbers.")
    normalized_dir = v_norm(dir_tuple)  # type: ignore[arg-type]
    if v_length(normalized_dir) < 1.0e-12:
        raise ValueError("User-defined axis direction vector cannot be zero.")
    return normalized_dir, point_tuple  # type: ignore[return-value]


def rotate_around_line(
    point: Tuple[float, float, float],
    axis_point: Tuple[float, float, float],
    axis_dir: Tuple[float, float, float],
    angle: float,
) -> Tuple[float, float, float]:
    """
    Rotate a 3D point around a line (axis_point + t * axis_dir) by 'angle'
    using Rodrigues' rotation formula. axis_dir need not be normalized.
    """
    ux, uy, uz = v_norm(axis_dir)
    px, py, pz = point
    ax, ay, az = axis_point
    rx, ry, rz = px - ax, py - ay, pz - az

    cos_t = math.cos(angle)
    sin_t = math.sin(angle)

    dot_ur = ux * rx + uy * ry + uz * rz
    cx = uy * rz - uz * ry
    cy = uz * rx - ux * rz
    cz = ux * ry - uy * rx

    rx_rot = rx * cos_t + cx * sin_t + ux * dot_ur * (1.0 - cos_t)
    ry_rot = ry * cos_t + cy * sin_t + uy * dot_ur * (1.0 - cos_t)
    rz_rot = rz * cos_t + cz * sin_t + uz * dot_ur * (1.0 - cos_t)

    return (ax + rx_rot, ay + ry_rot, az + rz_rot)


def normalize_angle(angle: float) -> float:
    """
    Wrap angle into [-pi, pi).
    """
    two_pi = 2.0 * math.pi
    return ((angle + math.pi) % two_pi) - math.pi


def perpendicular_unit_vector(axis_dir: Tuple[float, float, float]) -> Tuple[float, float, float]:
    """Return a deterministic unit vector perpendicular to axis_dir.

    Used for the optional 180-degree pre-flip of the moving helix when
    exploring the parallel vs. anti-parallel axis orientations under
    --axis_parallel y (V2.1 update).
    """
    u = v_norm(axis_dir)
    if v_length(u) < 1.0e-12:
        return (0.0, 0.0, 0.0)
    # Choose a base vector that is not close to colinear with u.
    if abs(u[0]) < 0.9:
        base = (1.0, 0.0, 0.0)
    else:
        base = (0.0, 1.0, 0.0)
    perp = v_cross(u, base)
    if v_length(perp) < 1.0e-8:
        # Very unlikely fallback, but keep it robust.
        base = (0.0, 0.0, 1.0)
        perp = v_cross(u, base)
    return v_norm(perp)


# ---------------------------------------------------------------------------
# PDB maps: chains, residues, P atoms
# ---------------------------------------------------------------------------

def build_nucleic_acid_maps(
    rec_list: List[pdb_atom_record],
):
    """
    Build mappings:
        chain_to_res_atoms[chain][resSeq] -> [atoms...]
        residue_to_P_atom[(chain, resSeq)] -> P atom
        chain_to_P_atoms[chain] -> list of P atoms
    """
    chain_to_res_atoms: Dict[str, Dict[int, List[pdb_atom_record]]] = defaultdict(
        lambda: defaultdict(list)
    )
    residue_to_P_atom: Dict[Tuple[str, int], pdb_atom_record] = {}
    chain_to_P_atoms: Dict[str, List[pdb_atom_record]] = defaultdict(list)

    for atom in rec_list:
        if atom.recordName not in ("ATOM", "HETATM"):
            continue
        chain = atom.chainID
        res = atom.resSeq
        chain_to_res_atoms[chain][res].append(atom)

        if atom.name.strip() == "P":
            key = (chain, res)
            if key not in residue_to_P_atom:
                residue_to_P_atom[key] = atom
            chain_to_P_atoms[chain].append(atom)

    return chain_to_res_atoms, residue_to_P_atom, chain_to_P_atoms


def median_min_distance(
    atoms_a: Iterable[pdb_atom_record],
    atoms_b: Iterable[pdb_atom_record],
) -> float:
    """
    For each atom in atoms_a, compute minimum distance to any atom in atoms_b,
    collect distances, return median.
    """
    atoms_a = list(atoms_a)
    atoms_b = list(atoms_b)
    if not atoms_a or not atoms_b:
        return float("inf")

    dists: List[float] = []
    for a in atoms_a:
        min_sq = float("inf")
        for b in atoms_b:
            dx = a.x - b.x
            dy = a.y - b.y
            dz = a.z - b.z
            dsq = dx * dx + dy * dy + dz * dz
            if dsq < min_sq:
                min_sq = dsq
        dists.append(math.sqrt(min_sq))

    dists.sort()
    return dists[len(dists) // 2]


def compute_chain_partner_map(
    chain_to_P_atoms: Dict[str, List[pdb_atom_record]],
) -> Dict[str, HelixID]:
    """
    Auto-detect helical partners.

    Standard behavior:
        pair chains whose P atoms are closest on average, returning
        chain_to_helix[chain] -> helix_id for automatically detected duplexes.

    Triplex behavior:
        if the input contains exactly three chains with P atoms, treat this as
        a potential triplex helix, prompt for the third-strand chain ID, and
        return one three-chain helix group. The helical axis for that group is
        then defined from the other two chains.
    """
    chains = sorted(chain_to_P_atoms.keys())
    if len(chains) < 2:
        raise ValueError("Need at least two chains with P atoms to define helices.")

    if len(chains) == 3:
        helix_id = _prompt_for_triplex_axis_override(
            tuple(chains),
            "automatic 3-chain input",
        )
        return {ch: helix_id for ch in chains}

    pair_dist: Dict[Tuple[str, str], float] = {}
    for i, ch1 in enumerate(chains):
        for j in range(i + 1, len(chains)):
            ch2 = chains[j]
            d = median_min_distance(chain_to_P_atoms[ch1], chain_to_P_atoms[ch2])
            pair_dist[(ch1, ch2)] = d

    def get_pair_dist(a: str, b: str) -> float:
        if a == b:
            return 0.0
        key = (a, b)
        if key in pair_dist:
            return pair_dist[key]
        key = (b, a)
        return pair_dist.get(key, float("inf"))

    chain_to_partner: Dict[str, str] = {}
    for ch in chains:
        best = None
        best_d = float("inf")
        for other in chains:
            if other == ch:
                continue
            d = get_pair_dist(ch, other)
            if d < best_d:
                best_d = d
                best = other
        if best is not None and best_d < float("inf"):
            chain_to_partner[ch] = best

    chain_to_helix: Dict[str, HelixID] = {}
    for ch, partner in chain_to_partner.items():
        helix_id = tuple(sorted((ch, partner)))
        chain_to_helix[ch] = helix_id

    return chain_to_helix


def build_chain_to_helix_from_defs(
    helix_defs: List[HelixID],
    chain_to_P_atoms: Dict[str, List[pdb_atom_record]],
) -> Dict[str, HelixID]:
    """
    Build chain->helix map from user-defined helix tokens like (AB), (CD),
    (ABMN), etc. All chains inside one parentheses group are treated as a
    single rigid helix group.

    If a helix group contains exactly three chains, treat it as a potential
    triplex helix, ask the user which chain is the third strand, and use the
    other two chains for helical-axis estimation.
    """
    chain_to_helix: Dict[str, HelixID] = {}

    for raw_h in helix_defs:
        # Sort to canonicalize the group (so (ABMN) and (MNBA) are treated identically)
        helix_id: HelixID = tuple(sorted(raw_h))
        if len(helix_id) == 3:
            _prompt_for_triplex_axis_override(helix_id, "user-defined helix")

        mapped_any_p_chain = False
        for ch in helix_id:
            if ch not in chain_to_P_atoms:
                continue
            if ch in chain_to_helix and chain_to_helix[ch] != helix_id:
                raise ValueError(
                    f"Chain '{ch}' appears in more than one helix definition: "
                    f"{helix_id_str(chain_to_helix[ch])} and {helix_id_str(helix_id)}."
                )
            chain_to_helix[ch] = helix_id
            mapped_any_p_chain = True

        if not mapped_any_p_chain:
            raise ValueError(
                f"User-defined helix {helix_id_str(helix_id)} contains no chains with P atoms; "
                "at least one P-bearing chain is needed to define or participate in a helix."
            )

    return chain_to_helix


def helix_id_str(h: HelixID) -> str:
    # e.g. ('A','B','M','N') -> "(ABMN)"
    return "(" + "".join(h) + ")"


# ---------------------------------------------------------------------------
# Helix axis estimation
# ---------------------------------------------------------------------------

def compute_helix_axis(
    chain_to_P_atoms: Dict[str, List[pdb_atom_record]],
    helix_id: HelixID,
    residue_ranges: Optional[ResolvedAxisRange] = None,
) -> Tuple[Tuple[float, float, float], Tuple[float, float, float]]:
    """
    Estimate the helical axis of a helix group.

    Normally all chains in helix_id contribute P atoms to the PCA-like axis
    estimate. For triplex helices, a registered axis-chain override is used so
    that only the duplex-forming strands define the axis. If residue_ranges is
    provided, only P atoms inside those chain-specific residue windows are used.

    Returns (axis_dir (unit), center point).
    """
    helix_key = _canonical_helix_id(helix_id)
    if residue_ranges:
        axis_chains = tuple(sorted(residue_ranges.keys()))
    elif helix_key in _HELIX_AXIS_RANGE_OVERRIDES:
        residue_ranges = _HELIX_AXIS_RANGE_OVERRIDES[helix_key]
        axis_chains = tuple(sorted(residue_ranges.keys()))
    else:
        axis_chains = get_axis_chains_for_helix(helix_key)

    coords: List[Tuple[float, float, float]] = []
    for ch in axis_chains:
        for atom in chain_to_P_atoms.get(ch, []):
            if residue_ranges and ch in residue_ranges:
                lo, hi = residue_ranges[ch]
                if atom.resSeq < lo or atom.resSeq > hi:
                    continue
            coords.append((atom.x, atom.y, atom.z))

    if len(coords) < 2:
        raise ValueError(
            f"Not enough P atoms to define axis for helix {helix_id_str(tuple(sorted(helix_id)))}."
        )

    n = len(coords)
    cx = sum(x for x, y, z in coords) / n
    cy = sum(y for x, y, z in coords) / n
    cz = sum(z for x, y, z in coords) / n

    cov = [[0.0, 0.0, 0.0] for _ in range(3)]
    for x, y, z in coords:
        dx = x - cx
        dy = y - cy
        dz = z - cz
        cov[0][0] += dx * dx
        cov[0][1] += dx * dy
        cov[0][2] += dx * dz
        cov[1][0] += dy * dx
        cov[1][1] += dy * dy
        cov[1][2] += dy * dz
        cov[2][0] += dz * dx
        cov[2][1] += dz * dy
        cov[2][2] += dz * dz

    for i in range(3):
        for j in range(3):
            cov[i][j] /= n

    v = v_norm((1.0, 1.0, 1.0))
    for _ in range(30):
        vx = cov[0][0] * v[0] + cov[0][1] * v[1] + cov[0][2] * v[2]
        vy = cov[1][0] * v[0] + cov[1][1] * v[1] + cov[1][2] * v[2]
        vz = cov[2][0] * v[0] + cov[2][1] * v[1] + cov[2][2] * v[2]
        norm = math.sqrt(vx * vx + vy * vy + vz * vz)
        if norm < 1.0e-12:
            break
        v = (vx / norm, vy / norm, vz / norm)

    axis_dir = v_norm(v)
    center = (cx, cy, cz)
    return axis_dir, center


def align_axes_for_pair(
    rec_list: List[pdb_atom_record],
    helix_fixed: HelixID,
    helix_moving: HelixID,
    chain_to_P_atoms: Dict[str, List[pdb_atom_record]],
    axis_dist: float,
    subset_fixed: Optional[Iterable[str]] = None,
    subset_moving: Optional[Iterable[str]] = None,
    axis_ranges_fixed: Optional[ResolvedAxisRange] = None,
    axis_ranges_moving: Optional[ResolvedAxisRange] = None,
) -> Tuple[Tuple[float, float, float], Tuple[float, float, float], Tuple[float, float, float]]:
    """
    Rigidly align the axis of helix_moving to that of helix_fixed, and
    translate helix_moving so that the two axes are separated by axis_dist.

    By default, the helical axes are estimated using ALL chains in the
    helix groups. If subset_fixed/subset_moving are provided, the axes
    are estimated using only those chains (intersected with the helix
    group), while the rigid transform is still applied to all chains in
    helix_moving.

    Returns:
        (axis_dir, center1, center2) after alignment, where center1/center2
        are the centers of the (possibly subset-based) axes.
    """
    # --- NEW: choose which chains define the axes (subset of the group) ---
    # For registered triplex helices, always use the duplex-forming chains for
    # axis estimation, regardless of which chains appear in the P-pair subset.
    fixed_axis_override = get_helix_axis_overrides().get(tuple(sorted(helix_fixed)))
    moving_axis_override = get_helix_axis_overrides().get(tuple(sorted(helix_moving)))
    axis_range_overrides = get_helix_axis_range_overrides()
    if axis_ranges_fixed is None:
        axis_ranges_fixed = axis_range_overrides.get(tuple(sorted(helix_fixed)))
    if axis_ranges_moving is None:
        axis_ranges_moving = axis_range_overrides.get(tuple(sorted(helix_moving)))

    if axis_ranges_fixed:
        fixed_axis_chains = tuple(sorted(axis_ranges_fixed.keys()))
    elif fixed_axis_override:
        fixed_axis_chains = fixed_axis_override
    elif subset_fixed:
        sf = set(subset_fixed)
        fixed_axis_chains = tuple(sorted(ch for ch in helix_fixed if ch in sf))
        if not fixed_axis_chains:
            fixed_axis_chains = helix_fixed
    else:
        fixed_axis_chains = helix_fixed

    if axis_ranges_moving:
        moving_axis_chains = tuple(sorted(axis_ranges_moving.keys()))
    elif moving_axis_override:
        moving_axis_chains = moving_axis_override
    elif subset_moving:
        sm = set(subset_moving)
        moving_axis_chains = tuple(sorted(ch for ch in helix_moving if ch in sm))
        if not moving_axis_chains:
            moving_axis_chains = helix_moving
    else:
        moving_axis_chains = helix_moving

    axis1_dir, center1 = compute_helix_axis(
        chain_to_P_atoms,
        fixed_axis_chains,
        residue_ranges=axis_ranges_fixed,
    )
    axis2_dir, center2 = compute_helix_axis(
        chain_to_P_atoms,
        moving_axis_chains,
        residue_ranges=axis_ranges_moving,
    )

    axis1_dir = v_norm(axis1_dir)
    axis2_dir = v_norm(axis2_dir)

    # Make them parallel, not anti-parallel
    if v_dot(axis1_dir, axis2_dir) < 0.0:
        axis2_dir = v_scale(axis2_dir, -1.0)

    cross_v = v_cross(axis2_dir, axis1_dir)
    norm_cross = v_length(cross_v)
    if norm_cross > 1.0e-6 and abs(v_dot(axis1_dir, axis2_dir)) < 1.0 - 1.0e-6:
        rot_axis = v_scale(cross_v, 1.0 / norm_cross)
        angle = math.acos(max(min(v_dot(axis2_dir, axis1_dir), 1.0), -1.0))

        for atom in rec_list:
            if atom.recordName not in ("ATOM", "HETATM"):
                continue
            if not atom_moves_with_helix(atom, helix_moving):
                continue
            x, y, z = rotate_around_line(
                (atom.x, atom.y, atom.z), center2, rot_axis, angle
            )
            atom.update_xyz(x, y, z)

        axis2_dir = axis1_dir

    # Re-estimate centers (using the same axis-defining chain subsets)
    axis1_dir, center1 = compute_helix_axis(
        chain_to_P_atoms,
        fixed_axis_chains,
        residue_ranges=axis_ranges_fixed,
    )
    axis2_dir, center2 = compute_helix_axis(
        chain_to_P_atoms,
        moving_axis_chains,
        residue_ranges=axis_ranges_moving,
    )
    axis1_dir = v_norm(axis1_dir)

    # Enforce axis_dist
    delta = v_sub(center2, center1)
    parallel = v_dot(delta, axis1_dir)
    perp_vec = v_sub(delta, v_scale(axis1_dir, parallel))
    perp_len = v_length(perp_vec)

    if perp_len < 1.0e-6:
        if abs(axis1_dir[0]) < 0.9:
            ref = (1.0, 0.0, 0.0)
        else:
            ref = (0.0, 1.0, 0.0)
        perp_unit = v_norm(v_cross(axis1_dir, ref))
    else:
        perp_unit = v_scale(perp_vec, 1.0 / perp_len)

    desired_delta = v_add(v_scale(axis1_dir, parallel), v_scale(perp_unit, axis_dist))
    translation = v_sub(v_add(center1, desired_delta), center2)

    # Apply the translation to ALL chains in the moving helix group
    for atom in rec_list:
        if atom.recordName not in ("ATOM", "HETATM"):
            continue
        if not atom_moves_with_helix(atom, helix_moving):
            continue
        new_xyz = v_add((atom.x, atom.y, atom.z), translation)
        atom.update_xyz(*new_xyz)

    # Final centers (still subset-based)
    _, center1 = compute_helix_axis(
        chain_to_P_atoms,
        fixed_axis_chains,
        residue_ranges=axis_ranges_fixed,
    )
    _, center2 = compute_helix_axis(
        chain_to_P_atoms,
        moving_axis_chains,
        residue_ranges=axis_ranges_moving,
    )

    return axis1_dir, center1, center2


# ---------------------------------------------------------------------------
# Objective function / optimisation for P-pair distances
# ---------------------------------------------------------------------------

def build_pair_objective(
    pairs: List[Tuple[Tuple[str, int], Tuple[str, int]]],
    residue_to_P_atom: Dict[Tuple[str, int], pdb_atom_record],
    axis_dir: Tuple[float, float, float],
    axis1_point: Tuple[float, float, float],
    axis2_point: Tuple[float, float, float],
    axis_parallel: bool,
    pre_flip: bool = False,
    fixed_beta: Optional[float] = None,
    fixed_rho: Optional[float] = None,
):
    """
    Build an objective f(params) that returns sum of squared distances between
    specified P-pairs, after transforming helix 2.

    Transform for helix 2 (in this order):

        1) rotate around helix 1's axis by phi,
        2) translate along axis_dir by d,
        3) rotate around helix 2's axis by tau,
        4) if axis_parallel is False, rotate around a P-pair optimised
           common perpendicular line (through anchor1) by beta.

    If axis_parallel is True:
        params = [d, tau, phi]
    Else if fixed_beta is None:
        params = [d, tau, phi, beta]
    Else:
        params = [d, tau, phi] and beta is held at fixed_beta.

    fixed_beta, when not None, is a user-defined beta angle in radians and is
    meaningful only when axis_parallel is False.

    Also returns d0, an initial guess for d, and anchor1 (the P-pair-based
    anchor point on helix 1's axis used for beta).

    Note: the 'pairs' passed in are already the subset chosen for alignment
    (either 1–2 P pairs or the two extreme P pairs along helix 1's axis).
    """
    if fixed_rho is not None:
        if fixed_beta is not None:
            raise ValueError("Use only one of fixed_beta or legacy fixed_rho.")
        note_legacy_rho_alias_once()
        fixed_beta = fixed_rho

    p1_list: List[Tuple[float, float, float]] = []
    p2_list: List[Tuple[float, float, float]] = []

    for posA, posB in pairs:
        atom1 = residue_to_P_atom.get(posA)
        atom2 = residue_to_P_atom.get(posB)
        if atom1 is None:
            raise ValueError(
                f"No P atom found for residue {posA[1]}{posA[0]} (helix 1)."
            )
        if atom2 is None:
            raise ValueError(
                f"No P atom found for residue {posB[1]}{posB[0]} (helix 2)."
            )
        p1_list.append((atom1.x, atom1.y, atom1.z))
        p2_list.append((atom2.x, atom2.y, atom2.z))

    axis_dir = v_norm(axis_dir)

    # Optional 180-degree pre-flip of helix-2 about a perpendicular axis through axis2_point.
    # This is used to explore parallel vs. anti-parallel axis orientations (V2.1 update).
    flip_axis_dir = perpendicular_unit_vector(axis_dir) if pre_flip else (0.0, 0.0, 0.0)

    # Choose anchor1: projection of midpoint/centroid of helix-1 P atoms onto helix-1 axis.
    # For exactly 2 pairs this is their midpoint; for 1 it's that point itself.
    if p1_list:
        cx = sum(p[0] for p in p1_list) / len(p1_list)
        cy = sum(p[1] for p in p1_list) / len(p1_list)
        cz = sum(p[2] for p in p1_list) / len(p1_list)
        P1_mean = (cx, cy, cz)

        delta_c1 = v_sub(P1_mean, axis1_point)
        s_bar = v_dot(delta_c1, axis_dir)
        anchor1 = v_add(axis1_point, v_scale(axis_dir, s_bar))
    else:
        anchor1 = axis1_point

    def objective(params: List[float]) -> float:
        if axis_parallel:
            if len(params) != 3:
                raise ValueError("Expected 3 params [d, tau, phi] for axis_parallel=True.")
            d, tau, phi = params
            beta = 0.0
        else:
            if fixed_beta is None:
                if len(params) != 4:
                    raise ValueError("Expected 4 params [d, tau, phi, beta] for axis_parallel=False.")
                d, tau, phi, beta = params
            else:
                if len(params) != 3:
                    raise ValueError("Expected 3 params [d, tau, phi] when fixed_beta is supplied.")
                d, tau, phi = params
                beta = fixed_beta

        # Center of helix-2 axis after phi and d
        C2_phi = rotate_around_line(axis2_point, axis1_point, axis_dir, phi)
        C2_phi_d = v_add(C2_phi, v_scale(axis_dir, d))
        axis2_base = C2_phi_d  # axis point for tau

        if axis_parallel:
            axis_perp_dir = None
            axis_perp_point = None
        else:
            # Point on helix-2 axis closest to anchor1
            delta_a2 = v_sub(anchor1, C2_phi_d)
            t = v_dot(delta_a2, axis_dir)
            C2_base = v_add(C2_phi_d, v_scale(axis_dir, t))
            r_perp = v_sub(C2_base, anchor1)
            r_len = v_length(r_perp)
            if r_len < 1.0e-6:
                # fallback direction ⟂ axis_dir
                if abs(axis_dir[0]) < 0.9:
                    base_vec = (1.0, 0.0, 0.0)
                else:
                    base_vec = (0.0, 1.0, 0.0)
                axis_perp_dir = v_norm(v_cross(axis_dir, base_vec))
            else:
                axis_perp_dir = v_scale(r_perp, 1.0 / r_len)
            axis_perp_point = anchor1

        total = 0.0
        for p1, p2_0 in zip(p1_list, p2_list):
            p = p2_0

            # Optional 180-degree pre-flip about a perpendicular axis through axis2_point
            # (explores parallel vs. anti-parallel axis orientations).
            if pre_flip:
                p = rotate_around_line(p, axis2_point, flip_axis_dir, math.pi)

            # phi: rotate around helix-1 axis
            if abs(phi) > 1.0e-10:
                p = rotate_around_line(p, axis1_point, axis_dir, phi)

            # d: translate along axis_dir
            if abs(d) > 1.0e-10:
                p = v_add(p, v_scale(axis_dir, d))

            # tau: rotate around helix-2 axis
            if abs(tau) > 1.0e-10:
                p = rotate_around_line(p, axis2_base, axis_dir, tau)

            # beta: tilt around common perpendicular through anchor1
            if (not axis_parallel) and abs(beta) > 1.0e-10:
                p = rotate_around_line(p, axis_perp_point, axis_perp_dir, beta)  # type: ignore[arg-type]

            dx = p1[0] - p[0]
            dy = p1[1] - p[1]
            dz = p1[2] - p[2]
            total += dx * dx + dy * dy + dz * dz

        return total

    # initial guess for d: align projections on axis_dir
    if not p1_list:
        d0 = 0.0
    else:
        diffs: List[float] = []
        for p1, p2 in zip(p1_list, p2_list):
            p2_eff = p2
            if pre_flip:
                p2_eff = rotate_around_line(p2_eff, axis2_point, flip_axis_dir, math.pi)
            diffs.append(v_dot(v_sub(p1, p2_eff), axis_dir))
        d0 = sum(diffs) / len(diffs)

    return objective, p1_list, p2_list, d0, anchor1


def build_user_axis_transform_objective(
    pairs: List[Tuple[Tuple[str, int], Tuple[str, int]]],
    residue_to_P_atom: Dict[Tuple[str, int], pdb_atom_record],
    axis_dir: Point3D,
    axis_point: Point3D,
):
    """
    Build an objective for user-defined-axis mode.

    The movable helix is rotated around the supplied line and then translated
    freely in XYZ. No fixed/moving helix axes are estimated in this mode.
    """
    p1_list: List[Point3D] = []
    p2_list: List[Point3D] = []

    for posA, posB in pairs:
        atom1 = residue_to_P_atom.get(posA)
        atom2 = residue_to_P_atom.get(posB)
        if atom1 is None:
            raise ValueError(
                f"No P atom found for residue {posA[1]}{posA[0]} (fixed helix)."
            )
        if atom2 is None:
            raise ValueError(
                f"No P atom found for residue {posB[1]}{posB[0]} (moving helix)."
            )
        p1_list.append((atom1.x, atom1.y, atom1.z))
        p2_list.append((atom2.x, atom2.y, atom2.z))

    axis_dir = v_norm(axis_dir)

    def translation_for_angle(angle: float) -> Point3D:
        if not p1_list:
            return (0.0, 0.0, 0.0)
        rotated_p2 = [
            rotate_around_line(p2, axis_point, axis_dir, angle)
            for p2 in p2_list
        ]
        return (
            sum(p1[0] - p2[0] for p1, p2 in zip(p1_list, rotated_p2)) / len(p1_list),
            sum(p1[1] - p2[1] for p1, p2 in zip(p1_list, rotated_p2)) / len(p1_list),
            sum(p1[2] - p2[2] for p1, p2 in zip(p1_list, rotated_p2)) / len(p1_list),
        )

    def objective(params: List[float]) -> float:
        if len(params) != 1:
            raise ValueError("Expected 1 param [angle] for user-defined-axis alignment.")
        angle = params[0]
        translation = translation_for_angle(angle)
        total = 0.0
        for p1, p2_0 in zip(p1_list, p2_list):
            p = rotate_around_line(p2_0, axis_point, axis_dir, angle)
            p = v_add(p, translation)
            dx = p1[0] - p[0]
            dy = p1[1] - p[1]
            dz = p1[2] - p[2]
            total += dx * dx + dy * dy + dz * dz
        return total

    return objective, p1_list, p2_list, translation_for_angle


def coordinate_descent(
    objective,
    initial_params: List[float],
    initial_steps: List[float],
    angle_indices: Set[int],
    max_iter: int = 80,
    shrink: float = 0.5,
    tol: float = 1.0e-3,
) -> Tuple[List[float], float]:
    """
    Simple coordinate-wise hill-climbing with shrinking step sizes.
    """
    x = list(initial_params)
    steps = list(initial_steps)
    best = objective(x)
    n = len(x)

    for _ in range(max_iter):
        improved = False
        for i in range(n):
            for sign in (+1.0, -1.0):
                trial = x.copy()
                trial[i] += sign * steps[i]
                if i in angle_indices:
                    trial[i] = normalize_angle(trial[i])
                val = objective(trial)
                if val < best - 1.0e-8:
                    x = trial
                    best = val
                    improved = True
                    break
        if not improved:
            if max(steps) < tol:
                break
            steps = [s * shrink for s in steps]

    return x, best


def apply_transform_to_helix(
    rec_list: List[pdb_atom_record],
    helix_moving: HelixID,
    axis_dir: Tuple[float, float, float],
    axis1_point: Tuple[float, float, float],
    axis2_point: Tuple[float, float, float],
    d: float,
    tau: float,
    phi: float,
    beta: float,
    axis_parallel: bool,
    anchor1: Tuple[float, float, float],
    pre_flip: bool = False,
) -> None:
    """
    Apply the optimised (d, tau, phi, beta) transform to all atoms of the
    moving helix group in-place.

    Order:
        1) rotate around helix 1's axis by phi,
        2) translate along axis_dir by d,
        3) rotate around helix 2's axis by tau,
        4) if axis_parallel is False, rotate around the P-pair-based common
           perpendicular through anchor1 by beta.
    """
    axis_dir = v_norm(axis_dir)

    # Optional 180-degree pre-flip (parallel vs. anti-parallel exploration).
    flip_axis_dir = perpendicular_unit_vector(axis_dir) if pre_flip else (0.0, 0.0, 0.0)

    # helix-2 axis center after phi and d
    C2_phi = rotate_around_line(axis2_point, axis1_point, axis_dir, phi)
    C2_phi_d = v_add(C2_phi, v_scale(axis_dir, d))
    axis2_base = C2_phi_d

    if axis_parallel:
        axis_perp_dir = None
        axis_perp_point = None
    else:
        delta_a2 = v_sub(anchor1, C2_phi_d)
        t = v_dot(delta_a2, axis_dir)
        C2_base = v_add(C2_phi_d, v_scale(axis_dir, t))
        r_perp = v_sub(C2_base, anchor1)
        r_len = v_length(r_perp)
        if r_len < 1.0e-6:
            if abs(axis_dir[0]) < 0.9:
                base_vec = (1.0, 0.0, 0.0)
            else:
                base_vec = (0.0, 1.0, 0.0)
            axis_perp_dir = v_norm(v_cross(axis_dir, base_vec))
        else:
            axis_perp_dir = v_scale(r_perp, 1.0 / r_len)
        axis_perp_point = anchor1

    for atom in rec_list:
        if atom.recordName not in ("ATOM", "HETATM"):
            continue
        if not atom_moves_with_helix(atom, helix_moving):
            continue

        p = (atom.x, atom.y, atom.z)

        # Optional 180-degree pre-flip about a perpendicular axis through axis2_point
        # (explores parallel vs. anti-parallel axis orientations).
        if pre_flip:
            p = rotate_around_line(p, axis2_point, flip_axis_dir, math.pi)

        # 1) phi around helix-1 axis
        if abs(phi) > 1.0e-10:
            p = rotate_around_line(p, axis1_point, axis_dir, phi)

        # 2) d along axis_dir
        if abs(d) > 1.0e-10:
            p = v_add(p, v_scale(axis_dir, d))

        # 3) tau around helix-2 axis
        if abs(tau) > 1.0e-10:
            p = rotate_around_line(p, axis2_base, axis_dir, tau)

        # 4) beta around common perpendicular through anchor1
        if (not axis_parallel) and abs(beta) > 1.0e-10:
            p = rotate_around_line(p, axis_perp_point, axis_perp_dir, beta)  # type: ignore[arg-type]

        atom.update_xyz(*p)


def apply_user_axis_transform_to_helix(
    rec_list: List[pdb_atom_record],
    helix_moving: HelixID,
    axis_dir: Point3D,
    axis_point: Point3D,
    translation: Point3D,
    angle: float,
) -> None:
    """Rotate the movable helix group around a user axis, then translate it in XYZ."""
    axis_dir = v_norm(axis_dir)
    for atom in rec_list:
        if atom.recordName not in ("ATOM", "HETATM"):
            continue
        if not atom_moves_with_helix(atom, helix_moving):
            continue
        p = rotate_around_line((atom.x, atom.y, atom.z), axis_point, axis_dir, angle)
        p = v_add(p, translation)
        atom.update_xyz(*p)


# ---------------------------------------------------------------------------
# Reciprocal-exchange helper (in-memory)
# ---------------------------------------------------------------------------

def apply_reciprocal_exchanges_in_memory(
    rec_list: List[pdb_atom_record],
    exchange_specs,
    cir_shift: int,
    return_metadata: bool = False,
    linker_phosphate_style=None,
):
    """Apply reciprocal exchanges to an in-memory PDB record list.

    This version uses `reciprocal_exchange_pdbV3_3.py` logic, which supports:

      - double / single reciprocal exchanges, and
      - bowtie exchanges (3'-3' + 5'-5' linkages with phosphate-only residues),

    and generates the appropriate LINK records (via edit_pdb_link).

    Returns
    -------
    (atom_rec_list, link_rec_list)
        atom_rec_list includes ATOM/HETATM records plus a TER record after each
        output chain.  link_rec_list contains LINK records for special/inverted
        bonds.
    """
    linker_phosphate_style = rex.coerce_linker_phosphate_style(linker_phosphate_style)

    # Convert legacy V2/V3.5-style exchange_specs into V3.3 dict specs.
    # Optional beta angles are alignment-only metadata and are ignored by
    # the reciprocal-exchange topology layer.
    specs_v3: List[Dict[str, object]] = []
    for spec in exchange_specs:
        pos1, pos2, kind, _beta_deg = unpack_exchange_spec(spec)
        c1, r1 = pos1
        c2, r2 = pos2
        specs_v3.append(
            {
                "pos1": (str(c1).upper(), int(r1)),
                "pos2": (str(c2).upper(), int(r2)),
                "kind": str(kind).lower(),
            }
        )

    bowtie_specs = [sp for sp in specs_v3 if sp.get("kind") == "bowtie"]

    # Extract atom records only (ignore TER for residue grouping).
    atom_recs = [
        r for r in rec_list
        if getattr(r, "recordName", "") in ("ATOM", "HETATM")
        and all(hasattr(r, attr) for attr in ("chainID", "resSeq", "name", "x", "y", "z"))
    ]

    nodes, label_to_idx = rex.build_residue_nodes(atom_recs)
    if not nodes:
        raise ValueError("No residues parsed from rec_list for reciprocal exchange.")

    orig_prev, orig_next = rex.build_original_prev_next(nodes)

    # Store/cut bowtie phosphates (pos2 only). Redirect stdout to stderr so
    # re_helix stdout stays relatively clean.
    with contextlib.redirect_stdout(sys.stderr):
        phos_store = rex.cut_and_store_bowtie_phosphates(nodes, label_to_idx, bowtie_specs)

    # Build original backbone graph and apply exchanges (order-independent).
    g = rex.build_original_graph(len(nodes), orig_next)
    junction_nodes, n_double, n_single, n_bowtie = rex.apply_exchanges_to_graph(
        g, nodes, label_to_idx, orig_prev, specs_v3
    )

    # Build ordered connected components and apply circular permutation.
    base_components = rex.build_ordered_components(g, nodes, junction_nodes, cir_shift)

    # Insert phosphate-only residues for 3'-3' edges.
    used_phos: Set[Tuple[str, int]] = set()
    final_components: List[Dict[str, object]] = []
    for comp in base_components:
        base_order = comp["order"]  # type: ignore[index]
        expanded_order = rex.insert_phosphate_nodes(
            base_order,
            g,
            nodes,
            phos_store,
            used_phos,
            linker_phosphate_style=linker_phosphate_style,
        )
        final_components.append({**comp, "order": expanded_order})

    # Assign new labels and collect atoms (this also appends TER records).
    output_atoms, _phos_new_label = rex.assign_new_labels_and_collect_atoms(final_components, nodes)

    # Generate LINK records.
    link_records, link_counts = rex.build_link_records(final_components, g, nodes, output_atoms)

    num_cycles = sum(1 for comp in final_components if comp.get("is_cycle"))
    num_paths = len(final_components) - num_cycles

    sys.stderr.write(
        f"[reciprocal_exchange] Applied {n_double} double, {n_single} single, {n_bowtie} bowtie exchanges; "
        f"chains: {len(final_components)} ({num_paths} linear, {num_cycles} circular); "
        f"LINK: {len(link_records)} (inv_backbone={link_counts.get('backbone_inverted', 0)}, "
        f"bowtie5to5={link_counts.get('bowtie_5to5', 0)}, bowtie3to3={link_counts.get('bowtie_3to3', 0)}); "
        f"3to3_linker={linker_phosphate_style.record_name}:{linker_phosphate_style.resname}.\n"
    )

    if return_metadata:
        metadata = {
            "specs_v3": specs_v3,
            "label_to_idx": label_to_idx,
            "orig_prev": orig_prev,
            "nodes": nodes,
            "phos_new_label": _phos_new_label,
            "final_components": final_components,
            "link_counts": link_counts,
            "n_double": n_double,
            "n_single": n_single,
            "n_bowtie": n_bowtie,
            "linker_phosphate_style": linker_phosphate_style,
        }
        return output_atoms, link_records, metadata

    return output_atoms, link_records


def _rex_dependency_remark_line() -> Optional[str]:
    if not hasattr(rex, "REMARK_PREFIX"):
        return None
    return (
        f"{rex.REMARK_PREFIX} SOFTWARE_DEPENDENCY "
        f"name={getattr(rex, 'SOFTWARE_NAME', 'reciprocal_exchange_pdb')} "
        f"version={getattr(rex, 'SOFTWARE_VERSION', 'unknown')} "
        f"developer={getattr(rex, 'SOFTWARE_DEVELOPER', 'DiLiuLab')}\n"
    )


def write_reciprocal_exchange_output(
    rec_list: List[pdb_atom_record],
    exchange_specs,
    cir_shift: int,
    command_text: str,
    output_path: str,
    output_stage: str,
    linker_phosphate_style=None,
) -> None:
    linker_phosphate_style = rex.coerce_linker_phosphate_style(linker_phosphate_style)
    rec_list_rex, link_rec_list, rex_metadata = apply_reciprocal_exchanges_in_memory(
        rec_list,
        exchange_specs,
        cir_shift=cir_shift,
        return_metadata=True,
        linker_phosphate_style=linker_phosphate_style,
    )

    rex_header = rex.build_re_script_header_lines(
        software_name=SOFTWARE_NAME,
        software_version=SOFTWARE_VERSION,
        developer=SOFTWARE_DEVELOPER,
        command=command_text,
        output_stage=output_stage,
        atom_rec_list=rec_list_rex,
        specs=rex_metadata.get("specs_v3"),
        label_to_idx=rex_metadata.get("label_to_idx"),
        orig_prev=rex_metadata.get("orig_prev"),
        nodes=rex_metadata.get("nodes"),
        phos_new_label=rex_metadata.get("phos_new_label"),
        component_orders=rex_metadata.get("final_components"),
        link_counts=rex_metadata.get("link_counts"),
        linker_phosphate_style=linker_phosphate_style,
    )

    dependency_line = _rex_dependency_remark_line()
    if dependency_line is not None:
        rex_header.append(dependency_line)

    with open(output_path, "w") as fout:
        rex.write_pdb_with_header(
            rec_list_rex,
            link_rec_list,
            fout,
            header_lines=rex_header,
            reorder_serial=True,
            linker_phosphate_style=linker_phosphate_style,
        )
# ---------------------------------------------------------------------------
# Replication helper (replicate all chains)
# ---------------------------------------------------------------------------

def replicate_all_chains(
    rec_list: List[pdb_atom_record],
    base_helices: List[HelixID],
    exchange_specs,
    explicit_helix_defs: Optional[List[HelixID]] = None,
) -> Tuple[List[ExchangeSpec], List[HelixID]]:
    """
    General replication logic:

    - Take ALL chains present in rec_list (ATOM/HETATM).
    - Rename them to consecutive letters A, B, C, ... according to the
      alphabetical order of the original chain IDs.
    - Replicate the entire set of chains as many times as needed so that
      all chain IDs appearing in exchange_specs exist as actual chains.
    - Propagate helix groupings across copies.

    Arguments
    ---------
    rec_list:
        Current list of pdb_atom_record / pdb_ter_record (input structure).
    base_helices:
        Helix IDs detected BEFORE replication (using either explicit
        helix definitions or auto-detection). Each HelixID is a tuple
        of chain IDs in the original PDB (before renaming).
    exchange_specs:
        List of ((chain1,res1),(chain2,res2),kind).
    explicit_helix_defs:
        If the user supplied (AB), (ABMN), etc., those definitions
        BEFORE replication. If provided, we prefer them over base_helices
        when building helix groupings, because they reflect the user's
        intended rigid groups.

    Returns
    -------
    new_exchange_specs:
        Exchange specs with chain IDs updated to the new scheme and
        including future copies as requested in the original arguments.
    helix_defs_repl:
        Helix definitions (HelixID groups) for ALL copies, to be used
        directly in align_helices_for_exchanges.
    """
    letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

    # Ensure any triplex helix templates have their axis-strand override registered
    # before chain IDs are remapped and copied.
    axis_override_templates = explicit_helix_defs if explicit_helix_defs else base_helices
    for helix_id in axis_override_templates:
        if len(tuple(sorted(helix_id))) == 3:
            _prompt_for_triplex_axis_override(tuple(sorted(helix_id)), "replication template")
    axis_overrides_old = get_helix_axis_overrides()
    axis_range_overrides_old = get_helix_axis_range_overrides()
    move_selections_old = get_helix_move_selections()

    # 1) Gather original chain IDs from the PDB (case-insensitive, but keep their actual chars)
    orig_chains: List[str] = sorted(
        {atom.chainID for atom in rec_list if atom.recordName in ("ATOM", "HETATM")}
    )
    if not orig_chains:
        raise ValueError("No ATOM/HETATM records with chain IDs found for replication.")

    if len(orig_chains) > len(letters):
        raise ValueError(
            f"Cannot remap {len(orig_chains)} chains to consecutive letters; "
            f"max supported is {len(letters)}."
        )

    # 2) Map original chain IDs to A, B, C, ... in alphabetical order
    #    mapping_upper_to_new: 'X'/'x' -> 'A' | 'B' | ...
    mapping_upper_to_new: Dict[str, str] = {}
    for i, ch in enumerate(orig_chains):
        new_ch = letters[i]
        mapping_upper_to_new[ch.upper()] = new_ch

    sys.stderr.write(
        "[re_helix] Replication requested: remapping original chains "
        + ", ".join(orig_chains)
        + " -> "
        + ", ".join(mapping_upper_to_new[ch.upper()] for ch in orig_chains)
        + ".\n"
    )

    # 3) Apply the chain-ID remapping to all atoms (including TER)
    for atom in rec_list:
        if not hasattr(atom, "chainID"):
            continue
        old = atom.chainID
        key = old.upper()
        if key in mapping_upper_to_new:
            atom.update_chainID(mapping_upper_to_new[key])

    # After remapping, base chains are just A, B, C, ...
    base_chains = sorted(
        {atom.chainID for atom in rec_list if atom.recordName in ("ATOM", "HETATM")}
    )
    n_base = len(base_chains)
    base_index: Dict[str, int] = {ch: i for i, ch in enumerate(base_chains)}

    # 4) Decide which helix groupings to use as the "template".
    #    If explicit helix defs were given, use those; otherwise use base_helices.
    if explicit_helix_defs:
        base_groups_raw = explicit_helix_defs
    else:
        base_groups_raw = base_helices

    # Map those groups onto the new chain IDs
    base_helices_new: List[HelixID] = []
    for h in base_groups_raw:
        try:
            new_group = tuple(sorted(mapping_upper_to_new[ch.upper()] for ch in h))
        except KeyError as e:
            raise ValueError(
                f"Helix definition {h} refers to chain '{e.args[0]}', "
                f"which is not present in the PDB."
            )
        base_helices_new.append(new_group)

    sys.stderr.write(
        "[re_helix]  Base helix groups after chain remap: "
        + ", ".join(helix_id_str(h) for h in base_helices_new)
        + "\n"
    )

    base_axis_overrides_new: Dict[HelixID, HelixID] = {}
    for helix_id, axis_chains in axis_overrides_old.items():
        try:
            new_helix = tuple(sorted(mapping_upper_to_new[ch.upper()] for ch in helix_id))
            new_axis = tuple(sorted(mapping_upper_to_new[ch.upper()] for ch in axis_chains))
        except KeyError:
            continue
        base_axis_overrides_new[new_helix] = new_axis

    # 5) Remap exchange specs: old base chains -> new ones; chains that don't
    #    exist in the original PDB (future copies) are left as-is BUT uppercased.
    new_exchange_specs = []
    chains_used: Set[str] = set()

    for spec in exchange_specs:
        (c1, r1), (c2, r2), kind, beta_deg = unpack_exchange_spec(spec)
        u1 = c1.upper()
        u2 = c2.upper()
        nc1 = mapping_upper_to_new.get(u1, u1)  # map if base chain, otherwise future copy
        nc2 = mapping_upper_to_new.get(u2, u2)

        if nc1 not in letters or nc2 not in letters:
            raise ValueError(
                f"Replication mode only supports chain IDs A-Z after remap; got '{nc1}'/'{nc2}'."
            )

        new_exchange_specs.append(((nc1, r1), (nc2, r2), kind, beta_deg))
        chains_used.add(nc1)
        chains_used.add(nc2)

    # 6) Determine how many copies are needed from the chain letters used
    idx_used = [letters.index(ch) for ch in chains_used]
    max_idx = max(idx_used)
    num_copies = max_idx // n_base + 1

    if num_copies * n_base > len(letters):
        raise ValueError(
            f"Need {num_copies} copies of {n_base} chains (= {num_copies * n_base} chains), "
            f"but only {len(letters)} chain IDs A-Z are available."
        )

    sys.stderr.write(
        f"[re_helix]  Replicating {n_base} base chains into {num_copies} "
        f"total copies to satisfy exchange chain IDs.\n"
    )

    # 7) Snapshot atoms from the base copy (copy index 0)
    base_atoms: List[pdb_atom_record] = [
        atom for atom in rec_list if atom.chainID in base_chains
    ]

    # 8) Create additional copies (copy indices 1..num_copies-1)
    for copy_idx in range(1, num_copies):
        offset = copy_idx * n_base
        sys.stderr.write(
            f"[re_helix]   Creating copy {copy_idx} "
            f"(chains {letters[offset]}..{letters[offset + n_base - 1]}).\n"
        )
        for atom in base_atoms:
            cls = type(atom)
            new_atom = cls(atom.string)
            i = base_index[atom.chainID]  # index in base_chains
            new_chain = letters[offset + i]
            new_atom.update_chainID(new_chain)
            rec_list.append(new_atom)

    # 9) Build helix definitions for ALL copies, based on base_helices_new template
    helix_defs_repl: List[HelixID] = []
    for copy_idx in range(num_copies):
        offset = copy_idx * n_base
        for h in base_helices_new:
            new_group: List[str] = []
            for ch in h:
                if ch not in base_index:
                    raise ValueError(
                        f"Helix group {h} references chain '{ch}' which is not in the base chain set."
                    )
                i = base_index[ch]
                new_ch = letters[offset + i]
                new_group.append(new_ch)
            helix_defs_repl.append(tuple(new_group))

    sys.stderr.write(
        "[re_helix]  Helix groups after replication: "
        + ", ".join(helix_id_str(h) for h in helix_defs_repl)
        + "\n"
    )

    axis_overrides_repl = _replicate_axis_overrides(
        base_axis_overrides_new,
        {ch: ch for ch in base_chains},
        base_index,
        num_copies,
        letters,
    )
    set_helix_axis_overrides(axis_overrides_repl)
    axis_range_overrides_repl = _replicate_axis_range_overrides(
        axis_range_overrides_old,
        mapping_upper_to_new,
        base_index,
        num_copies,
        letters,
    )
    set_helix_axis_range_overrides(axis_range_overrides_repl)
    move_selections_repl = _replicate_move_selections(
        move_selections_old,
        mapping_upper_to_new,
        base_index,
        num_copies,
        letters,
    )
    set_helix_move_selections(move_selections_repl)
    if axis_overrides_repl:
        sys.stderr.write(
            "[re_helix]  Axis-strand overrides after replication: "
            + ", ".join(
                f"{helix_id_str(h)} -> {helix_id_str(axis_overrides_repl[h])}"
                for h in sorted(axis_overrides_repl)
            )
            + "\n"
        )
    if axis_range_overrides_repl:
        sys.stderr.write(
            "[re_helix]  Axis-range overrides after replication: "
            + ", ".join(
                f"{helix_id_str(h)} -> {_format_axis_range_spec(axis_range_overrides_repl[h])}"
                for h in sorted(axis_range_overrides_repl)
            )
            + "\n"
        )
    if move_selections_repl:
        sys.stderr.write(
            "[re_helix]  Axis-move selections after replication: "
            + ", ".join(
                f"{helix_id_str(h)} -> {_format_axis_range_spec(move_selections_repl[h])}"
                for h in sorted(move_selections_repl)
            )
            + "\n"
        )

    return new_exchange_specs, helix_defs_repl


# ---------------------------------------------------------------------------
# Helix definition parsing: (AB), (CD), (ABMN), ...
# ---------------------------------------------------------------------------

def parse_helix_definition_tokens(
    tokens: List[str],
) -> Tuple[List[HelixID], List[str]]:
    """
    Extract helix definitions of the form '(AB)', '(CD)', '(ABMN)', etc.
    Returns (helix_defs, remaining_tokens).
    """
    helix_defs: List[HelixID] = []
    remaining: List[str] = []

    for t in tokens:
        if len(t) >= 4 and t[0] == "(" and t[-1] == ")":
            inside = t[1:-1]
            if len(inside) >= 2:
                h: HelixID = tuple(inside)
                helix_defs.append(h)
                continue
        remaining.append(t)

    return helix_defs, remaining


# ---------------------------------------------------------------------------
# Helix-graph building & alignment driver
# ---------------------------------------------------------------------------

def build_helix_pair_graph(
    exchange_specs,
    chain_to_helix: Dict[str, HelixID],
):
    """
    From exchange_specs, build:
        helix_pair_data: map frozenset({H1,H2}) ->
            {helix1, helix2, pairs, beta_angles}
        adjacency: helix -> set(neighbour helices)
    """
    helix_pair_data: Dict[frozenset, Dict[str, object]] = {}
    adjacency: Dict[HelixID, Set[HelixID]] = defaultdict(set)

    for spec in exchange_specs:
        pos1, pos2, _kind, beta_deg = unpack_exchange_spec(spec)
        chain1, res1 = pos1
        chain2, res2 = pos2

        if chain1 not in chain_to_helix:
            raise ValueError(
                f"Chain '{chain1}' from residue {res1}{chain1} is not assigned "
                f"to any helix (check explicit helix defs or P atoms)."
            )
        if chain2 not in chain_to_helix:
            raise ValueError(
                f"Chain '{chain2}' from residue {res2}{chain2} is not assigned "
                f"to any helix (check explicit helix defs or P atoms)."
            )

        helix1 = chain_to_helix[chain1]
        helix2 = chain_to_helix[chain2]

        if helix1 == helix2:
            # intra-helix; does not drive inter-helix alignment
            continue

        key = frozenset((helix1, helix2))
        if key not in helix_pair_data:
            helix_pair_data[key] = {
                "helix1": helix1,
                "helix2": helix2,
                "pairs": [],
                "beta_angles": [],
            }

        entry = helix_pair_data[key]
        h1 = entry["helix1"]  # type: ignore[assignment]
        h2 = entry["helix2"]  # type: ignore[assignment]

        if helix1 == h1 and helix2 == h2:
            pos_h1, pos_h2 = pos1, pos2
        elif helix1 == h2 and helix2 == h1:
            pos_h1, pos_h2 = pos2, pos1
        else:
            raise RuntimeError(
                f"Inconsistent helix orientation for pair {helix1} / {helix2}."
            )

        entry["pairs"].append((pos_h1, pos_h2))  # type: ignore[arg-type]
        entry["beta_angles"].append(beta_deg)  # type: ignore[index]

        adjacency[helix1].add(helix2)
        adjacency[helix2].add(helix1)

    return helix_pair_data, adjacency


def select_pairs_for_alignment(
    pairs_fixed_moving: List[Tuple[Tuple[str, int], Tuple[str, int]]],
    residue_to_P_atom: Dict[Tuple[str, int], pdb_atom_record],
    axis_dir: Tuple[float, float, float],
    axis1_point: Tuple[float, float, float],
) -> List[Tuple[Tuple[str, int], Tuple[str, int]]]:
    """
    Given all P-atom pairs for a helix pair (oriented as (fixed, moving)),
    choose which pairs to actually use for alignment:

      - If there are 1 or 2 valid P pairs (with available P atoms), use them all.
      - If there are >2, select the two pairs whose helix-1 P atoms are
        most separated ALONG THE HELIX-1 AXIS (axis_dir), i.e. the min and
        max projection indices. All other pairs are ignored.

    Returns the chosen subset of pairs in a stable order (by their position
    in the original list).
    """
    axis_dir = v_norm(axis_dir)

    # Filter out any pairs missing P atoms, while tracking their order.
    valid: List[Tuple[int, Tuple[Tuple[str, int], Tuple[str, int]], Tuple[float, float, float]]] = []
    for idx, (pos1, pos2) in enumerate(pairs_fixed_moving):
        atom1 = residue_to_P_atom.get(pos1)
        atom2 = residue_to_P_atom.get(pos2)
        if atom1 is None or atom2 is None:
            sys.stderr.write(
                f"[re_helix]   Warning: skipping pair {pos1[1]}{pos1[0]}-"
                f"{pos2[1]}{pos2[0]} due to missing P atom(s).\n"
            )
            continue
        p1 = (atom1.x, atom1.y, atom1.z)
        valid.append((idx, (pos1, pos2), p1))

    if len(valid) <= 2:
        # Use all remaining valid pairs
        return [pair for idx, pair, p1 in valid]

    # Compute projection along helix-1 axis for each helix-1 P atom
    proj_list: List[Tuple[int, float]] = []
    for j, (idx, pair, p1) in enumerate(valid):
        t = v_dot(v_sub(p1, axis1_point), axis_dir)
        proj_list.append((j, t))

    # Indices in 'valid' with min and max projection
    j_min, _ = min(proj_list, key=lambda x: x[1])
    j_max, _ = max(proj_list, key=lambda x: x[1])

    selected_indices = sorted({j_min, j_max})
    chosen_pairs = [valid[j][1] for j in selected_indices]

    return chosen_pairs


def align_helices_for_exchanges(
    rec_list: List[pdb_atom_record],
    exchange_specs,
    axis_dist: float,
    axis_parallel_flag: bool,
    explicit_helices: Optional[List[HelixID]] = None,
    fix_chain: Optional[str] = None,
    user_axis: Optional[UserAxisDefinition] = None,
) -> None:
    """
    Align helices (or multi-chain helix groups) in-place according to
    exchange_specs, using the P-pair objective.

    If fix_chain is provided (e.g. 'A'), then the helix group containing that
    chain is used as the fixed root in its connected component and is never
    moved; all other helices in that component are aligned relative to it.
    """
    chain_to_res_atoms, residue_to_P_atom, chain_to_P_atoms = build_nucleic_acid_maps(
        rec_list
    )

    if not chain_to_P_atoms:
        raise ValueError("No P atoms found in input PDB; alignment cannot proceed.")

    write_angle_definitions_once()
    if user_axis is not None:
        user_axis_dir, user_axis_point = user_axis
        sys.stderr.write(
            "[re_helix] Using user-defined alignment axis: "
            f"point=({user_axis_point[0]:.3f}, {user_axis_point[1]:.3f}, {user_axis_point[2]:.3f}), "
            f"direction=({user_axis_dir[0]:.6f}, {user_axis_dir[1]:.6f}, {user_axis_dir[2]:.6f}). "
            "Fixed/moving helix-axis estimation and axis_dist are skipped; each movable helix "
            "is optimized by rotation around this line plus a full XYZ translation.\n"
        )

    # Build chain -> helix map
    if explicit_helices:
        chain_to_helix = build_chain_to_helix_from_defs(explicit_helices, chain_to_P_atoms)
        sys.stderr.write(
            "[re_helix] Using user-defined helices: "
            + ", ".join(helix_id_str(tuple(sorted(h))) for h in explicit_helices)
            + "\n"
        )
    else:
        chain_to_helix = compute_chain_partner_map(chain_to_P_atoms)
        seen_helices = sorted({chain_to_helix[ch] for ch in chain_to_helix})
        sys.stderr.write(
            "[re_helix] Using automatic helix detection (P-atom proximity). "
            "Detected helices: "
            + ", ".join(helix_id_str(h) for h in seen_helices)
            + "\n"
        )

    axis_range_defs = get_helix_axis_range_definitions()
    axis_move_defs = get_helix_axis_move_definitions()
    if axis_range_defs:
        chain_to_helix = apply_axis_coupling_to_chain_map(
            chain_to_helix,
            axis_range_defs,
            axis_move_defs,
            chain_to_P_atoms,
        )
    if axis_range_defs:
        validate_axis_range_definitions(axis_range_defs, chain_to_P_atoms, chain_to_helix)
        sys.stderr.write(
            "[re_helix] Using axis residue-range overrides: "
            + ", ".join(_format_axis_range_spec(spec) for spec in axis_range_defs)
            + "\n"
        )
    if axis_move_defs:
        sys.stderr.write(
            "[re_helix] Using axis-coupled move selections: "
            + ", ".join(_format_axis_range_spec(spec) if spec else "(none)" for spec in axis_move_defs)
            + "\n"
        )

    # Optional fixed helix (by chain)
    fix_helix: Optional[HelixID] = None
    if fix_chain is not None:
        if fix_chain not in chain_to_helix:
            raise ValueError(
                f"Requested fixed chain '{fix_chain}' is not assigned to any helix "
                f"(check helix definitions or P atoms)."
            )
        fix_helix = chain_to_helix[fix_chain]
        sys.stderr.write(
            f"[re_helix] Helix {helix_id_str(fix_helix)} containing chain "
            f"'{fix_chain}' will be kept fixed in its component.\n"
        )

    # Build helix pair graph
    helix_pair_data, adjacency = build_helix_pair_graph(exchange_specs, chain_to_helix)

    if not helix_pair_data:
        sys.stderr.write(
            "[re_helix] No inter-helix P-atom pairs found in exchanges; "
            "nothing to align.\n"
        )
        return

    helix_nodes: Set[HelixID] = set(adjacency.keys())
    visited: Set[HelixID] = set()

    def process_component(root: HelixID, note: str = "") -> None:
        """
        Process one connected component of the helix graph, starting from 'root'.
        'root' is treated as fixed; every other helix in this component is moved
        at most once, relative to an already-fixed neighbour.
        """
        visited.add(root)
        msg_root = f"[re_helix] Processing helix component rooted at {helix_id_str(root)}"
        if note:
            msg_root += f" ({note})"
        msg_root += ".\n"
        sys.stderr.write(msg_root)

        queue: List[HelixID] = [root]

        while queue:
            fixed = queue.pop(0)
            for neighbour in sorted(adjacency[fixed]):
                if neighbour in visited:
                    # This edge would re-align two already aligned helices
                    sys.stderr.write(
                        f"[re_helix]  Skipping already-aligned helix pair "
                        f"{helix_id_str(fixed)}--{helix_id_str(neighbour)} "
                        f"(cycle / redundant constraint).\n"
                    )
                    continue

                key = frozenset((fixed, neighbour))
                entry = helix_pair_data[key]
                h1 = entry["helix1"]  # type: ignore[assignment]
                h2 = entry["helix2"]  # type: ignore[assignment]
                pairs_h1_h2: List[Tuple[Tuple[str, int], Tuple[str, int]]] = entry["pairs"]  # type: ignore[assignment]
                beta_angles_h1_h2: List[Optional[float]] = entry.get(
                    "beta_angles",
                    entry.get("rho_angles", []),
                )  # type: ignore[assignment]

                # Orient pairs so first pos is on 'fixed', second on 'neighbour'.
                # The beta-angle annotation follows the same site index. Its sign is
                # not changed when the fixed/moving roles are reversed because the
                # beta axis direction is reconstructed from the current fixed-to-moving
                # geometry.
                if fixed == h1 and neighbour == h2:
                    pairs_fixed_moving = pairs_h1_h2
                    beta_angles_fixed_moving = beta_angles_h1_h2
                elif fixed == h2 and neighbour == h1:
                    pairs_fixed_moving = [(b, a) for (a, b) in pairs_h1_h2]
                    beta_angles_fixed_moving = beta_angles_h1_h2
                else:
                    raise RuntimeError(
                        f"Internal error: helix {helix_id_str(fixed)} / "
                        f"{helix_id_str(neighbour)} not found in pair data."
                    )

                sys.stderr.write(
                    f"[re_helix]  Aligning helix {helix_id_str(fixed)} (fixed) "
                    f"with {helix_id_str(neighbour)} (moving); "
                    f"{len(pairs_fixed_moving)} P-atom pair(s) specified.\n"
                )

                fixed_beta_deg: Optional[float] = None
                fixed_beta_rad: Optional[float] = None
                defined_beta_angles = [angle for angle in beta_angles_fixed_moving if angle is not None]
                if defined_beta_angles:
                    angle_list = ", ".join(f"{angle:.3g}°" for angle in defined_beta_angles)
                    if user_axis is not None:
                        sys.stderr.write(
                            f"[re_helix]   Warning: fixed beta-angle definition(s) "
                            f"({angle_list}) for helix pair "
                            f"{helix_id_str(fixed)}/{helix_id_str(neighbour)} ignored because "
                            "user-defined-axis mode uses only rotation around the supplied "
                            "axis plus XYZ translation.\n"
                        )
                    elif axis_parallel_flag:
                        reason = "--axis_parallel y keeps helix axes parallel"
                        if len(pairs_fixed_moving) != 1:
                            reason += (
                                f" and {len(pairs_fixed_moving)} reciprocal-exchange sites "
                                f"connect this helix pair"
                            )
                        sys.stderr.write(
                            f"[re_helix]   Warning: fixed beta-angle definition(s) "
                            f"({angle_list}) for helix pair "
                            f"{helix_id_str(fixed)}/{helix_id_str(neighbour)} ignored because "
                            f"{reason}.\n"
                        )
                    elif len(pairs_fixed_moving) != 1:
                        sys.stderr.write(
                            f"[re_helix]   Warning: fixed beta-angle definition(s) "
                            f"({angle_list}) for helix pair "
                            f"{helix_id_str(fixed)}/{helix_id_str(neighbour)} ignored because "
                            f"{len(pairs_fixed_moving)} reciprocal-exchange sites connect this "
                            f"helix pair; fixed beta is used only for exactly one site.\n"
                        )
                    else:
                        fixed_beta_deg = defined_beta_angles[0]
                        fixed_beta_rad = fixed_beta_deg * math.pi / 180.0
                        sys.stderr.write(
                            f"[re_helix]    Using fixed single-site beta angle: "
                            f"{fixed_beta_deg:.2f}°.\n"
                        )

                # Identify which chains within each helix group actually
                # participate in these P-pairs; use them to define the axes.
                chains_fixed_for_axis: Set[str] = set()
                chains_moving_for_axis: Set[str] = set()
                for (c1, r1), (c2, r2) in pairs_fixed_moving:
                    chains_fixed_for_axis.add(c1)
                    chains_moving_for_axis.add(c2)

                fixed_positions = [(c1, r1) for (c1, r1), (_c2, _r2) in pairs_fixed_moving]
                moving_positions = [(c2, r2) for (_c1, _r1), (c2, r2) in pairs_fixed_moving]

                if user_axis is not None:
                    axis_dir, axis_point = user_axis
                    effective_pairs = select_pairs_for_alignment(
                        pairs_fixed_moving,
                        residue_to_P_atom,
                        axis_dir,
                        axis_point,
                    )
                    if not effective_pairs:
                        sys.stderr.write(
                            f"[re_helix]   Warning: no usable P-atom pairs "
                            f"for helix pair {helix_id_str(fixed)}/{helix_id_str(neighbour)}; "
                            f"skipping user-axis alignment.\n"
                        )
                        visited.add(neighbour)
                        queue.append(neighbour)
                        continue

                    pair_strs = []
                    for (c1, r1), (c2, r2) in effective_pairs:
                        pair_strs.append(f"{r1}{c1}-{r2}{c2}")
                    sys.stderr.write(
                        "[re_helix]    Using P pairs for user-axis alignment: "
                        + ", ".join(pair_strs)
                        + "\n"
                    )

                    objective, p1_list, _p2_list, translation_for_angle = (
                        build_user_axis_transform_objective(
                            effective_pairs,
                            residue_to_P_atom,
                            axis_dir,
                            axis_point,
                        )
                    )
                    if not p1_list:
                        sys.stderr.write(
                            f"[re_helix]   Warning: no valid P-atom coordinates "
                            f"for helix pair {helix_id_str(fixed)}/{helix_id_str(neighbour)}; "
                            f"skipping user-axis alignment.\n"
                        )
                        visited.add(neighbour)
                        queue.append(neighbour)
                        continue

                    best_params, best_val = coordinate_descent(
                        objective,
                        [0.0],
                        [math.pi / 2.0],
                        angle_indices={0},
                    )
                    angle_opt = best_params[0]
                    translation_opt = translation_for_angle(angle_opt)
                    angle_deg = angle_opt * 180.0 / math.pi
                    sys.stderr.write(
                        f"[re_helix]   Optimised user-axis alignment: translation = "
                        f"({translation_opt[0]:.3f}, {translation_opt[1]:.3f}, "
                        f"{translation_opt[2]:.3f}) A; angle = {angle_deg:.2f}°; "
                        f"sum(dist^2) = {best_val:.3f}.\n"
                    )

                    apply_user_axis_transform_to_helix(
                        rec_list,
                        neighbour,
                        axis_dir,
                        axis_point,
                        translation_opt,
                        angle_opt,
                    )

                    visited.add(neighbour)
                    queue.append(neighbour)
                    continue

                axis_ranges_fixed = find_axis_range_definition_for_positions(fixed, fixed_positions)
                axis_ranges_moving = find_axis_range_definition_for_positions(neighbour, moving_positions)
                axis_range_overrides = get_helix_axis_range_overrides()
                if axis_ranges_fixed is None:
                    axis_ranges_fixed = axis_range_overrides.get(tuple(sorted(fixed)))
                if axis_ranges_moving is None:
                    axis_ranges_moving = axis_range_overrides.get(tuple(sorted(neighbour)))

                if axis_ranges_fixed is not None:
                    sys.stderr.write(
                        f"[re_helix]    Fixed-axis residue range: {_format_axis_range_spec(axis_ranges_fixed)}\n"
                    )
                if axis_ranges_moving is not None:
                    sys.stderr.write(
                        f"[re_helix]    Moving-axis residue range: {_format_axis_range_spec(axis_ranges_moving)}\n"
                    )

                # 1) Align axes & set axis_dist (axes from participating chains,
                #    rigid transform applied to the full helix groups).
                axis_dir, center_fixed, center_moving = align_axes_for_pair(
                    rec_list,
                    fixed,
                    neighbour,
                    chain_to_P_atoms,
                    axis_dist,
                    subset_fixed=chains_fixed_for_axis,
                    subset_moving=chains_moving_for_axis,
                    axis_ranges_fixed=axis_ranges_fixed,
                    axis_ranges_moving=axis_ranges_moving,
                )

                # 1b) Choose which P pairs to actually use (possibly 2 extremes)
                effective_pairs = select_pairs_for_alignment(
                    pairs_fixed_moving,
                    residue_to_P_atom,
                    axis_dir,
                    center_fixed,
                )

                if not effective_pairs:
                    sys.stderr.write(
                        f"[re_helix]   Warning: no usable P-atom pairs "
                        f"for helix pair {helix_id_str(fixed)}/{helix_id_str(neighbour)}; "
                        f"skipping fine optimisation.\n"
                    )
                    visited.add(neighbour)
                    queue.append(neighbour)
                    continue

                # Report the specific P pairs used for alignment
                pair_strs = []
                for (c1, r1), (c2, r2) in effective_pairs:
                    pair_strs.append(f"{r1}{c1}-{r2}{c2}")
                sys.stderr.write(
                    "[re_helix]    Using P pairs for alignment: " + ", ".join(pair_strs) + "\n"
                )

                # 2) Build objective(s) + optimise parameters for this helix pair
                #    (using only effective_pairs).
                #
                # V2.1 update: With --axis_parallel y we interpret "parallel" as
                # "parallel OR anti-parallel" and explore both by optionally
                # applying a 180-degree pre-flip of helix 2 about a perpendicular
                # axis through its helix center.

                best_pre_flip = False

                if axis_parallel_flag:
                    sys.stderr.write(
                        "[re_helix]    axis_parallel=y: exploring parallel vs. anti-parallel axis orientations (0°/180° pre-flip).\n"
                    )

                    best_params: Optional[List[float]] = None
                    best_val: float = float("inf")
                    best_anchor1: Optional[Tuple[float, float, float]] = None

                    # Try both orientations: no flip (parallel) and 180-degree flip (anti-parallel).
                    for pre_flip_flag in (False, True):
                        objective, p1_list, _p2_list, d0, anchor1 = build_pair_objective(
                            effective_pairs,
                            residue_to_P_atom,
                            axis_dir,
                            center_fixed,
                            center_moving,
                            axis_parallel_flag,
                            pre_flip=pre_flip_flag,
                        )

                        if not p1_list:
                            continue

                        # params = [d, tau, phi]
                        x0 = [d0, 0.0, 0.0]
                        steps0 = [3.4, math.pi / 2.0, math.pi / 2.0]
                        angle_indices = {1, 2}

                        params, val = coordinate_descent(
                            objective,
                            x0,
                            steps0,
                            angle_indices=angle_indices,
                        )

                        orient = "parallel" if not pre_flip_flag else "anti-parallel"
                        sys.stderr.write(
                            f"[re_helix]     Trial ({orient}): sum(dist^2) = {val:.3f}.\n"
                        )

                        # Tie-break in favour of the first (no-flip) solution.
                        if val < best_val - 1.0e-8:
                            best_val = val
                            best_params = params
                            best_anchor1 = anchor1
                            best_pre_flip = pre_flip_flag

                    if best_params is None or best_anchor1 is None:
                        sys.stderr.write(
                            f"[re_helix]   Warning: no valid P-atom coordinates "
                            f"for helix pair {helix_id_str(fixed)}/{helix_id_str(neighbour)}; "
                            f"skipping fine optimisation.\n"
                        )
                        visited.add(neighbour)
                        queue.append(neighbour)
                        continue

                    # Unpack the chosen solution.
                    d_opt, tau_opt, phi_opt = best_params
                    beta_opt = 0.0
                    chosen_orient = "parallel" if not best_pre_flip else "anti-parallel"
                    best_val_final = best_val
                    anchor1_final = best_anchor1

                else:
                    # Original behaviour for --axis_parallel n, except that a valid
                    # single-site fixed-beta annotation holds beta constant and optimises
                    # only d/tau/phi.
                    objective, p1_list, _p2_list, d0, anchor1 = build_pair_objective(
                        effective_pairs,
                        residue_to_P_atom,
                        axis_dir,
                        center_fixed,
                        center_moving,
                        axis_parallel_flag,
                        fixed_beta=fixed_beta_rad,
                    )

                    if not p1_list:
                        sys.stderr.write(
                            f"[re_helix]   Warning: no valid P-atom coordinates "
                            f"for helix pair {helix_id_str(fixed)}/{helix_id_str(neighbour)}; "
                            f"skipping fine optimisation.\n"
                        )
                        visited.add(neighbour)
                        queue.append(neighbour)
                        continue

                    if fixed_beta_rad is None:
                        # params = [d, tau, phi, beta]
                        x0 = [d0, 0.0, 0.0, 0.0]
                        steps0 = [3.4, math.pi / 2.0, math.pi / 2.0, math.pi / 2.0]
                        angle_indices = {1, 2, 3}
                    else:
                        # params = [d, tau, phi]; beta is held at fixed_beta_rad
                        x0 = [d0, 0.0, 0.0]
                        steps0 = [3.4, math.pi / 2.0, math.pi / 2.0]
                        angle_indices = {1, 2}

                    best_params, best_val = coordinate_descent(
                        objective,
                        x0,
                        steps0,
                        angle_indices=angle_indices,
                    )

                    if fixed_beta_rad is None:
                        d_opt, tau_opt, phi_opt, beta_opt = best_params
                    else:
                        d_opt, tau_opt, phi_opt = best_params
                        beta_opt = fixed_beta_rad
                    chosen_orient = "n/a"
                    best_val_final = best_val
                    anchor1_final = anchor1

                tau_deg = tau_opt * 180.0 / math.pi
                phi_deg = phi_opt * 180.0 / math.pi
                beta_deg = beta_opt * 180.0 / math.pi

                msg = (
                    f"[re_helix]   Optimised: d = {d_opt:.3f} Å, "
                    f"tau = {tau_deg:.2f}°, phi = {phi_deg:.2f}°"
                )
                if not axis_parallel_flag:
                    msg += f", beta = {beta_deg:.2f}°"
                    if fixed_beta_rad is not None:
                        msg += " (fixed single-site angle)"
                else:
                    msg += f" (axis orientation: {chosen_orient})"
                msg += f"; sum(dist^2) = {best_val_final:.3f}.\n"
                sys.stderr.write(msg)

                # 3) Apply transform
                apply_transform_to_helix(
                    rec_list,
                    neighbour,
                    axis_dir,
                    center_fixed,
                    center_moving,
                    d_opt,
                    tau_opt,
                    phi_opt,
                    beta_opt,
                    axis_parallel_flag,
                    anchor1_final,
                    pre_flip=best_pre_flip,
                )

                visited.add(neighbour)
                queue.append(neighbour)

    # If a fixed helix is specified and participates in any inter-helix pairs,
    # process its component first.
    if fix_helix is not None and fix_helix in helix_nodes:
        process_component(fix_helix, note="fixed helix")
    elif fix_helix is not None and fix_helix not in helix_nodes:
        sys.stderr.write(
            f"[re_helix] Warning: fixed helix {helix_id_str(fix_helix)} "
            f"has no inter-helix P-atom pairs; it will not be involved in "
            f"alignment operations.\n"
        )

    # Process remaining components
    for root in sorted(helix_nodes):
        if root in visited:
            continue
        process_component(root)


# ---------------------------------------------------------------------------
# GUI helpers
# ---------------------------------------------------------------------------

def _default_output_base_from_pdb_path(pdb_path: str) -> str:
    """Return the GUI's default output base for an input PDB-like path."""
    path = pdb_path.strip()
    if path.lower().endswith(".pdb"):
        return path[:-4]
    return path


_GUI_HELP_TEXT: Dict[str, str] = {
    "pdb_in": """Input PDB file containing the nucleic-acid strands to align.

Example:
  TT_helixAB_33.pdb""",
    "helix_defs": """Optional explicit helix-group definitions. Each parenthesized token
defines one rigid helix group.

Examples:
  (AB) (CD)
  (ABMN)

Leave blank to use automatic helix detection.""",
    "output": """Base path for the output files. In normal mode the script writes:
  <base>_aligned.pdb
  <base>_aligned_rex.pdb

In RE-only mode the script writes:
  <base>_rex.pdb

Example:
  my_model""",
    "axis_dist": """Target distance between the two helix axes during the coarse axis
alignment step.

Units: Å
Default: 22.0""",
    "axis_parallel": """Choose whether the aligned helix axes must stay parallel.

y: keep the axes parallel (or anti-parallel) and optimise d/tau/phi.
n: also allow a tilt of helix 2 (beta), so the axes may become non-parallel.

For --axis_parallel n, a single-site pair may optionally define a fixed beta
angle in the pair row. That controlled angle is ignored with a warning when
more than one site connects the same helix pair, or when --axis_parallel y is
selected.""",
    "fix": """Optional chain ID whose helix group should remain fixed while the other
connected helices are aligned around it.

Example:
  A""",
    "replicate": """Replicate the full set of input chains when needed. The original chains
are first renamed to consecutive letters A, B, C, ... and extra full
copies are added as required by the exchange chain IDs.

When replication is active, --axis_range definitions are NOT copied
automatically to every copy; instead you can target final replicated
chain IDs directly.""",
    "re_only": """Apply reciprocal exchanges directly to the input structure
without running helix alignment first.

The output is <base>_rex.pdb. Explicit helix definitions and beta-angle values
may be left in the rows, but only the residue pairs and exchange kinds are
used for reciprocal exchange.""",
    "cir_shift": """Shift applied when writing circular strands after reciprocal exchange,
so the output nick is moved away from junction residues when possible.

Units: nt
Default: 8""",
    "linker_phosphate_resname": """Residue name used for phosphate-only residues inserted at bowtie 3'-3' linkages.

Options:
  X33: default custom linker residue, written as HETATM.
  custom 1-3 character name: written as HETATM with that residue name.
  DA or dA: written as regular ATOM DA, with only P, OP1, and OP2 kept, to avoid needing a custom residue definition during Phenix relaxation.

The command-line option is --linker_phosphate_resname.""",
    "pairs": """Exchange pairs are entered as residue positions joined across helices.
Blank or incomplete rows are ignored. Use at least one complete row with
pos1, pos2, and kind.

For many pairs, use the CLI-style single-line field below the rows. When that
line is filled, the individual rows above it are ignored.

In normal mode these pairs drive helix alignment and reciprocal exchange. In
RE-only mode they only drive reciprocal exchange.

pos1 / pos2 format: residue number + chain ID, such as 30A or A30
Units: nt + chain ID
beta angle: optional fixed beta angle in degrees for --axis_parallel n when
exactly one reciprocal-exchange site connects that helix pair; leave blank for
normal behavior. Ignored with a warning for multi-site helix pairs or when
--axis_parallel y is selected. This helps make controlled interhelical angles.
Legacy rho-angle wording means the same beta angle.
kind: d = double, s = single, b = bowtie

Example rows:
  30A   8D        d
  26A   9C   90   d
  13B  24C        s""",
    "pair_args": """Optional CLI-style exchange-pair argument line.

Use this when many exchange pairs are easier to paste as one command-line-style
string, such as:
  30A 8D d 26A 9C 90 d 13B 24C s

If this field is filled, the individual pair rows above it are ignored.""",
    "axis_range": """Optional chain or residue windows used to define a helical axis.

Format: comma-separated chain IDs or chain-specific ranges
Units: nt
Examples:
  A,B
  B26-B60,A1-A35

The adjacent "move with axis" field can list additional chains or residue
windows that should receive the same rigid transform, for example:
  C,D
  C1-C50,D

This can define a triplex-like group without stdin prompts: use axis A,B and
move C. The axis is fit from A/B while C moves together with that axis.
In replication mode, enter the final post-replication chain IDs you want
to target, or use a base-template row that covers the input chains before
replication to propagate it across copies. Blank rows are ignored.""",
    "user_axis": """Optional user-defined alignment axis.

Provide a direction vector and one point on the axis. When enabled, re_helix
does not estimate the fixed and moving helix axes from P atoms. Each movable
helix is optimized by rotating it around the supplied line and translating it
freely in XYZ.

Direction and point format:
  x y z

The direction vector is normalized automatically and cannot be zero.""",
}


def _build_equivalent_cli_command(
    script_path: str,
    pdb_in: str,
    helix_def_text: str,
    pair_rows: List[Dict[str, str]],
    pair_args_text: str,
    axis_range_rows: List[str],
    axis_move_rows: List[str],
    use_user_axis: bool,
    user_axis_dir: List[str],
    user_axis_point: List[str],
    output_base: str,
    axis_dist: str,
    axis_parallel: str,
    fix_chain: str,
    replicate: bool,
    re_only: bool,
    cir_shift: str,
    linker_phosphate_resname: str,
) -> List[str]:
    cmd = [sys.executable, script_path, pdb_in]

    helix_def_text = helix_def_text.strip()
    if helix_def_text:
        cmd.extend(helix_def_text.split())

    pair_args_text = pair_args_text.strip()
    if pair_args_text:
        try:
            pair_tokens = shlex.split(pair_args_text)
        except ValueError as exc:
            raise ValueError(f"Could not parse CLI-style exchange-pair line: {exc}") from exc
        if not pair_tokens:
            raise ValueError("Please provide at least one exchange pair in the CLI-style line.")
        cmd.extend(pair_tokens)
    else:
        n_pairs = 0
        for row in pair_rows:
            pos1 = row.get("pos1", "").strip()
            pos2 = row.get("pos2", "").strip()
            beta_angle = row.get("beta_angle", row.get("rho_angle", "")).strip()
            kind = row.get("kind", "").strip()
            if not pos1 or not pos2 or not kind:
                continue
            if beta_angle:
                cmd.extend([pos1, pos2, beta_angle, kind])
            else:
                cmd.extend([pos1, pos2, kind])
            n_pairs += 1

        if n_pairs == 0:
            raise ValueError("Please provide at least one complete exchange pair.")

    if not use_user_axis:
        for idx, spec in enumerate(axis_range_rows):
            spec = spec.strip()
            move_spec = axis_move_rows[idx].strip() if idx < len(axis_move_rows) else ""
            if move_spec and not spec:
                raise ValueError("A move-with-axis entry requires an axis definition in the same row.")
            if spec:
                cmd.extend(["--axis_range", spec])
                if move_spec:
                    cmd.extend(["--axis_move", move_spec])

    if use_user_axis:
        if len(user_axis_dir) != 3 or len(user_axis_point) != 3:
            raise ValueError("User axis direction and point must each have three values.")
        dir_values = [value.strip() for value in user_axis_dir]
        point_values = [value.strip() for value in user_axis_point]
        if not all(dir_values) or not all(point_values):
            raise ValueError("Please fill all user axis direction and point coordinates.")
        try:
            normalize_user_axis_definition(
                [float(value) for value in dir_values],
                [float(value) for value in point_values],
            )
        except ValueError as exc:
            raise ValueError(f"Invalid user-defined axis: {exc}") from exc
        cmd.extend(["--user_axis_dir", *dir_values])
        cmd.extend(["--user_axis_point", *point_values])

    if output_base.strip():
        cmd.extend(["-o", output_base.strip()])

    if axis_dist.strip():
        cmd.extend(["--axis_dist", axis_dist.strip()])

    axis_parallel = axis_parallel.strip().lower()
    if axis_parallel in ("y", "n"):
        cmd.extend(["--axis_parallel", axis_parallel])

    if fix_chain.strip():
        cmd.extend(["--fix", fix_chain.strip()])

    if replicate:
        cmd.append("--replicate")

    if re_only:
        cmd.append("--re_only")

    if cir_shift.strip():
        cmd.extend(["--cir_shift", cir_shift.strip()])

    linker_phosphate_resname = linker_phosphate_resname.strip()
    if linker_phosphate_resname and linker_phosphate_resname.upper() != "X33":
        cmd.extend(["--linker_phosphate_resname", linker_phosphate_resname])

    return cmd


def _launch_gui() -> None:
    try:
        import tkinter as tk
        from tkinter import filedialog, messagebox, scrolledtext, ttk
    except Exception as exc:
        raise RuntimeError(
            "GUI mode requires tkinter, but it could not be imported on this system."
        ) from exc

    root = tk.Tk()
    root.title(APP_TITLE)
    root.geometry("1320x930")

    gui_style = ttk.Style(root)
    gui_style.configure("Title.TLabel", font=("TkDefaultFont", 15, "bold"))
    gui_style.configure("Section.TLabelframe.Label", font=("TkDefaultFont", 10, "bold"))

    script_dir = Path(__file__).resolve().parent
    icon_candidates = [
        script_dir / "assets" / "icon.png",
        script_dir / "re_helix_lib" / "icon.png",
    ]
    icon_path = next((path for path in icon_candidates if path.exists()), None)
    if icon_path is not None:
        try:
            icon_image = tk.PhotoImage(file=str(icon_path))
            root.iconphoto(True, icon_image)
            root._re_helix_icon = icon_image
        except Exception:
            pass

    pair_count_var = tk.StringVar(value="3")
    axis_count_var = tk.StringVar(value="0")
    pdb_var = tk.StringVar()
    helix_defs_var = tk.StringVar()
    output_var = tk.StringVar()
    output_auto_state = {"last_default": ""}
    axis_dist_var = tk.StringVar(value="22.0")
    axis_parallel_var = tk.StringVar(value="y")
    user_axis_var = tk.BooleanVar(value=False)
    user_axis_dir_vars = [tk.StringVar(value=value) for value in ("0", "0", "1")]
    user_axis_point_vars = [tk.StringVar(value=value) for value in ("0", "0", "0")]
    fix_chain_var = tk.StringVar()
    replicate_var = tk.BooleanVar(value=False)
    re_only_var = tk.BooleanVar(value=False)
    cir_shift_var = tk.StringVar(value="8")
    linker_phosphate_resname_var = tk.StringVar(value="X33")
    pair_args_var = tk.StringVar()

    pair_widgets: List[Dict[str, object]] = []
    axis_widgets: List[object] = []
    user_axis_input_widgets: List[object] = []
    process_state: Dict[str, object] = {"proc": None, "queue": queue.Queue()}

    def show_help(title: str, message: str) -> None:
        messagebox.showinfo(title, message)

    def make_help_button(parent, title: str, key: str):
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
            command=lambda: show_help(title, _GUI_HELP_TEXT[key]),
        )

    def refresh_default_output_base(*_args) -> None:
        new_default = _default_output_base_from_pdb_path(pdb_var.get())
        if not new_default:
            output_auto_state["last_default"] = ""
            return
        current_output = output_var.get().strip()
        previous_default = output_auto_state["last_default"]
        if not current_output or current_output == previous_default:
            output_var.set(new_default)
        output_auto_state["last_default"] = new_default

    pdb_var.trace_add("write", refresh_default_output_base)

    root.rowconfigure(0, weight=1)
    root.columnconfigure(0, weight=1)

    canvas_host = ttk.Frame(root)
    canvas_host.grid(row=0, column=0, sticky="nsew")
    canvas_host.rowconfigure(0, weight=1)
    canvas_host.columnconfigure(0, weight=1)

    content_canvas = tk.Canvas(canvas_host, highlightthickness=0, borderwidth=0)
    content_canvas.grid(row=0, column=0, sticky="nsew")

    canvas_scrollbar = ttk.Scrollbar(
        canvas_host,
        orient="vertical",
        command=content_canvas.yview,
    )
    content_canvas.configure(yscrollcommand=canvas_scrollbar.set)

    scroll_root = ttk.Frame(content_canvas)
    scroll_window = content_canvas.create_window((0, 0), window=scroll_root, anchor="nw")

    scroll_state = {"visible": False, "pending": False}

    def refresh_scrollbar() -> None:
        scroll_state["pending"] = False
        if not root.winfo_exists():
            return
        root.update_idletasks()
        canvas_width = max(content_canvas.winfo_width(), 1)
        content_canvas.itemconfigure(scroll_window, width=canvas_width)
        content_width = max(scroll_root.winfo_reqwidth(), outer.winfo_reqwidth(), canvas_width)
        content_height = max(scroll_root.winfo_reqheight(), outer.winfo_reqheight())
        content_canvas.configure(scrollregion=(0, 0, content_width, content_height))
        viewport_height = max(content_canvas.winfo_height(), 0)
        need_scrollbar = content_height > viewport_height + 2
        if need_scrollbar and not scroll_state["visible"]:
            canvas_scrollbar.grid(row=0, column=1, sticky="ns")
            scroll_state["visible"] = True
        elif (not need_scrollbar) and scroll_state["visible"]:
            canvas_scrollbar.grid_remove()
            content_canvas.yview_moveto(0.0)
            scroll_state["visible"] = False

    def schedule_scrollbar_refresh(*_args) -> None:
        if scroll_state["pending"]:
            return
        scroll_state["pending"] = True
        root.after_idle(refresh_scrollbar)

    def _on_canvas_configure(event) -> None:
        content_canvas.itemconfigure(scroll_window, width=max(event.width, 1))
        schedule_scrollbar_refresh()

    content_canvas.bind("<Configure>", _on_canvas_configure)
    scroll_root.bind("<Configure>", schedule_scrollbar_refresh)

    def _coerce_gui_count(raw_value: str, minimum: int, maximum: int) -> Optional[int]:
        value_text = str(raw_value).strip()
        if not value_text:
            return None
        if not value_text.isdigit():
            return None
        return max(minimum, min(maximum, int(value_text)))

    def _validate_count_text(proposed: str) -> bool:
        if proposed == "":
            return True
        return proposed.isdigit()

    validate_count_cmd = (root.register(_validate_count_text), "%P")

    outer = ttk.Frame(scroll_root, padding=10)
    outer.pack(fill="both", expand=True)

    title_label = ttk.Label(
        outer,
        text=APP_TITLE,
        style="Title.TLabel",
        wraplength=1100,
        justify="left",
    )
    title_label.pack(fill="x", padx=2, pady=(0, 8))

    basics = ttk.LabelFrame(outer, text="Inputs", padding=10, style="Section.TLabelframe")
    basics.pack(fill="x", padx=2, pady=4)

    ttk.Label(basics, text="Input PDB file").grid(row=0, column=0, sticky="w")
    pdb_entry = ttk.Entry(basics, textvariable=pdb_var, width=84)
    pdb_entry.grid(row=0, column=1, sticky="ew", padx=4)

    def browse_pdb() -> None:
        path = filedialog.askopenfilename(
            title="Choose input PDB",
            filetypes=[("PDB-like files", "*.pdb *.txt *"), ("All files", "*")],
        )
        if path:
            pdb_var.set(path)

    ttk.Button(basics, text="Browse", command=browse_pdb).grid(row=0, column=2, padx=4)
    make_help_button(basics, "Input PDB", "pdb_in").grid(row=0, column=3, padx=(0, 4), sticky="w")

    ttk.Label(basics, text="Helix defs").grid(row=1, column=0, sticky="w", pady=(6, 0))
    ttk.Entry(basics, textvariable=helix_defs_var, width=84).grid(
        row=1, column=1, columnspan=2, sticky="ew", padx=4, pady=(6, 0)
    )
    make_help_button(basics, "Helix defs", "helix_defs").grid(
        row=1, column=3, padx=(0, 4), pady=(6, 0), sticky="w"
    )

    ttk.Label(basics, text="Output base").grid(row=2, column=0, sticky="w", pady=(6, 0))
    ttk.Entry(basics, textvariable=output_var, width=84).grid(
        row=2, column=1, columnspan=2, sticky="ew", padx=4, pady=(6, 0)
    )
    make_help_button(basics, "Output base", "output").grid(
        row=2, column=3, padx=(0, 4), pady=(6, 0), sticky="w"
    )

    basics.columnconfigure(1, weight=1)

    options = ttk.LabelFrame(outer, text="Options", padding=10, style="Section.TLabelframe")
    options.pack(fill="x", padx=2, pady=4)

    ttk.Label(options, text="axis_dist (Å)").grid(row=0, column=0, sticky="w")
    ttk.Entry(options, textvariable=axis_dist_var, width=12).grid(row=0, column=1, sticky="w", padx=4)
    make_help_button(options, "axis_dist", "axis_dist").grid(row=0, column=2, sticky="w")

    ttk.Label(options, text="axis_parallel").grid(row=0, column=3, sticky="w", padx=(12, 0))
    ttk.Combobox(
        options,
        textvariable=axis_parallel_var,
        values=("y", "n"),
        width=6,
        state="readonly",
    ).grid(row=0, column=4, sticky="w", padx=4)
    make_help_button(options, "axis_parallel", "axis_parallel").grid(row=0, column=5, sticky="w")

    ttk.Label(options, text="fix chain ID").grid(row=0, column=6, sticky="w", padx=(12, 0))
    ttk.Entry(options, textvariable=fix_chain_var, width=8).grid(row=0, column=7, sticky="w", padx=4)
    make_help_button(options, "fix chain", "fix").grid(row=0, column=8, sticky="w")

    ttk.Checkbutton(options, text="replicate", variable=replicate_var).grid(
        row=0, column=9, sticky="w", padx=(12, 0)
    )
    make_help_button(options, "replicate", "replicate").grid(row=0, column=10, sticky="w")

    ttk.Checkbutton(options, text="RE only", variable=re_only_var).grid(
        row=0, column=11, sticky="w", padx=(12, 0)
    )
    make_help_button(options, "RE only", "re_only").grid(row=0, column=12, sticky="w")

    ttk.Label(options, text="cir_shift (nt)").grid(row=0, column=13, sticky="w", padx=(12, 0))
    ttk.Entry(options, textvariable=cir_shift_var, width=8).grid(row=0, column=14, sticky="w", padx=4)
    make_help_button(options, "cir_shift", "cir_shift").grid(row=0, column=15, sticky="w")

    ttk.Label(options, text="3'-3' P resname").grid(row=1, column=0, sticky="w", pady=(6, 0))
    ttk.Combobox(
        options,
        textvariable=linker_phosphate_resname_var,
        values=("X33", "DA"),
        width=8,
        state="normal",
    ).grid(row=1, column=1, sticky="w", padx=4, pady=(6, 0))
    make_help_button(options, "3'-3' P resname", "linker_phosphate_resname").grid(
        row=1,
        column=2,
        sticky="w",
        pady=(6, 0),
    )

    pairs_box = ttk.LabelFrame(outer, text="Exchange pairs", padding=10, style="Section.TLabelframe")
    pairs_box.pack(fill="x", padx=2, pady=4)
    pairs_header = ttk.Frame(pairs_box)
    pairs_header.pack(fill="x")
    ttk.Label(pairs_header, text="Rows").pack(side="left")
    tk.Spinbox(
        pairs_header,
        from_=1,
        to=30,
        width=5,
        textvariable=pair_count_var,
        validate="key",
        validatecommand=validate_count_cmd,
        command=lambda: schedule_pair_rows(),
    ).pack(side="left", padx=6)
    make_help_button(pairs_header, "Exchange pairs", "pairs").pack(side="left")
    ttk.Label(
        pairs_header,
        text="pos1 / pos2 use residue+chain; optional beta angle is in degrees; kind = d / s / b.",
    ).pack(side="left", padx=8)
    pair_rows_frame = ttk.Frame(pairs_box)
    pair_rows_frame.pack(fill="x", pady=(8, 0))
    ttk.Label(pair_rows_frame, text="").grid(row=0, column=0, sticky="w")
    ttk.Label(pair_rows_frame, text="pos1 (nt+chain)").grid(row=0, column=1, sticky="w", padx=4)
    ttk.Label(pair_rows_frame, text="pos2 (nt+chain)").grid(row=0, column=2, sticky="w", padx=4)
    ttk.Label(pair_rows_frame, text="beta angle (deg, optional)").grid(row=0, column=3, sticky="w", padx=4)
    ttk.Label(pair_rows_frame, text="kind").grid(row=0, column=4, sticky="w", padx=4)

    pair_args_frame = ttk.Frame(pairs_box)
    pair_args_frame.pack(fill="x", pady=(8, 0))
    ttk.Label(pair_args_frame, text="CLI pair args").pack(side="left")
    ttk.Entry(pair_args_frame, textvariable=pair_args_var, width=100).pack(
        side="left",
        fill="x",
        expand=True,
        padx=4,
    )
    make_help_button(pair_args_frame, "CLI pair args", "pair_args").pack(side="left")

    axis_box = ttk.LabelFrame(outer, text="Axis definitions", padding=10, style="Section.TLabelframe")
    axis_box.pack(fill="x", padx=2, pady=4)
    axis_header = ttk.Frame(axis_box)
    axis_header.pack(fill="x")
    ttk.Label(axis_header, text="Rows").pack(side="left")
    axis_count_spinbox = tk.Spinbox(
        axis_header,
        from_=0,
        to=30,
        width=5,
        textvariable=axis_count_var,
        validate="key",
        validatecommand=validate_count_cmd,
        command=lambda: schedule_axis_rows(),
    )
    axis_count_spinbox.pack(side="left", padx=6)
    make_help_button(axis_header, "Axis definitions", "axis_range").pack(side="left")
    ttk.Label(axis_header, text="Axis: A,B or B26-B60,A1-A35; move: C,D or C1-C50,D").pack(side="left", padx=8)

    user_axis_frame = ttk.Frame(axis_box)
    user_axis_frame.pack(fill="x", pady=(8, 0))
    user_axis_check = ttk.Checkbutton(
        user_axis_frame,
        text="Use direction + point axis",
        variable=user_axis_var,
        command=lambda: refresh_user_axis_state(),
    )
    user_axis_check.pack(side="left")
    make_help_button(user_axis_frame, "Direction + point axis", "user_axis").pack(side="left", padx=(4, 8))
    ttk.Label(user_axis_frame, text="dir").pack(side="left")
    for idx, label_text in enumerate(("x", "y", "z")):
        ttk.Label(user_axis_frame, text=label_text).pack(side="left", padx=(6 if idx == 0 else 2, 0))
        entry = ttk.Entry(user_axis_frame, textvariable=user_axis_dir_vars[idx], width=9)
        entry.pack(side="left", padx=(2, 0))
        user_axis_input_widgets.append(entry)
    ttk.Label(user_axis_frame, text="point").pack(side="left", padx=(12, 0))
    for idx, label_text in enumerate(("x", "y", "z")):
        ttk.Label(user_axis_frame, text=label_text).pack(side="left", padx=(6 if idx == 0 else 2, 0))
        entry = ttk.Entry(user_axis_frame, textvariable=user_axis_point_vars[idx], width=9)
        entry.pack(side="left", padx=(2, 0))
        user_axis_input_widgets.append(entry)

    axis_rows_frame = ttk.Frame(axis_box)
    axis_rows_frame.pack(fill="x", pady=(8, 0))
    ttk.Label(axis_rows_frame, text="").grid(row=0, column=0, sticky="w")
    ttk.Label(axis_rows_frame, text="axis definition").grid(row=0, column=1, sticky="w", padx=4)
    ttk.Label(axis_rows_frame, text="move with axis").grid(row=0, column=2, sticky="w", padx=4)

    other_tools_box = ttk.LabelFrame(outer, text="Other tools", padding=10, style="Section.TLabelframe")
    other_tools_box.pack(fill="x", padx=2, pady=4)
    bend_helix_button = ttk.Button(other_tools_box, text="Bend Helix")
    bend_helix_button.pack(side="left")
    do_symmetry_button = ttk.Button(other_tools_box, text="Do Symmetry")
    do_symmetry_button.pack(side="left", padx=(6, 0))
    add_pdb_link_button = ttk.Button(other_tools_box, text="Add PDB LINK Record")
    add_pdb_link_button.pack(side="left", padx=(6, 0))
    insert_virtual_resi_button = ttk.Button(other_tools_box, text="Insert Virtual Resi")
    insert_virtual_resi_button.pack(side="left", padx=(6, 0))
    generate_lattice_button = ttk.Button(other_tools_box, text="Generate Lattice")
    generate_lattice_button.pack(side="left", padx=(6, 0))
    get_phenix_restraints_button = ttk.Button(other_tools_box, text="Get Phenix Restraints")
    get_phenix_restraints_button.pack(side="left", padx=(6, 0))

    buttons = ttk.Frame(outer)
    buttons.pack(fill="x", padx=2, pady=4)
    run_button = ttk.Button(buttons, text="Run")
    run_button.pack(side="left")
    ttk.Button(buttons, text="Close", command=root.destroy).pack(side="left", padx=6)

    log_box = ttk.LabelFrame(outer, text="Run log", padding=10, style="Section.TLabelframe")
    log_box.pack(fill="both", expand=True, padx=2, pady=4)
    log_widget = scrolledtext.ScrolledText(log_box, wrap="word", height=24)
    log_widget.pack(fill="both", expand=True)

    def append_log(message: str) -> None:
        log_widget.insert("end", message)
        log_widget.see("end")
        log_widget.update_idletasks()

    def launch_bundled_gui_tool(tool_label: str, script_name: str, extra_args: list) -> None:
        tool_script = Path(__file__).resolve().parent / "re_helix_lib" / script_name
        if not tool_script.exists():
            messagebox.showerror(
                APP_TITLE,
                f"Could not find {tool_label} tool at:\n{tool_script}",
            )
            return

        cmd = [sys.executable, str(tool_script), "--gui"]
        cmd.extend(extra_args)

        append_log(f"[GUI] Launching {tool_label} tool:\n")
        append_log("    " + " ".join(shlex.quote(part) for part in cmd) + "\n\n")

        try:
            proc = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
            )
        except OSError as exc:
            messagebox.showerror(APP_TITLE, f"Failed to launch {tool_label} tool:\n{exc}")
            append_log(f"[GUI] Failed to launch {tool_label} tool: {exc}\n")
            return

        append_log(f"[GUI] {tool_label} process started; stdout/stderr will be mirrored here.\n\n")

        def worker() -> None:
            assert proc.stdout is not None
            q = process_state["queue"]
            for line in proc.stdout:
                q.put(f"[{tool_label}] {line}")
            proc.stdout.close()
            returncode = proc.wait()
            q.put(f"[{tool_label}] Process exited with code {returncode}.\n")

        threading.Thread(target=worker, daemon=True).start()

    def launch_bend_helix_tool() -> None:
        current_pdb = pdb_var.get().strip()
        extra_args = []
        if current_pdb:
            extra_args.extend(["--input", current_pdb])
        launch_bundled_gui_tool("Bend Helix", "bend_helix.py", extra_args)

    def launch_do_symmetry_tool() -> None:
        current_pdb = pdb_var.get().strip()
        extra_args = [current_pdb] if current_pdb else []
        launch_bundled_gui_tool("Do Symmetry", "do_symmetry.py", extra_args)

    def launch_add_pdb_link_record_tool() -> None:
        current_pdb = pdb_var.get().strip()
        extra_args = [current_pdb] if current_pdb else []
        launch_bundled_gui_tool("Add PDB LINK Record", "add_pdb_link_record.py", extra_args)

    def launch_insert_virtual_resi_tool() -> None:
        current_pdb = pdb_var.get().strip()
        extra_args = [current_pdb] if current_pdb else []
        launch_bundled_gui_tool("Insert Virtual Resi", "insert_virtual_resi.py", extra_args)

    def launch_generate_lattice_tool() -> None:
        current_pdb = pdb_var.get().strip()
        extra_args = [current_pdb] if current_pdb else []
        launch_bundled_gui_tool("Generate Lattice", "generate_lattice.py", extra_args)

    def launch_get_phenix_restraints_tool() -> None:
        current_pdb = pdb_var.get().strip()
        extra_args = [current_pdb] if current_pdb else []
        launch_bundled_gui_tool("Get Phenix Restraints", "get_phenix_restraints.py", extra_args)

    row_targets = {"pair": 3, "axis": 0}
    render_state = {"pair_pending": False, "axis_pending": False}

    def current_pair_target() -> int:
        return int(row_targets["pair"])

    def current_axis_target() -> int:
        return int(row_targets["axis"])

    def refresh_user_axis_state() -> None:
        use_user_axis = bool(user_axis_var.get())
        user_entry_state = "normal" if use_user_axis else "disabled"
        range_entry_state = "disabled" if use_user_axis else "normal"
        axis_count_spinbox.configure(state=range_entry_state)
        for widget in user_axis_input_widgets:
            widget.configure(state=user_entry_state)
        for item in axis_widgets:
            _label, _axis_var, _move_var, axis_entry, move_entry = item
            axis_entry.configure(state=range_entry_state)
            move_entry.configure(state=range_entry_state)
        schedule_scrollbar_refresh()

    def render_pair_rows() -> None:
        target_opt = _coerce_gui_count(pair_count_var.get(), 1, 30)
        if target_opt is None:
            return
        target = target_opt
        row_targets["pair"] = target
        while len(pair_widgets) < target:
            row_index = len(pair_widgets)
            pos1_var = tk.StringVar()
            pos2_var = tk.StringVar()
            beta_var = tk.StringVar()
            kind_var = tk.StringVar(value="d")
            label = ttk.Label(pair_rows_frame, text=f"Pair {row_index + 1}")
            e1 = ttk.Entry(pair_rows_frame, textvariable=pos1_var, width=16)
            e2 = ttk.Entry(pair_rows_frame, textvariable=pos2_var, width=16)
            e_beta = ttk.Entry(pair_rows_frame, textvariable=beta_var, width=14)
            kind_box = ttk.Combobox(
                pair_rows_frame,
                textvariable=kind_var,
                values=("d", "s", "b"),
                width=6,
                state="readonly",
            )
            pair_widgets.append(
                {
                    "label_widget": label,
                    "pos1_var": pos1_var,
                    "pos2_var": pos2_var,
                    "beta_var": beta_var,
                    "kind_var": kind_var,
                    "widgets": [e1, e2, e_beta, kind_box],
                }
            )
        for idx, item in enumerate(pair_widgets):
            visible = idx < target
            label_widget = item["label_widget"]
            grid_row = idx + 1
            if visible:
                label_widget.configure(text=f"Pair {idx + 1}")
                label_widget.grid(row=grid_row, column=0, sticky="w", pady=2)
                for col, widget in enumerate(item["widgets"], start=1):
                    widget.grid(row=grid_row, column=col, sticky="w", padx=4, pady=2)
            else:
                label_widget.grid_remove()
                for widget in item["widgets"]:
                    widget.grid_remove()
        schedule_scrollbar_refresh()

    def render_axis_rows() -> None:
        target_opt = _coerce_gui_count(axis_count_var.get(), 0, 30)
        if target_opt is None:
            return
        target = target_opt
        row_targets["axis"] = target
        while len(axis_widgets) < target:
            row_index = len(axis_widgets)
            axis_var = tk.StringVar()
            move_var = tk.StringVar()
            label = ttk.Label(axis_rows_frame, text=f"Axis {row_index + 1}")
            axis_entry = ttk.Entry(axis_rows_frame, textvariable=axis_var, width=42)
            move_entry = ttk.Entry(axis_rows_frame, textvariable=move_var, width=42)
            axis_widgets.append((label, axis_var, move_var, axis_entry, move_entry))
        for idx, item in enumerate(axis_widgets):
            label, axis_var, move_var, axis_entry, move_entry = item
            grid_row = idx + 1
            if idx < target:
                label.configure(text=f"Axis {idx + 1}")
                label.grid(row=grid_row, column=0, sticky="w", pady=2)
                axis_entry.grid(row=grid_row, column=1, sticky="w", padx=4, pady=2)
                move_entry.grid(row=grid_row, column=2, sticky="w", padx=4, pady=2)
            else:
                label.grid_remove()
                axis_entry.grid_remove()
                move_entry.grid_remove()
        refresh_user_axis_state()
        schedule_scrollbar_refresh()

    def schedule_pair_rows(*_args) -> None:
        if render_state["pair_pending"]:
            return
        render_state["pair_pending"] = True

        def _run() -> None:
            render_state["pair_pending"] = False
            render_pair_rows()

        root.after_idle(_run)

    def schedule_axis_rows(*_args) -> None:
        if render_state["axis_pending"]:
            return
        render_state["axis_pending"] = True

        def _run() -> None:
            render_state["axis_pending"] = False
            render_axis_rows()

        root.after_idle(_run)

    def poll_queue() -> None:
        q = process_state["queue"]
        while True:
            try:
                item = q.get_nowait()
            except queue.Empty:
                break
            if isinstance(item, tuple) and item and item[0] == "done":
                returncode = int(item[1])
                if returncode == 0:
                    append_log("\n[GUI] Run finished successfully (exit code 0).\n")
                else:
                    append_log(
                        f"\n[GUI] Run finished with a non-zero exit code {returncode}.\n"
                    )
                run_button.state(["!disabled"])
                process_state["proc"] = None
            else:
                append_log(str(item))
        root.after(100, poll_queue)

    def start_run() -> None:
        if process_state.get("proc") is not None:
            messagebox.showinfo(APP_TITLE, "A run is already in progress.")
            return

        try:
            pair_rows: List[Dict[str, str]] = []
            for item in pair_widgets[: current_pair_target()]:
                pair_rows.append(
                    {
                        "pos1": item["pos1_var"].get(),
                        "pos2": item["pos2_var"].get(),
                        "beta_angle": item["beta_var"].get(),
                        "kind": item["kind_var"].get(),
                    }
                )
            axis_rows = [
                axis_var.get()
                for _label, axis_var, _move_var, _axis_entry, _move_entry
                in axis_widgets[: current_axis_target()]
            ]
            axis_move_rows = [
                move_var.get()
                for _label, _axis_var, move_var, _axis_entry, _move_entry
                in axis_widgets[: current_axis_target()]
            ]
            cmd = _build_equivalent_cli_command(
                str(Path(__file__).resolve()),
                pdb_var.get().strip(),
                helix_defs_var.get(),
                pair_rows,
                pair_args_var.get(),
                axis_rows,
                axis_move_rows,
                bool(user_axis_var.get()),
                [var.get() for var in user_axis_dir_vars],
                [var.get() for var in user_axis_point_vars],
                output_var.get(),
                axis_dist_var.get(),
                axis_parallel_var.get(),
                fix_chain_var.get(),
                bool(replicate_var.get()),
                bool(re_only_var.get()),
                cir_shift_var.get(),
                linker_phosphate_resname_var.get(),
            )
        except Exception as exc:
            append_log(f"[GUI] {exc}\n")
            return

        append_log("[GUI] Equivalent CLI command:\n")
        append_log("    " + " ".join(shlex.quote(part) for part in cmd) + "\n\n")

        try:
            proc = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
            )
        except OSError as exc:
            append_log(f"[GUI] Failed to start process: {exc}\n")
            return

        process_state["proc"] = proc
        run_button.state(["disabled"])
        q = process_state["queue"]

        def worker() -> None:
            assert proc.stdout is not None
            for line in proc.stdout:
                q.put(line)
            proc.stdout.close()
            returncode = proc.wait()
            q.put(("done", returncode))

        threading.Thread(target=worker, daemon=True).start()

    pair_count_var.trace_add("write", schedule_pair_rows)
    axis_count_var.trace_add("write", schedule_axis_rows)

    run_button.configure(command=start_run)
    bend_helix_button.configure(command=launch_bend_helix_tool)
    do_symmetry_button.configure(command=launch_do_symmetry_tool)
    add_pdb_link_button.configure(command=launch_add_pdb_link_record_tool)
    insert_virtual_resi_button.configure(command=launch_insert_virtual_resi_tool)
    generate_lattice_button.configure(command=launch_generate_lattice_tool)
    get_phenix_restraints_button.configure(command=launch_get_phenix_restraints_tool)
    render_pair_rows()
    render_axis_rows()
    schedule_scrollbar_refresh()
    root.after(100, poll_queue)
    root.mainloop()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    if len(sys.argv) == 1 or "--gui" in sys.argv:
        try:
            _launch_gui()
        except Exception as exc:
            print(f"Error launching GUI: {exc}", file=sys.stderr)
            sys.exit(1)
        return

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description=(
            "Align nucleic-acid helices (or multi-chain helix groups) in a PDB "
            "based on cross-helix P-atom pairs from a reciprocal-exchange-style "
            "specification, and then apply the reciprocal exchanges. Use "
            "--re_only to apply reciprocal exchanges without alignment."
        )
    )
    parser.add_argument(
        "pdb_in",
        help="Input PDB file.",
    )
    parser.add_argument(
        "ops",
        nargs="+",
        help=(
            "Tokens including optional helix definitions and exchange specs. "
            "Examples:\n"
            "  (AB) (CD) 30A 8D d 13B 24C s\n"
            "  (ABMN) 30A 8D d 13B 24C s\n"
            "  30A 8D d 13B 24C s\n"
            "  26A 9C 90 d     # fixed beta angle for one single-site helix pair\n"
            "Helix defs like (AB) or (ABMN) mean all chains inside the "
            "parentheses form one helix group that moves as a rigid block. "
            "If any are provided, helices are not auto-detected. "
            "Each exchange can also be written as <pos1> <pos2> <beta_deg> <kind>; "
            "the optional beta angle is used only for exactly one site between a "
            "helix pair under --axis_parallel n. Legacy docs called this rho_deg."
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="pdb_out_base",
        default=None,
        help=(
            "Base name for output files (extension optional). "
            "Normal outputs are <base>_aligned.pdb and <base>_aligned_rex.pdb; "
            "RE-only output is <base>_rex.pdb. "
            "Default base: input filename without extension."
        ),
    )
    parser.add_argument(
        "--gui",
        action="store_true",
        help=(
            "Launch the Tk GUI instead of running directly from the command line. "
            "If no command-line arguments are given, GUI mode is entered automatically."
        ),
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"{SOFTWARE_NAME} {SOFTWARE_VERSION}",
        help="Show the app version and exit.",
    )
    parser.add_argument(
        "--axis_range",
        action="append",
        default=[],
        help=(
            "Chain/range definition for helical-axis estimation. Repeat this option as needed. "
            "Examples: --axis_range B26-B60,A1-A35 or --axis_range A,B. Whole-chain letters "
            "use all P atoms on that chain for the axis. If paired with --axis_move, the same "
            "row defines an axis-coupled helix group and the listed move payload follows that axis. "
            "In replication mode, base-template axis groups can be propagated to copies, and final "
            "post-replication chain IDs can still be targeted explicitly."
        ),
    )
    parser.add_argument(
        "--axis_move",
        action="append",
        default=[],
        help=(
            "Additional chains or residue windows to move with the corresponding --axis_range row. "
            "Examples: --axis_move C,D or --axis_move C1-C50,D. This avoids triplex stdin prompts "
            "by letting an axis row such as --axis_range A,B move payload chain C with that axis."
        ),
    )
    parser.add_argument(
        "--user_axis_dir",
        "--user-axis-dir",
        nargs=3,
        type=float,
        metavar=("X", "Y", "Z"),
        default=None,
        help=(
            "Direction vector for a user-defined alignment axis. Must be supplied "
            "together with --user_axis_point. When used, fixed/moving helix-axis "
            "estimation is skipped and each moving helix is optimized by rotation "
            "around this axis plus a full XYZ translation."
        ),
    )
    parser.add_argument(
        "--user_axis_point",
        "--user-axis-point",
        nargs=3,
        type=float,
        metavar=("X", "Y", "Z"),
        default=None,
        help=(
            "Point on the user-defined alignment axis. Must be supplied together "
            "with --user_axis_dir."
        ),
    )
    parser.add_argument(
        "--axis_dist",
        type=float,
        default=22.0,
        help="Target distance between helical axes in Å (default: 22.0).",
    )
    parser.add_argument(
        "--axis_parallel",
        choices=["y", "Y", "n", "N"],
        default="y",
        help=(
            "If 'y' (default), helices are kept parallel and we optimise "
            "d (axial slide), tau (axial twist/spin of the moving helix) and "
            "phi (orbital azimuth around the fixed helix). "
            "If 'n', we additionally allow a tilt of helix 2 around a "
            "P-pair–dependent common perpendicular (beta), so the axes are "
            "no longer parallel. With a single reciprocal-exchange site between "
            "two helices, an operation can be written as <pos1> <pos2> "
            "<beta_deg> <kind> to hold beta at that angle while optimizing "
            "d/tau/phi. Legacy <rho_deg> examples are accepted as beta. "
            "Fixed beta definitions are ignored with a warning for --axis_parallel y "
            "or multi-site helix pairs."
        ),
    )
    parser.add_argument(
        "--fix",
        dest="fix_chain",
        default=None,
        help=(
            "Chain ID whose helix group should remain fixed during alignment. "
            "For example, '--fix A' keeps the helix group containing chain A "
            "as the immobile reference within its connected component."
        ),
    )
    parser.add_argument(
        "--replicate",
        action="store_true",
        help=(
            "Replicate the entire set of chains. All original chains are first "
            "renamed to consecutive letters A, B, C, ... (sorted by original "
            "chain ID), and then additional full copies of this chain set are "
            "created whenever those chain IDs are needed by the exchange "
            "specifications (e.g. to create copies with chains C/D, E/F, ...). "
            "If the input PDB appears to contain exactly one helix component, "
            "replicate mode is also enabled automatically even without this flag."
        ),
    )
    parser.add_argument(
        "--re_only",
        "--re-only",
        dest="re_only",
        action="store_true",
        help=(
            "Apply reciprocal exchanges directly to the input PDB without "
            "running helix alignment. Writes <base>_rex.pdb. If --replicate is "
            "also supplied, replication is performed first; otherwise RE-only "
            "mode behaves like reciprocal_exchange_pdb on the input structure."
        ),
    )
    parser.add_argument(
        "--re",
        action="store_true",
        help=(
            "Legacy flag: previously controlled whether reciprocal exchanges "
            "were applied. It is now ignored; exchanges are always applied and "
            "an *_aligned_rex.pdb file is always written."
        ),
    )
    parser.add_argument(
        "--cir_shift",
        type=int,
        default=8,
        help=(
            "Residue shift for circular strands when applying reciprocal "
            "exchanges (default: 8)."
        ),
    )
    parser.add_argument(
        "--linker_phosphate_resname",
        "--linker-phosphate-resname",
        default=None,
        help=(
            "Residue name for phosphate-only residues inserted at bowtie 3'-3' "
            "linkages. Default: X33 as HETATM. Custom 1-3 character names are "
            "written as HETATM by default. DA/dA writes regular ATOM DA unless "
            "--linker_phosphate_record is supplied."
        ),
    )
    parser.add_argument(
        "--linker_phosphate_record",
        "--linker-phosphate-record",
        choices=["ATOM", "HETATM", "atom", "hetatm"],
        default=None,
        help=(
            "Advanced: record type for inserted 3'-3' linker phosphates. "
            "Default is HETATM, except DA/dA defaults to ATOM."
        ),
    )

    args = parser.parse_args()
    try:
        linker_phosphate_style = rex.make_linker_phosphate_style(
            args.linker_phosphate_resname,
            args.linker_phosphate_record,
        )
    except ValueError as exc:
        parser.error(str(exc))

    reset_helix_axis_overrides()

    user_axis: Optional[UserAxisDefinition] = None
    if (args.user_axis_dir is None) != (args.user_axis_point is None):
        parser.error("--user_axis_dir and --user_axis_point must be supplied together.")
    if args.user_axis_dir is not None and args.user_axis_point is not None:
        try:
            user_axis = normalize_user_axis_definition(args.user_axis_dir, args.user_axis_point)
        except ValueError as exc:
            parser.error(str(exc))

    if args.re_only or user_axis is not None:
        axis_range_defs_input = []
        axis_move_defs_input = []
    else:
        try:
            axis_range_defs_input = parse_axis_range_specs(args.axis_range)
            axis_move_defs_input = parse_axis_move_specs(args.axis_move)
            pair_axis_move_definitions(axis_range_defs_input, axis_move_defs_input)
        except ValueError as exc:
            print(f"Error parsing --axis_range / --axis_move definitions: {exc}", file=sys.stderr)
            sys.exit(1)

    # Determine base name for outputs
    if args.pdb_out_base is None:
        if args.pdb_in.lower().endswith(".pdb"):
            base = args.pdb_in[:-4]
        else:
            base = args.pdb_in
    else:
        base = args.pdb_out_base
        if base.lower().endswith(".pdb"):
            base = base[:-4]

    pdb_out_aligned = base + "_aligned.pdb"
    pdb_out_aligned_rex = base + "_aligned_rex.pdb"
    pdb_out_rex_only = base + "_rex.pdb"

    # Flatten ops tokens (support quoted strings)
    if len(args.ops) == 1 and (" " in args.ops[0] or "\t" in args.ops[0]):
        tokens = args.ops[0].split()
    else:
        tokens = args.ops

    # Extract explicit helix defs like (AB), (CD), (ABMN) and leave the rest
    helix_defs, exch_tokens = parse_helix_definition_tokens(tokens)

    # Parse exchange specs
    try:
        exchange_specs = parse_exchange_specs(exch_tokens)
    except ValueError as exc:
        print(f"Error parsing exchange specifications: {exc}", file=sys.stderr)
        sys.exit(1)

    # Read PDB
    rec_list: List[pdb_atom_record] = []
    try:
        with open(args.pdb_in, "r") as fin:
            file2rec(fin, rec_list)
    except OSError as exc:
        print(f"Error reading PDB file '{args.pdb_in}': {exc}", file=sys.stderr)
        sys.exit(1)

    if not rec_list:
        print("No ATOM/HETATM/TER records found in input PDB.", file=sys.stderr)
        sys.exit(1)

    command_text = " ".join(shlex.quote(arg) for arg in [sys.executable] + sys.argv)

    if args.re_only and (args.axis_range or args.axis_move):
        sys.stderr.write("[re_helix] Warning: --axis_range/--axis_move is ignored in RE-only mode.\n")
    if user_axis is not None and (args.axis_range or args.axis_move):
        sys.stderr.write(
            "[re_helix] Warning: --axis_range/--axis_move is ignored when "
            "--user_axis_dir/--user_axis_point is supplied.\n"
        )
    if args.re_only and user_axis is not None:
        sys.stderr.write("[re_helix] Warning: user-defined alignment axis is ignored in RE-only mode.\n")

    if args.re_only and not args.replicate:
        try:
            write_reciprocal_exchange_output(
                rec_list,
                exchange_specs,
                cir_shift=args.cir_shift,
                command_text=command_text,
                output_path=pdb_out_rex_only,
                output_stage="reciprocal_exchange_only",
                linker_phosphate_style=linker_phosphate_style,
            )
        except OSError as exc:
            print(f"Error writing RE-only PDB '{pdb_out_rex_only}': {exc}", file=sys.stderr)
            sys.exit(1)
        except Exception as exc:
            print(f"Error during RE-only reciprocal exchange stage: {exc}", file=sys.stderr)
            sys.exit(1)

        print(f"Wrote RE-only PDB to '{pdb_out_rex_only}'.")
        return

    # Minimal maps for replication decision
    chain_to_res_atoms0, residue_to_P_atom0, chain_to_P_atoms0 = build_nucleic_acid_maps(rec_list)
    if not chain_to_P_atoms0:
        print("No P atoms found in input PDB; alignment cannot proceed.", file=sys.stderr)
        sys.exit(1)
    original_atom_chains0 = sorted(
        {atom.chainID for atom in rec_list if atom.recordName in ("ATOM", "HETATM")}
    )

    original_chain_lookup0 = {ch.upper(): ch for ch in original_atom_chains0}
    template_axis_defs_input: List[AxisSelection] = []
    template_move_defs_input: List[MoveSelection] = []
    for axis_def, move_def in pair_axis_move_definitions(axis_range_defs_input, axis_move_defs_input):
        row_chains = _selection_chains(axis_def) | _selection_chains(move_def)
        if row_chains and all(ch.upper() in original_chain_lookup0 for ch in row_chains):
            template_axis_defs_input.append(axis_def)
            template_move_defs_input.append(move_def)

    template_axis_defs_resolved: List[ResolvedAxisRange] = []
    template_move_defs_resolved: List[MoveSelection] = []
    axis_coupled_helices0: List[HelixID] = []
    if template_axis_defs_input:
        try:
            template_axis_defs_resolved = resolve_axis_range_definitions(
                template_axis_defs_input,
                chain_to_P_atoms0.keys(),
                chain_to_P_atoms0,
            )
            template_move_defs_resolved = resolve_axis_move_definitions(
                template_move_defs_input,
                original_atom_chains0,
            )
            axis_coupled_helices0 = build_axis_coupled_helix_defs(
                template_axis_defs_resolved,
                template_move_defs_resolved,
                helix_defs if helix_defs else None,
            )
            register_axis_coupling_overrides(
                template_axis_defs_resolved,
                template_move_defs_resolved,
                helix_defs if helix_defs else None,
            )
        except ValueError as exc:
            print(f"Error resolving base-template --axis_range / --axis_move definitions: {exc}", file=sys.stderr)
            sys.exit(1)

    # Determine how many helices the input PDB appears to contain
    if helix_defs:
        helices0 = sorted({tuple(sorted(h)) for h in helix_defs})
    elif axis_coupled_helices0 and set(chain_to_P_atoms0.keys()).issubset(
        set().union(*(set(h) for h in axis_coupled_helices0))
    ):
        helices0 = sorted(set(axis_coupled_helices0))
    else:
        chain_to_helix0 = compute_chain_partner_map(chain_to_P_atoms0)
        helices0 = sorted({chain_to_helix0[ch] for ch in chain_to_helix0})

    # Keep the old "auto-single-helix" behavior if you want.
    # If you prefer to replicate ONLY when --replicate is given, set auto_single_helix = False.
    auto_single_helix = (len(helices0) == 1)
    replicate_active = args.replicate or auto_single_helix

    if replicate_active:
        try:
            exchange_specs, helix_defs_repl = replicate_all_chains(
                rec_list,
                helices0,
                exchange_specs,
                helix_defs if helix_defs else (axis_coupled_helices0 if axis_coupled_helices0 else None),
            )
        except Exception as exc:
            print(f"Error during helix replication: {exc}", file=sys.stderr)
            sys.exit(1)
        helix_defs_for_align: Optional[List[HelixID]] = helix_defs_repl
    else:
        helix_defs_for_align = helix_defs if helix_defs else (axis_coupled_helices0 if axis_coupled_helices0 else None)

    if args.re_only:
        try:
            write_reciprocal_exchange_output(
                rec_list,
                exchange_specs,
                cir_shift=args.cir_shift,
                command_text=command_text,
                output_path=pdb_out_rex_only,
                output_stage="reciprocal_exchange_only",
                linker_phosphate_style=linker_phosphate_style,
            )
        except OSError as exc:
            print(f"Error writing RE-only PDB '{pdb_out_rex_only}': {exc}", file=sys.stderr)
            sys.exit(1)
        except Exception as exc:
            print(f"Error during RE-only reciprocal exchange stage: {exc}", file=sys.stderr)
            sys.exit(1)

        print(f"Wrote RE-only PDB to '{pdb_out_rex_only}'.")
        return

    _chain_to_res_atoms_curr, _residue_to_P_atom_curr, chain_to_P_atoms_curr = build_nucleic_acid_maps(rec_list)
    atom_chains_curr = sorted(
        {atom.chainID for atom in rec_list if atom.recordName in ("ATOM", "HETATM")}
    )
    if replicate_active:
        axis_range_defs_for_resolve = translate_axis_range_definitions_for_replication(
            axis_range_defs_input,
            chain_to_P_atoms0.keys(),
            chain_to_P_atoms_curr.keys(),
        )
        axis_move_defs_for_resolve = translate_axis_range_definitions_for_replication(
            axis_move_defs_input,
            original_atom_chains0,
            atom_chains_curr,
        )
    else:
        axis_range_defs_for_resolve = axis_range_defs_input
        axis_move_defs_for_resolve = axis_move_defs_input

    try:
        axis_range_defs_resolved = resolve_axis_range_definitions(
            axis_range_defs_for_resolve,
            chain_to_P_atoms_curr.keys(),
            chain_to_P_atoms_curr,
        )
        axis_move_defs_resolved = resolve_axis_move_definitions(
            axis_move_defs_for_resolve,
            atom_chains_curr,
        )
        pair_axis_move_definitions(axis_range_defs_resolved, axis_move_defs_resolved)
    except ValueError as exc:
        print(f"Error resolving --axis_range / --axis_move definitions: {exc}", file=sys.stderr)
        sys.exit(1)
    set_helix_axis_range_definitions(axis_range_defs_resolved)
    set_helix_axis_move_definitions(axis_move_defs_resolved)

    axis_parallel_flag = args.axis_parallel.lower() == "y"

    # 1) Align helices
    try:
        align_helices_for_exchanges(
            rec_list,
            exchange_specs,
            axis_dist=args.axis_dist,
            axis_parallel_flag=axis_parallel_flag,
            explicit_helices=helix_defs_for_align,
            fix_chain=args.fix_chain,
            user_axis=user_axis,
        )
    except Exception as exc:
        print(f"Error during helix alignment: {exc}", file=sys.stderr)
        sys.exit(1)

    command_text = " ".join(shlex.quote(arg) for arg in [sys.executable] + sys.argv)

    # 2) Write aligned PDB (pre-RE)
    try:
        pre_re_header = rex.build_re_script_header_lines(
            software_name=SOFTWARE_NAME,
            software_version=SOFTWARE_VERSION,
            developer=SOFTWARE_DEVELOPER,
            command=command_text,
            output_stage="aligned_pre_re",
            atom_rec_list=rec_list,
            extra_special_events=[
                f"reciprocal_exchange_pending count={len(exchange_specs)} output={pdb_out_aligned_rex}"
            ],
        )
        dependency_line = _rex_dependency_remark_line()
        if dependency_line is not None:
            pre_re_header.append(dependency_line)
        with open(pdb_out_aligned, "w") as fout:
            rex.write_pdb_with_header(rec_list, [], fout, header_lines=pre_re_header, reorder_serial=True)
    except OSError as exc:
        print(f"Error writing aligned PDB '{pdb_out_aligned}': {exc}", file=sys.stderr)
        sys.exit(1)
    except Exception as exc:
        print(f"Error building header for aligned PDB '{pdb_out_aligned}': {exc}", file=sys.stderr)
        sys.exit(1)

    print(f"Wrote aligned PDB (pre-RE) to '{pdb_out_aligned}'.")

    # 3) Always apply reciprocal exchanges to generate *_aligned_rex.pdb
    try:
        write_reciprocal_exchange_output(
            rec_list,
            exchange_specs,
            cir_shift=args.cir_shift,
            command_text=command_text,
            output_path=pdb_out_aligned_rex,
            output_stage="aligned_rex",
            linker_phosphate_style=linker_phosphate_style,
        )
    except OSError as exc:
        print(f"Error writing aligned+RE PDB '{pdb_out_aligned_rex}': {exc}", file=sys.stderr)
        sys.exit(1)
    except Exception as exc:
        print(f"Error during aligned+RE reciprocal exchange stage: {exc}", file=sys.stderr)
        sys.exit(1)

    print(f"Wrote aligned+RE PDB to '{pdb_out_aligned_rex}'.")


if __name__ == "__main__":
    main()
