#!/usr/bin/env python3
"""Reverse selected nucleic-acid strand serializations while preserving topology.

This tool reverses complete PDB residue blocks and continuously renumbers them.
It does not move atoms or make a reverse complement. Supported backbone LINK
topology is reconstructed before the edit and regenerated afterward so the
same covalent atom pairs remain connected in the reversed serialization.
"""

from __future__ import annotations

import argparse
import shlex
import sys
from dataclasses import dataclass, field
from io import StringIO
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

try:
    from re_helix_lib import edit_pdb_link
    from re_helix_lib import reciprocal_exchange_pdbV3_3 as rex
    from re_helix_lib.gui_icon import apply_optional_icon
    from re_helix_lib.permute_chain import (
        build_chain_inventory_remarks,
        collect_residue_blocks,
        coordinate_identity,
        format_resseq,
        is_chain_inventory_remark,
        replace_remark_residue_references,
        split_line_ending,
        update_coordinate_residue,
        update_het_line,
        update_ter_line,
    )
except ImportError:  # pragma: no cover - direct script execution fallback
    import edit_pdb_link
    import reciprocal_exchange_pdbV3_3 as rex
    from gui_icon import apply_optional_icon
    from permute_chain import (
        build_chain_inventory_remarks,
        collect_residue_blocks,
        coordinate_identity,
        format_resseq,
        is_chain_inventory_remark,
        replace_remark_residue_references,
        split_line_ending,
        update_coordinate_residue,
        update_het_line,
        update_ter_line,
    )


TOOL_NAME = "Reverse Strand Direction"
VERSION = "1.0"


@dataclass
class ReverseStats:
    chain_id: str
    residue_count: int
    is_cycle: bool
    numbering_start: int
    numbering_end: int
    old_first_resseq: int
    old_last_resseq: int
    new_first_source_resseq: int
    new_last_source_resseq: int
    coordinate_lines_changed: int = 0
    ter_lines_changed: int = 0
    het_lines_changed: int = 0
    remark_references_changed: int = 0
    remark_inventory_rebuilt: bool = False
    warnings: List[str] = field(default_factory=list)


@dataclass
class ReverseResult:
    output_path: Path
    stats: List[ReverseStats]
    topology_link_count: int
    preserved_link_count: int
    total_link_count: int


@dataclass(frozen=True)
class _InputLink:
    """Parsed LINK fields paired with the exact source line."""

    record: edit_pdb_link.pdb_link_record
    raw_line: str


@dataclass
class _ScopedInput:
    """Topology and LINK inventory needed for the selected strands only."""

    nodes: List[rex.ResidueNode]
    label_to_idx: Dict[rex.Label, int]
    topology: rex.InputBackboneTopology
    input_links: List[edit_pdb_link.pdb_link_record]
    retained_links: List[_InputLink]


def default_output_path(input_path: Path) -> Path:
    return input_path.with_name(f"{input_path.stem}_strand_reversed{input_path.suffix}")


def _record_name(line: str) -> str:
    return line[:6].strip()


def _normalize_chain_selector(chain_id: object) -> str:
    """Normalize a CLI/API selector, including a readable blank-chain alias."""
    raw = str(chain_id)
    stripped = raw.strip()
    if raw == " " or stripped.lower() in {"blank", "<blank>"} or stripped == "_":
        return " "
    if len(stripped) != 1:
        raise ValueError(
            "Each selected strand chain ID must be exactly one character; "
            "use 'blank' or '_' for a blank chain ID."
        )
    return stripped


def _display_chain_id(chain_id: str) -> str:
    return "<blank>" if chain_id == " " else chain_id


def _parse_link_line(line: str) -> edit_pdb_link.pdb_link_record:
    """Parse a LINK, including the valid PDB case of a blank chain ID."""
    body, ending = split_line_ending(line)
    padded = body.ljust(80)
    blank1 = padded[21] == " "
    blank2 = padded[51] == " "
    fixed_resseq_fields = False
    if _record_name(line) == "LINK":
        try:
            int(padded[22:26].strip())
            int(padded[52:56].strip())
            fixed_resseq_fields = True
        except ValueError:
            pass

    if fixed_resseq_fields and (blank1 or blank2):
        # The shared token-based LINK parser cannot distinguish an omitted
        # blank chain field and can silently shift endpoint-2 fields. Parse
        # through temporary one-character sentinels, then restore the blanks.
        proxy = list(padded)
        if blank1:
            proxy[21] = "_"
        if blank2:
            proxy[51] = "_"
        parsed = edit_pdb_link.pdb_link_record("".join(proxy) + (ending or "\n"))
        if blank1:
            parsed.update_chainID1(" ")
        if blank2:
            parsed.update_chainID2(" ")
        return parsed

    return edit_pdb_link.pdb_link_record(line)


def _parse_input(
    lines: Sequence[str],
    selected_chain_ids: Sequence[str],
) -> _ScopedInput:
    if any(_record_name(line) in {"MODEL", "ENDMDL"} for line in lines):
        raise ValueError("Multi-model PDB files are not supported by this tool.")

    selected = set(selected_chain_ids)
    last_coordinate_chain: Optional[str] = None
    completed_chains: Set[str] = set()
    terminated_chains: Set[str] = set()
    last_residue_by_chain: Dict[str, Tuple[int, str]] = {}
    for line in lines:
        identity = coordinate_identity(line)
        if identity is not None and identity[0] in selected and identity[2].strip():
            raise ValueError(
                "Insertion codes are not supported; remove them before reversing strand direction."
            )
        if (
            _record_name(line) == "TER"
            and len(line) > 21
            and line[21] in selected
        ):
            terminated_chains.add(line[21])
        if _record_name(line) not in {"ATOM", "HETATM"} or identity is None:
            continue
        chain_id, resseq, icode, _resname = identity
        if chain_id not in selected:
            continue
        if chain_id in terminated_chains:
            raise ValueError(
                f"Chain '{chain_id}' has coordinate records after its TER record."
            )
        if last_coordinate_chain is not None and chain_id != last_coordinate_chain:
            completed_chains.add(last_coordinate_chain)
            if chain_id in completed_chains:
                raise ValueError(f"Chain '{chain_id}' is interleaved in the coordinate section.")
        residue_key = (resseq, icode)
        previous = last_residue_by_chain.get(chain_id)
        if previous is not None and residue_key != previous and residue_key < previous:
            raise ValueError(
                f"Chain '{chain_id}' residue blocks are not serialized in ascending number order."
            )
        last_residue_by_chain[chain_id] = residue_key
        last_coordinate_chain = chain_id

    raw_link_lines = [line for line in lines if _record_name(line) == "LINK"]
    links = [_parse_link_line(line) for line in raw_link_lines]
    records: List[object] = []
    ignored_links: List[edit_pdb_link.pdb_link_record] = []
    coordinate_text = "".join(line for line in lines if _record_name(line) != "LINK")
    edit_pdb_link.file2rec_link(StringIO(coordinate_text), records, ignored_links)
    atom_records = [
        record
        for record in records
        if getattr(record, "recordName", "") in {"ATOM", "HETATM"}
        and (getattr(record, "chainID", "").strip() or " ") in selected
    ]
    if not atom_records:
        raise ValueError("No ATOM/HETATM records were found for the selected strand(s).")
    nodes, label_to_idx = rex.build_residue_nodes(atom_records)

    scoped_links: List[edit_pdb_link.pdb_link_record] = []
    retained_links: List[_InputLink] = []
    for link, raw_line in zip(links, raw_link_lines):
        label1 = rex._input_link_label(link.chainID1, link.resSeq1)
        label2 = rex._input_link_label(link.chainID2, link.resSeq2)
        endpoint1_selected = label1[0] in selected
        endpoint2_selected = label2[0] in selected
        both_selected = endpoint1_selected and endpoint2_selected
        classified = rex._classify_backbone_link(link)

        if classified is not None and (endpoint1_selected != endpoint2_selected):
            selected_label = label1 if endpoint1_selected else label2
            other_label = label2 if endpoint1_selected else label1
            raise ValueError(
                f"Chain '{selected_label[0]}' has a backbone connection to another "
                f"chain ({other_label[0]}{other_label[1]}); reverse the complete "
                "connected strand after assigning it one chain ID."
            )

        if both_selected:
            # Let the selected topology builder validate and consume supported
            # backbone links. Non-backbone and same-residue links remain raw.
            scoped_links.append(link)
            if classified is not None and label1 != label2:
                continue
        retained_links.append(_InputLink(record=link, raw_line=raw_line))

    topology = rex.build_input_backbone_topology(nodes, label_to_idx, scoped_links)
    return _ScopedInput(
        nodes=nodes,
        label_to_idx=label_to_idx,
        topology=topology,
        input_links=links,
        retained_links=retained_links,
    )


def _chain_node_order(
    chain_id: str,
    nodes: Sequence[rex.ResidueNode],
) -> List[int]:
    order = [index for index, node in enumerate(nodes) if node.orig_chain_id == chain_id]
    order.sort(key=lambda index: nodes[index].input_chain_rank)
    return order


def _reject_ambiguous_two_residue_topology(
    chain_id: str,
    node_order: Sequence[int],
    nodes: Sequence[rex.ResidueNode],
    input_links: Sequence[edit_pdb_link.pdb_link_record],
) -> None:
    """Reject the PDB representation that cannot distinguish a path from a cycle.

    A two-residue chain has only one unordered residue pair.  An identity-
    symmetry P(first)--O3'(second) LINK can therefore describe either an
    inverted one-bond path or the second bond of a covalently closed two-bond
    cycle.  ``BackboneGraph`` is intentionally a simple graph, so guessing
    either interpretation would risk changing topology during reversal.
    """
    if len(node_order) != 2:
        return

    first_label = nodes[node_order[0]].orig_label()
    second_label = nodes[node_order[1]].orig_label()
    selected_pair = frozenset((first_label, second_label))
    for link in input_links:
        classified = rex._classify_backbone_link(link)
        if classified is None:
            continue
        kind, endpoint1, endpoint2 = classified
        if kind != "std":
            continue
        label1 = (link.chainID1.strip() or " ", int(link.resSeq1))
        label2 = (link.chainID2.strip() or " ", int(link.resSeq2))
        if frozenset((label1, label2)) != selected_pair:
            continue
        endpoints = {label1: endpoint1, label2: endpoint2}
        if endpoints.get(first_label) == "P" and endpoints.get(second_label) == "O3":
            raise ValueError(
                f"Chain '{chain_id}' has ambiguous two-residue topology: its "
                "P(first)--O3'(second) LINK can represent either an inverted open "
                "path or a covalently closed two-residue cycle. Reversal was not "
                "performed."
            )


def _validate_chain_component(
    chain_id: str,
    node_order: Sequence[int],
    graph: rex.BackboneGraph,
    nodes: Sequence[rex.ResidueNode],
) -> bool:
    if not node_order:
        raise ValueError(f"No coordinate records found for chain '{chain_id}'.")
    selected = set(node_order)
    for index in node_order:
        external = [neighbor for neighbor in graph.neigh[index] if neighbor not in selected]
        if external:
            labels = ", ".join(
                f"{nodes[neighbor].orig_chain_id}{nodes[neighbor].orig_res_seq}"
                for neighbor in external
            )
            raise ValueError(
                f"Chain '{chain_id}' has backbone connections to another chain ({labels}); "
                "reverse the complete connected strand after assigning it one chain ID."
            )

    visited: Set[int] = set()
    stack = [node_order[0]]
    while stack:
        current = stack.pop()
        if current in visited:
            continue
        visited.add(current)
        stack.extend(neighbor for neighbor in graph.neigh[current] if neighbor in selected)
    if visited != selected:
        raise ValueError(
            f"Chain '{chain_id}' contains multiple disconnected backbone components."
        )

    degrees = [len(graph.neigh[index]) for index in node_order]
    is_cycle = len(node_order) > 1 and all(degree == 2 for degree in degrees)
    if not is_cycle:
        ends = sum(degree <= 1 for degree in degrees)
        if len(node_order) > 1 and ends != 2:
            raise ValueError(
                f"Chain '{chain_id}' is not a single unbranched path or cycle."
            )
    return is_cycle


def _structured_remark_update(
    line: str,
    chain_id: str,
    mapping: Dict[int, int],
) -> Tuple[str, int]:
    if not line.startswith("REMARK 950 RE_SCRIPT "):
        return line, 0
    if line.startswith("REMARK 950 RE_SCRIPT COMMAND "):
        return line, 0
    if line.startswith("REMARK 950 RE_SCRIPT SOFTWARE "):
        return line, 0
    return replace_remark_residue_references(line, chain_id, mapping)


def _reverse_chain_lines(
    lines: Sequence[str],
    chain_id: str,
    is_cycle: bool,
) -> Tuple[List[str], ReverseStats, Dict[int, int]]:
    blocks, first_coord_index, last_coord_index = collect_residue_blocks(lines, chain_id)
    if any(block.old_icode.strip() for block in blocks):
        raise ValueError(f"Chain '{chain_id}' uses insertion codes, which are unsupported.")
    for line in lines[first_coord_index:last_coord_index]:
        if _record_name(line) == "TER" and len(line) > 21 and line[21] == chain_id:
            raise ValueError(
                f"Chain '{chain_id}' occurs in multiple TER-separated coordinate sections."
            )

    if is_cycle and len(blocks) > 1:
        reordered = [blocks[0], *reversed(blocks[1:])]
    else:
        reordered = list(reversed(blocks))

    numbering_start = min(block.old_resseq for block in blocks)
    numbering_end = numbering_start + len(blocks) - 1
    format_resseq(numbering_start)
    format_resseq(numbering_end)
    mapping: Dict[int, int] = {}
    for offset, block in enumerate(reordered):
        block.new_resseq = numbering_start + offset
        mapping[block.old_resseq] = block.new_resseq

    stats = ReverseStats(
        chain_id=chain_id,
        residue_count=len(blocks),
        is_cycle=is_cycle,
        numbering_start=numbering_start,
        numbering_end=numbering_end,
        old_first_resseq=blocks[0].old_resseq,
        old_last_resseq=blocks[-1].old_resseq,
        new_first_source_resseq=reordered[0].old_resseq,
        new_last_source_resseq=reordered[-1].old_resseq,
    )

    inventory_indices = [
        index for index, line in enumerate(lines) if is_chain_inventory_remark(line, chain_id)
    ]
    inventory_first = inventory_indices[0] if inventory_indices else None
    inventory_set = set(inventory_indices)
    inventory_lines: List[str] = []
    if inventory_first is not None:
        _body, ending = split_line_ending(lines[inventory_first])
        inventory_lines = build_chain_inventory_remarks(
            chain_id,
            reordered,
            ending or "\n",
        )
        stats.remark_inventory_rebuilt = True

    transformed_blocks: List[str] = []
    for block in reordered:
        assert block.new_resseq is not None
        for line in block.lines:
            transformed = update_coordinate_residue(line, block.new_resseq)
            if transformed != line:
                stats.coordinate_lines_changed += 1
            transformed_blocks.append(transformed)

    output: List[str] = []
    emitted_coordinates = False
    emitted_inventory = False
    last_block = reordered[-1]
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

        updated = line
        record = _record_name(line)
        if record == "HET":
            updated, changed = update_het_line(updated, chain_id, mapping)
            if changed:
                stats.het_lines_changed += 1
        elif record == "TER":
            updated, changed = update_ter_line(
                updated,
                chain_id,
                last_block,
                numbering_end,
            )
            if changed:
                stats.ter_lines_changed += 1
        elif record == "REMARK":
            updated, changed_count = _structured_remark_update(
                updated,
                chain_id,
                mapping,
            )
            stats.remark_references_changed += changed_count
        output.append(updated)

    return output, stats, mapping


def _output_atom_records(lines: Sequence[str]):
    records: List[object] = []
    ignored_links: List[edit_pdb_link.pdb_link_record] = []
    edit_pdb_link.file2rec_link(StringIO("".join(lines)), records, ignored_links)
    return [
        record
        for record in records
        if getattr(record, "recordName", "") in {"ATOM", "HETATM"}
    ]


def _output_residue_positions(lines: Sequence[str]) -> Dict[Tuple[str, int], int]:
    positions: Dict[Tuple[str, int], int] = {}
    next_by_chain: Dict[str, int] = {}
    for line in lines:
        if _record_name(line) not in {"ATOM", "HETATM"}:
            continue
        identity = coordinate_identity(line)
        if identity is None:
            continue
        chain_id, resseq, _icode, _resname = identity
        key = (chain_id, resseq)
        if key not in positions:
            positions[key] = next_by_chain.get(chain_id, 0)
            next_by_chain[chain_id] = positions[key] + 1
    return positions


def _endpoint_atom(
    atom_index,
    node: rex.ResidueNode,
    endpoint: str,
):
    candidates = {
        "P": ["P"],
        "O3": ["O3'", "O3*", "O3"],
        "O5": ["O5'", "O5*", "O5"],
    }.get(endpoint)
    if candidates is None:
        raise ValueError(f"Unsupported backbone endpoint type: {endpoint}")
    atom = rex._find_required_atom(
        atom_index,
        node.new_chain_id,
        node.new_res_seq,
        candidates,
    )
    if atom is None:
        raise ValueError(
            f"Cannot preserve backbone topology for {node.new_chain_id}{node.new_res_seq}: "
            f"missing endpoint atom {endpoint}. This tool currently supports nucleic-acid strands."
        )
    return atom


def _edge_is_implicit_in_output(
    edge: rex.Edge,
    nodes: Sequence[rex.ResidueNode],
    positions: Dict[Tuple[str, int], int],
) -> bool:
    node_a = nodes[edge.a]
    node_b = nodes[edge.b]
    if node_a.new_chain_id != node_b.new_chain_id:
        return False
    if node_a.is_phos_bridge or node_b.is_phos_bridge:
        return False
    pos_a = positions[(node_a.new_chain_id, node_a.new_res_seq)]
    pos_b = positions[(node_b.new_chain_id, node_b.new_res_seq)]
    if pos_b == pos_a + 1:
        return edge.endpoints(edge.a, edge.b) == ("O3", "P")
    if pos_a == pos_b + 1:
        return edge.endpoints(edge.b, edge.a) == ("O3", "P")
    return False


def _format_edge_link(
    edge: rex.Edge,
    nodes: Sequence[rex.ResidueNode],
    atom_index,
) -> str:
    node_a = nodes[edge.a]
    node_b = nodes[edge.b]
    end_a, end_b = edge.end_a, edge.end_b

    if {end_a, end_b} == {"P", "O3"}:
        if end_a == "P":
            first_node, first_end, second_node, second_end = node_a, end_a, node_b, end_b
        else:
            first_node, first_end, second_node, second_end = node_b, end_b, node_a, end_a
    elif {end_a, end_b} == {"P", "O5"}:
        if end_a == "O5":
            first_node, first_end, second_node, second_end = node_a, end_a, node_b, end_b
        else:
            first_node, first_end, second_node, second_end = node_b, end_b, node_a, end_a
    else:
        first_node, first_end, second_node, second_end = node_a, end_a, node_b, end_b

    atom1 = _endpoint_atom(atom_index, first_node, first_end)
    atom2 = _endpoint_atom(atom_index, second_node, second_end)
    return rex._format_link_line(atom1, atom2)


def _regenerate_links(
    lines: Sequence[str],
    nodes: List[rex.ResidueNode],
    topology: rex.InputBackboneTopology,
    retained_links: Sequence[_InputLink],
    mappings: Dict[str, Dict[int, int]],
) -> Tuple[List[str], int, int]:
    for node in nodes:
        node.new_chain_id = node.orig_chain_id
        node.new_res_seq = mappings.get(node.orig_chain_id, {}).get(
            node.orig_res_seq,
            node.orig_res_seq,
        )

    output_atoms = _output_atom_records(lines)
    atom_index = edit_pdb_link.build_atom_index(output_atoms)
    positions = _output_residue_positions(lines)
    topology_links: List[str] = []
    for edge in topology.graph.edges.values():
        if _edge_is_implicit_in_output(edge, nodes, positions):
            continue
        topology_links.append(_format_edge_link(edge, nodes, atom_index))

    # LINKs outside the selected backbone graph stay byte-for-byte unchanged
    # unless one of their selected endpoint labels actually moves. Remapping
    # only the changed residue-number field avoids validating or rewriting an
    # unrelated, possibly non-nucleic-acid endpoint on an unselected chain.
    preserved_lines: List[str] = []
    for item in retained_links:
        source = item.record
        new_resseq1 = mappings.get(source.chainID1.strip() or " ", {}).get(
            source.resSeq1,
            source.resSeq1,
        )
        new_resseq2 = mappings.get(source.chainID2.strip() or " ", {}).get(
            source.resSeq2,
            source.resSeq2,
        )
        if new_resseq1 == source.resSeq1 and new_resseq2 == source.resSeq2:
            preserved_lines.append(item.raw_line)
            continue

        remapped = _parse_link_line(source.string)
        if new_resseq1 != source.resSeq1:
            remapped.update_resSeq1(new_resseq1)
        if new_resseq2 != source.resSeq2:
            remapped.update_resSeq2(new_resseq2)
        preserved_lines.append(remapped.string)

    link_lines = topology_links + preserved_lines
    return link_lines, len(topology_links), len(preserved_lines)


def _tool_header_lines(
    chains: Sequence[str],
    stats: Sequence[ReverseStats],
    topology_links: int,
    preserved_links: int,
    total_links: int,
    command_text: Optional[str],
    line_ending: str,
) -> List[str]:
    ending = line_ending or "\n"
    chain_tokens = ["blank" if chain == " " else chain for chain in chains]
    lines = [
        f"REMARK 950 REVERSE_STRAND SOFTWARE name=reverse_strand_direction version=V{VERSION}{ending}",
    ]
    if command_text:
        cleaned = " ".join(str(command_text).splitlines())
        lines.append(f"REMARK 950 REVERSE_STRAND COMMAND text={cleaned}{ending}")
    lines.append(
        "REMARK 950 REVERSE_STRAND EVENT action=reverse_serialization "
        f"chains={','.join(chain_tokens)} topology=preserved{ending}"
    )
    for item in stats:
        chain_token = "blank" if item.chain_id == " " else item.chain_id
        lines.append(
            "REMARK 950 REVERSE_STRAND CHAIN "
            f"chain={chain_token} type={'cycle' if item.is_cycle else 'path'} "
            f"count={item.residue_count} old_first={chain_token}:{item.old_first_resseq} "
            f"old_last={chain_token}:{item.old_last_resseq} "
            f"new_first_source={chain_token}:{item.new_first_source_resseq} "
            f"new_last_source={chain_token}:{item.new_last_source_resseq}{ending}"
        )
    lines.append(
        "REMARK 950 REVERSE_STRAND LINK_RECORDS "
        f"total={total_links} regenerated_backbone={topology_links} "
        f"preserved_other={preserved_links}{ending}"
    )
    return lines


def reverse_strand_lines(
    lines: Sequence[str],
    chain_ids: Sequence[str],
    command_text: Optional[str] = None,
) -> Tuple[List[str], List[ReverseStats], int, int, int]:
    original_lines = list(lines)

    present_order: List[str] = []
    for line in original_lines:
        if _record_name(line) not in {"ATOM", "HETATM"}:
            continue
        identity = coordinate_identity(line)
        if identity is not None and identity[0] not in present_order:
            present_order.append(identity[0])
    requested: Set[str] = set()
    for chain_id in chain_ids:
        requested.add(_normalize_chain_selector(chain_id))
    if not requested:
        raise ValueError("Select at least one strand chain ID to reverse.")
    missing = sorted(requested.difference(present_order))
    if missing:
        raise ValueError(
            "Selected chain(s) not found: "
            + ", ".join(_display_chain_id(chain) for chain in missing)
        )
    normalized_chains = [chain for chain in present_order if chain in requested]
    scoped = _parse_input(original_lines, normalized_chains)
    nodes = scoped.nodes
    topology = scoped.topology

    chain_cycles: Dict[str, bool] = {}
    for chain_id in normalized_chains:
        node_order = _chain_node_order(chain_id, nodes)
        _reject_ambiguous_two_residue_topology(
            chain_id,
            node_order,
            nodes,
            scoped.input_links,
        )
        chain_cycles[chain_id] = _validate_chain_component(
            chain_id,
            node_order,
            topology.graph,
            nodes,
        )

    output_lines = original_lines
    stats: List[ReverseStats] = []
    mappings: Dict[str, Dict[int, int]] = {}
    for chain_id in normalized_chains:
        output_lines, item, mapping = _reverse_chain_lines(
            output_lines,
            chain_id,
            chain_cycles[chain_id],
        )
        stats.append(item)
        mappings[chain_id] = mapping

    # The original RE_SCRIPT link-count inventory is no longer authoritative.
    body_without_links: List[str] = []
    for line in output_lines:
        if _record_name(line) == "LINK":
            continue
        if line.startswith("REMARK 950 RE_SCRIPT SPECIAL event=link_records "):
            continue
        body_without_links.append(line)

    new_links, topology_count, preserved_count = _regenerate_links(
        body_without_links,
        nodes,
        topology,
        scoped.retained_links,
        mappings,
    )
    ending = next(
        (split_line_ending(line)[1] for line in original_lines if split_line_ending(line)[1]),
        "\n",
    )
    header = _tool_header_lines(
        normalized_chains,
        stats,
        topology_count,
        preserved_count,
        len(new_links),
        command_text,
        ending,
    )
    insert_at = next(
        (
            index
            for index, line in enumerate(body_without_links)
            if _record_name(line) in {"ATOM", "HETATM", "ANISOU", "SIGATM", "SIGUIJ", "TER"}
        ),
        len(body_without_links),
    )
    link_lines = [line.rstrip("\r\n") + ending for line in new_links]
    final_lines = [
        *header,
        *body_without_links[:insert_at],
        *link_lines,
        *body_without_links[insert_at:],
    ]
    return final_lines, stats, topology_count, preserved_count, len(new_links)


def reverse_strands(
    input_pdb: Path,
    chain_ids: Sequence[str],
    output_pdb: Optional[Path] = None,
    command_text: Optional[str] = None,
    verbose: bool = True,
) -> ReverseResult:
    input_path = Path(input_pdb)
    if not input_path.is_file():
        raise ValueError(f"Input PDB does not exist: {input_path}")
    output_path = Path(output_pdb) if output_pdb is not None else default_output_path(input_path)
    try:
        same_path = output_path.resolve() == input_path.resolve()
    except OSError:
        same_path = output_path.absolute() == input_path.absolute()
    if same_path:
        raise ValueError("Input and output paths must differ; the input PDB will not be overwritten.")

    lines = input_path.read_text(errors="replace").splitlines(True)
    transformed, stats, topology_count, preserved_count, total_count = reverse_strand_lines(
        lines,
        chain_ids,
        command_text=command_text,
    )
    output_path.write_text("".join(transformed))
    result = ReverseResult(
        output_path=output_path,
        stats=stats,
        topology_link_count=topology_count,
        preserved_link_count=preserved_count,
        total_link_count=total_count,
    )
    if verbose:
        print(format_summary(result), flush=True)
    return result


def format_summary(result: ReverseResult) -> str:
    lines = [
        "Reversed strands: "
        + ", ".join(_display_chain_id(item.chain_id) for item in result.stats),
    ]
    for item in result.stats:
        displayed_chain = _display_chain_id(item.chain_id)
        lines.append(
            f"  {displayed_chain}: {item.residue_count} residues, "
            f"{'cycle anchored at its first residue' if item.is_cycle else 'path fully reversed'}, "
            f"source {displayed_chain}{item.new_first_source_resseq} is now first"
        )
    lines.extend(
        [
            f"LINK records: {result.total_link_count} "
            f"({result.topology_link_count} regenerated backbone, "
            f"{result.preserved_link_count} preserved other)",
            "Coordinates and atom identities were preserved; only residue serialization, numbering, and topology records changed.",
            f"Wrote: {result.output_path}",
        ]
    )
    return "\n".join(lines)


def build_cli_command(
    script_name: str,
    input_pdb: str,
    chain_ids: Sequence[str],
    output_pdb: str,
) -> str:
    parts = [sys.executable, script_name, input_pdb]
    for chain_id in chain_ids:
        parts.extend(["--strand", "blank" if chain_id == " " else chain_id])
    if output_pdb:
        parts.extend(["-o", output_pdb])
    return " ".join(shlex.quote(str(part)) for part in parts)


def _detected_chains(path: Path) -> List[Tuple[str, int, int, int]]:
    lines = path.read_text(errors="replace").splitlines(True)
    seen: Dict[str, List[Tuple[int, str]]] = {}
    for line in lines:
        if _record_name(line) not in {"ATOM", "HETATM"}:
            continue
        identity = coordinate_identity(line)
        if identity is None:
            continue
        chain_id, resseq, icode, _resname = identity
        key = (resseq, icode)
        if key not in seen.setdefault(chain_id, []):
            seen[chain_id].append(key)
    return [
        (chain_id, len(residues), residues[0][0], residues[-1][0])
        for chain_id, residues in seen.items()
    ]


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
    root.geometry("820x620")

    input_var = tk.StringVar(value=initial_path or "")
    output_var = tk.StringVar()
    status_var = tk.StringVar(value="Choose an input PDB and at least one strand.")
    chain_vars: Dict[str, tk.BooleanVar] = {}

    outer = ttk.Frame(root, padding=10)
    outer.pack(fill="both", expand=True)
    outer.columnconfigure(1, weight=1)
    outer.rowconfigure(5, weight=1)
    ttk.Label(outer, text=f"{TOOL_NAME} V{VERSION}", font=("TkDefaultFont", 12, "bold")).grid(
        row=0, column=0, columnspan=3, sticky="w", pady=(0, 8)
    )

    ttk.Label(outer, text="Input PDB").grid(row=1, column=0, sticky="e", padx=(0, 6), pady=4)
    ttk.Entry(outer, textvariable=input_var, width=72).grid(row=1, column=1, sticky="ew", pady=4)
    ttk.Label(outer, text="Output PDB").grid(row=2, column=0, sticky="e", padx=(0, 6), pady=4)
    ttk.Entry(outer, textvariable=output_var, width=72).grid(row=2, column=1, sticky="ew", pady=4)

    chain_box = ttk.LabelFrame(outer, text="Strands to reverse", padding=8)
    chain_box.grid(row=3, column=0, columnspan=3, sticky="ew", pady=8)
    explanation = (
        "Complete residue blocks are reversed and renumbered; coordinates are unchanged. "
        "This changes PDB serialization, not molecular geometry or base identity."
    )
    ttk.Label(outer, text=explanation, wraplength=780, justify="left").grid(
        row=4, column=0, columnspan=3, sticky="w", pady=(0, 6)
    )

    log = scrolledtext.ScrolledText(outer, wrap="word", height=17)
    log.grid(row=5, column=0, columnspan=3, sticky="nsew", pady=6)
    ttk.Label(outer, textvariable=status_var).grid(row=6, column=0, columnspan=3, sticky="w")
    buttons = ttk.Frame(outer)
    buttons.grid(row=7, column=0, columnspan=3, sticky="w", pady=(8, 0))

    def append(message: str) -> None:
        log.insert("end", message)
        log.see("end")
        log.update_idletasks()

    def refresh_chains() -> None:
        for child in chain_box.winfo_children():
            child.destroy()
        chain_vars.clear()
        path_text = input_var.get().strip()
        if not path_text:
            return
        path = Path(path_text)
        try:
            detected = _detected_chains(path)
        except Exception as exc:
            status_var.set(str(exc))
            return
        if not detected:
            status_var.set("No ATOM/HETATM chains detected.")
            return
        for index, (chain_id, count, first, last) in enumerate(detected):
            variable = tk.BooleanVar(value=False)
            chain_vars[chain_id] = variable
            displayed_chain = _display_chain_id(chain_id)
            label = (
                f"{displayed_chain}  ({count} residues; "
                f"{displayed_chain}{first} … {displayed_chain}{last})"
            )
            ttk.Checkbutton(chain_box, text=label, variable=variable).grid(
                row=index // 3,
                column=index % 3,
                sticky="w",
                padx=(0, 16),
                pady=2,
            )
        output_var.set(str(default_output_path(path)))
        status_var.set(f"Detected {len(detected)} strand(s).")

    def browse_input() -> None:
        selected = filedialog.askopenfilename(
            title="Choose input PDB",
            filetypes=[("PDB files", "*.pdb"), ("All files", "*")],
        )
        if selected:
            input_var.set(selected)
            refresh_chains()

    def browse_output() -> None:
        selected = filedialog.asksaveasfilename(
            title="Choose output PDB",
            defaultextension=".pdb",
            filetypes=[("PDB files", "*.pdb"), ("All files", "*")],
        )
        if selected:
            output_var.set(selected)

    ttk.Button(outer, text="Browse", command=browse_input).grid(row=1, column=2, padx=(6, 0))
    ttk.Button(outer, text="Browse", command=browse_output).grid(row=2, column=2, padx=(6, 0))

    def run() -> None:
        selected_chains = [chain for chain, variable in chain_vars.items() if variable.get()]
        if not selected_chains:
            messagebox.showerror(TOOL_NAME, "Select at least one strand to reverse.", parent=root)
            return
        command = build_cli_command(
            str(Path(__file__).resolve()),
            input_var.get().strip(),
            selected_chains,
            output_var.get().strip(),
        )
        append("Equivalent CLI command:\n    " + command + "\n\n")
        try:
            result = reverse_strands(
                Path(input_var.get().strip()),
                selected_chains,
                Path(output_var.get().strip()),
                command_text=command,
                verbose=False,
            )
        except Exception as exc:
            status_var.set(f"Error: {exc}")
            append(f"Error: {exc}\n")
            messagebox.showerror(TOOL_NAME, str(exc), parent=root)
            return
        summary = format_summary(result)
        append(summary + "\n")
        print(summary, flush=True)
        status_var.set(f"Wrote {result.output_path}")

    ttk.Button(buttons, text="Reverse selected strands", command=run).pack(side="left")
    ttk.Button(buttons, text="Refresh strands", command=refresh_chains).pack(side="left", padx=6)
    ttk.Button(buttons, text="Close", command=root.destroy).pack(side="left")

    if input_var.get().strip():
        refresh_chains()
    root.mainloop()
    return 0


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="reverse_strand_direction.py",
        description=(
            "Reverse complete residue blocks for selected nucleic-acid strands, "
            "renumber continuously, and regenerate LINK topology without moving atoms."
        ),
    )
    parser.add_argument("pdb", nargs="?", help="Input PDB file.")
    parser.add_argument(
        "--strand",
        "--reverse",
        action="append",
        dest="strands",
        metavar="CHAIN",
        help=(
            "One-character chain ID to reverse; repeat for multiple strands. "
            "Use 'blank' or '_' for a blank chain ID."
        ),
    )
    parser.add_argument("-o", "--output", help="Output PDB. Default: <input>_strand_reversed.pdb")
    parser.add_argument("--gui", action="store_true", help="Open the Tk GUI.")
    parser.add_argument("-v", "--version", action="version", version=f"{TOOL_NAME} V{VERSION}")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    if args.gui:
        return run_gui(args.pdb)
    if not args.pdb:
        parser.error("an input PDB is required unless --gui is used")
    if not args.strands:
        parser.error("at least one --strand CHAIN selection is required")
    command = " ".join(shlex.quote(arg) for arg in [sys.executable, *sys.argv])
    try:
        reverse_strands(
            Path(args.pdb),
            args.strands,
            Path(args.output) if args.output else None,
            command_text=command,
            verbose=True,
        )
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
