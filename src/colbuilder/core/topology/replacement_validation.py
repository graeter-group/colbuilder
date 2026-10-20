"""Verify model attachment from actual final ITP bonds and source residue maps."""

from collections import Counter
from pathlib import Path

from .crosslink_validation import canon, sections
from ..geometry.crosslink_network import CrosslinkIntegrityError
from ..geometry.replacement_policy import model_incidence, validate_attachment_plan


def validate_replacement_itps(network, groups):
    """Each group has an ITP, its source residue map and optional PDB residue numbers."""
    report = network.replacement_report
    if not report or report.get("mode") != "preserve_attachment":
        return None
    validate_attachment_plan(report, network.entities)
    observed, seen_models, seen_residues = Counter(), set(), set()
    files = []
    for group in groups:
        itp, residue_map, *number_lists = group
        itp = Path(itp)
        if not itp.is_file():
            raise CrosslinkIntegrityError(f"Missing final ITP for attachment validation: {itp}")
        if set(residue_map.values()) != set(range(len(residue_map))):
            raise CrosslinkIntegrityError(f"Invalid source residue map for {itp}")
        sources = {ordinal: key for key, ordinal in residue_map.items()}
        numbers = number_lists[0] if number_lists else [sources[i][1] for i in range(len(sources))]
        if len(numbers) != len(sources):
            raise CrosslinkIntegrityError(f"Invalid source residue numbers for {itp}")
        overlap = seen_models & {key[0] for key in residue_map}
        if overlap:
            raise CrosslinkIntegrityError(f"Models duplicated across final ITPs: {sorted(overlap)}")
        seen_models.update(key[0] for key in residue_map)
        if seen_residues & residue_map.keys():
            raise CrosslinkIntegrityError(f"Residues duplicated across final ITPs: {itp}")
        seen_residues.update(residue_map)
        atoms, residues, bonds = {}, [], []
        for section, fields, _ in sections(itp):
            if section == "atoms":
                aid, key, atom = int(fields[0]), (fields[2], fields[3]), fields[4]
                if aid in atoms:
                    raise CrosslinkIntegrityError(f"Duplicate atom index in {itp}: {aid}")
                if not residues or residues[-1][0] != key or atom in residues[-1][1]:
                    residues.append((key, set()))
                residues[-1][1].add(atom)
                ordinal = len(residues) - 1
                source = sources.get(ordinal)
                if source is None or (numbers[ordinal], source[3]) != key:
                    raise CrosslinkIntegrityError(f"PDB/ITP residue identity mismatch in {itp}: {key}")
                atoms[aid] = source + (atom,)
            elif section == "bonds":
                bonds.append(fields)
        if len(residues) != len(residue_map) or not atoms:
            raise CrosslinkIntegrityError(f"PDB/ITP residue inventory mismatch in {itp}")
        for fields in bonds:
            a, b = map(int, fields[:2])
            if a not in atoms or b not in atoms:
                raise CrosslinkIntegrityError(f"Bond outside atom inventory in {itp}: {a}, {b}")
            if atoms[a][0] != atoms[b][0]:
                if int(fields[2]) != 1:
                    raise CrosslinkIntegrityError(f"Unexpected inter-model bond function in {itp}: {fields}")
                observed[canon((atoms[a], atoms[b]))] += 1
        files.append(itp.name)
    expected = Counter(canon(bond) for entity in network.entities for bond in entity.bonds)
    if observed != expected:
        raise CrosslinkIntegrityError(
            "Final ITP crosslinks differ from replacement plan: "
            f"missing={list((expected-observed).elements())[:6]}, "
            f"unexpected/duplicated={list((observed-expected).elements())[:6]}"
        )
    # A complete PYD counts only when BOTH forming bonds are present once.
    confirmed = [entity for entity in network.entities
                 if all(observed[canon(bond)] == 1 for bond in entity.bonds)]
    degree = model_incidence(confirmed)
    protected = report["attachment"]["protected_models"]
    missing = sorted(mid for mid in protected if mid not in seen_models or not degree[mid])
    if missing:
        raise CrosslinkIntegrityError(f"Final ITPs leave protected models without crosslinks: {missing}")
    return {
        "status": "passed",
        "mode": "preserve_attachment",
        "itp_files": files,
        "protected_models": protected,
        "observed_complete_crosslinks": len(confirmed),
        "observed_intermodel_bonds": sum(observed.values()),
        "degree_from_itp_bonds": {str(mid): degree[mid] for mid in protected},
        "newly_unlinked_models": [],
        "mechanical_stability_tested": False,
    }
