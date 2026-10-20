"""Report glucosepane overlaps in the complete exported AMBER structure."""

from collections import defaultdict
from pathlib import Path
import json

import numpy as np
from scipy.spatial import cKDTree

from colbuilder.core.topology.crosslink_validation import sections
from colbuilder.core.utils.exceptions import TopologyGenerationError
from colbuilder.core.utils.logger import setup_logger
from colbuilder.core.utils.ring_geometry import RingScreen, is_hydrogen, residue_ring_indices
from colbuilder.core.utils.steric_policy import (
    CATASTROPHIC_HEAVY_A,
    CATASTROPHIC_HYDROGEN_A,
    WARN_14_A,
    WARN_HEAVY_A,
    WARN_HYDROGEN_A,
)

LOG = setup_logger(__name__)


def _atom_label(index, atom):
    residue, name, filename, local_id, residue_id, atom_type = atom
    return (
        f"global atom {index + 1} {residue}{residue_id}:{name} "
        f"(type {atom_type}, {filename}, local atom {local_id})"
    )


def validate_glucosepane_contacts(gro_path, itp_paths):
    """Check actual GRO coordinates, including H and contacts between ITPs.

    Finite overlaps only warn; invalid coordinates or GRO/ITP mappings still fail.
    This nonperiodic distance screening is not an energy/convergence test.
    The intended solvated simulation box must be checked separately after setup.
    """
    atoms, adjacency = [], defaultdict(set)
    ring_residues, ring_indices = [], {}
    for path in itp_paths:
        offset, local = len(atoms), []
        previous, residue_ordinal, names = None, -1, set()
        for section, fields, _ in sections(path):
            if section == "atoms":
                if int(fields[0]) != len(local) + 1:
                    raise TopologyGenerationError(
                        f"Noncontiguous atom numbering in {path}"
                    )
                local.append((fields[3], fields[4], Path(path).name, int(fields[0]),
                              fields[2], fields[1]))
                key = (fields[2], fields[3])
                if key != previous or fields[4] in names:
                    residue_ordinal += 1
                    names, previous = set(), key
                    ring_residues.append(((offset, residue_ordinal), fields[3]))
                names.add(fields[4])
                ring_indices[(offset, residue_ordinal, fields[4])] = offset + len(local) - 1
            elif section == "bonds":
                a, b = [int(x) - 1 + offset for x in fields[:2]]
                adjacency[a].add(b)
                adjacency[b].add(a)
        atoms.extend(local)
    selected = [i for i, atom in enumerate(atoms) if atom[0] in {"AGS", "LGX"}]
    if not selected:
        return {"glucosepane_atoms": 0, "gross_overlaps": 0, "compressed_contacts": 0}
    with Path(gro_path).open() as handle:
        handle.readline()
        count = int(handle.readline())
        if count != len(atoms):
            raise TopologyGenerationError(
                "GRO/ITP atom count mismatch during steric validation"
            )
        coordinates = []
        for expected in atoms:
            line = handle.readline()
            if (line[5:10].strip(), line[10:15].strip()) != expected[:2]:
                raise TopologyGenerationError(
                    "GRO/ITP atom order mismatch during steric validation"
                )
            coordinates.append([float(line[j : j + 8]) * 10 for j in (20, 28, 36)])
    coordinates = np.array(coordinates)
    if not np.isfinite(coordinates).all():
        raise TopologyGenerationError("Nonfinite coordinates in exported GRO")
    tree = cKDTree(coordinates)
    clashes, compressed, seen = [], [], set()
    for i in selected:
        excluded = {i} | adjacency[i] | {k for j in adjacency[i] for k in adjacency[j]}
        pairs14 = {k for j in excluded for k in adjacency[j]} - excluded
        for j in tree.query_ball_point(coordinates[i], WARN_HEAVY_A):
            pair = tuple(sorted((i, j)))
            if j in excluded or pair in seen:
                continue
            seen.add(pair)
            hydrogen = any(
                atoms[k][1].lstrip("0123456789").startswith("H") for k in pair
            )
            d = float(np.linalg.norm(coordinates[i] - coordinates[j]))
            severe_limit = CATASTROPHIC_HYDROGEN_A if hydrogen else CATASTROPHIC_HEAVY_A
            warning_limit = WARN_HYDROGEN_A if hydrogen else (
                WARN_14_A if j in pairs14 else WARN_HEAVY_A
            )
            if d < severe_limit:
                clashes.append((d, i, j))
            elif d < warning_limit:
                compressed.append((d, i, j))
    if clashes:
        LOG.warning(
            "Glucosepane geometry warning: %d gross nonbonded overlaps. "
            "Export continues; finite input coordinates do not guarantee finite forces "
            "or MD stability. Minimize and verify convergence before MD; "
            "do not add exclusions.",
            len(clashes),
        )
        for d, i, j in sorted(clashes):
            LOG.warning(
                "Glucosepane overlap warning: %.4f A between %s and %s.",
                d, _atom_label(i, atoms[i]), _atom_label(j, atoms[j]),
            )
    if compressed:
        d, i, j = min(compressed)
        LOG.warning(
            "Glucosepane geometry warning: %d compressed nonbonded contacts; "
            "closest %.4f A between %s and %s. Export continues. "
            "Minimize and check convergence before MD; do not add exclusions.",
            len(compressed), d, _atom_label(i, atoms[i]), _atom_label(j, atoms[j]),
        )
    rings, incomplete = residue_ring_indices(ring_residues, ring_indices)
    bonds = [(a, b) for a, neighbors in adjacency.items() for b in neighbors
             if a < b and not is_hydrogen(atoms[a][1]) and not is_hydrogen(atoms[b][1])]
    hits = RingScreen(coordinates, bonds, rings,
                      labels=[_atom_label(i, a) for i, a in enumerate(atoms)]).check(coordinates)
    ring_report = {"schema": 1, "scope": "actual exported GRO and ITP bonds; nonperiodic",
                   "rings_checked": len(rings), "incomplete_rings": incomplete, "contacts": hits,
                   "penetrations": sum(hit["status"] == "penetration" for hit in hits)}
    Path(gro_path).with_suffix(".rings.json").write_text(json.dumps(ring_report, indent=2) + "\n")
    for hit in hits:
        LOG.warning("Exported GRO ring contact (%s): ring=%s, crossing bond=%s. Export continues; inspect geometry.",
                    hit["status"], hit["ring_atoms"], hit["bond_atoms"])
    return {"glucosepane_atoms": len(selected), "gross_overlaps": len(clashes),
            "compressed_contacts": len(compressed), "ring_geometry": ring_report}
