"""Screen composed caps; repair penetrations with native sidechain torsions only."""
import hashlib
import json
from collections import defaultdict
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.optimize import differential_evolution
from scipy.spatial import cKDTree

from .backbone import BACKBONE, _xyz, read_pdb_residues
from .crosslink_network import CrosslinkIntegrityError
from colbuilder.core.sequence.steric_optimization import RADII, SidechainTorsions
from colbuilder.core.utils.logger import setup_logger
from colbuilder.core.utils.ring_geometry import (
    RingScreen, is_hydrogen, residue_ring_indices, residue_templates,
)

LOG = setup_logger(__name__)
FLEXIBLE_NATIVE = {"SER", "THR", "LYS", "ARG", "GLU", "GLN", "ASP", "ASN", "MET", "LEU", "ILE", "VAL"}


class CapsGeometry:
    def __init__(self, caps, network=None):
        self.documents, self.residues, self.indices, self.labels = {}, {}, {}, []
        self.atom_residues, self.adjacency, xyz = [], defaultdict(set), []
        templates = residue_templates()
        self.unsupported = set()
        for mid, path in sorted(caps.items()):
            lines = Path(path).read_text().splitlines(True)
            residues = read_pdb_residues(lines, str(path))
            self.documents[mid] = (Path(path), lines)
            for residue in residues:
                key = (mid, *residue.address)
                self.residues[key] = residue
                for name, line in residue.atoms.items():
                    if is_hydrogen(name):
                        continue
                    self.indices[(*key, name)] = len(xyz)
                    xyz.append(_xyz(line))
                    self.atom_residues.append(key)
                    self.labels.append(dict(model=int(mid), chain=residue.address[1],
                                            residue=residue.name + residue.address[2] + residue.address[3], atom=name))
            for ri, residue in enumerate(residues):
                template = templates.get("HIE" if residue.name == "HIS" else residue.name)
                if template is None:
                    self.unsupported.add(residue.name)
                    continue
                for pair in template["bonds"]:
                    ids = []
                    for name in pair:
                        other = ri
                        if name.startswith(("-", "+")):
                            other += -1 if name[0] == "-" else 1
                            name = name[1:]
                            if not 0 <= other < len(residues) or residues[other].segment != residue.segment:
                                break
                        index = self.indices.get((mid, *residues[other].address, name))
                        if index is None:
                            break
                        ids.append(index)
                    if len(ids) == 2:
                        self.adjacency[ids[0]].add(ids[1])
                        self.adjacency[ids[1]].add(ids[0])
        self.coordinates = np.array(xyz, dtype=float)
        if network is not None:
            marker_indices = {(k[0], k[3]+k[4], k[2], self.residues[k[:-1]].name, k[-1]): i
                              for k, i in self.indices.items()}
            for entity in network.entities:
                for a, b in entity.bonds:
                    # Marker port names are carried by the chemical network.
                    ia, ib = marker_indices[a], marker_indices[b]
                    self.adjacency[ia].add(ib)
                    self.adjacency[ib].add(ia)
        self.bonds = [(a, b) for a, neighbors in self.adjacency.items() for b in neighbors if a < b]
        self.rings, self.incomplete = residue_ring_indices(
            [(key, r.name) for key, r in self.residues.items()], self.indices)

    def screen(self, moving=None):
        return RingScreen(self.coordinates, self.bonds, self.rings, moving=moving, labels=self.labels)

    def repair(self, before):
        """Never move markers or backbone; accept only locally improving moves."""
        changes, repairs = {}, []
        keys = sorted({self.atom_residues[i] for hit in before if hit["status"] == "penetration" for i in hit["bond"]})
        for key in keys:
            residue = self.residues[key]
            if residue.name not in FLEXIBLE_NATIVE:
                continue
            atoms = [SimpleNamespace(name=n, coord=np.array(_xyz(row))) for n, row in residue.atoms.items()]
            torsion = SidechainTorsions(atoms, residue_templates()[residue.name])
            if not torsion.rotations:
                continue
            names = [n for n in torsion.names if n not in BACKBONE and not is_hydrogen(n)]
            moving = np.array([self.indices[(*key, n)] for n in names], dtype=int)
            mapped = [torsion.index[n] for n in names]
            fixed = np.array(sorted(set(range(len(self.coordinates))) - set(moving)), dtype=int)
            tree = cKDTree(self.coordinates[fixed])
            screen = self.screen(moving)
            baseline_hits = screen.check(self.coordinates)
            identity = lambda hit: (tuple(hit["ring"]), tuple(hit["bond"]))
            original = {identity(h) for h in baseline_hits if h["status"] == "penetration"}
            if not original:
                continue
            excluded = {i: {i} | self.adjacency[i] | {k for j in self.adjacency[i] for k in self.adjacency[j]} for i in moving}
            radii = np.array([RADII.get(label["atom"][0], 1.7) for label in self.labels])

            def evaluate(angles):
                local = torsion.coordinates(angles)
                coords = self.coordinates.copy()
                coords[moving] = local[mapped]
                all_hits = screen.check(coords)
                hits = [h for h in all_hits if h["status"] == "penetration"]
                pairs = [(i, int(fixed[j])) for i, near in zip(moving, tree.query_ball_point(coords[moving], 3.5))
                         for j in near if fixed[j] not in excluded[i]]
                # Internal native-residue geometry is invariant under these torsions.
                pairs = np.array(pairs, dtype=int).reshape(-1, 2)
                distances = np.linalg.norm(coords[pairs[:, 0]]-coords[pairs[:, 1]], axis=1)
                overlap = np.maximum(0, .85 * radii[pairs].sum(axis=1) - distances)
                minimum = float(distances.min()) if len(distances) else 100.
                score = (1e7 * len(hits) + 1e3 * (len(all_hits)-len(hits))
                         + 20 * float(np.sum(overlap**2)) + .01 * float(np.sum(angles**2)))
                return score, hits, minimum, local

            start = np.zeros(len(torsion.rotations))
            initial = evaluate(start)
            best = []

            def objective(angles):
                metrics = evaluate(angles)
                new_ids = {identity(h) for h in metrics[1]}
                if (len(new_ids) < len(original) and new_ids <= original
                        and metrics[2] >= min(.25, initial[2]) - 1e-8
                        and (not best or metrics[0] < best[0][0])):
                    best[:] = [metrics]
                return metrics[0]

            if len(start) == 1:
                for angle in np.linspace(-np.pi, np.pi, 361):
                    objective(np.array([angle]))
            else:
                seed = int.from_bytes(hashlib.sha256(repr(key).encode()).digest()[:4], "little")
                differential_evolution(objective, [(-np.pi, np.pi)] * len(start), maxiter=200,
                                       popsize=8, tol=.002, polish=True, seed=seed, x0=start)
            if not best:
                continue
            score, hits, minimum, local = best[0]
            self.coordinates[moving] = local[mapped]
            changes[key] = {name: local[torsion.index[name]] for name in torsion.names if name not in BACKBONE}
            repairs.append(dict(model=key[0], chain=residue.address[1], residue=residue.name+residue.address[2],
                                penetrations_before=len(original), penetrations_after=len(hits),
                                minimum_nonbonded_A=minimum,
                                max_displacement_A=float(np.linalg.norm(local-torsion.coords, axis=1).max())))
        return changes, repairs

    def write_changes(self, changes):
        for mid, (path, lines) in self.documents.items():
            output = list(lines)
            for key, coordinates in changes.items():
                if key[0] != mid:
                    continue
                residue = self.residues[key]
                for index in residue.line_indices:
                    line = lines[index]
                    name = line[12:16].strip()
                    if name in coordinates:
                        assert name not in BACKBONE
                        block = "".join(f"{v:8.3f}" for v in coordinates[name])
                        if len(block) != 24 or not np.isfinite(coordinates[name]).all():
                            raise CrosslinkIntegrityError("Invalid ring-repair coordinates")
                        output[index] = line[:30] + block + line[54:]
            if output != lines:
                path.write_text("".join(output))


def validate_caps_rings(caps, network=None, repair=False, policy="warn", report_path=None):
    if policy not in {"warn", "error"}:
        raise ValueError("Ring penetration policy must be warn or error")
    geometry = CapsGeometry(caps, network)
    before = geometry.screen().check(geometry.coordinates)
    changes, repairs = geometry.repair(before) if repair else ({}, [])
    # Validate the rounded coordinates that will actually be written to PDB.
    for key, values in changes.items():
        for name, xyz in values.items():
            index = geometry.indices.get((*key, name))
            if index is not None:
                geometry.coordinates[index] = [float(f"{v:.3f}") for v in xyz]
    after = geometry.screen().check(geometry.coordinates)
    identity = lambda h: (tuple(h["ring"]), tuple(h["bond"]))
    existing = {identity(h) for h in before if h["status"] == "penetration"}
    remaining = [h for h in after if h["status"] == "penetration"]
    if any(identity(h) not in existing for h in remaining):
        raise CrosslinkIntegrityError("Ring repair introduced a new penetration after PDB rounding")
    report = dict(schema=1, policy=policy, rings_checked=len(geometry.rings),
                  unsupported_residues=sorted(geometry.unsupported), incomplete_rings=geometry.incomplete,
                  penetrations_before=len(existing), penetrations_after=len(remaining),
                  before=before, after=after, repairs=repairs,
                  status=("unresolved" if remaining else "incomplete" if geometry.incomplete or geometry.unsupported
                          else "ambiguous" if after else "clear"))
    if report_path is not None:
        Path(report_path).write_text(json.dumps(report, indent=2) + "\n")
    for hit in remaining:
        LOG.warning("Unresolved ring penetration: ring=%s, crossing bond=%s", hit["ring_atoms"], hit["bond_atoms"])
    if geometry.incomplete or geometry.unsupported:
        LOG.warning("Ring screen incomplete: %d incomplete rings; unsupported residues=%s", len(geometry.incomplete), sorted(geometry.unsupported))
    if remaining and policy == "error":
        raise CrosslinkIntegrityError(f"{len(remaining)} unresolved bond/ring penetrations; see ring_geometry_report.json")
    geometry.write_changes(changes)
    LOG.info("Ring screen: %d confirmed penetrations before, %d after; %d native sidechains repaired",
             len(existing), len(remaining), len(repairs))
    return report
