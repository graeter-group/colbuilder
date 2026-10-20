"""Sterically screened glucosepane closure using bond-preserving torsions.

This is a heavy-atom geometry preconditioner, not a force-field minimizer.
It never changes the topology or adds distance-based exclusions.
"""

from collections import defaultdict, deque

import numpy as np
from scipy.optimize import differential_evolution
from scipy.spatial import cKDTree

from colbuilder.core.utils.constants import MAX_DIVALENT_DISTANCE
from colbuilder.core.utils.logger import setup_logger
from colbuilder.core.utils.ring_geometry import (
    RTP, RingScreen, is_hydrogen, residue_ring_indices, residue_templates,
)
from colbuilder.core.utils.steric_policy import (
    CATASTROPHIC_HEAVY_A,
    WARN_14_A,
    WARN_HEAVY_A,
)

LOG = setup_logger(__name__)
BACKBONE = {"N", "CA", "C", "O", "OXT", "H", "H1", "H2", "H3", "HA"}
RADII = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "H": 1.20}
TARGET_MAX_FORMING_A = 2.5
MAX_FORMING_A = MAX_DIVALENT_DISTANCE  # Existing divalent connectivity limit: 5 A.


class StericOptimizationError(ValueError):
    """No closed, gross-clash-free candidate was found; do not export it."""


def get_residue(structure, chain, number):
    matches = [r for r in structure[0][chain] if r.id[1] == int(number)]
    if len(matches) != 1:
        raise StericOptimizationError(f"Ambiguous or absent residue {number}.{chain}")
    return matches[0]


class SidechainTorsions:
    """Cut only acyclic, aliphatic sidechain bonds; retain rings and backbone."""

    def __init__(self, residue, template):
        self.names = [a.name for a in residue]
        self.coords = np.array([a.coord for a in residue], dtype=float)
        self.index = {name: i for i, name in enumerate(self.names)}
        adjacency = defaultdict(set)
        for a, b in template["bonds"]:
            if a in self.index and b in self.index:
                adjacency[a].add(b)
                adjacency[b].add(a)
        distances = {"CA": 0}
        queue = deque(["CA"])
        while queue:
            a = queue.popleft()
            for b in sorted(adjacency[a]):
                if b not in distances:
                    distances[b] = distances[a] + 1
                    queue.append(b)
        self.rotations = []
        for a, b in sorted(
            template["bonds"], key=lambda pair: min(distances.get(n, 999) for n in pair)
        ):
            if (
                a not in self.index
                or b not in self.index
                or is_hydrogen(a)
                or is_hydrogen(b)
            ):
                continue
            if distances.get(a, 999) > distances.get(b, 999):
                a, b = b, a
            if (a != "CA" and a in BACKBONE) or b in BACKBONE:
                continue
            # Do not rotate amide/conjugated bonds or any bond within a ring.
            if "CT" not in (template["atoms"][a], template["atoms"][b]):
                continue
            distal, queue = {b}, [b]
            while queue:
                node = queue.pop()
                for neighbor in adjacency[node]:
                    if {node, neighbor} == {a, b} or neighbor in distal:
                        continue
                    distal.add(neighbor)
                    queue.append(neighbor)
            if a in distal or distal & BACKBONE:
                continue
            moving = sorted(distal - {b})
            if not any(not is_hydrogen(n) for n in moving):
                continue
            self.rotations.append(
                (self.index[a], self.index[b], [self.index[n] for n in moving])
            )

    def coordinates(self, angles):
        coords = self.coords.copy()
        for (a, b, moving), angle in zip(self.rotations, angles):
            axis = coords[b] - coords[a]
            length = np.linalg.norm(axis)
            if length < 1e-8:
                raise StericOptimizationError("Zero-length torsion axis")
            axis /= length
            vectors = coords[moving] - coords[b]
            cosine, sine = np.cos(angle), np.sin(angle)
            coords[moving] = (
                coords[b]
                + vectors * cosine
                + np.cross(axis, vectors) * sine
                + np.outer(vectors @ axis, axis) * (1 - cosine)
            )
        return coords


class GlucosepaneGeometry:
    """Score the selected residues against both complete translated helices.

    Every occurrence of a changed residue moves identically: this matches the
    later replay into the seed and avoids optimizing inconsistent copies.
    """

    def __init__(self, structures, crosslink, rtp_path=RTP, flexible_neighbors=()):
        self.structures = structures
        self.crosslink = crosslink
        self.roles = [crosslink["R1"], crosslink["R2"]]
        self.templates = residue_templates(rtp_path)
        if {r["type"] for r in self.roles} != {"AGS", "LGX"}:
            raise StericOptimizationError(
                "Steric closure currently supports AGS/LGX only"
            )
        for sid, chain, number in flexible_neighbors:
            residue = get_residue(structures[sid], chain, number)
            self.roles.append(
                dict(
                    structure_id=sid,
                    chain=chain,
                    position=str(number),
                    type=residue.resname,
                )
            )
        identities = [(r["chain"], int(r["position"])) for r in self.roles]
        if len(set(identities)) != len(identities):
            raise StericOptimizationError(
                "Repeated flexible residue in glucosepane optimization"
            )
        self.torsions = []
        for role in self.roles:
            residue = get_residue(
                structures[role["structure_id"]], role["chain"], role["position"]
            )
            self.torsions.append(
                SidechainTorsions(residue, self.templates[residue.resname])
            )
        self.sizes = [len(t.rotations) for t in self.torsions]
        self.dimension = sum(self.sizes)
        if not self.dimension:
            raise StericOptimizationError("No valid sidechain torsions available")

        self.labels, coordinates, self.indices = [], [], {}
        ring_residues = []
        self.adjacency = defaultdict(set)
        for sid in ("copy1", "copy2"):
            for chain in structures[sid][0]:
                residues = list(chain)
                for ri, residue in enumerate(residues):
                    ring_residues.append(((sid, chain.id, residue.id[1]), residue.resname))
                    for atom in residue:
                        if is_hydrogen(atom.name):
                            continue
                        key = (sid, chain.id, residue.id[1], atom.name)
                        self.indices[key] = len(coordinates)
                        self.labels.append(key)
                        coordinates.append(atom.coord)
                for ri, residue in enumerate(residues):
                    # HIS is protonated later by pdb2gmx; its heavy-atom graph
                    # is shared by the AMBER histidine protonation templates.
                    template = self.templates.get(
                        "HIE" if residue.resname == "HIS" else residue.resname
                    )
                    if template is None:
                        raise StericOptimizationError(
                            f"Missing residue template: {residue.resname}"
                        )
                    for names in template["bonds"]:
                        ends = []
                        for name in names:
                            other = ri
                            if name.startswith(("-", "+")):
                                other += -1 if name[0] == "-" else 1
                                name = name[1:]
                                if not 0 <= other < len(residues):
                                    break
                                if abs(residues[other].id[1] - residue.id[1]) != 1:
                                    break
                            index = self.indices.get(
                                (sid, chain.id, residues[other].id[1], name)
                            )
                            if index is None:
                                break
                            ends.append(index)
                        if len(ends) == 2:
                            a, b = ends
                            self.adjacency[a].add(b)
                            self.adjacency[b].add(a)
        self.original = np.array(coordinates, dtype=float)
        self.ports = [
            self.indices[(r["structure_id"], r["chain"], int(r["position"]), r["atom"])]
            for r in self.roles[:2]
        ]
        a, b = self.ports
        self.adjacency[a].add(b)
        self.adjacency[b].add(a)

        self.instances = []
        moving = set()
        for role_index, (role, torsion) in enumerate(zip(self.roles, self.torsions)):
            for sid in ("copy1", "copy2"):
                residue = get_residue(structures[sid], role["chain"], role["position"])
                shift = (
                    np.array(residue["CA"].coord) - torsion.coords[torsion.index["CA"]]
                )
                names = [n for n in torsion.names if not is_hydrogen(n)]
                indices = [
                    self.indices[(sid, role["chain"], int(role["position"]), n)]
                    for n in names
                ]
                mapped = [torsion.index[n] for n in names]
                if not np.allclose(
                    self.original[indices], torsion.coords[mapped] + shift, atol=0.025
                ):
                    raise StericOptimizationError(
                        "Crosslink copies are not consistent translations"
                    )
                self.instances.append((role_index, indices, mapped, shift))
                moving.update(i for n, i in zip(names, indices) if n not in BACKBONE)
        self.moving = np.array(sorted(moving), dtype=int)
        self.fixed = np.array(sorted(set(range(len(coordinates))) - moving), dtype=int)
        self.fixed_tree = cKDTree(self.original[self.fixed])
        self.excluded, self.pairs14 = {}, {}
        for i in self.moving:
            visited, frontier = {i}, {i}
            shells = []
            for _ in range(3):
                frontier = {k for j in frontier for k in self.adjacency[j]} - visited
                visited |= frontier
                shells.append(frontier)
            self.excluded[i] = {i} | shells[0] | shells[1]
            self.pairs14[i] = shells[2]
        self.radii = np.array(
            [RADII.get(label[3].lstrip("0123456789")[0], 1.7) for label in self.labels]
        )
        self.dynamic_pairs = np.array(
            [
                (i, j)
                for n, i in enumerate(self.moving)
                for j in self.moving[n + 1 :]
                if j not in self.excluded[i]
            ],
            dtype=int,
        )
        self.joining_angles = [(n, a, b) for n in self.adjacency[a] - {b}]
        self.joining_angles += [(a, b, n) for n in self.adjacency[b] - {a}]
        rings, self.incomplete_rings = residue_ring_indices(ring_residues, self.indices, rtp_path)
        self.ring_screen = RingScreen(
            self.original, [(i, j) for i, neighbors in self.adjacency.items() for j in neighbors if i < j],
            rings, moving=self.moving, labels=self.labels,
        )

    def coordinates(self, angles):
        split = np.cumsum([0] + self.sizes)
        residues = [
            t.coordinates(angles[split[k] : split[k + 1]])
            for k, t in enumerate(self.torsions)
        ]
        coords = self.original.copy()
        for k, indices, mapped, shift in self.instances:
            coords[indices] = residues[k][mapped] + shift
        return coords, residues

    def evaluate(self, angles):
        coords, _ = self.coordinates(angles)
        candidates = self.fixed_tree.query_ball_point(coords[self.moving], 3.5)
        pairs = [
            (i, int(self.fixed[j]))
            for i, neighbors in zip(self.moving, candidates)
            for j in neighbors
            if int(self.fixed[j]) not in self.excluded[i]
        ]
        pairs.extend(map(tuple, self.dynamic_pairs))
        pairs = np.array(pairs, dtype=int).reshape(-1, 2)
        distances = np.linalg.norm(coords[pairs[:, 0]] - coords[pairs[:, 1]], axis=1)
        is14 = np.array([j in self.pairs14[i] for i, j in pairs])
        limits = (self.radii[pairs[:, 0]] + self.radii[pairs[:, 1]]) * np.where(
            is14, 0.72, 0.85
        )
        overlap = np.maximum(limits - distances, 0)
        steric = float(np.sum(overlap**2 * np.where(is14, 0.25, 1.0)))
        forming = float(np.linalg.norm(coords[self.ports[0]] - coords[self.ports[1]]))
        angle_penalty = 0.0
        for a, b, c in self.joining_angles:
            v, w = coords[a] - coords[b], coords[c] - coords[b]
            cosine = np.dot(v, w) / max(np.linalg.norm(v) * np.linalg.norm(w), 1e-10)
            angle_penalty += (cosine + 1 / 3) ** 2
        ordinary = np.flatnonzero(~is14)
        worst = ordinary[np.argmin(distances[ordinary])] if len(ordinary) else None
        minimum = float(distances[worst]) if worst is not None else float("inf")
        # Keep closure as a search target, but reserve the hard steric penalty
        # for nearly coincident atoms. Compressed contacts retain the soft score.
        closure_violation = (
            max(forming - (TARGET_MAX_FORMING_A - 0.05), 0) ** 2
            + max(1.25 - forming, 0) ** 2
        )
        contact_violation = float(
            np.sum(np.maximum(CATASTROPHIC_HEAVY_A + 0.05 - distances, 0) ** 2)
        )
        ring_hits = self.ring_screen.check(coords)
        penetrations = [hit for hit in ring_hits if hit["status"] == "penetration"]
        return {
            "score": (
                25 * (forming - 1.5) ** 2
                + 20 * steric
                + 2 * angle_penalty
                + 10000 * (closure_violation + contact_violation)
                + 1e7 * len(penetrations)
                + 1e3 * (len(ring_hits) - len(penetrations))
            ),
            "forming_A": forming,
            "minimum_nonbonded_A": minimum,
            "minimum_14_A": (
                float(distances[is14].min()) if np.any(is14) else float("inf")
            ),
            "steric_penalty": steric,
            "ring_penetrations": len(penetrations),
            "ring_contacts": ring_hits,
            "worst_contact": (
                [self.labels[i] for i in pairs[worst]] if worst is not None else []
            ),
        }

    @staticmethod
    def acceptable(result):
        """Permit strained pre-minimization candidates, not catastrophic overlaps."""
        return bool(
            1.2 <= result["forming_A"] <= MAX_FORMING_A
            and result["minimum_nonbonded_A"] >= CATASTROPHIC_HEAVY_A
            and result["minimum_14_A"] >= CATASTROPHIC_HEAVY_A
            and np.isfinite(result["score"])
            and result.get("ring_penetrations", 0) == 0
        )

    @staticmethod
    def quality_warnings(result):
        warnings = []
        if any(hit["status"] == "ambiguous" for hit in result.get("ring_contacts", [])):
            warnings.append("ambiguous bond/ring contact; inspect ring_contacts in geometry report")
        if result["forming_A"] > TARGET_MAX_FORMING_A:
            warnings.append(
                f"forming distance {result['forming_A']:.3f} A exceeds the "
                f"{TARGET_MAX_FORMING_A:.1f} A optimization target "
                f"(connectivity limit {MAX_FORMING_A:.1f} A)"
            )
        if result["minimum_nonbonded_A"] < WARN_HEAVY_A:
            warnings.append(
                f"compressed nonbonded heavy contact {result['minimum_nonbonded_A']:.3f} A"
            )
        if result["minimum_14_A"] < WARN_14_A:
            warnings.append(f"compressed 1-4 heavy contact {result['minimum_14_A']:.3f} A")
        return warnings

    def apply(self, angles):
        _, residues = self.coordinates(angles)
        for role, torsion, coordinates in zip(self.roles, self.torsions, residues):
            for structure in self.structures.values():
                residue = get_residue(structure, role["chain"], role["position"])
                shift = (
                    np.array(residue["CA"].coord) - torsion.coords[torsion.index["CA"]]
                )
                for name, coordinate in zip(torsion.names, coordinates):
                    # Preserve backbone coordinates exactly, including roundoff.
                    if name not in BACKBONE:
                        residue[name].coord = coordinate + shift


def optimize_glucosepane(structures, crosslink, maxiter=500):
    neighbors, previous_angles = [], np.array([])
    for expansion in range(3):
        geometry = GlucosepaneGeometry(
            structures, crosslink, flexible_neighbors=neighbors
        )
        start = np.zeros(geometry.dimension)
        start[: len(previous_angles)] = previous_angles
        LOG.info(
            "Glucosepane torsional/steric optimization: %s; flexible neighbors=%s",
            geometry.evaluate(start),
            neighbors,
        )
        feasible = []

        def objective(angles):
            metrics = geometry.evaluate(angles)
            if geometry.acceptable(metrics) and (
                not feasible or metrics["score"] < feasible[0][1]["score"]
            ):
                feasible[:] = [(angles.copy(), metrics)]
            return metrics["score"]

        result = differential_evolution(
            objective,
            [(-np.pi, np.pi)] * geometry.dimension,
            maxiter=maxiter,
            popsize=10,
            tol=0.002,
            polish=True,
            x0=start,
            seed=int(np.random.randint(0, 2**31 - 1)),
        )
        final = geometry.evaluate(result.x)
        if feasible and (
            not geometry.acceptable(final) or feasible[0][1]["score"] < final["score"]
        ):
            result.x, final = feasible[0]
        LOG.info("Glucosepane torsional/steric result: %s", final)
        if geometry.acceptable(final):
            break
        existing = {(r["chain"], int(r["position"])) for r in geometry.roles}
        extra = None
        ring_neighbors = [label for hit in final.get("ring_contacts", [])
                          if hit["status"] == "penetration" for label in hit["bond_atoms"]]
        for sid, chain, number, atom in ring_neighbors + list(final.get("worst_contact", [])):
            residue = get_residue(structures[sid], chain, number)
            if (chain, number) in existing or atom in BACKBONE:
                continue
            if residue.resname not in {
                "LYS",
                "ARG",
                "GLU",
                "GLN",
                "ASP",
                "ASN",
                "SER",
                "THR",
                "MET",
                "LEU",
                "ILE",
                "VAL",
            }:
                continue
            if SidechainTorsions(
                residue, geometry.templates[residue.resname]
            ).rotations:
                extra = (sid, chain, number)
                break
        if extra is None or expansion == 2:
            raise StericOptimizationError(
                f"Glucosepane geometry rejected: {final}; flexible neighbors={neighbors}"
            )
        neighbors.append(extra)
        previous_angles = result.x
    snapshot = [
        (atom, atom.coord.copy())
        for structure in structures.values()
        for role in geometry.roles
        for atom in get_residue(structure, role["chain"], role["position"])
    ]
    try:
        geometry.apply(result.x)
        # Test the coordinates that will actually be copied back to the seed.
        exported = GlucosepaneGeometry(
            structures, crosslink, flexible_neighbors=neighbors
        )
        check = exported.evaluate(np.zeros(exported.dimension))
        if not exported.acceptable(check):
            raise StericOptimizationError(
                f"Glucosepane replay validation failed: {check}"
            )
    except Exception:
        for atom, coordinate in snapshot:
            atom.coord = coordinate
        raise
    final["quality_warnings"] = exported.quality_warnings(check)
    if final["quality_warnings"]:
        LOG.warning(
            "Glucosepane geometry accepted with warnings: %s. "
            "Inspect and minimize before MD; this is not a converged structure.",
            "; ".join(final["quality_warnings"]),
        )
    final["flexible_neighbors"] = neighbors
    return final
