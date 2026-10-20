"""Nonperiodic heavy-bond/ring surface screening from residue bond graphs.

A confirmed penetration crosses both a centroid-fan surface and its best-fit
plane polygon, away from endpoints and the ring perimeter. Disagreement and
grazing contacts are reported as ambiguous, not silently called penetrations.
This is a geometric screen, not a topological linking invariant or MD engine.
"""
from collections import defaultdict
from functools import lru_cache
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.spatial import cKDTree

RTP = Path(__file__).resolve().parents[2] / "data/topology/amber99sb-star-ildnp.ff/aminoacids.rtp"


@lru_cache(maxsize=4)
def residue_templates(path=RTP):
    templates = defaultdict(lambda: {"atoms": {}, "bonds": []})
    name = section = None
    for line in Path(path).read_text().splitlines():
        text = line.split(";", 1)[0].strip()
        if not text or text.startswith("#"):
            continue
        if text.startswith("["):
            key = text.strip("[] ")
            if key in {"atoms", "bonds", "angles", "dihedrals", "impropers", "exclusions"}:
                section = key
            else:
                name, section = key, None
        elif section == "atoms":
            fields = text.split()
            templates[name]["atoms"][fields[0]] = fields[1]
        elif section == "bonds":
            templates[name]["bonds"].append(tuple(text.split()[:2]))
    return dict(templates)


def is_hydrogen(name):
    return name.lstrip("0123456789").startswith("H")


@lru_cache(maxsize=256)
def template_rings(name, path=RTP):
    template = residue_templates(path).get("HIE" if name == "HIS" else name)
    if template is None:
        return ()
    graph = nx.Graph((a, b) for a, b in template["bonds"]
                     if a in template["atoms"] and b in template["atoms"]
                     and not is_hydrogen(a) and not is_hydrogen(b))
    rings = []
    for nodes in nx.minimum_cycle_basis(graph):
        subgraph = graph.subgraph(nodes)
        if any(subgraph.degree(n) != 2 for n in nodes):
            raise ValueError(f"Non-simple minimum ring in RTP residue {name}: {nodes}")
        ordered, previous = [min(nodes)], None
        while len(ordered) < len(nodes):
            following = next(n for n in sorted(subgraph[ordered[-1]])
                             if n != previous and n not in ordered)
            previous = ordered[-1]
            ordered.append(following)
        rings.append(tuple(ordered))
    return tuple(sorted(rings))


def residue_ring_indices(residues, indices, path=RTP):
    """Residues are (identity, resname); indices map (*identity, atom) to IDs."""
    rings, incomplete = [], []
    for identity, name in residues:
        for names in template_rings(name, path):
            keys = [(*identity, atom) for atom in names]
            if all(key in indices for key in keys):
                rings.append(tuple(indices[key] for key in keys))
            else:
                incomplete.append(dict(residue=list(identity), resname=name,
                                       missing=[key[-1] for key in keys if key not in indices]))
    return rings, incomplete


def _surface_hits(poly, segments, tolerance):
    """Vectorized segment/triangle intersections, followed by plane confirmation."""
    center = poly.mean(axis=0)
    _, singular, axes = np.linalg.svd(poly - center)
    if singular[1] < 1e-8:
        return [(i, {"status": "ambiguous", "reason": "degenerate_ring"})
                for i in range(len(segments))]
    normal = axes[-1]
    start, end = segments[:, 0], segments[:, 1]
    direction = end - start
    e1, e2 = poly - center, np.roll(poly, -1, axis=0) - center
    h = np.cross(direction[:, None], e2[None])
    determinant = np.sum(e1[None] * h, axis=2)
    inv = np.divide(1., determinant, out=np.zeros_like(determinant),
                    where=abs(determinant) > 1e-10)
    s = start - center
    u = np.sum(s[:, None] * h, axis=2) * inv
    q = np.cross(s[:, None], e1[None])
    v = np.sum(direction[:, None] * q, axis=2) * inv
    t = np.sum(e2[None] * q, axis=2) * inv
    valid = ((abs(determinant) > 1e-10) & (u >= -1e-8) & (v >= -1e-8)
             & (u + v <= 1 + 1e-8) & (t >= 0) & (t <= 1))
    signed = (segments - center) @ normal
    denom = signed[:, 0] - signed[:, 1]
    plane_t = np.divide(signed[:, 0], denom, out=np.full(len(start), -1.), where=abs(denom) > 1e-10)
    coplanar = np.max(abs(signed), axis=1) <= tolerance
    plane_t[coplanar] = np.clip(
        np.sum((center-start[coplanar])*direction[coplanar], axis=1)
        / np.maximum(np.sum(direction[coplanar]**2, axis=1), 1e-12), 0, 1)
    candidates = np.flatnonzero(valid.any(axis=1) | ((plane_t >= 0) & (plane_t <= 1)))
    answer = []
    for i in candidates:
        plane_point = start[i] + plane_t[i] * direction[i]
        projected = (poly - plane_point) @ axes[:2].T
        following = np.roll(projected, -1, axis=0)
        winding = np.arctan2(projected[:, 0] * following[:, 1] - projected[:, 1] * following[:, 0],
                             np.sum(projected * following, axis=1)).sum()
        inside = 0 <= plane_t[i] <= 1 and abs(winding) > np.pi
        hits = sorted(t[i, valid[i]])
        unique = [x for j, x in enumerate(hits) if j == 0 or x - hits[j-1] > 1e-6]
        if not unique and not inside:
            continue
        hit_t = unique[0] if unique else plane_t[i]
        point = start[i] + hit_t * direction[i]
        edges = np.roll(poly, -1, axis=0) - poly
        along = np.clip(np.sum((point - poly) * edges, axis=1)
                        / np.maximum(np.sum(edges * edges, axis=1), 1e-12), 0, 1)
        margin = float(np.linalg.norm(point - poly - along[:, None] * edges, axis=1).min())
        endpoint_margin = float(min(hit_t, 1-hit_t) * np.linalg.norm(direction[i]))
        confirmed = (len(unique) == 1 and inside and margin > tolerance
                     and endpoint_margin > tolerance and min(abs(signed[i])) > tolerance)
        answer.append((int(i), dict(status="penetration" if confirmed else "ambiguous",
                                   fan_intersections=[float(x) for x in unique],
                                   plane_inside=bool(inside), edge_clearance_A=margin,
                                   endpoint_signed_A=signed[i].tolist(),
                                   ring_max_out_of_plane_A=float(abs((poly-center) @ normal).max()))))
    return answer


class RingScreen:
    """Cache static spatial indices; only rescreen surfaces affected by a move."""
    def __init__(self, coordinates, bonds, rings, moving=None, labels=None, tolerance=0.02):
        coordinates = np.asarray(coordinates, dtype=float)
        if coordinates.ndim != 2 or coordinates.shape[1] != 3:
            raise ValueError("Ring coordinates must have shape (n, 3)")
        if not np.isfinite(coordinates).all():
            raise ValueError("Nonfinite coordinates in ring screen")
        self.bonds = np.array(sorted({tuple(sorted(b)) for b in bonds}), dtype=int).reshape(-1, 2)
        self.rings = tuple(tuple(r) for r in rings)
        self.labels, self.tolerance = labels, tolerance
        mask = np.zeros(len(coordinates), dtype=bool)
        if moving is not None:
            mask[list(moving)] = True
        self.dynamic_bonds = np.flatnonzero(mask[self.bonds].any(axis=1))
        self.fixed_bonds = np.flatnonzero(~mask[self.bonds].any(axis=1))
        segments = coordinates[self.bonds[self.fixed_bonds]]
        self.bond_tree = cKDTree(segments.mean(axis=1))
        self.max_half_bond = float(np.linalg.norm(segments[:, 1]-segments[:, 0], axis=1).max()/2) if len(segments) else 0
        self.dynamic_rings = [i for i, ring in enumerate(self.rings) if moving is None or mask[list(ring)].any()]
        self.dynamic_ring_set = set(self.dynamic_rings)
        self.fixed_rings = [i for i in range(len(self.rings)) if i not in self.dynamic_ring_set]
        centers = [coordinates[list(self.rings[i])].mean(axis=0) for i in self.fixed_rings]
        self.ring_tree = cKDTree(np.array(centers).reshape(-1, 3))
        self.max_ring_radius = max((float(np.linalg.norm(coordinates[list(self.rings[i])] - c, axis=1).max())
                                    for i, c in zip(self.fixed_rings, centers)), default=0)

    def check(self, coordinates):
        coordinates = np.asarray(coordinates)
        if not np.isfinite(coordinates).all():
            raise ValueError("Nonfinite coordinates in ring screen")
        work = {i: set(self.dynamic_bonds) for i in self.dynamic_rings}
        for bond in self.dynamic_bonds:
            segment = coordinates[self.bonds[bond]]
            radius = self.max_ring_radius + np.linalg.norm(segment[1]-segment[0])/2 + self.tolerance
            for j in self.ring_tree.query_ball_point(segment.mean(axis=0), radius):
                work.setdefault(self.fixed_rings[j], set()).add(bond)
        hits = []
        for i, candidates in sorted(work.items()):
            ring = self.rings[i]
            poly = coordinates[list(ring)]
            center = poly.mean(axis=0)
            radius = np.linalg.norm(poly-center, axis=1).max() + self.max_half_bond + self.tolerance
            if i in self.dynamic_ring_set:
                candidates.update(self.fixed_bonds[self.bond_tree.query_ball_point(center, radius)])
            ring_set = set(ring)
            candidates = sorted(b for b in candidates if not ring_set.intersection(self.bonds[b]))
            if not candidates:
                continue
            for j, result in _surface_hits(poly, coordinates[self.bonds[candidates]], self.tolerance):
                bond = tuple(int(x) for x in self.bonds[candidates[j]])
                result.update(ring=list(ring), bond=list(bond))
                if self.labels is not None:
                    result.update(ring_atoms=[self.labels[a] for a in ring], bond_atoms=[self.labels[a] for a in bond])
                hits.append(result)
        return hits
