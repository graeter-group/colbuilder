"""Chemical crosslink identities shared by replacement and final topology.

Lattice-growth contacts are not a covalent graph. Resolve complete chemical
entities from the active caps, then retain their identities through mutation.
Long terminal bonds require an unambiguous configured recipe in a local caps
neighborhood; proximity alone never creates a long bond.
"""

from collections import Counter
from dataclasses import dataclass
import json
import math
from pathlib import Path
import random
import time

import numpy as np

from .backbone import read_pdb_residues
from .crosslink import Crosslink, read_crosslink
from ..utils.logger import setup_logger

LOG = setup_logger(__name__)
TOPOLOGY_CUTOFF = 5.0
TRIVALENT = {
    "PYD": (("LYX", "C12", "C13"), ("LY2", "CB"), ("LY3", "CG")),
    "DPD": (("LXX", "C12", "C13"), ("LX2", "CB"), ("LX3", "CG")),
    "PYL": (("LXY", "C12", "C13"), ("L2Y", "CG"), ("L3Y", "CG")),
    "DPL": (("LYY", "C12", "C13"), ("L2X", "CG"), ("L3X", "CG")),
}
DIVALENT = {
    "GCP": (("LGX", "CE"), ("AGS", "NZ")),
    "MOLD": (("LZD", "CE"), ("LZS", "NZ1")),
    "PENT": (("LPS", "CE"), ("APD", "NZ")),
    "HLKNL": (("L4Y", "CE"), ("L5Y", "NZ")),
    "LKNL": (("L4X", "CE"), ("L5X", "NZ")),
    "HLNL": (("LY4", "CE"), ("LY5", "NZ")),
    "LNL": (("LX4", "CE"), ("LX5", "NZ")),
}
MARKER_NAMES = {spec[0] for rule in (*TRIVALENT.values(), *DIVALENT.values()) for spec in rule}
AGE_NAMES = {"GCP", "MOLD", "PENT"}


class CrosslinkIntegrityError(ValueError):
    """The intended crosslink chemistry is ambiguous or incomplete."""


@dataclass
class Marker:
    key: tuple  # model id, residue id, chain, residue name
    atoms: dict

    def port(self, atom):
        if atom not in self.atoms:
            raise CrosslinkIntegrityError(f"Missing forming atom {self.key}:{atom}")
        return self.key + (atom,)


@dataclass(frozen=True)
class Entity:
    kind: str
    scope: str
    markers: tuple
    bonds: tuple

    @property
    def identity(self):
        return self.kind, tuple(sorted(self.markers))


@dataclass
class Network:
    markers: dict
    entities: tuple
    replacement_report: dict = None

    @property
    def counts(self):
        return dict(Counter(f"{e.kind}:{e.scope}" for e in self.entities))

    def pairs(self, model_ids):
        ids = set(model_ids)
        pairs = []
        for entity in self.entities:
            for a, b in entity.bonds:
                if a[0] not in ids and b[0] not in ids:
                    continue
                if not {a[0], b[0]} <= ids:
                    raise CrosslinkIntegrityError(
                        f"Crosslink split across topology groups: {a}, {b}"
                    )
                ends = []
                for port in (a, b):
                    mid, resid, chain, name, atom = port
                    ends.append(
                        Crosslink(
                            resid,
                            name,
                            chain,
                            self.markers[port[:4]].atoms[atom],
                            "T" if entity.kind in TRIVALENT else "D",
                            mid,
                            atom,
                        )
                    )
                pairs.append(tuple(ends))
        return pairs

    def write(self, path):
        payload = {
            "schema": 1,
            "counts": self.counts,
            "entities": [
                {
                    "kind": e.kind,
                    "scope": e.scope,
                    "markers": e.markers,
                    "bonds": e.bonds,
                    "distances_A": [distance(self.markers, a, b) for a, b in e.bonds],
                }
                for e in self.entities
            ],
        }
        if self.replacement_report is not None:
            payload["replacement"] = self.replacement_report
        Path(path).write_text(json.dumps(payload, indent=2) + "\n")


def load_markers(caps):
    markers = {}
    for mid, path in sorted(caps.items()):
        for residue in read_pdb_residues(
            Path(path).read_text().splitlines(keepends=True), str(path)
        ):
            if residue.name not in MARKER_NAMES:
                continue
            _, chain, number, insertion = residue.address
            key = (int(mid), number + insertion, chain, residue.name)
            if key in markers:
                raise CrosslinkIntegrityError(f"Ambiguous marker identity in {path}: {key}")
            atoms = {
                name: tuple(float(line[i : i + 8]) for i in (30, 38, 46))
                for name, line in residue.atoms.items()
            }
            if not all(math.isfinite(x) for xyz in atoms.values() for x in xyz):
                raise CrosslinkIntegrityError(f"Non-finite marker coordinates: {key}")
            markers[key] = Marker(key, atoms)
    return markers


def active_caps(system, root):
    """Select exactly one typed caps file for each active model; never stale ids."""
    result = {}
    for mid in system.get_models():
        model = system.get_model(mid)
        path = Path(root) / str(model.type) / f"{int(mid)}.caps.pdb"
        if not path.is_file():
            raise CrosslinkIntegrityError(f"Missing active caps file: {path}")
        result[int(mid)] = path
    return result


def distance(markers, a, b):
    return math.dist(markers[a[:4]].atoms[a[4]], markers[b[:4]].atoms[b[4]])


def _recipe(value):
    if not value:
        return frozenset()
    try:
        return frozenset(tuple(token.strip().rsplit(".", 1)) for token in value.split(" - "))
    except (AttributeError, ValueError) as exc:
        raise CrosslinkIntegrityError(f"Invalid terminal recipe: {value!r}") from exc


def resolve_network(caps, config=None, strict=True):
    markers = load_markers(caps)
    recipes = {
        f"enzymatic_{end}": _recipe(getattr(config, f"{end}_term_combination", None))
        for end in ("n", "c")
    }
    clouds = {}
    for marker in markers.values():
        clouds.setdefault(marker.key[0], []).extend(marker.atoms.values())
    neighbors = {mid: {mid} for mid in caps}
    mids = sorted(clouds)
    for i, a in enumerate(mids):
        for b in mids[i + 1 :]:
            delta = np.asarray(clouds[a])[:, None, :] - np.asarray(clouds[b])[None, :, :]
            if np.any(np.sum(delta * delta, axis=2) < TOPOLOGY_CUTOFF**2):
                neighbors[a].add(b)
                neighbors[b].add(a)
    entities, occupied = {}, {}

    def scope_for(kind, keys):
        if kind in AGE_NAMES:
            return "non_enzymatic"
        positions = frozenset((k[1], k[2]) for k in keys)
        hits = [scope for scope, recipe in recipes.items() if recipe and positions == recipe]
        if len(hits) > 1:
            raise CrosslinkIntegrityError(f"Overlapping N/C recipes: {positions}")
        return hits[0] if hits else "enzymatic"

    def add(kind, selected, rule):
        keys = tuple(m.key for m in selected)
        if kind in TRIVALENT:
            center, arm12, arm13 = selected
            bonds = (
                (center.port(rule[0][1]), arm12.port(rule[1][1])),
                (center.port(rule[0][2]), arm13.port(rule[2][1])),
            )
        else:
            bonds = ((selected[0].port(rule[0][1]), selected[1].port(rule[1][1])),)
        if any(a[0] == b[0] for a, b in bonds):
            raise CrosslinkIntegrityError(
                f"Expected inter-model crosslink, found intra-model: {keys}"
            )
        entity = Entity(kind, scope_for(kind, keys), keys, bonds)
        if entity.identity in entities:
            return
        conflicts = [key for key in keys if key in occupied]
        if conflicts:
            raise CrosslinkIntegrityError(f"Markers assigned to multiple crosslinks: {conflicts}")
        entities[entity.identity] = entity
        occupied.update((key, entity.identity) for key in keys)

    # A named terminal recipe in an unambiguous local caps neighborhood records
    # chemical intent even when optimization has not shortened the forming bond.
    for kind, rule in TRIVALENT.items():
        for center in markers.values():
            if center.key[3] != rule[0][0]:
                continue
            for recipe in recipes.values():
                if len(recipe) != 3 or (center.key[1], center.key[2]) not in recipe:
                    continue
                candidates = [
                    [
                        m
                        for m in markers.values()
                        if m.key[3] == spec[0]
                        and all(atom in m.atoms for atom in spec[1:])
                        and m.key[0] in neighbors[center.key[0]]
                        and (m.key[1], m.key[2]) in recipe
                    ]
                    for spec in rule
                ]
                if all(len(items) == 1 for items in candidates) and candidates[0][0] is center:
                    selected = [items[0] for items in candidates]
                    if frozenset((m.key[1], m.key[2]) for m in selected) == recipe:
                        add(kind, selected, rule)

    # Without recipe provenance, require unique atom-specific partners <5 A.
    available = [m for m in markers.values() if m.key not in occupied]
    for kind, rule in {**TRIVALENT, **DIVALENT}.items():
        for center in available:
            if center.key[3] != rule[0][0] or not all(atom in center.atoms for atom in rule[0][1:]):
                continue
            candidates = []
            for index, spec in enumerate(rule[1:], 1):
                a = center.port(rule[0][index] if kind in TRIVALENT else rule[0][1])
                matches = [
                    m
                    for m in available
                    if m.key[3] == spec[0]
                    and spec[1] in m.atoms
                    and m.key[0] != center.key[0]
                    and distance(markers, a, m.port(spec[1])) < TOPOLOGY_CUTOFF
                ]
                if len(matches) > 1:
                    raise CrosslinkIntegrityError(f"Ambiguous atom-specific partners for {a}")
                candidates.append(matches)
            if all(len(items) == 1 for items in candidates):
                add(kind, [center] + [items[0] for items in candidates], rule)
    unpaired = set(markers) - set(occupied)
    if strict and unpaired:
        raise CrosslinkIntegrityError(
            f"{len(unpaired)} unpaired markers in final caps; no complete, unambiguous "
            f"crosslink could be assigned: {sorted(unpaired)[:12]}"
        )
    network = Network(markers, tuple(sorted(entities.values(), key=lambda e: e.identity)))
    for e in network.entities:
        for a, b in e.bonds:
            d = distance(markers, a, b)
            if d >= TOPOLOGY_CUTOFF:
                LOG.warning(
                    "Configured %s bond retained at %.3f A: %s -- %s; geometry needs relaxation",
                    e.kind,
                    d,
                    a,
                    b,
                )
    return network


def select_replacements(network, ratio, scope, seed=None, mode="random"):
    if not math.isfinite(ratio) or not 0 <= ratio <= 100:
        raise CrosslinkIntegrityError("ratio_replace must be between 0 and 100")
    if scope not in {"enzymatic", "non_enzymatic", "all"}:
        raise CrosslinkIntegrityError(f"Unknown replacement scope: {scope}")
    if mode not in {"random", "preserve_attachment"}:
        raise CrosslinkIntegrityError(f"Unknown ratio_replace_mode: {mode}")
    seed = int(time.time()) if seed is None else seed
    if len({entity.identity for entity in network.entities}) != len(network.entities):
        raise CrosslinkIntegrityError("Duplicate chemical entities in replacement input")
    rng, buckets = random.Random(seed), {}
    for entity in network.entities:
        age = entity.scope == "non_enzymatic"
        if scope == "all" or (scope == "non_enzymatic") == age:
            buckets.setdefault(entity.scope, []).append(entity)
    selected, targets = [], {}
    for key, bucket in sorted(buckets.items()):
        bucket = sorted(bucket, key=lambda e: e.identity)
        if mode == "random":
            rng.shuffle(bucket)
        target = round(len(bucket) * ratio / 100)
        if mode == "random":
            selected.extend(bucket[:target])
        targets[key] = {
            "available": len(bucket),
            "removed": target,
            "remaining": len(bucket) - target,
        }
    from .replacement_policy import (
        attachment_report, constrained_selection, entity_record, validate_attachment_plan,
    )
    if mode == "preserve_attachment":
        selected = constrained_selection(network, buckets, targets, rng)
    selected_ids = {entity.identity for entity in selected}
    retained = [entity for entity in network.entities if entity.identity not in selected_ids]
    instructions = sorted(
        {
            f'{key[0]}.caps.pdb {"ARG" if key[3] in {"AGS", "APD"} else "LYS"} '
            f"{key[1]} {key[2]}"
            for e in selected
            for key in e.markers
        }
    )
    report = {
        "schema": 1,
        "mode": mode,
        "seed": seed,
        "ratio": ratio,
        "scope": scope,
        "targets": targets,
        "removed_entities": len(selected),
        "instructions": len(instructions),
        "removed_crosslinks": [entity_record(e) for e in selected],
        "retained_crosslinks": [entity_record(e) for e in retained],
        "attachment": attachment_report(network.entities, retained),
    }
    validate_attachment_plan(report, retained)
    return instructions, report


def replacement_survivors(network, instructions):
    """Resolve manual chain case exactly as the mutation guard, before mutation."""
    targets = set()
    for line in instructions:
        filename, target, number, chain = line.strip().strip("\"'").split()
        mid = int(Path(filename).name.split(".")[0])
        matches = [
            key
            for key in network.markers
            if key[:2] == (mid, number) and key[2].lower() == chain.lower()
        ]
        exact = [key for key in matches if key[2] == chain]
        matches = exact or matches
        if len(matches) > 1:
            raise CrosslinkIntegrityError(
                f"Ambiguous mutation target: {number}.{chain} in {filename}"
            )
        targets.update(matches)
    retained = []
    for entity in network.entities:
        affected = [key for key in entity.markers if key in targets]
        if affected and len(affected) != len(entity.markers):
            raise CrosslinkIntegrityError(f"Partial removal of {entity.kind}: {affected}")
        if not affected:
            retained.append(entity)
    return retained, set(network.markers) - targets


def after_replacement(network, instructions, caps):
    retained, expected = replacement_survivors(network, instructions)
    markers = load_markers(caps)
    if set(markers) != expected:
        raise CrosslinkIntegrityError("Replacement changed unexpected marker identities")
    occupied = {key for e in retained for key in e.markers}
    if set(markers) != occupied:
        raise CrosslinkIntegrityError(
            f"Unpaired markers remain after replacement: {sorted(set(markers) - occupied)}"
        )
    from .replacement_policy import validate_attachment_plan
    validate_attachment_plan(network.replacement_report, retained)
    return Network(markers, tuple(retained), network.replacement_report)


def refresh_network(network, caps):
    markers = load_markers(caps)
    if set(markers) != set(network.markers):
        raise CrosslinkIntegrityError("Final caps differ from the validated crosslink inventory")
    from .replacement_policy import validate_attachment_plan
    validate_attachment_plan(network.replacement_report, network.entities)
    return Network(markers, network.entities, network.replacement_report)


def install_network(system, network, caps, output=None):
    """Use one final graph for System, connect file and the AMBER bond ledger."""
    from .replacement_policy import validate_attachment_plan
    validate_attachment_plan(network.replacement_report, network.entities)
    adjacency = {int(mid): {int(mid)} for mid in system.get_models()}
    for entity in network.entities:
        for a, b in entity.bonds:
            adjacency[a[0]].add(b[0])
            adjacency[b[0]].add(a[0])
    remaining, groups = set(adjacency), []
    while remaining:
        stack, group = [min(remaining)], set()
        while stack:
            mid = stack.pop()
            if mid in group:
                continue
            group.add(mid)
            stack.extend(adjacency[mid] - group)
        remaining -= group
        groups.append(sorted(group))
    for group in groups:
        for mid in group:
            model = system.get_model(mid)
            model.connect = list(group)
            model.crosslink = read_crosslink(caps[mid])
            for marker in model.crosslink:
                marker.model_id = mid
    system.crosslink_network = network
    system.connect = {int(mid): list(system.get_model(mid).connect) for mid in system.get_models()}
    if output is not None:
        output = Path(output)
        output.mkdir(parents=True, exist_ok=True)
        network.write(output / "crosslink_network.json")
        if network.replacement_report is not None:
            (output / "replacement_report.json").write_text(
                json.dumps(network.replacement_report, indent=2) + "\n"
            )
        lines = [
            " ".join(f"{mid}.caps.pdb" for mid in group) + f" ; {system.get_model(group[0]).type}"
            for group in groups
        ]
        (output / "connect_from_colbuilder.txt").write_text("\n".join(lines) + "\n")
    return groups
