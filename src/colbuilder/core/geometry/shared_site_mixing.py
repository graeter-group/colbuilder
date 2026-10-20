"""Mix alternative chemical loci, retaining the union of additional crosslinks.

Ratio units are complete crosslinks within each positional family, not models
or marker atoms. A model's backbone comes from one variant; extra crosslinks
are transferred as rigid, complete entities onto that backbone.
"""

from collections import Counter
from dataclasses import dataclass
import json
import math
from pathlib import Path
import shutil
import tempfile

import numpy as np

from .backbone import BACKBONE, _xyz, read_pdb_residues
from .crosslink_network import (
    CrosslinkIntegrityError, Network, distance, install_network, load_markers,
    resolve_network,
)
from .unpaired_crosslinks import RESIDUE_TO_MUTATION
from .ring_validation import validate_caps_rings

OUTPUT_TYPE = "_mixed_sites"


@dataclass
class Locus:
    options: dict

    @property
    def positions(self):
        return frozenset(k[:3] for e in self.options.values() for k in e.markers)

    @property
    def family(self):
        return tuple(
            (label, entity.kind, tuple(sorted(k[1:] for k in entity.markers)))
            for label, entity in self.options.items()
        )


def match_loci(networks):
    """Match physical residue identities, never just nearby coordinates."""
    if len(networks) != 2:
        raise CrosslinkIntegrityError("shared_sites requires exactly two variants")
    a, b = networks
    left, right = networks[a].entities, networks[b].entities
    positions = lambda e: frozenset(k[:3] for k in e.markers)
    common, alternative, extra, used = [], [], [], set()
    for entity in left:
        hits = [other for other in right if positions(entity) & positions(other)]
        if not hits:
            extra.append(Locus({a: entity}))
            continue
        if len(hits) != 1 or hits[0].identity in used:
            raise CrosslinkIntegrityError(f"Ambiguous competing crosslinks at {positions(entity)}")
        other = hits[0]
        used.add(other.identity)
        locus = Locus({a: entity, b: other})
        if entity.identity == other.identity:
            if {frozenset(bond) for bond in entity.bonds} != {
                frozenset(bond) for bond in other.bonds
            }:
                raise CrosslinkIntegrityError("Common locus has conflicting forming bonds")
            common.append(locus)
        elif (entity.kind != other.kind
              and len(positions(entity) & positions(other)) >= 2
              and (positions(entity) <= positions(other) or positions(other) <= positions(entity))):
            alternative.append(locus)
        else:
            raise CrosslinkIntegrityError(
                f"Partially overlapping/incompatible crosslink sites: {locus.positions}"
            )
    extra.extend(Locus({b: e}) for e in right if e.identity not in used)
    return common, alternative, extra


def choose_alternatives(loci, ratios):
    """Exact per-family quotas, subject to one backbone variant per model."""
    labels = list(ratios)
    weights = [float(ratios[label]) for label in labels]
    if (len(labels) != 2 or not all(math.isfinite(v) and v >= 0 for v in weights)
            or not math.isclose(sum(weights), 100)):
        raise CrosslinkIntegrityError("shared_sites ratios must be two nonnegative percentages summing to 100")
    families = sorted({locus.family for locus in loci})
    family_index = {family: i for i, family in enumerate(families)}
    totals = Counter(locus.family for locus in loci)
    targets = []
    for family in families:
        raw = [totals[family] * w / 100 for w in weights]
        quota = [math.floor(v) for v in raw]
        for i in sorted(range(2), key=lambda i: (-(raw[i] - quota[i]), i))[:totals[family] - sum(quota)]:
            quota[i] += 1
        targets.append(quota[0])

    # Only competing loci constrain backbone ownership. Additional GCP edges
    # must not force all connected models to use the same variant.
    components = []
    for index, locus in enumerate(loci):
        models = {p[0] for p in locus.positions}
        hits = [i for i, (mids, _) in enumerate(components) if models & mids]
        members = {index}
        for i in reversed(hits):
            mids, indices = components.pop(i)
            models |= mids
            members |= indices
        components.append((models, members))
    components.sort(key=lambda c: tuple(sorted(c[0])))
    zero = (0,) * len(families)
    states = {zero: ()}
    for i, (_, members) in enumerate(components):
        counts = Counter(family_index[loci[j].family] for j in members)
        for state, chosen in list(states.items()):
            new = tuple(v + counts[k] for k, v in enumerate(state))
            if all(v <= targets[k] for k, v in enumerate(new)):
                states.setdefault(new, chosen + (i,))
        if len(states) > 200000:
            raise CrosslinkIntegrityError("Shared-site quota search exceeds 200000 states; split the mix into smaller site families")
    target = tuple(targets)
    if target not in states:
        raise CrosslinkIntegrityError(
            "Requested per-site ratios cannot be realized without splitting a backbone "
            f"component; first-variant targets={target}"
        )
    selected_components = set(states[target])
    selected, model_sources = {}, {}
    for i, (mids, indices) in enumerate(components):
        label = labels[0] if i in selected_components else labels[1]
        model_sources.update((mid, label) for mid in mids)
        selected.update((index, label) for index in indices)
    report = []
    for k, family in enumerate(families):
        achieved = Counter(selected[i] for i, locus in enumerate(loci) if locus.family == family)
        report.append({
            "recipe": family, "available": totals[family],
            "target": {labels[0]: targets[k], labels[1]: totals[family] - targets[k]},
            "achieved": {label: achieved[label] for label in labels},
        })
    return selected, model_sources, report


def _read_caps(variants):
    documents, lookup = {}, {}
    reference = next(iter(variants))
    ids = set(variants[reference])
    for label, caps in variants.items():
        if set(caps) != ids:
            raise CrosslinkIntegrityError("Mix variants have different clipped model IDs")
        for mid, path in sorted(caps.items()):
            lines = Path(path).read_text().splitlines(keepends=True)
            residues = read_pdb_residues(lines, str(path))
            keys = [(mid, r.address[2] + r.address[3], r.address[1]) for r in residues]
            if len(set(keys)) != len(keys) or len({r.address[0] for r in residues}) != 1:
                raise CrosslinkIntegrityError(f"Ambiguous residue IDs in {path}")
            documents[label, mid] = lines, residues
            lookup.update(((label, key), r) for key, r in zip(keys, residues))
            if label != reference:
                original = documents[reference, mid][1]
                if [(r.address, r.segment) for r in original] != [(r.address, r.segment) for r in residues]:
                    raise CrosslinkIntegrityError(f"Different backbone residue order/TER boundaries in model {mid}")
    return documents, lookup


def _rigid_fit(source, target, tolerance, description):
    source, target = np.asarray(source), np.asarray(target)
    source_center, target_center = source.mean(axis=0), target.mean(axis=0)
    u, singular, vt = np.linalg.svd((source - source_center).T @ (target - target_center))
    if singular[1] < 1e-8:
        raise CrosslinkIntegrityError(f"Degenerate backbone alignment for {description}")
    correction = np.eye(3)
    correction[2, 2] = np.linalg.det(u @ vt)
    rotation = u @ correction @ vt
    translation = target_center - source_center @ rotation
    residual = float(np.linalg.norm(source @ rotation + translation - target, axis=1).max())
    if residual > tolerance:
        raise CrosslinkIntegrityError(
            f"Cannot transfer {description} without deforming its geometry: backbone "
            f"alignment residual {residual:.4f} A exceeds {tolerance:.4f} A"
        )
    return rotation, translation, residual


def _render(lines, residues, replacements):
    """Preserve polymer records; remap serial-dependent records after grafting."""
    blocks, dropped = {}, set()
    for residue in residues:
        if residue.address in replacements:
            blocks[residue.line_indices[0]] = replacements[residue.address]
            dropped.update(residue.line_indices)
    output, serial, remap = [], 0, {}
    last_atom = None
    for index, line in enumerate(lines):
        record = line[:6].strip()
        if record in {"CONECT", "ANISOU", "MASTER"}:
            continue
        if index in blocks:
            residue = next(r for r in residues if r.line_indices[0] == index)
            for replacement in blocks[index]:
                serial += 1
                name = replacement[12:16].strip()
                if name in BACKBONE and name in residue.atoms:
                    remap[int(residue.atoms[name][6:11])] = serial
                output.append(replacement[:6] + f"{serial:5d}" + replacement[11:])
                last_atom = output[-1]
        if index in dropped:
            continue
        if record in {"ATOM", "HETATM", "TER"}:
            serial += 1
            if line[6:11].strip():
                remap[int(line[6:11])] = serial
            if record == "TER" and len(line) < 27:
                line = f"TER   {serial:5d}\n"
            else:
                line = line[:6] + f"{serial:5d}" + line[11:]
                if record == "TER" and last_atom is not None:
                    line = line[:17] + last_atom[17:27] + line[27:]
        if record in {"ATOM", "HETATM"}:
            last_atom = line
        elif record in {"MODEL", "ENDMDL"}:
            last_atom = None
        output.append(line)
    if serial > 99999:
        raise CrosslinkIntegrityError("Composed caps exceeds PDB atom serial capacity")
    for line in lines:
        if line.startswith("CONECT"):
            values = [int(line[i:i+5]) for i in range(6, len(line.rstrip()), 5) if line[i:i+5].strip()]
            if values and values[0] in remap:
                mapped = [remap[v] for v in values if v in remap]
                if len(mapped) > 1:
                    # Insert before END rather than after the PDB terminator.
                    pos = next((i for i, row in enumerate(output) if row[:6].strip() == "END"), len(output))
                    output.insert(pos, "CONECT" + "".join(f"{v:5d}" for v in mapped) + "\n")
        elif line.startswith("ANISOU") and int(line[6:11]) in remap:
            atom_serial = remap[int(line[6:11])]
            pos = next(i for i, row in enumerate(output)
                       if row[:6].strip() in {"ATOM", "HETATM"} and int(row[6:11]) == atom_serial)
            atom_line = output[pos]
            output.insert(pos + 1, line[:6] + f"{atom_serial:5d}" + atom_line[11:27] + line[27:])
    result = "".join(output)
    after = read_pdb_residues(result.splitlines(keepends=True), "shared_sites result")
    if [(r.address, r.segment) for r in residues] != [(r.address, r.segment) for r in after]:
        raise CrosslinkIntegrityError("Composition changed backbone residue order/TER boundaries")
    for before, new in zip(residues, after):
        for name in BACKBONE & before.atoms.keys():
            if name not in new.atoms or before.atoms[name][30:54] != new.atoms[name][30:54]:
                raise CrosslinkIntegrityError(f"Composition moved backbone atom {before.label}:{name}")
    return result


def compose_shared_sites(system, variants, ratios, root, config=None):
    """Compose validated caps, then atomically publish and install their graph.

    Input caps are never mutated. Failure before publication leaves System and
    any existing output untouched. Reusing an existing output is an error.
    """
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    destination = root / OUTPUT_TYPE
    if destination.exists():
        raise CrosslinkIntegrityError(f"Refusing to overwrite existing composed caps: {destination}")
    if list(variants) != list(ratios):
        raise CrosslinkIntegrityError("Variant order must match ratio_mix labels")
    reference = next(iter(variants))
    if set(map(int, system.get_models())) != set(variants[reference]):
        raise CrosslinkIntegrityError("Active system models differ from mixing caps")
    networks = {label: resolve_network(caps, config) for label, caps in variants.items()}
    common, alternatives, extras = match_loci(networks)
    chosen, sources, quotas = choose_alternatives(alternatives, ratios)
    sources = {mid: sources.get(mid, reference) for mid in variants[reference]}
    selected = [(chosen[i], locus.options[chosen[i]]) for i, locus in enumerate(alternatives)]
    for locus in common + extras:
        owners = {sources[p[0]] for p in locus.positions}
        label = next(iter(owners)) if len(owners) == 1 else reference
        if label not in locus.options:
            label = next(iter(locus.options))
        selected.append((label, locus.options[label]))

    documents, lookup = _read_caps(variants)
    footprints = {k[:3] for network in networks.values() for k in network.markers}
    marker_names = {k[3] for network in networks.values() for k in network.markers}
    native = lambda name: RESIDUE_TO_MUTATION.get(name, "LYS") if name in marker_names else name
    for (label, key), residue in lookup.items():
        if native(residue.name) != native(lookup[reference, key].name):
            raise CrosslinkIntegrityError(f"Different native residues in mix templates at {key}: {residue.name}, {lookup[reference, key].name}")
        if key not in footprints and residue.name != lookup[reference, key].name:
            raise CrosslinkIntegrityError(f"Different non-crosslink residue at {key}")

    replacements, alignments, expected = {}, [], set()
    tolerance = float(getattr(config, "mix_alignment_tolerance", 0.25))
    if not math.isfinite(tolerance) or not 0 < tolerance <= 1:
        raise CrosslinkIntegrityError("mix_alignment_tolerance must be >0 and <=1 A")
    for label, entity in selected:
        if expected & set(entity.markers):
            raise CrosslinkIntegrityError("Selected entities reuse a marker residue")
        expected.update(entity.markers)
        if all(sources[k[0]] == label for k in entity.markers):
            continue
        donor, recipient = [], []
        for key in entity.markers:
            for atom in ("N", "CA", "C"):
                try:
                    donor.append(_xyz(lookup[label, key[:3]].atoms[atom]))
                    recipient.append(_xyz(lookup[sources[key[0]], key[:3]].atoms[atom]))
                except KeyError as exc:
                    raise CrosslinkIntegrityError(f"Missing backbone atom {key}:{atom}") from exc
        rotation, translation, residual = _rigid_fit(donor, recipient, tolerance, entity.identity)
        alignments.append({"kind": entity.kind, "markers": entity.markers, "donor": label, "max_residual_A": residual})
        for key in entity.markers:
            mid = key[0]
            before, source = lookup[sources[mid], key[:3]], lookup[label, key[:3]]
            block = [line[:17] + source.name + line[20:] for name, line in before.atoms.items() if name in BACKBONE]
            for name, line in source.atoms.items():
                if name in BACKBONE:
                    continue
                xyz = np.asarray(_xyz(line)) @ rotation + translation
                coordinates = "".join(f"{v:8.3f}" for v in xyz)
                if len(coordinates) != 24 or not np.isfinite(xyz).all():
                    raise CrosslinkIntegrityError(f"Invalid transformed coordinates at {key}:{name}")
                anchor = next(iter(before.atoms.values()))
                block.append(line[:21] + anchor[21:27] + line[27:30] + coordinates + line[54:])
            replacements.setdefault(mid, {})[before.address] = block

    with tempfile.TemporaryDirectory(prefix=".shared-sites-", dir=root) as temporary:
        stage = Path(temporary) / OUTPUT_TYPE
        stage.mkdir()
        caps = {}
        for mid, label in sorted(sources.items()):
            path = stage / f"{mid}.caps.pdb"
            if mid in replacements:
                lines, residues = documents[label, mid]
                path.write_text(_render(lines, residues, replacements[mid]))
            else:
                shutil.copy2(variants[label][mid], path)
            caps[mid] = path
        markers = load_markers(caps)
        if set(markers) != expected:
            raise CrosslinkIntegrityError(f"Composed crosslink inventory differs: missing={expected-set(markers)}, extra={set(markers)-expected}")
        network = Network(markers, tuple(e for _, e in selected))
        # Re-resolving catches accidental alternative partners after grafting.
        resolved = resolve_network(caps, config)
        if {e.identity for e in resolved.entities} != {e.identity for e in network.entities}:
            raise CrosslinkIntegrityError("Composition changed crosslink partner identities")
        for label, entity in selected:
            for a, b in entity.bonds:
                if abs(distance(markers, a, b) - distance(networks[label].markers, a, b)) > 0.004:
                    raise CrosslinkIntegrityError(f"Composition changed forming-bond length: {a}, {b}")
        ring_report = validate_caps_rings(
            caps, network, repair=getattr(config, "mix_ring_repair", True),
            policy=getattr(config, "ring_penetration_policy", "warn"),
            report_path=root / "ring_geometry_report.json",
        )
        report = {
            "schema": 1, "strategy": "shared_sites", "ratio_mix": ratios,
            "families": quotas, "common_loci": len(common), "additional_loci": len(extras),
            "counts": network.counts, "model_sources": sources, "transfers": alignments,
            "selection": [{"source": label, "kind": e.kind, "markers": e.markers} for label, e in selected],
            "backbone": "N/CA/C/O and all other backbone coordinates unchanged from selected model variant",
            "ring_geometry": ring_report,
        }
        (stage / "crosslink_mix_report.json").write_text(json.dumps(report, indent=2) + "\n")
        stage.rename(destination)
    caps = {mid: destination / f"{mid}.caps.pdb" for mid in sources}
    for mid, label in sources.items():
        model = system.get_model(mid)
        model.mix_source = label
        model.type = OUTPUT_TYPE
        model.crosslink_type = OUTPUT_TYPE
    install_network(system, network, caps, root)
    shutil.copy2(destination / "crosslink_mix_report.json", root / "crosslink_mix_report.json")
    return report
