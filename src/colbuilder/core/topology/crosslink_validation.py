"""Fail-closed validation of the final crosslink graph and bonded terms."""

from collections import Counter, defaultdict
from pathlib import Path

from ..geometry.crosslink_network import MARKER_NAMES, TRIVALENT, DIVALENT, CrosslinkIntegrityError


def canon(path):
    path = tuple(path)
    return min(path, path[::-1])


def sections(path):
    section = None
    for number, line in enumerate(Path(path).read_text().splitlines(), 1):
        text = line.split(";", 1)[0].strip()
        if not text or text.startswith("#"):
            continue
        if text.startswith("["):
            section = text.strip("[] ").lower()
        else:
            yield section, text.split(), number


def validate_crosslink_terms(path, expected_bonds, rtp_path=None):
    atoms, residues = {}, []
    rows = defaultdict(list)
    for section, fields, line in sections(path):
        if section == "atoms":
            key, name, index = (fields[2], fields[3]), fields[4], int(fields[0])
            if not residues or residues[-1][0] != key or name in residues[-1][1]:
                residues.append((key, {}))
            residues[-1][1][name] = index
            atoms[index] = (len(residues) - 1, fields[3], name, fields[1], float(fields[6]))
        else:
            rows[section].append(fields)
    if not atoms:
        raise CrosslinkIntegrityError(f"No atoms in {path}")
    expected = {canon(pair) for pair in expected_bonds}
    bonds = Counter(canon(map(int, fields[:2])) for fields in rows["bonds"])
    adjacency = defaultdict(set)
    observed = set()
    for (a, b), count in bonds.items():
        if a not in atoms or b not in atoms:
            raise CrosslinkIntegrityError(f"Bond points outside atom inventory in {path}: {(a, b)}")
        adjacency[a].add(b)
        adjacency[b].add(a)
        aa, bb = atoms[a], atoms[b]
        marker_link = aa[0] != bb[0] and (aa[1] in MARKER_NAMES or bb[1] in MARKER_NAMES)
        if marker_link and {aa[2], bb[2]} != {"C", "N"}:
            observed.add((a, b))
        if (a, b) in expected and count != 1:
            raise CrosslinkIntegrityError(f"Duplicated crosslink bond in {path}: {(a, b)}")
    if observed != expected:
        raise CrosslinkIntegrityError(
            f"Crosslink bonds differ in {path}: missing={sorted(expected-observed)}, "
            f"unexpected={sorted(observed-expected)}"
        )

    ports = defaultdict(set)
    allowed_roles = set()
    for rule in TRIVALENT.values():
        ports[rule[0][0]].update(rule[0][1:])
        for i, (name, atom) in enumerate(rule[1:], 1):
            ports[name].add(atom)
            allowed_roles.add(frozenset(((rule[0][0], rule[0][i]), (name, atom))))
    for rule in DIVALENT.values():
        allowed_roles.add(frozenset(rule))
        for name, atom in rule:
            ports[name].add(atom)
    for a, b in expected:
        if frozenset((atoms[i][1], atoms[i][2]) for i in (a, b)) not in allowed_roles:
            raise CrosslinkIntegrityError(f"Chemically wrong marker bond in {path}: {(a, b)}")
    for (_, name), mapped in residues:
        if not ports[name] <= mapped.keys():
            raise CrosslinkIntegrityError(f"Missing forming atoms for {name} in {path}")
    for index, atom in atoms.items():
        if atom[2] in ports[atom[1]]:
            degree = sum(index in bond for bond in expected)
            if degree != 1:
                raise CrosslinkIntegrityError(
                    f"Unbonded or multiply bonded forming atom in {path}: "
                    f"{atom[1]}:{atom[2]} (index {index}, degree {degree})"
                )
    have_angles = {canon(map(int, fields[:3])) for fields in rows["angles"]}
    have_proper = {
        canon(map(int, fields[:4])) for fields in rows["dihedrals"] if int(fields[4]) not in (2, 4)
    }
    have_improper = {
        canon(map(int, fields[:4])) for fields in rows["dihedrals"] if int(fields[4]) in (2, 4)
    }
    have_pairs = Counter(canon(map(int, fields[:2])) for fields in rows["pairs"])
    angles, proper, pairs = set(), set(), set()
    # Independently enumerate simple paths from both sides of every bond end.
    starts = set()
    for a, b in expected:
        starts.update((a, b))
        starts.update(adjacency[a])
        starts.update(adjacency[b])
        for x in adjacency[a] | adjacency[b]:
            starts.update(adjacency[x])
    for start in starts:
        stack = [(start,)]
        while stack:
            walk = stack.pop()
            if len(walk) >= 3 and any(canon(edge) in expected for edge in zip(walk, walk[1:])):
                (angles if len(walk) == 3 else proper).add(canon(walk))
            if len(walk) < 4:
                stack.extend(
                    walk + (neighbor,) for neighbor in adjacency[walk[-1]] if neighbor not in walk
                )
    for walk in proper:
        a, b = walk[0], walk[-1]
        hydrogen_pair = all(atoms[i][2].lstrip("0123456789").startswith("H") for i in (a, b))
        if b not in adjacency[a] and not adjacency[a] & adjacency[b] and not hydrogen_pair:
            pairs.add(canon((a, b)))
    missing = {
        "angles": angles - have_angles,
        "proper_dihedrals": proper - have_proper,
        "pairs_1_4": pairs - set(have_pairs),
    }
    if any(missing.values()):
        details = {key: sorted(values)[:8] for key, values in missing.items() if values}
        raise CrosslinkIntegrityError(f"Missing crosslink terms in {path}: {details}")
    if any(have_pairs[pair] != 1 for pair in pairs):
        raise CrosslinkIntegrityError(f"Duplicated crosslink 1-4 pairs in {path}")
    if rows["exclusions"] and expected:
        raise CrosslinkIntegrityError(
            f"Unexpected explicit exclusions in crosslinked AMBER ITP {path}"
        )

    if rtp_path is not None:
        templates = defaultdict(lambda: defaultdict(list))
        residue, subsection = None, None
        known = {"atoms", "bonds", "angles", "dihedrals", "impropers", "exclusions"}
        for line in Path(rtp_path).read_text().splitlines():
            text = line.split(";", 1)[0].strip()
            if not text or text.startswith("#"):
                continue
            if text.startswith("["):
                key = text.strip("[] ")
                if key in known:
                    subsection = key
                else:
                    residue, subsection = key, None
            elif residue in MARKER_NAMES and subsection:
                templates[residue][subsection].append(text.split())
        for (_, name), mapped in residues:
            if name not in MARKER_NAMES:
                continue
            if name not in templates:
                raise CrosslinkIntegrityError(f"Missing RTP marker template: {name}")
            template = templates[name]
            for fields in template["atoms"]:
                atom = atoms.get(mapped.get(fields[0]))
                if atom is None or atom[3] != fields[1] or abs(atom[4] - float(fields[2])) > 2e-5:
                    raise CrosslinkIntegrityError(
                        f"RTP atom/type/charge mismatch: {name}:{fields[0]}"
                    )
            for section, size, present in [
                ("bonds", 2, set(bonds)),
                ("impropers", 4, have_improper),
            ]:
                for fields in template[section]:
                    names = fields[:size]
                    if any(n.startswith(("+", "-")) for n in names):
                        continue
                    if (
                        not all(n in mapped for n in names)
                        or canon(mapped[n] for n in names) not in present
                    ):
                        raise CrosslinkIntegrityError(f"Missing RTP {section}: {name}:{names}")
    return {
        "bonds": len(expected),
        "angles": len(angles),
        "proper_dihedrals": len(proper),
        "pairs_1_4": len(pairs),
    }
