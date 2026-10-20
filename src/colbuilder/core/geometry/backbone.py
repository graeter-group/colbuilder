"""Preserve polymer identity through sidechain mutation and validate AMBER output.

The input PDB's chain/TER boundaries, not a distance cutoff, define continuity.
Chimera contributes only the requested sidechains; it does not redefine polymers.
"""

from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Tuple
import math
import os
import tempfile

from ..utils.logger import setup_logger

LOG = setup_logger(__name__)
BACKBONE = frozenset(
    {
        "N",
        "CA",
        "C",
        "O",
        "OXT",
        "OC1",
        "OC2",
        "H",
        "HN",
        "H1",
        "H2",
        "H3",
        "HA",
        "HA1",
        "HA2",
        "HXT",
    }
)
SIDECHAIN = {
    "LYS": {"CB", "CG", "CD", "CE", "NZ"},
    "ARG": {"CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"},
}


class BackboneIntegrityError(ValueError):
    """A mutation or topology would change the input polymer connectivity."""


@dataclass
class Residue:
    address: Tuple[int, str, str, str]
    name: str
    segment: int
    atoms: Dict[str, str] = field(default_factory=dict)
    line_indices: List[int] = field(default_factory=list)

    @property
    def label(self):
        _, chain, number, insertion = self.address
        return f"{self.name}{number}{insertion}.{chain or '<blank>'}"


def read_pdb_residues(lines: List[str], source: str) -> List[Residue]:
    """Read ordered residue blocks without discarding TER or MODEL boundaries."""
    residues = []
    model, segment, boundary = 0, 0, True
    for index, line in enumerate(lines):
        record = line[:6].strip()
        if record in {"MODEL", "ENDMDL", "END", "TER"}:
            boundary = True
            if record == "MODEL":
                model += 1
            continue
        if record not in {"ATOM", "HETATM"}:
            continue
        if len(line) < 54:
            raise BackboneIntegrityError(f"{source}:{index + 1}: truncated atom record")
        address = (model, line[21], line[22:26].strip(), line[26].strip())
        name, atom = line[17:20].strip(), line[12:16].strip()
        if boundary or not residues or residues[-1].address != address:
            if boundary or not residues or residues[-1].address[:2] != address[:2]:
                segment += 1
            residues.append(Residue(address, name, segment))
        residue = residues[-1]
        if residue.name != name or atom in residue.atoms or line[16].strip():
            raise BackboneIntegrityError(
                f"{source}:{index + 1}: ambiguous residue/alternate location at {residue.label}"
            )
        residue.atoms[atom] = line
        residue.line_indices.append(index)
        boundary = False
    if not residues:
        raise BackboneIntegrityError(f"{source}: no PDB atoms")
    return residues


def _xyz(line):
    values = tuple(float(line[start : start + 8]) for start in (30, 38, 46))
    if not all(math.isfinite(value) for value in values):
        raise BackboneIntegrityError("Non-finite atom coordinates in mutation output")
    return values


def _atomic_write(path, content):
    mode = path.stat().st_mode
    with tempfile.NamedTemporaryFile(
        dir=path.parent, prefix=path.name + ".backbone-", delete=False
    ) as handle:
        temporary = Path(handle.name)
        try:
            handle.write(content)
            handle.flush()
            os.fchmod(handle.fileno(), mode)
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise
    try:
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _targets(residues, instructions, source):
    targets = {}
    addresses = [r.address for r in residues]
    if len(set(addresses)) != len(addresses) or len({a[0] for a in addresses}) != 1:
        raise BackboneIntegrityError(
            f"{source}: repeated residue identities or multiple models; mutation is ambiguous"
        )
    for target, number, chain in instructions:
        if target not in SIDECHAIN:
            raise BackboneIntegrityError(f"{source}: unsupported replacement target {target}")
        matches = [
            i
            for i, r in enumerate(residues)
            if r.address[2] + r.address[3] == number and r.address[1].lower() == chain.lower()
        ]
        exact = [i for i in matches if residues[i].address[1] == chain]
        matches = exact or matches
        if len(matches) != 1:
            raise BackboneIntegrityError(
                f"{source}: mutation target {number}.{chain} is missing or ambiguous"
            )
        i = matches[0]
        if i in targets and targets[i] != target:
            raise BackboneIntegrityError(
                f"{source}: conflicting replacements for {residues[i].label}"
            )
        if not {"N", "CA", "C", "O"}.issubset(residues[i].atoms):
            raise BackboneIntegrityError(
                f"{source}: incomplete input backbone at {residues[i].label}"
            )
        targets[i] = target
    return targets


def _restore_sidechains(original, residues, targets, mutated, source):
    after = read_pdb_residues(mutated.splitlines(keepends=True), source + " (Chimera)")
    if [r.address for r in residues] != [r.address for r in after]:
        raise BackboneIntegrityError(
            f"{source}: Chimera removed, reordered or renamed residue identities"
        )

    replacements = {}
    for i, (before, result) in enumerate(zip(residues, after)):
        expected_name = targets.get(i, before.name)
        if result.name != expected_name:
            raise BackboneIntegrityError(
                f"{source}: expected {expected_name} at {before.label}, found {result.name}"
            )
        if i not in targets:
            continue
        for atom in {"N", "CA", "C", "O"} & result.atoms.keys():
            if math.dist(_xyz(before.atoms[atom]), _xyz(result.atoms[atom])) > 0.01:
                raise BackboneIntegrityError(
                    f"{source}: Chimera moved backbone atom {before.label}:{atom}"
                )
        sidechain = {name: line for name, line in result.atoms.items() if name not in BACKBONE}
        heavy = {
            name
            for name, line in sidechain.items()
            if not (line[76:78].strip() == "H" or name.lstrip("0123456789").startswith("H"))
        }
        if heavy != SIDECHAIN[expected_name]:
            raise BackboneIntegrityError(
                f"{source}: incomplete/unexpected {expected_name} sidechain at "
                f"{before.label}: {sorted(heavy)}"
            )
        # Original backbone records carry the original coordinates and identity.
        rows = [
            (line[:17] + expected_name.rjust(3) + line[20:], int(line[6:11]))
            for name, line in before.atoms.items()
            if name in BACKBONE
        ]
        for line in sidechain.values():
            _xyz(line)
            rows.append(
                (line[:17] + expected_name.rjust(3) + before.atoms["CA"][20:27] + line[27:], None)
            )
        replacements[i] = rows

    by_line = {line: i for i, residue in enumerate(residues) for line in residue.line_indices}
    target_serials = {int(line[6:11]) for i in targets for line in residues[i].atoms.values()}
    rows, previous = [], None
    for index, line in enumerate(original):
        record = line[:6].strip()
        if index in by_line:
            i = by_line[index]
            previous = i
            if i in replacements:
                if index == residues[i].line_indices[0]:
                    rows.extend(replacements[i])
            else:
                rows.append((line, int(line[6:11])))
        elif record == "ANISOU" and int(line[6:11]) in target_serials:
            continue
        elif record == "MASTER":
            continue  # Its counts no longer describe the mutated atom inventory.
        else:
            if record == "TER" and previous in targets and len(line) >= 27:
                line = line[:17] + targets[previous].rjust(3) + line[20:]
            rows.append((line, None))

    serial, mapping, output = 0, {}, []
    for line, old_serial in rows:
        record = line[:6].strip()
        if record in {"ATOM", "HETATM", "TER"}:
            serial += 1
            if serial > 99999:
                raise BackboneIntegrityError(f"{source}: PDB serial capacity exceeded")
            if old_serial is not None:
                mapping[old_serial] = serial
            line = (
                line[:6] + f"{serial:5d}" + line[11:] if len(line) >= 11 else f"TER   {serial:5d}\n"
            )
        elif record == "ANISOU":
            line = line[:6] + f"{mapping[int(line[6:11])]:5d}" + line[11:]
        output.append(line)

    # Drop CONECT edges involving replaced sidechains; RTP defines the new ones.
    # Keep and renumber untouched edges and explicit backbone connectivity.
    final = []
    for line in output:
        if line.startswith("CONECT"):
            ids = [
                int(line[start : start + 5])
                for start in range(6, len(line.rstrip()), 5)
                if line[start : start + 5].strip()
            ]
            if not ids or ids[0] not in mapping:
                continue
            partners = [mapping[n] for n in ids[1:] if n in mapping]
            if not partners:
                continue
            line = "CONECT" + "".join(f"{n:5d}" for n in [mapping[ids[0]]] + partners) + "\n"
        final.append(line)
    restored = read_pdb_residues(final, source + " (restored)")
    if [(r.address, r.segment) for r in restored] != [(r.address, r.segment) for r in residues]:
        raise BackboneIntegrityError(f"{source}: failed to preserve polymer boundaries")
    return "".join(final).encode("utf-8")


@contextmanager
def preserve_backbone(type_dir: Path, instructions: List[str]):
    """Apply only requested sidechains, or roll back every affected PDB on failure.

    All original polymer boundaries and non-target residues are authoritative.
    New TER records from Chimera never enter the resulting PDB. Existing breaks
    are retained, including breaks within one chain ID and at numbering gaps.
    """
    grouped = {}
    root = Path(type_dir).resolve()
    for instruction in instructions:
        fields = str(instruction).strip().strip("\"'").split()
        if len(fields) != 4:
            raise BackboneIntegrityError(f"Invalid mutation instruction: {instruction!r}")
        filename, target, number, chain = fields
        path = (root / filename).resolve()
        if not path.is_relative_to(root):
            raise BackboneIntegrityError(f"Mutation file is outside the caps directory: {filename}")
        grouped.setdefault(path, []).append((target.upper(), number, chain))
    snapshots = {}
    for path, mutations in grouped.items():
        content = path.read_bytes()
        lines = content.decode("utf-8").splitlines(keepends=True)
        residues = read_pdb_residues(lines, str(path))
        targets = _targets(residues, mutations, str(path))
        snapshots[path] = (content, lines, residues, targets)
    try:
        yield
        repaired = {
            path: _restore_sidechains(lines, residues, targets, path.read_text(), str(path))
            for path, (_, lines, residues, targets) in snapshots.items()
        }
        for path, content in repaired.items():
            _atomic_write(path, content)
        LOG.info(
            "Preserved backbone atoms and polymer boundaries in %d mutated PDB file(s)",
            len(repaired),
        )
    except BaseException:
        for path, (content, _, _, _) in snapshots.items():
            # Chimera can remove a file before failing; recreate it for rollback.
            if not path.exists():
                path.touch()
            _atomic_write(path, content)
        raise


def validate_backbone_topology(pdb_file: Path, itp_file: Path) -> None:
    """Fail on missing backbone bonds or peptide bonds across true PDB boundaries."""
    pdb = read_pdb_residues(Path(pdb_file).read_text().splitlines(keepends=True), str(pdb_file))
    residues, bonds, atoms = [], set(), {}
    section = None
    with Path(itp_file).open() as handle:
        for line in handle:
            line = line.split(";", 1)[0].strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith("["):
                section = line.strip("[] ").lower()
                continue
            fields = line.split()
            if section == "atoms":
                atom_id, number, name, atom = int(fields[0]), fields[2], fields[3], fields[4]
                if not residues or residues[-1][0] != (number, name) or atom in residues[-1][1]:
                    residues.append(((number, name), {}))
                residues[-1][1][atom] = atom_id
                atoms[atom_id] = (len(residues) - 1, atom)
            elif section == "bonds" and int(fields[2]) in {1, 2, 3, 4, 5, 7, 8}:
                bonds.add(frozenset(map(int, fields[:2])))
    if len(residues) != len(pdb):
        raise BackboneIntegrityError(
            f"{itp_file}: PDB/ITP residue count differs ({len(pdb)} vs {len(residues)})"
        )
    expected_peptides = set()
    for i, source in enumerate(pdb):
        if residues[i][0][0] != source.address[2]:
            raise BackboneIntegrityError(f"{itp_file}: residue ordering differs at {source.label}")
        if (
            i
            and source.segment == pdb[i - 1].segment
            and "C" in pdb[i - 1].atoms
            and "N" in source.atoms
        ):
            left, right = residues[i - 1][1].get("C"), residues[i][1].get("N")
            pair = frozenset((left, right))
            if None in pair or pair not in bonds:
                raise BackboneIntegrityError(
                    f"{itp_file}: missing peptide bond {pdb[i - 1].label}:C--{source.label}:N"
                )
            expected_peptides.add(pair)
    for i, (source, ((number, name), mapped)) in enumerate(zip(pdb, residues)):
        for a, b in (("N", "CA"), ("CA", "C"), ("C", "O")):
            if {a, b}.issubset(source.atoms):
                terminal = i + 1 == len(pdb) or source.segment != pdb[i + 1].segment
                mapped_b = mapped.get(b)
                if b == "O" and mapped_b is None and terminal:
                    mapped_b = mapped.get("OC1", mapped.get("O1"))
                if (
                    a not in mapped
                    or mapped_b is None
                    or frozenset((mapped[a], mapped_b)) not in bonds
                ):
                    raise BackboneIntegrityError(
                        f"{itp_file}: missing backbone bond {source.label}:{a}--{b}"
                    )
    for pair in bonds:
        a, b = (atoms[index] for index in pair)
        if a[0] != b[0] and {a[1], b[1]} == {"C", "N"} and pair not in expected_peptides:
            raise BackboneIntegrityError(
                f"{itp_file}: unexpected peptide bond across input polymer boundaries: "
                f"{sorted(pair)}"
            )
