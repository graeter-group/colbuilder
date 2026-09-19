"""
Repair backbone bonds that martinize2 omits at strained crosslink geometry.

Why this is needed
------------------
martinize2 infers backbone connectivity from *distance*: it writes a
BB(i)-BB(i+1) bond only when the all-atom peptide C-N distance is
below its detection threshold (~1.9 A with ``-bonds-fudge 1.4``). 
At a trivalent crosslink the three triple helices are pulled close together 
and the backbone is locally stretched, so the peptide bond feeding the crosslink 
residue's neighbour can measure > 1.9 A (we observe 2.2-2.4 A). 
martinize2 then silently drops that BB-BB bond and the coarse-grained chain fragments 
into disconnected pieces.

Fix
---
After the Martini topology for a model has been written, add back exactly the
consecutive-residue BB-BB backbone bonds that are missing, reusing the force
constant / length of a neighbouring backbone bond (no geometry is changed, and
no bond is invented that martinize2 would not itself have written for
unstretched geometry).

a bond is added only when ALL hold:
  * only consecutive residues in the *same* chain (resid i -> i+1) are considered;
  * a bond is added only when the BB-BB pair is absent from BOTH ``[ bonds ]``
    and ``[ constraints ]`` (helical backbone may be represented as a
    constraint, which must not be duplicated);
  * adjacency never spans a ``TER`` record (chain / connection boundary).
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

from colbuilder.core.utils.logger import setup_logger

LOG = setup_logger(__name__)

# Extended-collagen Martini backbone bond (function type, length nm, force const).
# Used only as a fallback when no neighbouring BB-BB bond is available to copy.
DEFAULT_BB_PARAMS: Tuple[str, ...] = ("1", "0.356", "18000")

# The atom name martinize2 gives the backbone bead.
_BACKBONE_BEAD = "BB"

# Maximum BB-BB separation (Angstrom) for which a missing backbone bond will be
# re-added. Real consecutive backbone beads (including crosslink-strained ones)
# sit around 3.2-4.2 A; two separated strands' termini are far beyond this. This
# is the geometric guard that makes the repair independent of TER / chain-id
# labelling. Kept well above the observed maximum (~4.2 A) so genuinely strained
# bonds are still captured, and far below any inter-strand gap.
MAX_BB_BOND_DIST: float = 6.0

# A BB bead: (atom_index, resname, chain, resid, x, y, z).
_Bead = Tuple[int, str, str, int, float, float, float]


def _read_backbone_beads(merge_pdb: Path) -> List[Optional[_Bead]]:
    """Read backbone (BB) beads from the merged CG PDB, in atom order.

    The merge PDB is written in the same atom order as the corresponding
    ``col_N.itp`` ``[ atoms ]`` section, so the 1-based running index of an
    ATOM/HETATM record equals its itp atom index. A ``None`` entry is inserted
    wherever a ``TER`` record breaks the chain.
    """
    beads: List[Optional[_Bead]] = []
    idx = 0
    with open(merge_pdb) as fh:
        for line in fh:
            tag = line[0:6]
            if tag in ("ATOM  ", "HETATM"):
                idx += 1
                if line[12:16].strip() == _BACKBONE_BEAD:
                    try:
                        beads.append(
                            (
                                idx,
                                line[17:20].strip(),
                                line[21:22],
                                int(line[22:26]),
                                float(line[30:38]),
                                float(line[38:46]),
                                float(line[46:54]),
                            )
                        )
                    except ValueError:
                        continue
            elif tag.startswith("TER"):
                if beads and beads[-1] is not None:
                    beads.append(None)
    return beads


def _segments(beads: List[Optional[_Bead]]) -> List[List[_Bead]]:
    """Split the bead stream into per-chain contiguous segments (split on None)."""
    segs: List[List[_Bead]] = []
    cur: List[_Bead] = []
    for x in beads:
        if x is None:
            if cur:
                segs.append(cur)
                cur = []
        else:
            cur.append(x)
    if cur:
        segs.append(cur)
    return segs


def _find_section(lines: List[str], header: str) -> Tuple[int, int]:
    """Return (start, end) line indices of a ``[ header ]`` section.

    ``start`` is the line *after* the header; ``end`` is the first subsequent
    section header (``[ ... ]``) or end-of-file. Returns (-1, -1) if absent.
    """
    start = -1
    for i, line in enumerate(lines):
        if line.strip() == header:
            start = i + 1
            break
    if start < 0:
        return -1, -1
    end = len(lines)
    for i in range(start, len(lines)):
        if lines[i].lstrip().startswith("["):
            end = i
            break
    return start, end


def _rows(lines: List[str], start: int, end: int) -> List[List[str]]:
    """Return tokenised data rows in a section (skipping comments/directives)."""
    out: List[List[str]] = []
    if start < 0:
        return out
    for line in lines[start:end]:
        s = line.strip()
        if not s or s.startswith((";", "#")):
            continue
        out.append(s.split())
    return out


def _is_int(tok: str) -> bool:
    return tok.lstrip("-").isdigit()


def repair_backbone_bonds(
    itp_path: Union[str, Path], merge_pdb: Union[str, Path]
) -> List[str]:
    """Restore backbone bonds/angles/dihedrals martinize2 omitted at crosslinks.

    Parameters
    ----------
    itp_path
        Path to the ``col_N.itp`` topology written for the model.
    merge_pdb
        Path to the matching ``N.merge.pdb`` (same atom order as the itp).

    Returns
    -------
    list of str
        Human-readable descriptions of the joints that were repaired (empty if
        the topology was already complete).
    """
    itp_path = Path(itp_path)
    merge_pdb = Path(merge_pdb)

    if not itp_path.exists() or not merge_pdb.exists():
        LOG.debug(
            f"Backbone repair skipped (missing file): itp={itp_path.exists()}, "
            f"merge={merge_pdb.exists()}"
        )
        return []

    beads = _read_backbone_beads(merge_pdb)
    if not beads:
        return []
    segs = _segments(beads)

    lines = itp_path.read_text().splitlines(keepends=True)
    b_start, b_end = _find_section(lines, "[ bonds ]")
    if b_start < 0:
        LOG.warning(f"No [ bonds ] section in {itp_path.name}; cannot repair backbone.")
        return []
    c_start, c_end = _find_section(lines, "[ constraints ]")
    a_start, a_end = _find_section(lines, "[ angles ]")
    d_start, d_end = _find_section(lines, "[ dihedrals ]")

    bond_rows = _rows(lines, b_start, b_end)
    cons_rows = _rows(lines, c_start, c_end)
    angle_rows = _rows(lines, a_start, a_end)
    dih_rows = _rows(lines, d_start, d_end)

    connected = set()
    bond_params: Dict[frozenset, Tuple[str, ...]] = {}
    for t in bond_rows:
        if len(t) >= 2 and _is_int(t[0]) and _is_int(t[1]):
            key = frozenset((int(t[0]), int(t[1])))
            connected.add(key)
            bond_params[key] = tuple(t[2:])
    for t in cons_rows:
        if len(t) >= 2 and _is_int(t[0]) and _is_int(t[1]):
            connected.add(frozenset((int(t[0]), int(t[1]))))

    # Position lookup: bead atom_index -> (segment idx, position in segment).
    pos: Dict[int, Tuple[int, int]] = {}
    for si, seg in enumerate(segs):
        for i, bead in enumerate(seg):
            pos[bead[0]] = (si, i)

    # Existing backbone angle params keyed by middle bead, and dihedral params
    # keyed by central bond -- but only for terms made of *consecutive* BB beads,
    # so we copy true backbone parameters (never a side-chain term).
    angle_tmpl: Dict[int, Tuple[str, ...]] = {}
    for t in angle_rows:
        if len(t) < 4 or not all(_is_int(t[i]) for i in range(3)):
            continue
        a, b, c = int(t[0]), int(t[1]), int(t[2])
        if a in pos and b in pos and c in pos:
            (sa, ia), (sb, ib), (sc, ic) = pos[a], pos[b], pos[c]
            if sa == sb == sc and ib == ia + 1 and ic == ib + 1:
                angle_tmpl.setdefault(b, tuple(t[3:]))
    dih_tmpl: Dict[Tuple[int, int], Tuple[str, ...]] = {}
    for t in dih_rows:
        if len(t) < 5 or not all(_is_int(t[i]) for i in range(4)):
            continue
        ids = [int(t[i]) for i in range(4)]
        if all(x in pos for x in ids):
            ps = [pos[x] for x in ids]
            segset = {s for s, _ in ps}
            order = [i for _, i in ps]
            if len(segset) == 1 and order == list(range(order[0], order[0] + 4)):
                # params start AFTER the 4 atom indices (t[4:]); slicing t[3:]
                # would leave the 4th atom index in the function-type column,
                # producing malformed dihedrals ("Invalid dihedral type <atom>").
                dih_tmpl.setdefault((ids[1], ids[2]), tuple(t[4:]))

    def nearest_angle_params(mid: int) -> Optional[Tuple[str, ...]]:
        if not angle_tmpl:
            return None
        return angle_tmpl[min(angle_tmpl, key=lambda m: abs(m - mid))]

    def nearest_dih_params(bc: Tuple[int, int]) -> Optional[Tuple[str, ...]]:
        if not dih_tmpl:
            return None
        mid = (bc[0] + bc[1]) / 2.0
        return dih_tmpl[min(dih_tmpl, key=lambda k: abs((k[0] + k[1]) / 2.0 - mid))]

    have_angle = set()
    for t in angle_rows:
        if len(t) >= 3 and all(_is_int(t[i]) for i in range(3)):
            tri = (int(t[0]), int(t[1]), int(t[2]))
            have_angle.add(tri)
            have_angle.add(tri[::-1])
    have_dih = set()
    for t in dih_rows:
        if len(t) >= 4 and all(_is_int(t[i]) for i in range(4)):
            q = (int(t[0]), int(t[1]), int(t[2]), int(t[3]))
            have_dih.add(q)
            have_dih.add(q[::-1])

    new_bonds: List[str] = []
    new_angles: List[str] = []
    new_dih: List[str] = []
    report: List[str] = []

    for seg in segs:
        idxs = [b[0] for b in seg]
        for k in range(len(seg) - 1):
            a_bead, b_bead = seg[k], seg[k + 1]
            if b_bead[3] != a_bead[3] + 1:
                continue
            u, v = a_bead[0], b_bead[0]
            if frozenset((u, v)) in connected:
                continue
            dist = math.dist(a_bead[4:7], b_bead[4:7])
            if dist > MAX_BB_BOND_DIST:
                LOG.debug(
                    f"{itp_path.name}: {a_bead[1]}{a_bead[3]}-{b_bead[1]}{b_bead[3]} "
                    f"chain {a_bead[2]} left unbonded (BB-BB {dist:.1f} A > "
                    f"{MAX_BB_BOND_DIST} A; treated as a chain break)."
                )
                continue

            # --- bond ---
            params = bond_params.get(frozenset((idxs[k - 1], u))) if k >= 1 else None
            if not params:
                params = bond_params.get(frozenset((v, idxs[k + 2]))) if k + 2 < len(seg) else None
            if not params or len(params) < 3:
                params = DEFAULT_BB_PARAMS
            new_bonds.append(f"{u:>5} {v:>5} {params[0]} {params[1]} {params[2]}\n")
            connected.add(frozenset((u, v)))

            # --- angles spanning the joint: (k-1,k,k+1) and (k,k+1,k+2) ---
            n_ang = 0
            for trip in ((k - 1, k, k + 1), (k, k + 1, k + 2)):
                if not all(0 <= j < len(seg) for j in trip):
                    continue
                x, y, z = idxs[trip[0]], idxs[trip[1]], idxs[trip[2]]
                if (x, y, z) in have_angle:
                    continue
                ap = nearest_angle_params(y)
                if ap is None:
                    continue
                new_angles.append(f"{x:>5} {y:>5} {z:>5} {' '.join(ap)}\n")
                have_angle.add((x, y, z))
                have_angle.add((z, y, x))
                n_ang += 1

            # --- dihedrals whose central bond is the joint: windows k-2..k ---
            n_dih = 0
            for st in (k - 2, k - 1, k):
                quad = (st, st + 1, st + 2, st + 3)
                if not all(0 <= j < len(seg) for j in quad):
                    continue
                w, x, y, z = (idxs[j] for j in quad)
                if (w, x, y, z) in have_dih:
                    continue
                dp = nearest_dih_params((x, y))
                if dp is None:
                    continue
                new_dih.append(f"{w:>5} {x:>5} {y:>5} {z:>5} {' '.join(dp)}\n")
                have_dih.add((w, x, y, z))
                have_dih.add((z, y, x, w))
                n_dih += 1

            report.append(
                f"{a_bead[1]}{a_bead[3]}-{b_bead[1]}{b_bead[3]} chain {a_bead[2]} "
                f"(BB {u}-{v}, {dist:.1f} A): +1 bond, +{n_ang} angles, +{n_dih} dihedrals"
            )

    if not report:
        return []

    # Insert each block at the end of its section. Do it bottom-up so earlier
    # insertions do not shift the indices of later ones.
    def _block(kind: str, rows: List[str]) -> List[str]:
        return [f"; backbone {kind} re-added by colbuilder (martinize2 omitted at crosslink strain)\n"] + rows

    inserts = []
    if new_dih and d_end >= 0:
        inserts.append((d_end, _block("dihedrals", new_dih)))
    if new_angles and a_end >= 0:
        inserts.append((a_end, _block("angles", new_angles)))
    if new_bonds:
        inserts.append((b_end, _block("bonds", new_bonds)))
    for at, block in sorted(inserts, key=lambda p: p[0], reverse=True):
        lines[at:at] = block

    itp_path.write_text("".join(lines))
    LOG.info(
        f"Repaired {len(report)} backbone joint(s) omitted by martinize2 in "
        f"{itp_path.name}: {'; '.join(report)}"
    )
    return report
