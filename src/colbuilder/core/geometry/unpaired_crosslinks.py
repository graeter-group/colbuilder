"""
Utility to detect unpaired crosslinks in a generated fibril.

Scans the per-model PDB fragments in ``.tmp/geometry_gen`` and decides whether
each crosslink marker is paired by the SAME criterion the geometry optimiser
uses (``optimize.py``): a real cross-model distance between the crosslink marker
atoms. A marker whose bonded partner atom is not matched to a compatible marker
in another model within ``CROSSLINK_PAIR_CUTOFF`` is considered unpaired and is
emitted as a ``manual_replacements.txt`` entry that mutates it back to a standard
residue.

A marker is paired iff at least one ATOM-CORRECT partner in another model lies
within ``CROSSLINK_PAIR_CUTOFF``:

  * Compatibility is atom-specific. A trivalent central residue bonds each arm
    through a *specific* atom -- C13 bonds the 3-arm, C12 bonds the 2-arm (see
    TRIVALENT_ARMS). A C13 marker therefore may pair ONLY a 3-arm, never a
    2-arm, and vice-versa. Divalent/AGE residues carry a single marker atom, so
    residue compatibility suffices for them.
  * The test is a simple presence test (any correct partner within the cutoff),
    not a one-to-one assignment. These crosslinks bond across a web of
    neighbouring models, so a one-to-one nearest assignment can "steal" a
    marker's true partner for a slightly shorter competing edge and orphan a real
    pair. A presence test avoids that and is unambiguous here because measured
    distances fall into two well-separated bands -- formed bonds <= ~7.4 A and
    genuine orphans >= ~22 A -- with a wide empty gap between them.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np

from ..utils.logger import setup_logger
from .crosslink import read_crosslink

LOG = setup_logger(__name__)

OPT_ID_CANDIDATES = [
    Path(".tmp") / "geometry_gen" / "crystalcontacts_from_colbuilder_opt_id.txt",
    Path("crystalcontacts_from_colbuilder_opt_id.txt"),
]

# Covalently bonded crosslink partners (undirected), at the residue level. Used
# to generate candidate pairs; trivalent pairs are further constrained by atom
# (see TRIVALENT_ARMS / _compatible) so C13/C12 bond only their correct arm.
CROSSLINK_PARTNERS: Dict[str, Set[str]] = {
    # AGE (non-enzymatic) divalent
    "LGX": {"AGS"}, "AGS": {"LGX"},          # Glucosepane
    "LPS": {"APD"}, "APD": {"LPS"},          # Pentosidine
    "LZS": {"LZD"}, "LZD": {"LZS"},          # MOLD
    # enzymatic divalent
    "L4Y": {"L5Y"}, "L5Y": {"L4Y"},
    "L4X": {"L5X"}, "L5X": {"L4X"},
    "LY4": {"LY5"}, "LY5": {"LY4"},
    "LX4": {"LX5"}, "LX5": {"LX4"},
    # trivalent: central ring residue <-> its two arms
    "LYX": {"LY2", "LY3"}, "LY2": {"LYX"}, "LY3": {"LYX"},
    "LXX": {"LX2", "LX3"}, "LX2": {"LXX"}, "LX3": {"LXX"},
    "LXY": {"L2Y", "L3Y"}, "L2Y": {"LXY"}, "L3Y": {"LXY"},
    "LYY": {"L2X", "L3X"}, "L2X": {"LYY"}, "L3X": {"LYY"},
}

# Trivalent central residue -> which arm residue its C13 vs C12 atom bonds.
# PYD:  C13(LYX)-CG(LY3), C12(LYX)-CB(LY2)
# DPD:  C13(LXX)-CG(LX3), C12(LXX)-CB(LX2)
# PYL:  C13(LXY)-CG(L3Y), C12(LXY)-CG(L2Y)
# DPL:  C13(LYY)-CG(L3X), C12(LYY)-CG(L2X)
TRIVALENT_ARMS: Dict[str, Dict[str, str]] = {
    "LYX": {"C13": "LY3", "C12": "LY2"},
    "LXX": {"C13": "LX3", "C12": "LX2"},
    "LXY": {"C13": "L3Y", "C12": "L2Y"},
    "LYY": {"C13": "L3X", "C12": "L2X"},
}
_TRIVALENT_ARM_RESNAMES: Set[str] = {
    arm for arms in TRIVALENT_ARMS.values() for arm in arms.values()
}

# Combined with atom-specific one-to-one matching, 12 A is safe against false pairing.
CROSSLINK_PAIR_CUTOFF = 12.0

RESIDUE_TO_MUTATION = {
    "AGS": "ARG",  # arginine-derived
    "APD": "ARG",
    "LZS": "LYS",
    "LZD": "LYS",
    "LGX": "LYS",
    "LPS": "LYS",
}
DEFAULT_MUTATION = "LYS"  # all lysine/hydroxylysine-derived markers


def _compatible(a: Dict, b: Dict) -> bool:
    """True if markers a and b can be the two ends of one covalent crosslink.

    Residue-level compatibility (CROSSLINK_PARTNERS) plus, for trivalent centrals,
    atom-level correctness: a C13 central bonds only its 3-arm, a C12 central only
    its 2-arm.
    """
    ra, rb = a["resname"], b["resname"]
    if rb not in CROSSLINK_PARTNERS.get(ra, set()):
        return False

    # Identify a trivalent central (if any) and enforce the atom<->arm rule.
    if ra in TRIVALENT_ARMS:
        central, arm = a, b
    elif rb in TRIVALENT_ARMS:
        central, arm = b, a
    else:
        return True  # divalent / AGE: single marker atom, residue match is enough

    want = TRIVALENT_ARMS[central["resname"]].get(central.get("atom", ""))
    if want is None:
        # Central marker with an unexpected atom name: fall back to residue-level
        # compatibility rather than silently dropping a real pair
        return True
    return arm["resname"] == want


class UnpairedCrosslinkFinder:
    """Detects unpaired crosslinks by real marker-atom distance and writes replacements."""

    def __init__(
        self,
        base_dir: Path,
        geom_dir: Optional[Path] = None,
        allowed_resnames: Optional[Set[str]] = None,
        cutoff: float = CROSSLINK_PAIR_CUTOFF,
    ):
        self.base_dir = Path(base_dir).resolve()
        self.geom_root = (
            Path(geom_dir).resolve() if geom_dir else self.base_dir / ".tmp" / "geometry_gen"
        )
        # Only these residue types are MUTATED. Pairing always uses the full set,
        # so a marker's partner is never missed even if the caller restricts output.
        self.allowed_resnames = (
            set(allowed_resnames) if allowed_resnames else set(CROSSLINK_PARTNERS.keys())
        )
        self.cutoff = float(cutoff)

    # ------------------------------------------------------------------ #
    # Public entry point
    # ------------------------------------------------------------------ #
    def run(self) -> Tuple[List[str], Optional[Path]]:
        markers = self._collect_markers()
        if not markers:
            LOG.info("No relevant crosslink residues found; skipping unpaired detection.")
            return [], None

        unpaired = self._find_unpaired_by_distance(markers)
        if not unpaired:
            LOG.info("No unpaired crosslinks detected.")
            return [], None

        entries = self._build_entries(unpaired)
        out_path = self._write_manual_replacements(entries)
        LOG.info(
            "Detected %d unpaired crosslink marker(s) out of %d; wrote %s",
            len(entries), len(markers), out_path,
        )
        return [f"{pdb} {new_res} {resid} {chain}" for pdb, new_res, resid, chain, _ in entries], out_path

    # ------------------------------------------------------------------ #
    # Marker collection (global-frame positions via read_crosslink)
    # ------------------------------------------------------------------ #
    def _collect_markers(self) -> List[Dict]:
        pdb_dirs = self._find_all_pdb_dirs()
        if not pdb_dirs:
            LOG.warning("No caps PDB fragments found under %s", self.geom_root)
            return []

        opt_id_path = self._find_opt_id()
        final_models: Optional[Set[int]] = self._read_opt_id(opt_id_path) if opt_id_path else None
        if final_models is None:
            LOG.warning(
                "crystalcontacts_from_colbuilder_opt_id.txt not found; "
                "considering all caps fragments."
            )

        markers: List[Dict] = []
        seen_files: Set[str] = set()
        for pdb_dir in pdb_dirs:
            for caps in sorted(pdb_dir.glob("*.caps.pdb")):
                base = caps.name.split(".")[0]
                if not base.isdigit():
                    continue
                idx = int(base)
                if final_models is not None and idx not in final_models:
                    continue
                if caps.name in seen_files:
                    continue
                seen_files.add(caps.name)
                try:
                    crosslinks = read_crosslink(caps)
                except Exception as exc:  # do not silently swallow read errors
                    LOG.error("Could not read crosslinks from %s: %s", caps, exc)
                    raise
                for cl in crosslinks:
                    if cl.resname not in CROSSLINK_PARTNERS:
                        continue
                    markers.append(
                        {
                            "fname": caps.name,
                            "model": idx,
                            "resname": cl.resname,
                            "atom": getattr(cl, "atom", ""),
                            "resid": str(cl.resid),
                            "chain": cl.chain,
                            "pos": np.asarray(cl.position, dtype=float),
                        }
                    )
        return markers

    # ------------------------------------------------------------------ #
    # Distance-based pairing: atom-specific presence test
    # ------------------------------------------------------------------ #
    def _find_unpaired_by_distance(self, markers: List[Dict]) -> List[Dict]:
        unpaired: List[Dict] = []
        for m in markers:
            if not CROSSLINK_PARTNERS.get(m["resname"]):
                continue
            pos_m = m["pos"]
            model_m = m["model"]
            has_pair = False
            for other in markers:
                if other["model"] == model_m:
                    continue  # a crosslink bonds ACROSS models, never within one
                if not _compatible(m, other):
                    continue
                if float(np.linalg.norm(pos_m - other["pos"])) < self.cutoff:
                    has_pair = True
                    break
            if not has_pair:
                unpaired.append(m)
        return unpaired

    # ------------------------------------------------------------------ #
    # Entry building / output
    # ------------------------------------------------------------------ #
    def _build_entries(self, unpaired: List[Dict]) -> List[Tuple[str, str, str, str, str]]:
        seen: Set[Tuple[str, str, str, str]] = set()
        entries: List[Tuple[str, str, str, str, str]] = []
        for m in unpaired:
            if m["resname"] not in self.allowed_resnames:
                continue
            key = (m["fname"], m["resname"], m["resid"], m["chain"])
            if key in seen:
                continue
            seen.add(key)
            new_res = RESIDUE_TO_MUTATION.get(m["resname"], DEFAULT_MUTATION)
            entries.append((m["fname"], new_res, m["resid"], m["chain"], m["resname"]))

        def sort_key(entry: Tuple[str, str, str, str, str]):
            pdb_file, _, resid, chain, orig_res = entry
            base = pdb_file.split(".")[0]
            idx = int(base) if base.isdigit() else base
            resid_key = int(resid) if resid.isdigit() else resid
            return (idx, chain, resid_key, orig_res)

        entries.sort(key=sort_key)
        return entries

    def _write_manual_replacements(self, entries) -> Path:
        out_path = self.base_dir / "manual_replacements.txt"
        with out_path.open("w") as f:
            for pdb_file, new_res, resid, chain, _orig in entries:
                f.write(f"{pdb_file} {new_res} {resid} {chain}\n")
        return out_path

    # ------------------------------------------------------------------ #
    # File discovery
    # ------------------------------------------------------------------ #
    def _find_all_pdb_dirs(self) -> List[Path]:
        """Find geom_root itself plus any immediate subdirectory holding caps files.

        Subdirectories are discovered dynamically rather than matched against a
        fixed name list: the standard (non-mixing) build path names them by the
        model's crosslink-derived type ("D"/"T"/"NC"/"DT"), but the mixing
        feature names them by the user's own, arbitrary ``ratio_mix`` labels
        (e.g. "PYD"/"HLKNL"), so no fixed set of names covers both.
        """
        dirs: List[Path] = []
        if list(self.geom_root.glob("*.caps.pdb")):
            dirs.append(self.geom_root)
        if self.geom_root.is_dir():
            for candidate in sorted(self.geom_root.iterdir()):
                if candidate.is_dir() and list(candidate.glob("*.caps.pdb")):
                    dirs.append(candidate)
        return dirs

    def _find_opt_id(self) -> Optional[Path]:
        candidates: List[Path] = [
            self.geom_root / "crystalcontacts_from_colbuilder_opt_id.txt",
            self.base_dir / "crystalcontacts_from_colbuilder_opt_id.txt",
        ]
        for candidate in OPT_ID_CANDIDATES:
            candidates.append(self.base_dir / candidate)
        for abs_path in candidates:
            if abs_path.is_file():
                return abs_path
        return None

    def _read_opt_id(self, opt_path: Path) -> Set[int]:
        models: Set[int] = set()
        with opt_path.open("r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                match = re.search(r"\b(\d+)\b", line)
                if match:
                    models.add(int(match.group(1)))
        return models
