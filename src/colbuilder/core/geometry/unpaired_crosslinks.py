"""
Utility to detect unpaired crosslinks in a generated fibril.

Scans the per-model PDB fragments in ``.tmp/geometry_gen`` and decides whether
each crosslink marker is paired by the SAME criterion the geometry optimiser
uses (``optimize.py``): a real cross-model distance between the crosslink marker
atoms. A marker whose bonded partner atom is not within ``CROSSLINK_PAIR_CUTOFF``
of a compatible marker in another model is considered unpaired and is emitted as
a ``manual_replacements.txt`` entry that mutates it back to a standard residue.

This replaces the previous connect-row-membership heuristic, which credited a
marker as "paired" merely because its model shared a connect row with an
unrelated partner-type model — silently leaving genuinely unpaired markers
(notably AGE crosslinks) un-mutated.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np

from ..utils.logger import setup_logger
from .crosslink import read_crosslink

LOG = setup_logger(__name__)

PDB_SUBDIR_NAMES = ["D", "T", "NC", "DT"]
OPT_ID_CANDIDATES = [
    Path(".tmp") / "geometry_gen" / "crystalcontacts_from_colbuilder_opt_id.txt",
    Path("crystalcontacts_from_colbuilder_opt_id.txt"),
]

# Covalently bonded crosslink partners (undirected). A marker is paired iff a
# marker of one of its partner types is within CROSSLINK_PAIR_CUTOFF in another
# model. Trivalent central residues bond BOTH arms; the arms bond the central.
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

# Angstrom. Matches the paper's 0.3 nm contact criterion and optimize.py's 3.0 A
# unpaired-crosslink test. Measured formed crosslink marker distances are
# 1.7-2.4 A (AGE and enzymatic alike), so 3.0 A cleanly separates paired from
# unpaired (genuinely unpaired markers sit >20 A away).
CROSSLINK_PAIR_CUTOFF = 3.0

RESIDUE_TO_MUTATION = {
    "AGS": "ARG",  # arginine-derived
    "APD": "ARG",
    "LZS": "LYS",
    "LZD": "LYS",
    "LGX": "LYS",
    "LPS": "LYS",
}
DEFAULT_MUTATION = "LYS"  # all lysine/hydroxylysine-derived markers


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
                            "resid": str(cl.resid),
                            "chain": cl.chain,
                            "pos": np.asarray(cl.position, dtype=float),
                        }
                    )
        return markers

    # ------------------------------------------------------------------ #
    # Distance-based pairing (mirrors optimize.py._get_unpaired_crosslinks)
    # ------------------------------------------------------------------ #
    def _find_unpaired_by_distance(self, markers: List[Dict]) -> List[Dict]:
        unpaired: List[Dict] = []
        for m in markers:
            partners = CROSSLINK_PARTNERS.get(m["resname"], set())
            has_pair = False
            for other in markers:
                if other["model"] == m["model"]:
                    continue  # a crosslink bonds ACROSS models, never within one
                if other["resname"] not in partners:
                    continue
                if float(np.linalg.norm(m["pos"] - other["pos"])) < self.cutoff:
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
        dirs: List[Path] = []
        if list(self.geom_root.glob("*.caps.pdb")):
            dirs.append(self.geom_root)
        for sub in PDB_SUBDIR_NAMES:
            candidate = self.geom_root / sub
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
