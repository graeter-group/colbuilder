"""
Utility to detect unpaired enzymatic crosslinks in a generated fibril.

This mirrors the standalone ``find_unpaired.py`` helper but packages the
logic so it can be invoked directly from the geometry pipeline. It scans
the PDB fragments in ``.tmp/geometry_gen``, reconstructs the pairing
information, and produces ``manual_replacements.txt`` entries for any
crosslinks that do not find a partner on the same connect row.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

from ..utils.logger import setup_logger

LOG = setup_logger(__name__)

# Directory names and filenames used in the geometry stage
PDB_SUBDIR_NAMES = ["D", "T", "NC", "DT"]
CONNECT_SUBDIRS = ["", "D", "T", "NC", "DT"]
OPT_ID_CANDIDATES = [
    Path(".tmp") / "geometry_gen" / "crystalcontacts_from_colbuilder_opt_id.txt",
    Path("crystalcontacts_from_colbuilder_opt_id.txt"),
]

# Crosslink groupings: pairs and triplets
PAIR_TRIPLE_GROUPS = [
    ("LGX", "AGS"),
    ("L4Y", "L5Y"),
    ("L4X", "L5X"),
    ("LY4", "LY5"),
    ("LX4", "LX5"),
    ("APD", "LPS"),
    ("LZD", "LZS"),
    ("LYX", "LY3", "LY2"),
    ("LXX", "LX2", "LX3"),
    ("LXY", "L2Y", "L3Y"),
    ("LYY", "L2X", "L3X"),
]

ALL_RESIDUE_TYPES = sorted({name for group in PAIR_TRIPLE_GROUPS for name in group})

RESIDUE_TO_MUTATION = {
    "AGS": "ARG",
    "APD": "ARG",
    "LZS": "LYS",
    "LZD": "LYS",
    "LGX": "LYS",
    "LPS": "LYS",
}
DEFAULT_MUTATION = "LYS"


class UnpairedCrosslinkFinder:
    """Detects unpaired enzymatic crosslinks and writes manual replacements."""

    def __init__(
        self,
        base_dir: Path,
        geom_dir: Optional[Path] = None,
        allowed_resnames: Optional[Set[str]] = None,
    ):
        self.base_dir = Path(base_dir).resolve()
        self.geom_root = (
            Path(geom_dir).resolve() if geom_dir else self.base_dir / ".tmp" / "geometry_gen"
        )
        self.allowed_resnames = set(allowed_resnames) if allowed_resnames else set(ALL_RESIDUE_TYPES)

    def run(self) -> Tuple[List[str], Optional[Path]]:
        """
        Run the full detection flow.

        Returns:
            (entries, path): entries are replacement lines
                "<pdb_file> <NEW_RES> <RESID> <CHAIN>".
        """
        crosslinks_by_type, types_by_model = self._build_crosslink_data()
        if not crosslinks_by_type:
            LOG.info("No relevant crosslink residues found; skipping unpaired detection.")
            return [], None

        connect_path = self._find_connect_file()
        if not connect_path:
            LOG.warning("connect_from_colbuilder.txt not found; skipping unpaired detection.")
            return [], None

        rows = self._parse_connect_rows(connect_path)
        paired_models = self._analyse_pairings(crosslinks_by_type, types_by_model, rows)
        entries = self._build_unpaired_entries(crosslinks_by_type, paired_models)

        if not entries:
            LOG.info("No unpaired crosslinks detected.")
            return [], None

        out_path = self._write_manual_replacements(entries)
        return [f"{pdb} {new_res} {resid} {chain}" for pdb, new_res, resid, chain, _ in entries], out_path

    # ------------------------------------------------------------------ #
    # File discovery helpers
    # ------------------------------------------------------------------ #
    def _find_opt_id(self) -> Optional[Path]:
        """Locate crystalcontacts_from_colbuilder_opt_id.txt."""
        candidates: List[Path] = [
            self.geom_root / "crystalcontacts_from_colbuilder_opt_id.txt",
            self.base_dir / "crystalcontacts_from_colbuilder_opt_id.txt",
        ]
        # fall back to legacy locations under base_dir
        for candidate in OPT_ID_CANDIDATES:
            candidates.append(self.base_dir / candidate)
        for abs_path in candidates:
            if abs_path.is_file():
                return abs_path
        return None

    def _read_opt_id(self, opt_path: Path) -> Set[int]:
        """Return model indices that belong to the final PDB."""
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

    def _find_pdb_dir(self) -> Optional[Path]:
        """Locate the folder containing *.caps.pdb fragments."""
        # If caps are directly under geom_root, use it
        if list(self.geom_root.glob("*.caps.pdb")):
            return self.geom_root

        found: List[Path] = []
        for sub in PDB_SUBDIR_NAMES:
            candidate = self.geom_root / sub
            if candidate.is_dir():
                found.append(candidate)

        if not found:
            return None

        if len(found) > 1:
            LOG.warning(
                "Multiple geometry subfolders found (%s); using %s",
                ", ".join(p.name for p in found),
                found[0],
            )
        return found[0]

    def _find_connect_file(self) -> Optional[Path]:
        """Locate connect_from_colbuilder.txt."""
        candidates: List[Path] = [
            self.geom_root / "connect_from_colbuilder.txt",
            self.base_dir / "connect_from_colbuilder.txt",
        ]
        for sub in CONNECT_SUBDIRS:
            candidates.append(self.geom_root / sub / "connect_from_colbuilder.txt")
            candidates.append(self.base_dir / sub / "connect_from_colbuilder.txt")

        for cand in candidates:
            if cand.is_file():
                return cand
        return None

    # ------------------------------------------------------------------ #
    # Parsing helpers
    # ------------------------------------------------------------------ #
    @staticmethod
    def _parse_connect_rows(connect_path: Path) -> List[List[str]]:
        rows: List[List[str]] = []
        with connect_path.open("r") as f:
            for line in f:
                tokens = line.strip().split()
                models = [t for t in tokens if t.endswith(".caps.pdb") or t.endswith(".pdb")]
                if models:
                    rows.append(models)
        return rows

    @staticmethod
    def _scan_pdb_for_residues(pdb_path: Path, residue_types: Iterable[str]) -> Dict[str, List[Tuple[str, str]]]:
        found = {res: set() for res in residue_types}
        with pdb_path.open("r") as f:
            for line in f:
                if not (line.startswith("ATOM") or line.startswith("HETATM")):
                    continue
                res_name = line[17:20].strip()
                if res_name not in residue_types:
                    continue
                chain_id = line[21].strip() or "A"
                resid = line[22:26].strip()
                found[res_name].add((chain_id, resid))

        def sort_key(item: Tuple[str, str]):
            chain, resid = item
            return (chain, int(resid)) if resid.isdigit() else (chain, resid)

        return {k: sorted(v, key=sort_key) for k, v in found.items()}

    def _build_crosslink_data(self) -> Tuple[Dict[str, Dict[str, str]], Dict[str, Set[str]]]:
        pdb_dir = self._find_pdb_dir()
        if not pdb_dir:
            LOG.warning("No geometry subfolders found under %s", self.geom_root)
            return {}, {}

        opt_id_path = self._find_opt_id()
        if not opt_id_path:
            LOG.warning("crystalcontacts_from_colbuilder_opt_id.txt not found; skipping unpaired detection.")
            return {}, {}

        final_models = self._read_opt_id(opt_id_path)
        if not final_models:
            LOG.warning("No valid model indices found in %s", opt_id_path)
            return {}, {}

        crosslinks_by_type: Dict[str, Dict[str, str]] = {res: {} for res in self.allowed_resnames}
        types_by_model: Dict[str, Set[str]] = {}

        for fname in sorted(p.name for p in pdb_dir.glob("*.caps.pdb")):
            base = fname.split(".")[0]
            if not base.isdigit():
                continue
            idx = int(base)
            if idx not in final_models:
                continue

            pdb_path = pdb_dir / fname
            found = self._scan_pdb_for_residues(pdb_path, self.allowed_resnames)
            for res_name, residues in found.items():
                if not residues:
                    continue
                if res_name not in crosslinks_by_type:
                    continue
                chain_id, resid = residues[-1]
                line = f"{fname} {res_name} {resid} {chain_id}"
                crosslinks_by_type[res_name][fname] = line
                types_by_model.setdefault(fname, set()).add(res_name)

        crosslinks_by_type = {r: m for r, m in crosslinks_by_type.items() if m}
        return crosslinks_by_type, types_by_model

    @staticmethod
    def _analyse_pairings(
        crosslinks_by_type: Dict[str, Dict[str, str]],
        types_by_model: Dict[str, Set[str]],
        rows: List[List[str]],
    ) -> Dict[str, Set[str]]:
        all_types = set(crosslinks_by_type.keys())
        paired_models: Dict[str, Set[str]] = {t: set() for t in all_types}

        for row in rows:
            if len(row) < 2:
                # If a single model carries all residues in a group, count them as paired
                row_types_by_model = {
                    m: (types_by_model.get(m, set()) & all_types) for m in row
                }
                for group in PAIR_TRIPLE_GROUPS:
                    group_set = set(group)
                    if group_set.issubset(all_types):
                        for model_name, res_set in row_types_by_model.items():
                            if group_set.issubset(res_set):
                                for res_type in group:
                                    paired_models[res_type].add(model_name)
                continue

            row_types_by_model = {
                m: (types_by_model.get(m, set()) & all_types) for m in row
            }

            for group in PAIR_TRIPLE_GROUPS:
                group_set = set(group)
                if not group_set.issubset(all_types):
                    continue

                models_by_type: Dict[str, List[str]] = {}
                for res_type in group:
                    models_with_type = [m for m in row if res_type in row_types_by_model.get(m, set())]
                    if not models_with_type:
                        models_by_type = {}
                        break
                    models_by_type[res_type] = models_with_type

                if not models_by_type:
                    continue

                all_models_in_group = {m for lst in models_by_type.values() for m in lst}
                if len(all_models_in_group) < 2:
                    continue

                for res_type, models_list in models_by_type.items():
                    for m in models_list:
                        paired_models[res_type].add(m)

        return paired_models

    # ------------------------------------------------------------------ #
    # Output helpers
    # ------------------------------------------------------------------ #
    def _build_unpaired_entries(
        self,
        crosslinks_by_type: Dict[str, Dict[str, str]],
        paired_models: Dict[str, Set[str]],
    ) -> List[Tuple[str, str, str, str, str]]:
        entries: List[Tuple[str, str, str, str, str]] = []

        for res_type, mapping in crosslinks_by_type.items():
            all_models = set(mapping.keys())
            unpaired = all_models - paired_models.get(res_type, set())
            for model in unpaired:
                parts = mapping[model].split()
                if len(parts) < 4:
                    LOG.warning("Malformed mapping for %s: %s", model, mapping[model])
                    continue
                pdb_file, orig_res, resid, chain = parts[:4]
                new_res = RESIDUE_TO_MUTATION.get(orig_res, DEFAULT_MUTATION)
                entries.append((pdb_file, new_res, resid, chain, orig_res))

        def sort_key(entry: Tuple[str, str, str, str, str]):
            pdb_file, _, resid, chain, orig_res = entry
            base = pdb_file.split(".")[0]
            try:
                idx = int(base)
            except ValueError:
                idx = base
            try:
                resid_key = int(resid)
            except ValueError:
                resid_key = resid
            return (idx, chain, resid_key, orig_res)

        entries.sort(key=sort_key)
        return entries

    def _write_manual_replacements(
        self, entries: List[Tuple[str, str, str, str, str]]
    ) -> Path:
        out_path = self.base_dir / "manual_replacements.txt"
        with out_path.open("w") as f:
            for pdb_file, new_res, resid, chain, _ in entries:
                f.write(f"{pdb_file} {new_res} {resid} {chain}\n")
        LOG.info("Found %d unpaired residues; wrote %s", len(entries), out_path)
        return out_path