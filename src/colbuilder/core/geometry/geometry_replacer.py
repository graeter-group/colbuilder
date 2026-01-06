"""
Colbuilder Manual Crosslink Replacement Module (Simplified & Robust)

This module is strictly for MANUAL replacement of residues based on a provided list.
It includes robust logic to copy manual_replacements.txt from the launch directory 
and ensures consistent absolute paths. Crucially, it uses AGGRESSIVE syncing to 
ensure mutated files overwrite ALL copies in geometry_gen.
"""

import os
import traceback
import subprocess
import shutil
import random
import math
import time
from pathlib import Path
from typing import Optional, Any, Tuple, List, Dict, Set
from colorama import Fore, Style

# Import safe_load from yaml to parse raw config files
try:
    import yaml
except ImportError:
    yaml = None

from ..utils.exceptions import GeometryGenerationError
from ..utils.logger import setup_logger
from ..utils.config import ColbuilderConfig
from .system import System
from .connect import Connect
from .model import Model

LOG = setup_logger(__name__)

# Residue mapping for paired crosslink replacement
PAIRING_RULES = [
    (["LGX"], ["AGS"]),  # Glucosylpane
    (["LPS"], ["APD"]),  # Pentosidine
    (["LZS"], ["LZD"]),  # MOLD
    # Enzymatic divalent crosslinks
    (["L5Y"], ["L4Y"]),  # HLKNL
    (["L5X"], ["L4Y"]),  # LKNL
    (["LY5"], ["LY4"]),  # deH-HLNL
    (["LX5"], ["LX4"]),  # deH-LNL
]

REPLACEMENT_MAP: Dict[str, str] = {
    "AGS": "ARG",  # Glucosylpane
    "APD": "ARG",
    "LGX": "LYS",
    "LPS": "LYS",
    "LZS": "LYS",
    "LZD": "LYS",
    "LX3": "LYS",
    "LX2": "LYS",
    "LXX": "LYS",
    "L3Y": "LYS",
    "L2Y": "LYS",
    "LXY": "LYS",
    "L3X": "LYS",
    "L2X": "LYS",
    "LYY": "LYS",
    "L5Y": "LYS",
    "L4Y": "LYS",
    "L5X": "LYS",
    "LY5": "LYS",
    "LX5": "LYS",
    "LY4": "LYS",
    "LX4": "LYS",
    # Enzymatic C/N-term markers (PYD)
    "LYX": "LYS",
    "LY3": "LYS",
    "LY2": "LYS",
}
DEFAULT_REPLACEMENT = "LYS"
PAIR_DISTANCE_CUTOFF = 5.0

PAIRED_RESIDUES: Set[str] = {"AGS", "APD", "LGX", "LPS", "LZS", "LZD", "L5Y", "L4Y", "L5X", "LY5", "LX5", "LY4", "LX4"}
NON_ENZYMATIC_PAIRED: Set[str] = {"AGS", "APD", "LGX", "LPS", "LZS", "LZD"}
ENZYMATIC_SINGLETONS: Set[str] = {
    "LYX",
    "LY3",
    "LY2",
    "LX3",
    "LX2",
    "LXX",
    "L3Y",
    "L2Y",
    "LXY",
    "L3X",
    "L2X",
    "LYY",
    "L5Y",
    "L4Y",
    "L5X",
    "LY5",
    "LX5",
    "LY4",
    "LX4",
}
ENZYMATIC_TRIOS = [
    ("LYX", "LY2", "LY3"),  # PYD
    ("LX3", "LX2", "LXX"),  # DPD
    ("L3Y", "L2Y", "LXY"),  # PYL
    ("L3X", "L2X", "LYY"),  # DPL
]
PYD_THRESHOLD = 5.0  # distance cutoff to group enzymatic trios


class CrosslinkReplacer:
    """Handles ONLY manual replacement of residues in collagen systems."""

    def __init__(self):
        self.file_manager = None
        # Residues that are considered crosslink markers when parsing raw PDBs
        self._crosslink_resnames: Set[str] = set(REPLACEMENT_MAP.keys())

    async def replace(
        self, system: Optional[System], config: ColbuilderConfig, temp_dir: Path
    ) -> Tuple[Optional[System], Optional[Path]]:
        
        LOG.info(f"{Fore.CYAN}Starting MANUAL replacement mode.{Style.RESET_ALL}")
        
        if system is None:
            LOG.warning("System is None. Replacement not supported in this simplified mode.")
            return None, None

        await self.replace_in_system(system, config, temp_dir)
        
        # Ricalcolo path per output finale
        working_dir_root = Path(config.working_directory).resolve()
        # Fix path corruption check
        str_wd = str(working_dir_root)
        if ".tmp" in str_wd and "geometry_gen" in str_wd:
             working_dir_root = working_dir_root.parent.parent.parent 

        actual_temp_dir = working_dir_root / ".tmp" / "replace_manual"
        output_pdb = actual_temp_dir / f"{config.output or 'output'}.pdb"
        
        LOG.info(f"Writing final PDB structure to: {output_pdb}")
        
        system.write_pdb(
            pdb_out=output_pdb,
            fibril_length=config.fibril_length,
            temp_dir=actual_temp_dir,
        )

        return system, output_pdb

    async def replace_in_system(
        self, system: Any, config: ColbuilderConfig, temp_dir: Optional[Path] = None
    ) -> Any:
        try:
            # =================================================================================
            # PATH SETUP
            # =================================================================================
            working_dir_root = Path(config.working_directory).resolve()
            
            # Path corruption fix
            str_wd = str(working_dir_root)
            if ".tmp" in str_wd and "geometry_gen" in str_wd:
                LOG.debug("Detected geometry_gen in working directory; using parent as root.")
                parts = list(working_dir_root.parts)
                if ".tmp" in parts:
                    idx = parts.index(".tmp")
                    working_dir_root = Path(*parts[:idx])
                LOG.debug(f"Corrected root path: {working_dir_root}")

            if temp_dir is None:
                temp_dir = working_dir_root / ".tmp" / "replace_crosslinks"
            temp_dir.mkdir(parents=True, exist_ok=True)
            stale_manual = temp_dir / "manual_replacements.txt"
            if stale_manual.exists():
                try:
                    stale_manual.unlink()
                except Exception:
                    pass

            geometry_gen_dir = working_dir_root / ".tmp" / "geometry_gen"
            replace_manual_dir = working_dir_root / ".tmp" / "replace_manual"

            # =================================================================================
            # STEP 1: SOURCE SELECTION + FILE SYNC
            # =================================================================================

            source_dir = geometry_gen_dir
            if getattr(config, "auto_fix_unpaired", False) and getattr(config, "manual_replacements", None):
                try:
                    await self._apply_manual_replacements_to_dir(
                        source_dir=geometry_gen_dir,
                        dest_dir=replace_manual_dir,
                        manual_list=[
                            str(instr).strip()
                            for instr in config.manual_replacements
                            if str(instr).strip()
                        ],
                        working_dir_root=working_dir_root,
                        config=config,
                        system=system,
                    )
                    if getattr(config, "_replacement_verbose", True):
                        LOG.section("Running crosslinks replacement")
                    source_dir = replace_manual_dir
                except Exception as e:
                    LOG.warning("Auto-fix unpaired replacement failed: %s", e)
                    source_dir = geometry_gen_dir
            elif not geometry_gen_dir.exists() and replace_manual_dir.exists():
                source_dir = replace_manual_dir

            if not source_dir.exists():
                LOG.error(
                    f"{Fore.RED}ERROR: Source directory not found at {source_dir}.{Style.RESET_ALL}"
                )
                source_dir = Path.cwd()

            type_dir_candidates = [source_dir]
            for sub in ["DT", "TD", "D", "T", "NC", "DY"]:
                cand = source_dir / sub
                if cand.exists():
                    type_dir_candidates.append(cand)

            model_zero = system.get_model(model_id=0.0)
            system_type = model_zero.type if hasattr(model_zero, "type") else "D"
            type_dir = temp_dir / system_type
            type_dir.mkdir(parents=True, exist_ok=True)

            connect_file = None
            for candidate in [
                source_dir / "connect_from_colbuilder.txt",
                source_dir / "connect_from_colbuilder",
            ]:
                if candidate.exists():
                    connect_file = candidate
                    break
            if connect_file and connect_file.exists():
                if connect_file.suffix != ".txt":
                    connect_file_txt = connect_file.with_suffix(".txt")
                    if connect_file_txt.exists():
                        connect_file = connect_file_txt
                shutil.copy2(connect_file, temp_dir / "connect_from_colbuilder.txt")
                connect_file = temp_dir / "connect_from_colbuilder.txt"
            else:
                try:
                    connector = Connect(system=system)
                    connect_file = temp_dir / "connect_from_colbuilder"
                    connector.write_connect(system=system, connect_file=connect_file)
                    connect_file = connect_file.with_suffix(".txt")
                    LOG.info("Rebuilt connectivity file at %s", connect_file)
                except Exception as e:
                    LOG.warning("Could not rebuild connectivity file: %s", e)
                    connect_file = None

            caps_by_model: Dict[int, Path] = {}
            for search_dir in type_dir_candidates:
                for pdb_file in search_dir.glob("*.caps.pdb"):
                    try:
                        model_id = int(pdb_file.stem.split(".")[0])
                    except ValueError:
                        continue
                    if model_id not in caps_by_model:
                        caps_by_model[model_id] = pdb_file

            connect_groups = (
                self._load_connect_groups(connect_file) if connect_file else []
            )

            for model_id, source_caps in sorted(caps_by_model.items()):
                dest_caps = type_dir / source_caps.name
                if source_caps.resolve() != dest_caps.resolve():
                    shutil.copy2(source_caps, dest_caps)

            # =================================================================================
            # STEP 3: RECUPERO LISTA E PREPARAZIONE
            # =================================================================================
            manual_list: List[str] = []
            ratio_requested = (
                config.ratio_replace is not None
                and float(config.ratio_replace) > 0
            )
            generated_from_ratio = False

            if not ratio_requested and getattr(config, "manual_replacements", None):
                manual_list = [
                    str(instr).strip()
                    for instr in config.manual_replacements
                    if str(instr).strip()
                ]

            if ratio_requested and connect_groups:
                records = self._load_crosslinks_from_models(type_dir)
                manual_list = self._build_ratio_replacements_from_connect(
                    records=records,
                    connect_groups=connect_groups,
                    ratio_replace=float(config.ratio_replace),
                    scope=getattr(config, "ratio_replace_scope", "non_enzymatic"),
                    config=config,
                )
                generated_from_ratio = bool(manual_list)

            if not manual_list:
                # Quietly skip when nothing to do
                if ratio_requested:
                    LOG.warning(
                        "No replacement instructions generated (ratio=%s%%, scope=%s). "
                        "No matching markers found in connected caps.",
                        config.ratio_replace,
                        getattr(config, "ratio_replace_scope", "non_enzymatic"),
                    )
                try:
                    config._replacement_skipped = True  # type: ignore[attr-defined]
                except Exception:
                    pass
                return system

            if generated_from_ratio:
                try:
                    target_file_path = temp_dir / "manual_replacements.txt"
                    with open(target_file_path, "w") as f:
                        for instruction in manual_list:
                            f.write(f"{instruction}\n")
                except Exception as e:
                    LOG.warning("Could not persist ratio-based instructions: %s", e)

            replace_file = temp_dir / "replace.txt"
            with open(replace_file, "w") as f:
                for instruction in manual_list:
                    clean_instr = instruction.strip().strip('"').strip("'")
                    f.write(f"{clean_instr}\n")

            # If the expected type dir is empty, fall back to any caps dir we detected
            if not list(type_dir.glob("*.caps.pdb")):
                fallback_dirs = []
                for t in ["DT", "TD", "D", "T", "NC"]:
                    cand = temp_dir / t
                    if cand not in fallback_dirs:
                        fallback_dirs.append(cand)
                if list(temp_dir.glob("*.caps.pdb")):
                    fallback_dirs.append(temp_dir)

                for cand in fallback_dirs:
                    try:
                        if cand.exists() and list(Path(cand).glob("*.caps.pdb")):
                            type_dir = Path(cand)
                            break
                    except Exception:
                        continue

            # Ensure the working type_dir contains caps (copy from root if needed)
            caps_in_root = list(temp_dir.glob("*.caps.pdb"))
            for f in caps_in_root:
                dest = type_dir / f.name
                if not dest.exists():
                    shutil.copy2(f, dest)
            
            # =================================================================================
            # STEP 4: ESECUZIONE CHIMERA
            # =================================================================================
            success = await self._run_chimera_command(
                config, 
                str(replace_file), 
                type_dir, 
                working_dir_root
            )
            
            if not success:
                LOG.error("Chimera execution failed.")
                return system

            LOG.debug("Chimera replacements executed.")

            # =================================================================================
            # STEP 5: AGGRESSIVE BACK-PROPAGATION (FIX TOPOLOGIA)
            # Sovrascriviamo QUALSIASI copia del file in geometry_gen
            # =================================================================================
            source_dir = type_dir # Dove Chimera ha salvato i file mutati (es. replace_manual/D)
            updated_count = 0
            
            if source_dir.exists():
                mutated_files = list(source_dir.glob("*.caps.pdb"))
                
                for mutated_file in mutated_files:
                    filename = mutated_file.name
                    
                    # 1. Sovrascrivi nella root di geometry_gen
                    root_target = geometry_gen_dir / filename
                    if root_target.exists():
                        shutil.copy2(mutated_file, root_target)
                    
                    # 2. Sovrascrivi in TUTTE le sottocartelle di geometry_gen (es. D/, NC/, T/)
                    # main_topology.py copia ricorsivamente, quindi dobbiamo pulire tutto.
                    for subdir in geometry_gen_dir.iterdir():
                        if subdir.is_dir():
                            sub_target = subdir / filename
                            if sub_target.exists():
                                shutil.copy2(mutated_file, sub_target)
                                # LOG.debug(f"Overwrote {filename} in {subdir.name}")
                    
                    updated_count += 1

            # Sync-back details suppressed for cleaner logs
            # =================================================================================

            # =================================================================================
            # STEP 6: RECOMPUTE CONNECTIVITY AFTER MANUAL CHANGES
            # Updates in-memory connects and rewrites connect_from_colbuilder.txt
            # =================================================================================
            try:
                connector = Connect(system=system)
                new_connect = connector.run_connect(system=system)

                if new_connect:
                    for mid, conns in new_connect.items():
                        model_obj = system.get_model(model_id=mid)
                        if model_obj:
                            model_obj.connect = conns

                    connect_file_path = temp_dir / "connect_from_colbuilder.txt"
                    connector.write_connect(system=system, connect_file=connect_file_path)

                    try:
                        geom_connect_path = geometry_gen_dir / "connect_from_colbuilder.txt"
                        shutil.copy2(connect_file_path, geom_connect_path)
                    except Exception as e:
                        LOG.warning(f"Could not sync connectivity to geometry_gen: {e}")
                else:
                    LOG.warning("Connectivity recomputation returned empty; keeping previous connects.")
            except Exception as e:
                LOG.warning(f"Failed to recompute connectivity after manual replacement: {e}")
            # =================================================================================

            return system

        except Exception as e:
            LOG.error(f"Error in manual replacement: {str(e)}")
            traceback.print_exc()
            return system

    def _build_ratio_replacements(
        self, system: Any, ratio_replace: float, fibril_length: float, scope: str = "non_enzymatic"
    ) -> List[str]:
        """
        Select crosslinks to mutate based on a density (ratio_replace) value.

        Returns a list of swap instructions in the format expected by swapaa.py,
        mapping markers to standard residues according to REPLACEMENT_MAP and scope.
        """
        if system is None:
            LOG.warning("Cannot build ratio-based replacements without a system object.")
            return []

        if ratio_replace <= 0:
            return []

        try:
            scope = scope or "non_enzymatic"
            scope = str(scope).strip().lower()
            if scope not in {"enzymatic", "non_enzymatic", "all"}:
                LOG.warning("Unknown ratio_replace scope '%s'; defaulting to 'non_enzymatic'", scope)
                scope = "non_enzymatic"

            z_bounds = self._calculate_fibril_bounds(system, fibril_length)
        except Exception as e:
            LOG.warning("Failed to calculate fibril bounds, using all crosslinks: %s", e)
            z_bounds = (-math.inf, math.inf)

        if scope == "enzymatic":
            target_resnames = ENZYMATIC_SINGLETONS.copy()
        elif scope == "non_enzymatic":
            target_resnames = NON_ENZYMATIC_PAIRED.copy()
        else:  # "all"
            target_resnames = NON_ENZYMATIC_PAIRED.union(ENZYMATIC_SINGLETONS)

        crosslinks_by_type: Dict[str, List[Dict[str, Any]]] = {res: [] for res in target_resnames}

        for model_id in system.get_models():
            model = system.get_model(model_id=model_id)
            if not hasattr(model, "crosslink") or not model.crosslink:
                continue

            for crosslink in model.crosslink:
                resname = getattr(crosslink, "resname", None)
                if resname not in target_resnames:
                    continue
                if not hasattr(crosslink, "position"):
                    continue

                try:
                    z_pos = self._get_component(crosslink.position, 2)
                except Exception as e:
                    LOG.debug("Skipping crosslink without valid z-position: %s", e)
                    continue

                if not (z_bounds[0] <= z_pos <= z_bounds[1]):
                    continue

                crosslinks_by_type.setdefault(resname, []).append(
                    {"model_id": model_id, "crosslink": crosslink, "model": model}
                )

        pairs: List[List[Dict[str, Any]]] = []
        singles: List[Dict[str, Any]] = []
        pyd_trios: List[List[Dict[str, Any]]] = []

        for donors, acceptors in PAIRING_RULES:
            # Skip pairs outside scope
            if not (target_resnames.issuperset(set(donors)) and target_resnames.issuperset(set(acceptors))):
                continue

            donor_list: List[Dict[str, Any]] = []
            acceptor_list: List[Dict[str, Any]] = []

            for resname in donors:
                donor_list.extend(crosslinks_by_type.get(resname, []))
            for resname in acceptors:
                acceptor_list.extend(crosslinks_by_type.get(resname, []))

            if not donor_list or not acceptor_list:
                continue

            used_donors = set()
            used_acceptors = set()

            for i, donor_ref in enumerate(donor_list):
                if i in used_donors:
                    continue

                donor = donor_ref["crosslink"]
                best_idx = -1
                best_dist = float("inf")

                for j, acc_ref in enumerate(acceptor_list):
                    if j in used_acceptors:
                        continue

                    acceptor = acc_ref["crosslink"]

                    try:
                        dist = self._calculate_distance(donor.position, acceptor.position)
                    except Exception as e:
                        LOG.debug("Distance calculation failed, skipping pair: %s", e)
                        continue

                    if dist <= PAIR_DISTANCE_CUTOFF and dist < best_dist:
                        best_idx = j
                        best_dist = dist

                if best_idx != -1:
                    pairs.append([donor_ref, acceptor_list[best_idx]])
                    used_donors.add(i)
                    used_acceptors.add(best_idx)

        # Group enzymatic markers into PYD trios (LYX + LY2 + LY3) when included in scope
        if scope in {"enzymatic", "all"}:
            lyx_list = crosslinks_by_type.get("LYX", [])
            ly2_list = crosslinks_by_type.get("LY2", [])
            ly3_list = crosslinks_by_type.get("LY3", [])
            pyd_trios = self._build_pyd_trios(lyx_list, ly2_list, ly3_list)

        if not pairs:
            pass
        num_pair_to_replace = (
            min(len(pairs), max(1, math.ceil(len(pairs) * ratio_replace / 100.0)))
            if pairs
            else 0
        )

        num_trio_to_replace = (
            min(len(pyd_trios), max(1, math.ceil(len(pyd_trios) * ratio_replace / 100.0)))
            if pyd_trios
            else 0
        )

        # Singleton enzymatic markers (e.g., LYX/LY3/LY2) included in scope
        if scope in {"enzymatic", "all"}:
            for resname in target_resnames:
                if resname in PAIRED_RESIDUES:
                    continue
                singles.extend(crosslinks_by_type.get(resname, []))
            # Remove any singleton that is already part of a PYD trio
            trio_members = {
                (getattr(x["crosslink"], "resid", ""), getattr(x["crosslink"], "chain", ""), x.get("model_id", 0.0))
                for trio in pyd_trios
                for x in trio
            }
            singles = [
                s
                for s in singles
                if (getattr(s["crosslink"], "resid", ""), getattr(s["crosslink"], "chain", ""), s.get("model_id", 0.0))
                not in trio_members
            ]

        num_single_to_replace = (
            min(len(singles), max(1, math.ceil(len(singles) * ratio_replace / 100.0)))
            if singles
            else 0
        )

        random.seed(int(time.time()))
        random.shuffle(pairs)
        selected_pairs = pairs[:num_pair_to_replace] if num_pair_to_replace else []
        random.shuffle(pyd_trios)
        selected_trios = pyd_trios[:num_trio_to_replace] if num_trio_to_replace else []
        random.shuffle(singles)
        selected_singles = singles[:num_single_to_replace] if num_single_to_replace else []

        instructions: List[str] = []
        seen = set()

        for pair in selected_pairs:
            for cross_ref in pair:
                crosslink = cross_ref["crosslink"]
                resname = getattr(crosslink, "resname", "")
                resid = getattr(crosslink, "resid", "").strip()
                chain = getattr(crosslink, "chain", "").strip() or "A"
                model_id = cross_ref.get("model_id", 0.0)

                if not resid:
                    continue

                new_res = REPLACEMENT_MAP.get(resname, DEFAULT_REPLACEMENT)
                model_file = f"{int(float(model_id))}.caps.pdb"
                instruction = f"{model_file} {new_res} {resid} {chain}"

                if instruction not in seen:
                    seen.add(instruction)
                    instructions.append(instruction)

        for trio in selected_trios:
            for cross_ref in trio:
                crosslink = cross_ref["crosslink"]
                resname = getattr(crosslink, "resname", "")
                resid = getattr(crosslink, "resid", "").strip()
                chain = getattr(crosslink, "chain", "").strip() or "A"
                model_id = cross_ref.get("model_id", 0.0)

                if not resid:
                    continue

                new_res = REPLACEMENT_MAP.get(resname, DEFAULT_REPLACEMENT)
                model_file = f"{int(float(model_id))}.caps.pdb"
                instruction = f"{model_file} {new_res} {resid} {chain}"

                if instruction not in seen:
                    seen.add(instruction)
                    instructions.append(instruction)

        for single in selected_singles:
            crosslink = single["crosslink"]
            resname = getattr(crosslink, "resname", "")
            resid = getattr(crosslink, "resid", "").strip()
            chain = getattr(crosslink, "chain", "").strip() or "A"
            model_id = single.get("model_id", 0.0)

            if not resid:
                continue

            new_res = REPLACEMENT_MAP.get(resname, DEFAULT_REPLACEMENT)
            model_file = f"{int(float(model_id))}.caps.pdb"
            instruction = f"{model_file} {new_res} {resid} {chain}"

            if instruction not in seen:
                seen.add(instruction)
                instructions.append(instruction)

        def sort_key(line: str) -> Tuple:
            parts = line.split()
            model_part = parts[0] if parts else ""
            try:
                model_idx = int(model_part.split(".")[0])
            except (ValueError, IndexError):
                model_idx = model_part

            resid_part = parts[2] if len(parts) > 2 else ""
            try:
                resid_idx = int(resid_part)
            except ValueError:
                resid_idx = resid_part

            chain_part = parts[3] if len(parts) > 3 else ""
            return (model_idx, chain_part, resid_idx, parts[1] if len(parts) > 1 else "")

        instructions.sort(key=sort_key)
        total_candidates = len(pairs) * 2 + len(singles) + len(pyd_trios) * 3
        selected_pairs_res = len(selected_pairs) * 2
        selected_singles_res = len(selected_singles)
        selected_trio_res = len(selected_trios) * 3
        # Keep silent when no replacements are selected to avoid noisy output
        return instructions

    def _calculate_fibril_bounds(
        self, system: Any, fibril_length: float
    ) -> Tuple[float, float]:
        try:
            z_values = []

            for model_id in system.get_models():
                model = system.get_model(model_id=model_id)
                try:
                    cog = model.get_cog()
                    z_values.append(self._get_z_position(cog))
                except Exception as e:
                    LOG.debug("Could not read COG for model %s: %s", model_id, e)

            if not z_values:
                return (-math.inf, math.inf)

            min_z = min(z_values) - 500.0
            max_z = max(z_values) + 500.0

            if max_z - min_z > 5000.0:
                model = system.get_model(model_id=0.0)
                cog = model.get_cog()
                z_center = self._get_z_position(cog)
                fibril_length_angstroms = fibril_length * 10.0 if fibril_length else 0.0
                return (
                    z_center - fibril_length_angstroms * 5,
                    z_center + fibril_length_angstroms * 5,
                )

            return (min_z, max_z)
        except Exception as e:
            LOG.warning("Failed to calculate fibril bounds: %s", e)
            return (-math.inf, math.inf)

    def _get_z_position(self, position: Any) -> float:
        if hasattr(position, "shape"):
            if len(position.shape) == 1 and position.shape[0] >= 3:
                return float(position[2])
            if len(position.shape) == 2 and position.shape[1] >= 3:
                return float(position[0][2])

        if hasattr(position, "__getitem__"):
            z_val = position[2]
            if hasattr(z_val, "__getitem__") and hasattr(z_val, "__len__") and len(z_val) > 0:
                return float(z_val[0])
            return float(z_val)

        return float(position)

    def _calculate_distance(self, pos1: Any, pos2: Any) -> float:
        x1 = self._get_component(pos1, 0)
        y1 = self._get_component(pos1, 1)
        z1 = self._get_component(pos1, 2)

        x2 = self._get_component(pos2, 0)
        y2 = self._get_component(pos2, 1)
        z2 = self._get_component(pos2, 2)

        dx = x2 - x1
        dy = y2 - y1
        dz = z2 - z1

        return math.sqrt(dx * dx + dy * dy + dz * dz)

    def _get_component(self, position: Any, index: int) -> float:
        if hasattr(position, "shape"):
            if len(position.shape) == 1 and position.shape[0] > index:
                return float(position[index])
            if len(position.shape) == 2 and position.shape[1] > index:
                return float(position[0][index])

        if hasattr(position, "__getitem__"):
            val = position[index]
            if hasattr(val, "__getitem__") and hasattr(val, "__len__") and len(val) > 0:
                return float(val[0])
            return float(val)

        return float(position)

    async def _run_chimera_command(
        self,
        config: ColbuilderConfig,
        replace_file_path: str,
        type_dir_path: Path,
        root_dir: Path
    ) -> bool:
        try:
            swapaa_script = None
            search_paths = [
                Path(config.CHIMERA_SCRIPTS_DIR) / "swapaa.py" if config.CHIMERA_SCRIPTS_DIR else None,
                root_dir / "chimera_scripts" / "swapaa.py",
                Path("/home/guido/miniforge3/envs/colbuilder/lib/python3.9/site-packages/colbuilder/chimera_scripts/swapaa.py")
            ]
            
            for p in search_paths:
                if p and p.exists():
                    swapaa_script = p
                    break

            if not swapaa_script:
                LOG.error("Cannot find swapaa.py script.")
                return False

            # Keep execution string intact for Chimera, but log a concise command line
            cmd = f'chimera --nogui --silent --script "{swapaa_script} {replace_file_path} {type_dir_path}"'

            if getattr(config, "_replacement_verbose", True):
                LOG.info(
                    "Running command: chimera --nogui --silent --script %s",
                    swapaa_script,
                )
            LOG.debug("swapaa args: %s %s", replace_file_path, type_dir_path)
            result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            
            if result.returncode != 0:
                LOG.error(f"Chimera Error: {result.stderr}")
                return False

            if getattr(config, "_replacement_verbose", True):
                LOG.info("Chimera swapaa executed via %s", swapaa_script)
            return True

        except Exception as e:
            LOG.error(f"Chimera execution exception: {e}")
            return False

    async def replace_direct(self, config: ColbuilderConfig, temp_dir: Path):
        """
        Direct replacement workflow using an existing PDB (geometry_generation=False).

        This allows running ratio_replace/replace_file without regenerating geometry.
        """
        input_pdb = Path(config.replace_file) if config.replace_file else None
        if not input_pdb or not input_pdb.exists():
            raise GeometryGenerationError(
                message=f"Input PDB file not found: {input_pdb}",
                error_code="GEO_ERR_004",
            )

        working_dir_root = Path(config.working_directory).resolve()
        temp_dir.mkdir(parents=True, exist_ok=True)

        local_pdb = temp_dir / input_pdb.name
        if input_pdb.resolve() != local_pdb.resolve():
            shutil.copy2(input_pdb, local_pdb)

        type_dir = temp_dir / "NC"
        type_dir.mkdir(parents=True, exist_ok=True)

        # Split the input PDB into per-model caps files
        self._split_pdb_into_models(local_pdb, type_dir)
        system_from_caps: Optional[System] = None

        # Copy manual_replacements.txt if present
        candidate_path = working_dir_root / "manual_replacements.txt"
        if not candidate_path.exists():
            candidate_path = Path.cwd() / "manual_replacements.txt"
        target_file_path = temp_dir / "manual_replacements.txt"
        if candidate_path.exists() and candidate_path.resolve() != target_file_path.resolve():
            shutil.copy2(candidate_path, target_file_path)

        manual_list: List[str] = []
        ratio_requested = config.ratio_replace is not None and float(config.ratio_replace) > 0
        generated_from_ratio = False

        if getattr(config, "manual_replacements", None):
            manual_list = [
                str(instr).strip()
                for instr in config.manual_replacements
                if str(instr).strip()
            ]

        if not manual_list and ratio_requested:
            records = self._load_crosslinks_from_models(type_dir)
            manual_list = self._build_ratio_replacements_from_records(
                records=records,
                ratio_replace=float(config.ratio_replace),
                fibril_length=getattr(config, "fibril_length", 0.0),
                scope=getattr(config, "ratio_replace_scope", "non_enzymatic"),
            )
            generated_from_ratio = bool(manual_list)

        if not manual_list and target_file_path.exists() and not ratio_requested:
            with open(target_file_path, "r") as f:
                manual_list = [
                    line.strip()
                    for line in f
                    if line.strip() and not line.startswith("#")
                ]

        if not manual_list and hasattr(config, "__dict__"):
            manual_list = config.__dict__.get("manual_replacements") or []

        if not manual_list:
            LOG.warning(
                f"{Fore.YELLOW}No replacement instructions generated. Copying input to output without changes.{Style.RESET_ALL}"
            )
            output_pdb = temp_dir / f"{config.output or 'output'}.pdb"
            if output_pdb.resolve() != local_pdb.resolve():
                shutil.copy2(local_pdb, output_pdb)
            else:
                output_pdb = local_pdb
            return None, output_pdb

        if generated_from_ratio:
            try:
                with open(target_file_path, "w") as f:
                    for instruction in manual_list:
                        f.write(f"{instruction}\n")
                LOG.info(
                    "Saved ratio-based replacement instructions to %s for transparency.",
                    target_file_path,
                )
            except Exception as e:
                LOG.warning("Could not persist ratio-based instructions: %s", e)

        replace_file = temp_dir / "replace.txt"
        with open(replace_file, "w") as f:
            for instruction in manual_list:
                clean_instr = str(instruction).strip().strip('"').strip("'")
                f.write(f"{clean_instr}\n")

        success = await self._run_chimera_command(
            config,
            str(replace_file),
            type_dir,
            working_dir_root,
        )

        if not success:
            raise GeometryGenerationError(
                message="Chimera replacement failed in direct mode",
                error_code="GEO_ERR_004",
            )

        output_pdb = temp_dir / f"{config.output or 'output'}.pdb"
        self._combine_caps_to_pdb(
            source_dir=type_dir,
            output_pdb=output_pdb,
            template_pdb=local_pdb,
        )

        try:
            system_from_caps = self._categorize_caps_and_build_system(
                base_dir=temp_dir, source_dir=type_dir, reference_pdb=local_pdb
            )
            if system_from_caps:
                self._recompute_connectivity_from_caps(system_from_caps, temp_dir)
        except Exception as e:
            LOG.warning(f"Failed to build system/connectivity from replaced caps: {e}")

        return system_from_caps, output_pdb

    def _load_crosslinks_from_models(self, type_dir: Path) -> List[Dict[str, Any]]:
        """Parse caps models to extract crosslink residues and their positions."""
        records: List[Dict[str, Any]] = []
        for pdb_file in type_dir.glob("*.caps.pdb"):
            try:
                model_id = int(pdb_file.stem.split(".")[0])
            except ValueError:
                continue

            residue_atoms: Dict[Tuple[str, str], List[List[float]]] = {}

            with open(pdb_file, "r") as fh:
                for line in fh:
                    if not line.startswith(("ATOM", "HETATM")) or len(line) < 54:
                        continue
                    resname = line[17:20].strip()
                    if resname not in self._crosslink_resnames:
                        continue
                    resid = line[22:26].strip()
                    chain = line[21].strip() or "A"
                    key = (resname, resid, chain)
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                    except ValueError:
                        continue
                    residue_atoms.setdefault(key, []).append([x, y, z])

            for (resname, resid, chain), atoms in residue_atoms.items():
                if not atoms:
                    continue
                x_avg = sum(a[0] for a in atoms) / len(atoms)
                y_avg = sum(a[1] for a in atoms) / len(atoms)
                z_avg = sum(a[2] for a in atoms) / len(atoms)
                records.append(
                    {
                        "model_id": model_id,
                        "resname": resname,
                        "resid": resid,
                        "chain": chain,
                        "position": [x_avg, y_avg, z_avg],
                    }
                )

        return records

    def _load_connect_groups(self, connect_file: Path) -> List[List[int]]:
        """Parse connect_from_colbuilder.txt into groups of connected model ids."""
        if not connect_file or not connect_file.exists():
            return []

        groups: List[List[int]] = []
        try:
            with open(connect_file, "r") as fh:
                for line in fh:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    parts = line.split(";")[0].strip().split()
                    if not parts:
                        continue
                    model_ids: List[int] = []
                    for token in parts:
                        token = token.strip()
                        if not token:
                            continue
                        stem = token.split(".")[0]
                        try:
                            model_ids.append(int(stem))
                        except ValueError:
                            continue
                    if model_ids:
                        groups.append(sorted(set(model_ids)))
        except Exception as e:
            LOG.warning("Failed to parse connect file %s: %s", connect_file, e)
            return []

        return groups

    def _parse_term_combination(
        self, combo: Optional[str]
    ) -> Set[Tuple[str, str]]:
        """Parse a term combination like '9.C - 5.B - 944.B' into resid/chain pairs."""
        if not combo:
            return set()
        pairs: Set[Tuple[str, str]] = set()
        for token in combo.split("-"):
            token = token.strip()
            if not token or "." not in token:
                continue
            resid, chain = token.split(".", 1)
            resid = resid.strip()
            chain = chain.strip()
            if resid and chain:
                pairs.add((resid, chain))
        return pairs

    async def _apply_manual_replacements_to_dir(
        self,
        source_dir: Path,
        dest_dir: Path,
        manual_list: List[str],
        working_dir_root: Path,
        config: ColbuilderConfig,
        system: Any,
    ) -> bool:
        """Apply manual replacements to caps in a dedicated directory."""
        if not manual_list:
            return False

        dest_dir.mkdir(parents=True, exist_ok=True)
        stale_manual = dest_dir / "manual_replacements.txt"
        if stale_manual.exists():
            try:
                stale_manual.unlink()
            except Exception:
                pass

        type_dir_candidates = [source_dir]
        for sub in ["DT", "TD", "D", "T", "NC", "DY"]:
            cand = source_dir / sub
            if cand.exists():
                type_dir_candidates.append(cand)

        try:
            model_zero = system.get_model(model_id=0.0)
            system_type = model_zero.type if hasattr(model_zero, "type") else "D"
        except Exception:
            system_type = "D"
        type_dir = dest_dir / system_type
        type_dir.mkdir(parents=True, exist_ok=True)

        caps_by_model: Dict[int, Path] = {}
        for search_dir in type_dir_candidates:
            for pdb_file in search_dir.glob("*.caps.pdb"):
                try:
                    model_id = int(pdb_file.stem.split(".")[0])
                except ValueError:
                    continue
                if model_id not in caps_by_model:
                    caps_by_model[model_id] = pdb_file

        for model_id, source_caps in sorted(caps_by_model.items()):
            dest_caps = type_dir / source_caps.name
            if source_caps.resolve() != dest_caps.resolve():
                shutil.copy2(source_caps, dest_caps)

        connect_file = None
        for candidate in [
            source_dir / "connect_from_colbuilder.txt",
            source_dir / "connect_from_colbuilder",
        ]:
            if candidate.exists():
                connect_file = candidate
                break
        if connect_file and connect_file.exists():
            if connect_file.suffix != ".txt":
                connect_file_txt = connect_file.with_suffix(".txt")
                if connect_file_txt.exists():
                    connect_file = connect_file_txt
            shutil.copy2(connect_file, dest_dir / "connect_from_colbuilder.txt")
        else:
            try:
                connector = Connect(system=system)
                connect_file = dest_dir / "connect_from_colbuilder"
                connector.write_connect(system=system, connect_file=connect_file)
            except Exception as e:
                LOG.warning("Could not rebuild connectivity for auto-fix: %s", e)

        replace_file = dest_dir / "replace.txt"
        with open(replace_file, "w") as f:
            for instruction in manual_list:
                clean_instr = instruction.strip().strip('"').strip("'")
                f.write(f"{clean_instr}\n")

        success = await self._run_chimera_command(
            config=config,
            replace_file_path=str(replace_file),
            type_dir_path=type_dir,
            root_dir=working_dir_root,
        )

        if not success:
            LOG.warning("Chimera execution failed for auto-fix manual replacements.")
            return False

        return True

    def _build_ratio_replacements_from_connect(
        self,
        records: List[Dict[str, Any]],
        connect_groups: List[List[int]],
        ratio_replace: float,
        scope: str = "non_enzymatic",
        config: Optional[ColbuilderConfig] = None,
    ) -> List[str]:
        """Generate replacement instructions using connect groups and parsed caps records."""
        if not records or not connect_groups or ratio_replace <= 0:
            return []

        scope = (scope or "non_enzymatic").strip().lower()
        if scope not in {"enzymatic", "non_enzymatic", "all"}:
            LOG.warning(
                "Unknown ratio_replace scope '%s'; defaulting to 'non_enzymatic'", scope
            )
            scope = "non_enzymatic"

        records_by_model: Dict[int, List[Dict[str, Any]]] = {}
        for rec in records:
            try:
                model_id = int(rec.get("model_id", -1))
            except (ValueError, TypeError):
                continue
            records_by_model.setdefault(model_id, []).append(rec)

        n_term_pairs = self._parse_term_combination(
            getattr(config, "n_term_combination", None) if config else None
        )
        c_term_pairs = self._parse_term_combination(
            getattr(config, "c_term_combination", None) if config else None
        )

        entity_pool: List[Tuple[str, List[Dict[str, Any]]]] = []
        for group in connect_groups:
            group_records: List[Dict[str, Any]] = []
            for model_id in group:
                group_records.extend(records_by_model.get(model_id, []))
            if not group_records:
                continue

            by_type: Dict[str, List[Dict[str, Any]]] = {}
            for rec in group_records:
                by_type.setdefault(rec.get("resname", ""), []).append(rec)

            pairs: List[List[Dict[str, Any]]] = []
            used_records: Set[Tuple[int, str, str, str]] = set()

            for donors, acceptors in PAIRING_RULES:
                donor_list: List[Dict[str, Any]] = []
                acceptor_list: List[Dict[str, Any]] = []
                for res in donors:
                    donor_list.extend(by_type.get(res, []))
                for res in acceptors:
                    acceptor_list.extend(by_type.get(res, []))

                if not donor_list or not acceptor_list:
                    continue

                used_donors: Set[int] = set()
                used_acceptors: Set[int] = set()

                for i, donor in enumerate(donor_list):
                    if i in used_donors:
                        continue
                    best_idx = -1
                    best_dist = float("inf")
                    for j, acceptor in enumerate(acceptor_list):
                        if j in used_acceptors:
                            continue
                        dist = self._calculate_distance(
                            donor["position"], acceptor["position"]
                        )
                        if dist <= 10.0 and dist < best_dist:
                            best_idx = j
                            best_dist = dist
                    if best_idx != -1:
                        pair = [donor, acceptor_list[best_idx]]
                        pairs.append(pair)
                        used_donors.add(i)
                        used_acceptors.add(best_idx)

            for pair in pairs:
                for rec in pair:
                    used_records.add(
                        (
                            int(float(rec.get("model_id", -1))),
                            rec.get("resid", ""),
                            rec.get("chain", ""),
                            rec.get("resname", ""),
                        )
                    )

            for pair in pairs:
                resnames = {rec.get("resname", "") for rec in pair}
                if resnames.issubset(NON_ENZYMATIC_PAIRED):
                    scope_tag = "non_enzymatic"
                else:
                    scope_tag = None
                    pair_pairs = {(r.get("resid", ""), r.get("chain", "")) for r in pair}
                    if n_term_pairs and pair_pairs.issubset(n_term_pairs):
                        scope_tag = "enzymatic_n"
                    elif c_term_pairs and pair_pairs.issubset(c_term_pairs):
                        scope_tag = "enzymatic_c"
                if scope_tag:
                    entity_pool.append((scope_tag, pair))

            for scope_tag, term_pairs in (
                ("enzymatic_n", n_term_pairs),
                ("enzymatic_c", c_term_pairs),
            ):
                if not term_pairs:
                    continue
                term_records = [
                    rec
                    for rec in group_records
                    if rec.get("resname", "") in ENZYMATIC_SINGLETONS
                    and (rec.get("resid", ""), rec.get("chain", "")) in term_pairs
                ]
                term_found = {(r.get("resid", ""), r.get("chain", "")) for r in term_records}
                if term_found and term_found.issuperset(term_pairs):
                    entity_pool.append((scope_tag, term_records))

        if not entity_pool:
            return []

        if scope == "enzymatic":
            eligible_entities = [
                e
                for e in entity_pool
                if e[0] in {"enzymatic_n", "enzymatic_c"}
            ]
        elif scope == "non_enzymatic":
            eligible_entities = [e for e in entity_pool if e[0] == "non_enzymatic"]
        else:
            eligible_entities = entity_pool

        if not eligible_entities:
            return []

        random.seed(int(time.time()))

        selected_entities: List[Tuple[str, List[Dict[str, Any]]]] = []
        if scope == "all":
            buckets = {
                "non_enzymatic": [e for e in eligible_entities if e[0] == "non_enzymatic"],
                "enzymatic_c": [e for e in eligible_entities if e[0] == "enzymatic_c"],
                "enzymatic_n": [e for e in eligible_entities if e[0] == "enzymatic_n"],
            }
            for key, bucket in buckets.items():
                if not bucket:
                    continue
                random.shuffle(bucket)
                target = max(1, math.ceil(len(bucket) * ratio_replace / 100.0))
                selected_entities.extend(bucket[:target])
        elif scope == "enzymatic":
            buckets = {
                "enzymatic_c": [e for e in eligible_entities if e[0] == "enzymatic_c"],
                "enzymatic_n": [e for e in eligible_entities if e[0] == "enzymatic_n"],
            }
            for key, bucket in buckets.items():
                if not bucket:
                    continue
                random.shuffle(bucket)
                target = max(1, math.ceil(len(bucket) * ratio_replace / 100.0))
                selected_entities.extend(bucket[:target])
        else:
            random.shuffle(eligible_entities)
            target = max(1, math.ceil(len(eligible_entities) * ratio_replace / 100.0))
            selected_entities = eligible_entities[:target]

        seen: Set[str] = set()
        instructions: List[str] = []
        for _scope_tag, entity_records in selected_entities:
            for rec in entity_records:
                new_res = REPLACEMENT_MAP.get(rec["resname"], DEFAULT_REPLACEMENT)
                instruction = (
                    f"{int(float(rec['model_id']))}.caps.pdb {new_res} {rec['resid']} {rec['chain']}"
                )
                if instruction not in seen:
                    seen.add(instruction)
                    instructions.append(instruction)

        def sort_key(line: str) -> Tuple:
            parts = line.split()
            model_part = parts[0] if parts else ""
            try:
                model_idx = int(model_part.split(".")[0])
            except (ValueError, IndexError):
                model_idx = model_part
            resid_part = parts[2] if len(parts) > 2 else ""
            try:
                resid_idx = int(resid_part)
            except ValueError:
                resid_idx = resid_part
            chain_part = parts[3] if len(parts) > 3 else ""
            return (model_idx, chain_part, resid_idx, parts[1] if len(parts) > 1 else "")

        instructions.sort(key=sort_key)
        return instructions

    def _build_ratio_replacements_from_records(
        self,
        records: List[Dict[str, Any]],
        ratio_replace: float,
        fibril_length: float,
        scope: str = "non_enzymatic",
    ) -> List[str]:
        """Generate replacement instructions using parsed crosslink records."""
        if not records or ratio_replace <= 0:
            return []

        scope = (scope or "non_enzymatic").strip().lower()
        if scope not in {"enzymatic", "non_enzymatic", "all"}:
            LOG.warning(
                "Unknown ratio_replace scope '%s'; defaulting to 'non_enzymatic'", scope
            )
            scope = "non_enzymatic"

        if scope == "enzymatic":
            target_resnames = ENZYMATIC_SINGLETONS.copy()
        elif scope == "non_enzymatic":
            target_resnames = NON_ENZYMATIC_PAIRED.copy()
        else:
            target_resnames = NON_ENZYMATIC_PAIRED.union(ENZYMATIC_SINGLETONS)

        z_bounds = self._calculate_bounds_from_records(records, fibril_length)
        filtered = [
            r
            for r in records
            if r["resname"] in target_resnames
            and z_bounds[0] <= self._get_component(r["position"], 2) <= z_bounds[1]
        ]

        by_type: Dict[str, List[Dict[str, Any]]] = {res: [] for res in target_resnames}
        for rec in filtered:
            by_type.setdefault(rec["resname"], []).append(rec)

        pairs: List[List[Dict[str, Any]]] = []
        singles: List[Dict[str, Any]] = []

        for donors, acceptors in PAIRING_RULES:
            if not (target_resnames.issuperset(set(donors)) and target_resnames.issuperset(set(acceptors))):
                continue

            donor_list: List[Dict[str, Any]] = []
            acceptor_list: List[Dict[str, Any]] = []
            for res in donors:
                donor_list.extend(by_type.get(res, []))
            for res in acceptors:
                acceptor_list.extend(by_type.get(res, []))

            if not donor_list or not acceptor_list:
                continue

            used_donors: Set[int] = set()
            used_acceptors: Set[int] = set()

            for i, donor in enumerate(donor_list):
                if i in used_donors:
                    continue
                best_idx = -1
                best_dist = float("inf")
                for j, acceptor in enumerate(acceptor_list):
                    if j in used_acceptors:
                        continue
                    # Use a looser cutoff here (10 Å) matching legacy direct ratio replacement
                    dist = self._calculate_distance(donor["position"], acceptor["position"])
                    if dist <= 10.0 and dist < best_dist:
                        best_idx = j
                        best_dist = dist
                if best_idx != -1:
                    pairs.append([donor, acceptor_list[best_idx]])
                    used_donors.add(i)
                    used_acceptors.add(best_idx)

        # PYD trios (LYX/LY2/LY3) when enzymatic scope is active
        if scope in {"enzymatic", "all"}:
            lyx_list = by_type.get("LYX", [])
            ly2_list = by_type.get("LY2", [])
            ly3_list = by_type.get("LY3", [])
            pyd_trios = self._build_pyd_trios(lyx_list, ly2_list, ly3_list)

        num_pair_to_replace = (
            min(len(pairs), max(1, math.ceil(len(pairs) * ratio_replace / 100.0)))
            if pairs
            else 0
        )

        num_trio_to_replace = (
            min(len(pyd_trios), max(1, math.ceil(len(pyd_trios) * ratio_replace / 100.0)))
            if pyd_trios
            else 0
        )

        if scope in {"enzymatic", "all"}:
            for resname in target_resnames:
                if resname in PAIRED_RESIDUES:
                    continue
                singles.extend(by_type.get(resname, []))
            trio_members = {
                (r.get("resid", ""), r.get("chain", ""), r.get("model_id", 0.0))
                for trio in pyd_trios
                for r in trio
            }
            singles = [
                s
                for s in singles
                if (s.get("resid", ""), s.get("chain", ""), s.get("model_id", 0.0)) not in trio_members
            ]
        # Fallback: if we found no valid pairs or trios but still have candidates, treat them as singles
        if not pairs and not pyd_trios and filtered:
            singles = filtered.copy()

        num_single_to_replace = (
            min(len(singles), max(1, math.ceil(len(singles) * ratio_replace / 100.0)))
            if singles
            else 0
        )

        random.seed(int(time.time()))
        random.shuffle(pairs)
        random.shuffle(singles)
        random.shuffle(pyd_trios)
        selected_pairs = pairs[:num_pair_to_replace]
        selected_trios = pyd_trios[:num_trio_to_replace]
        selected_singles = singles[:num_single_to_replace]

        seen: Set[str] = set()
        instructions: List[str] = []

        for pair in selected_pairs:
            for rec in pair:
                new_res = REPLACEMENT_MAP.get(rec["resname"], DEFAULT_REPLACEMENT)
                instruction = f"{int(float(rec['model_id']))}.caps.pdb {new_res} {rec['resid']} {rec['chain']}"
                if instruction not in seen:
                    seen.add(instruction)
                    instructions.append(instruction)

        for trio in selected_trios:
            for rec in trio:
                new_res = REPLACEMENT_MAP.get(rec["resname"], DEFAULT_REPLACEMENT)
                instruction = f"{int(float(rec['model_id']))}.caps.pdb {new_res} {rec['resid']} {rec['chain']}"
                if instruction not in seen:
                    seen.add(instruction)
                    instructions.append(instruction)

        for rec in selected_singles:
            new_res = REPLACEMENT_MAP.get(rec["resname"], DEFAULT_REPLACEMENT)
            instruction = f"{int(float(rec['model_id']))}.caps.pdb {new_res} {rec['resid']} {rec['chain']}"
            if instruction not in seen:
                seen.add(instruction)
                instructions.append(instruction)

        def sort_key(line: str) -> Tuple:
            parts = line.split()
            model_part = parts[0] if parts else ""
            try:
                model_idx = int(model_part.split(".")[0])
            except (ValueError, IndexError):
                model_idx = model_part
            resid_part = parts[2] if len(parts) > 2 else ""
            try:
                resid_idx = int(resid_part)
            except ValueError:
                resid_idx = resid_part
            chain_part = parts[3] if len(parts) > 3 else ""
            return (model_idx, chain_part, resid_idx, parts[1] if len(parts) > 1 else "")

        instructions.sort(key=sort_key)
        LOG.info(
            "Selected %d residues for replacement (pairs: %d residues, PYD trios: %d residues, singles: %d) out of %d candidates (scope=%s, ratio=%s%%).",
            len(instructions),
            len(selected_pairs) * 2,
            len(selected_trios) * 3,
            len(selected_singles),
            len(pairs) * 2 + len(pyd_trios) * 3 + len(singles),
            scope,
            ratio_replace,
        )
        return instructions

    def _categorize_caps_and_build_system(
        self, base_dir: Path, source_dir: Path, reference_pdb: Path
    ) -> Optional[System]:
        """
        Classify caps models by remaining marker type and build a lightweight System.

        The ratio_replace direct workflow initially stores all fragments under NC/.
        To enable downstream connectivity/topology, we scan the caps files, infer
        their marker composition (D/T/DT/NC), and place the entire set into a
        single folder whose name captures every marker type present (NC is only
        used when no markers remain).
        """
        search_dirs = [source_dir]
        for sub in ["D", "T", "NC", "DT", "TD", "DY"]:
            sub_dir = base_dir / sub
            if sub_dir.exists() and sub_dir not in search_dirs:
                search_dirs.append(sub_dir)

        caps_candidates: List[Path] = []
        seen: Set[Path] = set()
        for directory in search_dirs:
            for pdb_file in sorted(directory.glob("*.caps.pdb")):
                resolved = pdb_file.resolve()
                if resolved in seen:
                    continue
                seen.add(resolved)
                caps_candidates.append(pdb_file)

        if not caps_candidates:
            LOG.warning("No caps files found to categorise in %s", base_dir)
            return None

        system = System()
        try:
            from .crystal import Crystal

            system.crystal = Crystal(pdb=str(reference_pdb))
        except Exception as e:
            LOG.warning("Could not attach crystal to replacement system: %s", e)

        temp_models: List[Tuple[Path, Model]] = []
        found_types: Set[str] = set()
        for caps_file in caps_candidates:
            try:
                model_id = float(caps_file.stem.split(".")[0])
            except ValueError:
                LOG.debug("Skipping caps file without numeric id: %s", caps_file)
                continue

            model_obj = Model(
                id=model_id, transformation=[0.0, 0.0, 0.0], pdb_file=str(caps_file)
            )
            temp_models.append((caps_file, model_obj))
            if model_obj.type:
                found_types.add(model_obj.type)

        canonical_types = {t for t in found_types if t and t != "NC"}
        if not canonical_types:
            aggregated_type = "NC"
        else:
            aggregated_type = "".join(sorted(canonical_types))

        dest_dir = base_dir / aggregated_type
        dest_dir.mkdir(parents=True, exist_ok=True)

        for caps_file, model_obj in temp_models:
            dest_path = dest_dir / caps_file.name
            if caps_file.resolve() != dest_path.resolve():
                try:
                    shutil.move(caps_file, dest_path)
                except Exception:
                    shutil.copy2(caps_file, dest_path)

            model_obj.pdb_file = str(dest_path)
            model_obj.type = aggregated_type
            system.add_model(model_obj)

        system.get_models()

        # Clean up any stale type directories that are now empty
        for sub in ["D", "T", "NC", "DT", "TD", "DY"]:
            sub_dir = base_dir / sub
            if sub_dir == dest_dir or not sub_dir.exists():
                continue
            try:
                if next(sub_dir.iterdir(), None) is None:
                    shutil.rmtree(sub_dir)
            except Exception:
                continue

        return system

    def _recompute_connectivity_from_caps(self, system: System, output_dir: Path) -> None:
        """Rebuild connect_from_colbuilder.txt from the categorised caps files."""
        connector = Connect(system=system)
        new_connect = connector.run_connect(system=system)

        if new_connect:
            for mid, conns in new_connect.items():
                try:
                    model_obj = system.get_model(model_id=mid)
                except Exception:
                    model_obj = None
                if model_obj:
                    model_obj.connect = conns

            connect_file_path = output_dir / "connect_from_colbuilder.txt"
            connector.write_connect(system=system, connect_file=connect_file_path)
            LOG.info(f"Connectivity written to {connect_file_path}")
        else:
            LOG.warning("Connectivity recomputation returned empty; forcing self-connect entries.")
            fallback_connect = {mid: [mid] for mid in system.get_models()}
            for mid, conns in fallback_connect.items():
                try:
                    model_obj = system.get_model(model_id=mid)
                except Exception:
                    continue
                model_obj.connect = conns
            connector.connect = fallback_connect
            connect_file_path = output_dir / "connect_from_colbuilder.txt"
            connector.write_connect(system=system, connect_file=connect_file_path)
            LOG.info(f"Connectivity (self-only) written to {connect_file_path}")

    def _calculate_bounds_from_records(
        self, records: List[Dict[str, Any]], fibril_length: float
    ) -> Tuple[float, float]:
        """Rudimentary z-bounds based on parsed crosslink positions."""
        if not records:
            return (-math.inf, math.inf)

        z_vals = [self._get_component(r["position"], 2) for r in records if "position" in r]
        if not z_vals:
            return (-math.inf, math.inf)

        z_min = min(z_vals) - 500.0
        z_max = max(z_vals) + 500.0

        if z_max - z_min > 5000.0 and fibril_length:
            center = (z_min + z_max) / 2.0
            fibril_len_ang = fibril_length * 10.0
            return (center - fibril_len_ang * 5, center + fibril_len_ang * 5)

        return (z_min, z_max)

    def _build_pyd_trios(
        self,
        lyx_list: List[Dict[str, Any]],
        ly2_list: List[Dict[str, Any]],
        ly3_list: List[Dict[str, Any]],
    ) -> List[List[Dict[str, Any]]]:
        """
        Group LYX/LY2/LY3 into PYD trios based on proximity (legacy cutoff 15 Å).
        """
        trios: List[List[Dict[str, Any]]] = []
        used_ly2: Set[int] = set()
        used_ly3: Set[int] = set()

        for lyx in lyx_list:
            best_ly2 = None
            best_ly2_dist = float("inf")
            for idx, ly2 in enumerate(ly2_list):
                if idx in used_ly2:
                    continue
                try:
                    dist = self._calculate_distance(lyx["crosslink"].position if "crosslink" in lyx else lyx["position"],
                                                    ly2["crosslink"].position if "crosslink" in ly2 else ly2["position"])
                except Exception:
                    continue
                if dist <= PYD_THRESHOLD and dist < best_ly2_dist:
                    best_ly2 = (idx, ly2)
                    best_ly2_dist = dist

            best_ly3 = None
            best_ly3_dist = float("inf")
            for jdx, ly3 in enumerate(ly3_list):
                if jdx in used_ly3:
                    continue
                try:
                    dist = self._calculate_distance(lyx["crosslink"].position if "crosslink" in lyx else lyx["position"],
                                                    ly3["crosslink"].position if "crosslink" in ly3 else ly3["position"])
                except Exception:
                    continue
                if dist <= PYD_THRESHOLD and dist < best_ly3_dist:
                    best_ly3 = (jdx, ly3)
                    best_ly3_dist = dist

            if best_ly2 and best_ly3:
                used_ly2.add(best_ly2[0])
                used_ly3.add(best_ly3[0])
                trios.append([lyx, best_ly2[1], best_ly3[1]])

        return trios

    def _combine_caps_to_pdb(self, source_dir: Path, output_pdb: Path, template_pdb: Path) -> None:
        """Combine modified caps files into a single PDB."""
        with open(output_pdb, "w") as out:
            try:
                with open(template_pdb, "r") as tmpl:
                    first_line = tmpl.readline().strip()
                    if first_line.startswith("CRYST1"):
                        out.write(first_line + "\n")
                    else:
                        out.write("REMARK   Generated by colbuilder replacement\n")
            except Exception:
                out.write("REMARK   Generated by colbuilder replacement\n")

            pdb_files = sorted(
                [f for f in source_dir.glob("*.caps.pdb")],
                key=lambda x: int(x.stem.split(".")[0]),
            )

            for pdb_file in pdb_files:
                with open(pdb_file, "r") as fh:
                    for line in fh:
                        if line.startswith(("ATOM", "HETATM", "TER")):
                            if line.startswith("HETATM"):
                                line = "ATOM  " + line[6:]
                            out.write(line)
            out.write("END\n")

    def _split_pdb_into_models(self, input_pdb: Path, type_dir: Path) -> int:
        """
        Split a Colbuilder-generated PDB into per-model caps files.
        Mirrors the legacy behaviour used by ratio replacement in standard_geo.
        """
        with open(input_pdb, "r") as fh:
            all_lines = fh.readlines()

        cryst_line = None
        for line in all_lines:
            if line.startswith("CRYST1"):
                cryst_line = line
                break

        atom_lines = [ln for ln in all_lines if ln.startswith(("ATOM", "HETATM", "TER"))]
        models: List[List[str]] = []
        current: List[str] = []
        chain_sequence: List[str] = []

        for line in atom_lines:
            if line.startswith(("ATOM", "HETATM")):
                chain = line[21]
                if not chain_sequence or chain != chain_sequence[-1]:
                    chain_sequence.append(chain)
                    if len(chain_sequence) > 3 and chain_sequence[-4:] == ["A", "B", "C", "A"]:
                        models.append(current)
                        current = []
                        chain_sequence = ["A"]
            current.append(line)

        if current:
            models.append(current)

        if not models:
            models = [atom_lines]

        type_dir.mkdir(parents=True, exist_ok=True)
        for idx, model_lines in enumerate(models):
            caps_file = type_dir / f"{idx}.caps.pdb"
            with open(caps_file, "w") as fh:
                if cryst_line:
                    fh.write(cryst_line)
                fh.writelines(model_lines)
                if not model_lines[-1].startswith("TER"):
                    fh.write("TER\n")

        return len(models)


async def replace_in_system(system: Any, config: ColbuilderConfig) -> Any:
    replacer = CrosslinkReplacer()
    return await replacer.replace_in_system(system, config)


async def direct_replace_geometry(config: ColbuilderConfig) -> None:
    replacer = CrosslinkReplacer()
    temp_dir = Path(config.working_directory) / ".tmp" / "replacement_direct"
    temp_dir.mkdir(parents=True, exist_ok=True)
    await replacer.replace_direct(config, temp_dir)
