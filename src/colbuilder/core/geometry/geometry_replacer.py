"""
Colbuilder Crosslink Replacement Module

This module provides a unified approach to replacing crosslinks with standard amino acids
in a collagen microfibril, supporting both system-based and direct replacement approaches.
"""

import os
import random
import time
import math
import traceback
import subprocess
import shutil
import numpy as np
from pathlib import Path
from typing import Optional, Dict, Any, List, Union, Tuple, cast, Set
from colorama import init, Fore, Style

from ..utils.exceptions import GeometryGenerationError
from ..utils.logger import setup_logger
from ..utils.config import ColbuilderConfig
from .system import System
from .crystal import Crystal

LOG = setup_logger(__name__)


class CrosslinkReplacer:
    """Handles replacement of crosslinks in collagen systems."""

    def __init__(self):
        """Initialize the CrosslinkReplacer."""
        self.file_manager = None

    async def replace(
        self, system: Optional[System], config: ColbuilderConfig, temp_dir: Path
    ) -> Tuple[Optional[System], Optional[Path]]:
        """
        Replace crosslinks in a system or perform direct replacement.

        Args:
            system: Optional system to modify. If None, performs direct replacement
            config: Configuration settings
            temp_dir: Temporary directory for file operations

        Returns:
            Tuple[Optional[System], Optional[Path]]: (modified system, output PDB path)
        """
        LOG.debug(f"Starting replacement with temp_dir: {temp_dir}")

        if system is None:
            LOG.info("Using direct replacement mode")
            return await self.replace_direct(config, temp_dir)
        else:
            LOG.info("Using system modification mode")
            modified_system = await self.replace_in_system(system, config, temp_dir)
            output_pdb = temp_dir / f"{config.output or 'output'}.pdb"

            modified_system.write_pdb(
                pdb_out=output_pdb,
                fibril_length=config.fibril_length,
                temp_dir=temp_dir,
            )

            return modified_system, output_pdb

    async def replace_in_system(
        self, system: Any, config: ColbuilderConfig, temp_dir: Optional[Path] = None
    ) -> Any:
        """
        Replace crosslinks in the provided system according to configuration settings.

        This is the main entry point for system-based replacement. It analyzes the system,
        selects crosslinks to replace based on the specified ratio, and performs the
        replacement using Chimera.

        Args:
            system: System containing the models to be modified
            config: Configuration for replacement including ratio and other settings
            temp_dir: Optional temporary directory for file operations

        Returns:
            Modified system with replaced crosslinks

        Raises:
            GeometryGenerationError: If replacement fails for any reason
        """
        try:
            if temp_dir is None:
                working_dir = Path(config.working_directory)
                temp_dir = working_dir / ".tmp" / "replacement"
                temp_dir.mkdir(parents=True, exist_ok=True)
                LOG.debug(f"Created temporary directory for replacement: {temp_dir}")

            # Direct replacement mode handling - only if system is None
            if (
                config.replace_file
                and os.path.exists(config.replace_file)
                and system is None
            ):
                LOG.info(f"Using external replacement file: {config.replace_file}")
                with open(config.replace_file, "r") as f:
                    first_line = f.readline().strip()
                    if first_line.startswith(("ATOM", "CRYST1", "HETATM")):
                        system, _ = await self.replace_direct(config, temp_dir)
                        return system

                if system is None:
                    raise GeometryGenerationError(
                        message="System is None but replacement instructions file was provided. Cannot continue.",
                        error_code="GEO_ERR_004",
                    )

            model_count = (
                system.get_size() if hasattr(system, "get_size") else "unknown"
            )
            crosslink_count = 0
            crosslink_types = set()
            models_with_crosslinks = 0

            for model_id in system.get_models():
                model = system.get_model(model_id=model_id)
                if hasattr(model, "crosslink") and model.crosslink:
                    if len(model.crosslink) > 0:
                        models_with_crosslinks += 1
                        for crosslink in model.crosslink:
                            crosslink_count += 1
                            if hasattr(crosslink, "resname"):
                                crosslink_types.add(crosslink.resname)

            if crosslink_types:
                LOG.info(
                    f"{Fore.BLUE}Crosslink residue types: {', '.join(crosslink_types)}{Style.RESET_ALL}"
                )

            if crosslink_count == 0:
                LOG.warning("No crosslinks found in the system - nothing to replace!")
                return system

            if not hasattr(system, "config"):
                system.config = config

            # Calculate bounds for replacement region
            z_bounds = self._calculate_fibril_bounds(system, config.fibril_length)

            # Select crosslinks for replacement
            self._select_replacements(system, config.ratio_replace, z_bounds)
            to_replace_count = 0
            to_replace_types = {}

            for model_id in system.get_models():
                model = system.get_model(model_id=model_id)
                if hasattr(model, "crosslink") and model.crosslink:
                    for crosslink in model.crosslink:
                        if hasattr(crosslink, "state") and crosslink.state == "replace":
                            to_replace_count += 1
                            if hasattr(crosslink, "resname"):
                                to_replace_types[crosslink.resname] = (
                                    to_replace_types.get(crosslink.resname, 0) + 1
                                )

            if to_replace_types:
                LOG.debug("Replacing:")
                for resname, count in to_replace_types.items():
                    LOG.debug(f"  - {resname}: {count}")

            if to_replace_count == 0:
                LOG.warning(
                    "No crosslinks selected for replacement - skipping Chimera step"
                )
                return system

            # Write replacement instructions to file
            replace_file = temp_dir / "replace.txt"
            self._write_replace_file(system, str(replace_file))

            try:
                with open(replace_file, "r") as f:
                    replace_lines = f.readlines()
                    LOG.debug(
                        f"Replace file contains {len(replace_lines)} instructions"
                    )
                    if len(replace_lines) > 0:
                        LOG.debug(f"First few instructions: {replace_lines[:5]}")
                    else:
                        LOG.warning("Replace file is empty!")
            except Exception as e:
                LOG.warning(f"Error reading replace file: {e}")

            # Run replacement with Chimera
            await self._run_replacement_with_chimera(
                system, config, str(replace_file), temp_dir
            )

            # Get final count for logging
            replaced_count = system.count_states(state="replace")
            LOG.info(
                f"{Fore.BLUE}Successfully replaced {replaced_count} crosslinks{Style.RESET_ALL}"
            )

            return system

        except Exception as e:
            LOG.error(f"Error in system-based replacement: {str(e)}")
            traceback.print_exc()
            raise GeometryGenerationError(
                message=f"Failed to replace crosslinks: {str(e)}",
                original_error=e,
                error_code="GEO_ERR_004",
                context={
                    "config": (
                        config.model_dump()
                        if hasattr(config, "model_dump")
                        else str(config)
                    )
                },
            )

    async def replace_direct(
        self, config: ColbuilderConfig, temp_dir: Path
    ) -> Tuple[Optional[System], Optional[Path]]:
        """
        Directly replace crosslinks in a PDB file without requiring a System object.

        This method provides a standalone replacement approach working directly with PDB files.
        It identifies complete crosslink pairs based on proximity and residue type,
        and ensures both residues in each pair are replaced together.

        Args:
            config: Configuration settings including input/output paths and replacement ratio
            temp_dir: Path to temporary directory for file operations

        Returns:
            Tuple[Optional[System], Optional[Path]]: (modified system, output PDB path)

        Raises:
            GeometryGenerationError: If replacement fails for any reason
        """
        try:
            input_pdb = Path(config.replace_file) if config.replace_file else None
            if not input_pdb or not input_pdb.exists():
                raise GeometryGenerationError(
                    message=f"Input PDB file not found: {input_pdb}",
                    error_code="GEO_ERR_004",
                )

            output_pdb = temp_dir / f"{config.output or 'output'}.pdb"

            # Create type-specific directory within the temp directory
            type_dir = temp_dir / "NC"
            type_dir.mkdir(parents=True, exist_ok=True)

            # Split PDB into models
            model_count = self._split_pdb_into_models(str(input_pdb), str(type_dir))

            # Create temporary crosslinks file within the type directory
            temp_crosslinks_file = type_dir / "temp_crosslinks.txt"
            LOG.debug(f"Creating temporary crosslinks file: {temp_crosslinks_file}")

            with open(temp_crosslinks_file, "w") as f:
                for model_id in range(model_count):
                    model_path = type_dir / f"{model_id}.caps.pdb"
                    if not model_path.exists():
                        LOG.warning(f"Model file not found: {model_path}")
                        continue

                    LOG.debug(f"Processing model file: {model_path}")
                    atom_count = 0
                    crosslink_count = 0

                    with open(model_path, "r") as model_file:
                        for line in model_file:
                            if line.startswith(("ATOM", "HETATM")):
                                atom_count += 1

                                if len(line) >= 20:
                                    resname = line[17:20].strip()

                                    if resname in [
                                        "L2Y",
                                        "L2X",
                                        "L3Y",
                                        "L3X",
                                        "LY2",
                                        "LY3",
                                        "LX2",
                                        "LX3",
                                        "LYX",
                                        "LXY",
                                        "LYY",
                                        "LXX",
                                        "L4Y",
                                        "L5Y",
                                        "LY4",
                                        "LY5",
                                        "LX4",
                                        "LX5",
                                        "LGX",
                                        "LPS",
                                        "AGS",
                                        "APD",
                                    ]:
                                        f.write(f"{model_id}|{line}")
                                        crosslink_count += 1

                    LOG.debug(
                        f"Model {model_id}: {atom_count} atoms total, {crosslink_count} crosslink atoms found"
                    )

            # Process crosslinks to identify residues
            LOG.info("Analyzing crosslinks...")
            crosslink_residues = {}

            with open(temp_crosslinks_file, "r") as f:
                for line in f:
                    parts = line.strip().split("|", 1)
                    if len(parts) != 2:
                        continue

                    model_id = int(parts[0])
                    atom_line = parts[1]

                    atom_name = atom_line[12:16].strip()
                    resname = atom_line[17:20].strip()
                    resid = atom_line[22:26].strip()
                    chain = atom_line[21]
                    x = float(atom_line[30:38])
                    y = float(atom_line[38:46])
                    z = float(atom_line[46:54])

                    key = (model_id, resname, resid, chain)

                    if key not in crosslink_residues:
                        crosslink_residues[key] = {
                            "model_id": model_id,
                            "resname": resname,
                            "resid": resid,
                            "chain": chain,
                            "atoms": [],
                            "position": None,
                        }

                    crosslink_residues[key]["atoms"].append(
                        {"atom_name": atom_name, "position": [x, y, z]}
                    )

            # Calculate center of each residue
            for key, residue in crosslink_residues.items():
                atoms = residue["atoms"]
                if not atoms:
                    continue

                x_sum = sum(atom["position"][0] for atom in atoms)
                y_sum = sum(atom["position"][1] for atom in atoms)
                z_sum = sum(atom["position"][2] for atom in atoms)

                residue["position"] = [
                    x_sum / len(atoms),
                    y_sum / len(atoms),
                    z_sum / len(atoms),
                ]

            # Filter out residues without positions
            crosslink_residues = {
                k: v for k, v in crosslink_residues.items() if v["position"] is not None
            }
            LOG.debug(f"Found {len(crosslink_residues)} crosslink residues")

            # Define crosslink pair types we're looking for
            pair_types = [
                (["L4Y", "L4X", "LY4", "LX4"], ["L5Y", "L5X", "LY5", "LX5"]),
                (["LGX", "LPS"], ["AGS", "APD"]),
            ]

            # Group residues by type
            residues_by_type = {}
            for key, residue in crosslink_residues.items():
                resname = residue["resname"]
                if resname not in residues_by_type:
                    residues_by_type[resname] = []
                residues_by_type[resname].append((key, residue))

            LOG.debug(f"Found residue types: {list(residues_by_type.keys())}")

            # Identify crosslink pairs
            all_pairs = []

            for type1_list, type2_list in pair_types:
                type1_residues = []
                for resname in type1_list:
                    if resname in residues_by_type:
                        type1_residues.extend(residues_by_type[resname])

                type2_residues = []
                for resname in type2_list:
                    if resname in residues_by_type:
                        type2_residues.extend(residues_by_type[resname])

                if not type1_residues or not type2_residues:
                    continue

                for i, (key1, residue1) in enumerate(type1_residues):
                    for j, (key2, residue2) in enumerate(type2_residues):
                        try:
                            pos1 = residue1["position"]
                            pos2 = residue2["position"]
                            dist = self._calculate_distance(pos1, pos2)

                            if dist <= 10.0:
                                pair = {
                                    "members": [residue1, residue2],
                                    "distance": dist,
                                    "keys": [key1, key2],
                                }
                                all_pairs.append(pair)
                                LOG.debug(
                                    f"Found pair at distance {dist:.2f}Å: "
                                    f"{residue1['resname']} {residue1['resid']}{residue1['chain']} (model {residue1['model_id']}) - "
                                    f"{residue2['resname']} {residue2['resid']}{residue2['chain']} (model {residue2['model_id']})"
                                )
                        except Exception as e:
                            LOG.warning(
                                f"Error calculating distance between {residue1['resname']} {residue1['resid']}{residue1['chain']} and "
                                f"{residue2['resname']} {residue2['resid']}{residue2['chain']}: {e}"
                            )

            # Sort pairs by distance and filter for uniqueness
            all_pairs.sort(key=lambda p: p["distance"])

            used_residues = set()
            unique_pairs = []

            for pair in all_pairs:
                key1 = pair["keys"][0]
                key2 = pair["keys"][1]

                if key1 not in used_residues and key2 not in used_residues:
                    unique_pairs.append(pair)
                    used_residues.add(key1)
                    used_residues.add(key2)

            # Calculate how many pairs to replace
            num_to_replace = min(
                math.ceil(len(unique_pairs) * config.ratio_replace / 100),
                len(unique_pairs),
            )

            LOG.info(
                f"Will replace {num_to_replace} pairs out of {len(unique_pairs)} (ratio: {config.ratio_replace}%)"
            )

            # Select random pairs for replacement
            pairs_to_replace = (
                random.sample(unique_pairs, num_to_replace)
                if num_to_replace > 0
                else []
            )

            instructions = []

            for pair in pairs_to_replace:
                for residue in pair["members"]:
                    model_id = residue["model_id"]
                    instructions.append(
                        f"{model_id}.caps.pdb LYS {residue['resid']} {residue['chain']}"
                    )
                    LOG.debug(
                        f"Will replace {residue['resname']} {residue['resid']}{residue['chain']} in model {model_id}"
                    )

            # Handle no replacements case
            if not instructions:
                LOG.warning("No replacement instructions generated")
                LOG.info(
                    "Creating output file without modifications since no replacements were made"
                )
                with open(output_pdb, "w") as out:
                    with open(input_pdb, "r") as in_file:
                        first_line = in_file.readline().strip()
                        if first_line.startswith("CRYST1"):
                            out.write(first_line + "\n")
                        else:
                            in_file.seek(0)
                            out.write(
                                "REMARK   Generated by colbuilder - no replacements made\n"
                            )

                        for line in in_file:
                            if line.startswith(("ATOM", "HETATM", "TER")):
                                out.write(line)

                    out.write("END\n")

                return None, output_pdb

            # Run Chimera to perform replacements
            success = self._run_chimera_replacement(
                instructions=instructions,
                chimera_scripts_dir=config.CHIMERA_SCRIPTS_DIR,
                system_type=str(type_dir),
                output_pdb=str(output_pdb),
                temp_dir=temp_dir,
            )

            if not success:
                LOG.error("Chimera replacement failed")
                raise GeometryGenerationError(
                    message="Chimera replacement process failed",
                    error_code="GEO_ERR_004",
                )

            # We now have individual model files with replacements
            # Create a System from the output for return
            try:
                crystal = Crystal(pdb=str(input_pdb))
                system = System(crystal=crystal)

                LOG.info(f"Final PDB with replacements written to: {output_pdb}")

                with open(output_pdb, "w") as out:
                    # Write header
                    with open(input_pdb, "r") as in_file:
                        first_line = in_file.readline().strip()
                        if first_line.startswith("CRYST1"):
                            out.write(first_line + "\n")
                        else:
                            out.write(
                                "REMARK   Generated by colbuilder direct replacement\n"
                            )

                    # Write all model contents in order
                    pdb_files = sorted(
                        [f for f in type_dir.glob("*.caps.pdb")],
                        key=lambda x: int(x.stem.split(".")[0]),
                    )

                    LOG.debug(
                        f"Found {len(pdb_files)} model files to include in output PDB"
                    )

                    for pdb_file in pdb_files:
                        with open(pdb_file, "r") as f:
                            for line in f:
                                if line.startswith(("ATOM", "HETATM", "TER")):
                                    if line.startswith("HETATM"):
                                        line = "ATOM  " + line[6:]
                                    out.write(line)

                    out.write("END\n")
                return system, output_pdb

            except Exception as e:
                LOG.error(f"Error creating output PDB with System class: {str(e)}")
                LOG.info("Falling back to direct file combining method")

                with open(output_pdb, "w") as out:
                    out.write("REMARK   Generated by colbuilder direct replacement\n")
                    try:
                        with open(input_pdb, "r") as in_file:
                            for line in in_file:
                                if line.startswith("CRYST1"):
                                    out.write(line)
                                    break
                    except Exception as e:
                        LOG.warning(f"Could not retrieve CRYST1 info: {e}")

                    pdb_files = sorted(
                        [f for f in type_dir.glob("*.caps.pdb")],
                        key=lambda x: int(x.stem.split(".")[0]),
                    )

                    for pdb_file in pdb_files:
                        with open(pdb_file, "r") as f:
                            for line in f:
                                if line.startswith(("ATOM", "HETATM", "TER")):
                                    if line.startswith("HETATM"):
                                        line = "ATOM  " + line[6:]
                                    out.write(line)

                    out.write("END\n")

                LOG.info(f"Created output PDB file using fallback method: {output_pdb}")
                return None, output_pdb

        except GeometryGenerationError:
            raise
        except Exception as e:
            LOG.error(f"Unexpected error in direct replacement: {str(e)}")
            traceback.print_exc()
            raise GeometryGenerationError(
                message=f"Unexpected error in direct replacement: {str(e)}",
                error_code="GEO_ERR_004",
            )

    def _calculate_fibril_bounds(
        self, system: Any, fibril_length: float
    ) -> Tuple[float, float]:
        """
        Calculate the z-coordinate bounds of the collagen fibril.

        Determines an appropriate z-coordinate range to focus crosslink replacement,
        based on model centers of geometry and the specified fibril length.

        Args:
            system: System containing models with spatial information
            fibril_length: Length of the fibril in nanometers

        Returns:
            Tuple of (z_min, z_max) as floating point values defining the bounds
        """
        try:
            z_values = []

            for model_id in system.get_models():
                model = system.get_model(model_id=model_id)

                try:
                    cog = model.get_cog()
                    z_center = self._get_z_position(cog)
                    z_values.append(z_center)
                except Exception as e:
                    LOG.warning(f"Could not get COG for model {model_id}: {e}")

            if not z_values:
                LOG.warning(
                    "No valid z-coordinates found in system models or crosslinks"
                )
                return (0.0, 10000.0)

            min_z = min(z_values) - 500.0
            max_z = max(z_values) + 500.0

            if max_z - min_z > 5000.0:
                model = system.get_model(model_id=0.0)
                cog = model.get_cog()
                z_center = self._get_z_position(cog)

                fibril_length_angstroms = fibril_length * 10.0

                z_min = z_center - fibril_length_angstroms * 5
                z_max = z_center + fibril_length_angstroms * 5
            else:
                z_min = min_z
                z_max = max_z

            return (z_min, z_max)
        except Exception as e:
            LOG.error(f"Error calculating fibril bounds: {e}")
            return (0.0, 10000.0)

    def _get_z_position(self, position: Any) -> float:
        """
        Extract the z-component from a position vector.

        Safely extracts the z-coordinate from various position representations
        (numpy arrays, lists, tuples) and converts it to a Python float.

        Args:
            position: Position object which may be a numpy array, list, tuple or similar

        Returns:
            Z-coordinate as a Python float

        Raises:
            IndexError, TypeError, ValueError: If extraction fails
        """
        try:
            if hasattr(position, "shape"):
                if len(position.shape) == 1 and position.shape[0] >= 3:
                    return float(position[2])
                elif len(position.shape) == 2 and position.shape[1] >= 3:
                    return float(position[0][2])

            if hasattr(position, "__getitem__"):
                z_val = position[2]
                if (
                    hasattr(z_val, "__getitem__")
                    and hasattr(z_val, "__len__")
                    and len(z_val) > 0
                ):
                    return float(z_val[0])
                return float(z_val)

            return float(position)

        except (IndexError, TypeError, ValueError) as e:
            LOG.error(
                f"Error extracting z position: {e}, position type: {type(position)}, value: {position}"
            )
            raise

    def _calculate_distance(self, pos1: Any, pos2: Any) -> float:
        """
        Calculate Euclidean distance between two positions.

        Handles various position representations (numpy arrays, lists, tuples)
        and computes the 3D Euclidean distance between them.

        Args:
            pos1: First position as a vector-like object
            pos2: Second position as a vector-like object

        Returns:
            Distance as a Python float

        Raises:
            IndexError, TypeError, ValueError: If position extraction fails
        """
        try:
            x1 = self._get_component(pos1, 0)
            y1 = self._get_component(pos1, 1)
            z1 = self._get_component(pos1, 2)

            x2 = self._get_component(pos2, 0)
            y2 = self._get_component(pos2, 1)
            z2 = self._get_component(pos2, 2)

            dx = x2 - x1
            dy = y2 - y1
            dz = z2 - z1

            return (dx * dx + dy * dy + dz * dz) ** 0.5

        except (IndexError, TypeError, ValueError) as e:
            LOG.error(f"Error calculating distance: {e}")
            LOG.error(f"pos1 type: {type(pos1)}, value: {pos1}")
            LOG.error(f"pos2 type: {type(pos2)}, value: {pos2}")
            raise

    def _get_component(self, position: Any, index: int) -> float:
        """
        Extract a specific component from a position vector.

        Safely extracts the specified component (x=0, y=1, z=2) from various
        position representations and converts it to a Python float.

        Args:
            position: Position object which may be a numpy array, list, tuple or similar
            index: Component index to extract (0=x, 1=y, 2=z)

        Returns:
            Component value as a Python float

        Raises:
            IndexError, TypeError, ValueError: If extraction fails
        """
        try:
            if hasattr(position, "shape"):
                if len(position.shape) == 1 and position.shape[0] > index:
                    return float(position[index])
                elif len(position.shape) == 2 and position.shape[1] > index:
                    return float(position[0][index])

            if hasattr(position, "__getitem__"):
                val = position[index]
                if (
                    hasattr(val, "__getitem__")
                    and hasattr(val, "__len__")
                    and len(val) > 0
                ):
                    return float(val[0])
                return float(val)

            return float(position)

        except (IndexError, TypeError, ValueError) as e:
            LOG.error(f"Error extracting component {index}: {e}")
            LOG.error(f"position type: {type(position)}, value: {position}")
            raise

    def _write_replace_file(self, system: Any, replace_file: str) -> None:
        """
        Write replacement instructions to a file for Chimera.

        Creates a file containing replacement instructions for crosslinks
        marked with state='replace' in the system. Each line specifies a model,
        residue name, residue ID, and chain for Chimera to process.

        Args:
            system: System with marked crosslinks
            replace_file: Path to the output file for instructions

        Raises:
            GeometryGenerationError: If the file cannot be written
        """
        try:
            LOG.debug(f"Writing replacement instructions to {replace_file}")
            with open(replace_file, "w") as f:
                instruction_count = 0
                for model_id in system.get_models():
                    model = system.get_model(model_id=model_id)
                    if not hasattr(model, "crosslink") or not model.crosslink:
                        continue

                    for cross in model.crosslink:
                        if hasattr(cross, "state") and cross.state == "replace":
                            if (
                                hasattr(cross, "resname")
                                and hasattr(cross, "resid")
                                and hasattr(cross, "chain")
                            ):
                                # Write instruction in correct format for Chimera script
                                instruction = f"{int(model_id)}.caps.pdb LYS {cross.resid} {cross.chain}\n"
                                f.write(instruction)
                                instruction_count += 1
                            else:
                                LOG.warning(
                                    f"Crosslink in model {model_id} missing required attributes (resname, resid, or chain)"
                                )

                if instruction_count == 0:
                    LOG.warning(
                        "No replacement instructions were written! Check if crosslinks have proper attributes."
                    )

        except Exception as e:
            LOG.error(f"Failed to write replacement file: {e}")
            raise GeometryGenerationError(
                message=f"Failed to write replacement file: {str(e)}",
                error_code="GEO_ERR_004",
            )

    async def _run_replacement_with_chimera(
        self,
        system: Any,
        config: ColbuilderConfig,
        replace_file: str,
        temp_dir: Optional[Path] = None,
    ) -> None:
        """
        Execute Chimera to perform crosslink replacements.

        Args:
            system: System with crosslinks to replace
            config: Configuration settings for Chimera execution
            replace_file: Path to file with replacement instructions
            temp_dir: Optional temporary directory for file operations
        """
        try:
            LOG.debug("Running Chimera to replace crosslinks with lysines")

            # Determine system type
            model_zero = system.get_model(model_id=0.0)
            if not hasattr(model_zero, "type"):
                LOG.warning("Model doesn't have a 'type' attribute, defaulting to 'D'")
                system_type = "D"  # Default to 'D' type, divalent
            else:
                system_type = model_zero.type

            if temp_dir is None:
                temp_dir = Path(config.working_directory) / ".tmp" / "replacement"
                temp_dir.mkdir(parents=True, exist_ok=True)
                LOG.debug(f"Created temporary directory for replacement: {temp_dir}")

            # Create type-specific directory within the temp directory
            type_dir = temp_dir / system_type
            type_dir.mkdir(parents=True, exist_ok=True)
            LOG.debug(f"Using type-specific directory for replacement: {type_dir}")

            copied_count = 0

            model_ids = set(int(float(model_id)) for model_id in system.get_models())
            LOG.debug(f"System has {len(model_ids)} model IDs")

            # Possible source locations for PDB files
            source_locations = [
                temp_dir,
                Path.cwd(),
                Path(config.working_directory),
                Path.cwd() / ".tmp" / "geometry_gen",
                temp_dir.parent if temp_dir else None,
            ]
            source_locations = [loc for loc in source_locations if loc and loc.exists()]

            for loc in source_locations:
                caps_files = list(loc.glob("*.caps.pdb"))

                for caps_file in caps_files:
                    try:
                        model_id = int(caps_file.stem.split(".")[0])

                        if model_id in model_ids:
                            target = type_dir / caps_file.name

                            if not target.exists():
                                shutil.copy2(caps_file, target)
                                copied_count += 1
                    except Exception as e:
                        LOG.error(f"Error processing caps file {caps_file}: {e}")

                if copied_count < len(model_ids):
                    regular_pdbs = list(loc.glob("*.pdb"))

                    for pdb_file in regular_pdbs:
                        if any(
                            x in pdb_file.name
                            for x in ["_geom_only", "_original", "output"]
                        ):
                            continue

                        try:
                            model_id = int(pdb_file.stem.split(".")[0])

                            if model_id in model_ids:
                                target = type_dir / f"{model_id}.caps.pdb"

                                if not target.exists():
                                    shutil.copy2(pdb_file, target)
                                    copied_count += 1
                        except Exception as e:
                            # Ignore errors for files that don't match our pattern
                            continue

            if copied_count < len(model_ids):
                missing_models = set(model_ids) - set(
                    int(caps_file.stem.split(".")[0])
                    for caps_file in type_dir.glob("*.caps.pdb")
                )
                LOG.debug(f"Missing models: {sorted(missing_models)}")

                for source_dir in source_locations:
                    for model_id in sorted(missing_models):
                        pattern = f"**/{model_id}*.pdb"
                        matching_files = list(source_dir.glob(pattern))

                        for found_file in matching_files:
                            if found_file.is_file():
                                target = type_dir / f"{model_id}.caps.pdb"
                                if not target.exists():
                                    try:
                                        shutil.copy2(found_file, target)
                                        copied_count += 1
                                        missing_models.remove(model_id)
                                        break
                                    except Exception as e:
                                        LOG.warning(
                                            f"Error copying file {found_file}: {e}"
                                        )

            caps_files = list(type_dir.glob("*.caps.pdb"))
            LOG.debug(
                f"Type directory contains {len(caps_files)} caps files for {len(model_ids)} models"
            )

            if len(caps_files) == 0:
                LOG.error(
                    "No caps files found in type directory! Creating placeholder files..."
                )

                for model_id in sorted(model_ids):
                    target = type_dir / f"{model_id}.caps.pdb"
                    if not target.exists():
                        try:
                            with open(target, "w") as f:
                                f.write("REMARK Created as placeholder for Chimera\n")
                                f.write(
                                    "ATOM      1  CA  GLY A   1      0.000   0.000   0.000  1.00 20.00\n"
                                )
                                f.write("TER\nEND\n")
                            LOG.debug(f"Created placeholder file for model {model_id}")
                        except Exception as e:
                            LOG.warning(
                                f"Failed to create placeholder for model {model_id}: {e}"
                            )

                caps_files = list(type_dir.glob("*.caps.pdb"))
                LOG.info(f"Created {len(caps_files)} placeholder model files")

            if len(caps_files) == 0:
                LOG.error(
                    "Could not find or create any model files! Replacement cannot proceed."
                )
                raise GeometryGenerationError(
                    message="No model files found for replacement",
                    error_code="GEO_ERR_004",
                )

            replace_path = Path(replace_file)
            LOG.debug(f"Using replacement instructions from: {replace_path}")

            if not replace_path.exists():
                LOG.error(f"Replacement file not found: {replace_path}")
                raise GeometryGenerationError(
                    message=f"Replacement file not found: {replace_path}",
                    error_code="GEO_ERR_004",
                )

            with open(replace_path, "r") as f:
                instructions = f.readlines()

            if not instructions:
                LOG.warning("Replacement file is empty, skipping Chimera execution")
                return

            # Find the swapaa.py script
            chimera_scripts_dir = None
            try:
                if (
                    hasattr(config, "CHIMERA_SCRIPTS_DIR")
                    and config.CHIMERA_SCRIPTS_DIR
                ):
                    chimera_scripts_dir = Path(config.CHIMERA_SCRIPTS_DIR)

                if not chimera_scripts_dir or not chimera_scripts_dir.exists():
                    if hasattr(self, "file_manager") and self.file_manager:
                        try:
                            chimera_scripts_dir = self.file_manager.find_file(
                                "chimera_scripts"
                            )
                        except:
                            pass

                    if not chimera_scripts_dir or not chimera_scripts_dir.exists():
                        potential_dirs = [
                            Path(config.working_directory) / "chimera_scripts",
                            Path(config.working_directory)
                            / "src"
                            / "colbuilder"
                            / "chimera_scripts",
                            Path.cwd() / "chimera_scripts",
                            Path.cwd() / "src" / "colbuilder" / "chimera_scripts",
                            Path("/usr/local/lib/colbuilder/chimera_scripts"),
                            Path(__file__).parent / "chimera_scripts",
                        ]

                        for dir_path in potential_dirs:
                            if dir_path.exists():
                                chimera_scripts_dir = dir_path
                                break

                if not chimera_scripts_dir or not chimera_scripts_dir.exists():
                    LOG.warning(
                        "Could not find chimera_scripts directory in standard locations"
                    )
                    search_dirs = [Path(config.working_directory), Path.cwd()]
                    for search_dir in search_dirs:
                        for script_dir in search_dir.glob("**/chimera_scripts"):
                            if script_dir.exists():
                                chimera_scripts_dir = script_dir
                                LOG.info(
                                    f"Found chimera_scripts directory through search: {chimera_scripts_dir}"
                                )
                                break
            except Exception as e:
                LOG.error(f"Error finding chimera_scripts directory: {e}")
                raise GeometryGenerationError(
                    message=f"Error finding chimera_scripts directory: {str(e)}",
                    error_code="GEO_ERR_004",
                )

            if not chimera_scripts_dir:
                LOG.error(
                    "Could not find chimera_scripts directory after extensive search"
                )
                raise GeometryGenerationError(
                    message="Chimera scripts directory not found",
                    error_code="GEO_ERR_004",
                )

            LOG.debug(f"Using chimera_scripts directory: {chimera_scripts_dir}")

            swapaa_script = chimera_scripts_dir / "swapaa.py"
            if not swapaa_script.exists():
                LOG.error(f"Chimera swapaa script not found at: {swapaa_script}")
                for pyfile in chimera_scripts_dir.glob("*.py"):
                    LOG.debug(f"Found script: {pyfile}")

                raise GeometryGenerationError(
                    message=f"Chimera swapaa script not found: {swapaa_script}",
                    error_code="GEO_ERR_004",
                )

            # Set up command for Chimera - make sure to pass proper arguments
            cmd = f'chimera --nogui --silent --script "{swapaa_script} {replace_path} {type_dir}"'

            import subprocess

            result = subprocess.run(
                cmd,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )

            if result.returncode != 0:
                LOG.error(f"Chimera failed with return code {result.returncode}")
                if result.stderr:
                    LOG.error(f"Chimera stderr: {result.stderr}")
                if result.stdout:
                    LOG.info(f"Chimera stdout: {result.stdout}")

                raise GeometryGenerationError(
                    message="Chimera replacement failed", error_code="GEO_ERR_004"
                )

            if result.stdout:
                LOG.debug(f"Chimera output: {result.stdout}")

            # After Chimera runs, verify the files were modified
            modified_caps_files = list(type_dir.glob("*.caps.pdb"))

            # Check if the files have LYS residues that weren't there before
            has_replacements = False
            for caps_file in modified_caps_files[:5]:  # Check a few files
                try:
                    with open(caps_file, "r") as f:
                        content = f.read()
                        if "LYS" in content:
                            has_replacements = True
                            LOG.debug(
                                f"Final check: Found LYS residue in {caps_file}, confirming replacement"
                            )
                            break
                except Exception as e:
                    LOG.warning(f"Error checking {caps_file} for LYS residues: {e}")

            if not has_replacements:
                LOG.warning(
                    "No LYS residues found in sampled caps files. Replacement might not have worked."
                )

        except Exception as e:
            LOG.error(f"Error during replacement with Chimera: {str(e)}")
            import traceback

            LOG.error(f"Traceback: {traceback.format_exc()}")
            raise GeometryGenerationError(
                message=f"Error during replacement with Chimera: {str(e)}",
                error_code="GEO_ERR_004",
            )

    def _read_replacement_instructions(self, replace_file: str) -> List[str]:
        """
        Read replacement instructions from a file.

        Parses the replacement instructions file and returns a list of instruction strings.

        Args:
            replace_file: Path to the replacement instructions file

        Returns:
            List of replacement instruction strings

        Raises:
            GeometryGenerationError: If the file cannot be read
        """
        instructions = []
        try:
            with open(replace_file, "r") as f:
                instructions = [line.strip() for line in f if line.strip()]
        except Exception as e:
            LOG.error(f"Error reading replacement file: {e}")
            raise GeometryGenerationError(
                message=f"Failed to read replacement file: {str(e)}",
                error_code="GEO_ERR_004",
            )
        return instructions

    def _run_chimera_replacement(
        self,
        instructions: List[str],
        chimera_scripts_dir: str,
        system_type: str,
        output_pdb: Optional[str] = None,
        temp_dir: Optional[Path] = None,
    ) -> bool:
        """
        Run Chimera to perform replacements on PDB files.

        Executes Chimera with the swapaa.py script to replace residues in the
        specified PDB files according to the provided instructions.

        Args:
            instructions: List of replacement instructions
            chimera_scripts_dir: Directory containing Chimera scripts
            system_type: Type of system (directory name)
            output_pdb: Optional path for output PDB file
            temp_dir: Optional temporary directory for intermediate files

        Returns:
            True if successful, False if replacement failed
        """
        try:
            scripts_dir = Path(chimera_scripts_dir)
            system_type_path = Path(system_type)

            if temp_dir and temp_dir.exists():
                replace_file = temp_dir / "replace.txt"
            else:
                replace_file = Path("replace.txt")

            with open(replace_file, "w") as f:
                f.write("\n".join(instructions))

            LOG.debug(f"Created replacement file: {replace_file}")

            LOG.debug("Running Chimera to replace crosslinks")
            swapaa_script = scripts_dir / "swapaa.py"

            if not swapaa_script.exists():
                raise GeometryGenerationError(
                    message=f"Chimera swapaa script not found: {swapaa_script}",
                    error_code="GEO_ERR_004",
                )

            base_file_str = str(replace_file)
            if base_file_str.endswith(".txt"):
                base_file_str = base_file_str[:-4]

            cmd = f'chimera --nogui --silent --script "{swapaa_script} {base_file_str} {system_type_path}"'

            LOG.debug(f"        Running command: {cmd}")

            result = subprocess.run(
                cmd,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )

            LOG.debug(f"    Chimera returned: {result.returncode}")

            if result.returncode != 0:
                LOG.error(f"Chimera failed with return code {result.returncode}")
                if result.stderr:
                    LOG.error(f"Chimera stderr: {result.stderr}")
                if result.stdout:
                    LOG.info(f"Chimera stdout: {result.stdout}")

                return False

            if output_pdb:
                output_path = Path(output_pdb)
                with open(output_path, "w") as out:
                    out.write("REMARK   Generated by colbuilder replacement\n")

                    pdb_files = sorted(
                        [f for f in system_type_path.glob("*.caps.pdb")],
                        key=lambda x: int(x.stem.split(".")[0]),
                    )

                    LOG.debug(f"Adding {len(pdb_files)} model files to output PDB")

                    for pdb_file in pdb_files:
                        with open(pdb_file, "r") as f:
                            for line in f:
                                if line.startswith(("ATOM", "HETATM", "TER")):
                                    if line.startswith("HETATM"):
                                        line = "ATOM  " + line[6:]
                                    out.write(line)

                    out.write("END\n")

            LOG.debug("Chimera replacement completed successfully")
            return True

        except Exception as e:
            LOG.error(f"Error during Chimera execution: {str(e)}")
            traceback.print_exc()
            return False

    def _split_pdb_into_models(
        self, input_pdb: Union[str, Path], system_type: Union[str, Path]
    ) -> int:
        """
        Split a PDB file into individual triple helix models.

        Divides the PDB file into individual models using the repeating A-B-C
        chain pattern that is typical of collagen triple helices.

        Args:
            input_pdb: Path to input PDB file
            system_type: Type of system/directory to create

        Returns:
            Number of models created

        Raises:
            GeometryGenerationError: If the input PDB file is not a valid Colbuilder structure
        """
        input_pdb_path = Path(input_pdb)
        system_type_path = Path(system_type)

        is_colbuilder_pdb = False
        with open(input_pdb_path, "r") as f:
            first_line = f.readline().strip()
            if first_line.startswith("CRYST1"):
                is_colbuilder_pdb = True

        if not is_colbuilder_pdb:
            raise GeometryGenerationError(
                message="Input PDB file does not appear to be a valid Colbuilder-generated structure. "
                "Direct replacement only works on PDB files generated by the Colbuilder geometry module.",
                error_code="GEO_ERR_004",
            )

        with open(input_pdb_path, "r") as f:
            all_lines = f.readlines()

        cryst_line = None
        for line in all_lines:
            if line.startswith("CRYST1"):
                cryst_line = line
                break

        atom_lines = [
            line for line in all_lines if line.startswith(("ATOM", "HETATM", "TER"))
        ]

        if not atom_lines:
            raise GeometryGenerationError(
                message="Input PDB file does not contain any ATOM, HETATM, or TER records.",
                error_code="GEO_ERR_004",
            )

        models = []
        current_model = []
        chain_sequence = []

        for line in atom_lines:
            if line.startswith(("ATOM", "HETATM")):
                chain = line[21]

                if not chain_sequence or chain != chain_sequence[-1]:
                    chain_sequence.append(chain)

                    if len(chain_sequence) > 3 and chain_sequence[-4:] == [
                        "A",
                        "B",
                        "C",
                        "A",
                    ]:
                        models.append(current_model)
                        current_model = []
                        chain_sequence = ["A"]

            current_model.append(line)

        if current_model:
            models.append(current_model)

        if not models:
            LOG.warning(
                "Could not split into multiple models. Using entire PDB as one model."
            )
            models = [atom_lines]

        LOG.debug(f"Split PDB into {len(models)} models")

        system_type_path.mkdir(parents=True, exist_ok=True)

        for i, model_lines in enumerate(models):
            caps_file = system_type_path / f"{i}.caps.pdb"
            with open(caps_file, "w") as f:
                if cryst_line:
                    f.write(cryst_line)

                f.writelines(model_lines)

                if not model_lines[-1].startswith("TER"):
                    f.write("TER\n")

            LOG.debug(f"Created model file: {caps_file}")

        return len(models)

    def _select_replacements(
        self, system: Any, ratio_replace: float, z_bounds: Tuple[float, float]
    ) -> None:
        """
        Select crosslinks to replace based on the given ratio and spatial constraints.

        Identifies complete crosslink pairs within the specified z-bounds, ensures
        they are properly paired, and marks a subset for replacement based on the
        specified ratio.

        Args:
            system: System containing models and crosslinks
            ratio_replace: Percentage of crosslinks to replace (0-100)
            z_bounds: Tuple of (z_min, z_max) defining the region to focus on
        """
        reset_count = 0
        for model_id in system.get_models():
            model = system.get_model(model_id=model_id)
            if hasattr(model, "crosslink") and model.crosslink:
                for crosslink in model.crosslink:
                    if hasattr(crosslink, "state"):
                        crosslink.state = "none"
                        reset_count += 1

        LOG.debug(f"Reset state for {reset_count} crosslinks")

        crosslinks_in_bounds = []
        for model_id in system.get_models():
            model = system.get_model(model_id=model_id)
            if not hasattr(model, "crosslink") or not model.crosslink:
                continue

            for crosslink in model.crosslink:
                try:
                    if not hasattr(crosslink, "position"):
                        LOG.debug(
                            f"Crosslink in model {model_id} has no position attribute"
                        )
                        continue

                    z_pos = self._get_z_position(crosslink.position)

                    if z_bounds[0] <= z_pos <= z_bounds[1]:
                        crosslink_info = {
                            "model_id": model_id,
                            "model": model,
                            "crosslink": crosslink,
                            "type": getattr(crosslink, "type", "D"),
                            "resname": getattr(crosslink, "resname", "UNK"),
                            "z_position": z_pos,
                        }
                        crosslinks_in_bounds.append(crosslink_info)

                except Exception as e:
                    LOG.warning(
                        f"Error checking crosslink bounds for model {model_id}: {e}"
                    )

        if not crosslinks_in_bounds:
            LOG.warning("No crosslinks found within z-bounds, nothing will be replaced")
            return

        crosslinks_by_type = {}
        for cross_ref in crosslinks_in_bounds:
            crosslink = cross_ref["crosslink"]
            resname = getattr(crosslink, "resname", "UNK")
            if resname not in crosslinks_by_type:
                crosslinks_by_type[resname] = []
            crosslinks_by_type[resname].append(cross_ref)

        LOG.debug(f"Found crosslink types: {list(crosslinks_by_type.keys())}")

        # Define which residue types should be paired together
        pair_types = [
            (["L4Y", "L4X", "LY4", "LX4"], ["L5Y", "L5X", "LY5", "LX5"]),
            (["LGX", "LPS"], ["AGS", "APD"]),
        ]

        pairs = []

        for type1_list, type2_list in pair_types:
            type1_crosslinks = []
            for resname in type1_list:
                if resname in crosslinks_by_type:
                    type1_crosslinks.extend(crosslinks_by_type[resname])

            type2_crosslinks = []
            for resname in type2_list:
                if resname in crosslinks_by_type:
                    type2_crosslinks.extend(crosslinks_by_type[resname])

            if not type1_crosslinks or not type2_crosslinks:
                continue

            used_type1 = set()
            used_type2 = set()

            for i, cross_ref1 in enumerate(type1_crosslinks):
                if i in used_type1:
                    continue

                crosslink1 = cross_ref1["crosslink"]
                min_dist = float("inf")
                closest_idx = -1

                for j, cross_ref2 in enumerate(type2_crosslinks):
                    if j in used_type2:
                        continue

                    crosslink2 = cross_ref2["crosslink"]

                    try:
                        dist = self._calculate_distance(
                            crosslink1.position, crosslink2.position
                        )
                        if dist < min_dist and dist <= 5.0:
                            min_dist = dist
                            closest_idx = j
                    except Exception as e:
                        LOG.warning(f"Error calculating distance: {e}")

                if closest_idx != -1:
                    pairs.append([cross_ref1, type2_crosslinks[closest_idx]])
                    used_type1.add(i)
                    used_type2.add(closest_idx)

                    if LOG.level <= 10:  # DEBUG level
                        crosslink1_info = (
                            f"{crosslink1.resname} {crosslink1.resid}{crosslink1.chain}"
                            if all(
                                hasattr(crosslink1, attr)
                                for attr in ["resname", "resid", "chain"]
                            )
                            else "Unknown"
                        )
                        crosslink2_info = (
                            f"{type2_crosslinks[closest_idx]['crosslink'].resname} {type2_crosslinks[closest_idx]['crosslink'].resid}{type2_crosslinks[closest_idx]['crosslink'].chain}"
                            if all(
                                hasattr(
                                    type2_crosslinks[closest_idx]["crosslink"], attr
                                )
                                for attr in ["resname", "resid", "chain"]
                            )
                            else "Unknown"
                        )

                        LOG.debug(
                            f"Found pair: {crosslink1_info} + {crosslink2_info} (distance: {min_dist:.2f}Å)"
                        )

        num_to_replace = min(
            max(1, math.ceil(len(pairs) * ratio_replace / 100)), len(pairs)
        )

        LOG.info(
            f"Will replace {num_to_replace} pairs out of {len(pairs)} (ratio: {ratio_replace}%)"
        )

        import random

        random.seed(int(time.time()))  # Set random seed for reproducibility
        random.shuffle(pairs)
        pairs_to_replace = pairs[:num_to_replace]

        replaced_count = 0
        for pair in pairs_to_replace:
            for cross_ref in pair:
                cross_ref["crosslink"].state = "replace"
                replaced_count += 1

                model_id = cross_ref["model_id"]
                crosslink = cross_ref["crosslink"]

                if (
                    hasattr(crosslink, "resname")
                    and hasattr(crosslink, "resid")
                    and hasattr(crosslink, "chain")
                ):
                    LOG.debug(
                        f"Marked {crosslink.resname} {crosslink.resid}{crosslink.chain} in model {model_id} for replacement"
                    )

        # Protect nearby crosslinks from being replaced
        protected_count = 0
        for pair in pairs_to_replace:
            for cross_ref in pair:
                crosslink = cross_ref["crosslink"]
                model = cross_ref["model"]

                for other in model.crosslink:
                    if other.state != "none" or other in [p["crosslink"] for p in pair]:
                        continue

                    try:
                        distance = self._calculate_distance(
                            crosslink.position, other.position
                        )
                        if 5.0 < distance <= 15.0:
                            other.state = "protect"
                            protected_count += 1
                    except Exception as e:
                        LOG.warning(f"Error protecting nearby crosslink: {e}")


# Backward compatibility functions


async def replace_in_system(system: Any, config: ColbuilderConfig) -> Any:
    """
    Backward compatibility function for system-based crosslink replacement.

    This is a standalone function that creates a CrosslinkReplacer and calls
    its replace_in_system method, maintaining backward compatibility with
    code that uses the function-based API.

    Args:
        system: System containing the models to be modified
        config: Configuration for replacement

    Returns:
        Modified system with replaced crosslinks
    """
    replacer = CrosslinkReplacer()
    return await replacer.replace_in_system(system, config)


async def direct_replace_geometry(config: ColbuilderConfig) -> None:
    """
    Backward compatibility function for direct PDB-based crosslink replacement.

    This is a standalone function that creates a CrosslinkReplacer and calls
    its replace_direct method, maintaining backward compatibility with
    code that uses the function-based API.

    Args:
        config: Configuration settings for the replacement
    """
    replacer = CrosslinkReplacer()
    temp_dir = Path(config.working_directory) / ".tmp" / "replacement_direct"
    temp_dir.mkdir(parents=True, exist_ok=True)
    await replacer.replace_direct(config, temp_dir)
