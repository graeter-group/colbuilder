"""
Colbuilder Main Geometry Module

Main entry point for geometry operations, coordinating the crystal,
mixing, and replacement services.
"""

import os
import time
import shutil
import traceback
import asyncio
import logging
import copy
from pathlib import Path
from typing import Optional, Set, List, Union, Tuple, Dict

from colbuilder.core.utils.exceptions import GeometryGenerationError
from colbuilder.core.utils.config import ColbuilderConfig
from colbuilder.core.utils.logger import setup_logger
from colbuilder.core.utils.files import FileManager, managed_resources
from .crystal_builder import CrystalBuilder
from .crosslink_mixer import CrosslinkMixer
from .geometry_replacer import CrosslinkReplacer
from .system import System
from .crystal import Crystal
from colorama import init, Fore, Style

LOG = setup_logger(__name__)

STANDARD_TEMP_FILES = {"replace.txt"}
STANDARD_TEMP_DIRS = {}


def cleanup_temp_files(
    temp_files: Optional[Set[str]] = None,
    temp_dirs: Optional[Set[str]] = None,
    include_standard: bool = True,
) -> None:
    """
    Clean up temporary files and directories.

    Args:
        temp_files: Optional set of specific temporary files to clean up
        temp_dirs: Optional set of specific temporary directories to clean up
        include_standard: Whether to also clean up standard temporary files/directories
    """
    files_to_clean = set(temp_files) if temp_files else set()
    dirs_to_clean = set(temp_dirs) if temp_dirs else set()

    if include_standard:
        files_to_clean.update(STANDARD_TEMP_FILES)
        dirs_to_clean.update(STANDARD_TEMP_DIRS)

    for dir_path in dirs_to_clean:
        if os.path.exists(dir_path) and os.path.isdir(dir_path):
            try:
                shutil.rmtree(dir_path)
                LOG.info(f"Removed temporary directory: {dir_path}")
            except Exception as e:
                LOG.warning(
                    f"Failed to remove temporary directory {dir_path}: {str(e)}"
                )

    for file_path in files_to_clean:
        if os.path.exists(file_path) and os.path.isfile(file_path):
            try:
                os.remove(file_path)
                LOG.info(f"Removed temporary file: {file_path}")
            except Exception as e:
                LOG.warning(f"Failed to remove temporary file {file_path}: {str(e)}")


class GeometryService:
    """
    Main service coordinating geometry operations.

    This class orchestrates the different components involved in collagen
    geometry generation, mixing, and crosslink replacement.
    """

    def __init__(
        self, config: ColbuilderConfig, file_manager: Optional[FileManager] = None
    ) -> None:
        """
        Initialize geometry service with configuration.

        Args:
            config: Configuration settings for geometry operations
            file_manager: Optional file manager for consistent file handling
        """
        self.config = config
        self.crystal_service = CrystalBuilder()
        self.mixer_service = CrosslinkMixer()
        self.replacer_service = CrosslinkReplacer()
        self.temp_files: Set[str] = set()
        self.temp_dirs: Set[str] = set()
        self.file_manager = file_manager or FileManager(config)
        self.original_dir = Path.cwd()
        self.temp_dir = None

        self.crystal_service.set_file_manager(self.file_manager)

    def _track_temp_resources(
        self, files: Optional[List[str]] = None, dirs: Optional[List[str]] = None
    ) -> None:
        """
        Track temporary files and directories for later cleanup.

        Args:
            files: List of file paths to track
            dirs: List of directory paths to track
        """
        if files:
            self.temp_files.update(files)
            for file_path in files:
                LOG.info(f"Tracking temporary file: {file_path}")

        if dirs:
            self.temp_dirs.update(dirs)
            for dir_path in dirs:
                LOG.info(f"Tracking temporary directory: {dir_path}")

    def _cleanup(self) -> None:
        """Clean up all tracked temporary resources."""
        LOG.info("Cleaning up temporary files...")
        cleanup_temp_files(self.temp_files, self.temp_dirs)

    async def _handle_mixing_only(self) -> Tuple[Optional[System], Optional[Path]]:
        """
        Handle mixing-only operation without geometry generation.

        This method performs crosslink mixing.
        The process includes:
        1. Validating mixing ratios
        2. Copying source files to a dedicated mixing directory
        3. Creating a basic system from the source files
        4. Executing the mixing service to blend different crosslink types
        5. Writing the final mixed system to a PDB file

        Returns:
            Tuple[Optional[System], Optional[Path]]: The mixed system object and output PDB file path

        Raises:
            GeometryGenerationError: If mixing fails, input validation fails, or output generation fails
        """
        if not self.config.files_mix:
            raise GeometryGenerationError(
                message="No files_mix provided for mixing operation",
                error_code="GEO_ERR_003",
            )

        original_dir: Path = Path.cwd()
        original_files_mix = self.config.files_mix

        try:
            mixing_dir: Path = self.file_manager.ensure_mixing_dir()
            LOG.debug(
                f"{Fore.BLUE}Using mixing directory: {mixing_dir}{Style.RESET_ALL}"
            )

            # Process ratio_mix format
            if isinstance(self.config.ratio_mix, str):
                ratio_dict: Dict[str, int] = {}
                for part in self.config.ratio_mix.split():
                    if ":" in part:
                        key, value = part.split(":")
                        try:
                            ratio_dict[key] = int(value)
                        except ValueError:
                            LOG.error(f"Invalid ratio value in {part}")
                            ratio_dict[key] = 0
                self.config.ratio_mix = ratio_dict
                LOG.info(f"Converted ratio_mix: {self.config.ratio_mix}")

            # Validate ratio_mix
            if not isinstance(self.config.ratio_mix, dict):
                raise GeometryGenerationError(
                    message=f"Invalid ratio_mix format: {self.config.ratio_mix}",
                    error_code="GEO_ERR_009",
                )

            # Copy mix files to mixing directory
            updated_mix_files: List[Path] = []
            for mix_file in self.config.files_mix:
                dest_path = self.file_manager.copy_to_directory(
                    mix_file, dest_dir=mixing_dir
                )
                updated_mix_files.append(dest_path)

            # Validate file count
            if len(updated_mix_files) < len(self.config.ratio_mix):
                LOG.warning(
                    f"Not enough files ({len(updated_mix_files)}) for ratio_mix types ({len(self.config.ratio_mix)})"
                )

            # Update config with local paths
            self.config.files_mix = updated_mix_files

            # Create base system from first input file
            system: System = System()
            crystal: Crystal = Crystal(pdb=str(updated_mix_files[0]))
            system.crystal = crystal

            # Perform mixing operation
            system, output_pdb = await self.mixer_service.mix(
                system, self.config, mixing_dir
            )

            # Verify the system has models
            if system:
                model_count = len(list(system.get_models()))
                LOG.info(f"Mixed system contains {model_count} models")
                
                type_counts = {}
                for model_id in system.get_models():
                    model = system.get_model(model_id=model_id)
                    if hasattr(model, 'type'):
                        model_type = model.type
                        type_counts[model_type] = type_counts.get(model_type, 0) + 1
                LOG.info(f"Type distribution: {type_counts}")
            else:
                LOG.warning("System object is None after mixing")

            # Copy to final output location
            if output_pdb and output_pdb.exists():
                final_pdb: Path = self.file_manager.get_output_path(
                    self.config.output, ".pdb"
                )

                self.file_manager.copy_to_output(output_pdb, dest_name=final_pdb.name)
                LOG.info(f"Final mixed PDB copied to: {final_pdb}")
                
                # Return both the system object and the PDB path
                return system, final_pdb
            else:
                raise GeometryGenerationError(
                    message="Mixing operation did not produce an output file",
                    error_code="GEO_ERR_004",
                )

        except GeometryGenerationError:
            raise
        except Exception as e:
            error_msg = f"Mixing operation failed: {str(e)}"
            LOG.error(error_msg)
            LOG.debug(f"Traceback: {traceback.format_exc()}")
            raise GeometryGenerationError(
                message=error_msg, error_code="GEO_ERR_005", original_error=e
            )

        finally:
            self.config.files_mix = original_files_mix
            os.chdir(original_dir)

    async def _handle_direct_replacement(
        self,
    ) -> Tuple[Optional[System], Optional[Path]]:
        """
        Handle direct replacement operation without geometry generation.

        This method performs crosslink replacement directly on an input PDB file when
        no geometry generation is required. The process includes:
        1. Validating the replacement file
        2. Copying the file to a dedicated replacement directory
        3. Running the replacement service to modify crosslinks
        4. Copying the resulting PDB to the output location

        Returns:
            Tuple[Optional[System], Optional[Path]]: Modified system (which might be None) and output PDB path

        Raises:
            GeometryGenerationError: If replacement file is missing or invalid, or if replacement fails
        """
        if not self.config.replace_file:
            raise GeometryGenerationError(
                message="No replace_file provided for direct replacement operation",
                error_code="GEO_ERR_004",
            )

        original_dir: Path = Path.cwd()

        try:
            replacement_dir: Path = self.file_manager.ensure_replacement_dir()
            LOG.debug(f"Using replacement directory: {replacement_dir}")

            replace_path: Path = Path(self.config.replace_file)
            if not replace_path.exists():
                raise GeometryGenerationError(
                    message=f"Replacement file not found: {replace_path}",
                    error_code="GEO_ERR_004",
                )

            local_replace: Path = replacement_dir / replace_path.name
            shutil.copy2(replace_path, local_replace)
            LOG.debug(f"Copied replacement file to: {local_replace}")

            temp_config: ColbuilderConfig = self.config.copy()
            temp_config.replace_file = local_replace
            temp_config.working_directory = replacement_dir

            self.replacer_service.file_manager = self.file_manager

            LOG.info(f"Performing direct replacement using file: {replace_path.name}")
            system: Optional[System]
            output_pdb: Optional[Path]
            system, output_pdb = await self.replacer_service.replace_direct(
                temp_config, replacement_dir
            )

            if output_pdb and output_pdb.exists():
                final_pdb: Path = self.file_manager.get_output_path(
                    self.config.output, ".pdb"
                )

                self.file_manager.copy_to_output(output_pdb, dest_name=final_pdb.name)
                return system, final_pdb
            else:
                raise GeometryGenerationError(
                    message="Replacement operation did not produce an output file",
                    error_code="GEO_ERR_004",
                )

        except GeometryGenerationError:
            raise
        except Exception as e:
            error_msg = f"Replacement operation failed: {str(e)}"
            LOG.error(error_msg)
            LOG.debug(f"Traceback: {traceback.format_exc()}")
            raise GeometryGenerationError(
                message=error_msg, error_code="GEO_ERR_004", original_error=e
            )

        finally:
            os.chdir(original_dir)

    async def _handle_full_generation(self) -> Tuple[Optional[System], Optional[Path]]:
        """
        Handle full geometry generation with optional mixing and replacement.
        Avoids redundant PDB writing when possible.

        Returns:
            Tuple[Optional[System], Optional[Path]]: Generated system and output PDB path

        Raises:
            GeometryGenerationError: If any step fails
        """
        original_dir = Path.cwd()

        try:
            geometry_dir = self.file_manager.ensure_geometry_dir()
            LOG.debug(f"Using geometry directory for generation: {geometry_dir}")

            self.temp_dir = geometry_dir

            os.chdir(geometry_dir)

            temp_config = self.config.copy()
            if self.config.pdb_file:
                pdb_path = Path(self.config.pdb_file)
                if pdb_path.exists():
                    local_pdb = geometry_dir / pdb_path.name
                    shutil.copy2(pdb_path, local_pdb)
                    temp_config.pdb_file = local_pdb
                else:
                    raise GeometryGenerationError(
                        message=f"PDB file not found: {pdb_path}",
                        error_code="GEO_ERR_001",
                    )

            temp_config.working_directory = geometry_dir

            system = None
            output_pdb_path = None
            final_temp_dir = geometry_dir

            # Perform geometry generation
            if temp_config.geometry_generator:
                LOG.info(
                    f"{Fore.MAGENTA}Geometry mode: building system from crystal{Style.RESET_ALL}"
                )
                self.crystal_service.set_file_manager(self.file_manager)
                system = await self.crystal_service.build(temp_config)

                geometry_only_system = copy.deepcopy(system)

                if not (temp_config.mix_bool or temp_config.replace_bool):
                    output_prefix = temp_config.output or temp_config.species
                    output_pdb_path = geometry_dir / f"{output_prefix}.pdb"

                    import logging

                    original_level = LOG.level
                    try:
                        LOG.setLevel(logging.ERROR)  # Suppress non-error messages
                        system.write_pdb(
                            pdb_out=output_pdb_path,
                            fibril_length=temp_config.fibril_length,
                            cleanup=False,
                            temp_dir=geometry_dir,
                        )
                    finally:
                        LOG.setLevel(original_level)

                    final_pdb = self.file_manager.copy_to_output(output_pdb_path)
                    return system, final_pdb

            # Perform mixing if required
            if temp_config.mix_bool and system:
                LOG.section("Mixing geometry...")
                mixing_dir = self.file_manager.ensure_mixing_dir()
                system = await self.mixer_service.mix(system, temp_config, mixing_dir)
                LOG.info("Mixing completed.")

            # Perform replacement if required
            if temp_config.replace_bool and system:
                LOG.section("Replacing crosslinks...")

                model_type = (
                    system.get_model(model_id=0.0).type
                    if hasattr(system.get_model(model_id=0.0), "type")
                    else "D"
                )

                geometry_type_dir = geometry_dir / model_type
                if geometry_type_dir.exists():
                    pre_replacement_counts = {}
                    atom_count = 0
                    for caps_file in geometry_type_dir.glob("*.caps.pdb"):
                        with open(caps_file, "r") as f:
                            for line in f:
                                if (
                                    line.startswith(("ATOM", "HETATM"))
                                    and len(line) >= 20
                                ):
                                    atom_count += 1
                                    resname = line[17:20].strip()
                                    pre_replacement_counts[resname] = (
                                        pre_replacement_counts.get(resname, 0) + 1
                                    )

                replacement_dir = self.file_manager.ensure_replacement_dir()
                LOG.debug(f"Using replacement directory: {replacement_dir}")

                self.replacer_service.file_manager = self.file_manager

                system = await self.replacer_service.replace_in_system(
                    system, temp_config, replacement_dir
                )

                replacement_type_dir = replacement_dir / model_type

                if replacement_type_dir.exists():
                    caps_files = list(replacement_type_dir.glob("*.caps.pdb"))

                    if caps_files:
                        output_prefix = temp_config.output or temp_config.species
                        output_pdb_path = replacement_dir / f"{output_prefix}.pdb"
                        final_temp_dir = replacement_dir

                        system.write_pdb(
                            pdb_out=output_pdb_path,
                            fibril_length=temp_config.fibril_length,
                            cleanup=False,
                            temp_dir=replacement_dir,
                        )

                        if output_pdb_path.exists():
                            try:
                                post_replacement_counts = {}
                                with open(output_pdb_path, "r") as f:
                                    for line in f:
                                        if (
                                            line.startswith(("ATOM", "HETATM"))
                                            and len(line) >= 20
                                        ):
                                            resname = line[17:20].strip()
                                            post_replacement_counts[resname] = (
                                                post_replacement_counts.get(resname, 0)
                                                + 1
                                            )

                            except Exception as e:
                                LOG.warning(f"Error analyzing output PDB: {e}")

                            final_pdb = self.file_manager.copy_to_output(
                                output_pdb_path
                            )
                            LOG.info(
                                f"Final PDB with replacements written to: {final_pdb}"
                            )
                            return system, final_pdb
                    else:
                        LOG.warning(
                            "No caps files found in replacement directory. Falling back to geometry directory."
                        )
                else:
                    LOG.warning(
                        f"Replacement type directory not found: {replacement_type_dir}"
                    )

                # If we get here, replacement was requested but replacement directory didn't have caps
                # Write to the geometry directory instead
                output_prefix = temp_config.output or temp_config.species
                output_pdb_path = geometry_dir / f"{output_prefix}.pdb"

                LOG.info(
                    f"Writing system with replacements to geometry directory: {output_pdb_path}"
                )
                system.write_pdb(
                    pdb_out=output_pdb_path,
                    fibril_length=temp_config.fibril_length,
                    cleanup=False,
                    temp_dir=geometry_dir,
                )

                if output_pdb_path.exists():
                    final_pdb = self.file_manager.copy_to_output(output_pdb_path)
                    LOG.info(f"Final PDB with replacements copied to: {final_pdb}")
                    return system, final_pdb

            # Write the final system PDB if we've done mixing (since replacement handles its own writing)
            if system and temp_config.mix_bool:
                output_prefix = temp_config.output or temp_config.species
                output_pdb_path = final_temp_dir / f"{output_prefix}.pdb"

                LOG.info(f"Writing final mixed system PDB to {output_pdb_path}")
                system.write_pdb(
                    pdb_out=output_prefix,
                    fibril_length=temp_config.fibril_length,
                    cleanup=False,
                    temp_dir=final_temp_dir,
                )

                if output_pdb_path.exists():
                    final_pdb = self.file_manager.copy_to_output(output_pdb_path)
                    LOG.info(f"Final mixed PDB copied to: {final_pdb}")
                    return system, final_pdb
                else:
                    LOG.error(f"Output PDB not found at: {output_pdb_path}")
                    raise GeometryGenerationError(
                        message="Failed to create output PDB file",
                        error_code="GEO_ERR_012",
                    )

            if system:
                output_prefix = temp_config.output or temp_config.species
                output_pdb_path = geometry_dir / f"{output_prefix}.pdb"

                LOG.warning(
                    f"No PDB was written earlier, writing fallback PDB to {output_pdb_path}"
                )
                system.write_pdb(
                    pdb_out=output_pdb_path,
                    fibril_length=temp_config.fibril_length,
                    cleanup=False,
                    temp_dir=geometry_dir,
                )

                if output_pdb_path.exists():
                    final_pdb = self.file_manager.copy_to_output(output_pdb_path)
                    LOG.info(f"Fallback PDB copied to: {final_pdb}")
                    return system, final_pdb

            # Final fallback in case the system exists but no PDB could be written
            LOG.error("Failed to generate a valid PDB file")
            return system, None

        except Exception as e:
            LOG.error(f"Error during full generation: {str(e)}")
            import traceback

            LOG.error(f"Traceback: {traceback.format_exc()}")
            raise GeometryGenerationError(
                message=f"Failed to complete geometry generation: {str(e)}",
                original_error=e,
                error_code="GEO_ERR_001",
            )
        finally:
            os.chdir(original_dir)
            LOG.debug(f"Returned to original directory: {original_dir}")

            if (
                hasattr(self, "config")
                and hasattr(self.config, "debug")
                and self.config.debug
            ):
                LOG.info(f"Debug mode enabled, preserving temporary directories")
                for debug_dir in [
                    self.temp_dir,
                    self.file_manager.replacement_dir,
                    self.file_manager.mixing_dir,
                ]:
                    if debug_dir and debug_dir.exists():
                        try:
                            with open(debug_dir / "_DEBUG_INFO.txt", "w") as f:
                                f.write(f"Debug information for {debug_dir.name}\n")
                                f.write(f"Created at: {time.ctime()}\n")
                                if hasattr(self, "config"):
                                    f.write(f"Configuration:\n")
                                    for key, value in vars(self.config).items():
                                        if not key.startswith("_"):
                                            f.write(f"  {key}: {value}\n")
                        except Exception as e:
                            LOG.warning(
                                f"Failed to create debug info file in {debug_dir}: {e}"
                            )

    async def build_geometry(self) -> Tuple[Optional[System], Optional[Path]]:
        """
        Main entry point for geometry generation.

        This method orchestrates the entire geometry generation process based on
        the configuration settings, handling crystal building, mixing, and
        crosslink replacement.

        Returns:
            Generated system or None if only direct replacement was performed

        Raises:
            GeometryGenerationError: If any step in the process fails
        """
        try:
            if self.config.replace_bool and not self.config.geometry_generator:
                return await self._handle_direct_replacement()

            if self.config.mix_bool and not self.config.geometry_generator:
                return await self._handle_mixing_only()

            if (
                not self.config.geometry_generator
                and not self.config.replace_bool
                and not self.config.mix_bool
            ):
                LOG.info(
                    "Set -geometry flag to generate microfibrillar structure PDB file"
                )
                return None, None

            return await self._handle_full_generation()

        except GeometryGenerationError:
            raise
        except Exception as e:
            LOG.error(f"Unexpected error in geometry generation: {str(e)}")
            traceback.print_exc()
            raise GeometryGenerationError(
                message=f"Unexpected error in geometry generation: {str(e)}",
                original_error=e,
                error_code="GEO_ERR_001",
                context={"config": self.config.model_dump()},
            )


async def build_geometry(
    config: ColbuilderConfig, file_manager: Optional[FileManager] = None
) -> Optional[System]:
    """
    Build geometry from configuration.

    Standalone function that creates a GeometryService and invokes its
    build_geometry method.

    Args:
        config: Configuration for geometry operations
        file_manager: Optional file manager for consistent file handling

    Returns:
        Generated system or None
    """
    service = GeometryService(config, file_manager)
    system, _ = await service.build_geometry()
    return system


async def build_geometry_anywhere(
    config: ColbuilderConfig, file_manager: Optional[FileManager] = None
) -> Tuple[Path, Path]:
    """
    Main entry point for standalone geometry generation.

    This function serves as an alternative entry point for geometry generation outside
    of the usual pipeline. It creates a clean environment, handles file copying, and
    delegates to the CrystalBuilder to generate the structure.

    Args:
        config: Configuration settings with geometry parameters
        file_manager: Optional file manager to use for consistent file handling

    Returns:
        Tuple[Path, Path]: Tuple containing paths to geometry directory and final PDB file

    Raises:
        GeometryGenerationError: If input validation fails or geometry generation fails
    """
    original_dir: Path = Path.cwd()

    try:
        if file_manager is None:
            file_manager = FileManager(config)

        if not config.pdb_file:
            raise GeometryGenerationError(
                message="PDB file not specified for geometry generation",
                error_code="GEO_ERR_005",
            )

        if not config.contact_distance and not config.crystalcontacts_file:
            raise GeometryGenerationError(
                message="Either contact_distance or crystalcontacts_file must be provided",
                error_code="GEO_ERR_001",
                context={
                    "contact_distance": config.contact_distance,
                    "crystalcontacts_file": config.crystalcontacts_file,
                },
            )

        geometry_dir: Path = file_manager.ensure_geometry_dir()
        os.chdir(geometry_dir)

        try:
            pdb_path: Path = file_manager.copy_to_directory(
                config.pdb_file, dest_dir=geometry_dir
            )
            config.pdb_file = pdb_path

            crystal_builder: CrystalBuilder = CrystalBuilder(file_manager)

            system: System = await crystal_builder.build(config)
            LOG.info("Crystal structure built successfully")

            output_pdb_path: Path = file_manager.get_output_path(config.output, ".pdb")
            if not output_pdb_path.exists():
                raise GeometryGenerationError(
                    message=f"Expected output file not found: {output_pdb_path}",
                    error_code="GEO_ERR_001",
                )

            return geometry_dir, output_pdb_path

        finally:
            os.chdir(original_dir)
            LOG.debug(f"Returned to original directory: {original_dir}")

    except GeometryGenerationError:
        raise
    except Exception as e:
        LOG.error(
            f"Unexpected error during geometry generation: {str(e)}", exc_info=True
        )
        raise GeometryGenerationError(
            message=f"Failed to complete geometry generation: {str(e)}",
            original_error=e,
            error_code="GEO_ERR_001",
        )


async def mix_geometry(system: System, config: ColbuilderConfig) -> System:
    """
    Mix geometry types in system.

    Standalone function that creates a CrosslinkMixer and invokes its
    mix method.

    Args:
        system: System containing models to be mixed
        config: Configuration for mixing operation

    Returns:
        Mixed system
    """
    mixer = CrosslinkMixer()
    mixing_dir = Path(config.working_directory) / ".tmp" / "mixing_crosslinks"
    mixing_dir.mkdir(parents=True, exist_ok=True)

    return await mixer.mix(system, config, mixing_dir)


async def replace_geometry(
    system: System, config: ColbuilderConfig
) -> Optional[System]:
    """
    Replace crosslinks in system according to replacement ratio.

    Standalone function that handles both direct replacement (when no system
    is provided) and system-based replacement.

    Args:
        system: System containing models with crosslinks, or None for direct replacement
        config: Configuration for replacement

    Returns:
        System with replaced crosslinks, or None for direct replacement
    """
    file_manager = FileManager(config)
    original_dir = Path.cwd()

    try:
        replace_dir = Path(config.working_directory) / ".tmp" / "replace_crosslinks"
        replace_dir.mkdir(parents=True, exist_ok=True)
        file_manager.temp_dirs.add(replace_dir)

        LOG.info(f"Created clean replacement directory: {replace_dir}")

        os.chdir(replace_dir)

        replacer = CrosslinkReplacer()
        replacer.file_manager = file_manager

        temp_config = config.copy()
        temp_config.working_directory = replace_dir

        is_direct_replacement = (
            system is None
            and config.replace_file
            and os.path.exists(config.replace_file)
            and _is_pdb_file(config.replace_file)
        )

        if is_direct_replacement:
            replace_path = Path(config.replace_file)
            local_replace = replace_dir / replace_path.name
            shutil.copy2(replace_path, local_replace)
            temp_config.replace_file = local_replace

            LOG.info("Using direct replacement approach")
            system, output_pdb = await replacer.replace_direct(temp_config, replace_dir)

            if output_pdb and output_pdb.exists():
                output_prefix = temp_config.species or "output"
                if temp_config.output:
                    output_prefix = temp_config.output

                final_pdb = file_manager.get_output_path(output_prefix, ".pdb")
                shutil.copy2(output_pdb, final_pdb)
                LOG.info(f"Copied final output to: {final_pdb}")

                if final_pdb.exists():
                    try:
                        crystal = Crystal(pdb=str(final_pdb))
                        return System(crystal=crystal)
                    except Exception as e:
                        LOG.warning(
                            f"Could not create minimal system after direct replacement: {e}"
                        )

            return None
        else:
            if not system:
                raise GeometryGenerationError(
                    message="No system provided for replacement and input file is not a valid PDB",
                    error_code="GEO_ERR_004",
                )

            system = await replacer.replace_in_system(system, temp_config, replace_dir)

            output_prefix = temp_config.species or "output"
            if temp_config.output:
                output_prefix = temp_config.output

            output_pdb_path = replace_dir / f"{output_prefix}.pdb"
            LOG.info(f"Writing final system PDB to {output_pdb} - main_geometry l781")
            system.write_pdb(
                pdb_out=output_prefix,
                fibril_length=temp_config.fibril_length,
                cleanup=False,
                temp_dir=replace_dir,
            )
            LOG.info(f"PDB file written to {output_pdb_path}")

            if output_pdb_path.exists():
                final_pdb = file_manager.copy_to_output(output_pdb_path)
            else:
                LOG.warning(
                    f"Output PDB file not found at expected path: {output_pdb_path}"
                )

            return system
    finally:
        os.chdir(original_dir)
        LOG.debug(f"Returned to original directory: {original_dir}")


def _is_pdb_file(file_path: str) -> bool:
    """
    Check if a file is a PDB file based on its contents.

    Args:
        file_path: Path to the file to check

    Returns:
        True if the file appears to be a PDB file, False otherwise
    """
    try:
        with open(file_path, "r") as f:
            first_line = f.readline().strip()
            return first_line.startswith(("ATOM", "CRYST1", "HETATM"))
    except Exception:
        return False
