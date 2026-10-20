"""
Colbuilder Main Geometry Module

Main entry point for geometry operations, coordinating the crystal,
mixing, and replacement services.
"""

import os
import time
import shutil
import traceback
import logging
from pathlib import Path
from typing import Optional, List, Tuple

from colbuilder.core.utils.exceptions import GeometryGenerationError
from colbuilder.core.utils.config import ColbuilderConfig
from colbuilder.core.utils.logger import setup_logger
from colbuilder.core.utils.files import FileManager
from .crystal_builder import CrystalBuilder
from .crosslink_mixer import CrosslinkMixer
from .geometry_replacer import CrosslinkReplacer
from .unpaired_crosslinks import UnpairedCrosslinkFinder
from .system import System
from .crystal import Crystal
from colorama import Fore, Style

LOG = setup_logger(__name__)


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
        self.file_manager = file_manager or FileManager(config)
        self.original_dir = Path.cwd()
        self.temp_dir = None

        # Set file manager on all services
        self.crystal_service.set_file_manager(self.file_manager)

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
            
            working_dir_root = Path(self.config.working_directory).resolve()
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

                # Auto-detect unpaired crosslinks
                auto_manual_replacements: List[str] = []
                auto_manual_file: Optional[Path] = None
                if system:
                    try:
                        finder = UnpairedCrosslinkFinder(base_dir=working_dir_root)
                        auto_manual_replacements, auto_manual_file = finder.run()
                        if auto_manual_replacements:
                            user_requested_replacement = bool(
                                (
                                    temp_config.ratio_replace is not None
                                    and temp_config.ratio_replace > 0
                                )
                                or temp_config.replace_file
                                or temp_config.manual_replacements
                            )
                            if getattr(temp_config, "ratio_replace_mode", "random") == "preserve_attachment":
                                # Preserve the local cleanup phase without turning it into
                                # a manual override of the user's exact ratio selection.
                                temp_config._auto_unpaired_replacements = auto_manual_replacements
                                temp_config.replace_bool = True
                                temp_config.auto_fix_unpaired = True
                                self.config.auto_fix_unpaired = True
                                if temp_config.ratio_replace is None and temp_config.replace_file is None:
                                    temp_config.ratio_replace = 0
                                LOG.info(
                                    "Auto-detected %d unpaired marker(s); cleaning them before attachment-preserving replacement.",
                                    len(auto_manual_replacements),
                                )
                            elif user_requested_replacement:
                                # Setting manual_replacements here would silently
                                # discard the user's own ratio_replace/replace_file/
                                # manual_replacements request: manual_replacements
                                # takes precedence over ratio_replace by design (see
                                # geometry_replacer.py), so auto-populating it with
                                # just the auto-detected unpaired markers would mean
                                # the user's explicit request never runs at all.
                                # Leave their config alone instead -- the unpaired
                                # markers only get fixed if the user's own selection
                                # happens to include them.
                                LOG.warning(
                                    "Auto-detected %d unpaired crosslink marker(s), "
                                    "but not auto-fixing them: a replacement was "
                                    "already requested (ratio_replace/replace_file/"
                                    "manual_replacements), which takes precedence. "
                                    "These markers stay unpaired unless your own "
                                    "replacement selection covers them.",
                                    len(auto_manual_replacements),
                                )
                            else:
                                temp_config.manual_replacements = auto_manual_replacements
                                temp_config.replace_bool = True
                                temp_config.auto_fix_unpaired = True

                                # Propagate to shared config so downstream (topology) uses replace_manual
                                try:
                                    self.config.auto_fix_unpaired = True
                                except Exception:
                                    pass

                                if (
                                    temp_config.ratio_replace is None
                                    and temp_config.replace_file is None
                                ):
                                    temp_config.ratio_replace = 0

                                LOG.info(
                                    "Auto-detected %d unpaired crosslink marker(s); enabling replacement.",
                                    len(auto_manual_replacements),
                                )
                        else:
                            LOG.debug("No unpaired crosslinks detected")
                    except Exception as e:
                        LOG.error(
                            "Automatic unpaired crosslink detection failed: %s", e, exc_info=True
                        )
                        raise
                else:
                    temp_config.auto_fix_unpaired = False

                # Only write geometry-only PDB if no mixing or replacement will follow
                if not (temp_config.mix_bool or temp_config.replace_bool):
                    output_prefix = temp_config.output or temp_config.species
                    output_pdb_path = geometry_dir / f"{output_prefix}.pdb"

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
                system, _ = await self.mixer_service.mix(system, temp_config, mixing_dir)
                final_temp_dir = mixing_dir
                LOG.info("Mixing completed.")

            # Determine if replacement is requested
            replacement_requested = bool(
                temp_config.replace_bool
                or getattr(temp_config, "auto_fix_unpaired", False)
                or (
                    temp_config.ratio_replace is not None
                    and float(temp_config.ratio_replace) > 0
                )
                or temp_config.replace_file
                or getattr(temp_config, "manual_replacements", None)
            )
            
            # Set verbosity flag for replacement logging
            try:
                temp_config._replacement_verbose = replacement_requested  # type: ignore[attr-defined]
            except Exception:
                pass

            # (Unpaired-crosslink auto-detection now runs earlier, right after the
            # system is built and before the geometry-only early return, so that a
            # detected unpaired crosslink can enable the replacement stage.)

            # Enable replacement stage when a ratio was provided without explicit replace_bool
            if (
                not temp_config.replace_bool
                and temp_config.ratio_replace is not None
                and float(temp_config.ratio_replace) > 0
            ):
                temp_config.replace_bool = True

            # For replacement we need the original working dir (not the geometry subdir)
            temp_config.working_directory = working_dir_root

            # Perform replacement if required
            if temp_config.replace_bool and system:
                model_type = (
                    system.get_model(model_id=0.0).type
                    if hasattr(system.get_model(model_id=0.0), "type")
                    else "D"
                )

                replacement_dir = self.file_manager.ensure_replacement_dir()
                LOG.debug(f"Using replacement directory: {replacement_dir}")

                self.replacer_service.file_manager = self.file_manager

                system = await self.replacer_service.replace_in_system(
                    system, temp_config, replacement_dir
                )

                # Check if replacement was skipped
                if getattr(temp_config, "_replacement_skipped", False):
                    temp_config.replace_bool = False
                    LOG.debug("Replacement was skipped - no matching crosslinks found")
                else:
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
                                final_pdb = self.file_manager.copy_to_output(
                                    output_pdb_path
                                )
                                if getattr(temp_config, "_replacement_verbose", True):
                                    LOG.info(
                                        f"Final PDB with replacements written to: {final_pdb}"
                                    )
                                return system, final_pdb
                        else:
                            LOG.warning(
                                "No caps files found in replacement directory. Falling back to geometry directory."
                            )
                    else:
                        LOG.debug(
                            f"Replacement type directory not found: {replacement_type_dir}"
                        )

                    # If we get here, replacement was requested but replacement directory didn't have caps
                    # Write to the geometry directory instead
                    output_prefix = temp_config.output or temp_config.species
                    output_pdb_path = geometry_dir / f"{output_prefix}.pdb"

                    if getattr(temp_config, "_replacement_verbose", True):
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
                        if getattr(temp_config, "_replacement_verbose", True):
                            LOG.info(f"Final PDB with replacements copied to: {final_pdb}")
                        return system, final_pdb

            # Write the final system PDB if we've done mixing (since replacement handles its own writing)
            if system and temp_config.mix_bool:
                output_prefix = temp_config.output or temp_config.species
                output_pdb_path = final_temp_dir / f"{output_prefix}.pdb"

                LOG.info(f"Writing final mixed system PDB to {output_pdb_path}")
                system.write_pdb(
                    pdb_out=output_pdb_path,
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

            # Fallback: If we reach this point, we have a system but no PDB was written
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
            Tuple[Optional[System], Optional[Path]]: Generated system and output PDB path

        Raises:
            GeometryGenerationError: If any step in the process fails
        """
        try:
            # Clear stale caps/type-subdirs from previous runs so a re-run cannot
            # silently reuse old geometry (different fibril length, contact distance, etc.).
            try:
                if getattr(self, "file_manager", None) is not None:
                    self.file_manager.clear_generation_dirs()
            except Exception as e:
                LOG.debug(f"Could not clear generation dirs at run start: {e}")

            # Allow direct replacement when a PDB is provided, even if replace_bool was not set explicitly
            if self.config.replace_file and not self.config.geometry_generator:
                self.config.replace_bool = True

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
