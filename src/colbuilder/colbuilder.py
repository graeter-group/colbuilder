"""
Colbuilder Main Module

This module serves as the main entry point for the Colbuilder system,
coordinating sequence, geometry, and topology generation operations.

The module manages the complete pipeline of collagen microfibril generation,
including sequence processing, geometry generation, and topology creation.

Key Features:
    - Configuration management
    - Operation coordination
    - Resource management
    - Progress tracking
    - Error handling

Example Usage:
    colbuilder --config_file config.yaml --sequence_generator --geometry_generator

Dependencies:
    - asyncio: For asynchronous operations
    - click: For command line interface
    - colorama: For terminal coloring
    - pydantic: For configuration management
"""

import sys
import os
import time
import logging
import asyncio
import traceback 
from pathlib import Path
from typing import Dict, Any, Tuple, Optional, Union, List
import click
import yaml
from colorama import init, Fore, Style
import shutil
from datetime import datetime

init()

# Package version
VERSION = "0.0.0"

from colbuilder.core.utils.logger import setup_logger
from colbuilder.core.utils.dec import timeit
from colbuilder.core.utils.files import FileManager
from colbuilder.core.utils.config import (
    ColbuilderConfig,
    get_config,
    OperationMode,
    load_yaml_config,
    resolve_relative_paths,
    validate_config,
)
from colbuilder.core.utils.exceptions import (
    ColbuilderError,
    ColbuilderErrorDetail,
    ConfigurationError,
    SystemError,
    SequenceGenerationError,
    GeometryGenerationError,
    TopologyGenerationError,
    ErrorCategory,
    ErrorSeverity,
)
from colbuilder.core.geometry.system import System

ConfigDict = Dict[str, Any]
RatioDict = Dict[str, int]
PathLike = Union[str, Path]

LOG = setup_logger(__name__)


def print_version(ctx: click.Context, param: click.Parameter, value: bool) -> None:
    """
    Print the version number and exit.

    Args:
        ctx: Click context
        param: Click parameter
        value: Flag value
    """
    if not value or ctx.resilient_parsing:
        return
    click.echo(f"colbuilder version {VERSION}")
    ctx.exit()


from colbuilder.core.sequence.main_sequence import build_sequence
from colbuilder.core.geometry.main_geometry import build_geometry_anywhere
from colbuilder.core.topology.main_topology import build_topology


def configure_loggers():
    """Configure external loggers to prevent duplicated output."""
    problematic_loggers = [
        # Add specific logger names that are causing duplicates
        # 'colbuilder.core.geometry.specific_module_causing_duplicates'
    ]

    for logger_name in problematic_loggers:
        logger = logging.getLogger(logger_name)
        logger.propagate = False

    for logger_name in logging.root.manager.loggerDict:
        if logger_name not in problematic_loggers and not logger_name.startswith(
            "colbuilder"
        ):
            logger = logging.getLogger(logger_name)
            if logger.level == logging.NOTSET:
                logger.setLevel(logging.WARNING)


def copy_config_to_tmp(config_file_path: Path, tmp_dir: Path) -> Optional[Path]:
    """
    Copy the configuration YAML file to the .tmp directory.
    Always copies the config regardless of debug mode.

    Args:
        config_file_path: Path to the configuration file
        tmp_dir: Path to the .tmp directory

    Returns:
        Path to the copied file or None if unsuccessful
    """
    if not config_file_path or not config_file_path.exists():
        LOG.debug(f"No config file to copy: {config_file_path}")
        return None

    tmp_dir.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    config_copy_path = tmp_dir / f"config_{timestamp}.yaml"

    try:
        shutil.copy2(config_file_path, config_copy_path)
        LOG.info(f"Config file saved: {config_copy_path}")
        return config_copy_path
    except Exception as e:
        LOG.error(f"Could not copy config file: {e}")
        return None


def parse_ratio_mix(ratio_str: str) -> RatioDict:
    """
    Parse mixing ratio string into a dictionary.

    Converts a string representation of mixing ratios into a
    dictionary mapping types to percentages.

    Args:
        ratio_str: String in format "Type:percentage Type:percentage"

    Returns:
        Dictionary mapping types to percentages

    Raises:
        GeometryGenerationError: If parsing fails or ratios invalid
    """
    try:
        ratio_mix = dict(item.split(":") for item in ratio_str.split())
        ratio_mix = {k: int(v) for k, v in ratio_mix.items()}
        if sum(ratio_mix.values()) != 100:
            raise ValueError("Mix ratios must sum to 100%")
        return ratio_mix
    except (ValueError, IndexError) as e:
        raise GeometryGenerationError(
            message="Invalid mixing ratio format",
            original_error=e,
            error_code="GEO_ERR_003",
            context={"ratio_string": ratio_str, "error_details": str(e)},
        )


def display_title() -> None:
    """Display the application title."""
    LOG.title("ColBuilder: Collagen Microfibril Builder")
    LOG.info(f"{Fore.CYAN}Version: {VERSION}{Style.RESET_ALL}")
    LOG.info("Copyright (c) 2024, ColBuilder Development Team")
    LOG.info("Distributed under the terms of the Apache License 2.0")
    LOG.info("")


@timeit
async def run_sequence_generation(config: ColbuilderConfig) -> Tuple[Optional[Path], Path]:
    """
    Generate coordinates for collagen molecule from sequence information.

    This function handles the sequence generation step of the pipeline,
    including homology modeling and structure optimization.

    Args:
        config: Configuration settings

    Returns:
        Tuple containing the MSA path (None for mutated PDB workflow) and the final PDB path

    Raises:
        SequenceGenerationError: If sequence generation fails
    """
    try:
        LOG.subsection("Generating Sequence")
        return await build_sequence(config)
    except Exception as e:
        LOG.error(f"Sequence generation failed: {str(e)}")
        if not isinstance(e, SequenceGenerationError):
            raise SequenceGenerationError(
                message=f"Sequence generation failed: {str(e)}",
                original_error=e,
                error_code="SEQ_ERR_001",
            )
        raise


@timeit
async def run_geometry_generation(
    config: ColbuilderConfig, file_manager: Optional[FileManager] = None
) -> Tuple[Optional[System], Path]:
    """
    Generate fibril geometry or handle mixing/replacement operations.

    Args:
        config: Configuration settings
        file_manager: Optional file manager for consistent file handling

    Returns:
        Tuple containing the generated system (which might be None) and the output PDB path

    Raises:
        GeometryGenerationError: If geometry generation or mixing fails
    """
    try:
        LOG.subsection("Building Geometry or Mixing")

        current_file_manager = file_manager or FileManager(config)

        # mixing-only
        if config.mix_bool and not config.geometry_generator:
            from colbuilder.core.geometry.main_geometry import GeometryService

            geometry_service = GeometryService(config, current_file_manager)
            system, pdb_path = await geometry_service._handle_mixing_only()
            return system, pdb_path

        # geometry generation (or combined geometry + mixing/replacement)
        if not config.pdb_file:
            raise GeometryGenerationError(
                message="PDB file not specified for geometry generation",
                error_code="GEO_ERR_005",
            )

        pdb_path = Path(config.pdb_file).resolve()
        if not pdb_path.exists():
            raise GeometryGenerationError(
                message=f"PDB file not found: {pdb_path}", error_code="GEO_ERR_005"
            )

        if config.contact_distance is None and not config.crystalcontacts_file:
            raise GeometryGenerationError(
                message="Either contact_distance or crystalcontacts_file must be provided",
                error_code="GEO_ERR_001",
                context={
                    "contact_distance": config.contact_distance,
                    "crystalcontacts_file": config.crystalcontacts_file,
                },
            )

        output_path, pdb_path = await build_geometry_anywhere(
            config, current_file_manager
        )

        if not pdb_path.exists():
            raise GeometryGenerationError(
                message=f"Expected output file not found: {pdb_path}",
                error_code="GEO_ERR_001",
            )

        try:
            from colbuilder.core.geometry.crystal import Crystal

            crystal = Crystal(pdb=str(pdb_path))
            system = System(crystal=crystal)
        except Exception as e:
            LOG.warning(f"Could not create system object from output PDB: {e}")
            system = None

        return system, pdb_path

    except Exception as e:
        LOG.error(f"Geometry generation failed: {str(e)}")
        if not isinstance(e, GeometryGenerationError):
            raise GeometryGenerationError(
                message=f"Geometry generation failed: {str(e)}",
                original_error=e,
                error_code="GEO_ERR_001",
            )
        raise


@timeit
#TODO: double check this!!
async def run_topology_generation(
    config: ColbuilderConfig,
    system_path: Path,
    existing_system: Optional[System] = None,
    file_manager: Optional[FileManager] = None,
) -> Tuple[Path, Path]:
    """
    Generate topology files for molecular dynamics simulations.
    
    Args:
        config: Configuration settings
        system_path: Path to the fibril PDB file
        existing_system: Optional pre-loaded System object
        file_manager: Optional file manager for handling directories
        
    Returns:
        Tuple of (topology_directory, system_path)
        
    Raises:
        TopologyGenerationError: If topology generation fails
    """
    try:
        if file_manager is None:
            file_manager = FileManager(config)

        # Check for existing cap files from geometry/mixing steps
        mixing_dir = Path(".tmp") / "mixing_crosslinks"
        if not mixing_dir.exists():
            mixing_dir = Path.cwd() / ".tmp" / "mixing_crosslinks"

        geometry_dir = Path(".tmp") / "geometry_gen"
        if not geometry_dir.exists():
            geometry_dir = Path.cwd() / ".tmp" / "geometry_gen"

        topology_dir_path = file_manager.get_temp_dir("topology_gen")

        cap_files_found = False
        
        # Copy cap files from mixing directory if available
        if mixing_dir.exists():
            LOG.debug(f"Searching for cap files in mixing directory: {mixing_dir}")
            cap_files = list(mixing_dir.glob("**/*.caps.pdb"))
            LOG.debug(f"Found {len(cap_files)} cap files in mixing directory")
            
            if cap_files:
                cap_files_found = True
                for cap_file in cap_files:
                    relative_path = cap_file.relative_to(mixing_dir)
                    dest_file = topology_dir_path / relative_path
                    dest_file.parent.mkdir(parents=True, exist_ok=True)
                    
                    try:
                        shutil.copy2(cap_file, dest_file)
                        LOG.debug(f"Copied cap file: {relative_path}")
                    except Exception as e:
                        LOG.warning(f"Failed to copy {cap_file}: {e}")
        
        # Copy cap files from geometry directory if available
        if geometry_dir.exists():
            LOG.debug(f"Searching for cap files in geometry directory: {geometry_dir}")
            cap_files = list(geometry_dir.glob("**/*.caps.pdb"))
            LOG.debug(f"Found {len(cap_files)} cap files in geometry directory")
            
            if cap_files:
                cap_files_found = True
                for cap_file in cap_files:
                    relative_path = cap_file.relative_to(geometry_dir)
                    dest_file = topology_dir_path / relative_path
                    dest_file.parent.mkdir(parents=True, exist_ok=True)
                    
                    try:
                        shutil.copy2(cap_file, dest_file)
                        LOG.debug(f"Copied cap file: {relative_path}")
                    except Exception as e:
                        LOG.warning(f"Failed to copy {cap_file}: {e}")

        if not cap_files_found and existing_system and system_path.exists():
            LOG.info("No cap files found, extracting models from fibril PDB")
            await extract_and_cap_models_from_pdb(
                system_path, existing_system, config, topology_dir_path
            )

        if existing_system:
            LOG.info("Using existing system from geometry generation")
            system = existing_system

            LOG.debug(f"System has {len(list(system.get_models()))} models")
            for model_id in system.get_models():
                model = system.get_model(model_id=model_id)
                if model:
                    LOG.debug(
                        f"Model {model_id} - Type: {model.type if hasattr(model, 'type') else 'Unknown'}"
                    )
        else:
            LOG.info("Creating new system from PDB file")
            from colbuilder.core.geometry.crystal import Crystal
            from colbuilder.core.geometry.system import System

            crystal = Crystal(pdb=str(system_path))
            system = System(crystal=crystal)

        # Generate topology files
        await build_topology(system, config, file_manager)

        # Topology files are created in [species]_topology_files directory
        topology_dir = Path(f"{config.species}_topology_files")
        if not topology_dir.exists():
            LOG.warning(f"Topology directory not found: {topology_dir}")
            topology_dir = Path()

        return topology_dir, system_path

    except Exception as e:
        LOG.error(f"Topology generation failed: {str(e)}")
        if not isinstance(e, TopologyGenerationError):
            raise TopologyGenerationError(
                message=f"Topology generation failed: {str(e)}",
                original_error=e,
                error_code="TOP_ERR_001",
            )
        raise


async def extract_and_cap_models_from_pdb(
    pdb_path: Path,
    system: System,
    config: ColbuilderConfig,
    output_dir: Path
) -> None:
    """
    Prepare fibril PDB for topology generation.
    
    Extracts individual models from a multi-model PDB, detects crosslink-based
    connectivity, and organizes them by crosslink type. Works for both Amber99 
    and Martini3 force fields.
    
    Args:
        pdb_path: Path to the fibril PDB file
        system: System object with model information
        config: Configuration object
        output_dir: Directory for output files
    """
    from colbuilder.core.utils.crosslink_detector import CrosslinkDetector
    
    LOG.info(f"Preparing fibril structure for topology generation: {pdb_path}")
    
    detector = CrosslinkDetector()
    type_dir = detector.prepare_for_topology(pdb_path, output_dir)
    
    LOG.info("Analyzing crosslink connectivity between models...")
    connections = detector.find_crosslink_connections(type_dir, cutoff=5.0)
    
    # Update system models with connectivity information
    for model_id in system.get_models():
        model = system.get_model(model_id=model_id)
        if model:
            int_model_id = int(model_id)
            if int_model_id in connections:
                model.connect = sorted(connections[int_model_id])
                LOG.debug(f"Model {int_model_id} connected to: {model.connect}")
            else:
                model.connect = [model_id]
                LOG.debug(f"Model {int_model_id} has no connections (self only)")
    
    connected_count = sum(1 for m_id in system.get_models() 
                         if len(system.get_model(model_id=m_id).connect) > 1)
    LOG.info(f"Connectivity analysis complete: {connected_count}/{len(list(system.get_models()))} models have connections")
    
    LOG.info(f"Ready for {config.force_field} topology generation")


async def run_pipeline(config: ColbuilderConfig) -> Dict[str, Path]:
    """
    Run the complete Colbuilder pipeline based on configuration.
    
    Args:
        config: Configuration settings
        
    Returns:
        Dictionary containing paths to all generated outputs
        
    Raises:
        Various ColbuilderError subclasses depending on which step fails
    """
    results = {
        "sequence_msa": None,
        "sequence_pdb": None,
        "geometry_pdb": None,
        "topology_dir": None,
    }

    try:
        file_manager = FileManager(config)
        current_system = None

        # Sequence Generation
        if config.sequence_generator:
            LOG.section("Running sequence generation")
            sequence_msa, sequence_pdb = await run_sequence_generation(config)
            results["sequence_msa"] = sequence_msa
            results["sequence_pdb"] = sequence_pdb

            # Update PDB file for next steps if needed
            if sequence_pdb and not config.pdb_file:
                config.pdb_file = sequence_pdb
                LOG.info(
                    f"Using generated sequence PDB for further processing: {sequence_pdb}"
                )

        # Topology-only mode
        if (config.topology_generator and 
            not config.geometry_generator and 
            not config.mix_bool and 
            not config.replace_bool):
            
            LOG.section("Running topology-only mode")
            
            if not config.pdb_file:
                raise TopologyGenerationError(
                    message="PDB file required for topology-only mode",
                    error_code="TOP_ERR_008",
                )
            
            pdb_path = Path(config.pdb_file).resolve()
            if not pdb_path.exists():
                raise TopologyGenerationError(
                    message=f"PDB file not found: {pdb_path}",
                    error_code="TOP_ERR_008",
                )
            
            LOG.info(f"Using existing fibril PDB: {pdb_path}")
            results["geometry_pdb"] = pdb_path
            
            from colbuilder.core.geometry.crystal import Crystal
            from colbuilder.core.geometry.system import System
            from colbuilder.core.geometry.model import Model
            from colbuilder.core.utils.crosslink_detector import CrosslinkDetector
            
            detector = CrosslinkDetector()
            
            structure_type = detector.detect_structure_type(pdb_path)
            LOG.info(f"Detected structure type: {structure_type}")
            
            crystal = Crystal(pdb=str(pdb_path))
            current_system = System(crystal=crystal)
            
            model_ids = detector.get_model_ids(pdb_path)
            
            if not model_ids:
                # Count by TER records (3 TERs = 1 triple helix)
                LOG.info("No MODEL records found, counting TER-separated triple helices")
                ter_count = 0
                with open(pdb_path, 'r') as f:
                    for line in f:
                        if line.startswith('TER'):
                            ter_count += 1
                
                num_models = ter_count // 3
                model_ids = list(range(num_models))
                LOG.info(f"Detected {num_models} triple helix models from TER records")
            
            for model_id in model_ids:
                model = Model(
                    id=float(model_id),
                    transformation=crystal.get_default_transformation(),
                    unit_cell=crystal.get_s_matrix(
                        t_matrix=crystal.get_default_transformation()
                    ),
                    pdb_file=str(pdb_path),
                )
                model.type = structure_type
                model.crosslink_type = structure_type
                current_system.add_model(model=model)
                LOG.debug(f"Added model {model_id} with type {structure_type}")
            
            LOG.info(f"Created system with {len(model_ids)} models of type {structure_type}")
            results["geometry_system"] = current_system

        # Handle direct replacement without geometry generation
        elif config.replace_bool and not config.geometry_generator:
            LOG.section("Running direct replacement without geometry generation")
            from colbuilder.core.geometry.main_geometry import GeometryService

            geometry_service = GeometryService(config, file_manager)
            current_system, pdb_path = (
                await geometry_service._handle_direct_replacement()
            )
            results["geometry_pdb"] = pdb_path
            results["geometry_system"] = current_system
            LOG.info(f"Direct replacement completed, output PDB: {pdb_path}")

        # Mix-only 
        elif config.mix_bool and not config.geometry_generator:
            LOG.section("Running crosslinks mixing mode")
            current_system, pdb_path = await run_geometry_generation(
                config, file_manager
            )
            results["geometry_pdb"] = pdb_path
            results["geometry_system"] = current_system
            LOG.info(f"Mixing completed, output PDB: {pdb_path}")

        # Full geometry generation (optionally with mixing/replacement)
        elif config.geometry_generator:
            LOG.section("Running geometry generation")
            from colbuilder.core.geometry.main_geometry import GeometryService

            geometry_service = GeometryService(config, file_manager)
            current_system, pdb_path = await geometry_service._handle_full_generation()
            results["geometry_pdb"] = pdb_path
            results["geometry_system"] = current_system
            LOG.info(f"Geometry generation completed, output PDB: {pdb_path}")

        # Topology Generation
        if config.topology_generator and results["geometry_pdb"]:
            LOG.section("Running topology generation")

            if "geometry_system" in results and results["geometry_system"]:
                topology_dir, system_path = await run_topology_generation(
                    config,
                    results["geometry_pdb"],
                    existing_system=results["geometry_system"],
                    file_manager=file_manager,
                )
            else:
                topology_dir, system_path = await run_topology_generation(
                    config, results["geometry_pdb"], file_manager=file_manager
                )

            results["topology_dir"] = topology_dir

        # Cleanup temporary files unless in debug mode
        if not config.debug:
            file_manager.cleanup()

        return results

    except SequenceGenerationError as e:
        LOG.error(f"Sequence generation error: {str(e)}")
        e.log_error()
        raise
    except GeometryGenerationError as e:
        LOG.error(f"Geometry generation error: {str(e)}")
        e.log_error()
        raise
    except TopologyGenerationError as e:
        LOG.error(f"Topology generation error: {str(e)}")
        e.log_error()
        raise
    except ColbuilderError as e:
        LOG.error(f"Colbuilder error: {str(e)}")
        e.log_error()
        raise
    except Exception as e:
        LOG.critical(f"Unhandled exception in pipeline: {str(e)}")
        LOG.info(f"Exception details: {traceback.format_exc()}")
        raise SystemError(
            message=f"Unhandled exception in pipeline: {str(e)}",
            original_error=e,
            error_code="SYS_ERR_001",
            context={
                "config": (
                    config.model_dump()
                    if hasattr(config, "model_dump")
                    else str(config)
                )
            },
        )


def log_configuration_summary(cfg: ColbuilderConfig) -> None:
    """
    Log a summary of the configuration.

    This function provides a summary of the current
    configuration settings for user verification.

    Args:
        cfg: Configuration to summarize
    """
    # Detect topology-only mode
    topology_only = (
        cfg.topology_generator and 
        not cfg.geometry_generator and 
        not cfg.sequence_generator and
        not cfg.mix_bool and 
        not cfg.replace_bool
    )
    
    if topology_only:
        # Minimal configuration for topology-only mode
        sections = {
            "Input Configuration": lambda: [
                f"Species: {cfg.species}",
                f"Force Field: {cfg.force_field}",
                f"Input PDB: {cfg.pdb_file}" if cfg.pdb_file else None,
            ],
            "Operation Modes": lambda: [
                "Topology Generation \u2713",
            ],
        }
    else:
        # Full configuration for geometry/sequence generation
        sections = {
            "Fibril Parameters": lambda: [
                f"Species: {cfg.species}",
                (
                    f"Contact Distance: {cfg.contact_distance}"
                    if cfg.contact_distance
                    else None
                ),
                f"Fibril Length: {cfg.fibril_length}",
                "Crosslinks:",
                f"    Mix Ratio: {cfg.ratio_mix}" if cfg.mix_bool else None,
                f"    Mix Files: {cfg.files_mix}" if cfg.mix_bool else None,
                f"    Replace Ratio: {cfg.ratio_replace}%" if cfg.replace_bool else None,
                (
                    f"    N-terminal: {cfg.n_term_type}, {cfg.n_term_combination}"
                    if (cfg.crosslink and not cfg.mix_bool)
                    else f"    N-terminal: No additional crosslinks"
                ),
                (
                    f"    C-terminal: {cfg.c_term_type}, {cfg.c_term_combination}"
                    if (cfg.crosslink and not cfg.mix_bool)
                    else f"    C-terminal: No additional crosslinks"
                ),
            ],
            "Operation Modes": lambda: [
                "Sequence Generation \u2713" if cfg.sequence_generator else None,
                "Geometry Generation \u2713" if cfg.geometry_generator else None,
                "Mix Crosslinks \u2713" if cfg.mix_bool else None,
                "Replace Crosslinks \u2713" if cfg.replace_bool else None,
                "Topology Generation \u2713" if cfg.topology_generator else None,
            ],
        }

    for section, get_items in sections.items():
        LOG.subsection(section)
        for item in filter(None, get_items()):
            LOG.info(f"- {item}")

    if cfg.config_file:
        LOG.info(f"Config File: {cfg.config_file}")
    if cfg.pdb_file and (cfg.geometry_generator or topology_only):
        LOG.info(f"Input File: {cfg.pdb_file}")
    if not topology_only:
        LOG.info(f"Output File: {cfg.output}.pdb")
    LOG.info(f"Working Directory: {cfg.working_directory}")


def initialize_logging(debug=False, working_dir=None, config_file=None):
    """
    Initialize logging system for the application.
    
    Args:
        debug: Enable debug logging if True
        working_dir: Working directory for log files
        config_file: Configuration file to copy to temp directory
        
    Returns:
        Configured logger instance
    """
    tmp_dir = working_dir / ".tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    if config_file:
        copy_config_to_tmp(config_file, tmp_dir)

    if debug:
        os.environ["COLBUILDER_DEBUG"] = "1"

        for logger_name in logging.root.manager.loggerDict:
            if logger_name.startswith("colbuilder"):
                logger = logging.getLogger(logger_name)
                logger.setLevel(logging.DEBUG)

    from colbuilder.core.utils.logger import initialize_root_logger

    return initialize_root_logger(debug=debug, log_dir=tmp_dir)


@click.command()
@click.option(
    "--config_file",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="YAML configuration file",
)
@click.option("--species", type=str, help="Species name (e.g., homo_sapiens)")
@click.option("--sequence_generator", is_flag=True, help="Run sequence generation")
@click.option("--geometry_generator", is_flag=True, help="Run geometry generation")
@click.option(
    "-fasta",
    "--fasta_file",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Fasta-input file for collagen triple helix sequence",
)
@click.option(
    "--crosslink_copies",
    nargs=2,
    type=click.Choice(["D0", "D1", "D2", "D3", "D4", "D5"]),
    default=["D0", "D5"],
    help="Pair of unit cell translations for crosslink optimization (default: D0 D5)"
)
@click.option(
    "--mutated_pdb",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Pre-mutated PDB file to skip homology modeling and apply additional crosslinks"
)
@click.option(
    "--additional_1_type",
    type=str,
    help="First additional crosslink type to apply to mutated PDB"
)
@click.option(
    "--additional_2_type", 
    type=str,
    help="Second additional crosslink type to apply to mutated PDB"
)
@click.option(
    "--additional_1_combination",
    type=str,
    help="First additional crosslink position (e.g., '9.C - 947.A')"
)
@click.option(
    "--additional_2_combination",
    type=str,
    help="Second additional crosslink position (e.g., '1047.C - 104.C')"
)
@click.option(
    "-pdb",
    "--pdb_file",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="PDB-input file for single triple helix or template fibril",
)
@click.option(
    "-wd",
    "--working_directory",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    default=Path.cwd(),
    help="Set working directory",
)
@click.option(
    "-dc",
    "--contact_distance",
    type=float,
    help="Contact distance as input for radial size of microfibril",
)
@click.option("-length", "--fibril_length", type=float, help="Length of microfibril")
@click.option(
    "-contacts",
    "--crystalcontacts_file",
    type=click.Path(path_type=Path),
    help="Read crystalcontacts from file",
)
@click.option(
    "-connect",
    "--connect_file",
    type=click.Path(path_type=Path),
    help="Read connect between contacts from file",
)
@click.option(
    "-optimize",
    "--crystalcontacts_optimize",
    is_flag=True,
    help="Optimize crystalcontacts",
)
@click.option(
    "-space",
    "--solution_space",
    nargs=3,
    type=float,
    default=[1, 1, 1],
    help="Solution space of optimisation problem [ d_x d_y d_z ]",
)
@click.option(
    "-mix", "--mix_bool", is_flag=True, help="Generate a mixed crosslinked microfibril"
)
@click.option(
    "-ratio_mix",
    "--ratio_mix",
    type=str,
    help='Ratio for mix-crosslink setup in format "Type:percentage Type:percentage"',
)
@click.option(
    "-files_mix",
    "--files_mix",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    multiple=True,
    help="PDB-files with different crosslink-types",
)
@click.option(
    "-replace",
    "--replace_bool",
    is_flag=True,
    help="Generate a microfibril with less crosslinks",
)
@click.option(
    "-ratio_replace",
    "--ratio_replace",
    type=float,
    help="Ratio of crosslinks to be replaced with Lysines",
)
@click.option(
    "-replace_file",
    "--replace_file",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="File with information about crosslinks to be replaced with Lysine",
)
@click.option(
    "-topology", "--topology_generator", is_flag=True, help="Generate topology files"
)
@click.option(
    "-ff",
    "--force_field",
    help="Specify force field to be used, e.g. -ff amber99 OR -ff martini3",
)
@click.option(
    "-top_debug",
    "--topology_debug",
    is_flag=True,
    help="Keep intermediate topology files for debugging",
)
@click.option("--debug", is_flag=True, help="Enable debug logging")
@click.option(
    "--version",
    is_flag=True,
    callback=print_version,
    expose_value=False,
    is_eager=True,
    help="Show the version and exit.",
)

@timeit
def main(**kwargs: Any) -> int:
    """
    Main entry point for Colbuilder.

    This function initializes the system, sets up configuration,
    and coordinates the execution of all requested operations.

    Args:
        **kwargs: Command line arguments and options

    Returns:
        Exit code (0 for success, 1 for failure)
    """
    try:
        debug = kwargs.get("debug", False)
        working_directory = kwargs.get("working_directory", Path.cwd())
        config_file = kwargs.get("config_file")

        global LOG
        LOG = initialize_logging(
            debug=debug, working_dir=working_directory, config_file=config_file
        )

        configure_loggers()

        display_title()

        raw_files_mix = None

        if config_file:
            config_path = Path(config_file)
            config_dir = config_path.parent

            try:
                LOG.info(f"Loading configuration from: {config_path}")
                file_config = load_yaml_config(config_path)
                file_config = resolve_relative_paths(file_config, config_dir)

                if "files_mix" in file_config and file_config["files_mix"]:
                    LOG.debug(f"{file_config['files_mix']}")
                    raw_files_mix = file_config["files_mix"]

                config_data = file_config

            except Exception as e:
                LOG.error(f"Error loading configuration file: {str(e)}")
                if debug:
                    LOG.debug(f"Exception details: {traceback.format_exc()}")
                return 1
        else:
            config_data = {}

        # Override with command-line arguments if they are explicitly provided
        for key, value in kwargs.items():
            if key == "config_file":
                continue

            elif key == "files_mix" and value:
                LOG.info(f"Command line files_mix: {value}")
                config_data[key] = value
                raw_files_mix = value
                continue

        LOG.debug(f"Final configuration data before validation: {config_data}")

        if config_data.get("mix_bool", False) and "files_mix" in config_data and config_data["files_mix"]:
            files = []
            for file_path in config_data["files_mix"]:
                if isinstance(file_path, str):
                    file_path = Path(file_path)

                if not file_path.is_absolute():
                    base_dir = config_dir if config_file else Path.cwd()
                    file_path = (base_dir / file_path).resolve()

                if file_path.exists():
                    files.append(file_path)
                    LOG.info(f"Found mix file: {file_path}")
                else:
                    LOG.warning(f"Mix file not found: {file_path}")

            config_data["files_mix"] = tuple(files)
            LOG.debug(
                f"Final files_mix before validation: {config_data.get('files_mix')}"
            )
        elif config_data.get("mix_bool", False):
            LOG.debug("mix_bool is True but no files_mix provided")
        else:
            if "files_mix" in config_data:
                config_data["files_mix"] = None
            LOG.debug("mix_bool is False, ignoring files_mix")

        try:
            config = validate_config(config_data)

            if config.mix_bool and (not config.files_mix or len(config.files_mix) == 0):
                if raw_files_mix:
                    LOG.debug(f"Using raw files_mix: {raw_files_mix}")
                    files = []

                    for file_path in raw_files_mix:
                        if isinstance(file_path, str):
                            file_path = Path(file_path)

                        if not file_path.is_absolute():
                            base_dir = config_dir if config_file else Path.cwd()
                            file_path = (base_dir / file_path).resolve()

                        files.append(file_path)

                    config.files_mix = tuple(files)

            global_config = get_config(existing_config=config)

        except Exception as e:
            LOG.error(f"Configuration setup failed: {str(e)}")
            LOG.info(f"Exception details: {traceback.format_exc()}")
            return 1

        log_configuration_summary(config)

        # Only validate/create mix files if mixing is actually enabled
        if config.mix_bool:
            LOG.debug(f"Files Mix: {config.files_mix}")

            if not config.files_mix or len(config.files_mix) == 0:
                LOG.warning("No mix files found, but mixing is enabled!")

                if raw_files_mix:
                    LOG.info("Creating empty mix files for testing...")
                    files = []

                    for file_name in raw_files_mix:
                        if isinstance(file_name, Path):
                            file_path = file_name
                        else:
                            file_path = Path(file_name)

                        if not file_path.is_absolute():
                            file_path = (
                                config.working_directory / file_path.name
                            ).resolve()

                        if not file_path.exists():
                            try:
                                LOG.info(f"Creating empty file: {file_path}")
                                with open(file_path, "w") as f:
                                    f.write(
                                        "REMARK This is an empty PDB file created for testing\n"
                                    )
                                    f.write("END\n")
                                files.append(file_path)
                            except Exception as e:
                                LOG.error(f"Failed to create test file: {e}")
                        else:
                            files.append(file_path)

                    if files:
                        config.files_mix = tuple(files)
                        LOG.info(f"Created test files: {config.files_mix}")
                    files = []

                    for file_name in raw_files_mix:
                        if isinstance(file_name, Path):
                            file_path = file_name
                        else:
                            file_path = Path(file_name)

                        if not file_path.is_absolute():
                            file_path = (
                                config.working_directory / file_path.name
                            ).resolve()

                        if not file_path.exists():
                            try:
                                LOG.info(f"Creating empty file: {file_path}")
                                with open(file_path, "w") as f:
                                    f.write(
                                        "REMARK This is an empty PDB file created for testing\n"
                                    )
                                    f.write("END\n")
                                files.append(file_path)
                            except Exception as e:
                                LOG.error(f"Failed to create test file: {e}")
                        else:
                            files.append(file_path)

                    if files:
                        config.files_mix = tuple(files)
                        LOG.info(f"Created test files: {config.files_mix}")

        LOG.subsection("Starting ColBuilder Pipeline")
        results = asyncio.run(run_pipeline(config))

        LOG.section("ColBuilder Pipeline Complete")
        LOG.info(
            f"{Fore.MAGENTA}Done! Colbuilder completed successfully.{Style.RESET_ALL}"
        )
        return 0

    except ColbuilderError as e:
        if "LOG" in globals():
            LOG.error(f"Colbuilder error: {str(e)}")
            e.log_error()
        else:
            print(f"ERROR: {str(e)}")
        return 1

    except Exception as e:
        if "LOG" in globals():
            LOG.critical(f"Unhandled exception: {str(e)}")
            from colbuilder.core.utils.logger import log_exception

            log_exception(LOG, e)
        else:
            print(f"CRITICAL ERROR: {str(e)}")
            import traceback

            traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())