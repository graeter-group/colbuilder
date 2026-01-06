"""
Module for generating collagen structures from input sequences.

The SequenceGenerator class orchestrates the entire process of generating collagen 
triple helical structures, including sequence alignment, structure modeling, 
crosslink application and optimization.

Key Features:
- Sequence alignment using MUSCLE
- 3D structure modeling with MODELLER
- Crosslink application with support for terminal and non-terminal crosslinks
- Crosslink optimization using Monte Carlo methods
- Support for pre-mutated PDB workflow with additional crosslinks (for AGE crosslinks)
"""

# Copyright (c) 2024, ColBuilder Development Team
# Distributed under the terms of the Apache License 2.0

import asyncio
import os
import shutil
from pathlib import Path
from typing import Tuple, List, Optional, Dict, Any
from contextlib import asynccontextmanager
from dataclasses import asdict

import numpy as np
import pandas as pd

from colbuilder.core.utils.constants import (
    TEMP_FILE_SUFFIX,
    DISORIENTED_SUFFIX,
    PDB_EXTENSION,
    FASTA_EXTENSION,
)
from colbuilder.core.utils.data_structures import CrosslinkPair, OptimizationState
from colbuilder.core.utils.crosslinks import (
    CrosslinkOptimizer,
    extract_crosslinks_from_dataframe,
)
from colbuilder.core.utils.files import update_pdb_header, suppress_output
from colbuilder.core.sequence.alignment import align_sequences
from colbuilder.core.sequence.modeller import run_modeller
from colbuilder.core.sequence.mutate_crosslinks import apply_crosslinks
from colbuilder.core.utils.config import ColbuilderConfig
from colbuilder.core.utils.files import FileManager
from colbuilder.core.utils.exceptions import SequenceGenerationError, SystemError
from colbuilder.core.utils.logger import setup_logger

LOG = setup_logger(__name__)


class SequenceGenerator:
    """
    Manages the generation of collagen structures from sequences.
    
    This class coordinates sequence alignment, structure modeling, crosslink application
    and optimization, and file management throughout the generation process.
    """

    def __init__(
        self, 
        config: ColbuilderConfig, 
        file_manager: Optional[FileManager] = None
    ) -> None:
        """
        Initialize the sequence generator.

        Args:
            config: Configuration object for structure generation
            file_manager: Optional file manager instance
        """
        self.config = config
        self.file_manager = file_manager or FileManager(config)
        self._temp_dir: Optional[Path] = None
        self._crosslinks: List[CrosslinkPair] = []
        self._additional_crosslinks: List[CrosslinkPair] = []
        self._additional_crosslinks_1: List[CrosslinkPair] = []
        self._additional_crosslinks_2: List[CrosslinkPair] = []
        self._original_crosslinks: List[CrosslinkPair] = []
        self._state: Dict[str, Any] = {}

        if config.mutated_pdb:
            self.steps = 0
            if config.additional_1_type or config.additional_2_type:
                self.steps += 1  # Apply additional crosslinks
                self.steps += 1  # Optimize additional crosslinks
        else:
            self.steps = 4 if config.crosslink else 2

    @asynccontextmanager
    async def _manage_resources(self):
        """
        Context manager for resource lifecycle management.
        
        Handles working directory setup, temporary file management,
        and cleanup operations.
        
        Yields:
            None
            
        Raises:
            SequenceGenerationError: If working environment setup fails
        """
        try:
            original_dir = Path.cwd()
            working_dir = Path(self.config.working_directory)

            if not working_dir.exists():
                LOG.error(f"Working directory doesn't exist: {working_dir}")
                raise SequenceGenerationError(
                    "Working directory does not exist",
                    error_code="SEQ_ERR_004",
                    context={"working_directory": str(working_dir)},
                )

            temp_dir = self.file_manager.get_temp_path("sequence_gen", create_dir=True)
            os.chdir(temp_dir)

            yield

        except Exception as e:
            LOG.error(f"Error setting up directories: {e}")
            raise SequenceGenerationError(
                "Failed to set up working environment",
                error_code="SEQ_ERR_004",
                context={
                    "working_dir": str(working_dir),
                    "temp_dir": str(temp_dir) if "temp_dir" in locals() else None,
                    "current_dir": str(Path.cwd()),
                    "error": str(e),
                },
            )
        finally:
            try:
                os.chdir(original_dir)
            except Exception as e:
                LOG.warning(f"Failed to restore original directory {original_dir}: {e}")

            if not self.config.debug and self._temp_dir and self._temp_dir.exists():
                try:
                    shutil.rmtree(self._temp_dir)
                except Exception as e:
                    LOG.warning(f"Failed to clean up temporary directory {self._temp_dir}: {e}")

            self._temp_dir = None

    async def generate(self) -> Tuple[Optional[Path], Path]:
        """
        Generate collagen structure from sequence or apply additional mutations.
        
        Returns:
            Tuple[Optional[Path], Path]: MSA path (None for mutated PDB workflow) and final PDB
            
        Raises:
            SequenceGenerationError: If generation fails
        """
        LOG.debug(f"Working directory: {self.config.working_directory}")
        LOG.debug(f"Debug mode: {self.config.debug}")
        
        if self.config.mutated_pdb:
            LOG.info("Using pre-mutated PDB workflow")
            return await self._process_mutated_pdb()
        else:
            return await self._process_from_sequence()

    async def _process_mutated_pdb(self) -> Tuple[Optional[Path], Path]:
        """
        Process a pre-mutated PDB file with additional crosslinks.
        
        Returns:
            Tuple[Optional[Path], Path]: MSA path (None for mutated PDB workflow) and final PDB
            
        Raises:
            SequenceGenerationError: If processing fails
        """
        async with self._manage_resources():
            try:
                mutated_pdb_path = Path(self.config.mutated_pdb).resolve()
                
                if not mutated_pdb_path.exists():
                    raise SequenceGenerationError(
                        "Mutated PDB file not found",
                        error_code="SEQ_ERR_004",
                        context={"file_path": str(mutated_pdb_path)}
                    )
                
                LOG.info(f"Loading pre-mutated PDB: {mutated_pdb_path}")
                
                file_prefix = mutated_pdb_path.stem
                temp_pdb = Path.cwd() / mutated_pdb_path.name
                shutil.copy2(mutated_pdb_path, temp_pdb)
                
                LOG.info("Skipping original crosslink optimization (already in mutated PDB)")
                
                if self.config.additional_1_type or self.config.additional_2_type:
                    await self._load_additional_crosslinks()
                    
                    if not self._additional_crosslinks:
                        LOG.error("No additional crosslinks were loaded from the database")
                        LOG.error("Please check:")
                        LOG.error(f"  1. Crosslink type names are correct (case-sensitive)")
                        LOG.error(f"  2. Species '{self.config.species}' has these crosslink types")
                        LOG.error(f"  3. Crosslink combinations match the format in the database")
                        
                        raise SequenceGenerationError(
                            "Additional crosslinks not found in database",
                            error_code="SEQ_ERR_002",
                            context={
                                "species": self.config.species,
                                "additional_1_type": self.config.additional_1_type,
                                "additional_2_type": self.config.additional_2_type,
                                "additional_1_combination": self.config.additional_1_combination,
                                "additional_2_combination": self.config.additional_2_combination
                            }
                        )
                    
                    if self._additional_crosslinks:
                        current_step = 1
                        LOG.info(f"Step {current_step}/{self.steps} Applying additional crosslinks")
                        LOG.info(f"  Number of crosslinks to apply: {len(self._additional_crosslinks)}")
                        with suppress_output():
                            temp_pdb = await self._apply_additional_crosslinks(temp_pdb, file_prefix)
                        LOG.debug(f"Additional crosslinks applied - Output: {temp_pdb}")
                
                final_output = await self._finalize_output_mutated(temp_pdb, file_prefix)
                final_pdb = self.file_manager.copy_to_output(final_output)
                
                return None, final_pdb
                
            except Exception as e:
                LOG.error(f"Error processing mutated PDB: {str(e)}", exc_info=True)
                raise SequenceGenerationError(
                    "Error processing mutated PDB",
                    error_code="SEQ_ERR_001",
                    context={
                        "mutated_pdb": str(self.config.mutated_pdb),
                        "error": str(e)
                    }
                )

    async def _process_from_sequence(self) -> Tuple[Path, Path]:
        """
        Process collagen structure generation from sequence file.
        
        Returns:
            Tuple[Path, Path]: Paths to final MSA and PDB files
            
        Raises:
            SequenceGenerationError: If processing fails
        """
        async with self._manage_resources():
            try:
                fasta_file = Path(self.config.fasta_file).resolve()
                
                if not fasta_file.exists():
                    LOG.error(f"FASTA file not found: {fasta_file}")
                    raise SequenceGenerationError(
                        "FASTA file not found",
                        error_code="SEQ_ERR_004",
                        context={"file_path": str(fasta_file)}
                    )
                
                file_prefix = fasta_file.stem
                temp_fasta = Path.cwd() / fasta_file.name
                shutil.copy2(fasta_file, temp_fasta)
                
                msa_output, modeller_output = await self._run_alignment(temp_fasta)
                output_pdb = await self._run_modelling(modeller_output, file_prefix)
                
                if self.config.crosslink:
                    await self._load_crosslinks()
                    with suppress_output():
                        output_pdb = await self._apply_crosslinks(output_pdb, file_prefix)
                
                final_output = await self._finalize_output(output_pdb, file_prefix)
                
                final_msa = self.file_manager.copy_to_output(msa_output)
                final_pdb = self.file_manager.copy_to_output(final_output)
                
                return final_msa, final_pdb
                
            except Exception as e:
                LOG.error(f"Error during generation: {str(e)}", exc_info=True)
                raise

    async def _run_alignment(self, fasta_path: Path) -> Tuple[Path, Path]:
        """
        Run sequence alignment using MUSCLE.
        
        Args:
            fasta_path: Path to input FASTA file
            
        Returns:
            Tuple[Path, Path]: Paths to MSA output and MODELLER input files
            
        Raises:
            SequenceGenerationError: If alignment fails
        """
        try:
            LOG.info(f"Step 1/{self.steps} Performing sequence alignment")

            fasta_path = fasta_path.resolve()
            template_fasta = self.file_manager.find_file("sequence/template.fasta").resolve()
            template_pdb = self.file_manager.find_file("sequence/template.pdb").resolve()

            file_prefix = fasta_path.stem
            LOG.debug(f"Template FASTA: {template_fasta}")
            LOG.debug(f"Template PDB: {template_pdb}")

            for path, name in [
                (fasta_path, "Input FASTA"),
                (template_fasta, "Template FASTA"),
                (template_pdb, "Template PDB")
            ]:
                if not path.exists():
                    raise SequenceGenerationError(
                        f"{name} file not found",
                        error_code="SEQ_ERR_004",
                        context={f"{name.lower().replace(' ', '_')}_path": str(path)},
                    )

            with suppress_output():
                msa_output_path, modeller_output = align_sequences(
                    fasta_path, template_fasta, file_prefix, template_pdb
                )

            msa_output = Path(msa_output_path).resolve()
            modeller_out = Path(modeller_output).resolve()

            LOG.debug(f"Alignment completed. MSA output: {msa_output}")
            LOG.debug(f"Modeller output: {modeller_out}")

            for path, name in [
                (msa_output, "MSA output"),
                (modeller_out, "Modeller output")
            ]:
                if not path.exists():
                    raise SequenceGenerationError(
                        f"{name} file not created",
                        error_code="SEQ_ERR_001",
                        context={f"{name.lower().replace(' ', '_')}": str(path)},
                    )

            return msa_output, modeller_out

        except Exception as e:
            LOG.error(f"Sequence alignment failed: {str(e)}", exc_info=True)
            raise SequenceGenerationError(
                "Sequence alignment failed",
                error_code="SEQ_ERR_001",
                context={
                    "fasta_file": str(fasta_path),
                    "template_fasta": str(self.config.TEMPLATE_FASTA_PATH),
                    "current_dir": str(Path.cwd()),
                    "error": str(e),
                },
            )

    async def _run_modelling(self, modeller_output: Path, file_prefix: str) -> Path:
        """
        Run MODELLER for 3D structure generation.
        
        Args:
            modeller_output: Path to MODELLER input file
            file_prefix: Prefix for output files
            
        Returns:
            Path: Path to generated PDB file
            
        Raises:
            SequenceGenerationError: If modeling fails
        """
        try:
            LOG.info(f"Step 2/{self.steps} Generating structure with MODELLER")

            modeller_output = modeller_output.resolve()
            template_pdb = self.file_manager.find_file("sequence/template.pdb").resolve()
            restyp_lib = self.file_manager.find_file("sequence/modeller/restyp_mod.lib").resolve()
            top_heav_lib = self.file_manager.find_file("sequence/modeller/top_heav_mod.lib").resolve()
            par_mod_lib = self.file_manager.find_file("sequence/modeller/par_mod.lib").resolve()

            LOG.debug("Using libraries:")
            LOG.debug(f"  Template PDB: {template_pdb}")
            LOG.debug(f"  RESTYP: {restyp_lib}")
            LOG.debug(f"  TOP_HEAV: {top_heav_lib}")
            LOG.debug(f"  PAR_MOD: {par_mod_lib}")

            for path, name in [
                (modeller_output, "Modeller input"),
                (template_pdb, "Template PDB"),
                (restyp_lib, "RESTYP library"),
                (top_heav_lib, "TOP_HEAV library"),
                (par_mod_lib, "PAR_MOD library"),
            ]:
                if not path.exists():
                    raise SequenceGenerationError(
                        f"{name} file not found",
                        error_code="SEQ_ERR_004",
                        context={"missing_file": str(path)},
                    )

            with suppress_output():
                output_pdb = run_modeller(
                    aligned_file=str(modeller_output),
                    template_pdb=str(template_pdb),
                    output_prefix=file_prefix,
                    restyp_lib=str(restyp_lib),
                    top_heav_lib=str(top_heav_lib),
                    par_mod_lib=str(par_mod_lib),
                )

            output_path = Path(output_pdb).resolve()
            LOG.debug(f"MODELLER completed. Output PDB: {output_path}")

            if not output_path.exists():
                raise SequenceGenerationError(
                    "MODELLER output file not created",
                    error_code="SEQ_ERR_001",
                    context={"output_pdb": str(output_path)},
                )

            return output_path

        except Exception as e:
            LOG.error(f"MODELLER failed: {str(e)}", exc_info=True)
            raise SequenceGenerationError(
                "MODELLER structure generation failed",
                error_code="SEQ_ERR_001",
                context={
                    "modeller_output": str(modeller_output),
                    "template_pdb": str(self.config.TEMPLATE_PDB_PATH),
                    "current_dir": str(Path.cwd()),
                    "error": str(e),
                },
            )

    async def _load_crosslinks(self) -> None:
        """
        Load crosslink information from CSV file for standard workflow.
        
        Raises:
            SequenceGenerationError: If crosslink loading fails
        """
        if not self.config.crosslink:
            return

        try:
            crosslinks_file = self.file_manager.find_file("sequence/crosslinks.csv")
            crosslinks_df = pd.read_csv(crosslinks_file)
            species_crosslinks = crosslinks_df[
                crosslinks_df["species"] == self.config.species
            ]

            if species_crosslinks.empty:
                raise SequenceGenerationError(
                    f"No crosslinks found for species: {self.config.species}",
                    error_code="SEQ_ERR_002",
                    context={"species": self.config.species},
                )

            self._crosslinks = []

            if self.config.n_term_type:
                LOG.debug(f"Loading N-terminal crosslinks of type: {self.config.n_term_type}")
                n_crosslinks = extract_crosslinks_from_dataframe(
                    species_crosslinks,
                    "N",
                    self.config.n_term_type,
                    self.config.n_term_combination,
                )
                if n_crosslinks:
                    self._crosslinks.extend(n_crosslinks)
                    LOG.debug(f"Loaded {len(n_crosslinks)} N-terminal crosslinks")

            if self.config.c_term_type:
                LOG.debug(f"Loading C-terminal crosslinks of type: {self.config.c_term_type}")
                c_crosslinks = extract_crosslinks_from_dataframe(
                    species_crosslinks,
                    "C",
                    self.config.c_term_type,
                    self.config.c_term_combination,
                )
                if c_crosslinks:
                    self._crosslinks.extend(c_crosslinks)
                    LOG.debug(f"Loaded {len(c_crosslinks)} C-terminal crosslinks")

            if not self._crosslinks:
                LOG.warning("No crosslinks were loaded despite crosslinks being enabled")

        except SequenceGenerationError:
            raise
        except Exception as e:
            raise SequenceGenerationError(
                "Error loading crosslinks",
                original_error=e,
                error_code="SEQ_ERR_002",
                context={
                    "crosslinks_file": str(self.config.CROSSLINKS_FILE),
                    "species": self.config.species,
                    "n_term_type": self.config.n_term_type,
                    "c_term_type": self.config.c_term_type,
                },
            )

    async def _load_additional_crosslinks(self) -> None:
        """
        Load additional crosslink information for mutated PDB workflow.
        
        Raises:
            SequenceGenerationError: If crosslink loading fails
        """
        try:
            crosslinks_file = self.file_manager.find_file("sequence/crosslinks.csv")
            crosslinks_df = pd.read_csv(crosslinks_file)
            
            LOG.debug(f"Available species: {crosslinks_df['species'].unique()}")
            
            species_crosslinks = crosslinks_df[
                crosslinks_df["species"] == self.config.species
            ]
            
            if species_crosslinks.empty:
                raise SequenceGenerationError(
                    f"No crosslinks found for species: {self.config.species}",
                    error_code="SEQ_ERR_002",
                    context={"species": self.config.species}
                )
            
            LOG.debug(f"Available crosslink types for {self.config.species}: {species_crosslinks['type'].unique()}")
            
            self._additional_crosslinks = []
            self._additional_crosslinks_1 = []
            self._additional_crosslinks_2 = []
            
            if self.config.additional_1_type:
                LOG.debug(f"Loading additional crosslinks type 1: {self.config.additional_1_type}")
                LOG.debug(f"Looking for combination: {self.config.additional_1_combination}")
                
                terminal_1 = getattr(self.config, 'additional_1_terminal', None)
                
                if terminal_1:
                    crosslinks_1 = extract_crosslinks_from_dataframe(
                        species_crosslinks,
                        terminal_1,
                        self.config.additional_1_type,
                        self.config.additional_1_combination
                    )
                else:
                    crosslinks_1 = extract_crosslinks_from_dataframe(
                        species_crosslinks,
                        "N",
                        self.config.additional_1_type,
                        self.config.additional_1_combination
                    )
                    if not crosslinks_1:
                        crosslinks_1 = extract_crosslinks_from_dataframe(
                            species_crosslinks,
                            "C",
                            self.config.additional_1_type,
                            self.config.additional_1_combination
                        )
                    
                    if not crosslinks_1:
                        LOG.debug("Trying to find crosslink without terminal restriction")
                        type_crosslinks = species_crosslinks[
                            species_crosslinks['type'] == self.config.additional_1_type
                        ]
                        if not type_crosslinks.empty:
                            LOG.debug(f"Found {len(type_crosslinks)} crosslinks of type {self.config.additional_1_type}")
                            LOG.debug(f"Available combinations: {type_crosslinks['combination'].values}")
                            LOG.debug(f"Available terminals: {type_crosslinks['terminal'].values}")
                        else:
                            LOG.debug(f"No crosslinks found with type: {self.config.additional_1_type}")
                            LOG.debug(f"Available types in CSV: {species_crosslinks['type'].unique()}")
                
                if crosslinks_1:
                    self._additional_crosslinks_1.extend(crosslinks_1)
                    self._additional_crosslinks.extend(crosslinks_1)
                    LOG.info(f"Loaded {len(crosslinks_1)} additional crosslinks type 1: {self.config.additional_1_type}")
                else:
                    LOG.warning(f"No crosslinks found for type: {self.config.additional_1_type} with combination: {self.config.additional_1_combination}")
            
            if self.config.additional_2_type:
                LOG.debug(f"Loading additional crosslinks type 2: {self.config.additional_2_type}")
                LOG.debug(f"Looking for combination: {self.config.additional_2_combination}")
                terminal_2 = getattr(self.config, 'additional_2_terminal', None)
                
                if terminal_2:
                    crosslinks_2 = extract_crosslinks_from_dataframe(
                        species_crosslinks,
                        terminal_2,
                        self.config.additional_2_type,
                        self.config.additional_2_combination
                    )
                else:
                    crosslinks_2 = extract_crosslinks_from_dataframe(
                        species_crosslinks,
                        "N",
                        self.config.additional_2_type,
                        self.config.additional_2_combination
                    )
                    if not crosslinks_2:
                        crosslinks_2 = extract_crosslinks_from_dataframe(
                            species_crosslinks,
                            "C",
                            self.config.additional_2_type,
                            self.config.additional_2_combination
                        )
                        
                    if not crosslinks_2:
                        LOG.debug("Trying to find crosslink without terminal restriction")
                        type_crosslinks = species_crosslinks[
                            species_crosslinks['type'] == self.config.additional_2_type
                        ]
                        if not type_crosslinks.empty:
                            LOG.debug(f"Found {len(type_crosslinks)} crosslinks of type {self.config.additional_2_type}")
                            LOG.debug(f"Available combinations: {type_crosslinks['combination'].values}")
                            LOG.debug(f"Available terminals: {type_crosslinks['terminal'].values}")
                        else:
                            LOG.debug(f"No crosslinks found with type: {self.config.additional_2_type}")
                            LOG.debug(f"Available types in CSV: {species_crosslinks['type'].unique()}")
                
                if crosslinks_2:
                    self._additional_crosslinks_2.extend(crosslinks_2)
                    self._additional_crosslinks.extend(crosslinks_2)
                    LOG.info(f"Loaded {len(crosslinks_2)} additional crosslinks type 2: {self.config.additional_2_type}")
                else:
                    LOG.warning(f"No crosslinks found for type: {self.config.additional_2_type} with combination: {self.config.additional_2_combination}")
                    
            if not self._additional_crosslinks:
                LOG.warning("No additional crosslinks were loaded")
                    
        except Exception as e:
            raise SequenceGenerationError(
                "Error loading additional crosslinks",
                original_error=e,
                error_code="SEQ_ERR_002",
                context={
                    "species": self.config.species,
                    "additional_1_type": self.config.additional_1_type,
                    "additional_2_type": self.config.additional_2_type,
                    "additional_1_combination": self.config.additional_1_combination,
                    "additional_2_combination": self.config.additional_2_combination
                }
            )

    async def _apply_crosslinks(self, input_pdb: Path, file_prefix: str) -> Path:
        """
        Apply crosslinks to PDB structure for standard workflow.
        
        Args:
            input_pdb: Path to input PDB file
            file_prefix: Prefix for output files
            
        Returns:
            Path: Path to PDB file with crosslinks applied
            
        Raises:
            SequenceGenerationError: If crosslink application fails
        """
        if not self.config.crosslink or not self._crosslinks:
            return input_pdb

        try:
            n_suffix = f"N_{self.config.n_term_type}" if self.config.n_term_type else "N_NONE"
            c_suffix = f"C_{self.config.c_term_type}" if self.config.c_term_type else "C_NONE"
            output_pdb_crosslinked = f"{file_prefix}_{n_suffix}_{c_suffix}_temp.pdb"

            LOG.info(f"Step 3/{self.steps} Applying crosslinks to structure")

            n_crosslinks = [c for c in self._crosslinks if c.terminal_type == "N"]
            c_crosslinks = [c for c in self._crosslinks if c.terminal_type == "C"]

            def get_crosslink_data(crosslink: CrosslinkPair) -> Dict[str, str]:
                """Extract crosslink data into dictionary format."""
                data = {
                    "P1": crosslink.position1.position_str,
                    "R1": crosslink.position1.residue_type,
                    "A1": crosslink.position1.atom_name,
                    "P2": crosslink.position2.position_str,
                    "R2": crosslink.position2.residue_type,
                    "A2": crosslink.position2.atom_name,
                }
                if hasattr(crosslink, "position3") and crosslink.position3 is not None:
                    data.update({
                        "P3": crosslink.position3.position_str,
                        "R3": crosslink.position3.residue_type,
                        "A3": crosslink.position3.atom_name,
                    })
                return data

            n_crosslink_data = None if not n_crosslinks else get_crosslink_data(n_crosslinks[0])
            c_crosslink_data = None if not c_crosslinks else get_crosslink_data(c_crosslinks[0])

            LOG.debug(f"N-terminal crosslink data: {n_crosslink_data}")
            LOG.debug(f"C-terminal crosslink data: {c_crosslink_data}")

            result_path = Path(
                apply_crosslinks(
                    str(input_pdb),
                    output_pdb_crosslinked,
                    n_crosslink_data,
                    c_crosslink_data,
                    self.config,
                )
            )

            if not result_path.exists():
                raise SequenceGenerationError(
                    "Crosslink application failed to create output file",
                    error_code="SEQ_ERR_002",
                    context={"output_path": str(result_path)},
                )

            return result_path

        except AttributeError as e:
            LOG.error(f"Invalid crosslink data structure: {str(e)}")
            raise SequenceGenerationError(
                "Invalid crosslink data structure",
                original_error=e,
                error_code="SEQ_ERR_002",
                context={
                    "n_crosslinks": str(n_crosslinks) if "n_crosslinks" in locals() else None,
                    "c_crosslinks": str(c_crosslinks) if "c_crosslinks" in locals() else None,
                },
            )
        except SequenceGenerationError:
            raise
        except Exception as e:
            LOG.error(f"Error applying crosslinks: {str(e)}")
            raise SequenceGenerationError(
                "Error applying crosslinks",
                original_error=e,
                error_code="SEQ_ERR_002",
                context={
                    "input_pdb": str(input_pdb),
                    "file_prefix": file_prefix,
                    "n_terminal_type": self.config.n_term_type,
                    "c_terminal_type": self.config.c_term_type,
                },
            )

    async def _apply_additional_crosslinks(self, input_pdb: Path, file_prefix: str) -> Path:
        """
        Apply additional crosslinks to mutated PDB.
        
        Args:
            input_pdb: Path to input PDB file
            file_prefix: Prefix for output files
            
        Returns:
            Path: Path to PDB file with additional crosslinks applied
            
        Raises:
            SequenceGenerationError: If crosslink application fails
        """
        try:
            suffixes = []
            
            if self.config.additional_1_type:
                suffixes.append(f"ADD1_{self.config.additional_1_type}")
            if self.config.additional_2_type:
                suffixes.append(f"ADD2_{self.config.additional_2_type}")
                
            suffix = "+".join(suffixes) if suffixes else "NO_ADDITIONAL"
            output_pdb_crosslinked = f"{file_prefix}+{suffix}_temp.pdb"
            
            add_crosslinks_1 = self._additional_crosslinks_1
            add_crosslinks_2 = self._additional_crosslinks_2
            
            def get_crosslink_data(crosslink: CrosslinkPair) -> Dict[str, str]:
                """Extract crosslink data into dictionary format."""
                data = {
                    "P1": crosslink.position1.position_str,
                    "R1": crosslink.position1.residue_type,
                    "A1": crosslink.position1.atom_name,
                    "P2": crosslink.position2.position_str,
                    "R2": crosslink.position2.residue_type,
                    "A2": crosslink.position2.atom_name,
                }
                if hasattr(crosslink, "position3") and crosslink.position3 is not None:
                    data.update({
                        "P3": crosslink.position3.position_str,
                        "R3": crosslink.position3.residue_type,
                        "A3": crosslink.position3.atom_name,
                    })
                return data
            
            add_1_data = None if not add_crosslinks_1 else get_crosslink_data(add_crosslinks_1[0])
            add_2_data = None if not add_crosslinks_2 else get_crosslink_data(add_crosslinks_2[0])
            
            LOG.debug(f"Additional crosslink 1 data: {add_1_data}")
            LOG.debug(f"Additional crosslink 2 data: {add_2_data}")
            
            result_path = Path(
                apply_crosslinks(
                    str(input_pdb),
                    output_pdb_crosslinked,
                    add_1_data,
                    add_2_data,
                    self.config,
                )
            )
            
            if not result_path.exists():
                raise SequenceGenerationError(
                    "Additional crosslink application failed to create output file",
                    error_code="SEQ_ERR_002",
                    context={"output_path": str(result_path)}
                )
            
            return result_path
            
        except Exception as e:
            raise SequenceGenerationError(
                "Error applying additional crosslinks",
                original_error=e,
                error_code="SEQ_ERR_002",
                context={
                    "input_pdb": str(input_pdb),
                    "additional_1_type": self.config.additional_1_type,
                    "additional_2_type": self.config.additional_2_type
                }
            )

    async def _finalize_output(self, input_pdb: Path, file_prefix: str) -> Path:
        """
        Finalize output files and perform crosslink optimization for standard workflow.
        
        Args:
            input_pdb: Path to input PDB file
            file_prefix: Prefix for output files
            
        Returns:
            Path: Path to final optimized PDB file
            
        Raises:
            SequenceGenerationError: If finalization fails
        """
        try:
            if not self.config.crosslink:
                n_suffix = "N_NONE"
                c_suffix = "C_NONE"
                file_prefix = f"{file_prefix}_{n_suffix}_{c_suffix}"
            elif self.config.crosslink and self._crosslinks:
                n_suffix = f"N_{self.config.n_term_type}" if self.config.n_term_type else "N_NONE"
                c_suffix = f"C_{self.config.c_term_type}" if self.config.c_term_type else "C_NONE"
                file_prefix = f"{file_prefix}_{n_suffix}_{c_suffix}"

            dis_output_name = f"{file_prefix}_disoriented.pdb"
            dis_output_path = Path.cwd() / dis_output_name

            shutil.copy2(input_pdb, dis_output_path)
            LOG.debug(f"Created disoriented output: {dis_output_path}")

            final_name = f"{file_prefix}.pdb"
            final_output_path = Path.cwd() / final_name

            if self.config.crosslink and self._crosslinks:
                LOG.info(f"Step 4/{self.steps} Optimizing crosslink positions")
                chimera_scripts_dir = self.file_manager.find_file("chimera_scripts")

                crosslink_copies = getattr(self.config, 'crosslink_copies', ["D0", "D5"])
                LOG.info(f"Using crosslink translation pair: {crosslink_copies[0]} and {crosslink_copies[1]}")

                optimizer = CrosslinkOptimizer(
                    crosslink_pairs=self._crosslinks,
                    chimera_scripts_dir=chimera_scripts_dir,
                    crosslink_copies=crosslink_copies
                )

                final_distance, final_output = await optimizer.optimize(
                    input_pdb=dis_output_path, output_pdb=final_output_path
                )

                update_pdb_header(final_output_path, str(self.config.pdb_first_line))

                self._state.update({
                    "optimization_distance": final_distance,
                    "optimization_attempts": optimizer.state.attempt_number,
                    "crosslink_copies_used": crosslink_copies 
                })

                LOG.info(f"Crosslink optimization complete:")
                LOG.info(f"  Translation pair: {crosslink_copies[0]} - {crosslink_copies[1]}")
                LOG.info(f"  Optimization attempts: {optimizer.state.attempt_number}")

            else:
                LOG.info("No crosslinks optimization needed, preparing final output...")
                shutil.copy2(dis_output_path, final_output_path)
                update_pdb_header(final_output_path, str(self.config.pdb_first_line))

            try:
                if dis_output_path.exists() and not self.config.debug:
                    dis_output_path.unlink()
                    LOG.debug(f"Cleaned up temporary file: {dis_output_path}")
            except Exception as e:
                LOG.warning(f"Failed to delete temporary file {dis_output_path}: {str(e)}")

            if not final_output_path.exists():
                raise SequenceGenerationError(
                    "Final output file not created",
                    error_code="SEQ_ERR_004",
                    context={"final_output": str(final_output_path)},
                )

            return final_output_path

        except SequenceGenerationError:
            raise
        except Exception as e:
            LOG.error(f"Error in finalize_output: {str(e)}", exc_info=True)
            raise SequenceGenerationError(
                "Error finalizing output",
                error_code="SEQ_ERR_001",
                context={
                    "input_pdb": str(input_pdb),
                    "file_prefix": file_prefix,
                    "state": self._state,
                    "current_dir": str(Path.cwd()),
                    "working_dir": str(self.config.working_directory),
                },
            )

    async def _finalize_output_mutated(self, input_pdb: Path, file_prefix: str) -> Path:
        """
        Finalize output for mutated PDB workflow with optimization of additional crosslinks only.
        
        Args:
            input_pdb: Path to input PDB file
            file_prefix: Prefix for output files
            
        Returns:
            Path: Path to final optimized PDB file
            
        Raises:
            SequenceGenerationError: If finalization fails
        """
        try:
            original_prefix = file_prefix
            
            suffixes = []
            if self.config.additional_1_type:
                suffixes.append(f"ADD1_{self.config.additional_1_type}")
            if self.config.additional_2_type:
                suffixes.append(f"ADD2_{self.config.additional_2_type}")
            
            if suffixes:
                file_prefix = f"{original_prefix}+{'+'.join(suffixes)}"
            
            dis_output_name = f"{file_prefix}_disoriented.pdb"
            dis_output_path = Path.cwd() / dis_output_name
            
            shutil.copy2(input_pdb, dis_output_path)
            LOG.debug(f"Created disoriented output: {dis_output_path}")
            
            final_output_path = Path.cwd() / f"{file_prefix}.pdb"
            
            if self._additional_crosslinks:
                current_step = 2
                LOG.info(f"Step {current_step}/{self.steps} Optimizing additional crosslink positions")
                LOG.info(f"  Crosslinks to optimize: {len(self._additional_crosslinks)}")
                
                chimera_scripts_dir = self.file_manager.find_file("chimera_scripts")
                crosslink_copies = getattr(self.config, 'crosslink_copies', ["D0", "D5"])
                LOG.info(f"  Using translation pair: {crosslink_copies[0]} and {crosslink_copies[1]}")
                
                optimizer = CrosslinkOptimizer(
                    crosslink_pairs=self._additional_crosslinks,
                    chimera_scripts_dir=chimera_scripts_dir,
                    crosslink_copies=crosslink_copies
                )
                
                additional_distance, additional_output = await optimizer.optimize(
                    input_pdb=dis_output_path,
                    output_pdb=final_output_path
                )
                
                self._state["additional_optimization"] = {
                    "final_distance": additional_distance,
                    "attempts": optimizer.state.attempt_number,
                    "crosslink_copies": crosslink_copies
                }
                
                LOG.info(f"Additional crosslinks optimization complete:")
                LOG.info(f"  Final distance: {additional_distance:.2f} Å")
                LOG.info(f"  Attempts: {optimizer.state.attempt_number}")
            else:
                LOG.info("No additional crosslinks to optimize, preparing final output...")
                shutil.copy2(dis_output_path, final_output_path)
            
            update_pdb_header(final_output_path, str(self.config.pdb_first_line))
            
            if not self.config.debug:
                for temp_file in [dis_output_path]:
                    if temp_file.exists() and temp_file != final_output_path:
                        try:
                            temp_file.unlink()
                            LOG.debug(f"Cleaned up temporary file: {temp_file}")
                        except Exception as e:
                            LOG.warning(f"Failed to delete temporary file {temp_file}: {str(e)}")
            
            return final_output_path
            
        except Exception as e:
            raise SequenceGenerationError(
                "Error finalizing mutated PDB output",
                error_code="SEQ_ERR_001",
                context={
                    "input_pdb": str(input_pdb),
                    "error": str(e)
                }
            )

    @property
    def state(self) -> Dict[str, Any]:
        """
        Get current generation state.
        
        Returns:
            Dict[str, Any]: Copy of current state dictionary
        """
        return self._state.copy()

    @property
    def has_crosslinks(self) -> bool:
        """
        Check if structure has any crosslinks.
        
        Returns:
            bool: True if crosslinks exist
        """
        return bool(self._crosslinks) or bool(self._additional_crosslinks)

    def get_crosslink_info(self) -> Dict[str, List[Dict[str, Any]]]:
        """
        Get information about current crosslinks.
        
        Returns:
            Dict[str, List[Dict[str, Any]]]: Dictionary containing crosslink information
                - 'additional': List of all additional crosslinks
                - 'additional_1': List of additional crosslinks type 1
                - 'additional_2': List of additional crosslinks type 2
                - 'original': List of original crosslinks (if not in mutated PDB mode)
                - 'all': Combined list of all crosslinks
        """
        info = {
            "additional": [asdict(crosslink) for crosslink in self._additional_crosslinks],
            "additional_1": [asdict(crosslink) for crosslink in self._additional_crosslinks_1],
            "additional_2": [asdict(crosslink) for crosslink in self._additional_crosslinks_2]
        }
        
        if not self.config.mutated_pdb:
            info["original"] = [asdict(crosslink) for crosslink in self._crosslinks]
            info["all"] = [asdict(crosslink) for crosslink in (self._crosslinks + self._additional_crosslinks)]
        else:
            info["all"] = info["additional"]
            
        return info
