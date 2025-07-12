"""
This module provides tools for generating collagen structures from input sequences.

The `SequenceGenerator` class orchestrates the entire process of generating collagen structures,
including sequence alignment, structure modeling, crosslink application and optimization.

Key Features:
--------------
1. **Sequence Alignment**:
   - Align input sequences with template sequences using MUSCLE.
   - Generate multiple sequence alignments (MSA) and prepare input for structure modeling.

2. **Structure Modeling**:
   - Use MODELLER to generate 3D structures from aligned sequences.
   - Validate and handle required input files for modeling.

3. **Crosslink Application**:
   - Apply crosslinks to the generated structure based on user-defined configurations.
   - Support for N-terminal and C-terminal crosslinks with flexible residue and atom specifications.
   - Support for additional crosslinks that can be placed anywhere in the molecule.

4. **Crosslink Optimization**:
   - Optimize crosslink positions using Monte Carlo method and Chimera scripts.
   - Minimize distances between crosslinked residues to satisfy geometric constraints.
   - For mutated PDB workflow: only optimize newly added crosslinks, not existing ones.

5. **Mutated PDB Workflow**:
   - Load pre-mutated PDB files that already contain optimized crosslinks.
   - Apply only additional crosslinks without re-optimizing existing ones.
   - Optimize only the newly added crosslinks.

Usage:
------
This module is designed to be used as part of a pipeline for generating collagen structures.
The main entry point is the `generate` method of the `SequenceGenerator` class, which performs
all steps in sequence and outputs the final structure.

Example:
--------
```python
from colbuilder.core.sequence.sequence_generator import SequenceGenerator
from colbuilder.core.utils.config import ColbuilderConfig

# Normal workflow - generate from sequence
config = ColbuilderConfig(
    working_directory="/path/to/working_dir",
    fasta_file="/path/to/input.fasta",
    crosslink=True,
    debug=False
)

# Mutated PDB workflow - add crosslinks to existing structure
config = ColbuilderConfig(
    working_directory="/path/to/working_dir",
    mutated_pdb="/path/to/existing.pdb",
    additional_1_type="Glucosepane",
    additional_1_combination="523.A - 286.C",
    debug=False
)

# Initialize the sequence generator
generator = SequenceGenerator(config)

# Generate the structure
final_msa, final_pdb = await generator.generate()

print(f"Final MSA file: {final_msa}")
print(f"Final PDB file: {final_pdb}")
```
"""

# Copyright (c) 2024, ColBuilder Development Team
# Distributed under the terms of the Apache License 2.0

import asyncio
import os
import tempfile
import shutil
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Tuple, List, Optional, Dict, Any
from contextlib import asynccontextmanager
from dataclasses import asdict

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
    Manages the generation of collagen structure from sequences.

    This class coordinates the entire structure from sequence generation process, including:
    - Sequence alignment
    - Structure modeling
    - Crosslink application and optimization
    - File management and cleanup

    Attributes:
        config (ColbuilderConfig): Configuration for sequence generation
        _temp_dir (Optional[Path]): Path to temporary working directory
        _crosslinks (List[CrosslinkPair]): List of crosslinks to apply
        _state (Dict[str, Any]): Current generation state
    """

    def __init__(
        self, config: ColbuilderConfig, file_manager: Optional[FileManager] = None
    ):
        """
        Initialize the structure from sequence generator.

        Args:
            config: Configuration object for structure generation
        """
        self.config = config
        self.file_manager = file_manager or FileManager(config)
        self._temp_dir: Optional[Path] = None
        self._crosslinks: List[CrosslinkPair] = []
        self._additional_crosslinks: List[CrosslinkPair] = []
        self._original_crosslinks: List[CrosslinkPair] = []
        self._state: Dict[str, Any] = {}

        # Adjust steps based on workflow
        if config.mutated_pdb:
            # Mutated PDB workflow: apply additional crosslinks + optimize only additional
            self.steps = 0
            if config.additional_1_type or config.additional_2_type:
                self.steps += 1  # Apply additional crosslinks
                self.steps += 1  # Optimize additional crosslinks
        else:
            # Normal workflow: align + model + crosslink + optimize
            self.steps = 4 if config.crosslink else 2

    @asynccontextmanager
    async def _manage_resources(self):
        """Context manager for resource lifecycle management."""
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
                    LOG.warning(
                        f"Failed to clean up temporary directory {self._temp_dir}: {e}"
                    )

            self._temp_dir = None

    async def _run_alignment(self, fasta_path: Path) -> Tuple[Path, Path]:
        """Run sequence alignment."""
        try:
            LOG.info(f"Step 1/{self.steps} Performing sequence alignment")

            fasta_path = fasta_path.resolve()
            template_fasta = self.file_manager.find_file(
                "sequence/template.fasta"
            ).resolve()
            template_pdb = self.file_manager.find_file(
                "sequence/template.pdb"
            ).resolve()

            file_prefix = fasta_path.stem
            LOG.debug(f"Template FASTA: {template_fasta}")
            LOG.debug(f"Template PDB: {template_pdb}")

            if not fasta_path.exists():
                raise SequenceGenerationError(
                    "Input FASTA file not found",
                    error_code="SEQ_ERR_004",
                    context={"fasta_path": str(fasta_path)},
                )

            if not template_fasta.exists():
                raise SequenceGenerationError(
                    "Template FASTA file not found",
                    error_code="SEQ_ERR_004",
                    context={"template_fasta": str(template_fasta)},
                )

            if not template_pdb.exists():
                raise SequenceGenerationError(
                    "Template PDB file not found",
                    error_code="SEQ_ERR_004",
                    context={"template_pdb": str(template_pdb)},
                )

            with suppress_output():
                msa_output_path, modeller_output = align_sequences(
                    fasta_path, template_fasta, file_prefix, template_pdb
                )

            msa_output = Path(msa_output_path).resolve()
            modeller_out = Path(modeller_output).resolve()

            LOG.debug(f"Alignment completed. MSA output: {msa_output}")
            LOG.debug(f"Modeller output: {modeller_out}")

            if not msa_output.exists():
                raise SequenceGenerationError(
                    "MSA output file not created",
                    error_code="SEQ_ERR_001",
                    context={"msa_output": str(msa_output)},
                )

            if not modeller_out.exists():
                raise SequenceGenerationError(
                    "Modeller output file not created",
                    error_code="SEQ_ERR_001",
                    context={"modeller_output": str(modeller_out)},
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
        """Run MODELLER for structure generation."""
        try:
            LOG.info(f"Step 2/{self.steps} Generating structure with MODELLER")

            modeller_output = modeller_output.resolve()
            template_pdb = self.file_manager.find_file(
                "sequence/template.pdb"
            ).resolve()
            restyp_lib = self.file_manager.find_file(
                "sequence/modeller/restyp_mod.lib"
            ).resolve()
            top_heav_lib = self.file_manager.find_file(
                "sequence/modeller/top_heav_mod.lib"
            ).resolve()
            par_mod_lib = self.file_manager.find_file(
                "sequence/modeller/par_mod.lib"
            ).resolve()

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

    async def generate(self) -> Tuple[Path, Path]:
        """Generate collagen structure from sequence or apply additional mutations."""
        LOG.debug(f"Working directory: {self.config.working_directory}")
        LOG.debug(f"Debug mode: {self.config.debug}")
        
        # Check if we're using mutated PDB workflow
        if self.config.mutated_pdb:
            LOG.info("Using pre-mutated PDB workflow")
            return await self._process_mutated_pdb()
        else:
            # Original workflow
            return await self._process_from_sequence()
    
    async def _process_mutated_pdb(self) -> Tuple[Path, Path]:
        """Process a pre-mutated PDB file with additional crosslinks."""
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
                
                # For mutated PDB workflow, we assume original crosslinks are already optimized
                # We do NOT load or optimize them again
                LOG.info("Skipping original crosslink optimization (already in mutated PDB)")
                
                # Apply additional crosslinks
                if self.config.additional_1_type or self.config.additional_2_type:
                    await self._load_additional_crosslinks()
                    if self._additional_crosslinks:
                        current_step = 1
                        LOG.info(f"Step {current_step}/{self.steps} Applying additional crosslinks")
                        with suppress_output():
                            temp_pdb = await self._apply_additional_crosslinks(temp_pdb, file_prefix)
                        LOG.debug(f"Additional crosslinks applied - Output: {temp_pdb}")
                
                final_output = await self._finalize_output_mutated(temp_pdb, file_prefix)
                
                # For mutated PDB workflow, we don't have an MSA file
                # Create a dummy MSA path or use the original if it exists
                msa_path = Path.cwd() / f"{file_prefix}.msa"
                if not msa_path.exists():
                    # Create a simple MSA file indicating this was from mutated PDB
                    with open(msa_path, 'w') as f:
                        f.write(f"# MSA file for mutated PDB workflow\n")
                        f.write(f"# Source: {mutated_pdb_path}\n")
                
                final_msa = self.file_manager.copy_to_output(msa_path)
                final_pdb = self.file_manager.copy_to_output(final_output)
                
                return final_msa, final_pdb
                
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
        """Original workflow: process from sequence."""
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

    async def _apply_crosslinks(self, input_pdb: Path, file_prefix: str) -> Path:
        """
        Apply crosslinks to the PDB structure if enabled.

        Args:
            input_pdb: Path to input PDB file
            file_prefix: Prefix for output files

        Returns:
            Path to PDB file with crosslinks applied (or original if disabled)

        Raises:
            SequenceGenerationError: If crosslink application fails
        """
        if not self.config.crosslink or not self._crosslinks:
            return input_pdb

        try:
            n_suffix = (
                f"N_{self.config.n_term_type}" if self.config.n_term_type else "N_NONE"
            )
            c_suffix = (
                f"C_{self.config.c_term_type}" if self.config.c_term_type else "C_NONE"
            )
            output_pdb_crosslinked = f"{file_prefix}_{n_suffix}_{c_suffix}_temp.pdb"

            LOG.info(f"Step 3/{self.steps} Applying crosslinks to structure")

            n_crosslinks = [c for c in self._crosslinks if c.terminal_type == "N"]
            c_crosslinks = [c for c in self._crosslinks if c.terminal_type == "C"]

            def get_crosslink_data(crosslink):
                data = {
                    "P1": crosslink.position1.position_str,
                    "R1": crosslink.position1.residue_type,
                    "A1": crosslink.position1.atom_name,
                    "P2": crosslink.position2.position_str,
                    "R2": crosslink.position2.residue_type,
                    "A2": crosslink.position2.atom_name,
                }
                if hasattr(crosslink, "position3") and crosslink.position3 is not None:
                    data.update(
                        {
                            "P3": crosslink.position3.position_str,
                            "R3": crosslink.position3.residue_type,
                            "A3": crosslink.position3.atom_name,
                        }
                    )
                return data

            n_crosslink_data = (
                None if not n_crosslinks else get_crosslink_data(n_crosslinks[0])
            )
            c_crosslink_data = (
                None if not c_crosslinks else get_crosslink_data(c_crosslinks[0])
            )

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
                    "n_crosslinks": (
                        str(n_crosslinks) if "n_crosslinks" in locals() else None
                    ),
                    "c_crosslinks": (
                        str(c_crosslinks) if "c_crosslinks" in locals() else None
                    ),
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

    async def _load_crosslinks(self) -> None:
        """
        Load crosslink information from configuration.
        Only called when crosslinks are enabled.

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
                LOG.debug(
                    f"Loading N-terminal crosslinks of type: {self.config.n_term_type}"
                )
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
                LOG.debug(
                    f"Loading C-terminal crosslinks of type: {self.config.c_term_type}"
                )
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
                LOG.warning(
                    "No crosslinks were loaded despite crosslinks being enabled"
                )

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
        """Load additional crosslink information for mutated PDB workflow."""
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
                    context={"species": self.config.species}
                )
            
            self._additional_crosslinks = []
            
            # Load additional crosslinks type 1
            if self.config.additional_1_type:
                LOG.debug(f"Loading additional crosslinks type 1: {self.config.additional_1_type}")
                # Check if there's a specific terminal indicator in the config
                # If not provided, we'll search both N and C terminals
                terminal_1 = getattr(self.config, 'additional_1_terminal', None)
                
                if terminal_1:
                    crosslinks_1 = extract_crosslinks_from_dataframe(
                        species_crosslinks,
                        terminal_1,
                        self.config.additional_1_type,
                        self.config.additional_1_combination
                    )
                else:
                    # Try both terminals and use whichever has the crosslink
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
                
                if crosslinks_1:
                    # Mark these as additional crosslinks
                    for cl in crosslinks_1:
                        cl.crosslink_type = self.config.additional_1_type
                    self._additional_crosslinks.extend(crosslinks_1)
                    LOG.debug(f"Loaded {len(crosslinks_1)} additional crosslinks type 1")
            
            # Load additional crosslinks type 2
            if self.config.additional_2_type:
                LOG.debug(f"Loading additional crosslinks type 2: {self.config.additional_2_type}")
                terminal_2 = getattr(self.config, 'additional_2_terminal', None)
                
                if terminal_2:
                    crosslinks_2 = extract_crosslinks_from_dataframe(
                        species_crosslinks,
                        terminal_2,
                        self.config.additional_2_type,
                        self.config.additional_2_combination
                    )
                else:
                    # Try both terminals
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
                
                if crosslinks_2:
                    # Mark these as additional crosslinks
                    for cl in crosslinks_2:
                        cl.crosslink_type = self.config.additional_2_type
                    self._additional_crosslinks.extend(crosslinks_2)
                    LOG.debug(f"Loaded {len(crosslinks_2)} additional crosslinks type 2")
                    
        except Exception as e:
            raise SequenceGenerationError(
                "Error loading additional crosslinks",
                original_error=e,
                error_code="SEQ_ERR_002",
                context={
                    "species": self.config.species,
                    "additional_1_type": self.config.additional_1_type,
                    "additional_2_type": self.config.additional_2_type
                }
            )
    
    async def _apply_additional_crosslinks(self, input_pdb: Path, file_prefix: str) -> Path:
        """Apply additional crosslinks to mutated PDB."""
        try:
            # Build suffix based only on additional crosslinks
            suffixes = []
            
            if self.config.additional_1_type:
                suffixes.append(f"ADD1_{self.config.additional_1_type}")
            if self.config.additional_2_type:
                suffixes.append(f"ADD2_{self.config.additional_2_type}")
                
            suffix = "_".join(suffixes) if suffixes else "NO_ADDITIONAL"
            output_pdb_crosslinked = f"{file_prefix}_{suffix}_temp.pdb"
            
            # Separate additional crosslinks by type
            add_crosslinks_1 = [c for c in self._additional_crosslinks 
                               if hasattr(c, 'crosslink_type') and c.crosslink_type == self.config.additional_1_type]
            add_crosslinks_2 = [c for c in self._additional_crosslinks 
                               if hasattr(c, 'crosslink_type') and c.crosslink_type == self.config.additional_2_type]
            
            # If crosslink_type is not set, fall back to terminal_type
            if not add_crosslinks_1 and not add_crosslinks_2:
                add_crosslinks_1 = [c for c in self._additional_crosslinks if c.terminal_type == "N"]
                add_crosslinks_2 = [c for c in self._additional_crosslinks if c.terminal_type == "C"]
            
            # Helper function to extract crosslink data
            def get_crosslink_data(crosslink):
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
            
            # Prepare crosslink data
            add_1_data = None if not add_crosslinks_1 else get_crosslink_data(add_crosslinks_1[0])
            add_2_data = None if not add_crosslinks_2 else get_crosslink_data(add_crosslinks_2[0])
            
            LOG.debug(f"Additional crosslink 1 data: {add_1_data}")
            LOG.debug(f"Additional crosslink 2 data: {add_2_data}")
            
            # Apply crosslinks using the existing function
            # Note: We're passing additional crosslinks in place of N/C terminal crosslinks
            result_path = Path(
                apply_crosslinks(
                    str(input_pdb),
                    output_pdb_crosslinked,
                    add_1_data,  # Using position 1 for first additional crosslink
                    add_2_data,  # Using position 2 for second additional crosslink
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
        """Finalize output files and perform crosslink optimization if needed."""
        try:
            if not self.config.crosslink:
                n_suffix = f"N_NONE"
                c_suffix = f"C_NONE"
                file_prefix = f"{file_prefix}_{n_suffix}_{c_suffix}"
            elif self.config.crosslink and self._crosslinks:
                n_suffix = (
                    f"N_{self.config.n_term_type}"
                    if self.config.n_term_type
                    else "N_NONE"
                )
                c_suffix = (
                    f"C_{self.config.c_term_type}"
                    if self.config.c_term_type
                    else "C_NONE"
                )
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

                self._state.update(
                    {
                        "optimization_distance": final_distance,
                        "optimization_attempts": optimizer.state.attempt_number,
                        "crosslink_copies_used": crosslink_copies 
                    }
                )

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
                LOG.warning(
                    f"Failed to delete temporary file {dis_output_path}: {str(e)}"
                )

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
        """Finalize output for mutated PDB workflow with separate optimizations."""
        try:
            # Extract original crosslinks from filename if possible
            original_prefix = file_prefix
            
            # Build suffix for additional crosslinks only
            suffixes = []
            if self.config.additional_1_type:
                suffixes.append(f"ADD1_{self.config.additional_1_type}")
            if self.config.additional_2_type:
                suffixes.append(f"ADD2_{self.config.additional_2_type}")
            
            # If we have additional crosslinks, append them to the original prefix
            if suffixes:
                file_prefix = f"{original_prefix}_{'+'.join(suffixes)}"
            
            dis_output_name = f"{file_prefix}_disoriented.pdb"
            dis_output_path = Path.cwd() / dis_output_name
            
            shutil.copy2(input_pdb, dis_output_path)
            LOG.debug(f"Created disoriented output: {dis_output_path}")
            
            final_output_path = Path.cwd() / f"{file_prefix}.pdb"
            
            # Only optimize additional crosslinks (if they exist)
            if self._additional_crosslinks:
                current_step = 2  # Step 1 was applying crosslinks
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
                # No additional crosslinks to optimize
                LOG.info("No additional crosslinks to optimize, preparing final output...")
                shutil.copy2(dis_output_path, final_output_path)
            
            # Update PDB header
            update_pdb_header(final_output_path, str(self.config.pdb_first_line))
            
            # Cleanup
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
        """Get current generation state."""
        return self._state.copy()

    @property
    def has_crosslinks(self) -> bool:
        """Check if structure has crosslinks."""
        return bool(self._crosslinks) or bool(self._additional_crosslinks)

    def get_crosslink_info(self) -> Dict[str, List[Dict[str, Any]]]:
        """Get information about current crosslinks."""
        info = {
            "additional": [asdict(crosslink) for crosslink in self._additional_crosslinks]
        }
        
        # Only include original crosslinks if we're not in mutated PDB mode
        if not self.config.mutated_pdb:
            info["original"] = [asdict(crosslink) for crosslink in self._crosslinks]
            info["all"] = [asdict(crosslink) for crosslink in (self._crosslinks + self._additional_crosslinks)]
        else:
            info["all"] = info["additional"]
            
        return info