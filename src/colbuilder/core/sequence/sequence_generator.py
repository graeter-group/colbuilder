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

4. **Crosslink Optimization**:
   - Optimize crosslink positions using Monte Carlo method and Chimera scripts.
   - Minimize distances between crosslinked residues to satisfy geometric constraints.

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

# Load configuration
config = ColbuilderConfig(
    working_directory="/path/to/working_dir",
    fasta_file="/path/to/input.fasta",
    crosslink=True,
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
        self._state: Dict[str, Any] = {}
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
        """Generate collagen structure from sequence."""
        LOG.debug(f"Working directory: {self.config.working_directory}")
        LOG.debug(f"Input FASTA: {self.config.fasta_file}")
        LOG.debug(f"Debug mode: {self.config.debug}")

        original_dir = Path.cwd().resolve()

        async with self._manage_resources():
            try:
                fasta_file = Path(self.config.fasta_file).resolve()

                if not fasta_file.exists():
                    LOG.error(f"FASTA file not found: {fasta_file}")
                    raise SequenceGenerationError(
                        "FASTA file not found",
                        error_code="SEQ_ERR_004",
                        context={"file_path": str(fasta_file)},
                    )

                file_prefix = fasta_file.stem
                temp_fasta = Path.cwd() / fasta_file.name
                shutil.copy2(fasta_file, temp_fasta)

                msa_output, modeller_output = await self._run_alignment(temp_fasta)
                LOG.debug(
                    f"Alignment complete - MSA: {msa_output}, Modeller: {modeller_output}"
                )

                output_pdb = await self._run_modelling(modeller_output, file_prefix)
                LOG.debug(f"MODELLER complete - Output: {output_pdb}")

                if self.config.crosslink:
                    await self._load_crosslinks()
                    with suppress_output():
                        output_pdb = await self._apply_crosslinks(
                            output_pdb, file_prefix
                        )

                final_output = await self._finalize_output(output_pdb, file_prefix)

                final_msa = self.file_manager.copy_to_output(msa_output)
                final_pdb = self.file_manager.copy_to_output(final_output)

                return final_msa, final_pdb

            except (SequenceGenerationError, SystemError):
                LOG.error("Known error occurred during generation", exc_info=True)
                raise
            except Exception as e:
                LOG.error(
                    f"Unexpected error during generation: {str(e)}", exc_info=True
                )
                raise SequenceGenerationError(
                    "Unexpected error during sequence generation",
                    error_code="SEQ_ERR_001",
                    context={
                        "state": self._state,
                        "error": str(e),
                        "error_type": type(e).__name__,
                    },
                )

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

                optimizer = CrosslinkOptimizer(
                    crosslink_pairs=self._crosslinks,
                    chimera_scripts_dir=chimera_scripts_dir,
                )

                final_distance, final_output = await optimizer.optimize(
                    input_pdb=dis_output_path, output_pdb=final_output_path
                )

                update_pdb_header(final_output_path, str(self.config.pdb_first_line))

                self._state.update(
                    {
                        "optimization_distance": final_distance,
                        "optimization_attempts": optimizer.state.attempt_number,
                    }
                )

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

    @property
    def state(self) -> Dict[str, Any]:
        """Get current generation state."""
        return self._state.copy()

    @property
    def has_crosslinks(self) -> bool:
        """Check if structure has crosslinks."""
        return bool(self._crosslinks)

    def get_crosslink_info(self) -> List[Dict[str, Any]]:
        """Get information about current crosslinks."""
        return [asdict(crosslink) for crosslink in self._crosslinks]
