"""
This module serves as the main entry point for the sequence generation process in the ColBuilder pipeline.

It provides a high-level function, `build_sequence`, which orchestrates the generation of collagen
structures from input configurations. This includes sequence alignment, structure modeling, and
crosslink application, leveraging the `SequenceGenerator` class.

Key Features:
--------------
1. **Sequence Generation**:
   - Align input sequences with templates.
   - Generate multiple sequence alignments (MSA) and 3D protein structures.

2. **Error Handling**:
   - Raise detailed exceptions for sequence generation or system operation failures.

3. **Integration with ColBuilder**:
   - Designed to work seamlessly within the ColBuilder pipeline for collagen modeling.

Usage:
------
This module is intended to be used as part of the ColBuilder pipeline. The main function,
`build_sequence`, takes a configuration object as input and outputs the paths to the generated
MSA and PDB files.

Example:
--------
```python
from colbuilder.core.sequence.main_sequence import build_sequence
from colbuilder.core.utils.config import ColbuilderConfig

# Define configuration
config = ColbuilderConfig(
    working_directory="/path/to/working_dir",
    fasta_file="/path/to/input.fasta",
    crosslink=True,
    debug=False
)

# Generate sequence
msa_output, pdb_output = await build_sequence(config)

print(f"MSA file saved to: {msa_output}")
print(f"PDB file saved to: {pdb_output}")
```
"""

# Copyright (c) 2024, ColBuilder Development Team
# Distributed under the terms of the Apache License 2.0

from pathlib import Path
from typing import Optional, Tuple
from colorama import Fore, Style
from colbuilder.core.utils.config import ColbuilderConfig
from .sequence_generator import SequenceGenerator

from ..utils.logger import setup_logger

LOG = setup_logger(__name__)


# In colbuilder/core/sequence/main_sequence.py

async def build_sequence(config: ColbuilderConfig) -> Tuple[Optional[Path], Path]:
    """
    Build collagen sequence with structure generation and optional crosslink optimization.
    
    Args:
        config: Configuration object containing all sequence generation parameters
        
    Returns:
        Tuple[Optional[Path], Path]: MSA path (None for mutated PDB workflow) and final PDB structure
    """
    try:
        # Log sequence generation parameters
        LOG.info("Starting sequence generation with parameters:")
        LOG.info(f"  Species: {config.species}")
        LOG.info(f"  Crosslinks enabled: {config.crosslink}")
        
        if config.crosslink:
            LOG.info(f"  N-terminal type: {config.n_term_type or 'None'}")
            LOG.info(f"  C-terminal type: {config.c_term_type or 'None'}")
            LOG.info(f"  Crosslink optimization translations: {getattr(config, 'crosslink_copies', ['D0', 'D5'])}")
        
        # Initialize the sequence generator
        generator = SequenceGenerator(config)
        
        # Generate the structure
        final_msa, final_pdb = await generator.generate()
        
        # Log results
        LOG.info("Sequence generation completed successfully")
        if final_msa is None:
            LOG.info("  MSA output: None (mutated PDB workflow)")
        else:
            LOG.info(f"  MSA output: {final_msa}")
        LOG.info(f"  PDB output: {final_pdb}")
        
        # If crosslinks were optimized, log the optimization results
        if config.crosslink and generator.has_crosslinks:
            state = generator.state
            if 'optimization_distance' in state:
                LOG.info(f"Crosslink optimization results:")
                LOG.info(f"  Attempts: {state['optimization_attempts']}")
                if 'crosslink_copies_used' in state:
                    LOG.info(f"  Translations used: {state['crosslink_copies_used']}")
        
        return final_msa, final_pdb
        
    except Exception as e:
        LOG.error(f"Sequence generation failed: {str(e)}")
        raise


__all__ = ["build_sequence"]
