"""
This module provides a wrapper for MODELLER®, a program for comparative protein structure modeling.

MODELLER® is a trademark of the Regents of the University of California, developed by Andrej Sali
and colleagues at the University of California, San Francisco. For more information, visit:
https://salilab.org/modeller/

The `ModellerWrapper` class simplifies the integration of MODELLER into the ColBuilder project,
allowing users to generate protein structure models from aligned sequences and template structures.
It supports custom residue type, topology, and parameter libraries for advanced modeling.

Key Features:
--------------
1. **Protein Structure Modeling**:
   - Generate 3D protein structures using aligned sequences and template PDB files.
   - Support for custom residue types, topology, and parameter libraries.

2. **Custom Segment Renaming**:
   - Rename chain segments and renumber residues for consistency in output models.

3. **Integration with ColBuilder**:
   - Designed to work seamlessly within the ColBuilder pipeline for collagen modeling.

Usage:
------
This module is designed to be used as part of a pipeline for generating collagen structures.
The main entry point is the `run_modeller` function or the `ModellerWrapper` class, which
handles the entire modeling process and outputs the final PDB file.

Example:
--------
```python
from colbuilder.core.sequence.modeller import run_modeller

# Define input files and parameters
aligned_file = "aligned_sequences.ali"
template_pdb = "template_structure.pdb"
output_prefix = "model_output"
restyp_lib = "/path/to/restyp.lib"
top_heav_lib = "/path/to/top_heav.lib"
par_mod_lib = "/path/to/par.lib"

# Run MODELLER
output_pdb = run_modeller(
    aligned_file=aligned_file,
    template_pdb=template_pdb,
    output_prefix=output_prefix,
    restyp_lib=restyp_lib,
    top_heav_lib=top_heav_lib,
    par_mod_lib=par_mod_lib
)

print(f"Generated structure saved to: {output_pdb}")
```
"""

# Copyright (c) 2024, ColBuilder Development Team
# Distributed under the terms of the Apache License 2.0

import os
from typing import Optional
from modeller import Environ
from modeller.automodel import AutoModel

from colbuilder.core.utils.logger import setup_logger
from colbuilder.core.utils.dec import timeit

LOG = setup_logger(__name__)


class ModellerWrapper:
    """
    Attributes:
        aligned_file (str): Path to the aligned sequence file.
        template_pdb (str): Path to the template PDB file.
        output_prefix (str): Prefix for output files.
        restyp_lib (str): Path to the residue type library file.
        top_heav_lib (str): Path to the topology library file.
        par_mod_lib (str): Path to the parameter library file.
        output_pdb (Optional[str]): Path to the output PDB file.
    """

    def __init__(
        self,
        aligned_file: str,
        template_pdb: str,
        output_prefix: str,
        restyp_lib: str,
        top_heav_lib: str,
        par_mod_lib: str,
    ):
        """
        Initialize the ModellerWrapper.

        Args:
            aligned_file (str): Path to the aligned sequence file.
            template_pdb (str): Path to the template PDB file.
            output_prefix (str): Prefix for output files.
            restyp_lib (str): Path to custom residue type library file.
            top_heav_lib (str): Path to custom topology library file.
            par_mod_lib (str): Path to custom parameter library file.
        """
        self.aligned_file = aligned_file
        self.template_pdb = template_pdb
        self.output_prefix = output_prefix
        self.restyp_lib = restyp_lib
        self.top_heav_lib = top_heav_lib
        self.par_mod_lib = par_mod_lib
        self.output_pdb: Optional[str] = None

    @property
    def output_pdb(self) -> Optional[str]:
        return self._output_pdb

    @output_pdb.setter
    def output_pdb(self, value: Optional[str]):
        self._output_pdb = value

    def run_modeller(self) -> None:
        """
        Run MODELLER to generate a protein structure model.

        This method sets up the MODELLER environment, runs the AutoModel,
        and renames the output file.

        Raises:
            Exception: If an error occurs during the MODELLER process.
        """

        class MyModel(AutoModel):
            def special_patches(self, aln):
                self.rename_segments(
                    segment_ids=["A", "B", "C"], renumber_residues=[1, 1, 1]
                )

        try:
            env = Environ(
                rand_seed=-8123,
                restyp_lib_file=str(self.restyp_lib),
                copy=None,
            )

            env.io.hetatm = True
            env.libs.topology.read(str(self.top_heav_lib))
            env.libs.parameters.read(str(self.par_mod_lib))
            template_dir = os.path.dirname(self.template_pdb)
            env.io.atom_files_directory = [".", template_dir]

            a = MyModel(
                env,
                alnfile=str(self.aligned_file),
                knowns="template",
                sequence="target",
            )

            a.very_fast()
            a.starting_model = 1
            a.ending_model = 1
            a.make()

            self.output_pdb = f"{self.output_prefix}_final_model.pdb"
            try:
                os.rename(a.outputs[0]["name"], self.output_pdb)
            except OSError as e:
                LOG.error(
                    f"Could not rename output file to {self.output_pdb}: {str(e)}"
                )
                raise
        except Exception as e:
            LOG.error(f"An error occurred during the MODELLER process: {str(e)}")
            raise

    @timeit
    def execute_modeller(self) -> None:
        """Time and delegate to run_modeller, which already logs and raises on failure."""
        self.run_modeller()


@timeit
def run_modeller(
    aligned_file: str,
    template_pdb: str,
    output_prefix: str,
    restyp_lib: str,
    top_heav_lib: str,
    par_mod_lib: str,
) -> str:
    """
    Run MODELLER to generate a protein structure model.

    Args:
        aligned_file (str): Path to the aligned sequence file.
        template_pdb (str): Path to the template PDB file.
        output_prefix (str): Prefix for output files.
        restyp_lib (str): Path to the residue type library file.
        top_heav_lib (str): Path to the topology library file.
        par_mod_lib (str): Path to the parameter library file.

    Returns:
        str: Path to the output PDB file.

    Raises:
        Exception: If an error occurs during the MODELLER process.
    """
    modeller = ModellerWrapper(
        str(aligned_file),
        str(template_pdb),
        str(output_prefix),
        str(restyp_lib),
        str(top_heav_lib),
        str(par_mod_lib),
    )
    modeller.execute_modeller()
    return modeller.output_pdb
