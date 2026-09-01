# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0

from pymol import cmd, editor
import subprocess
from typing import List, Dict, Any, Optional
import os
import shutil
from pathlib import Path

from colbuilder.core.utils.logger import setup_logger

LOG = setup_logger(__name__)


class Caps:
    """
    Adding Caps to a single triple helix.
    Original version implemented by Agnieszka Obarska-Kosinska taken from

    Obarska-Kosinska A, Rennekamp B, Ünal A, Gräter F.
    ColBuilder: A server to build collagen fibril models.
    Biophys J. 2021 Sep 7;120(17):3544-3549.
    doi: 10.1016/j.bpj.2021.07.009.
    Epub 2021 Jul 13. PMID: 34265261; PMCID: PMC8456305.

    Attributes:
        system (Any): The system object containing PDB information.
        system_size (int): Size of the system.
        chains (List[str]): List of chain identifiers.
        caps (List[str]): List of cap types.
        is_line (tuple): Tuple of valid line types in PDB file.
        chain_length (Dict[str, int]): Dictionary to store chain lengths.
        model (Dict[str, List[int]]): Dictionary to store residue numbers for each chain.
    """

    def __init__(self, system: Any):
        self.system = system
        self.system_size = system.size
        self.chains: List[str] = ["A", "B", "C"]
        self.caps: List[str] = ["N", "C"]
        self.is_line: tuple = ("ATOM  ", "HETATM", "ANISOU", "TER   ")
        self.chain_length: Dict[str, int] = {k: 0 for k in self.chains}
        self.model: Dict[str, List[int]] = {k: [] for k in self.chains}
        self.get_chain_length(system=system)

    def read_residues(self, pdb_id: int) -> None:
        """
        Reads pdb-file for each chain in the triple helix to obtain
        exact location of where to place the cap.

        Args:
            pdb_id (int): PDB file identifier for single triple helix.

        Raises:
            FileNotFoundError: If the PDB file is not found.
        """
        self.model = {k: [] for k in self.chains}
        pdb_file = f"{pdb_id}.pdb"

        if not os.path.exists(pdb_file):
            raise FileNotFoundError(f"PDB file not found: {pdb_file}")

        with open(pdb_file, "r") as f:
            for line in f:
                if line[0:6] in self.is_line and line[13:15] == "CA":
                    chain = line[21]
                    if chain in self.chains:
                        self.model[chain].append(int(line[22:26]))

    def get_line(self, cap: str, chain_id: str) -> str:
        """
        Writes command to be used in pymol to add a cap.

        Args:
            cap (str): Cap type ('N' or 'C').
            chain_id (str): Chain identifier.

        Returns:
            str: PyMOL command to add cap.

        Raises:
            ValueError: If an invalid cap type is provided.
        """
        if cap not in self.caps:
            raise ValueError(f"Invalid cap type: {cap}. Must be 'N' or 'C'.")

        index = 0 if cap == "N" else -1
        return f"resi {self.model[chain_id][index]} and chain {chain_id} and name {cap}"

    def get_chain_length(self, system: Any) -> None:
        """
        Gets the chain length for the initial model.

        Args:
            system (Any): System object containing PDB file information.
        """
        self.read_residues(pdb_id=system.crystal.pdb_file)
        for chain in self.model:
            self.chain_length[chain] = len(self.model[chain])

    def add_caps(
        self,
        pdb_id: int,
        crosslink_type: Optional[str] = None,
        temp_dir: Optional[Path] = None,
    ) -> str:
        """
        Adds caps to both ends of each model.

        Args:
            pdb_id (int): PDB identifier
            crosslink_type (Optional[str]): Type of crosslink. Uses 'NC' for non-crosslinked structures
            temp_dir (Optional[Path]): Directory for output files

        Returns:
            str: Path to the new PDB file with caps
        """
        cwd = Path.cwd()
        pdb_file = f"{pdb_id}.pdb"

        if not os.path.exists(pdb_file):
            LOG.info(f"PDB file not found: {pdb_file}")
            raise FileNotFoundError(f"PDB file not found: {pdb_file}")

        if temp_dir is None:
            output_dir = cwd.resolve()
        else:
            temp_dir = Path(temp_dir).resolve()

            # Check if we're already in a mix context (D or T directory)
            if temp_dir.name in ["D", "T"] and "mix_crosslinks" in str(temp_dir):
                output_dir = temp_dir
            elif "geometry_gen" in str(temp_dir) and not "mix_crosslinks" in str(
                temp_dir
            ):
                model_type = crosslink_type or "NC"
                type_dir = temp_dir / model_type
                if not type_dir.exists():
                    type_dir.mkdir(parents=True, exist_ok=True)
                output_dir = type_dir
            else:
                output_dir = temp_dir

        if not output_dir.exists():
            output_dir.mkdir(parents=True, exist_ok=True)

        cmd.load(pdb_file)
        for cap in self.caps:
            for chain in self.chains:
                line_cap = self.get_line(cap=cap, chain_id=chain)
                resi = int(line_cap.split(" ")[1])
                if (cap == "N" and resi != 1) or (
                    cap == "C" and resi != int(self.chain_length[chain])
                ):
                    cmd.edit(line_cap)
                    editor.attach_amino_acid(
                        "pk1", "ace" if cap == "N" else "nme", ss=0
                    )

        cmd.save("tmp.pdb")
        cmd.delete(name=str(pdb_id))
        caps_output_file = output_dir / f"{pdb_id}.caps.pdb"

        result = self.write_caps(pdb="tmp.pdb", pdb_id=pdb_id, output_dir=output_dir)
        return result

    def write_caps(self, pdb: str, pdb_id: int, output_dir: Path) -> str:
        """
        Write PDB file with caps.

        Args:
            pdb (str): Path to temporary PDB file.
            pdb_id (int): PDB identifier.
            output_dir (Path): Directory to save the output file.

        Returns:
            str: Path to the new PDB file with caps.

        Raises:
            FileNotFoundError: If the temporary PDB file is not found.
        """
        if not os.path.exists(pdb):
            raise FileNotFoundError(f"Temporary PDB file not found: {pdb}")

        output_file = output_dir / f"{pdb_id}.caps.pdb"

        with open(pdb, "r") as f_in, open(output_file, "w") as f_out:
            for line in f_in:
                if line[0:6] in self.is_line and line[0:3] != "TER":
                    if line.startswith("HETATM"):
                        line = "ATOM  " + line[6:]
                    f_out.write(line)
                    if (line[17:20] == "NME" and line[12:16] == "3HH3") or (
                        line[13:16] == "OXT"
                    ):
                        f_out.write("TER \n")

        # Cleanup temporary files
        try:
            if os.path.exists(f"{pdb_id}.pdb"):
                os.remove(f"{pdb_id}.pdb")
            if os.path.exists(pdb):
                os.remove(pdb)
        except Exception as e:
            LOG.warning(f"Error cleaning up temporary files: {e}")

        return str(output_file)
