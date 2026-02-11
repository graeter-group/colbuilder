# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0

import subprocess
import os

from colbuilder.core.utils.logger import setup_logger
from colbuilder.core.utils.config import ColbuilderConfig

LOG = setup_logger(__name__)


class Chimera(object):
    """
    Generate collagen microfibril with the "crystal contacts" command from UCSF Chimera.

    Chimera needs to be installed beforehand from:
    https://www.cgl.ucsf.edu/chimera/download.html

    Important: Install Chimera (python 2.7 based) and not ChimeraX (python 3.)

    Attributes:
        pdb_file (str): Path to the PDB file.

    Methods:
        matrixget: Gets transformation matrices based on crystal contacts.
        matrixset: Sets PDB models based on updated, symmetrized transformation matrices.
        swapaa: Swaps amino acid defined by user for mutation.
    """

    def __init__(self, config, pdb=None):
        """
        Initialize Chimera object.

        Args:
            config (ColbuilderConfig): Configuration object containing settings.
            pdb (Optional[str]): Path to the PDB file. Defaults to None.
        """
        self.chimera_dir = config.CHIMERA_SCRIPTS_DIR
        self.pdb_file = pdb

    def matrixget(self, pdb=None, contact_distance=0, crystalcontacts=""):
        """
        Call Chimera via Python 2.7 script from terminal to get transformation matrices.

        Args:
            pdb (str, optional): Path to the PDB file. If None, uses self.pdb_file.
            contact_distance (int): Contact distance for crystal contacts.
            crystalcontacts (str): Output file name for crystal contacts.

        Returns:
            subprocess.CompletedProcess: Result of the subprocess run.
        """

        if pdb is None:
            pdb = self.pdb_file

        pdb_full_path = os.path.abspath(pdb)
        script_path = os.path.join(self.chimera_dir, "matrixget.py")

        if crystalcontacts.endswith(".txt"):
            crystalcontacts = crystalcontacts[:-4]

        env = os.environ.copy()
        env["PDB_FILE"] = pdb_full_path
        env["CONTACT_DISTANCE"] = str(contact_distance)
        env["CRYSTALCONTACTS_FILE"] = crystalcontacts

        chimera_command = ["chimera", "--nogui", "--silent", "--script", script_path]

        result = subprocess.run(
            " ".join(chimera_command),
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=env,
        )

        if result.returncode != 0:
            LOG.error(f"Chimera command failed with return code {result.returncode}")
            LOG.error(f"Chimera stdout: {result.stdout.decode()}")
            LOG.error(f"Chimera stderr: {result.stderr.decode()}")
            raise RuntimeError(f"Chimera command failed: {result.stderr.decode()}")

        return result

    def matrixset(self, pdb=None, crystalcontacts="", system_size=0, fibril_length=0.0):
        if pdb is None:
            pdb = self.pdb_file

        pdb_full_path = os.path.abspath(pdb)
        script_path = os.path.join(self.chimera_dir, "matrixset.py")

        crystalcontacts_str = str(crystalcontacts)
        if crystalcontacts_str.endswith(".txt"):
            crystalcontacts_str = crystalcontacts_str[:-4]

        env = os.environ.copy()
        env["PDB_FILE"] = pdb_full_path
        env["CRYSTALCONTACTS_FILE"] = crystalcontacts_str
        env["SYSTEM_SIZE"] = str(system_size)
        env["FIBRIL_LENGTH"] = str(fibril_length)

        chimera_command = ["chimera", "--nogui", "--silent", "--script", script_path]

        LOG.debug(f"     Running command: {' '.join(chimera_command)}")

        result = subprocess.run(
            " ".join(chimera_command),
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=env,
        )

        if result.returncode != 0:
            LOG.error(f"Chimera command failed with return code {result.returncode}")
            LOG.error(f"Chimera stdout: {result.stdout.decode()}")
            LOG.error(f"Chimera stderr: {result.stderr.decode()}")
            raise RuntimeError(f"Chimera command failed: {result.stderr.decode()}")

        expected_file = f"{crystalcontacts_str}_id.txt"
        if os.path.exists(expected_file):
            LOG.debug(f"        File created successfully: {expected_file}")
            with open(expected_file, "r") as f:
                pass
        else:
            LOG.error(f"    File not created: {expected_file}")
            raise FileNotFoundError(f"Expected file not created: {expected_file}")

        return result

    def swapaa(self, replace: str, system_type: str) -> subprocess.CompletedProcess:
        """
        Call Chimera via Python 2.7 script from terminal to swap amino acids.

        Args:
            replace (str): Path to the file containing amino acid replacement instructions (without .txt extension).
            system_type (str): Type of the system (directory name).

        Returns:
            subprocess.CompletedProcess: Result of the subprocess run.
        """
        swapaa_script = os.path.join(self.chimera_dir, "swapaa.py")

        if not os.path.exists(swapaa_script):
            LOG.error(f"Chimera swapaa script not found: {swapaa_script}")
            raise FileNotFoundError(f"Chimera swapaa script not found: {swapaa_script}")

        # Construct the command - our swapaa.py expects two arguments: replace_file and system_type
        # Use proper quoting for the script and arguments
        cmd = f'chimera --nogui --silent --script "{swapaa_script} {replace} {system_type}"'

        LOG.debug(f"Running Chimera command: {cmd}")

        try:
            result = subprocess.run(
                cmd,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )

            LOG.debug(
                f"Chimera command completed with return code: {result.returncode}"
            )
            if result.stdout:
                # LOG.debug(f"Chimera stdout: {result.stdout}")
                pass
            if result.stderr:
                # LOG.debug(f"Chimera stderr: {result.stderr}")
                pass
            return result
        except Exception as e:
            LOG.error(f"Error executing Chimera command: {str(e)}")
            return subprocess.CompletedProcess(
                args=cmd,
                returncode=1,
                stdout=b"",
                stderr=str(e).encode() if hasattr(str(e), "encode") else str(e),
            )
