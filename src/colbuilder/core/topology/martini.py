"""
Martini topology generation module.

This module implements the Martini 3 force field topology generation for molecular systems,
with a focus on collagen microfibrils. It provides functionality for:

- PDB file processing and manipulation
- Coarse-graining using Martinize2
- GO-like potential implementation
- System topology generation
- File organization and management

The module requires the Martinize2 tool and custom contact map utilities.
"""

import cmd
import os
import subprocess
import shutil
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple, Union
import asyncio
from tqdm import tqdm
from colorama import Fore, Style

from colbuilder.core.topology.itp import Itp
from colbuilder.core.topology.crosslink import Crosslink
from colbuilder.core.geometry.system import System
from colbuilder.core.utils.dec import timeit
from colbuilder.core.utils.config import ColbuilderConfig
from colbuilder.core.utils.exceptions import TopologyGenerationError
from colbuilder.core.utils.martinize_finder import (
    get_active_conda_env,
    find_and_install_custom_force_field,
    get_conda_command_with_path,
)
from colbuilder.core.utils.logger import setup_logger
from colbuilder.core.utils.files import FileManager

LOG = setup_logger(__name__)


class Martini:
    """
    Martini 3 force field topology generator.

    This class handles the generation and manipulation of molecular topologies using
    the Martini 3 force field, specifically designed for coarse-grained simulations
    of collagen systems.

    Attributes
    ----------
    system : Any
        Molecular system being processed
    ff : str
        Force field name
    is_line : tuple[str, ...]
        Valid PDB line identifiers
    is_chain : tuple[str, ...]
        Valid chain identifiers
    """

    def __init__(self, system: Any = None, ff: Optional[str] = None):
        """
        Initialize Martini topology generator.

        Parameters
        ----------
        system : Any
            Molecular system to process
        ff : Optional[str]
            Force field name
        """
        self.system = system
        self.ff = ff
        self.is_line = ("ATOM  ", "HETATM", "ANISOU", "TER   ")
        self.is_chain = ("A", "B", "C")

    def merge_pdbs(self, model_id: Optional[int] = None, cnt_model: Optional[int] = None) -> Optional[str]:
        """
        Merge multiple PDB files based on system connectivity.
        
        Merges all CG PDB files created for a model's connections into a single file.
        This includes both self-connections and cross-connections to ensure proper
        crosslink detection.
        """
        model = self.system.get_model(model_id=model_id)
        if model is None or model.connect is None:
            LOG.error(f"No model or connections for model_id {model_id}")
            return None

        if cnt_model is None:
            LOG.error("cnt_model is None in merge_pdbs")
            return None
            
        output_file = f"{int(cnt_model)}.merge.pdb"
        merged_count = 0

        try:
            with open(output_file, "w") as f:
                for connect_id in model.connect:
                    input_file = f"{int(model_id)}.{int(connect_id)}.CG.pdb"
                    
                    if not os.path.exists(input_file):
                        LOG.warning(f"CG PDB file not found: {input_file}")
                        continue

                    with open(input_file, "r") as infile:
                        lines_written = 0
                        for line in infile:
                            if line[0:6] in self.is_line:
                                f.write(line)
                                lines_written += 1
                        LOG.debug(f"Merged {lines_written} lines from {input_file}")
                        merged_count += 1

                f.write("END\n")
            
            if merged_count == 0:
                LOG.error(f"No CG files were merged for model {model_id}")
                if os.path.exists(output_file):
                    os.remove(output_file)
                return None
                
            LOG.debug(f"Successfully merged {merged_count} files into {output_file}")
            
            if len(model.connect) > 1:
                LOG.debug(f"Model {model_id} merge file contains {merged_count} structures from connections: {model.connect}")
            
            return output_file
        except Exception as e:
            LOG.error(f"Error merging PDBs: {str(e)}")
            return None

    def read_pdb(self, pdb_id: Optional[int] = None) -> List[str]:
        """
        Read PDB for Martinize2 processing.

        Parameters
        ----------
        pdb_id : Optional[int]
            The identifier for the PDB in the system

        Returns
        -------
        List[str]
            List of PDB file lines
        """
        pdb = []
        model = self.system.get_model(model_id=pdb_id)
        
        if model is None:
            LOG.error(f"Model not found for pdb_id: {pdb_id}")
            return pdb

        if model.type:
            if pdb_id is not None:
                file_path = Path(model.type) / f"{int(pdb_id)}.caps.pdb"
                if file_path.exists():
                    try:
                        with open(file_path, "r") as file:
                            pdb = [line for line in file if line[0:6] in self.is_line]
                        return pdb
                    except Exception as e:
                        LOG.error(f"Error reading PDB file {file_path}: {str(e)}")

        alt_paths = []
        if pdb_id is not None:
            alt_paths.append(Path(f"{int(pdb_id)}.caps.pdb")) 
            if model.type:
                alt_paths.append(Path(f"./{model.type}/{int(pdb_id)}.caps.pdb")) 

        for path in alt_paths:
            if path and path.exists():
                try:
                    with open(path, "r") as file:
                        pdb = [line for line in file if line[0:6] in self.is_line]
                    return pdb
                except Exception as e:
                    LOG.error(f"Error reading PDB file {path}: {str(e)}")

        LOG.error(f"PDB file not found for pdb_id: {pdb_id}")
        return pdb

    def set_pdb(self, pdb: Optional[List[str]] = None) -> Tuple[List[str], List[str]]:
        """
        Prepare PDBs for Martinize2 by renumbering residues and adding chain terminators.

        Parameters
        ----------
        pdb : Optional[List[str]]
            List of PDB file lines

        Returns
        -------
        Tuple[List[str], List[str]]
            Tuple containing (order, map) where order is the renumbered PDB and map is for contact mapping
        """
        if not pdb:
            LOG.warning("Empty PDB provided to set_pdb")
            return [], []

        try:
            first_cnt = int(pdb[1][22:26]) - 1
            cnt, cnt_map = 0, 0
            order, map = [], []
            chain_store = "A"

            for line in pdb:
                if line[0:3] == "TER":
                    continue

                if line[21:22] != chain_store:
                    order.append("TER\n")
                    if chain_store == "A":
                        chain_store = "B"
                    elif chain_store == "B":
                        chain_store = "C"

                if first_cnt < int(line[22:26]):
                    first_cnt = int(line[22:26])
                    cnt += 1
                    cnt_map += 1
                if first_cnt > int(line[22:26]) and str(line[21:22]) in self.is_chain:
                    first_cnt = int(line[22:26])
                    cnt = 1
                    cnt_map += 1

                if cnt < 10:
                    order.append(line[:22] + "   " + str(int(cnt)) + line[26:])
                elif 10 <= cnt < 100:
                    order.append(line[:22] + "  " + str(int(cnt)) + line[26:])
                elif 100 <= cnt < 1000:
                    order.append(line[:22] + " " + str(int(cnt)) + line[26:])
                elif 1000 <= cnt < 10000:
                    order.append(line[:22] + str(int(cnt)) + line[26:])

                if cnt_map < 10:
                    map.append(line[:22] + "   " + str(int(cnt_map)) + line[26:])
                elif 10 <= cnt_map < 100:
                    map.append(line[:22] + "  " + str(int(cnt_map)) + line[26:])
                elif 100 <= cnt_map < 1000:
                    map.append(line[:22] + " " + str(int(cnt_map)) + line[26:])
                elif 1000 <= cnt_map < 10000:
                    map.append(line[:22] + str(int(cnt_map)) + line[26:])

            return order, map
        except Exception as e:
            LOG.error(f"Error in set_pdb: {str(e)}")
            return [], []

    def cap_pdb(self, pdb: Optional[List[str]] = None) -> Tuple[List[str], str, str]:
        """
        Decide martinize2 -nter/-cter flags robustly.
        1) Rename terminal ALA -> CLA in the PDB *before* flag decisions.
        2) If any chain starts/ends on a special block, use 'none' for that side.
        Special blocks include ACE/CLA and crosslink blocks (LY2/LY3/L4Y/L5Y/LYX),
        and also NME (since it's already an explicit cap).
        """
        if not pdb:
            LOG.warning("Empty PDB provided to cap_pdb")
            return [], "none", "none"

        try:
            # 1) Rename terminal ALA -> CLA (bookkeeping only; never pass 'CLA' as a CLI mod)
            chain_length = self.get_chain_length(pdb)  # e.g. {"A": " 178", "B": "...", "C": "..."}
            for i in range(len(pdb)):
                line = pdb[i]
                if not line.startswith(("ATOM  ", "HETATM")):
                    continue
                if line[17:20] == "ALA":
                    tag = line[21:26]  # e.g., "A 178"
                    if tag in {f"A{chain_length['A']}", f"B{chain_length['B']}", f"C{chain_length['C']}"}:
                        pdb[i] = line[:17] + "CLA " + line[21:]

            # 2) Inspect first/last residue per chain *after* the rename above
            first_res = {}   # chain -> resname
            last_res  = {}   # chain -> (resid, resname)
            for line in pdb:
                if not line.startswith(("ATOM  ", "HETATM")):
                    continue
                ch = line[21:22]
                if ch not in self.is_chain:
                    continue
                resname = line[17:20].strip()
                try:
                    resid = int(line[22:26])
                except ValueError:
                    continue
                if ch not in first_res:
                    first_res[ch] = resname
                if ch not in last_res or resid >= last_res[ch][0]:
                    last_res[ch] = (resid, resname)

            firsts = set(first_res.values()) if first_res else set()
            lasts  = {name for _, name in last_res.values()} if last_res else set()

            # Anything here means "don't ask martinize2 to apply a terminal mod":
            special_first = {"ACE", "CLA", "LY2", "LY3", "L4Y", "L5Y", "LYX"}
            special_last  = {"ACE", "CLA", "LY2", "LY3", "L4Y", "L5Y", "LYX", "NME"}

            # N-terminus: be conservative (ACE often causes issues on GLN etc. in some FF)
            nter_flag = "none" if (not firsts or (firsts & special_first)) else "none"

            # C-terminus: allow NME only if all chains end on standard residues (not special, not NME)
            cter_flag = "NME"
            if not lasts or (lasts & special_last):
                cter_flag = "none"

            # LOG.debug(f"cap_pdb decided: -nter {nter_flag}, -cter {cter_flag} (first={firsts}, last={lasts})")
            return pdb, cter_flag, nter_flag

        except Exception as e:
            LOG.error(f"Error in cap_pdb: {str(e)}")
            return pdb or [], "none", "none"

    def get_chain_length(self, pdb: Optional[List[str]] = None) -> Dict[str, str]:
        """
        Get length for each chain of the triple helix.

        Parameters
        ----------
        pdb : Optional[List[str]]
            List of PDB file lines

        Returns
        -------
        Dict[str, str]
            Dictionary with chain IDs as keys and lengths as values
        """
        if not pdb:
            LOG.warning("Empty PDB provided to get_chain_length")
            return {"A": "", "B": "", "C": ""}

        try:
            chain_length = {key: "" for key in ["A", "B", "C"]}

            for line_it in range(len(pdb) - 1):
                if pdb[line_it][21:22] == "A" and pdb[line_it + 1][21:22] == "B":
                    chain_length["A"] = pdb[line_it][22:26]
                if pdb[line_it][21:22] == "B" and pdb[line_it + 1][21:22] == "C":
                    chain_length["B"] = pdb[line_it][22:26]

            if len(pdb) > 1 and pdb[-1][21:22] == "C":
                chain_length["C"] = pdb[-1][22:26]

            return chain_length
        except Exception as e:
            LOG.error(f"Error in get_chain_length: {str(e)}")
            return {"A": "", "B": "", "C": ""}

    def write_pdb(
        self, pdb: Optional[List[str]] = None, file: Optional[str] = None
    ) -> None:
        """
        Write PDB lines to a file.

        Parameters
        ----------
        pdb : Optional[List[str]]
            List of PDB file lines
        file : Optional[str]
            Path to the output file
        """
        if not pdb or not file:
            LOG.warning(f"Missing parameters in write_pdb")
            return

        try:
            with open(file, "w") as f:
                for line in pdb:
                    if line[0:3] != "END":
                        f.write(line)
                f.write("END")
        except PermissionError:
            LOG.error(f"Permission denied when writing to file: {file}")
            raise
        except Exception as e:
            LOG.error(f"Error writing PDB to file: {str(e)}")

    def get_system_pdb(self, size: Optional[int] = None) -> List[str]:
        """
        Combine all model PDBs into a system PDB.

        Parameters
        ----------
        size : Optional[int]
            Number of models to include

        Returns
        -------
        List[str]
            Combined PDB lines for the system
        """
        if not size:
            LOG.warning("No size provided to get_system_pdb")
            return []

        pdb = []
        found_models = 0

        try:
            for cnt_model in range(size):
                merge_pdb_path = f"{int(cnt_model)}.merge.pdb"
                if not os.path.exists(merge_pdb_path):
                    LOG.warning(f"Merged PDB file not found: {merge_pdb_path}")
                    continue

                with open(merge_pdb_path, "r") as f:
                    model_lines = f.readlines()
                    pdb.extend(model_lines)
                    found_models += 1

            LOG.debug(f"Found {found_models} out of {size} expected merge files")
            return pdb
        except Exception as e:
            LOG.error(f"Error getting system PDB: {str(e)}")
            return []

    def write_system_topology(
        self, topology_file: str = "system.top", size: Optional[int] = None
    ) -> None:
        """
        Write final topology file for the system.

        Parameters
        ----------
        topology_file : str
            Path to the output topology file
        size : Optional[int]
            Number of models to include
        """
        if not size:
            LOG.warning("No size provided to write_system_topology")
            return

        try:
            with open(topology_file, "w") as f:
                f.write("; This is the topology for the collagen microfibril\n")
                f.write("#define GO_VIRT\n")
                f.write('#include "martini_v3.0.0.itp"\n')
                f.write('#include "go-sites.itp"\n\n')

                for m in range(size):
                    f.write(f'#include "col_{m}.itp"\n')
                    f.write(f'#include "col_{m}_go-excl.itp"\n')

                f.write('\n#include "martini_v3.0.0_solvents_v1.itp"\n')
                f.write('#include "martini_v3.0.0_ions_v1.itp"\n')

                f.write("\n[ system ]\n")
                f.write("Collagen, Martini 3 and Go-Potentials \n")
                f.write("\n[ molecules ]\n")

                for t in range(size):
                    f.write(f"col_{t}     1\n")

            self.write_go_topology(name_type="sites.itp")

        except PermissionError:
            LOG.error(f"Permission denied when writing to system topology file: {topology_file}")
            raise
        except Exception as e:
            LOG.error(f"Error writing system topology: {str(e)}")

    def write_go_topology(self, name_type: Optional[str] = None) -> None:
        """
        Write topology for go-like potentials.

        Parameters
        ----------
        name_type : Optional[str]
            Type of GO topology file to write
        """
        if not name_type:
            LOG.warning("No name_type provided to write_go_topology")
            return

        try:
            with open(f"go-{name_type}", "w") as f:
                f.write(f'#include "col_go-{name_type}"\n')

            with open(f"col_go-{name_type}", "w") as f:
                f.write("[ atomtypes ]\n")
                f.write("; protein BB virtual particle\n")
                f.write("col 0.0 0.000 A 0.0 0.0 \n")

        except PermissionError:
            LOG.error("Permission denied when writing GO topology files")
            raise
        except Exception as e:
            LOG.error(f"Error writing GO topology: {str(e)}")

    def translate_pdb(self, pdb: Optional[List[str]] = None) -> List[str]:
        if not pdb:
            LOG.warning("Empty PDB provided to translate_pdb")
            return []

        translated = []
        for line in pdb:
            if line[0:6] in self.is_line:
                try:
                    # Columns: x[30:38], y[38:46], z[46:54]
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    new_line = f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}"
                    translated.append(new_line)
                except ValueError:
                    # If parse fails, keep original to avoid data loss but log it
                    LOG.warning(f"Bad coord formatting, keeping original: {line.rstrip()}")
                    translated.append(line)
            else:
                translated.append(line)
        return translated


    def write_gro(
        self,
        system: Optional[Any] = None,
        gro_file: Optional[str] = None,
        processed_models: Optional[List[int]] = None,
    ) -> None:
        """
        Write a GRO (Gromos87) file for the processed models.

        This is a placeholder method for API compatibility with the Amber class.
        Martini uses PDB files primarily, but this method could be implemented
        to convert PDB to GRO format if needed.
        """
        LOG.debug("Write_gro called but not implemented for Martini - using PDB format instead")
        return None


def check_output_files(topology_dir: Path, expected_models: int):
    """Diagnostic function to check which files were created"""
    cg_files = list(topology_dir.glob("*.CG.pdb"))
    merge_files = list(topology_dir.glob("*.merge.pdb"))
    itp_files = list(topology_dir.glob("col_*.itp"))
    
    LOG.debug(f"CG files found: {len(cg_files)}")
    LOG.debug(f"Merge files found: {len(merge_files)} (expected: {expected_models})")
    LOG.debug(f"ITP files found: {len(itp_files)}")
    
    for i in range(expected_models):
        if not (topology_dir / f"{i}.merge.pdb").exists():
            LOG.warning(f"Missing merge file: {i}.merge.pdb")
    
    return len(merge_files)


@timeit
async def build_martini3(
    system: System, config: ColbuilderConfig, file_manager: Optional[FileManager] = None
) -> Martini:
    """
    Build a Martini 3 topology for the given molecular system.
    """
    ff = f"{config.force_field}"
    go_epsilon = getattr(config, "go_epsilon", 9.414)

    if file_manager is None:
        file_manager = FileManager(config)

    topology_dir = file_manager.get_temp_dir("topology_gen")
    original_dir = Path.cwd()

    martini = Martini(system=system, ff=ff)
    steps = 4
    cnt_model = 0

    try:
        os.chdir(topology_dir)

        source_ff_dir = config.FORCE_FIELD_DIR
        source_contactmap_dir = source_ff_dir / "contactmap"

        local_contactmap_dir = Path("contactmap")
        local_contactmap_dir.mkdir(exist_ok=True)

        LOG.info(f"Step 1/{steps} Using coordinates from geometry generation")
        try:
            LOG.debug("Skipping system translation - using original coordinates from geometry generation")
        except Exception as e:
            LOG.error(f"Error in Step 1: {str(e)}")
            raise TopologyGenerationError(
                message="Failed in coordinate setup for Martini topology",
                original_error=e,
                error_code="TOP_MART_001",
                context={"crystal": str(system.crystal)},
            )

        LOG.info(f"Step 2/{steps} Setting up Martini force field files")
        try:
            if not source_ff_dir.exists():
                raise TopologyGenerationError(
                    message=f"Martini force field directory not found: {source_ff_dir}",
                    error_code="TOP_MART_005",
                    context={"force_field_dir": str(source_ff_dir)},
                )

            go_script = source_ff_dir / "create_goVirt.py"
            if go_script.exists():
                shutil.copy2(go_script, Path("create_goVirt.py"))
            else:
                LOG.warning(f"create_goVirt.py not found in force field directory: {go_script}")

            contact_map_path = local_contactmap_dir / "contact_map"
            if not contact_map_path.exists():
                source_executable = source_contactmap_dir / "contact_map"
                if source_executable.exists():
                    shutil.copy2(source_executable, contact_map_path)
                    contact_map_path.chmod(contact_map_path.stat().st_mode | 0o111)
                else:
                    for src_file in source_contactmap_dir.glob("*.c"):
                        shutil.copy2(src_file, local_contactmap_dir / src_file.name)
                    for header_file in source_contactmap_dir.glob("*.h"):
                        shutil.copy2(header_file, local_contactmap_dir / header_file.name)

                    source_makefile = source_contactmap_dir / "makefile"
                    if source_makefile.exists():
                        shutil.copy2(source_makefile, local_contactmap_dir / "makefile")
                        make_process = await asyncio.create_subprocess_shell(
                            "make",
                            cwd=str(local_contactmap_dir),
                            stdout=asyncio.subprocess.PIPE,
                            stderr=asyncio.subprocess.PIPE,
                        )
                        stdout, stderr = await make_process.communicate()
                        if make_process.returncode != 0:
                            LOG.error(f"Failed to build contact map tool: {stderr.decode()}")
                        else:
                            LOG.debug("Successfully built contact map tool")
                    else:
                        source_makefile = source_ff_dir / "makefile"
                        if source_makefile.exists():
                            shutil.copy2(source_makefile, Path("makefile"))
                            make_process = await asyncio.create_subprocess_shell(
                                "make",
                                stdout=asyncio.subprocess.PIPE,
                                stderr=asyncio.subprocess.PIPE,
                            )
                            stdout, stderr = await make_process.communicate()
                            if make_process.returncode != 0:
                                LOG.error(f"Failed to build contact map tool: {stderr.decode()}")
                            else:
                                LOG.debug("Successfully built contact map tool")
                        else:
                            LOG.warning("No Makefile found to build contact map tool")

            if not contact_map_path.exists():
                LOG.warning("Contact map tool could not be found or built, may cause issues")

        except Exception as e:
            LOG.error(f"Error setting up Martini force field files: {str(e)}")
            raise TopologyGenerationError(
                message="Failed to set up Martini force field files",
                original_error=e,
                error_code="TOP_MART_005",
                context={"force_field_dir": str(config.FORCE_FIELD_DIR)},
            )

        LOG.info(f"Step 3/{steps} Processing models with Martinize2")
        processed_models = []

        try:
            martinize2_command = getattr(config, "martinize2_command", None)
            if not martinize2_command:
                from shutil import which
                martinize2_cmd = which("martinize2")
                if martinize2_cmd:
                    config.martinize2_command = "martinize2"
                else:
                    raise TopologyGenerationError(
                        message="Martinize2 command not found in configuration or PATH",
                        error_code="TOP_MART_005",
                        context={"force_field": ff},
                    )
        except Exception as e:
            LOG.error(f"Error checking for Martinize2 command: {str(e)}")
            raise TopologyGenerationError(
                message="Failed to configure Martinize2 command",
                original_error=e,
                error_code="TOP_MART_005",
                context={"force_field": ff},
            )

        LOG.info(f"{Fore.BLUE}Building coarse-grained topology:{Style.RESET_ALL}")

        if len(list(system.get_models())) > 0:
            first_model = system.get_model(model_id=list(system.get_models())[0])
            model_type = first_model.type
            type_dir = topology_dir / model_type
            type_dir.mkdir(exist_ok=True, parents=True)

        models_list = [model_id for model_id in system.get_models()]
        model_status = {}
        failed_martinize, failed_contact, failed_go, failed_itp = [], [], [], []
        processed_in_topology: set = set()

        for model_id in tqdm(models_list, desc="Building topology", unit="%"):
            model = system.get_model(model_id=model_id)
            if model is None or model.connect is None:
                LOG.warning(f"Skipping model {model_id}: No connections found")
                model_status[model_id] = "no_connections"
                continue

            if model_id in processed_in_topology:
                LOG.debug(f"Skipping model {model_id}: Already included in another model's topology file")
                model_status[model_id] = "already_processed"
                continue

            try:
                for connect_id in model.connect:
                    try:
                        pdb = martini.read_pdb(pdb_id=connect_id)
                        if not pdb:
                            LOG.warning(f"Empty PDB for connect_id: {connect_id}")
                            continue

                        pdb = martini.translate_pdb(pdb=pdb)
                        pdb_capped, cter, nter = martini.cap_pdb(pdb=pdb)
                        order, map_pdb = martini.set_pdb(pdb=pdb_capped)

                        martini.write_pdb(pdb=order, file="tmp.pdb")
                        martini.write_pdb(pdb=map_pdb, file="map.pdb")

                        def _martinize_cmd(nter_flag: str, cter_flag: str) -> str:
                            args = (
                                f"-f tmp.pdb -sep -merge A,B,C "
                                f"-collagen -from amber99 -o topol.top -bonds-fudge 1.4 -p backbone "
                                f"-ff {martini.ff}00C -x {int(model_id)}.{int(connect_id)}.CG.pdb "
                                f"-nter {nter_flag} -cter {cter_flag} "
                                f"-govs-include -govs-moltype col_{int(model_id)}.{int(connect_id)} -maxwarn 4"
                            )
                            return get_conda_command_with_path("martinize2", args)

                        cmd = _martinize_cmd(nter, cter)
                        process = await asyncio.create_subprocess_shell(
                            cmd, stdout=asyncio.subprocess.PIPE, stderr=asyncio.subprocess.PIPE
                        )
                        stdout, stderr = await process.communicate()
                        stderr_txt = stderr.decode(errors="ignore")
                        stdout_txt = stdout.decode(errors="ignore")

                        if process.returncode != 0:
                            LOG.error(f"Martinize2 failed for model {model_id} -> connection {connect_id}")
                            LOG.error(f"martinize2 stderr:\n{stderr_txt}")
                            LOG.debug(f"martinize2 stdout:\n{stdout_txt}")
                            failed_martinize.append((model_id, connect_id))
                            continue

                        cg_file = f"{int(model_id)}.{int(connect_id)}.CG.pdb"
                        if not os.path.exists(cg_file):
                            LOG.error(f"CG file not created: {cg_file}")
                        else:
                            LOG.debug(f"CG file created successfully: {cg_file}")

                        contact_cmd = f"./contact_map ../map.pdb > ../map.out"
                        contact_process = await asyncio.create_subprocess_shell(
                            contact_cmd,
                            cwd=str(local_contactmap_dir),
                            stdout=asyncio.subprocess.PIPE,
                            stderr=asyncio.subprocess.PIPE,
                        )
                        contact_stdout, contact_stderr = await contact_process.communicate()
                        if contact_process.returncode != 0:
                            LOG.error(f"Contact map failed for model {model_id} -> connection {connect_id}")
                            failed_contact.append((model_id, connect_id))
                            continue

                        go_cmd = (
                            f"python create_goVirt.py -s {int(model_id)}.{int(connect_id)}.CG.pdb "
                            f"-f map.out --moltype col_{int(model_id)}.{int(connect_id)} --go_eps {go_epsilon}"
                        )
                        go_process = await asyncio.create_subprocess_shell(
                            go_cmd,
                            stdout=asyncio.subprocess.PIPE,
                            stderr=asyncio.subprocess.PIPE,
                        )
                        go_stdout, go_stderr = await go_process.communicate()
                        if go_process.returncode != 0:
                            LOG.error(f"GO virtual sites creation failed for model {model_id} -> connection {connect_id}")
                            failed_go.append((model_id, connect_id))
                            continue

                    except Exception as e:
                        LOG.error(f"Error processing model {model_id} -> connection {connect_id}: {str(e)}")
                        continue

                merged_pdb = martini.merge_pdbs(model_id=int(model_id), cnt_model=cnt_model)
                if merged_pdb:
                    model_status[model_id] = "merged"
                    try:
                        itp_ = Itp(system=system, model_id=int(model_id))
                        itp_.read_model(model_id=int(model_id))
                        itp_.go_to_pairs(model_id=int(model_id))
                        itp_.make_topology(model_id=int(model_id), cnt_model=cnt_model)
                        processed_models.append(model_id)
                        for connect_id in model.connect:
                            processed_in_topology.add(connect_id)
                    except Exception as e:
                        LOG.error(f"Error processing ITP for model {model_id}: {str(e)}")
                        failed_itp.append(model_id)
                else:
                    model_status[model_id] = "no_merge_file"
                    LOG.error(f"Model {model_id} failed: No merged PDB created")

                cnt_model += 1

            except Exception as e:
                LOG.error(f"Error processing model {model_id}: {str(e)}")
                model_status[model_id] = f"error: {str(e)}"

        LOG.info(f"Successfully processed: {len(processed_models)} models")

        if failed_martinize:
            LOG.warning(f"Failed at Martinize2: {failed_martinize}")
        if failed_contact:
            LOG.warning(f"Failed at contact map: {failed_contact}")
        if failed_go:
            LOG.warning(f"Failed at GO creation: {failed_go}")
        if failed_itp:
            LOG.warning(f"Failed at ITP processing: {failed_itp}")

        check_output_files(Path.cwd(), cnt_model)

        if not processed_models:
            raise TopologyGenerationError(
                message="No models were successfully processed with Martinize2",
                error_code="TOP_MART_002",
            )

        LOG.info(f"Step 4/{steps} Creating system PDB and topology")

        output_topology_dir = file_manager.ensure_dir(f"{config.species}_{ff}_topology_files")

        try:
            system_pdb = martini.get_system_pdb(size=cnt_model)
            pdb_file_path = Path(f"collagen_fibril_CG_{config.species}.pdb")
            martini.write_pdb(pdb=system_pdb, file=str(pdb_file_path))
            if pdb_file_path.exists():
                file_manager.copy_to_directory(pdb_file_path, dest_dir=output_topology_dir)
        except Exception as e:
            raise TopologyGenerationError(
                message="Failed to create system PDB file",
                original_error=e,
                error_code="TOP_MART_003",
                context={"output": config.species},
            )

        try:
            final_topology_file = f"collagen_fibril_{config.species}.top"
            martini.write_system_topology(topology_file=final_topology_file, size=cnt_model)

            topology_file_path = Path(final_topology_file)
            if topology_file_path.exists():
                file_manager.copy_to_directory(topology_file_path, dest_dir=output_topology_dir)

            for itp_file in Path().glob("col_[0-9]*.itp"):
                file_manager.copy_to_directory(itp_file, dest_dir=output_topology_dir)

            for excl_file in Path().glob("col_[0-9]*_go-excl.itp"):
                file_manager.copy_to_directory(excl_file, dest_dir=output_topology_dir)

            for go_site_file in Path().glob("*go-sites.itp"):
                file_manager.copy_to_directory(go_site_file, dest_dir=output_topology_dir)

            source_martini_files = list(source_ff_dir.glob("martini_v3.0.0*"))
            for source_file in source_martini_files:
                file_manager.copy_to_directory(source_file, dest_dir=output_topology_dir)

            if not source_martini_files:
                LOG.warning("No Martini force field files found in source directory")

        except Exception as e:
            raise TopologyGenerationError(
                message="Failed to create system topology files",
                original_error=e,
                error_code="TOP_MART_004",
                context={"output": config.species},
            )

        try:
            if not config.debug:
                subprocess.run(r"rm \#*", shell=True, check=False)
        except Exception as e:
            LOG.warning(f"Error cleaning up temporary files: {str(e)}")

        LOG.info(f"{Fore.BLUE}Martini topology generated successfully for {len(processed_models)} models.{Style.RESET_ALL}")

        os.chdir(original_dir)
        return martini

    except TopologyGenerationError:
        os.chdir(original_dir)
        raise
    except Exception as e:
        os.chdir(original_dir)
        raise TopologyGenerationError(
            message="Unexpected error in Martini topology generation",
            original_error=e,
            error_code="TOP_MART_001",
            context={
                "force_field": ff,
                "error_details": str(e),
                "location": "build_martini3",
            },
        )
