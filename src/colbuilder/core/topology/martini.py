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

import os
import re
import sys
import subprocess
import shutil
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple, Union
import asyncio
from tqdm import tqdm
from colorama import Fore, Style

from colbuilder.core.topology.itp import Itp
from colbuilder.core.topology.backbone_repair import repair_backbone_bonds
from colbuilder.core.geometry.system import System
from colbuilder.core.utils.dec import timeit
from colbuilder.core.utils.config import ColbuilderConfig
from colbuilder.core.utils.exceptions import TopologyGenerationError
from colbuilder.core.utils.martinize_finder import (
    get_active_conda_env,
    find_and_install_custom_force_field,
    get_conda_command_with_path,
    find_martinize2_executable,
)
from colbuilder.core.utils.logger import setup_logger
from colbuilder.core.utils.files import FileManager

LOG = setup_logger(__name__)

# martinize2 -bonds-fudge used for every model (see _martinize_cmd in build_martini3).
MARTINIZE2_BONDS_FUDGE = 1.4
# vermouth (MakeBonds) guesses an inter-residue bond only below 0.5 * (r_vdW(C) + r_vdW(N)) * fudge,
# with the Bondi radii C 1.70 A and N 1.55 A, i.e. 2.275 A at fudge 1.4. Used for logging only.
VERMOUTH_CN_GUESS_LIMIT = 0.5 * (1.70 + 1.55) * MARTINIZE2_BONDS_FUDGE
# Consecutive residues of a chain whose peptide C-N is longer than this are not connected
# explicitly (a real gap, not a strained bond); the chain checks below then stop the build.
MAX_PEPTIDE_CN_ANGSTROM = 4.0
# Crosslink residues whose backbone must stay connected to its chain neighbours.
CROSSLINK_RESNAMES = frozenset({"LY2", "LY3", "L4Y", "L5Y", "LYX", "AGS", "LGX", "LZS", "LZD"})


def _pdb_residue_runs(lines: List[str]) -> List[Tuple[str, List[Tuple[str, str, str]]]]:
    """Return the chains of PDB lines in file order as (chain, [(resname, resid, icode), ...]).

    A new chain starts at every TER record and wherever the chain identifier changes.
    """
    runs: List[Tuple[str, List[Tuple[str, str, str]]]] = []
    cur_chain: Optional[str] = None
    cur: List[Tuple[str, str, str]] = []
    last_key = None
    for line in lines:
        if line.startswith("TER"):
            if cur:
                runs.append((cur_chain, cur))
            cur, cur_chain, last_key = [], None, None
            continue
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        chain = line[21:22]
        key = (line[17:20].strip(), line[22:26].strip(), line[26:27].strip())
        if chain != cur_chain:
            if cur:
                runs.append((cur_chain, cur))
            cur, cur_chain, last_key = [], chain, None
        if key != last_key:
            cur.append(key)
            last_key = key
    if cur:
        runs.append((cur_chain, cur))
    return runs


def _cg_backbone_beads(cg_pdb: Union[str, Path]) -> List[Optional[Tuple[int, str, str, int]]]:
    """BB beads of a CG PDB in atom order as (atom_index, resname, chain, resid); None marks a TER record.

    atom_index is the 1-based running index of the ATOM/HETATM record, which equals the
    atom number of the corresponding (merged) itp.
    """
    beads: List[Optional[Tuple[int, str, str, int]]] = []
    idx = 0
    with open(cg_pdb) as fh:
        for line in fh:
            if line.startswith(("ATOM  ", "HETATM")):
                idx += 1
                if line[12:16].strip() == "BB":
                    beads.append((idx, line[17:20].strip(), line[21:22], int(line[22:26])))
            elif line.startswith("TER") and beads and beads[-1] is not None:
                beads.append(None)
    return beads


def _bead_chain_runs(beads: List[Optional[Tuple[int, str, str, int]]]) -> List[List[Tuple[int, str, str, int]]]:
    """Split BB beads into chains: a chain ends at a TER record or where the chain identifier changes."""
    runs: List[List[Tuple[int, str, str, int]]] = []
    cur: List[Tuple[int, str, str, int]] = []
    for bead in beads:
        if bead is None or (cur and bead[2] != cur[-1][2]):
            if cur:
                runs.append(cur)
            cur = []
        if bead is not None:
            cur.append(bead)
    if cur:
        runs.append(cur)
    return runs


def check_cg_residue_order(aa_pdb: List[str], cg_pdb: Union[str, Path], label: str) -> None:
    """Raise if martinize2 did not keep the residues of the input PDB, chain by chain and in order.

    create_goVirt.py maps the contacts of the all-atom map (residues numbered in input order) onto
    the CG backbone beads by position, and backbone_repair only joins consecutive residues, so any
    reordering or splitting of a chain by martinize2 silently corrupts the Go model and the backbone.
    """
    aa_runs = [(chain, [res[0] for res in residues]) for chain, residues in _pdb_residue_runs(aa_pdb)]
    cg_runs = [(run[0][2], [bead[1] for bead in run]) for run in _bead_chain_runs(_cg_backbone_beads(cg_pdb))]
    if aa_runs == cg_runs:
        return
    problems = []
    if [c for c, _ in aa_runs] != [c for c, _ in cg_runs]:
        problems.append(
            f"chains in input {[f'{c}:{len(r)}' for c, r in aa_runs]} but in martinize2 output "
            f"{[f'{c}:{len(r)}' for c, r in cg_runs]} (chain:residues)"
        )
    for (ca, ra), (cc, rc) in zip(aa_runs, cg_runs):
        if (ca, ra) != (cc, rc):
            first = next((i for i, (x, y) in enumerate(zip(ra, rc)) if x != y), min(len(ra), len(rc)))
            problems.append(
                f"chain {ca}: {len(ra)} residues in input, chain {cc}: {len(rc)} BB beads in output, "
                f"first difference at position {first + 1} (input {ra[first:first + 3]}, output {rc[first:first + 3]})"
            )
            break
    raise TopologyGenerationError(
        message=(
            f"{label}: martinize2 changed the residue order or chain composition of the input "
            f"({'; '.join(problems)}). The Go contacts and backbone of this model would be wrong."
        ),
        error_code="TOP_MART_008",
        context={"model": label, "cg_pdb": str(cg_pdb)},
    )


def check_go_contacts_per_chain(go_table: Union[str, Path], cg_pdb: Union[str, Path], moltype: str) -> Dict[str, int]:
    """Count the Go contacts touching each chain of one martinize2 model; raise if a chain has none.

    The Go virtual site ``{moltype}_{k}`` belongs to the k-th BB bead of the CG structure.
    """
    runs = _bead_chain_runs(_cg_backbone_beads(cg_pdb))
    chain_of_k: Dict[int, int] = {}
    k = 0
    for ci, run in enumerate(runs):
        for _ in run:
            k += 1
            chain_of_k[k] = ci
    counts = [0] * len(runs)
    prefix = f"{moltype}_"
    with open(go_table) as fh:
        for line in fh:
            tokens = line.split()
            if len(tokens) < 2 or not tokens[0].startswith(prefix) or not tokens[1].startswith(prefix):
                continue
            for name in tokens[:2]:
                suffix = name[len(prefix):]
                ci = chain_of_k.get(int(suffix)) if suffix.isdigit() else None
                if ci is not None:
                    counts[ci] += 1
    labels = [f"{run[0][2]}{ci + 1}" for ci, run in enumerate(runs)]
    empty = [label for label, n in zip(labels, counts) if n == 0]
    if empty:
        raise TopologyGenerationError(
            message=(
                f"{moltype}: no Go contacts for chain(s) {empty} (contacts per chain: "
                f"{dict(zip(labels, counts))}). The contact map and the CG backbone do not match."
            ),
            error_code="TOP_MART_008",
            context={"moltype": moltype, "go_table": str(go_table)},
        )
    return dict(zip(labels, counts))


def check_backbone_continuity(itp_path: Union[str, Path], merge_pdb: Union[str, Path]) -> None:
    """Raise if a chain of a written col_N.itp is not one continuous BB-BB bonded backbone.

    Every pair of consecutive BB beads of a chain (TER records and chain identifiers of the merged
    CG PDB delimit the chains) must be joined by a bond or a constraint.
    """
    connected = set()
    section = None
    with open(itp_path) as fh:
        for line in fh:
            s = line.split(";")[0].strip()
            if not s or s.startswith("#"):
                continue
            if s.startswith("["):
                section = s.strip("[] ")
                continue
            if section in ("bonds", "constraints"):
                tokens = s.split()
                if len(tokens) >= 2 and tokens[0].isdigit() and tokens[1].isdigit():
                    connected.add(frozenset((int(tokens[0]), int(tokens[1]))))
    missing = []
    for run in _bead_chain_runs(_cg_backbone_beads(merge_pdb)):
        for a, b in zip(run, run[1:]):
            if frozenset((a[0], b[0])) not in connected:
                tag = " (crosslink residue)" if a[1] in CROSSLINK_RESNAMES or b[1] in CROSSLINK_RESNAMES else ""
                missing.append(f"chain {a[2]} {a[1]}{a[3]}-{b[1]}{b[3]} (BB {a[0]}-{b[0]}){tag}")
    if missing:
        raise TopologyGenerationError(
            message=(
                f"{Path(itp_path).name}: {len(missing)} chain(s) not continuous, missing backbone bond(s): "
                f"{'; '.join(missing[:10])}{' ...' if len(missing) > 10 else ''}"
            ),
            error_code="TOP_MART_008",
            context={"itp": str(itp_path), "merge_pdb": str(merge_pdb)},
        )


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

        Primary path: use System connectivity (model.connect).
        Fallback: if no connections are present in the System, infer connect IDs
        from files on disk matching '{model_id}.[0-9]*.CG.pdb'.
        """
        if cnt_model is None:
            LOG.error("cnt_model is None in merge_pdbs")
            return None

        output_file = f"{int(cnt_model)}.merge.pdb"
        merged_count = 0

        # Try to read connections from the System
        sys_connect_ids: Optional[list] = None
        if model_id is not None and self.system is not None:
            try:
                model = self.system.get_model(model_id=model_id)
                if model is not None and getattr(model, "connect", None):
                    sys_connect_ids = list(model.connect)
            except Exception as e:
                LOG.debug(f"Could not read connect list from System for model {model_id}: {e}")

        # Fallback: infer from files
        if not sys_connect_ids:
            inferred_files = sorted(Path().glob(f"{int(model_id)}.[0-9]*.CG.pdb")) if model_id is not None else []
            if not inferred_files:
                LOG.error(f"No model or connections for model_id {model_id} and no CG files found on disk")
                if os.path.exists(output_file):
                    os.remove(output_file)
                return None
            LOG.debug(f"Model {model_id}: no connections in System; merging from disk ({len(inferred_files)} CG files).")

            try:
                with open(output_file, "w") as f:
                    for cg_path in inferred_files:
                        with open(cg_path, "r") as infile:
                            lines_written = 0
                            for line in infile:
                                if line[0:6] in self.is_line or line.startswith("TER"):
                                    f.write(line)
                                    lines_written += 1
                            LOG.debug(f"Merged {lines_written} lines from {cg_path.name}")
                            merged_count += 1
                    f.write("END\n")
                LOG.debug(f"Successfully merged {merged_count} files into {output_file}")
                return output_file
            except Exception as e:
                LOG.error(f"Error merging PDBs (fallback): {e}")
                return None

        # Normal path using System connect IDs
        try:
            with open(output_file, "w") as f:
                for connect_id in sys_connect_ids:
                    input_file = f"{int(model_id)}.{int(connect_id)}.CG.pdb"
                    if not os.path.exists(input_file):
                        LOG.warning(f"CG PDB file not found: {input_file}")
                        continue
                    with open(input_file, "r") as infile:
                        lines_written = 0
                        for line in infile:
                            if line[0:6] in self.is_line or line.startswith("TER"):
                                f.write(line)
                                lines_written += 1
                        LOG.debug(f"Merged {lines_written} lines from {input_file}")
                        merged_count += 1
                f.write("END\n")

            if merged_count < len(sys_connect_ids):
                # A partial merge silently drops the failed connect_id's entire
                # triple helix from the output (e.g. a crosslink partner whose
                # martinize2 run failed) while still looking like a complete,
                # successful group. Treat any missing member as a full failure
                # so it surfaces as one instead of shipping an incomplete topology.
                LOG.error(
                    f"Incomplete merge for model {model_id}: {merged_count}/"
                    f"{len(sys_connect_ids)} connection(s) succeeded"
                )
                if os.path.exists(output_file):
                    os.remove(output_file)
                return None

            LOG.debug(f"Successfully merged {merged_count} files into {output_file}")
            if sys_connect_ids and len(sys_connect_ids) > 1:
                LOG.debug(f"Model {model_id} merge contains {merged_count} structures from connections: {sys_connect_ids}")
            return output_file
        except Exception as e:
            LOG.error(f"Error merging PDBs: {e}")
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
                            pdb = [line for line in file if line[0:6] in self.is_line or line.startswith("TER")]
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
                        pdb = [line for line in file if line[0:6] in self.is_line or line.startswith("TER")]
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

    def backbone_conect(self, pdb: Optional[List[str]] = None) -> Tuple[List[str], List[str]]:
        """
        CONECT records for the peptide bonds of the input, so that martinize2 does not guess them.

        martinize2 (vermouth MakeBonds) bonds residues by distance only. At crosslink sites the
        all-atom peptide C-N can be stretched beyond its limit (2.3 A at -bonds-fudge 1.4; the
        N-terminal PYD LY3-SER bond of chain C measures 3.1 A): the chain is then cut in two
        fragments, and vermouth orders the fragments by their lowest residue, so a fragment
        attached by a steric contact to an earlier chain comes first and the N-terminal piece
        is written after the end of its chain. Listing every peptide bond as CONECT keeps each
        chain one fragment in input order; bonds that vermouth finds anyway are unchanged.

        Only residues adjacent in the file, in the same chain (no TER in between) and with
        consecutive residue numbers are connected, and only up to MAX_PEPTIDE_CN_ANGSTROM.

        Parameters
        ----------
        pdb : Optional[List[str]]
            PDB lines of the model with their original residue numbers (set_pdb only renumbers
            residues, so the atom serials are those of the martinize2 input).

        Returns
        -------
        Tuple[List[str], List[str]]
            CONECT lines and a description of every bond longer than vermouth's guess limit.
        """
        if not pdb:
            return [], []
        residues: List[Dict[str, Any]] = []
        cur = None
        for line in pdb:
            if line.startswith("TER"):
                cur = None
                residues.append(None)
                continue
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            key = (line[21:22], line[22:27])
            if cur is None or cur["key"] != key:
                cur = {"key": key, "resname": line[17:20].strip(), "atoms": {}}
                residues.append(cur)
            try:
                cur["atoms"][line[12:16].strip()] = (
                    int(line[6:11]),
                    (float(line[30:38]), float(line[38:46]), float(line[46:54])),
                )
            except ValueError:
                cur["atoms"][line[12:16].strip()] = None

        conect, strained = [], []
        for prev, nxt in zip(residues, residues[1:]):
            if prev is None or nxt is None or prev["key"][0] != nxt["key"][0]:
                continue
            try:
                consecutive = int(nxt["key"][1][:4]) == int(prev["key"][1][:4]) + 1
            except ValueError:
                consecutive = False
            if not consecutive or prev["key"][1][4:].strip() or nxt["key"][1][4:].strip():
                continue
            c_atom, n_atom = prev["atoms"].get("C"), nxt["atoms"].get("N")
            if not c_atom or not n_atom:
                continue
            dist = sum((a - b) ** 2 for a, b in zip(c_atom[1], n_atom[1])) ** 0.5
            label = (
                f"chain {prev['key'][0]} {prev['resname']} {prev['key'][1].strip()}-"
                f"{nxt['resname']} {nxt['key'][1].strip()} (C-N {dist:.2f} A)"
            )
            if dist > MAX_PEPTIDE_CN_ANGSTROM:
                LOG.warning(f"Peptide C-N of {label} is longer than {MAX_PEPTIDE_CN_ANGSTROM} A; not connected")
                continue
            conect.append(f"CONECT{c_atom[0]:5d}{n_atom[0]:5d}\n")
            if dist > VERMOUTH_CN_GUESS_LIMIT:
                strained.append(label)
        return conect, strained

    def cap_pdb(self, pdb: Optional[List[str]] = None) -> Tuple[List[str], str, str]:
        """
        Decide martinize2 -nter/-cter flags robustly.
        1) Rename terminal ALA -> CLA in the PDB *before* flag decisions.
        2) If any chain starts/ends on a special block, use 'none' for that side.
        Special blocks include ACE/CLA and crosslink blocks
        (LY2/LY3/L4Y/L5Y/LYX/AGS/LGX/LZS/LZD), and also NME (since it is
        already an explicit cap).
        """
        if not pdb:
            LOG.warning("Empty PDB provided to cap_pdb")
            return [], "none", "none"

        try:
            # 1) Rename terminal ALA -> CLA
            chain_length = self.get_chain_length(pdb)  # e.g. {"A": " 178", "B": "...", "C": "..."}
            for i in range(len(pdb)):
                line = pdb[i]
                if not line.startswith(("ATOM  ", "HETATM")):
                    continue
                if line[17:20] == "ALA":
                    tag = line[21:26]  # e.g. "A 178"
                    if tag in {f"A{chain_length['A']}", f"B{chain_length['B']}", f"C{chain_length['C']}"}:
                        pdb[i] = line[:17] + "CLA " + line[21:]

            # 2) Inspect first/last residue per chain *after* the rename above
            first_res = {}   # chain -> resname
            last_res  = {}   # chain -> (resid, resname)
            last_oxt  = {}   # chain -> whether the current last residue has an OXT atom
            for line in pdb:
                if not line.startswith(("ATOM  ", "HETATM")):
                    continue
                ch = line[21:22]
                if ch not in self.is_chain:
                    continue
                resname = line[17:20].strip()
                atomname = line[12:16].strip()
                try:
                    resid = int(line[22:26])
                except ValueError:
                    continue
                if ch not in first_res:
                    first_res[ch] = resname
                if ch not in last_res or resid > last_res[ch][0]:
                    last_res[ch] = (resid, resname)
                    last_oxt[ch] = atomname == "OXT"
                elif resid == last_res[ch][0] and atomname == "OXT":
                    last_oxt[ch] = True

            firsts = set(first_res.values()) if first_res else set()
            lasts  = {name for _, name in last_res.values()} if last_res else set()
            any_last_has_oxt = any(last_oxt.values())

            # Residues meaning "this terminus is already capped / is a crosslink
            # block", so martinize2 must NOT add another terminal modification.
            special_first = {
                "ACE", "CLA", "LY2", "LY3", "L4Y", "L5Y", "LYX", "AGS", "LGX",
                "LZS", "LZD"
            }
            special_last = {
                "ACE", "CLA", "LY2", "LY3", "L4Y", "L5Y", "LYX", "AGS", "LGX",
                "LZS", "LZD", "NME"
            }

            # N-terminus: the geometry stage (caps.py) pre-caps every chain with an
            # ACE residue, so the first residue is normally ACE and needs no -nter
            # modification. We always use 'none' (martinize2's ACE modification is
            # also known to fail on some first residues, e.g. GLN). Warn if a chain
            # is unexpectedly NOT pre-capped so it isn't silently left uncapped.
            if firsts and not (firsts & special_first):
                LOG.warning(
                    f"N-terminus not pre-capped (first residues={sorted(firsts)}); "
                    f"using -nter none. Check the geometry capping stage if unexpected."
                )
            nter_flag = "none"

            # C-terminus: apply the NME cap unless the chain already ends in a
            # cap/crosslink block, or its last residue already carries its own
            # OXT atom -- a complete, naturally-terminated carboxylic acid.
            # Applying NME on top of an existing OXT is a modification conflict
            # that crashes vermouth's modification-patching (ValueError: Cannot
            # apply modification to block), confirmed by direct reproduction:
            # martinize2 succeeds on the identical input once -cter is 'none'
            # instead of 'NME' for such a residue.
            cter_flag = (
                "none"
                if (not lasts or (lasts & special_last) or any_last_has_oxt)
                else "NME"
            )

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

            # Find each chain's last residue directly, by keeping the resid of
            # the most recently seen ATOM/HETATM line per chain (atoms are in
            # increasing resid order within a chain, so the last one seen is
            # the chain's last residue). This must not rely on line adjacency
            # (comparing pdb[i] to pdb[i+1]) -- real per-model PDBs carry a TER
            # record after every chain, which breaks any adjacency-based
            # transition detection and previously left chain_length empty for
            # every chain, silently disabling the ALA -> CLA rename below.
            for line in pdb:
                if not line.startswith(("ATOM  ", "HETATM")):
                    continue
                ch = line[21:22]
                if ch in chain_length:
                    chain_length[ch] = line[22:26]

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
        self,
        topology_file: str = "system.top",
        size: Optional[int] = None,
        processed_ids: Optional[List[int]] = None,
        glucosepane_timestep: Optional[int] = None,
    ) -> None:
        """
        Write final topology file for the system.

        Parameters
        ----------
        topology_file : str
            Path to the output topology file
        size : Optional[int]
            Number of models considered (fallback range when processed_ids is None).
        processed_ids : Optional[List[int]]
            col_{id} indices that actually produced an .itp. When given, only these
            are #included, preventing GROMACS errors from missing files when some
            models failed or were skipped.
        glucosepane_timestep : Optional[int]
            Glucosepane Martini 3 parameter set used (20 or 10 fs). If given and an included
            col_{id}.itp contains glucosepane (AGS), a comment with the set and the time step it
            is meant for is written at the top of the file.
        """
        if processed_ids is not None:
            ids_for_includes = sorted({int(m) for m in processed_ids})
        elif size:
            ids_for_includes = list(range(size))
        else:
            LOG.warning("No size or processed_ids provided to write_system_topology")
            return

        try:
            with open(topology_file, "w") as f:
                f.write("; This is the topology for the collagen microfibril\n")
                if glucosepane_timestep is not None and any(
                    Path(f"col_{int(m)}.itp").exists()
                    and re.search(r"^\s*\d+\s+\S+\s+\d+\s+AGS\s", Path(f"col_{int(m)}.itp").read_text(), re.MULTILINE)
                    for m in ids_for_includes
                ):
                    kind = "model A_b2, two virtual sites" if int(glucosepane_timestep) == 20 else "no virtual sites"
                    f.write(
                        f"; Glucosepane (LGX/AGS): Martini 3 parameter set for a {int(glucosepane_timestep)} fs "
                        f"time step ({kind}), glucosepane_martini_timestep: {int(glucosepane_timestep)}\n"
                    )
                f.write("#define GO_VIRT\n")
                f.write('#include "martini_v3.0.0.itp"\n')
                f.write('#include "go-sites.itp"\n\n')

                included = []
                for m in ids_for_includes:
                    # Belt-and-suspenders: never include a file that isn't on disk.
                    if not Path(f"col_{int(m)}.itp").exists():
                        LOG.warning(f"Skipping missing topology include: col_{int(m)}.itp")
                        continue
                    f.write(f'#include "col_{int(m)}.itp"\n')
                    f.write(f'#include "col_{int(m)}_go-excl.itp"\n')
                    included.append(int(m))

                f.write('\n#include "martini_v3.0.0_solvents_v1.itp"\n')
                f.write('#include "martini_v3.0.0_ions_v1.itp"\n')

                f.write("\n[ system ]\n")
                f.write("Collagen, Martini 3 and Go-Potentials \n")
                f.write("\n[ molecules ]\n")

                for t in included:
                    f.write(f"col_{int(t)}     1\n")

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
            # Only ATOM/HETATM carry coordinates. TER/ANISOU lines have no
            # parseable x/y/z, so don't try (and don't emit a spurious
            # "Bad coord formatting" warning for every TER record).
            if line[0:6] in ("ATOM  ", "HETATM"):
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


def _build_connected_groups(system: System) -> List[List[Any]]:
    """
    Collapse crosslink-connected models into transitive groups.

    Mirrors Amber.get_connected_groups so that each physically connected set
    of triple helices is emitted exactly once. model.connect is symmetric and
    self-inclusive (e.g. {self} + direct neighbors); a BFS gives the full group.
    """
    all_models = list(system.get_models())
    processed: set = set()
    groups: List[List[Any]] = []

    for model_id in all_models:
        if model_id in processed:
            continue

        model = system.get_model(model_id=model_id)
        if not model or not getattr(model, "connect", None):
            groups.append([model_id])
            processed.add(model_id)
            continue

        group: set = set()
        stack = [model_id]
        while stack:
            cur = stack.pop()
            if cur in group:
                continue
            group.add(cur)
            cur_model = system.get_model(model_id=cur)
            if cur_model and getattr(cur_model, "connect", None):
                for nb in cur_model.connect:
                    if nb not in group:
                        stack.append(nb)

        groups.append(sorted(group))
        processed.update(group)

    return groups


def check_glucosepane_parameter_set(itp_ids: List[int], glucosepane_timestep: int) -> None:
    """Raise if the AGS beads of the written col_{id}.itp do not belong to the requested glucosepane set."""
    want_virtual = int(glucosepane_timestep) == 20
    for m in itp_ids or []:
        path = Path(f"col_{int(m)}.itp")
        if not path.exists():
            continue
        section = None
        for line in path.read_text().splitlines():
            s = line.split(";")[0].strip()
            if s.startswith("["):
                section = s.strip("[] ")
                continue
            p = s.split()
            if section != "atoms" or len(p) < 5 or p[3] != "AGS" or p[4] not in ("SC3", "SC6"):
                continue
            mass = float(p[7]) if len(p) > 7 else None
            is_virtual = mass == 0.0
            if is_virtual != want_virtual:
                raise TopologyGenerationError(
                    message=(
                        f"{path}: AGS {p[4]} has mass {mass}, which does not belong to the requested "
                        f"glucosepane_martini_timestep {int(glucosepane_timestep)} set; the installed martini300C "
                        "force field may have been changed by another run in the same environment"
                    ),
                    error_code="TOP_MART_007",
                    context={"itp": str(path), "glucosepane_martini_timestep": int(glucosepane_timestep)},
                )


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
        processed_models = []             # fibril model_ids that succeeded — for counts/logging ONLY
                                          # (col_N.itp file indices live in written_itp_ids)
        written_itp_ids: List[int] = []   # col_{cnt_model} indices actually written

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
        LOG.info(f"{Fore.BLUE}Building coarse-grained topology:{Style.RESET_ALL}")

        if len(list(system.get_models())) > 0:
            first_model = system.get_model(model_id=list(system.get_models())[0])
            model_type = first_model.type
            type_dir = topology_dir / model_type
            type_dir.mkdir(exist_ok=True, parents=True)

        models_list = [model_id for model_id in system.get_models()]

        # Group crosslink-connected models so each physical triple helix is
        # written exactly once. Without this, every model re-merges its whole
        # neighborhood and connected helices are emitted multiple times,
        # producing duplicate beads (0 A overlaps) and inflated particle counts.
        connected_groups = _build_connected_groups(system)
        LOG.debug(
            f"Martini: collapsed {len(models_list)} models into "
            f"{len(connected_groups)} connected group(s)"
        )

        model_status = {}
        failed_martinize, failed_contact, failed_go, failed_itp = [], [], [], []
        processed_in_topology: set = set()

        for group in tqdm(connected_groups, desc="Building topology", unit="grp"):
            # Representative model id drives all file naming for this group.
            model_id = group[0]
            model = system.get_model(model_id=model_id)
            if model is None:
                LOG.warning(f"Skipping group {group}: representative model not found in System")
                model_status[model_id] = "no_model"
                continue

            # Process every member of the group exactly once, and make sure the
            # downstream helpers that read connectivity from the System
            # (merge_pdbs, Itp.allocate/read_model) see the full group.
            connect_ids = list(group)
            model.connect = list(group)

            if not connect_ids:
                LOG.warning(f"Skipping group {group}: no members")
                model_status[model_id] = "no_connections"
                continue

            try:
                for connect_id in connect_ids:
                    try:
                        pdb = martini.read_pdb(pdb_id=connect_id)
                        if not pdb:
                            LOG.warning(f"Empty PDB for connect_id: {connect_id}")
                            continue

                        pdb = martini.translate_pdb(pdb=pdb)
                        pdb_capped, cter, nter = martini.cap_pdb(pdb=pdb)
                        order, map_pdb = martini.set_pdb(pdb=pdb_capped)
                        conect, strained = martini.backbone_conect(pdb=pdb_capped)
                        if strained:
                            LOG.info(
                                f"Model {model_id} -> connection {connect_id}: peptide bond(s) longer than "
                                f"martinize2's distance guess passed explicitly: {'; '.join(strained)}"
                            )

                        martini.write_pdb(pdb=order + conect, file="tmp.pdb")
                        martini.write_pdb(pdb=map_pdb, file="map.pdb")

                        def _martinize_cmd(nter_flag: str, cter_flag: str) -> str:
                            args = (
                                f"-f tmp.pdb -sep -merge A,B,C "
                                f"-collagen -from amber99 -o topol.top -bonds-fudge {MARTINIZE2_BONDS_FUDGE} -p backbone "
                                f"-ff {martini.ff}00C -x {int(model_id)}.{int(connect_id)}.CG.pdb "
                                f"-nter {nter_flag} -cter {cter_flag} "
                                f"-govs-include -govs-moltype col_{int(model_id)}.{int(connect_id)} -maxwarn 4"
                            )
                            # An explicit user override wins.
                            configured = getattr(config, "martinize2_command", None)
                            if configured and configured not in ("martinize2", "auto"):
                                return f"{configured} {args}"
                            # Prefer the martinize2 that lives in the SAME environment
                            # as the running interpreter: that is the vermouth into
                            # which the custom martini300C force field was just
                            # installed. This avoids picking a martinize2 from another
                            # Python (e.g. a pyenv shim ahead on PATH) whose vermouth
                            # does not know the custom force field ("Unknown force
                            # field martini300C").
                            local_exe = Path(sys.executable).parent / "martinize2"
                            if local_exe.exists():
                                return f"{local_exe} {args}"
                            # Otherwise fall back to a direct martinize2 on PATH, and
                            # only wrap in `conda run` when it is reachable solely via
                            # a conda env.
                            try:
                                exe, use_conda, _env = find_martinize2_executable()
                            except Exception:
                                exe, use_conda = "martinize2", True
                            if use_conda:
                                return get_conda_command_with_path("martinize2", args)
                            return f"{exe} {args}"

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
                            check_cg_residue_order(order, cg_file, f"col_{int(model_id)}.{int(connect_id)}")

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
                            f"{sys.executable} create_goVirt.py -s {int(model_id)}.{int(connect_id)}.CG.pdb "
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
                        check_go_contacts_per_chain(
                            f"col_{int(model_id)}.{int(connect_id)}_go-table.itp",
                            cg_file,
                            f"col_{int(model_id)}.{int(connect_id)}",
                        )

                    except TopologyGenerationError:
                        raise
                    except Exception as e:
                        LOG.error(f"Error processing model {model_id} -> connection {connect_id}: {str(e)}")
                        continue

                merged_pdb = martini.merge_pdbs(model_id=int(model_id), cnt_model=cnt_model)
                if merged_pdb:
                    model_status[model_id] = "merged"
                    try:
                        itp_ = Itp(
                            system=system,
                            model_id=int(model_id),
                            glucosepane_timestep=config.glucosepane_martini_timestep,
                        )
                        itp_.read_model(model_id=int(model_id))
                        itp_.go_to_pairs(model_id=int(model_id))
                        itp_.make_topology(model_id=int(model_id), cnt_model=cnt_model)

                        # martinize2 infers backbone bonds by distance and drops
                        # the BB-BB peptide bond wherever crosslink strain pushes
                        # the all-atom C-N beyond its ~1.9 A threshold, fragmenting
                        # the CG chain. Add back the missing consecutive-residue
                        # backbone bonds (the amber path is unaffected because
                        # pdb2gmx bonds from templates, not geometry).
                        try:
                            repair_backbone_bonds(
                                itp_path=f"col_{int(cnt_model)}.itp",
                                merge_pdb=f"{int(cnt_model)}.merge.pdb",
                            )
                        except Exception as e:
                            LOG.warning(
                                f"Backbone-bond repair failed for col_{int(cnt_model)}.itp: {e}"
                            )
                        check_backbone_continuity(
                            itp_path=f"col_{int(cnt_model)}.itp",
                            merge_pdb=f"{int(cnt_model)}.merge.pdb",
                        )

                        processed_models.append(model_id)
                        written_itp_ids.append(cnt_model)   # this col_{cnt_model}.itp now exists

                        # Track what we processed
                        for connect_id in connect_ids:
                            processed_in_topology.add(connect_id)

                    except TopologyGenerationError:
                        raise
                    except Exception as e:
                        LOG.error(f"Error processing ITP for model {model_id}: {str(e)}")
                        failed_itp.append(model_id)
                else:
                    model_status[model_id] = "no_merge_file"
                    LOG.error(f"Model {model_id} failed: No merged PDB created")

                cnt_model += 1

            except TopologyGenerationError:
                raise
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

        # The glucosepane blocks come from the martini300C force field installed in vermouth; check that
        # the written topologies carry the requested set (R2 = AGS:SC6 and R5 = AGS:SC3 are massless
        # virtual sites in the 20 fs set and massive beads in the 10 fs set), e.g. against an install
        # changed by another run in the same environment.
        check_glucosepane_parameter_set(written_itp_ids, config.glucosepane_martini_timestep)

        try:
            final_topology_file = f"collagen_fibril_{config.species}.top"
            martini.write_system_topology(
                topology_file=final_topology_file,
                size=cnt_model,
                processed_ids=written_itp_ids,
                glucosepane_timestep=config.glucosepane_martini_timestep,
            )

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
            # Gate intermediate-file cleanup on topology_debug (the flag meant for
            # keeping topology intermediates), not the general logging `debug` flag.
            if not getattr(config, "topology_debug", False):
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
