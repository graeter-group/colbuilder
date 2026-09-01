"""
Amber topology generation module.

This module provides functionality for generating AMBER99 force field topology files
for molecular systems, particularly focused on collagen microfibrils with proper
crosslink handling between models.
"""

import os
import shutil
from pathlib import Path
import asyncio
import numpy as np
from colorama import Fore, Style
from typing import List, Any, Optional, Dict, Union, Tuple, Set
import re

from colbuilder.core.geometry.system import System
from colbuilder.core.geometry.crosslink import read_crosslink, Crosslink
from colbuilder.core.utils.dec import timeit
from colbuilder.core.utils.files import FileManager
from colbuilder.core.utils.config import ColbuilderConfig
from colbuilder.core.utils.exceptions import TopologyGenerationError
from colbuilder.core.utils.logger import setup_logger

LOG = setup_logger(__name__)

REQUIRED_FF_FILES = ['residuetypes.dat', 'specbond.dat']

class Amber:
    """AMBER99 force field topology generator with crosslink handling."""

    def __init__(self, system: Optional[System] = None, ff: Optional[str] = None) -> None:
        self.system = system
        self.ff = ff + '.ff' if ff else None
        self.pdb_line_types = ('ATOM  ', 'HETATM', 'ANISOU', 'TER   ')

    def get_connected_groups(self) -> List[List[int]]:
        """Group models that are connected together."""
        if not self.system:
            return []

        all_models = list(self.system.get_models())
        processed = set()
        groups = []

        for model_id in all_models:
            if model_id in processed:
                continue

            model = self.system.get_model(model_id=model_id)
            if not model or not model.connect:
                groups.append([model_id])
                processed.add(model_id)
                continue

            group = set()
            to_process = [model_id]

            while to_process:
                current_id = to_process.pop()
                if current_id in group:
                    continue

                group.add(current_id)
                current_model = self.system.get_model(model_id=current_id)

                if current_model and current_model.connect:
                    for connected_id in current_model.connect:
                        if connected_id not in group:
                            to_process.append(connected_id)

            if group:
                groups.append(sorted(list(group)))
                processed.update(group)

        return groups

    def merge_connected_models(
        self,
        model_group: List[int]
    ) -> Optional[Tuple[str, str, List[Tuple[Crosslink, Crosslink]]]]:
        """
        Merge connected models and detect crosslink pairs from individual model files.

        Returns:
            Tuple of (model_type, group_id, crosslink_pairs) or None
        """
        if not model_group or not self.system:
            return None

        first_model = self.system.get_model(model_id=float(model_group[0]))
        if not first_model or not first_model.type:
            return None

        model_type = first_model.type
        os.makedirs(model_type, exist_ok=True)

        group_id = "_".join(str(int(mid)) for mid in sorted(model_group))
        output_file = os.path.join(model_type, f"{group_id}.merge.pdb")

        # DETECT CROSSLINKS FROM INDIVIDUAL MODEL FILES BEFORE MERGING
        crosslink_pairs = []
        if len(model_group) > 1:
            all_crosslinks = []

            for mid in model_group:
                caps_file = os.path.join(model_type, f"{int(mid)}.caps.pdb")
                if os.path.exists(caps_file):
                    try:
                        cls = read_crosslink(caps_file)
                        for cl in cls:
                            cl.model_id = int(mid) # tag
                            all_crosslinks.append(cl)
                        LOG.debug(f"Found {len(cls)} crosslinks in model {int(mid)}")
                    except Exception as e:
                        LOG.warning(f"Could not read crosslinks from {caps_file}: {e}")

            candidate_pairs = []
            for i, cl1 in enumerate(all_crosslinks):
                for j in range(i + 1, len(all_crosslinks)):
                    cl2 = all_crosslinks[j]
                    if cl1.model_id == cl2.model_id:
                        continue
                    if not self._are_compatible_crosslinks(cl1, cl2):
                        continue
                    distance = float(np.linalg.norm(cl1.position - cl2.position))
                    if distance < 5.0:
                        candidate_pairs.append((distance, i, j))

            candidate_pairs.sort(key=lambda t: t[0])
            used_markers: Set[int] = set()
            for distance, i, j in candidate_pairs:
                if i in used_markers or j in used_markers:
                    continue
                used_markers.add(i)
                used_markers.add(j)
                cl1, cl2 = all_crosslinks[i], all_crosslinks[j]
                crosslink_pairs.append((cl1, cl2))
                LOG.debug(
                    f"  Crosslink pair: model {cl1.model_id} ({cl1.resname}{getattr(cl1, 'atom', '')}) - "
                    f"model {cl2.model_id} ({cl2.resname}{getattr(cl2, 'atom', '')}), distance: {distance:.2f} Å"
                )

            if crosslink_pairs:
                LOG.info(f"Found {len(crosslink_pairs)} crosslink pairs for group {group_id}")

        def write_caps(path, out):
            if os.path.exists(path):
                with open(path, "r") as f_in:
                    for line in f_in:
                        if line.startswith(self.pdb_line_types):
                            out.write(line)
            else:
                LOG.debug(f"Caps file not found: {path}")

        with open(output_file, "w") as f_out:
            for mid in model_group:
                caps = os.path.join(model_type, f"{int(mid)}.caps.pdb")
                write_caps(caps, f_out)

            for mid in model_group:
                model = self.system.get_model(model_id=float(mid))
                if model and model.connect:
                    for cid in model.connect:
                        if int(cid) not in model_group:
                            caps = os.path.join(model_type, f"{int(cid)}.caps.pdb")
                            write_caps(caps, f_out)

            f_out.write("END\n")

        if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
            return (model_type, group_id, crosslink_pairs)

        LOG.error(f"Failed to create merged PDB for group {model_group}")
        return None

    def _are_compatible_crosslinks(self, cl1: Crosslink, cl2: Crosslink) -> bool:
        """
        Check if two crosslinks are compatible for bonding.

        Args:
            cl1: First crosslink
            cl2: Second crosslink

        Returns:
            True if crosslinks can form a bond
        """
        # Divalent crosslinks (D-type)
        divalent_aldehyde = {'L4Y', 'L4X', 'LY4', 'LX4', 'LGX', 'LPS', 'LZD'}
        divalent_amine = {'L5Y', 'L5X', 'LY5', 'LX5', 'AGS', 'APD', 'LZS'}

        # Trivalent crosslinks (T-type). Each aldehyde (ring) residue bonds two
        # specific arms — one via ring carbon C13, one via C12 — fixed per type
        # (data/sequence/crosslinks.csv). Blocks cross-arm and cross-type bonds.
        trivalent_bond_map = {
            "LYX": {"C13": "LY3", "C12": "LY2"},  # PYD
            "LXX": {"C13": "LX3", "C12": "LX2"},  # DPD
            "LXY": {"C13": "L3Y", "C12": "L2Y"},  # PYL
            "LYY": {"C13": "L3X", "C12": "L2X"},  # DPL
        }

        # Divalent pairs: aldehyde + amine
        if cl1.resname in divalent_aldehyde and cl2.resname in divalent_amine:
            return True
        if cl1.resname in divalent_amine and cl2.resname in divalent_aldehyde:
            return True

        # Trivalent pairs: the aldehyde ring carbon must bond its designated arm.
        def _trivalent_ok(aldehyde: Crosslink, amine: Crosslink) -> bool:
            arms = trivalent_bond_map.get(aldehyde.resname)
            if arms is None:
                return False
            atom = (getattr(aldehyde, "atom", "") or "").strip()
            if atom in arms:
                return arms[atom] == amine.resname
            # Atom name unknown: restrict to this central's two valid arms (still
            # blocks cross-type; the 1-to-1 match below resolves C13 vs C12).
            return amine.resname in arms.values()

        if _trivalent_ok(cl1, cl2):
            return True
        if _trivalent_ok(cl2, cl1):
            return True

        return False

    def find_atom_indices_for_crosslinks(
        self,
        itp_file: str,
        crosslink_pairs: List[Tuple[Crosslink, Crosslink]],
        merged_pdb_file: str,
    ) -> List[Tuple[Tuple[Crosslink, Crosslink], Tuple[int, int]]]:
        """Map crosslink markers to atom indices, skipping invalid same-side T pairs."""
        mapped: List[Tuple[Tuple[Crosslink, Crosslink], Tuple[int, int]]] = []

        with open(itp_file, 'r') as f:
            lines = f.readlines()

        atoms_section = False
        atom_data = []
        for line in lines:
            if line.strip().startswith('[ atoms ]'):
                atoms_section = True
                continue
            elif atoms_section and line.strip().startswith('['):
                break
            elif atoms_section and not line.strip().startswith(';') and line.strip():
                parts = line.split()
                if len(parts) >= 8:
                    atom_data.append({
                        'index': int(parts[0]),
                        'atom_name': parts[4],
                        'residue_nr': int(parts[2]),
                        'residue_name': parts[3],
                    })

        atom_by_index = {a['index']: a for a in atom_data}

        # Read coordinates from GRO file
        gro_file = str(Path(itp_file).with_suffix('.gro'))
        gro_coords: Dict[int, np.ndarray] = {}
        if os.path.exists(gro_file):
            with open(gro_file, 'r') as f:
                gro_lines = f.readlines()
            for line in gro_lines[2:-1]:
                if len(line) >= 44:
                    try:
                        atom_idx = int(line[15:20].strip())
                        x = float(line[20:28]) * 10  # nm to Å
                        y = float(line[28:36]) * 10
                        z = float(line[36:44]) * 10
                        gro_coords[atom_idx] = np.array([x, y, z])
                    except (ValueError, IndexError):
                        continue

        # Map each crosslink pair to atom indices
        for cl1, cl2 in crosslink_pairs:
            best_atom1_idx = None
            best_atom2_idx = None
            best_dist1 = float('inf')
            best_dist2 = float('inf')

            for atom in atom_data:
                atom_idx = atom['index']

                # Check if this atom matches cl1
                if (atom['residue_name'] == cl1.resname and
                    self._is_crosslink_atom(atom['atom_name'], cl1.resname, cl1.type)):
                    if atom_idx in gro_coords:
                        d = np.linalg.norm(gro_coords[atom_idx] - cl1.position)
                        if d < best_dist1 and d < 5.0:
                            best_dist1 = d
                            best_atom1_idx = atom_idx

                # Check if this atom matches cl2
                if (atom['residue_name'] == cl2.resname and
                    self._is_crosslink_atom(atom['atom_name'], cl2.resname, cl2.type)):
                    if atom_idx in gro_coords:
                        d = np.linalg.norm(gro_coords[atom_idx] - cl2.position)
                        if d < best_dist2 and d < 5.0:
                            best_dist2 = d
                            best_atom2_idx = atom_idx

            # Reject invalid same-side T–T selections
            if best_atom1_idx and best_atom2_idx:
                a1 = atom_by_index.get(best_atom1_idx)
                a2 = atom_by_index.get(best_atom2_idx)
                if cl1.type == "T" and cl2.type == "T" and a1 and a2:
                    a1_isA = self._is_crosslink_atom(a1['atom_name'], a1['residue_name'], "T", "A")
                    a1_isB = self._is_crosslink_atom(a1['atom_name'], a1['residue_name'], "T", "B")
                    a2_isA = self._is_crosslink_atom(a2['atom_name'], a2['residue_name'], "T", "A")
                    a2_isB = self._is_crosslink_atom(a2['atom_name'], a2['residue_name'], "T", "B")
                    if (a1_isA and a2_isA) or (a1_isB and a2_isB):
                        LOG.debug(
                            f"Skipping same-side T pair "
                            f"{a1['residue_name']}{a1['residue_nr']}–{a2['residue_name']}{a2['residue_nr']}"
                        )
                        continue

                mapped.append(((cl1, cl2), (best_atom1_idx, best_atom2_idx)))
            else:
                LOG.warning(
                    f"Could not find atoms for crosslink pair: "
                    f"{cl1.resname}{cl1.resid} - {cl2.resname}{cl2.resid}"
                )

        return mapped

    def _is_crosslink_atom(
        self,
        atom_name: str,
        resname: str,
        crosslink_type: str,
        want_side: Optional[str] = None,
    ) -> bool:
        """
        Check if this atom is the target bonding atom for the crosslink type.
        For T-type crosslinks we distinguish sides to block same-side bonds:
        - Side "A": aldehyde side (e.g., LYX carbonyl C13/C12)
        - Side "B": amine side (e.g., LY3 CG or LY2 CB)
        """
        if crosslink_type == "T":
            is_A = (resname in ("LYX", "LXY", "LYY", "LXX") and atom_name in ("C13", "C12"))
            is_B = (
                (resname in ("LY3", "LX3", "L3Y", "L2Y", "L3X", "L2X") and atom_name == "CG")
                or (resname in ("LY2", "LX2") and atom_name == "CB")
            )
            if want_side == "A":
                return is_A
            if want_side == "B":
                return is_B
            return is_A or is_B

        elif crosslink_type == "D":
            if resname in ("L4Y", "L4X", "LY4", "LX4") and atom_name == "CE":
                return True
            elif resname in ("L5Y", "L5X", "LY5", "LX5") and atom_name == "NZ":
                return True
            elif resname in ("LGX", "LPS") and atom_name == "CE":
                return True
            elif resname in ("AGS", "APD") and atom_name == "NZ":
                return True
            elif resname in ("LZD") and atom_name == "CE":
                return True
            elif resname in ("LZS") and atom_name == "NZ1":
                return True

        return False

    def add_crosslink_topology_to_itp(
        self,
        itp_file: str,
        crosslink_pairs: List[Tuple[Crosslink, Crosslink]],
        merged_pdb_file: str
    ) -> None:
        """Add complete crosslink topology (bonds, angles, dihedrals) to an existing ITP file."""
        if not crosslink_pairs:
            return

        LOG.debug(f"        Adding crosslink topology for {len(crosslink_pairs)} crosslink pairs to {itp_file}")

        try:
            mapped_pairs = self.find_atom_indices_for_crosslinks(itp_file, crosslink_pairs, merged_pdb_file)

            # Surface mapping failures: a dropped pair means NO bond and NO
            # exclusions for that crosslink -> close contacts that crash MD.
            total = len(crosslink_pairs)
            n_mapped = len(mapped_pairs)
            if n_mapped < total:
                LOG.error(
                    f"Crosslink mapping incomplete in {os.path.basename(itp_file)}: "
                    f"{n_mapped}/{total} pair(s) mapped, {total - n_mapped} dropped. "
                    f"Dropped crosslinks get no bond and no exclusions (likely clashes "
                    f"during equilibration). Check residue/atom-name coverage in "
                    f"_is_crosslink_atom and the 5 A marker cutoff."
                )

            if not mapped_pairs:
                LOG.warning(f"No atom indices found for crosslink topology in {itp_file}")
                return

            valid_bonds: Set[Tuple[int, int]] = set()
            valid_bond_data: List[Dict] = []

            for (cl1, cl2), (atom1_idx, atom2_idx) in mapped_pairs:
                if atom1_idx == atom2_idx:
                    continue
                bond_key = tuple(sorted((atom1_idx, atom2_idx)))
                if bond_key not in valid_bonds:
                    valid_bonds.add(bond_key)
                    valid_bond_data.append({
                        'atoms': bond_key,
                        'cl1': cl1,
                        'cl2': cl2,
                        'original_indices': (atom1_idx, atom2_idx),
                    })

            if not valid_bond_data:
                LOG.warning(f"No valid crosslink topology to add to {itp_file}")
                return

            # Add the crosslink bonds, then complete the bonded topology around
            # them in a single graph pass (angles, proper dihedrals, 1-4 pairs).
            self._add_crosslink_bonds(itp_file, valid_bond_data)

            try:
                self._complete_crosslink_topology(itp_file, valid_bond_data)
            except Exception as e:
                LOG.warning(f"Failed to complete crosslink angles/dihedrals/pairs, "
                            f"but bonds were added: {str(e)}")

            # No explicit [ exclusions ] are written: grompp derives all nonbonded
            # exclusions from the final bond graph of the moleculetype (nrexcl=3),
            # including the crosslink bonds added above, so an explicit block would
            # only duplicate what grompp already does.

            LOG.debug(f"    Successfully added crosslink topology to {itp_file}")

        except Exception as e:
            LOG.error(f"Failed to add crosslink topology to {itp_file}: {str(e)}")

    def _add_crosslink_bonds(self, itp_file: str, valid_bond_data: List[Dict]) -> None:
        """Add crosslink bonds, skipping any pdb2gmx already wrote via specbond.dat
        (adding it twice would double-count the bond)."""
        with open(itp_file, 'r') as f:
            lines = f.readlines()

        def _canon(a: int, b: int) -> Tuple[int, int]:
            return (a, b) if a < b else (b, a)

        # Bonds already present (e.g. written by pdb2gmx from specbond.dat).
        existing_bonds: Set[Tuple[int, int]] = set()
        in_bonds = False
        for line in lines:
            s = line.strip()
            if s.startswith('[ bonds ]'):
                in_bonds = True
                continue
            if in_bonds and s.startswith('['):
                break
            if in_bonds and s and not s.startswith(';'):
                p = s.split()
                if len(p) >= 2:
                    try:
                        existing_bonds.add(_canon(int(p[0]), int(p[1])))
                    except ValueError:
                        pass

        bonds_section_start = -1
        bonds_section_end = -1

        for i, line in enumerate(lines):
            if line.strip().startswith('[ bonds ]'):
                bonds_section_start = i
            elif bonds_section_start >= 0 and line.strip().startswith('[') and not line.strip().startswith('[ bonds ]'):
                bonds_section_end = i
                break

        if bonds_section_start >= 0:
            if bonds_section_end >= 0:
                last_content_line = bonds_section_end - 1
                while (last_content_line > bonds_section_start and
                       (not lines[last_content_line].strip() or
                        lines[last_content_line].strip().startswith(';'))):
                    last_content_line -= 1
                insert_pos = last_content_line + 1
            else:
                insert_pos = len(lines)
        else:
            insert_pos = -1
            for i, line in enumerate(lines):
                if line.strip().startswith('[ atoms ]'):
                    for j in range(i+1, len(lines)):
                        if lines[j].strip().startswith('['):
                            insert_pos = j
                            break
                    break

            if insert_pos >= 0:
                lines.insert(insert_pos, '\n[ bonds ]\n')
                lines.insert(insert_pos + 1, ';   ai    aj funct\n')
                insert_pos += 2

        if insert_pos >= 0:
            crosslink_entries = []
            skipped = 0
            for i, bond_data in enumerate(valid_bond_data):
                atom1_idx, atom2_idx = bond_data['atoms']
                cl1, cl2 = bond_data['cl1'], bond_data['cl2']

                if _canon(atom1_idx, atom2_idx) in existing_bonds:
                    # Already written by pdb2gmx (specbond.dat) -> don't duplicate.
                    skipped += 1
                    continue

                comment = f"; Crosslink bond {i+1}: {cl1.resname}{cl1.resid}{cl1.chain} - {cl2.resname}{cl2.resid}{cl2.chain} (Type: {cl1.type}-{cl2.type})\n"
                bond_entry = f"{atom1_idx} {atom2_idx}     1\n"

                crosslink_entries.append(comment)
                crosslink_entries.append(bond_entry)

            if skipped:
                LOG.debug(
                    f"    Skipped {skipped} crosslink bond(s) already present "
                    f"(pdb2gmx/specbond.dat) in {os.path.basename(itp_file)}"
                )

            for entry in reversed(crosslink_entries):
                lines.insert(insert_pos, entry)

            final_pos = insert_pos + len(crosslink_entries)
            if (final_pos < len(lines) and
                lines[final_pos].strip().startswith('[') and
                (final_pos == 0 or lines[final_pos - 1].strip())):
                lines.insert(final_pos, '\n')

        with open(itp_file, 'w') as f:
            f.writelines(lines)

    def _complete_crosslink_topology(self, itp_file: str, valid_bond_data: List[Dict]) -> None:
        """Complete the bonded topology around the crosslink bonds in one graph
        pass, basically same enumeration pdb2gmx performs from connectivity.

        From the final covalent graph (pdb2gmx's bonds plus the crosslink bonds
        added just before this), emit every term that uses a crosslink bond:
          * angles    : every 2-bond path  i-j-k     -> function 1
          * dihedrals : every 3-bond path  i-j-k-l   -> function 9 (proper)
          * 1-4 pairs : the end atoms of each such 3-bond path, at graph distance
                        exactly 3, excluding H-H       -> function 1
        Terms pdb2gmx already wrote are skipped. No impropers are inferred 
        and no explicit exclusions are written (grompp derives them from nrexcl). 
        Parameters are resolved by grompp from the atom types, as for pdb2gmx-generated terms.
        """
        if not valid_bond_data:
            return

        with open(itp_file, 'r') as f:
            lines = f.readlines()

        def canon2(a: int, b: int) -> Tuple[int, int]:
            return (a, b) if a < b else (b, a)

        def canon3(t: Tuple[int, int, int]) -> Tuple[int, int, int]:
            return t if t[0] <= t[2] else (t[2], t[1], t[0])

        def canon4(t: Tuple[int, int, int, int]) -> Tuple[int, int, int, int]:
            return t if t[0] <= t[3] else (t[3], t[2], t[1], t[0])

        atom_name: Dict[int, str] = {}
        bonds: List[Tuple[int, int]] = []
        have_angles: Set[Tuple[int, int, int]] = set()
        have_diheds: Set[Tuple[int, int, int, int]] = set()
        have_pairs: Set[Tuple[int, int]] = set()
        section = None
        for line in lines:
            s = line.strip()
            if s.startswith('['):
                toks = s.strip('[] ').split()
                section = toks[0].lower() if toks else None
                continue
            if not s or s.startswith((';', '#')):
                continue
            p = s.split()
            try:
                if section == 'atoms' and len(p) >= 5:
                    atom_name[int(p[0])] = p[4]
                elif section == 'bonds' and len(p) >= 2:
                    bonds.append((int(p[0]), int(p[1])))
                elif section == 'angles' and len(p) >= 3:
                    have_angles.add(canon3((int(p[0]), int(p[1]), int(p[2]))))
                elif section == 'dihedrals' and len(p) >= 4:
                    have_diheds.add(canon4((int(p[0]), int(p[1]), int(p[2]), int(p[3]))))
                elif section == 'pairs' and len(p) >= 2:
                    have_pairs.add(canon2(int(p[0]), int(p[1])))
            except ValueError:
                continue

        # Final covalent graph (crosslink bonds already present in [ bonds ]).
        adj: Dict[int, Set[int]] = {}
        for a, b in bonds:
            adj.setdefault(a, set()).add(b)
            adj.setdefault(b, set()).add(a)

        def is_H(i: int) -> bool:
            return atom_name.get(i, '').startswith('H')

        def graph_dist(u: int, v: int, cap: int = 3):
            if u == v:
                return 0
            seen = {u}
            frontier = [u]
            for d in range(1, cap + 1):
                nxt = []
                for a in frontier:
                    for w in adj.get(a, ()):
                        if w == v:
                            return d
                        if w not in seen:
                            seen.add(w)
                            nxt.append(w)
                frontier = nxt
            return None

        xbonds = {canon2(*bd['atoms']) for bd in valid_bond_data}

        angles: Set[Tuple[int, int, int]] = set()
        dihedrals: Set[Tuple[int, int, int, int]] = set()
        for (a1, a2) in xbonds:
            n1 = adj.get(a1, set())
            n2 = adj.get(a2, set())
            # angles using the crosslink bond as an edge
            for x in n1:
                if x != a2:
                    angles.add(canon3((x, a1, a2)))
            for y in n2:
                if y != a1:
                    angles.add(canon3((a1, a2, y)))
            # dihedrals using the crosslink bond as first, central or last bond
            for x in n1:
                if x == a2:
                    continue
                for y in n2:
                    if y != a1 and y != x:
                        dihedrals.add(canon4((x, a1, a2, y)))
            for y in n2:
                if y == a1:
                    continue
                for z in adj.get(y, set()):
                    if z != a1 and z != a2:
                        dihedrals.add(canon4((a1, a2, y, z)))
            for x in n1:
                if x == a2:
                    continue
                for w in adj.get(x, set()):
                    if w != a1 and w != a2:
                        dihedrals.add(canon4((a2, a1, x, w)))

        # 1-4 pairs = the two end atoms of each dihedral, when their shortest
        # graph distance is exactly 3 (a genuine 1-4, not a 1-2/1-3 via a shorter
        # route) and they are not both hydrogens.
        pairs: Set[Tuple[int, int]] = set()
        for (i, j, k, l) in dihedrals:
            if i == l or (is_H(i) and is_H(l)):
                continue
            if graph_dist(i, l) == 3:
                pairs.add(canon2(i, l))

        new_angles = sorted(a for a in angles if a not in have_angles)
        new_diheds = sorted(d for d in dihedrals if d not in have_diheds)
        new_pairs = sorted(pr for pr in pairs if pr not in have_pairs and pr not in xbonds)

        if not (new_angles or new_diheds or new_pairs):
            return

        def append_to_section(header: str, entries: List[str]) -> None:
            if not entries:
                return
            start = end = -1
            for idx, line in enumerate(lines):
                if line.strip().startswith(header):
                    start = idx
                elif start >= 0 and line.strip().startswith('[') and not line.strip().startswith(header):
                    end = idx
                    break
            if start >= 0:
                insert_at = end if end >= 0 else len(lines)
                while insert_at - 1 > start and (
                    not lines[insert_at - 1].strip() or lines[insert_at - 1].strip().startswith(';')
                ):
                    insert_at -= 1
                for entry in reversed(entries):
                    lines.insert(insert_at, entry)
            else:
                lines.append(f"\n{header}\n")
                lines.extend(entries)

        append_to_section('[ angles ]', [f"{a} {b} {c}     1\n" for (a, b, c) in new_angles])
        append_to_section('[ pairs ]', [f"{a} {b}     1\n" for (a, b) in new_pairs])
        append_to_section('[ dihedrals ]', [f"{a} {b} {c} {d}     9\n" for (a, b, c, d) in new_diheds])

        with open(itp_file, 'w') as f:
            f.writelines(lines)

        LOG.debug(
            f"    Completed crosslink topology in {os.path.basename(itp_file)}: "
            f"+{len(new_angles)} angles, +{len(new_diheds)} dihedrals, +{len(new_pairs)} pairs"
        )

    def ensure_posre_include(self, itp_path, group_id):
        """Normalize POSRES include placement."""
        from pathlib import Path
        itp_path = Path(itp_path)
        posre_name = f"posre_{group_id}.itp"
        posre_path = itp_path.parent / posre_name

        generic = itp_path.parent / "posre.itp"
        if not posre_path.exists() and generic.exists():
            try:
                shutil.copy2(generic, posre_path)
            except Exception as e:
                LOG.warning(f"Failed to normalize generic posre.itp -> {posre_name}: {e}")
        if not posre_path.exists():
            LOG.warning(f"Position restraint file not found for {group_id}: expected {posre_name}")
            return

        lines = itp_path.read_text().splitlines()

        new_lines = []
        skip_block = False
        for ln in lines:
            s = ln.strip()
            if s.startswith("#ifdef") and "POSRES" in s:
                if new_lines and new_lines[-1].strip().startswith(";") and (
                    "posre" in new_lines[-1].lower() or "position restraint" in new_lines[-1].lower()
                ):
                    new_lines.pop()
                    if new_lines and not new_lines[-1].strip():
                        new_lines.pop()
                skip_block = True
                continue
            if skip_block:
                if s.startswith("#endif"):
                    skip_block = False
                continue
            if s.startswith("#include") and "posre" in s:
                if new_lines and new_lines[-1].strip().startswith(";") and (
                    "posre" in new_lines[-1].lower() or "position restraint" in new_lines[-1].lower()
                ):
                    new_lines.pop()
                    if new_lines and not new_lines[-1].strip():
                        new_lines.pop()
                continue
            new_lines.append(ln)
        lines = new_lines

        block = [
            "",
            "; Include Position restraint file",
            "#ifdef POSRES",
            f'#include "{posre_name}"',
            "#endif",
            ""
        ]
        lines.extend(block)

        itp_path.write_text("\n".join(lines) + "\n")

    def write_itp(
        self,
        itp_file: Union[str, Path],
        molecule_name: str,
        merged_pdb_file: Optional[str] = None,
        crosslink_pairs: Optional[List[Tuple[Crosslink, Crosslink]]] = None
    ) -> None:
        """Process and write Include Topology (ITP) file with crosslink bonds."""
        itp_file = Path(itp_file)

        with open(itp_file, 'r') as f:
            itp_model = f.readlines()

        try:
            itp_file.unlink()
        except Exception:
            pass

        output_file = itp_file.with_suffix(".itp")

        # pdb2gmx writes a single moleculetype block. GROMACS <2023 names it
        # "Protein_chain_A"; GROMACS >=2023 names it just "Protein". Match either,
        # anchored so comments/other lines don't trigger.
        molname_re = re.compile(r'^\s*Protein(?:_chain_\w+)?\s+\d+\s*$')

        with open(output_file, 'w') as f:
            write = False
            for line in itp_model:
                if 'Include water topology' in line:
                    break
                if write:
                    f.write(line)
                elif molname_re.match(line):
                    f.write('[ moleculetype ]\n')
                    f.write(f'{molecule_name}  3\n')
                    write = True

        if crosslink_pairs and merged_pdb_file:
            LOG.info(f"Adding {len(crosslink_pairs)} crosslink pairs to {output_file.name}")
            self.add_crosslink_topology_to_itp(str(output_file), crosslink_pairs, merged_pdb_file)

    def write_topology(self, topology_file: str, processed_groups: List[Tuple[str, str]]) -> None:
        """Generate AMBER99-ILDNP-STAR force field topology file for connected model groups."""
        if not processed_groups:
            raise ValueError("processed_groups cannot be empty")
        if not self.ff:
            raise ValueError("Force field (self.ff) is not set")

        with open(topology_file, 'w') as f:
            f.write('; Topology for Collagen Microfibril from Colbuilder 2.0\n')
            f.write(f'#include "./{self.ff}/forcefield.itp"\n')

            for group_type, group_id in processed_groups:
                itp_file = f"col_{group_id}.itp"
                if os.path.exists(itp_file):
                    f.write(f'#include "{itp_file}"\n')

            f.write(f'#include "./{self.ff}/ions.itp"\n')
            f.write(f'#include "./{self.ff}/tip3p.itp"\n')
            f.write('\n\n[ system ]\n ;name\nCollagen Microfibril in Water\n\n[ molecules ]\n;name  number\n')

            for group_type, group_id in processed_groups:
                itp_file = f"col_{group_id}.itp"
                if os.path.exists(itp_file):
                    f.write(f'col_{group_id}   1\n')

    def write_gro(self, gro_file: str, processed_groups: List[Tuple[str, str]]) -> None:
        """Write GRO file for connected model groups."""
        if not processed_groups:
            raise ValueError("processed_groups cannot be empty")

        all_atom_lines = []
        last_box_line = "   1.00000   1.00000   1.00000\n"

        for group_type, group_id in processed_groups:
            group_gro = f"col_{group_id}.gro"
            if os.path.exists(group_gro):
                with open(group_gro, 'r') as gro_f:
                    gro_lines = gro_f.readlines()
                    all_atom_lines.extend(gro_lines[2:-1])
                    last_box_line = gro_lines[-1]
                os.remove(group_gro)
            else:
                LOG.warning(f"GRO file not found for group: {group_id}")

        with open(gro_file, 'w') as f:
            f.write("GROMACS GRO-FILE\n")
            f.write(f"{len(all_atom_lines)}\n")
            for line in all_atom_lines:
                f.write(line)
            f.write(last_box_line)

        LOG.info(f"GRO file written with {len(all_atom_lines)} atoms from {len(processed_groups)} groups")


@timeit
async def build_amber99(system: System, config: ColbuilderConfig, file_manager: Optional[FileManager] = None) -> Amber:
    """Build an AMBER99 topology for the given molecular system."""
    ff = f"{config.force_field}sb-star-ildnp"
    ff_name = f"{ff}.ff"
    source_ff_dir = config.FORCE_FIELD_DIR / ff_name
    working_dir = Path.cwd()
    copied_ff_dir = working_dir / ff_name

    amber = Amber(system=system, ff=ff)
    file_manager = file_manager or FileManager(config)
    steps = 3

    LOG.info(f'Step 1/{steps} Setting up Amber99 force field')
    try:
        if not copied_ff_dir.exists():
            if not source_ff_dir.exists():
                raise TopologyGenerationError(
                    message=f"Force field directory not found: {source_ff_dir}",
                    error_code="TOP_ERR_002",
                    context={"ff_dir": str(source_ff_dir)}
                )

            try:
                shutil.copytree(source_ff_dir, copied_ff_dir)
                for filename in REQUIRED_FF_FILES:
                    src_file = source_ff_dir / filename
                    if not src_file.exists():
                        raise TopologyGenerationError(
                            message=f"Required force field file not found: {filename}",
                            error_code="TOP_ERR_003",
                            context={"missing_file": filename}
                        )
                    dest_file = working_dir / filename
                    shutil.copy2(src_file, dest_file)

            except Exception as e:
                raise TopologyGenerationError(
                    message="Force field setup failed",
                    original_error=e,
                    error_code="TOP_ERR_002",
                    context={"ff_dir": str(source_ff_dir)}
                )

        LOG.info(f'Step 2/{steps} Grouping connected models and processing with GROMACS')

        connected_groups = amber.get_connected_groups()
        LOG.debug(f"    Found {len(connected_groups)} molecular groups: {connected_groups}")

        processed_groups = []

        for group in connected_groups:
            try:
                merge_result = amber.merge_connected_models(group)
                if merge_result is None:
                    LOG.warning(f"Skipping group {group} - merge failed")
                    continue

                model_type, group_id, crosslink_pairs = merge_result
                merge_pdb_path = working_dir / model_type / f"{group_id}.merge.pdb"

                if not merge_pdb_path.exists() or not os.path.getsize(merge_pdb_path):
                    LOG.error(f'Invalid merged PDB file: {merge_pdb_path}')
                    continue

                gmx_cmd = (f'export GMXLIB={working_dir} && gmx pdb2gmx -f {merge_pdb_path} '
                          f'-ignh -merge all -ff {ff} -water tip3p '
                          f'-p col_{group_id}.top -o col_{group_id}.gro '
                          f'-i posre_{group_id}.itp')

                result = await asyncio.create_subprocess_shell(
                    gmx_cmd,
                    stdout=asyncio.subprocess.PIPE,
                    stderr=asyncio.subprocess.PIPE
                )
                stdout, stderr = await result.communicate()

                if result.returncode != 0:
                    LOG.error(f"GROMACS stderr for group {group}: {stderr.decode()}")
                    raise TopologyGenerationError(
                        message=f'GROMACS pdb2gmx failed for group {group}',
                        error_code="TOP_ERR_005",
                        context={
                            "group": str(group),
                            "stderr": stderr.decode()
                        }
                    )

                molecule_name = f"col_{group_id}"
                amber.write_itp(
                    itp_file=working_dir / f'col_{group_id}.top',
                    molecule_name=f'col_{group_id}',
                    merged_pdb_file=str(merge_pdb_path),
                    crosslink_pairs=crosslink_pairs
                )
                amber.ensure_posre_include(working_dir / f'col_{group_id}.itp', group_id)

                processed_groups.append((model_type, group_id))
                LOG.debug(f"    Successfully processed connected group {group} as {molecule_name}")

            except TopologyGenerationError:
                raise
            except Exception as e:
                LOG.error(f'Group {group} processing failed: {str(e)}')

        if not processed_groups:
            raise TopologyGenerationError(
                message='No model groups were successfully processed',
                error_code="TOP_ERR_006"
            )

        LOG.info(f'Step 3/{steps} Generating system topology files')
        try:
            topology_file = str(working_dir / f"collagen_fibril_{config.species}.top")
            gro_file = str(working_dir / f"collagen_fibril_{config.species}.gro")

            amber.write_topology(topology_file=topology_file, processed_groups=processed_groups)
            amber.write_gro(gro_file=gro_file, processed_groups=processed_groups)

            LOG.info(f"Successfully generated topology for {len(processed_groups)} molecular groups")
            LOG.debug(f"    Groups processed: {[group_id for _, group_id in processed_groups]}")

        except Exception as e:
            raise TopologyGenerationError(
                message='Final topology file generation failed',
                original_error=e,
                error_code="TOP_ERR_007",
                context={"output": config.species}
            )

        return amber

    except TopologyGenerationError:
        raise
    except Exception as e:
        raise TopologyGenerationError(
            message="Amber topology generation failed",
            original_error=e,
            error_code="TOP_ERR_001",
            context={"force_field": ff}
        )
