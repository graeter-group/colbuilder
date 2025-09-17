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
    
    def merge_connected_models(self, model_group: List[int]) -> Optional[Tuple[str, str]]:
        if not model_group or not self.system:
            return None

        first_model = self.system.get_model(model_id=float(model_group[0]))
        if not first_model or not first_model.type:
            return None

        model_type = first_model.type
        os.makedirs(model_type, exist_ok=True)

        group_id = "_".join(str(int(mid)) for mid in sorted(model_group))
        output_file = os.path.join(model_type, f"{group_id}.merge.pdb")

        def write_caps(path, out):
            if os.path.exists(path):
                with open(path, "r") as f_in:
                    for line in f_in:
                        if line.startswith(self.pdb_line_types):
                            out.write(line)
            else:
                LOG.warning(f"Caps file not found: {path}")

        with open(output_file, "w") as f_out:
            # Always include each model’s own caps
            for mid in model_group:
                caps = os.path.join(model_type, f"{int(mid)}.caps.pdb")
                write_caps(caps, f_out)

            # Optionally also include connection partners’ caps (if different from self)
            for mid in model_group:
                model = self.system.get_model(model_id=float(mid))
                if model and model.connect:
                    for cid in model.connect:
                        if int(cid) not in model_group:
                            caps = os.path.join(model_type, f"{int(cid)}.caps.pdb")
                            write_caps(caps, f_out)

            f_out.write("END\n")

        if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
            return (model_type, group_id)
        LOG.error(f"Failed to create merged PDB for group {model_group}")
        return None


    def find_crosslink_pairs(self, merged_pdb_file: str, distance_cutoff: float = 5.0) -> List[Tuple[Crosslink, Crosslink]]:
        """Find crosslink pairs that should be bonded based on distance and type."""
        crosslinks = read_crosslink(merged_pdb_file)
        pairs = []
        
        for i, cl1 in enumerate(crosslinks):
            for j, cl2 in enumerate(crosslinks[i+1:], i+1):
                if (cl1.resname == cl2.resname and cl1.resid == cl2.resid and cl1.chain == cl2.chain):
                    continue
                
                if cl1.type == cl2.type:
                    distance = np.linalg.norm(cl1.position - cl2.position)
                    if distance <= distance_cutoff:
                        pairs.append((cl1, cl2))
                        LOG.debug(f"        Found crosslink pair: {cl1.resname}{cl1.resid}{cl1.chain} - {cl2.resname}{cl2.resid}{cl2.chain} (distance: {distance:.2f} Å, type: {cl1.type}-{cl2.type})")

        return pairs

    def find_atom_indices_for_crosslinks(self, itp_file: str, crosslink_pairs: List[Tuple[Crosslink, Crosslink]], merged_pdb_file: str) -> List[Tuple[int, int]]:
        """Find atom indices in the topology file for crosslink pairs."""
        atom_indices = []
        
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
                        'residue_name': parts[3]
                    })
        
        gro_file = itp_file.replace('.itp', '.gro')
        gro_coords = {}
        
        if os.path.exists(gro_file):
            with open(gro_file, 'r') as f:
                gro_lines = f.readlines()
                for line in gro_lines[2:-1]:
                    if len(line) >= 44:
                        try:
                            atom_idx = int(line[15:20].strip())
                            x = float(line[20:28]) * 10
                            y = float(line[28:36]) * 10
                            z = float(line[36:44]) * 10
                            gro_coords[atom_idx] = np.array([x, y, z])
                        except (ValueError, IndexError):
                            continue
        
        for cl1, cl2 in crosslink_pairs:
            best_atom1_idx = None
            best_atom2_idx = None
            best_dist1 = float('inf')
            best_dist2 = float('inf')
            
            for atom in atom_data:
                atom_idx = atom['index']
                
                if (atom['residue_name'] == cl1.resname and 
                    self._is_crosslink_atom(atom['atom_name'], cl1.resname, cl1.type)):
                    
                    if atom_idx in gro_coords:
                        dist = np.linalg.norm(gro_coords[atom_idx] - cl1.position)
                        if dist < best_dist1 and dist < 5.0:
                            best_dist1 = dist
                            best_atom1_idx = atom_idx
                
                if (atom['residue_name'] == cl2.resname and 
                    self._is_crosslink_atom(atom['atom_name'], cl2.resname, cl2.type)):
                    
                    if atom_idx in gro_coords:
                        dist = np.linalg.norm(gro_coords[atom_idx] - cl2.position)
                        if dist < best_dist2 and dist < 5.0:
                            best_dist2 = dist
                            best_atom2_idx = atom_idx
            
            if best_atom1_idx and best_atom2_idx:
                atom_indices.append((best_atom1_idx, best_atom2_idx))
            else:
                LOG.warning(f"Could not find atoms for crosslink pair: {cl1.resname}{cl1.resid} - {cl2.resname}{cl2.resid}")
        
        return atom_indices

    def _is_crosslink_atom(self, atom_name: str, resname: str, crosslink_type: str) -> bool:
        """Check if this atom is the target bonding atom for the crosslink type."""
        if crosslink_type == "T":
            if resname in ("LYX", "LXY", "LYY", "LXX") and atom_name in ("C13", "C12"):
                return True
            elif resname in ("LY3", "LX3", "L3Y", "L2Y", "L3X", "L2X") and atom_name == "CG":
                return True
            elif resname in ("LY2", "LX2") and atom_name == "CB":
                return True
        elif crosslink_type == "D":
            if resname in ("L4Y", "L4X", "LY4", "LX4") and atom_name == "CE":
                return True
            elif resname in ("L5Y", "L5X", "LY5", "LX5") and atom_name == "NZ":
                return True
            elif resname in ("LGX", "LPS") and atom_name == "CE":
                return True
            elif resname in ("AGS", "APD") and atom_name == "NZ":
                return True
        
        return False

    def parse_topology_sections(self, itp_file: str) -> Dict[str, List[List[int]]]:
        """Parse existing topology to get bonds, angles, and dihedrals."""
        topology = {'bonds': [], 'angles': [], 'dihedrals': []}
        
        try:
            with open(itp_file, 'r') as f:
                lines = f.readlines()
        except Exception as e:
            LOG.warning(f"Could not read {itp_file}: {e}")
            return topology
        
        current_section = None
        
        for line in lines:
            line = line.strip()
            if not line or line.startswith(';'):
                continue
                
            if line.startswith('[ bonds ]'):
                current_section = 'bonds'
            elif line.startswith('[ angles ]'):
                current_section = 'angles'
            elif line.startswith('[ dihedrals ]'):
                current_section = 'dihedrals'
            elif line.startswith('['):
                current_section = None
            elif current_section and not line.startswith(';'):
                parts = line.split()
                try:
                    if len(parts) >= 3 and current_section == 'bonds':
                        topology['bonds'].append([int(parts[0]), int(parts[1])])
                    elif len(parts) >= 4 and current_section == 'angles':
                        topology['angles'].append([int(parts[0]), int(parts[1]), int(parts[2])])
                    elif len(parts) >= 5 and current_section == 'dihedrals':
                        topology['dihedrals'].append([int(parts[0]), int(parts[1]), int(parts[2]), int(parts[3])])
                except (ValueError, IndexError):
                    continue
        
        return topology

    def build_connectivity_graph(self, bonds: List[List[int]]) -> Dict[int, Set[int]]:
        """Build a connectivity graph from bond list."""
        graph = {}
        for bond in bonds:
            atom1, atom2 = bond[0], bond[1]
            if atom1 not in graph:
                graph[atom1] = set()
            if atom2 not in graph:
                graph[atom2] = set()
            graph[atom1].add(atom2)
            graph[atom2].add(atom1)
        
        return graph

    def generate_crosslink_angles(self, crosslink_bonds: List[Tuple[int, int]], connectivity: Dict[int, Set[int]]) -> List[Tuple[int, int, int]]:
        """Generate angles involving crosslink bonds."""
        angles = []
        
        for atom1, atom2 in crosslink_bonds:
            if atom1 in connectivity:
                for x in connectivity[atom1]:
                    if x != atom2:
                        angles.append((x, atom1, atom2))
            
            if atom2 in connectivity:
                for y in connectivity[atom2]:
                    if y != atom1:
                        angles.append((atom1, atom2, y))
        
        return angles

    def generate_crosslink_dihedrals(self, crosslink_bonds: List[Tuple[int, int]], connectivity: Dict[int, Set[int]]) -> List[Tuple[int, int, int, int]]:
        """Generate dihedrals involving crosslink bonds."""
        dihedrals = []
        
        for atom1, atom2 in crosslink_bonds:
            if atom1 in connectivity and atom2 in connectivity:
                for x in connectivity[atom1]:
                    if x != atom2:
                        for y in connectivity[atom2]:
                            if y != atom1 and y != x:
                                dihedrals.append((x, atom1, atom2, y))
        
        return dihedrals

    def _is_backbone_atom(self, atom_name: str) -> bool:
        """Check if an atom is a backbone atom."""
        return atom_name in ('N', 'CA', 'C', 'O', 'H')

    def _get_atom_name_by_index(self, itp_file: str, atom_index: int) -> Optional[str]:
        """Get atom name by its index from the topology file."""
        try:
            with open(itp_file, 'r') as f:
                lines = f.readlines()
            
            atoms_section = False
            for line in lines:
                if line.strip().startswith('[ atoms ]'):
                    atoms_section = True
                    continue
                elif atoms_section and line.strip().startswith('['):
                    break
                elif atoms_section and not line.strip().startswith(';') and line.strip():
                    parts = line.split()
                    if len(parts) >= 5 and int(parts[0]) == atom_index:
                        return parts[4]
        except Exception:
            pass
        return None

    def _dihedral_involves_backbone(self, itp_file: str, dihedral: Tuple[int, int, int, int]) -> bool:
        """Check if a dihedral involves any backbone atoms."""
        for atom_idx in dihedral:
            atom_name = self._get_atom_name_by_index(itp_file, atom_idx)
            if atom_name and self._is_backbone_atom(atom_name):
                return True
        return False

    def add_crosslink_topology_to_itp(self, itp_file: str, crosslink_pairs: List[Tuple[Crosslink, Crosslink]], merged_pdb_file: str) -> None:
        """Add complete crosslink topology (bonds, angles, dihedrals) to an existing ITP file."""
        if not crosslink_pairs:
            return

        LOG.debug(f"        Adding crosslink topology for {len(crosslink_pairs)} crosslink pairs to {itp_file}")

        try:
            atom_indices = self.find_atom_indices_for_crosslinks(itp_file, crosslink_pairs, merged_pdb_file)
            
            if not atom_indices:
                LOG.warning(f"No atom indices found for crosslink topology in {itp_file}")
                return
            
            valid_bonds = set()
            valid_bond_data = []
            
            for (cl1, cl2), (atom1_idx, atom2_idx) in zip(crosslink_pairs, atom_indices):
                if atom1_idx == atom2_idx:
                    continue
                
                bond_key = tuple(sorted([atom1_idx, atom2_idx]))
                
                if bond_key not in valid_bonds:
                    valid_bonds.add(bond_key)
                    valid_bond_data.append({
                        'atoms': bond_key,
                        'cl1': cl1,
                        'cl2': cl2,
                        'original_indices': (atom1_idx, atom2_idx)
                    })
            
            if not valid_bond_data:
                LOG.warning(f"No valid crosslink topology to add to {itp_file}")
                return
            
            self._add_crosslink_bonds(itp_file, valid_bond_data)
            
            try:
                self._add_crosslink_angles_and_dihedrals(itp_file, valid_bond_data)
            except Exception as e:
                LOG.warning(f"Failed to add angles/dihedrals, but bonds were added successfully: {str(e)}")
            
            LOG.debug(f"    Successfully added crosslink topology to {itp_file}")
            
        except Exception as e:
            LOG.error(f"Failed to add crosslink topology to {itp_file}: {str(e)}")

    def _add_crosslink_bonds(self, itp_file: str, valid_bond_data: List[Dict]) -> None:
        """Add crosslink bonds using standard GROMACS format."""
        with open(itp_file, 'r') as f:
            lines = f.readlines()
        
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
            for i, bond_data in enumerate(valid_bond_data):
                atom1_idx, atom2_idx = bond_data['atoms']
                cl1, cl2 = bond_data['cl1'], bond_data['cl2']
                
                comment = f"; Crosslink bond {i+1}: {cl1.resname}{cl1.resid}{cl1.chain} - {cl2.resname}{cl2.resid}{cl2.chain} (Type: {cl1.type}-{cl2.type})\n"
                bond_entry = f"{atom1_idx} {atom2_idx}     1\n"
                
                crosslink_entries.append(comment)
                crosslink_entries.append(bond_entry)
            
            for entry in reversed(crosslink_entries):
                lines.insert(insert_pos, entry)
            
            final_pos = insert_pos + len(crosslink_entries)
            if (final_pos < len(lines) and 
                lines[final_pos].strip().startswith('[') and 
                (final_pos == 0 or lines[final_pos - 1].strip())):
                lines.insert(final_pos, '\n')
        
        with open(itp_file, 'w') as f:
            f.writelines(lines)

    def _add_crosslink_angles_and_dihedrals(self, itp_file: str, valid_bond_data: List[Dict]) -> None:
        """Add angles and dihedrals for crosslink bonds."""
        try:
            existing_topology = self.parse_topology_sections(itp_file)
            
            all_bonds = existing_topology['bonds'].copy()
            crosslink_bonds = [bond_data['atoms'] for bond_data in valid_bond_data]
            all_bonds.extend(crosslink_bonds)
            
            connectivity = self.build_connectivity_graph(all_bonds)
            
            crosslink_angles = self.generate_crosslink_angles(crosslink_bonds, connectivity)
            crosslink_dihedrals = self.generate_crosslink_dihedrals(crosslink_bonds, connectivity)

            if crosslink_angles or crosslink_dihedrals:
                with open(itp_file, 'r') as f:
                    lines = f.readlines()
                
                if crosslink_angles:
                    lines = self._add_angles_to_lines(lines, crosslink_angles)
                
                if crosslink_dihedrals:
                    lines = self._add_dihedrals_to_lines(lines, crosslink_dihedrals)
                
                with open(itp_file, 'w') as f:
                    f.writelines(lines)
            
        except Exception as e:
            LOG.warning(f"Could not add angles/dihedrals: {str(e)}")

    def _add_angles_to_lines(self, lines: List[str], crosslink_angles: List[Tuple[int, int, int]]) -> List[str]:
        """Add crosslink angles using standard GROMACS format."""
        if not crosslink_angles:
            return lines
            
        angles_section_start = -1
        angles_section_end = -1
        
        for i, line in enumerate(lines):
            if line.strip().startswith('[ angles ]'):
                angles_section_start = i
            elif angles_section_start >= 0 and line.strip().startswith('[') and not line.strip().startswith('[ angles ]'):
                angles_section_end = i
                break
        
        if angles_section_start >= 0:
            if angles_section_end >= 0:
                last_content_line = angles_section_end - 1
                while (last_content_line > angles_section_start and 
                       (not lines[last_content_line].strip() or 
                        lines[last_content_line].strip().startswith(';'))):
                    last_content_line -= 1
                insert_pos = last_content_line + 1
            else:
                insert_pos = len(lines)
        else:
            insert_pos = self._find_section_end(lines, '[ bonds ]')
            if insert_pos >= 0:
                lines.insert(insert_pos, '\n[ angles ]\n')
                lines.insert(insert_pos + 1, ';   ai    aj    ak funct\n')
                insert_pos += 2
        
        if insert_pos >= 0:
            angle_entries = []
            angle_entries.append("; Crosslink angles\n")
            for atom1, atom2, atom3 in crosslink_angles:
                angle_entry = f"{atom1} {atom2} {atom3}     1\n"
                angle_entries.append(angle_entry)
            
            for entry in reversed(angle_entries):
                lines.insert(insert_pos, entry)
            
            final_pos = insert_pos + len(angle_entries)
            if (final_pos < len(lines) and 
                lines[final_pos].strip().startswith('[') and 
                (final_pos == 0 or lines[final_pos - 1].strip())):
                lines.insert(final_pos, '\n')
        
        return lines

    def _add_dihedrals_to_lines(self, lines: List[str], crosslink_dihedrals: List[Tuple[int, int, int, int]]) -> List[str]:
        """Add crosslink dihedrals using standard GROMACS format with proper type assignment."""
        if not crosslink_dihedrals:
            return lines
        
        temp_itp_file = "temp_for_backbone_check.itp"
        with open(temp_itp_file, 'w') as f:
            f.writelines(lines)
            
        dihedrals_section_start = -1
        dihedrals_section_end = -1
        
        for i, line in enumerate(lines):
            if line.strip().startswith('[ dihedrals ]'):
                dihedrals_section_start = i
            elif dihedrals_section_start >= 0 and line.strip().startswith('[') and not line.strip().startswith('[ dihedrals ]'):
                dihedrals_section_end = i
                break
        
        if dihedrals_section_start >= 0:
            if dihedrals_section_end >= 0:
                last_content_line = dihedrals_section_end - 1
                while (last_content_line > dihedrals_section_start and 
                       (not lines[last_content_line].strip() or 
                        lines[last_content_line].strip().startswith(';'))):
                    last_content_line -= 1
                insert_pos = last_content_line + 1
            else:
                insert_pos = len(lines)
        else:
            insert_pos = self._find_section_end(lines, '[ angles ]')
            if insert_pos >= 0:
                lines.insert(insert_pos, '\n[ dihedrals ]\n')
                lines.insert(insert_pos + 1, ';   ai    aj    ak    al funct\n')
                insert_pos += 2
        
        if insert_pos >= 0:
            dihedral_entries = []
            dihedral_entries.append("; Crosslink dihedrals\n")
            
            for atom1, atom2, atom3, atom4 in crosslink_dihedrals:
                if self._dihedral_involves_backbone(temp_itp_file, (atom1, atom2, atom3, atom4)):
                    dihedral_type = 4
                else:
                    dihedral_type = 9
                
                dihedral_entry = f"{atom1} {atom2} {atom3} {atom4}     {dihedral_type}\n"
                dihedral_entries.append(dihedral_entry)
            
            for entry in reversed(dihedral_entries):
                lines.insert(insert_pos, entry)
            
            final_pos = insert_pos + len(dihedral_entries)
            if (final_pos < len(lines) and 
                lines[final_pos].strip().startswith('[') and 
                (final_pos == 0 or lines[final_pos - 1].strip())):
                lines.insert(final_pos, '\n')
        
        try:
            os.remove(temp_itp_file)
        except:
            pass
        
        return lines

    def _find_section_end(self, lines: List[str], section_header: str) -> int:
        """Find the end of a given section."""
        for i, line in enumerate(lines):
            if line.strip().startswith(section_header):
                for j in range(i+1, len(lines)):
                    if lines[j].strip().startswith('['):
                        return j
        return len(lines)

    def ensure_posre_include(self, itp_path, group_id):
        """
        Normalize POSRES include placement:
        - Remove existing POSRES blocks or lone #include lines (and their leading comment).
        - Append a clean POSRES block at the END of the .itp (valid: after [atoms] and any other sections).
        """
        from pathlib import Path
        itp_path = Path(itp_path)
        posre_name = f"posre_{group_id}.itp"
        posre_path = itp_path.parent / posre_name

        # Ensure a posre file exists; normalize from generic if needed
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

        # 1) Strip existing POSRES includes/blocks + any leading comment line
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

        # 2) Append clean POSRES block at EOF
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


    def write_itp(self, itp_file: Union[str, Path], molecule_name: str, merged_pdb_file: Optional[str] = None) -> None:
        """Process and write Include Topology (ITP) file with crosslink bonds."""
        itp_file = Path(itp_file)
        
        with open(itp_file, 'r') as f:
            itp_model = f.readlines()

        try:
            itp_file.unlink()
        except Exception:
            pass
        
        output_file = itp_file.parent / str(itp_file.name).replace("top", "itp")
        
        with open(output_file, 'w') as f:
            write = False
            for line in itp_model:
                if 'Include water topology' in line:
                    break
                if write:
                    f.write(line)
                elif 'Protein_chain_A' in line:
                    f.write('[ moleculetype ]\n')
                    f.write(f'{molecule_name}  3\n')
                    write = True
        
        if merged_pdb_file and os.path.exists(merged_pdb_file):
            try:
                crosslink_pairs = self.find_crosslink_pairs(merged_pdb_file)
                if crosslink_pairs:
                    self.add_crosslink_topology_to_itp(str(output_file), crosslink_pairs, merged_pdb_file)
            except Exception as e:
                LOG.warning(f"Failed to add crosslink topology: {str(e)}")

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
                
                model_type, group_id = merge_result
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
                    merged_pdb_file=str(merge_pdb_path)
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