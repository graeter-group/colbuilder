# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0

import numpy as np
from sklearn.metrics import pairwise_distances as pdist
from typing import List, Dict, Any, Optional, Tuple, Union
import os

from colbuilder.core.utils.logger import setup_logger

LOG = setup_logger(__name__)


class Crosslink:
    """
    Setup crosslink topology for collagen models.
    
    This class handles the identification and parameterization of crosslinks 
    in collagen molecular models. It processes merged PDB files to identify 
    crosslink sites and generates the necessary bonded parameters (bonds, 
    angles, dihedrals) for crosslinks in Martini coarse-grained models.
    
    The class supports two types of crosslinks:
    - Divalent HLKNL-crosslinks between L4Y and L5Y residues
    - Trivalent PYD-crosslinks between LYX, LY2, and LY3 residues
    """
    
    def __init__(self, cnt_model: Optional[int] = None) -> None:
        """
        Initialize the Crosslink object with model parameters.
        
        Parameters
        ----------
        cnt_model : Optional[int]
            Model counter used for file naming. If provided, the merged PDB
            file will be named '{cnt_model}.merge.pdb'
        """
        self.file: str = f"{int(cnt_model)}.merge.pdb" if cnt_model is not None else ""
        self.crosslink_coords: List[List[float]] = []
        self.crosslink_pdb: List[List[Any]] = []  # Mix of str and float
        self.crosslink_neighbors: List[Any] = []
        self.crosslink_connect: List[List[List[Any]]] = []
        self.crosslink_pairs: List[Tuple[List[Any], List[Any]]] = []  # Store valid pairs
        self.crosslink_bonded: Dict[str, List[List[Any]]] = {
            'bonds': [], 
            'angles': [], 
            'dihedrals': []
        }
        
        # Distance thresholds for different crosslink types (in Å)
        # Adjust these values based on your coarse-grained model
        self.crosslink_thresholds = {
            'LYX_LY2': 15.0,  # Maximum distance for LYX SC4 - LY2 SC1 crosslinks
            'LYX_LY3': 15.0,  # Maximum distance for LYX SC5 - LY3 SC1 crosslinks
            'L4Y_L5Y': 15.0   # Maximum distance for L4Y SC1 - L5Y SC2 crosslinks
        }
        
        # HLKNL-crosslink parameters
        self.dly45: str = '0.415'    # L4Y-L5Y bond equilibrium distance (nm)
        self.kly45: str = '7000'     # L4Y-L5Y bond force constant (kJ/mol/nm^2)
        self.al45y_1: str = '140'    # L4Y-L5Y SC1-SC2-SC1(L4Y) angle (degrees)
        self.al45y_2: str = '140'    # L4Y-L5Y SC2(L5Y)-SC1-BB angle (degrees)
        self.k_angle: str = '153'    # Universal angle force constant (kJ/mol/rad^2)
        
        # PYD-crosslink parameters
        self.klyxly2: str = '9000'   # LYX-LY2 bond force constant (kJ/mol/nm^2)
        self.klyxly3: str = '12000'  # LYX-LY3 bond force constant (kJ/mol/nm^2)
        self.dlyxly2: str = '0.290'  # LYX-LY2 bond equilibrium distance (nm)
        self.dlyxly3: str = '0.230'  # LYX-LY3 bond equilibrium distance (nm)
        
        # PYD-crosslink angle parameters (degrees)
        self.al2yx_1: str = '100'    # LY2-LYX TP1q-TC6q-TC4 angle
        self.al2yx_2: str = '60'     # LY2-LYX TQ2p-TP1q-TC4 angle
        self.al2yx_3: str = '130'    # LY2-LYX TP1q-TC4-SP2 angle
        self.al3yx_1: str = '100'    # LY3-LYX TP1q-TC6q-TC4 angle
        self.al3yx_2: str = '110'    # LY3-LYX TQ2p-TC6q-TC4 angle
        self.al3yx_3: str = '130'    # LY3-LYX TC6q-TC4-SP2 angle
        
        LOG.debug(f"Initialized Crosslink for model counter {cnt_model}")

    def get_crosslink_coords(self, cnt_model: Optional[int] = None) -> List[List[float]]:
        """
        Extract coordinates of crosslink sites from a PDB file.
        
        Following the working version pattern: stores coordinates as floats.
        """
        if cnt_model is None:
            file = self.file
        else:
            file = f"{int(cnt_model)}.merge.pdb"
            
        LOG.debug(f"Reading crosslink coordinates from {file}")
        
        if not os.path.exists(file):
            LOG.error(f"PDB file not found: {file}")
            return self.crosslink_coords
            
        self.crosslink_coords = []
        self.crosslink_pdb = []
        it_pdb = 0
        
        try:
            with open(file, 'r') as f:
                for line in f:
                    if line[0:4] == 'ATOM':
                        it_pdb += 1
                        
                    if ((line[17:20] == 'LYX' and line[12:15] in ['SC4', 'SC5']) or
                        (line[17:20] in ['LY2', 'LY3'] and line[12:15] == 'SC1') or
                        (line[17:20] == 'L4Y' and line[12:15] == 'SC1') or
                        (line[17:20] == 'L5Y' and line[12:15] == 'SC2')):
                        
                        coords = [float(line[29:38]), float(line[38:46]), float(line[46:56])]
                        
                        self.crosslink_pdb.append([
                            str(it_pdb),           # atom index as string
                            line[17:20],           # residue name
                            line[12:15],           # atom name  
                            line[21:26],           # chain info
                            coords[0],             # x coordinate as float
                            coords[1],             # y coordinate as float
                            coords[2]              # z coordinate as float
                        ])
                        self.crosslink_coords.append(coords)
                        
                        LOG.debug(f"Found crosslink site: {line[17:20].strip()} {line[12:15].strip()} (atom {it_pdb})")
            
            LOG.debug(f"Found {len(self.crosslink_coords)} crosslink sites")
                
        except Exception as e:
            LOG.error(f"Error reading PDB file {file}: {str(e)}")
            
        return self.crosslink_coords

    def get_crosslink_connect(self, cnt_model: Optional[int] = None) -> List[List[List[Any]]]:
        """
        Get nearest crosslinks to determine connections.
        
        This method is kept for backwards compatibility but now uses the improved
        pair-finding algorithm internally.
        """
        LOG.debug("Using improved crosslink pair detection algorithm")
        
        # Use the new pair-finding method and convert to old format for compatibility
        pairs = self.find_crosslink_pairs(cnt_model=cnt_model)
        
        self.crosslink_connect = []
        if pairs:
            # Group all atoms involved in pairs
            all_atoms = []
            for pair in pairs:
                if pair[0] not in all_atoms:
                    all_atoms.append(pair[0])
                if pair[1] not in all_atoms:
                    all_atoms.append(pair[1])
            
            if all_atoms:
                self.crosslink_connect.append(all_atoms)
                LOG.debug(f"Converted {len(pairs)} pairs to connection group with {len(all_atoms)} atoms")
        
        return self.crosslink_connect

    def find_crosslink_pairs(self, cnt_model: Optional[int] = None) -> List[Tuple[List[Any], List[Any]]]:
        """
        Find valid crosslink pairs based on distance and compatibility.
        
        This replaces the old distance-matrix approach with a more targeted method
        that specifically looks for compatible crosslink pairs.
        """
        self.get_crosslink_coords(cnt_model=cnt_model)
        
        if not self.crosslink_coords:
            LOG.warning("No crosslink coordinates found")
            return []
        
        self.crosslink_pairs = []
        
        # Separate atoms by type for targeted pairing
        lyx_sc4_atoms = []
        lyx_sc5_atoms = []
        ly2_sc1_atoms = []
        ly3_sc1_atoms = []
        l4y_sc1_atoms = []
        l5y_sc2_atoms = []
        
        for atom in self.crosslink_pdb:
            if atom[1] == 'LYX' and atom[2] == 'SC4':
                lyx_sc4_atoms.append(atom)
            elif atom[1] == 'LYX' and atom[2] == 'SC5':
                lyx_sc5_atoms.append(atom)
            elif atom[1] == 'LY2' and atom[2] == 'SC1':
                ly2_sc1_atoms.append(atom)
            elif atom[1] == 'LY3' and atom[2] == 'SC1':
                ly3_sc1_atoms.append(atom)
            elif atom[1] == 'L4Y' and atom[2] == 'SC1':
                l4y_sc1_atoms.append(atom)
            elif atom[1] == 'L5Y' and atom[2] == 'SC2':
                l5y_sc2_atoms.append(atom)
        
        LOG.debug(f"Crosslink atoms:")
        LOG.debug(f"  LYX SC4: {len(lyx_sc4_atoms)}")
        LOG.debug(f"  LYX SC5: {len(lyx_sc5_atoms)}")
        LOG.debug(f"  LY2 SC1: {len(ly2_sc1_atoms)}")
        LOG.debug(f"  LY3 SC1: {len(ly3_sc1_atoms)}")
        LOG.debug(f"  L4Y SC1: {len(l4y_sc1_atoms)}")
        LOG.debug(f"  L5Y SC2: {len(l5y_sc2_atoms)}")
        
        # Find LYX SC4 - LY2 SC1 pairs
        for lyx_atom in lyx_sc4_atoms:
            for ly2_atom in ly2_sc1_atoms:
                dist = np.linalg.norm(np.array(lyx_atom[-3:]) - np.array(ly2_atom[-3:]))
                LOG.debug(f"    Distance LYX SC4 (atom {lyx_atom[0]}) - LY2 SC1 (atom {ly2_atom[0]}): {dist:.3f} Å")
                
                if dist <= self.crosslink_thresholds['LYX_LY2']:
                    self.crosslink_pairs.append((lyx_atom, ly2_atom))
                    LOG.info(f" Added LYX-LY2 pair: atoms {lyx_atom[0]} - {ly2_atom[0]} (distance: {dist:.3f} Å)")
        
        # Find LYX SC5 - LY3 SC1 pairs
        for lyx_atom in lyx_sc5_atoms:
            for ly3_atom in ly3_sc1_atoms:
                dist = np.linalg.norm(np.array(lyx_atom[-3:]) - np.array(ly3_atom[-3:]))
                LOG.debug(f"    Distance LYX SC5 (atom {lyx_atom[0]}) - LY3 SC1 (atom {ly3_atom[0]}): {dist:.3f} Å")
                
                if dist <= self.crosslink_thresholds['LYX_LY3']:
                    self.crosslink_pairs.append((lyx_atom, ly3_atom))
                    LOG.info(f" Added LYX-LY3 pair: atoms {lyx_atom[0]} - {ly3_atom[0]} (distance: {dist:.3f} Å)")
        
        # Find L4Y SC1 - L5Y SC2 pairs
        for l4y_atom in l4y_sc1_atoms:
            for l5y_atom in l5y_sc2_atoms:
                dist = np.linalg.norm(np.array(l4y_atom[-3:]) - np.array(l5y_atom[-3:]))
                LOG.debug(f"    Distance L4Y SC1 (atom {l4y_atom[0]}) - L5Y SC2 (atom {l5y_atom[0]}): {dist:.3f} Å")
                
                if dist <= self.crosslink_thresholds['L4Y_L5Y']:
                    self.crosslink_pairs.append((l4y_atom, l5y_atom))
                    LOG.info(f" Added L4Y-L5Y pair: atoms {l4y_atom[0]} - {l5y_atom[0]} (distance: {dist:.3f} Å)")
        
        return self.crosslink_pairs
    
    def set_crosslink_bonded(self, cnt_model: Optional[int] = None, 
                           crosslink_connect: Optional[List[List[List[Any]]]] = None) -> Dict[str, List[List[Any]]]:
        """
        Setup topology for crosslink bonded parameters: bonds, angles and dihedrals.
        
        Now uses the improved pair-finding algorithm for better crosslink detection.
        """
        LOG.debug(f"Setting up crosslink bonded parameters for model {cnt_model}")
        
        self.crosslink_bonded = {'bonds': [], 'angles': [], 'dihedrals': []}
        
        self.find_crosslink_pairs(cnt_model=cnt_model)
        
        if not self.crosslink_pairs:
            LOG.debug("No crosslink pairs found, returning empty parameters")
            return self.crosslink_bonded
            
        connections_found = 0
        
        try:
            for clx, cly in self.crosslink_pairs:
                dist = np.linalg.norm(np.array(clx[-3:]) - np.array(cly[-3:]))
                
                LOG.debug(f"Processing crosslink between {clx[1]}{clx[2]} and {cly[1]}{cly[2]}: {dist:.3f} Å")
                
                # LYX SC4 - LY2 SC1 crosslinks
                if (clx[1] == 'LYX' and clx[2] == 'SC4' and 
                    cly[1] == 'LY2' and cly[2] == 'SC1'):
                    
                    self.crosslink_bonded['bonds'].append([
                        clx[0], cly[0], '1', self.dlyxly2, f"{self.klyxly2}\n"
                    ])
                    self.crosslink_bonded['angles'].append([
                        str(int(clx[0])+1), clx[0], cly[0], '1', self.al2yx_1, f"{self.k_angle}\n"
                    ])
                    self.crosslink_bonded['angles'].append([
                        str(int(clx[0])-1), clx[0], cly[0], '1', self.al2yx_2, f"{self.k_angle}\n"
                    ])
                    self.crosslink_bonded['angles'].append([
                        clx[0], cly[0], str(int(cly[0])-1), '1', self.al2yx_3, f"{self.k_angle}\n"
                    ])
                    connections_found += 1
                    LOG.info(f"Added LYX-LY2 crosslink between {clx[0]} and {cly[0]} (distance: {dist:.3f} Å)")
                    
                # LYX SC5 - LY3 SC1 crosslinks  
                elif (clx[1] == 'LYX' and clx[2] == 'SC5' and 
                      cly[1] == 'LY3' and cly[2] == 'SC1'):
                    
                    self.crosslink_bonded['bonds'].append([
                        clx[0], cly[0], '1', self.dlyxly3, f"{self.klyxly3}\n"
                    ])
                    self.crosslink_bonded['angles'].append([
                        str(int(clx[0])-1), clx[0], cly[0], '1', self.al3yx_1, f"{self.k_angle}\n"
                    ])
                    self.crosslink_bonded['angles'].append([
                        str(int(clx[0])-2), clx[0], cly[0], '1', self.al3yx_2, f"{self.k_angle}\n"
                    ])
                    self.crosslink_bonded['angles'].append([
                        clx[0], cly[0], str(int(cly[0])-1), '1', self.al3yx_3, f"{self.k_angle}\n"
                    ])
                    connections_found += 1
                    LOG.info(f" Added LYX-LY3 crosslink between {clx[0]} and {cly[0]} (distance: {dist:.3f} Å)")
                    
                # L4Y SC1 - L5Y SC2 crosslinks
                elif (clx[1] == 'L4Y' and clx[2] == 'SC1' and 
                      cly[1] == 'L5Y' and cly[2] == 'SC2'):
                    
                    self.crosslink_bonded['bonds'].append([
                        clx[0], cly[0], '1', self.dly45, f"{self.kly45}\n"
                    ])
                    self.crosslink_bonded['angles'].append([
                        clx[0], cly[0], str(int(cly[0])-1), '1', self.al45y_1, f"{self.k_angle}\n"
                    ])
                    self.crosslink_bonded['angles'].append([
                        str(int(clx[0])-1), clx[0], cly[0], '1', self.al45y_2, f"{self.k_angle}\n"
                    ])
                    connections_found += 1
                    LOG.info(f" Added L4Y-L5Y crosslink between {clx[0]} and {cly[0]} (distance: {dist:.3f} Å)")
                
                # Handle reverse order pairs (cly, clx instead of clx, cly)
                elif (cly[1] == 'LYX' and cly[2] == 'SC4' and 
                      clx[1] == 'LY2' and clx[2] == 'SC1'):
                    
                    self.crosslink_bonded['bonds'].append([
                        cly[0], clx[0], '1', self.dlyxly2, f"{self.klyxly2}\n"
                    ])
                    self.crosslink_bonded['angles'].append([
                        str(int(cly[0])+1), cly[0], clx[0], '1', self.al2yx_1, f"{self.k_angle}\n"
                    ])
                    self.crosslink_bonded['angles'].append([
                        str(int(cly[0])-1), cly[0], clx[0], '1', self.al2yx_2, f"{self.k_angle}\n"
                    ])
                    self.crosslink_bonded['angles'].append([
                        cly[0], clx[0], str(int(clx[0])-1), '1', self.al2yx_3, f"{self.k_angle}\n"
                    ])
                    connections_found += 1
                    LOG.info(f" Added LYX-LY2 crosslink between {cly[0]} and {clx[0]} (distance: {dist:.3f} Å)")
                    
                elif (cly[1] == 'LYX' and cly[2] == 'SC5' and 
                      clx[1] == 'LY3' and clx[2] == 'SC1'):
                    
                    self.crosslink_bonded['bonds'].append([
                        cly[0], clx[0], '1', self.dlyxly3, f"{self.klyxly3}\n"
                    ])
                    self.crosslink_bonded['angles'].append([
                        str(int(cly[0])-1), cly[0], clx[0], '1', self.al3yx_1, f"{self.k_angle}\n"
                    ])
                    self.crosslink_bonded['angles'].append([
                        str(int(cly[0])-2), cly[0], clx[0], '1', self.al3yx_2, f"{self.k_angle}\n"
                    ])
                    self.crosslink_bonded['angles'].append([
                        cly[0], clx[0], str(int(clx[0])-1), '1', self.al3yx_3, f"{self.k_angle}\n"
                    ])
                    connections_found += 1
                    LOG.info(f" Added LYX-LY3 crosslink between {cly[0]} and {clx[0]} (distance: {dist:.3f} Å)")
                    
                elif (cly[1] == 'L4Y' and cly[2] == 'SC1' and 
                      clx[1] == 'L5Y' and clx[2] == 'SC2'):
                    
                    self.crosslink_bonded['bonds'].append([
                        cly[0], clx[0], '1', self.dly45, f"{self.kly45}\n"
                    ])
                    self.crosslink_bonded['angles'].append([
                        cly[0], clx[0], str(int(clx[0])-1), '1', self.al45y_1, f"{self.k_angle}\n"
                    ])
                    self.crosslink_bonded['angles'].append([
                        str(int(cly[0])-1), cly[0], clx[0], '1', self.al45y_2, f"{self.k_angle}\n"
                    ])
                    connections_found += 1
                    LOG.info(f" Added L4Y-L5Y crosslink between {cly[0]} and {clx[0]} (distance: {dist:.3f} Å)")
                else:
                    LOG.debug(f"    Distance {dist:.3f} Å between {clx[1]}{clx[2]} and {cly[1]}{cly[2]} - not a recognized crosslink type")
                
            LOG.info(f"Created {len(self.crosslink_bonded['bonds'])} bonds and "
                    f"{len(self.crosslink_bonded['angles'])} angles from "
                    f"{connections_found} crosslink connections")
                
        except Exception as e:
            LOG.error(f"Error setting up crosslink bonded parameters: {str(e)}")
            
        return self.crosslink_bonded