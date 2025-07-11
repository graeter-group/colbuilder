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
    
    Attributes
    ----------
    file : str
        Path to the merged PDB file
    crosslink_coords : List[List[float]]
        Coordinates of identified crosslink sites
    crosslink_pdb : List[List[str]]
        PDB data for crosslink atoms
    crosslink_neighbors : List[Any]
        Neighboring atoms for each crosslink site
    crosslink_connect : List[List[List[str]]]
        Connected crosslink sites
    crosslink_bonded : Dict[str, List[List[Any]]]
        Bonded parameters for crosslinks
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
        self.crosslink_pdb: List[List[str]] = []
        self.crosslink_neighbors: List[Any] = []
        self.crosslink_connect: List[List[List[str]]] = []
        self.crosslink_bonded: Dict[str, List[List[Any]]] = {
            'bonds': [], 
            'angles': [], 
            'dihedrals': []
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
        
        Parses PDB file to identify atoms involved in crosslinks, including:
        - LYX SC4/SC5 atoms for PYD crosslinks
        - LY2/LY3 SC1 atoms for PYD crosslinks  
        - L4Y SC1 and L5Y SC2 atoms for HLKNL crosslinks
        
        Parameters
        ----------
        cnt_model : Optional[int]
            Model counter for PDB file naming. If None, uses instance file path
        
        Returns
        -------
        List[List[float]]
            List of [x,y,z] coordinates for each crosslink site
        
        Notes
        -----
        Also populates self.crosslink_pdb with atom metadata including:
        atom index, residue name, atom name, chain ID, and coordinates
        """
        file = self.file if cnt_model is None else f"{int(cnt_model)}.merge.pdb"
        LOG.debug(f"Reading crosslink coordinates from {file}")
        
        if not os.path.exists(file):
            LOG.error(f"PDB file not found: {file}")
            return self.crosslink_coords
            
        self.crosslink_coords = []
        self.crosslink_pdb = []
        atom_index = 0
        
        try:
            with open(file, 'r') as f:
                for line in f:
                    if not line.startswith('ATOM'):
                        continue
                        
                    atom_index += 1
                    residue = line[17:20]
                    atom = line[12:15]
                    
                    if ((residue == 'LYX' and atom in ['SC4', 'SC5']) or
                        (residue in ['LY2', 'LY3'] and atom == 'SC1') or
                        (residue == 'L4Y' and atom == 'SC1') or
                        (residue == 'L5Y' and atom == 'SC2')):
                        
                        coords = [
                            float(line[29:38]),
                            float(line[38:46]), 
                            float(line[46:56])
                        ]
                        
                        self.crosslink_pdb.append([
                            str(atom_index),
                            residue,
                            atom, 
                            line[21:26],
                            *[str(c) for c in coords]
                        ])
                        self.crosslink_coords.append(coords)
            
            LOG.debug(f"Found {len(self.crosslink_coords)} crosslink sites")
                
        except Exception as e:
            LOG.error(f"Error reading PDB file {file}: {str(e)}")
            
        return self.crosslink_coords

    def get_crosslink_connect(self, cnt_model: Optional[int] = None) -> List[List[List[str]]]:
        """
        Identify connections between crosslink sites based on proximity.
        
        Uses pairwise distance calculations to find nearest neighbor crosslink 
        sites that could potentially form bonds. For each site, up to 4 nearest
        neighbors are considered as potential connection points.
        
        Parameters
        ----------
        cnt_model : Optional[int]
            Model counter for PDB file naming. If None, uses instance file path
        
        Returns
        -------
        List[List[List[str]]]
            Nested list of crosslink site metadata for connected pairs
        """
        self.get_crosslink_coords(cnt_model=cnt_model)
        
        if not self.crosslink_coords:
            LOG.warning("No crosslink coordinates found")
            return []
        
        try:
            pairs = pdist(np.array(self.crosslink_coords))
            processed_indices = []
            self.crosslink_connect = []
            
            for distances in pairs:
                connections = []
                for idx in np.argsort(distances)[:4]:
                    if idx not in processed_indices:
                        connections.append(self.crosslink_pdb[idx])
                        processed_indices.append(idx)
                        
                if connections:
                    self.crosslink_connect.append(connections)
                    
            LOG.debug(f"Found {len(self.crosslink_connect)} potential connections")
                
        except Exception as e:
            LOG.error(f"Error finding crosslink connections: {str(e)}")
            
        return self.crosslink_connect
    
    def set_crosslink_bonded(self, cnt_model: Optional[int] = None, 
                           crosslink_connect: Optional[List[List[List[str]]]] = None) -> Dict[str, List[List[Any]]]:
        """
        Setup topology for crosslink bonded parameters: bonds, angles, and dihedrals.
        
        Identifies different types of crosslinks and sets up the
        appropriate bonded parameters for each.
        
        Parameters
        ----------
        cnt_model : Optional[int]
            The model counter for file naming if different from initialization
        crosslink_connect : Optional[List[List[List[str]]]]
            List of connected crosslink sites (optional)
            
        Returns
        -------
        Dict[str, List[List[Any]]]
            Dictionary of bonded parameters for crosslinks
        """
        LOG.debug(f"Setting up crosslink bonded parameters for model {cnt_model}")
        
        self.crosslink_bonded = {'bonds': [], 'angles': [], 'dihedrals': []}
        
        if crosslink_connect is None:
            crosslink_connect = self.get_crosslink_connect(cnt_model=cnt_model)
            
        if not crosslink_connect:
            LOG.debug("No crosslink connections found, returning empty parameters")
            return self.crosslink_bonded
            
        try:
            connections_found = 0
            for c in crosslink_connect:
                for i, clx in enumerate(c):
                    for j, cly in enumerate(c):
                        if i == j:
                            continue
                            
                        try:
                            dist = np.linalg.norm(np.array(clx[-3:]) - np.array(cly[-3:]))
                        except (ValueError, TypeError) as e:
                            LOG.warning(f"Error calculating distance: {str(e)}")
                            continue
                            
                        if dist < 10.0:
                            if clx[1] == 'LYX' and clx[2] == 'SC4' and cly[1] == 'LY2':
                                self._add_lyx_ly2_bonds(clx, cly)
                                connections_found += 1
                                LOG.debug(f"Added LYX-LY2 bonds between {clx[0]} and {cly[0]}")
                                
                            elif clx[1] == 'LYX' and clx[2] == 'SC5' and cly[1] == 'LY3':
                                self._add_lyx_ly3_bonds(clx, cly)
                                connections_found += 1
                                LOG.debug(f"Added LYX-LY3 bonds between {clx[0]} and {cly[0]}")
                                
                            if clx[1] == 'L4Y' and clx[2] == 'SC1' and cly[1] == 'L5Y' and cly[2] == 'SC2':
                                self._add_l4y_l5y_bonds(clx, cly)
                                connections_found += 1
                                LOG.debug(f"Added L4Y-L5Y bonds between {clx[0]} and {cly[0]}")
                    
            LOG.debug(f"Created {len(self.crosslink_bonded['bonds'])} bonds and "
                    f"{len(self.crosslink_bonded['angles'])} angles from "
                    f"{connections_found} crosslink connections")
                
        except Exception as e:
            LOG.error(f"Error setting up crosslink bonded parameters: {str(e)}")
            
        return self.crosslink_bonded
    
    def _add_lyx_ly2_bonds(self, clx: List[Any], cly: List[Any]) -> None:
        """
        Add bonded parameters for LYX-LY2 crosslink.
        
        Parameters
        ----------
        clx : List[Any]
            First crosslink atom data
        cly : List[Any]
            Second crosslink atom data
        """
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
    
    def _add_lyx_ly3_bonds(self, clx: List[Any], cly: List[Any]) -> None:
        """
        Add bonded parameters for LYX-LY3 crosslink.
        
        Parameters
        ----------
        clx : List[Any]
            First crosslink atom data
        cly : List[Any]
            Second crosslink atom data
        """
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
    
    def _add_l4y_l5y_bonds(self, clx: List[Any], cly: List[Any]) -> None:
        """
        Add bonded parameters for L4Y-L5Y crosslink.
        
        Parameters
        ----------
        clx : List[Any]
            First crosslink atom data
        cly : List[Any]
            Second crosslink atom data
        """
        self.crosslink_bonded['bonds'].append([
            clx[0], cly[0], '1', self.dly45, f"{self.kly45}\n"
        ])
        
        self.crosslink_bonded['angles'].append([
            clx[0], cly[0], str(int(cly[0])-1), '1', self.al45y_1, f"{self.k_angle}\n"
        ])
        self.crosslink_bonded['angles'].append([
            str(int(clx[0])-1), clx[0], cly[0], '1', self.al45y_2, f"{self.k_angle}\n"
        ])