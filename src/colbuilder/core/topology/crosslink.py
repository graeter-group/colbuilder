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
        
        Following the working version pattern.
        """
        self.get_crosslink_coords(cnt_model=cnt_model)
        
        if not self.crosslink_coords:
            LOG.warning("No crosslink coordinates found")
            return []
        
        if len(self.crosslink_coords) < 2:
            LOG.warning("Need at least 2 crosslink sites to form connections")
            return []
        
        try:
            pairs = pdist(self.crosslink_coords)
            out = []
            self.crosslink_connect = []
            
            for p in pairs:
                tmp = []
                for k in np.argsort(p)[0:4]:
                    if k not in out:
                        tmp.append(self.crosslink_pdb[k])
                        out.append(k)
                        
                if tmp:
                    self.crosslink_connect.append(tmp)
                    LOG.debug(f"Found connection group with {len(tmp)} atoms")
            
            LOG.debug(f"Found {len(self.crosslink_connect)} potential crosslink connections")
                
        except Exception as e:
            LOG.error(f"Error finding crosslink connections: {str(e)}")
            self.crosslink_connect = []
            
        return self.crosslink_connect
    
    def set_crosslink_bonded(self, cnt_model: Optional[int] = None, 
                           crosslink_connect: Optional[List[List[List[Any]]]] = None) -> Dict[str, List[List[Any]]]:
        """
        Setup topology for crosslink bonded parameters: bonds, angles and dihedrals.
        
        Following the working version pattern: uses atom indices directly.
        """
        LOG.debug(f"Setting up crosslink bonded parameters for model {cnt_model}")
        
        self.crosslink_bonded = {'bonds': [], 'angles': [], 'dihedrals': []}
        
        if crosslink_connect is None:
            crosslink_connect = self.get_crosslink_connect(cnt_model=cnt_model)
            
        if not crosslink_connect:
            LOG.debug("No crosslink connections found, returning empty parameters")
            return self.crosslink_bonded
            
        connections_found = 0
        
        try:
            for c in crosslink_connect:
                for clx in c:
                    for cly in c:
                        if clx == cly:  # Skip self-comparison
                            continue
                            
                        dist = np.linalg.norm(np.array(clx[-3:]) - np.array(cly[-3:]))
                        
                        LOG.debug(f"Distance between {clx[1]}{clx[2]} and {cly[1]}{cly[2]}: {dist:.3f} Å")
                        
                        if dist < 10.0:  # Within crosslinking distance
                            # LYX-LY2 crosslinks
                            if (clx[1] == 'LYX' and clx[2] == 'SC4' and cly[1] == 'LY2'):
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
                                LOG.debug(f"Added LYX-LY2 crosslink between {clx[0]} and {cly[0]} (distance: {dist:.3f} Å)")
                                
                            # LYX-LY3 crosslinks
                            elif (clx[1] == 'LYX' and clx[2] == 'SC5' and cly[1] == 'LY3'):
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
                                LOG.debug(f"Added LYX-LY3 crosslink between {clx[0]} and {cly[0]} (distance: {dist:.3f} Å)")
                                
                            # L4Y-L5Y crosslinks
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
                                LOG.debug(f"Added L4Y-L5Y crosslink between {clx[0]} and {cly[0]} (distance: {dist:.3f} Å)")
                        else:
                            LOG.debug(f"Distance {dist:.3f} Å too large for crosslinking between {clx[1]}{clx[2]} and {cly[1]}{cly[2]}")
                
            LOG.debug(f"Created {len(self.crosslink_bonded['bonds'])} bonds and "
                    f"{len(self.crosslink_bonded['angles'])} angles from "
                    f"{connections_found} crosslink connections")
                
        except Exception as e:
            LOG.error(f"Error setting up crosslink bonded parameters: {str(e)}")
            
        return self.crosslink_bonded