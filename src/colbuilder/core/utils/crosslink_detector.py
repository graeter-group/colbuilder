"""
Crosslink detection utilities for inferring structure types from PDB files.

This module provides functionality to detect crosslink types in collagen
structures by analyzing residue names in PDB files. It's used primarily
for topology-only mode to determine the appropriate force field parameters.
"""

from pathlib import Path
from typing import Set, List, Dict, Tuple, Optional
import shutil
import numpy as np

from colbuilder.core.utils.logger import setup_logger

LOG = setup_logger(__name__)


class CrosslinkDetector:
    """
    Detect and classify crosslink types in collagen structures.
    
    This class analyzes PDB files to identify crosslink residues and infer
    the crosslink type (D=divalent, T=trivalent, or M=mixed) based on the 
    residues present in the structure.
    """
    
    # Divalent crosslink residues (D-type)
    DIVALENT_RESIDUES = {
        'LY4', 'LX4', 'L4Y', 'L4X',
        'LY5', 'LX5', 'L5Y', 'L5X',
        'LGX', 'LPS',
        'AGS', 'APD',
    }
    
    # Trivalent crosslink residues (T-type)
    TRIVALENT_RESIDUES = {
        'LYX', 'LXY', 'LYY', 'LXX',
        'LY2', 'LX2',
        'LY3', 'LX3', 'L3Y', 'L2Y', 'L3X', 'L2X',
    }
    
    ALL_CROSSLINK_RESIDUES = DIVALENT_RESIDUES | TRIVALENT_RESIDUES
    
    @classmethod
    def detect_structure_type(cls, pdb_path: Path) -> str:
        """
        Detect the overall crosslink type from a PDB file.
        
        Analyzes the entire structure (regardless of MODEL records) to determine
        what type of crosslinks are present.
        
        Args:
            pdb_path: Path to PDB file
            
        Returns:
            Structure type: 'D' (divalent), 'T' (trivalent), 'M' (mixed), or 'D' (default)
        """
        residues = cls.get_all_residues(pdb_path)
        
        if not residues:
            LOG.warning(f"No residues found in {pdb_path}")
            return 'D'
        
        divalent_found = bool(residues & cls.DIVALENT_RESIDUES)
        trivalent_found = bool(residues & cls.TRIVALENT_RESIDUES)
        
        if divalent_found and trivalent_found:
            LOG.info("Mixed crosslinks detected (both divalent and trivalent)")
            return 'M'  # Mixed
        elif trivalent_found:
            LOG.info("Trivalent crosslinks detected")
            return 'T'
        elif divalent_found:
            LOG.info("Divalent crosslinks detected")
            return 'D'
        else:
            # No crosslinks found - could be a non-crosslinked structure
            LOG.info("No crosslinks detected, using default type D")
            return 'D'  # Default to divalent
    
    @classmethod
    def get_all_residues(cls, pdb_path: Path) -> Set[str]:
        """
        Extract all unique residue names from a PDB file.
        
        Scans the entire PDB file (all models, all chains) and returns
        the set of all unique residue names found.
        
        Args:
            pdb_path: Path to PDB file
            
        Returns:
            Set of unique residue names (3-letter codes)
        """
        residues = set()
        
        try:
            with open(pdb_path, 'r') as f:
                for line in f:
                    if line.startswith(('ATOM', 'HETATM')):
                        if len(line) >= 20:
                            resname = line[17:20].strip()
                            residues.add(resname)
            
            crosslinks_found = residues & cls.ALL_CROSSLINK_RESIDUES
            if crosslinks_found:
                LOG.debug(f"Crosslink residues present: {sorted(crosslinks_found)}")
            
            return residues
            
        except Exception as e:
            LOG.error(f"Error reading PDB file {pdb_path}: {e}")
            return set()
    
    @classmethod
    def get_model_ids(cls, pdb_path: Path) -> List[int]:
        """
        Extract all model IDs from a PDB file.
        
        Args:
            pdb_path: Path to PDB file
            
        Returns:
            List of model IDs found in the file (sorted)
        """
        model_ids = []
        
        try:
            with open(pdb_path, 'r') as f:
                for line in f:
                    if line.startswith('MODEL'):
                        try:
                            model_id = int(line.split()[1])
                            model_ids.append(model_id)
                        except (IndexError, ValueError):
                            continue
            
            return sorted(model_ids)
            
        except Exception as e:
            LOG.error(f"Error reading PDB file {pdb_path}: {e}")
            return []
    
    @classmethod
    def has_model_records(cls, pdb_path: Path) -> bool:
        """
        Check if a PDB file contains MODEL/ENDMDL records.
        
        Args:
            pdb_path: Path to PDB file
            
        Returns:
            True if MODEL records are present, False otherwise
        """
        try:
            with open(pdb_path, 'r') as f:
                for line in f:
                    if line.startswith('MODEL'):
                        return True
            return False
        except Exception as e:
            LOG.error(f"Error checking for MODEL records: {e}")
            return False
    
    @classmethod
    def detect_crosslinks_in_models(cls, model_dir: Path) -> Dict[int, List[Tuple[str, int, np.ndarray]]]:
        """
        Detect crosslink positions in all extracted model PDB files.
        
        Returns a dictionary mapping model_id to list of (resname, atom_num, position) tuples
        for crosslink residues.
        
        Args:
            model_dir: Directory containing extracted model PDB files
            
        Returns:
            Dictionary: {model_id: [(resname, atom_num, xyz_coords), ...]}
        """
        crosslink_data = {}
        
        model_files = sorted(model_dir.glob("*.pdb"))
        
        for model_file in model_files:
            if model_file.stem.endswith('.caps'):
                continue
            
            try:
                model_id = int(model_file.stem)
            except ValueError:
                continue
            
            crosslinks = []
            
            with open(model_file, 'r') as f:
                for line in f:
                    if not line.startswith(('ATOM', 'HETATM')):
                        continue
                    
                    resname = line[17:20].strip()
                    
                    if resname in cls.ALL_CROSSLINK_RESIDUES:
                        try:
                            atom_num = int(line[6:11])
                            x = float(line[30:38])
                            y = float(line[38:46])
                            z = float(line[46:54])
                            position = np.array([x, y, z])
                            
                            crosslinks.append((resname, atom_num, position))
                        except (ValueError, IndexError):
                            continue
            
            if crosslinks:
                crosslink_data[model_id] = crosslinks
                LOG.debug(f"Model {model_id}: found {len(crosslinks)} crosslink atoms")
        
        return crosslink_data

    @classmethod
    def find_crosslink_connections(
        cls, 
        model_dir: Path, 
        cutoff: float = 5.0
    ) -> Dict[int, Set[int]]:
        """
        Find connections between models based on crosslink proximity.
        
        Two models are connected if any of their crosslink atoms are within
        the cutoff distance.
        
        Args:
            model_dir: Directory containing extracted model PDB files
            cutoff: Distance cutoff in Angstroms for crosslink connection
            
        Returns:
            Dictionary mapping model_id to set of connected model_ids
        """
        crosslink_data = cls.detect_crosslinks_in_models(model_dir)
        
        connections = {model_id: set([model_id]) for model_id in crosslink_data.keys()}

        model_ids = sorted(crosslink_data.keys())
        
        # Find inter-model connections
        for i, model_a in enumerate(model_ids):
            for model_b in model_ids[i+1:]:
                for resname_a, atom_a, pos_a in crosslink_data[model_a]:
                    for resname_b, atom_b, pos_b in crosslink_data[model_b]:
                        distance = np.linalg.norm(pos_a - pos_b)
                        
                        if distance < cutoff:
                            connections[model_a].add(model_b)
                            connections[model_b].add(model_a)
                            LOG.debug(
                                f"Connection found: Model {model_a} ({resname_a}) <-> "
                                f"Model {model_b} ({resname_b}), distance: {distance:.2f} Å"
                            )
                            # Break after finding first connection between these models
                            break
                    else:
                        continue
                    break
        
        connected_count = len([c for c in connections.values() if len(c) > 1])
        LOG.info(f"Found connections for {connected_count} models")
        
        return connections

    @classmethod
    def split_by_ter_records(cls, pdb_path: Path, output_dir: Path) -> List[Path]:
        """
        Split a fibril PDB into individual models using TER records.
        
        Assumes each triple helix model consists of three chains (A, B, C)
        and ends with a TER record after chain C.
        
        Args:
            pdb_path: Path to fibril PDB file
            output_dir: Directory where extracted models will be saved
            
        Returns:
            List of paths to extracted model files
        """
        LOG.info("  Splitting PDB by TER records (triple helix pattern)")
        
        extracted_files = []
        model_num = 0
        current_model_lines = []
        ter_count = 0
        
        try:
            with open(pdb_path, 'r') as f:
                for line in f:
                    if line.startswith(('CRYST1', 'END', 'REMARK', 'TITLE', 'HEADER')):
                        continue
                    
                    if line.startswith(('ATOM', 'HETATM', 'TER', 'ANISOU')):
                        current_model_lines.append(line)
                        
                        if line.startswith('TER'):
                            ter_count += 1
                            
                            # Every 3 TERs = complete triple helix model
                            if ter_count % 3 == 0:
                                model_file = output_dir / f"{model_num}.pdb"
                                with open(model_file, 'w') as f_out:
                                    f_out.writelines(current_model_lines)
                                    f_out.write("END\n")
                                
                                extracted_files.append(model_file)
                                LOG.debug(f"Extracted model {model_num} ({len(current_model_lines)} lines)")
                                
                                current_model_lines = []
                                model_num += 1
                
                if current_model_lines:
                    model_file = output_dir / f"{model_num}.pdb"
                    with open(model_file, 'w') as f_out:
                        f_out.writelines(current_model_lines)
                        f_out.write("END\n")
                    extracted_files.append(model_file)
                    LOG.debug(f"Extracted final model {model_num} ({len(current_model_lines)} lines)")
            
            LOG.info(f"Successfully extracted {len(extracted_files)} models to {output_dir}")
            return extracted_files
            
        except Exception as e:
            LOG.error(f"Error splitting PDB by TER records: {e}")
            return []
    
    @classmethod
    def extract_models_by_model_records(cls, pdb_path: Path, output_dir: Path) -> List[Path]:
        """
        Extract individual models from a PDB file using MODEL/ENDMDL records.
        
        Args:
            pdb_path: Path to multi-model PDB file
            output_dir: Directory where extracted models will be saved
            
        Returns:
            List of paths to extracted model files
        """
        LOG.info(f"Extracting models using MODEL/ENDMDL records")
        
        extracted_files = []
        current_model = None
        current_lines = []
        
        try:
            with open(pdb_path, 'r') as f:
                for line in f:
                    if line.startswith('MODEL'):
                        try:
                            current_model = int(line.split()[1])
                            current_lines = []
                        except (IndexError, ValueError):
                            LOG.warning(f"Could not parse MODEL line: {line.strip()}")
                            continue
                    
                    elif line.startswith('ENDMDL'):
                        if current_model is not None and current_lines:
                            model_file = output_dir / f"{current_model}.pdb"
                            with open(model_file, 'w') as f_out:
                                f_out.writelines(current_lines)
                                f_out.write("END\n")
                            
                            extracted_files.append(model_file)
                            LOG.debug(f"Extracted model {current_model} to {model_file.name}")
                        
                        current_model = None
                        current_lines = []
                    
                    elif current_model is not None:
                        if line.startswith(('ATOM', 'HETATM', 'TER', 'ANISOU')):
                            current_lines.append(line)
            
            LOG.info(f"Successfully extracted {len(extracted_files)} models to {output_dir}")
            return extracted_files
            
        except Exception as e:
            LOG.error(f"Error extracting models from PDB: {e}")
            return []
    
    @classmethod
    def extract_models(
        cls, 
        pdb_path: Path, 
        structure_type: str,
        output_dir: Path
    ) -> List[Path]:
        """
        Extract individual models from a PDB file.
        
        Automatically detects whether to use MODEL/ENDMDL records or TER-based
        splitting. Results are saved in a type-specific subdirectory.
        
        Args:
            pdb_path: Path to PDB file
            structure_type: Pre-detected structure type (D/T/M)
            output_dir: Base directory where type-specific subdirectory will be created
            
        Returns:
            List of paths to extracted model files
        """
        # Create type-specific subdirectory
        type_dir = output_dir / structure_type
        type_dir.mkdir(parents=True, exist_ok=True)
        
        # Check for MODEL records
        has_models = cls.has_model_records(pdb_path)
        
        if has_models:
            return cls.extract_models_by_model_records(pdb_path, type_dir)
        else:
            LOG.info("No MODEL records found - using TER-based splitting")
            return cls.split_by_ter_records(pdb_path, type_dir)
    
    @classmethod
    def create_caps_files(cls, model_dir: Path) -> None:
        """
        Create .caps.pdb copies of all model files in a directory.
        
        These files are used by certain downstream processing steps.
        
        Args:
            model_dir: Directory containing model PDB files
        """
        model_files = sorted(model_dir.glob("[0-9]*.pdb"))
        
        for model_file in model_files:
            # Skip if it's already a caps file
            if model_file.stem.endswith('.caps'):
                continue
            
            caps_file = model_file.parent / f"{model_file.stem}.caps.pdb"
            
            try:
                shutil.copy2(model_file, caps_file)
                LOG.debug(f"Created caps file: {caps_file.name}")
            except Exception as e:
                LOG.warning(f"Could not create caps file for {model_file.name}: {e}")
    
    @classmethod
    def prepare_for_topology(
        cls, 
        pdb_path: Path, 
        output_dir: Path
    ) -> Path:
        """
        Prepare a PDB structure for topology generation.
        
        This is the main entry point for topology-only mode. It:
        1. Detects the crosslink type once
        2. Extracts individual models
        3. Creates caps files
        4. Returns the path to the prepared structure directory
        
        Args:
            pdb_path: Path to input PDB file
            output_dir: Directory where type-specific subdirectory will be created
            
        Returns:
            Path to the type-specific directory containing extracted models
        """
        structure_type = cls.detect_structure_type(pdb_path)
        
        extracted_files = cls.extract_models(pdb_path, structure_type, output_dir)
        
        if not extracted_files:
            LOG.error("No models were extracted from the PDB file")
            return output_dir / structure_type
        
        type_dir = output_dir / structure_type
        cls.create_caps_files(type_dir)
        
        LOG.info(f"Prepared {len(extracted_files)} models in {type_dir}")
        
        return type_dir