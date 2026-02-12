# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0

from __future__ import annotations
from typing import Dict, List, Optional, Tuple, Any, Union
from pathlib import Path
import numpy as np
import os
import shutil

from colbuilder.core.utils.logger import setup_logger

LOG = setup_logger(__name__)

class System:
    """
    Base class representing a system of models.

    The system is represented as a dictionary where keys are model IDs and values are model objects.

    Attributes
    ----------
    system : Dict[float, Any]
        Dictionary of models in the system.
    connect : Dict[float, List[float]]
        Dictionary of connections between models.
    models : List[float]
        List of model IDs in the system.
    crystal : Any
        Crystal information for the system.
    crystalcontacts : Any
        Crystal contacts information for the system.
    size_models : int
        Number of models in the system.
    is_line : Tuple[str, ...]
        Tuple of valid line types in PDB files.
    size : int
        Total number of models in the system.
    type : str
        Type of the system.
    pdb_fibril : Path
        Path to the PDB file of the fibril.
    """

    def __init__(self, crystal: Optional[Any] = None, crystalcontacts: Optional[Any] = None, pdb_fibril: Path = Path()):
        self.system: Dict[float, Any] = {}
        self.connect: Dict[float, List[float]] = {}
        self.models: List[float] = []
        self.crystal = crystal
        self.crystalcontacts = crystalcontacts
        self.size_models: int = 0
        self.is_line: Tuple[str, ...] = ('ATOM  ', 'HETATM', 'ANISOU', 'TER   ')
        self.size: int = 0
        self.type: str = ''
        self.pdb_fibril: Path = pdb_fibril

    def add_model(self, model: Any) -> None:
        """
        Add a model to the system.

        Parameters
        ----------
        model : Any
            The model to be added to the system.
        """
        if hasattr(model, 'id') and hasattr(model, 'type') and hasattr(model, 'connect'):
            self.system[model.id] = model
        else:
            LOG.error(f"Model missing required attributes: {model}")

    def set_crystal(self, crystal: Optional[Any] = None) -> None:
        """
        Set crystal information for all models in the system.

        Parameters
        ----------
        crystal : Optional[Any], default=None
            Crystal information to be set. If None, uses the system's crystal.
        """
        if crystal is None:
            crystal = self.crystal
        for model in self.system.values():
            model.crystal = crystal

    def get_size(self) -> int:
        """
        Get the number of models in the system.

        Returns
        -------
        int
            The number of models in the system.
        """
        self.size = len(self.system)
        return self.size

    def get_connect_size(self) -> int:
        """
        Get the number of connected models in the system.

        Returns
        -------
        int
            The number of connected models in the system.
        """
        return sum(1 for model in self.models if self.get_model(model_id=model).connect is not None)

    def get_model(self, model_id: float) -> Any:
        """
        Get a model from the system by its ID.

        Parameters
        ----------
        model_id : float
            The ID of the model to retrieve.

        Returns
        -------
        Any
            The model with the specified ID.
        """
        return self.system[model_id]

    def get_models(self) -> List[float]:
        """
        Get a list of all model IDs in the system.

        Returns
        -------
        List[float]
            A list of all model IDs in the system.
        """
        self.models = list(self.system.keys())
        return self.models

    def get_connect(self) -> Dict[float, List[float]]:
        """
        Get the connections for all models in the system.

        Returns
        -------
        Dict[float, List[float]]
            A dictionary of connections for each model.
        """
        self.connect = {model.id: model.connect for model in self.system.values()}
        return self.connect

    def delete_model(self, model_id: float) -> Dict[float, Any]:
        """
        Delete a model from the system.

        Parameters
        ----------
        model_id : float
            The ID of the model to delete.

        Returns
        -------
        Dict[float, Any]
            The updated system dictionary.
        """
        del self.system[model_id]
        return self.system

    def translate_system(self, crystal: Any, center: List[float]) -> None:
        """
        Translate the whole system to a certain position.

        Parameters
        ----------
        crystal : Any
            Crystal information for the translation.
        center : List[float]
            The target center position [x, y, z].
        """
        translate = [0, 0, center[2] - self.center_system(crystal=crystal)]
        for model in self.system.values():
            if model.connect is not None:
                for connect_id in model.connect:
                    crystal.translate_crystal(
                        pdb=Path(model.type) / f"{int(connect_id)}.caps.pdb",
                        translate=translate,
                        bool_system=True
                    )

    def center_system(self, crystal: Any) -> float:
        """
        Calculate the center of the system.

        Parameters
        ----------
        crystal : Any
            Crystal information for the calculation.

        Returns
        -------
        float
            The z-coordinate of the system's center.
        """
        cog = []
        for model_id in self.get_models():
            if self.get_model(model_id=model_id).connect is not None:
                for connect_id in self.get_model(model_id=model_id).connect:
                    full_path = f"{self.get_model(model_id=model_id).type}/{int(connect_id)}.caps.pdb"
                    try:
                        cog.append(crystal.get_cog(pdb=full_path))
                    except FileNotFoundError:
                        LOG.error(f"File not found: {full_path}.pdb")
                        continue
        return np.mean(cog)

    def count_states(self, state: str) -> int:
        """
        Count all models with a certain state.

        Parameters
        ----------
        state : str
            The state to count ('no', 'mut', or 'prot').

        Returns
        -------
        int
            The number of models with the specified state.
        """
        return sum(model.count_state(state=state) for model in self.system.values())
    
    def safe_remove_directory(self, directory: Union[str, Path]) -> None:
        """
        Safely remove a directory and all its contents.
        
        Parameters
        ----------
        directory : Union[str, Path]
            The directory to remove
        """
        try:
            directory_path = Path(directory).resolve() 
            
            if not directory_path.exists():
                LOG.debug(f"Directory does not exist, nothing to remove: {directory_path}")
                return
                
            if not any(directory_path.name == x for x in ['NC', 'T', 'D', 'TD', 'DT', 'M']):
                LOG.warning(f"Refusing to remove directory that isn't a type directory: {directory_path}")
                return
                
            shutil.rmtree(directory_path)
            LOG.debug(f"Successfully removed directory: {directory_path}")
        except Exception as e:
            LOG.warning(f"Failed to remove directory {directory_path}: {str(e)}")
            
    def write_pdb(self, pdb_out: Union[str, Path], fibril_length: float, cleanup: bool = True, temp_dir: Optional[Path] = None):
        """
        Write the system to a PDB file with proper type-specific caps handling.
        
        This method writes the system in a deterministic order: one caps file per model ID,
        respecting each model's assigned type. This ensures consistent output regardless of
        whether connectivity information is present.
        
        Args:
            pdb_out: Path to output PDB file
            fibril_length: Length of fibril in nm
            cleanup: Whether to clean up temporary files
            temp_dir: Directory containing caps files to use
            
        Raises:
            ValueError: If temp_dir is not provided
            FileNotFoundError: If required type directories are not found
        """
        pdb_out_path = Path(pdb_out)
        if not pdb_out_path.is_absolute():
            pdb_out_path = Path.cwd() / pdb_out_path
        pdb_out_path.parent.mkdir(exist_ok=True, parents=True)
        if pdb_out_path.suffix != '.pdb':
            pdb_out_path = pdb_out_path.with_suffix('.pdb')

        if temp_dir is None:
            raise ValueError("temp_dir must be provided - specify which directory to use for caps files")
        
        # Collect model types and organize models by type
        model_types = set()
        model_type_count = {}
        models_by_type: Dict[str, List[float]] = {}
        
        for model_id, model in self.system.items():
            if hasattr(model, 'type') and model.type:
                model_types.add(model.type)
                model_type_count[model.type] = model_type_count.get(model.type, 0) + 1
                models_by_type.setdefault(model.type, []).append(model_id)

        LOG.debug(f"Model type counts: {model_type_count}")
        
        # Find type directories and verify they exist
        type_dirs = {}
        for model_type in model_types:
            type_dir = temp_dir / str(model_type)
            if type_dir.exists():
                type_dirs[model_type] = type_dir
                caps_files_count = len(list(type_dir.glob("*.caps.pdb")))
                LOG.debug(f"Found {caps_files_count} caps files in {type_dir}")
            else:
                LOG.error(f"No {model_type} directory found in {temp_dir}!")
                raise FileNotFoundError(f"Required directory not found: {type_dir}")
        
        # Build index of available caps files
        caps_by_model = {}
        caps_by_type = {}
        for model_type, type_dir in type_dirs.items():
            caps_by_type[model_type] = []
            for caps_file in type_dir.glob("*.caps.pdb"):
                try:
                    connect_id = float(caps_file.stem.split('.')[0])
                    caps_by_model[(model_type, connect_id)] = caps_file
                    caps_by_type[model_type].append(connect_id)
                except (ValueError, IndexError):
                    LOG.warning(f"Could not parse caps file name: {caps_file}")
                    continue
        
        # Log connectivity status for debugging
        has_connections = False
        models_with_connections = 0
        models_by_type_with_connections = {}
        for model_id, model in self.system.items():
            if getattr(model, 'connect', None): 
                has_connections = True
                models_with_connections += 1
                model_type = getattr(model, 'type', 'unknown')
                models_by_type_with_connections[model_type] = models_by_type_with_connections.get(model_type, 0) + 1
        
        if has_connections:
            LOG.debug(f"System has connections: {models_with_connections} models with connectivity")
        else:
            LOG.debug("System has no connectivity information")
        
        try:
            # Get crystal header if available
            crystal_header = None
            if self.crystal and self.crystal.pdb_file:
                crystal_pdb = Path(self.crystal.pdb_file).with_suffix('.pdb')
                if crystal_pdb.exists():
                    with open(crystal_pdb, 'r') as crystal_file:
                        crystal_header = crystal_file.readline()
            
            written_models = 0
            written_by_type = {}
            written_connections = 0
            written_caps_files = set() 
            
            with open(pdb_out_path, 'w') as f:
                # Write crystal header
                if crystal_header:
                    f.write(crystal_header)
                
                # Write model files if they exist (typically for initial crystal structures)
                for model_id, model in self.system.items():
                    if hasattr(model, 'pdb_file') and Path(model.pdb_file).exists():
                        LOG.debug(f"Writing model {model_id} from pdb_file")
                        self._write_model_file_content(f, model.pdb_file)
                        written_models += 1
                        model_type = getattr(model, 'type', 'unknown')
                        written_by_type[model_type] = written_by_type.get(model_type, 0) + 1
                
                # Write caps files: one per model ID, respecting assigned type
                # This deterministic approach works for both connected and non-connected systems
                for model_type, ids in models_by_type.items():
                    for model_id in sorted(ids):
                        key = (model_type, float(model_id))
                        
                        # Skip if already written
                        if key in written_caps_files:
                            continue
                        
                        # Find the caps file
                        caps_file = caps_by_model.get(key)
                        if not caps_file:
                            LOG.warning(
                                f"No caps file found for model {model_id} type {model_type}"
                            )
                            continue
                        
                        # Write the caps file content
                        LOG.debug(f"Writing caps file {caps_file} for model {model_id} type {model_type}")
                        self._write_caps_file_content(f, caps_file)
                        written_connections += 1
                        written_by_type[model_type] = written_by_type.get(model_type, 0) + 1
                        written_caps_files.add(key)
                
                f.write("END\n")
            
            # Verify the output
            self._verify_output_file(pdb_out_path)
            
            LOG.debug(f"PDB writing complete: {written_connections} caps files written")
            LOG.debug(f"Written by type: {written_by_type}")
            
        except Exception as e:
            LOG.error(f"Error writing PDB file {pdb_out_path}: {str(e)}")
            import traceback
            LOG.debug(f"Traceback: {traceback.format_exc()}")
            raise

    def _write_caps_file_content(self, file_handle, caps_file):
        """
        Helper method to write caps file content to an open file handle.
        
        Filters out CRYST1, END, and ENDMDL lines, converts HETATM to ATOM,
        and ensures lines don't exceed 80 characters.
        
        Args:
            file_handle: Open file handle to write to
            caps_file: Path to caps file to read from
        """
        try:
            with open(caps_file, 'r') as caps_file_obj:
                content_written = False
                for line in caps_file_obj:
                    if line.startswith(self.is_line) or line.startswith("TER"):
                        content_written = True
                        # Convert HETATM to ATOM
                        if line.startswith('HETATM'):
                            line = 'ATOM  ' + line[6:]
                        # Ensure line has content and proper length
                        if len(line.rstrip()) > 0:
                            if len(line) > 81:
                                line = line[:80] + '\n'
                            file_handle.write(line)
                    elif line.startswith(("END", "ENDMDL")):
                        break
                
                if not content_written:
                    LOG.warning(f"No atom records found in caps file: {caps_file}")
        except Exception as e:
            LOG.error(f"Error reading caps file {caps_file}: {str(e)}")

    def _write_model_file_content(self, file_handle, model_file):
        """
        Helper method to write model file content to an open file handle.
        
        Similar to _write_caps_file_content but for initial model PDB files.
        
        Args:
            file_handle: Open file handle to write to
            model_file: Path to model file to read from
        """
        try:
            with open(model_file, 'r') as model_file_obj:
                for line in model_file_obj:
                    if line.startswith(self.is_line) or line.startswith("TER"):
                        # Convert HETATM to ATOM
                        if line.startswith('HETATM'):
                            line = 'ATOM  ' + line[6:]
                        # Ensure line has content and proper length
                        if len(line.rstrip()) > 0:
                            if len(line) > 81:
                                line = line[:80] + '\n'
                            file_handle.write(line)
                    elif line.startswith(("END", "ENDMDL")):
                        break
        except Exception as e:
            LOG.error(f"Error reading model file {model_file}: {str(e)}")

    def _verify_output_file(self, pdb_out_path):
        """
        Helper method to verify the output file has content and valid structure.
        
        Checks file size, counts atoms and TER records, and tallies residue types.
        Logs warnings for suspicious files but doesn't raise exceptions.
        
        Args:
            pdb_out_path: Path to the output PDB file to verify
        """
        try:
            if pdb_out_path.exists():
                file_size = pdb_out_path.stat().st_size
                
                if file_size == 0:
                    LOG.error("Output PDB file is empty!")
                elif file_size < 1000:
                    LOG.warning(f"Output PDB file is very small ({file_size} bytes)")
                
                atom_count = 0
                ter_count = 0
                resi_types = {}
                with open(pdb_out_path, 'r') as check_file:
                    for line in check_file:
                        if line.startswith(("ATOM", "HETATM")):
                            atom_count += 1
                            if len(line) >= 20:
                                resi_type = line[17:20].strip()
                                resi_types[resi_type] = resi_types.get(resi_type, 0) + 1
                        elif line.startswith("TER"):
                            ter_count += 1
                
                LOG.debug(f"PDB verification: {atom_count} atoms, {ter_count} TER records")
                if resi_types:
                    LOG.debug(f"Residue type summary: {dict(list(resi_types.items())[:10])}")
            else:
                LOG.error(f"Output PDB file not found: {pdb_out_path}")
        except Exception as e:
            LOG.warning(f"Error performing PDB verification: {e}")
            