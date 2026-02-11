# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0

from pathlib import Path
import numpy as np
from typing import Dict, List, Optional, Any, Union

from colbuilder.core.utils.logger import setup_logger

LOG = setup_logger(__name__)

class Mix:
    """
    Manages the mixing and type assignment of components in a collagen system.

    This class handles:
    - Assignment of crosslink types according to specified ratios
    - Preservation of spatial relationships between models
    - Management of connected components
    - Type distribution based on component connectivity

    Attributes
    ----------
    ratio_mix : Dict[str, int]
        Mapping of component types to their ratios (e.g., {"D": 80, "T": 20})
    system : Any
        System object containing models to be mixed
    connect_mix : Dict[float, str]
        Mapping of model IDs to their assigned types

    Methods
    -------
    add_mix(ratio_mix, system)
        Assigns types to models according to ratios while preserving connectivity
    get_mix(ratio_mix)
        Returns a random component type based on ratios
    get_mix_from_connect_file(system, connect_file)
        Assigns types based on a connect file specification

    Example
    -------
    >>> mixer = Mix(ratio_mix={"D": 80, "T": 20})
    >>> mixer.add_mix(system=collagen_system)
    >>> assigned_system = mixer.get_mix_from_connect_file("connects.txt")
    """

    def __init__(self, ratio_mix: Optional[Dict[str, int]] = None, system: Optional[Any] = None, connect_mix: Optional[Dict[float, str]] = None):
        """
        Initialize the Mix object.

        Args:
            ratio_mix (Optional[Dict[str, int]]): Initial ratio mix dictionary.
            system (Optional[Any]): Initial system object.
            connect_mix (Optional[Dict[float, str]]): Initial connect mix dictionary.
        """
        self.ratio_mix = ratio_mix or {}
        self.system = system
        self.connect_mix: Dict[float, str] = connect_mix or {}
    
    def add_mix(self, ratio_mix: Optional[Dict[str, int]] = None, system: Optional[Any] = None) -> Any:
        """
        Add or update the mix ratio and system.
        
        This method assigns crosslink types to models in the system according to the ratio mix,
        while preserving spatial relationships and connections between models.
        
        Args:
            ratio_mix: Optional dictionary of crosslink types and their ratios
            system: Optional system to update
            
        Returns:
            The updated system with crosslink types assigned
        """
        if ratio_mix:
            self.ratio_mix = ratio_mix
        if system:
            self.system = system

        if not self.system or not self.ratio_mix:
            raise ValueError("Both system and ratio_mix must be initialized")

        if isinstance(self.ratio_mix, str):
            ratio_dict = {}
            for part in self.ratio_mix.split():
                if ':' in part:
                    key, value = part.split(':')
                    try:
                        ratio_dict[key] = int(value)
                    except ValueError:
                        LOG.error(f"Invalid ratio value in {part}")
                        ratio_dict[key] = 0
            self.ratio_mix = ratio_dict
            LOG.debug(f"Converted ratio_mix from string: {self.ratio_mix}")

        total_models = self.system.get_size()
        
        model_ids = list(self.system.get_models())
        
        total_ratio = sum(self.ratio_mix.values())
        type_counts = {}
        remaining = total_models
        
        for type_name, ratio in sorted(self.ratio_mix.items(), key=lambda x: x[1], reverse=True):
            count = max(1, int((ratio / total_ratio) * total_models))
            type_counts[type_name] = count
            remaining -= count
        
        if remaining > 0:
            for type_name in sorted(type_counts.keys(), key=lambda x: self.ratio_mix[x], reverse=True):
                type_counts[type_name] += 1
                remaining -= 1
                if remaining == 0:
                    break
        
        connected_components = []
        visited = set()
        
        for model_id in model_ids:
            if model_id in visited:
                continue
                
            component = set()
            queue = [model_id]
            
            while queue:
                current_id = queue.pop(0)
                if current_id in component:
                    continue
                    
                component.add(current_id)
                visited.add(current_id)
                
                current_model = self.system.get_model(model_id=current_id)
                if hasattr(current_model, 'connect') and current_model.connect:
                    for connect_id in current_model.connect:
                        if connect_id not in component and connect_id in model_ids:
                            queue.append(connect_id)
            
            if component:
                connected_components.append(component)
        
        connected_components.sort(key=len, reverse=True)
        
        type_assignments = {t: 0 for t in type_counts.keys()}
        assigned_types = {}
        
        for component in connected_components:
            available_types = [t for t in type_counts.keys() if type_assignments[t] < type_counts[t]]
            
            if not available_types:
                chosen_type = max(self.ratio_mix.keys(), key=lambda t: self.ratio_mix[t])
            else:
                percentages = {
                    t: type_assignments[t] / type_counts[t]
                    for t in available_types
                }
                
                chosen_type = min(percentages.keys(), key=lambda t: percentages[t])
            
            for model_id in component:
                assigned_types[model_id] = chosen_type
            
            type_assignments[chosen_type] += len(component)
            
        # Apply the assigned types to each model
        for model_id in model_ids:
            if model_id in assigned_types:
                model = self.system.get_model(model_id=model_id)
                model_type = assigned_types[model_id]
                
                model.type = model_type
                model.crosslink_type = model_type
                
                LOG.debug(f"Assigned type {model.type} to model {model_id}")
                
                # Set the type for connected models too
                if hasattr(model, 'connect') and model.connect:
                    for connect_id in model.connect:
                        try:
                            connected_model = self.system.get_model(model_id=connect_id)
                            connected_model.type = model_type
                            connected_model.crosslink_type = model_type
                        except Exception as e:
                            LOG.warning(f"  - Could not set type for connected model {connect_id}: {e}")
            else:
                LOG.warning(f"Model {model_id} was not assigned a type")
        
        # Log the final type distribution
        actual_distribution = {}
        for model_id in model_ids:
            model = self.system.get_model(model_id=model_id)
            if hasattr(model, 'type'):
                actual_distribution[model.type] = actual_distribution.get(model.type, 0) + 1
        
        return self.system

    def get_mix(self, ratio_mix: Optional[List[int]] = None) -> str:
        """
        Get a random mix based on the provided or stored ratio mix.

        Args:
            ratio_mix (Optional[List[int]]): List of ratios to use for mixing.

        Returns:
            str: Randomly selected component based on the ratios.
        """
        if ratio_mix is None:
            ratio_mix = list(self.ratio_mix.values())
            
        return np.random.choice(list(self.ratio_mix.keys()), p=[i/100 for i in ratio_mix])

    def get_mix_from_connect_file(self, system: Optional[Any] = None, connect_file: Optional[str] = None) -> Any:
        """
        Update the system and assign types to models based on a connect file.

        Args:
            system (Optional[Any]): New system object to update.
            connect_file (Optional[str]): Path to the connect file.

        Returns:
            Any: The updated system object.
        """
        if system:
            self.system = system
        self.connect_mix = self.get_connect_mix(connect_file=connect_file)
        
        for model_id in self.system.get_models():
            model = self.system.get_model(model_id=model_id)
            if model.connect is not None and len(model.connect) > 1:
                for connect_id in model.connect:
                    self.system.get_model(model_id=connect_id).type = self.connect_mix[model_id]
        return self.system
    
    def get_connect_mix(self, connect_file: Optional[Union[str, Path]] = None) -> Dict[float, str]:
        """
        Read and parse the connect file to create a connect mix dictionary.

        Args:
            connect_file (Optional[Union[str, Path]]): Path to the connect file.

        Returns:
            Dict[float, str]: A dictionary mapping model IDs to their types.

        Raises:
            ValueError: If connect_file is not provided.
            FileNotFoundError: If the connect file is not found.
            IsADirectoryError: If the connect file path is a directory.
        """
        if connect_file is None:
            raise ValueError("connect_file must be provided")
        
        connect_file = Path(connect_file)
        
        if connect_file.exists():
            file_to_open = connect_file
        elif connect_file.with_suffix('.txt').exists():
            file_to_open = connect_file.with_suffix('.txt')
        else:
            file_to_open = connect_file
        
        try:
            with open(file_to_open, 'r') as f:
                return {float(l.split(';')[0].split(' ')[0].split('.')[0]): l.split(';')[1].strip() for l in f}
        except FileNotFoundError:
            raise FileNotFoundError(f"Connect file not found: {file_to_open}")
        except IsADirectoryError:
            raise IsADirectoryError(f"Connect file is a directory: {file_to_open}")

if __name__ == "__main__":
    pass