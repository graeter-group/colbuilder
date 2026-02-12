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
    
    def _build_components(self, model_ids: List[float]) -> List[List[float]]:
        """
        Build connected components from the system connectivity.
        
        Models with no connectivity information are treated as single-node components.
        Uses depth-first search to find all connected models.
        
        Args:
            model_ids: List of all model IDs in the system
            
        Returns:
            List of connected components, where each component is a list of model IDs.
            Components are sorted by size (largest first), then by first model ID.
        """
        components: List[List[float]] = []
        visited = set()

        for model_id in model_ids:
            if model_id in visited:
                continue
            
            # Use stack for depth-first search
            stack = [model_id]
            comp: List[float] = []
            
            while stack:
                current = stack.pop()
                if current in visited:
                    continue
                    
                visited.add(current)
                comp.append(current)
                
                # Get neighbors from connectivity
                try:
                    neighbors = self.system.get_model(model_id=current).connect or []
                except Exception:
                    neighbors = []
                
                for neigh in neighbors:
                    if neigh not in visited:
                        stack.append(neigh)
            
            components.append(sorted(comp))
        
        # Sort components: largest first, then by first model ID
        components.sort(key=lambda c: (-len(c), c[0]))
        return components

    def _compute_targets(self, total_models: int) -> Dict[str, int]:
        """
        Compute deterministic model counts per type from the requested ratios.
        
        Distributes models according to the ratio_mix proportions. Fractional
        models are distributed to types with the largest remainders to ensure
        the total equals total_models.
        
        Args:
            total_models: Total number of models to distribute
            
        Returns:
            Dictionary mapping each type to its target count of models
        """
        total_ratio = float(sum(self.ratio_mix.values())) or 1.0
        
        # Calculate raw (fractional) counts for each type
        raw = {
            t: (self.ratio_mix[t] / total_ratio) * total_models 
            for t in self.ratio_mix
        }

        # Floor all values to get integer base counts
        targets = {t: int(np.floor(v)) for t, v in raw.items()}
        
        # Distribute remaining models based on largest fractional remainders
        remainder = total_models - sum(targets.values())
        if remainder > 0:
            # Sort by fractional part (descending), then by ratio (descending)
            fractions = sorted(
                ((t, raw[t] - targets[t]) for t in self.ratio_mix),
                key=lambda x: (-x[1], -self.ratio_mix[x[0]]),
            )
            for i in range(remainder):
                t, _ = fractions[i % len(fractions)]
                targets[t] += 1
                
        return targets

    def add_mix(self, ratio_mix: Optional[Dict[str, int]] = None, system: Optional[Any] = None) -> Any:
        """
        Add or update the mix ratio and system.
        
        This method assigns crosslink types to models in the system according to the ratio mix,
        while preserving spatial relationships and connections between models. Connected
        components are kept together and assigned the same type.
        
        The algorithm:
        1. Computes target counts for each type based on ratios
        2. Identifies connected components in the system
        3. Assigns types to components using round-robin to avoid starving small ratios
        4. Ensures all models in a component get the same type
        
        Args:
            ratio_mix: Optional dictionary of crosslink types and their ratios
            system: Optional system to update
            
        Returns:
            The updated system with crosslink types assigned
            
        Raises:
            ValueError: If both system and ratio_mix are not initialized
        """
        if ratio_mix:
            self.ratio_mix = ratio_mix
        if system:
            self.system = system

        if not self.system or not self.ratio_mix:
            raise ValueError("Both system and ratio_mix must be initialized")

        # Handle string format for ratio_mix (e.g., "D:80 T:20")
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

        # Compute target counts and build connected components
        target_counts = self._compute_targets(total_models)
        components = self._build_components(model_ids)

        # Deterministic assignment: assign components to types using round-robin
        # This ensures small ratios aren't starved by early assignment of large components
        remaining = target_counts.copy()
        fallback_type = max(self.ratio_mix, key=self.ratio_mix.get)
        assigned_types: Dict[float, str] = {}

        # Create ordered list of types for round-robin assignment
        type_order = sorted(self.ratio_mix, key=self.ratio_mix.get, reverse=True)
        type_idx = 0

        for comp in components:
            # Round-robin across types with remaining quota
            chosen = None
            for _ in range(len(type_order)):
                candidate = type_order[type_idx % len(type_order)]
                type_idx += 1
                if remaining.get(candidate, 0) > 0:
                    chosen = candidate
                    break

            type_choice = chosen or fallback_type

            # Assign the chosen type to all models in this component
            for mid in comp:
                assigned_types[mid] = type_choice

            # Decrement the remaining count for this type
            remaining[type_choice] = remaining.get(type_choice, 0) - len(comp)

        # Apply the assigned types to each model
        for model_id in model_ids:
            model = self.system.get_model(model_id=model_id)
            model_type = assigned_types.get(model_id, fallback_type)
            
            # Set both type and crosslink_type attributes
            model.type = model_type
            model.crosslink_type = model_type
            
            LOG.debug(f"Assigned type {model.type} to model {model_id}")
            
            # Also set type for connected models to ensure consistency
            if hasattr(model, 'connect') and model.connect:
                for connect_id in model.connect:
                    try:
                        connected_model = self.system.get_model(model_id=connect_id)
                        connected_model.type = model_type
                        connected_model.crosslink_type = model_type
                    except Exception as e:
                        LOG.warning(f"Could not set type for connected model {connect_id}: {e}")

        # Log the final type distribution
        actual_distribution: Dict[str, int] = {}
        for model_id in model_ids:
            model = self.system.get_model(model_id=model_id)
            if hasattr(model, 'type'):
                actual_distribution[model.type] = actual_distribution.get(model.type, 0) + 1

        LOG.info(f"Final type distribution: {actual_distribution} (targets: {target_counts})")

        return self.system

    def get_mix(self, ratio_mix: Optional[List[int]] = None) -> str:
        """
        Get a random mix based on the provided or stored ratio mix.

        Uses weighted random selection according to the ratio values.

        Args:
            ratio_mix (Optional[List[int]]): List of ratios to use for mixing.

        Returns:
            str: Randomly selected component type based on the ratios.
        """
        if ratio_mix is None:
            ratio_mix = list(self.ratio_mix.values())
            
        return np.random.choice(list(self.ratio_mix.keys()), p=[i/100 for i in ratio_mix])

    def get_mix_from_connect_file(self, system: Optional[Any] = None, connect_file: Optional[str] = None) -> Any:
        """
        Update the system and assign types to models based on a connect file.

        The connect file should contain lines in the format:
        "model_id.caps.pdb ; TYPE"
        
        This method reads the file and assigns the specified type to each model
        and its connected neighbors.

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

        The connect file format is:
        "model_id.caps.pdb ; TYPE"
        
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
        
        # Try to find the file with or without .txt extension
        if connect_file.exists():
            file_to_open = connect_file
        elif connect_file.with_suffix('.txt').exists():
            file_to_open = connect_file.with_suffix('.txt')
        else:
            file_to_open = connect_file
        
        try:
            with open(file_to_open, 'r') as f:
                return {
                    float(l.split(';')[0].split(' ')[0].split('.')[0]): 
                    l.split(';')[1].strip() 
                    for l in f
                }
        except FileNotFoundError:
            raise FileNotFoundError(f"Connect file not found: {file_to_open}")
        except IsADirectoryError:
            raise IsADirectoryError(f"Connect file is a directory: {file_to_open}")

if __name__ == "__main__":
    pass
