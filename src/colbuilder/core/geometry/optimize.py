# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0

from __future__ import annotations
import numpy as np
from itertools import product
from typing import Dict, List, Optional, Any, Set, Tuple, Union
import os
from pathlib import Path

from colbuilder.core.geometry import model
from colbuilder.core.utils.logger import setup_logger

LOG = setup_logger(__name__)

class Optimizer:
    """
    Optimizer sets up a grid, depending on the z-position, and extends the current plane by
    adding nodes in integer space. Final check if added node is kept: Added node and 
    point-reflection node have to be connected in euclidean space. 
    As a result, the initial system is optimized.

    Attributes:
        system (Any): The system object to be optimized.
        solution_space (List[int]): The solution space dimensions [dx, dy, dz].
        t_matrix (Dict[float, List[float]]): Transformation matrix dictionary.
        s_matrix (Dict[float, List[int]]): Shift matrix dictionary.
        grid (List[List[float]]): The grid of nodes.
        paired_crosslinks (Set[str]): Set of crosslinks that have been paired.
    """

    def __init__(self, system: Any, solution_space: List[int]):
        self.system = system
        self.solution_space = solution_space
        self.t_matrix = system.crystalcontacts.read_t_matrix()
        self.s_matrix = {k: system.crystal.get_s_matrix(t_matrix=self.t_matrix[k]) for k in self.t_matrix}
        self.grid: List[List[float]] = []
        self.paired_crosslinks: Set[str] = set()

    def get_grid(self, z_grid: int, s_matrix: Optional[Dict[float, List[int]]] = None) -> np.ndarray:
        """
        Get the grid for a specific z position.

        Args:
            z_grid (int): The z-position of the grid.
            s_matrix (Optional[Dict[float, List[int]]]): The shift matrix to use.

        Returns:
            np.ndarray: The grid for the given z position.
        """
        s_matrix = s_matrix or self.s_matrix
        
        return np.array([
            [coords[0], coords[1], z_grid]
            for coords in s_matrix.values()
            if abs(coords[2] - z_grid) < 1e-6  # Compare with small tolerance
        ])

    def extend_grid(self, z_grid: int, s_matrix: Optional[Dict[float, List[int]]] = None) -> np.ndarray:
        """
        Update grid (x,y) at specific z-pos by adding nodes or filling up the grid.

        Args:
            z_grid (int): The z-position of the grid.
            s_matrix (Optional[Dict[float, List[int]]]): The shift matrix to use.

        Returns:
            np.ndarray: The extended grid.
        """
        s_matrix = s_matrix or self.s_matrix
        grid = self.get_grid(z_grid=z_grid, s_matrix=s_matrix)
        
        if grid.size == 0:
            x_min, x_max = -self.solution_space[0], self.solution_space[0]
            y_min, y_max = -self.solution_space[1], self.solution_space[1]
        else:
            x_max, x_min = np.max(grid[:, 0]), np.min(grid[:, 0])
            y_max, y_min = np.max(grid[:, 1]), np.min(grid[:, 1])
        
        d_x, d_y = int(self.solution_space[0]), int(self.solution_space[1])
        
        x_mesh = np.arange(x_min - d_x, x_max + d_x + 1)
        y_mesh = np.arange(y_min - d_y, y_max + d_y + 1)
        
        return np.array(np.meshgrid(x_mesh, y_mesh, [z_grid])).T.reshape(-1, 3)

    def get_nodes(self, z_grid: int, s_matrix: Optional[Dict[float, List[int]]] = None) -> List[List[float]]:
        """
        Compare nodes before and after grid extension step and update nodes in integer space.

        Args:
            z_grid (int): The z-position of the grid.
            s_matrix (Optional[Dict[float, List[int]]]): The shift matrix to use.

        Returns:
            List[List[float]]: List of new nodes.
        """
        s_matrix = s_matrix or self.s_matrix
        grid_init = set(map(tuple, self.get_grid(z_grid=z_grid, s_matrix=s_matrix)))
        grid_extend = set(map(tuple, self.extend_grid(z_grid=z_grid, s_matrix=s_matrix)))
        self.grid = [list(node) for node in grid_extend.difference(grid_init)]
        return self.grid

    def set_grid(self, z_grid: Optional[int] = None, s_matrix: Optional[Dict[float, List[int]]] = None) -> List[List[float]]:
        """
        Set up new meshgrid with more nodes to fill-up gaps at each z-pos.

        Args:
            z_grid (Optional[int]): The z-position of the grid.
            s_matrix (Optional[Dict[float, List[int]]]): The shift matrix to use.

        Returns:
            List[List[float]]: The new grid.
        """
        s_matrix = s_matrix or self.s_matrix
        if z_grid is None:
            z_grid = int(max(max(v) for v in s_matrix.values()))
        return self.get_nodes(z_grid, s_matrix)

    def run_optimize(self, s_matrix: Optional[Dict[float, List[int]]] = None, connect: Any = None, system: Optional[Any] = None) -> Any:
        """
        Wrapper to perform the grid-optimization in integer space with multiple passes.
        """
        s_matrix = s_matrix or self.s_matrix
        system = system or self.system
        
        LOG.debug("Starting initial optimization pass")
        system = self.optimize_crystalcontacts(s_matrix=s_matrix, connect=connect, system=system)
        
        if connect:
            system = self._update_system_connections(system, connect)
        
        LOG.debug("Starting targeted optimization for incomplete crosslinks")
        system = self.optimize_incomplete_crosslinks(s_matrix=s_matrix, connect=connect, system=system)
        
        if connect:
            system = self._update_system_connections(system, connect)
        
        LOG.debug("Starting final optimization pass")
        system = self.optimize_crystalcontacts(s_matrix=s_matrix, connect=connect, system=system)
        
        if connect:
            system = self._update_system_connections(system, connect)

        self._log_optimization_results(system)
        
        return system

    def _update_system_connections(self, system: Any, connect: Any) -> Any:
        """
        Update all connections in the system.
        
        Args:
            system (Any): The system to update.
            connect (Any): The connect object.
            
        Returns:
            Any: The updated system.
        """
        LOG.debug("     Updating system connections")
        contact_connect = connect.get_contact_connect(system=system)
        
        for model_id, connected_models in contact_connect.items():
            try:
                system.get_model(model_id=model_id).add_connect(
                    connect_id=model_id,
                    connect=connected_models
                )
            except Exception as e:
                LOG.warning(f"Failed to update connections for model {model_id}: {str(e)}")
                
        return system

    def _log_optimization_results(self, system: Any) -> None:
        """
        Log information about the optimized system.
        
        Args:
            system (Any): The optimized system.
        """
        total_models = len(list(system.get_models()))
        crosslinked_models = sum(1 for model_id in system.get_models() 
                                if hasattr(system.get_model(model_id=model_id), 'crosslink') and
                                system.get_model(model_id=model_id).crosslink)
        
        connected_models = sum(1 for model_id in system.get_models()
                             if hasattr(system.get_model(model_id=model_id), 'connect') and
                             system.get_model(model_id=model_id).connect)
        
        LOG.debug(f"Optimization completed: {total_models} total models, {crosslinked_models} with crosslinks, {connected_models} with connections")
        
    def check_node_connect(self, connect: Any, system: Any, node: List[float]) -> bool:
        """
        Check if the added node is connected to any other node of the grid.

        Args:
            connect (Any): The connect object.
            system (Any): The system object.
            node (List[float]): The node to check.

        Returns:
            bool: True if the node is connected, False otherwise.
        """
        return connect.run_connect(system=system, unit_cell=node)

    def optimize_crystalcontacts(self, s_matrix: Optional[Dict[float, List[int]]] = None, connect: Any = None, system: Optional[Any] = None) -> Any:
        """
        Algorithm adds models to the system if they are connected to any other model.
        Double check based on the point reflection of system: 
        Both added model and point reflection have to be connected to be added.

        Args:
            s_matrix (Optional[Dict[float, List[int]]]): The shift matrix to use.
            connect (Any): The connect object.
            system (Optional[Any]): The system to optimize.

        Returns:
            Any: The optimized system.
        """
        s_matrix = s_matrix or self.s_matrix
        system = system or self.system

        d_z = int(self.solution_space[2])
        z_grid = int(max(max(v) for v in s_matrix.values()))

        # Search in a wider z-range
        for plane in range(z_grid - d_z*2, z_grid + d_z*2 + 1):
            for node in self.set_grid(z_grid=plane, s_matrix=s_matrix):
                # Skip if we already have a model at this position
                if self._node_exists_in_system(node, system):
                    continue
                    
                if self.check_node_connect(connect=connect, system=system, node=node):
                    pr_node = [-i for i in node]
                    
                    # Skip if we already have a model at the reflection position
                    if self._node_exists_in_system(pr_node, system):
                        continue
                        
                    if self.check_node_connect(connect=connect, system=system, node=pr_node):
                        self._add_model_pair(system, node, pr_node)
                        
        return system

    def _node_exists_in_system(self, node: List[float], system: Any) -> bool:
        """
        Check if a node already exists in the system.
        
        Args:
            node (List[float]): The node coordinates to check.
            system (Any): The system to check against.
            
        Returns:
            bool: True if the node exists, False otherwise.
        """
        for model_id in system.get_models():
            model = system.get_model(model_id=model_id)
            if hasattr(model, 'unit_cell') and model.unit_cell:
                if np.allclose(node, model.unit_cell, rtol=1e-5, atol=1e-8):
                    return True
        return False

    def optimize_incomplete_crosslinks(self, s_matrix: Optional[Dict[float, List[int]]] = None, connect: Any = None, system: Optional[Any] = None) -> Any:
        """
        Special optimization pass focusing specifically on completing unpaired crosslinks.
        
        Args:
            s_matrix (Optional[Dict[float, List[int]]]): The shift matrix to use.
            connect (Any): The connect object.
            system (Optional[Any]): The system to optimize.
            
        Returns:
            Any: The optimized system.
        """
        s_matrix = s_matrix or self.s_matrix
        system = system or self.system
        
        incomplete_models = self._find_models_with_incomplete_crosslinks(system)
        LOG.debug(f"    Found {len(incomplete_models)} models with incomplete crosslinks")
        
        if not incomplete_models:
            return system
            
        models_fixed = 0
        for model_id in incomplete_models:
            model = system.get_model(model_id=model_id)
            
            unpaired_crosslinks = self._get_unpaired_crosslinks(model, system)
            for crosslink in unpaired_crosslinks:
                
                if self._create_complementary_model(crosslink, model, system, connect):
                    LOG.debug(f"        Successfully added complementary model for crosslink in model {model_id}")
                    models_fixed += 1
                    
        LOG.debug(f"    Fixed {models_fixed} crosslinks in targeted optimization pass")
        return system
        
    def _find_models_with_incomplete_crosslinks(self, system: Any) -> List[float]:
        """
        Find all models that have unpaired crosslinks.
        
        Args:
            system (Any): The system to analyze.
            
        Returns:
            List[float]: List of model IDs with unpaired crosslinks.
        """
        result = []
        
        for model_id in system.get_models():
            model = system.get_model(model_id=model_id)
            
            if not hasattr(model, 'crosslink') or not model.crosslink:
                continue
                
            if not hasattr(model, 'connect') or not model.connect:
                result.append(model_id)
                continue
            
            if self._get_unpaired_crosslinks(model, system):
                result.append(model_id)
                
        return result
        
    def _get_unpaired_crosslinks(self, model: Any, system: Any) -> List[Any]:
        """
        Identify which crosslinks in a model don't have matching pairs.
        
        Args:
            model (Any): The model to analyze.
            system (Any): The system containing all models.
            
        Returns:
            List[Any]: List of unpaired crosslinks.
        """
        if not hasattr(model, 'connect') or not model.connect:
            return model.crosslink
        
        unpaired = []
        for cross in model.crosslink:
            crosslink_id = f"{model.id}_{cross.resid}_{cross.type}"
            
            if crosslink_id in self.paired_crosslinks:
                continue
                
            has_pair = False
            
            for connect_id in model.connect:
                if connect_id == model.id:  # Skip self-connections
                    continue
                    
                try:
                    connected_model = system.get_model(model_id=connect_id)
                    
                    if not hasattr(connected_model, 'crosslink') or not connected_model.crosslink:
                        continue
                        
                    for connected_cross in connected_model.crosslink:
                        if np.linalg.norm(cross.position - connected_cross.position) < 3.0:
                            self.paired_crosslinks.add(crosslink_id)
                            self.paired_crosslinks.add(f"{connected_model.id}_{connected_cross.resid}_{connected_cross.type}")
                            has_pair = True
                            break
                            
                    if has_pair:
                        break
                        
                except Exception as e:
                    LOG.warning(f"Error checking connected model {connect_id}: {str(e)}")
                    continue
            
            if not has_pair:
                unpaired.append(cross)
                
        return unpaired
        
    def _create_complementary_model(self, crosslink: Any, source_model: Any, system: Any, connect: Any) -> bool:
        """
        Try to create a complementary model that would pair with the given crosslink,
        using a prioritized search based on crystal symmetry.
        
        Args:
            crosslink (Any): The crosslink needing a pair.
            source_model (Any): The model containing the crosslink.
            system (Any): The system to modify.
            connect (Any): The connect object.
            
        Returns:
            bool: True if successfully created a matching model, False otherwise.
        """
        source_pos = crosslink.position
        
        try:
            cs_matrix = system.crystal.read_cs_matrix()
            approx_unit_cell = system.crystal.get_s_matrix(cs_matrix=cs_matrix, t_matrix=source_pos.tolist())
            center_unit_cell = [round(x) for x in approx_unit_cell]
            
            priority_offsets = [
                # Primary nearest neighbors (face-adjacent)
                (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1),
                
                # Secondary nearest neighbors (edge-adjacent)
                (1, 1, 0), (1, -1, 0), (-1, 1, 0), (-1, -1, 0),
                (1, 0, 1), (1, 0, -1), (-1, 0, 1), (-1, 0, -1),
                (0, 1, 1), (0, 1, -1), (0, -1, 1), (0, -1, -1),
                
                # Tertiary nearest neighbors (corner-adjacent)
                (1, 1, 1), (1, 1, -1), (1, -1, 1), (1, -1, -1),
                (-1, 1, 1), (-1, 1, -1), (-1, -1, 1), (-1, -1, -1),
                
                # Further positions if needed
                (2, 0, 0), (-2, 0, 0), (0, 2, 0), (0, -2, 0), (0, 0, 2), (0, 0, -2),
            ]
            
            for offset in priority_offsets:
                candidate = [center_unit_cell[0] + offset[0], 
                            center_unit_cell[1] + offset[1], 
                            center_unit_cell[2] + offset[2]]
                
                if self._node_exists_in_system(candidate, system):
                    continue
                    
                transformation = system.crystal.get_t_matrix(s_matrix=candidate)
                
                temp_model = model.Model(
                    id=float(-1),
                    unit_cell=candidate,
                    transformation=transformation,
                    pdb_file=system.crystal.pdb_file
                )
                
                if connect.get_connect(ref_model=source_model, model=temp_model):
                    # LOG.debug(f"Found complementary position at unit cell {candidate}")
                    
                    pr_candidate = [-x for x in candidate]
                    pr_transformation = system.crystal.get_t_matrix(s_matrix=pr_candidate)
                    
                    new_id = float(system.get_size())
                    
                    new_model = model.Model(
                        id=new_id,
                        unit_cell=candidate,
                        transformation=transformation,
                        pdb_file=system.crystal.pdb_file
                    )
                    
                    system.add_model(model=new_model)
                    
                    if not self._node_exists_in_system(pr_candidate, system):
                        pr_model = model.Model(
                            id=new_id + 1,
                            unit_cell=pr_candidate,
                            transformation=pr_transformation,
                            pdb_file=system.crystal.pdb_file
                        )
                        system.add_model(model=pr_model)
                    
                    self._mark_model_crosslinks_paired(new_model)
                    
                    crosslink_id = f"{source_model.id}_{crosslink.resid}_{crosslink.type}"
                    self.paired_crosslinks.add(crosslink_id)
                    
                    return True
                    
            return False
                
        except Exception as e:
            LOG.warning(f"Error creating complementary model: {str(e)}")
            return False
    
    def _mark_model_crosslinks_paired(self, model: Any) -> None:
        """
        Mark all crosslinks in a model as paired.
        
        Args:
            model (Any): The model whose crosslinks should be marked as paired.
        """
        if not hasattr(model, 'crosslink') or not model.crosslink:
            return
            
        for cross in model.crosslink:
            crosslink_id = f"{model.id}_{cross.resid}_{cross.type}"
            self.paired_crosslinks.add(crosslink_id)
        
    def _add_model_pair(self, system: Any, node: List[float], pr_node: List[float]):
        """
        Add a pair of models to the system.

        Args:
            system (Any): The system to add models to.
            node (List[float]): The node coordinates.
            pr_node (List[float]): The point-reflected node coordinates.
        """
        for current_node in (node, pr_node):
            if self._node_exists_in_system(current_node, system):
                continue
                
            new_model = model.Model(
                id=float(system.get_size()),
                unit_cell=current_node,
                transformation=system.crystal.get_t_matrix(s_matrix=current_node),
                pdb_file=system.crystal.pdb_file
            )
            
            system.add_model(model=new_model)

            self._mark_model_crosslinks_paired(new_model)