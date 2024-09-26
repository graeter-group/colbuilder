# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0

from __future__ import annotations
import numpy as np
from typing import Dict, List, Optional, Any
import os

from colbuilder.core.geometry import model

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
    """

    def __init__(self, system: Any, solution_space: List[int]):
        self.system = system
        self.solution_space = solution_space
        self.t_matrix = system.crystalcontacts.read_t_matrix()
        self.s_matrix = {k: system.crystal.get_s_matrix(t_matrix=self.t_matrix[k]) for k in self.t_matrix}
        self.grid: List[List[float]] = []

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
        Wrapper to perform the grid-optimization in integer space.

        Args:
            s_matrix (Optional[Dict[float, List[int]]]): The shift matrix to use.
            connect (Any): The connect object.
            system (Optional[Any]): The system to optimize.

        Returns:
            Any: The optimized system.
        """
        s_matrix = s_matrix or self.s_matrix
        system = system or self.system
        return self.optimize_crystalcontacts(s_matrix=s_matrix, connect=connect, system=system)

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

        for plane in range(z_grid - d_z, z_grid + 1):
            for node in self.set_grid(z_grid=plane, s_matrix=s_matrix):
                if self.check_node_connect(connect=connect, system=system, node=node):
                    pr_node = [-i for i in node]
                    if self.check_node_connect(connect=connect, system=system, node=pr_node):
                        self._add_model_pair(system, node, pr_node)
        return system

    def _add_model_pair(self, system: Any, node: List[float], pr_node: List[float]):
        """
        Add a pair of models to the system.

        Args:
            system (Any): The system to add models to.
            node (List[float]): The node coordinates.
            pr_node (List[float]): The point-reflected node coordinates.
        """
        for current_node in (node, pr_node):
            system.add_model(model.Model(
                id=float(system.get_size()),
                unit_cell=current_node,
                transformation=system.crystal.get_t_matrix(s_matrix=current_node),
                pdb_file=system.crystal.pdb_file
            ))
