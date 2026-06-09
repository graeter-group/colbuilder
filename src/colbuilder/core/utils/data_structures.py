"""
Data Structures for Crosslinking and Optimization in ColBuilder

This module defines structured data classes to represent crosslinking information and track
optimization states in the ColBuilder pipeline. These data structures provide a clear and
consistent way to manage residue positions, crosslink pairs, and optimization progress.

Key Features:
--------------
1. **Crosslink Representation**:
   - `CrosslinkPosition`: Represents a residue and atom involved in a crosslink, including chain ID,
     residue type, and atom name.
   - `CrosslinkPair`: Represents a crosslink between two or three positions, supporting both
     divalent and trivalent crosslinks.

2. **Optimization State Tracking**:
   - `OptimizationState`: Tracks the state of crosslink optimization, including current and best
     distances, atomic coordinates, and optimization history.

3. **Validation**:
   - Ensures valid chain IDs (`A`, `B`, `C`) and terminal types (`N`, `C`) for crosslinks.
   - Raises `ValueError` for invalid inputs.

Classes:
---------
1. **CrosslinkPosition**:
   - Attributes:
     - `residue_number` (int): Residue sequence number.
     - `chain_id` (str): Chain identifier (`A`, `B`, `C`).
     - `residue_type` (str): Three-letter residue code.
     - `atom_name` (str): Atom identifier in the residue.
   - Methods:
     - `position_str`: Returns the position in the format `number.chain`.

2. **CrosslinkPair**:
   - Attributes:
     - `position1` (CrosslinkPosition): First position in the crosslink.
     - `position2` (CrosslinkPosition): Second position in the crosslink.
     - `position3` (Optional[CrosslinkPosition]): Optional third position for trivalent crosslinks.
     - `terminal_type` (str): Terminal type (`N` or `C`).
   - Methods:
     - `is_trivalent`: Returns `True` if the crosslink is trivalent.

3. **OptimizationState**:
   - Attributes:
     - `current_distance` (float): Current distance between crosslinked atoms.
     - `attempt_number` (int): Current optimization attempt number.
     - `best_distance` (float): Best distance achieved so far.
     - `coordinates` (Optional[np.ndarray]): Current atomic coordinates.
     - `optimization_history` (List[Dict[str, Any]]): History of optimization attempts.
   - Methods:
     - `update`: Updates the optimization state with new distance and coordinates.
     - `increment_attempt`: Increments the optimization attempt counter.

Usage:
------
These data structures are used throughout the ColBuilder pipeline to manage crosslinking and
optimization processes.

Example:
--------
```python
from colbuilder.core.utils.data_structures import CrosslinkPosition, CrosslinkPair, OptimizationState

# Define crosslink positions
pos1 = CrosslinkPosition(residue_number=10, chain_id="A", residue_type="LYS", atom_name="NZ")
pos2 = CrosslinkPosition(residue_number=25, chain_id="B", residue_type="ASP", atom_name="OD1")

# Create a crosslink pair
crosslink = CrosslinkPair(position1=pos1, position2=pos2, terminal_type="N")
print(crosslink.is_trivalent)  # Output: False

# Track optimization state
state = OptimizationState()
state.update(distance=5.0, coords=np.array([[1.0, 2.0, 3.0]]))
state.increment_attempt()
print(state.best_distance)  # Output: 5.0
```
"""

from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any
import numpy as np


@dataclass(frozen=True)
class CrosslinkPosition:
    """
    Represents a position in a crosslink, including residue and atom information.

    Attributes:
        residue_number (int): The residue sequence number
        chain_id (str): Chain identifier (A, B, C)
        residue_type (str): Three-letter residue code
        atom_name (str): Atom identifier in the residue

    Raises:
        ValueError: If chain_id is not A, B, or C
    """

    residue_number: int
    chain_id: str
    residue_type: str
    atom_name: str
    atom_name2: Optional[str] = None  # second atom (CSV A32) for the third residue of a trivalent crosslink

    def __post_init__(self):
        if self.chain_id not in {"A", "B", "C"}:
            raise ValueError(f"Chain ID must be A, B, or C, got {self.chain_id}")

    @property
    def position_str(self) -> str:
        """Returns the position in format 'number.chain'"""
        return f"{self.residue_number}.{self.chain_id}"


@dataclass(frozen=True)
class CrosslinkPair:
    """
    Represents a crosslink between two or three positions.

    Attributes:
        position1 (CrosslinkPosition): First position in crosslink
        position2 (CrosslinkPosition): Second position in crosslink
        position3 (Optional[CrosslinkPosition]): Optional third position for trivalent crosslinks
        terminal_type (str): Terminal type (N or C)
    """

    position1: CrosslinkPosition
    position2: CrosslinkPosition
    position3: Optional[CrosslinkPosition] = None
    terminal_type: str = field(default="N")

    def __post_init__(self):
        if self.terminal_type not in {"N", "C"}:
            raise ValueError(f"Terminal type must be N or C, got {self.terminal_type}")

    @property
    def is_trivalent(self) -> bool:
        """Returns True if this is trivalent crosslink."""
        return self.position3 is not None


@dataclass
class OptimizationState:
    """
    Tracks the state of crosslink optimization.

    Attributes:
        current_distance: Current distance between crosslinked atoms
        attempt_number: Current optimization attempt number
        best_distance: Best distance achieved so far
        coordinates: Numpy array of current atomic coordinates
    """

    current_distance: float = float("inf")
    attempt_number: int = 0
    best_distance: float = float("inf")
    coordinates: Optional[np.ndarray] = None
    optimization_history: List[Dict[str, Any]] = field(default_factory=list)

    def update(self, distance: float, coords: Optional[np.ndarray] = None) -> None:
        """
        Update optimization state with new values.

        Args:
            distance: Current distance between crosslinked atoms
            coords: Optional coordinates array. If None, coordinates won't be updated.
        """
        self.current_distance = distance
        if distance < self.best_distance:
            self.best_distance = distance
            if coords is not None:
                self.coordinates = coords.copy()
        self.optimization_history.append(
            {"attempt": self.attempt_number, "distance": distance}
        )

    def increment_attempt(self) -> None:
        """Increment the attempt counter."""
        self.attempt_number += 1
