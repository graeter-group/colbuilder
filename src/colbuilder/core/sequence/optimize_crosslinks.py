"""
This module provides tools for optimizing collagen structures to satisfy crosslinking constraints.

It implements methods to model and optimize crosslinked collagen structures by applying geometric
transformations to residues and minimizing the distances between crosslinked atoms.

Key Features:
--------------
1. **Crosslink Distance Calculation**:
   - Calculate Euclidean distances between crosslinked residues or atoms.
   - Support for both divalent (two residues) and trivalent (three residues) crosslinks.

2. **Geometric Transformations**:
   - Rotate residues around backbone torsion angles (phi and psi).
   - Rotate side chains around chi1 and chi2 angles.
   - Apply general rotation matrices to side chains.
   - Rotate side chains relative to the backbone plane.

3. **Monte Carlo Optimization**:
   - Use Monte Carlo methods to iteratively optimize crosslink distances.
   - Simulated annealing is employed to escape local minima and improve convergence.

4. **Transformation Tracking**:
   - Track all transformations applied to residues for reproducibility.
   - Apply transformations to residues in the original structure after optimization.

5. **PDB File Handling**:
   - Load and save PDB structures.
   - Extract atomic coordinates and manipulate residues.

6. **Crosslink Matching**:
   - Identify potential residue matches for crosslinking based on residue type and position.
   - Select the best matching crosslinks based on initial distances.

Usage:
------
This module is designed to be used as part of a pipeline for optimizing collagen structures with
crosslinking constraints. The main entry point is the `optimize_structure` function, which takes
PDB files and crosslink specifications as input and outputs an optimized structure.

Example:
--------
```python
from pathlib import Path

# Define input files and crosslink information
initial_pdb = "input_structure.pdb"
copy1_pdb = "translated_copy.pdb"
copy2_pdb = "original_copy.pdb"
optimized_pdb = "optimized_structure.pdb"
crosslink_info = [
    {
        "residue1_type": "L5Y",
        "residue1_position": "9",
        "atom1": "NZ",
        "residue2_type": "L4Y",
        "residue2_position": "947",
        "atom2": "CE",
        "residue3_type": "NONE",
        "residue3_position": "",
        "atom31": "",
        "atom32": ""
    }
]

# Optimize the structure
total_distance, tracker = optimize_structure(
    initial_pdb=initial_pdb,
    copy1_pdb=copy1_pdb,
    copy2_pdb=copy2_pdb,
    crosslink_info=crosslink_info,
    optimized_pdb=optimized_pdb
)

print(f"Total crosslink distance after optimization: {total_distance:.2f}")
```
"""

# Copyright (c) 2024, ColBuilder Development Team
# Distributed under the terms of the Apache License 2.0

import numpy as np
from scipy.spatial.transform import Rotation
from Bio import PDB
import random
import warnings
import math
from collections import defaultdict
from pathlib import Path
from typing import Dict, Tuple, Any, Optional, Union, List
import numpy.typing as npt
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure
import subprocess

from colbuilder.core.utils.logger import setup_logger

LOG = setup_logger(__name__)

warnings.filterwarnings("ignore", message="Ignoring unrecognized record.*")


def load_pdb(filename: Union[str, Path]) -> PDB.Structure.Structure:
    """Load PDB structure from file."""
    parser = PDB.PDBParser()
    return parser.get_structure("molecule", filename)


def save_pdb(structure: PDB.Structure.Structure, filename: Union[str, Path]) -> None:
    """Save structure to PDB file."""
    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(filename)


def distance(coord1: npt.NDArray[np.float64], coord2: npt.NDArray[np.float64]) -> float:
    """Calculate Euclidean distance between two coordinates."""
    return np.linalg.norm(coord1 - coord2)


def get_atom_coords(
    structure: PDB.Structure.Structure, chain_id: str, residue_id: str, atom_name: str
) -> npt.NDArray[np.float64]:
    """
    Get atomic coordinates from structure.

    Args:
        structure: PDB structure
        chain_id: Chain identifier
        residue_id: Residue number
        atom_name: Atom name

    Returns:
        Numpy array of coordinates

    Raises:
        ValueError: If atom not found in structure
    """
    chain = structure[0][chain_id]
    for residue in chain:
        if residue.id[1] == int(residue_id):
            return residue[atom_name].coord
    raise ValueError(
        f"Atom {atom_name} in residue {residue_id} not found in chain {chain_id}"
    )


def get_distances(
    structures: Dict[str, PDB.Structure.Structure], crosslink: Dict[str, Dict[str, Any]]
) -> Tuple[float, float]:
    """
    Calculate distances between crosslinked residues.

    Args:
        structures: Dictionary of PDB structures
        crosslink: Crosslink specification dictionary

    Returns:
        Tuple of distances between crosslinked atoms
    """
    r1_coord = get_atom_coords(
        structures[crosslink["R1"]["structure_id"]],
        crosslink["R1"]["chain"],
        crosslink["R1"]["position"],
        crosslink["R1"]["atom"],
    )
    r2_coord = get_atom_coords(
        structures[crosslink["R2"]["structure_id"]],
        crosslink["R2"]["chain"],
        crosslink["R2"]["position"],
        crosslink["R2"]["atom"],
    )

    if crosslink["R3"]["type"] == "NONE":
        dist = distance(r1_coord, r2_coord)
        return dist, dist

    r3_coord1 = get_atom_coords(
        structures[crosslink["R3"]["structure_id"]],
        crosslink["R3"]["chain"],
        crosslink["R3"]["position"],
        crosslink["R3"]["atom31"],
    )
    r3_coord2 = get_atom_coords(
        structures[crosslink["R3"]["structure_id"]],
        crosslink["R3"]["chain"],
        crosslink["R3"]["position"],
        crosslink["R3"]["atom32"],
    )
    return distance(r1_coord, r3_coord1), distance(r2_coord, r3_coord2)


def get_backbone_plane(residue: Residue) -> npt.NDArray[np.float64]:
    """Get backbone plane normal vector from residue N, CA, C atoms."""
    n_coord = residue["N"].coord
    ca_coord = residue["CA"].coord
    c_coord = residue["C"].coord

    v1 = ca_coord - n_coord
    v2 = c_coord - ca_coord
    normal = np.cross(v1, v2)
    return normal / np.linalg.norm(normal)


def get_phi_psi_atoms(
    residue: Residue,
) -> Optional[Dict[str, Union[Dict[str, npt.NDArray[np.float64]], Tuple[int, int]]]]:
    """
    Get atoms involved in phi/psi angles and valid residue range.

    Args:
        residue: Residue to get atoms for

    Returns:
        Dictionary with phi/psi atoms and valid residue range, or None if terminal residue
    """
    try:
        chain = residue.get_parent()
        res_id = residue.id[1]

        if res_id <= 1 or res_id >= len(chain) - 1:
            return None

        atoms = {
            "phi": {
                "prev_C": chain[res_id - 1]["C"].coord,
                "N": residue["N"].coord,
                "CA": residue["CA"].coord,
                "C": residue["C"].coord,
            },
            "psi": {
                "N": residue["N"].coord,
                "CA": residue["CA"].coord,
                "C": residue["C"].coord,
                "next_N": chain[res_id + 1]["N"].coord,
            },
            "range": (res_id, min(len(chain) - 1, res_id + 0)),
        }
        return atoms
    except Exception:
        return None


def get_chi1_axis(residue: Residue) -> npt.NDArray[np.float64]:
    """Get normalized axis vector for chi1 rotation."""
    ca_coord = residue["CA"].coord
    cb_coord = residue["CB"].coord
    return (cb_coord - ca_coord) / np.linalg.norm(cb_coord - ca_coord)


def get_chi2_axis(residue: Residue) -> Optional[npt.NDArray[np.float64]]:
    """Get normalized axis vector for chi2 rotation if CG exists."""
    if "CG" in residue:
        cb_coord = residue["CB"].coord
        cg_coord = residue["CG"].coord
        return (cg_coord - cb_coord) / np.linalg.norm(cg_coord - cb_coord)
    return None


def log_crosslink_info(crosslink: Dict[str, Dict[str, str]], index: int) -> None:
    """Log detailed crosslink information."""
    LOG.debug(f"\nCrosslink {index + 1} Details:")
    LOG.debug("=" * 50)

    for residue in ["R1", "R2", "R3"]:
        LOG.debug(f"\n{residue} information:")
        LOG.debug("-" * 20)
        res_info = crosslink[residue]
        LOG.debug(f"Structure ID: {res_info['structure_id']}")
        LOG.debug(f"Chain: {res_info['chain']}")
        LOG.debug(f"Position: {res_info['position']}")
        LOG.debug(f"Type: {res_info['type']}")
        if residue in ["R1", "R2"]:
            LOG.debug(f"Atom: {res_info['atom']}")
        else:
            LOG.debug(f"Atom 1: {res_info['atom31']}")
            LOG.debug(f"Atom 2: {res_info['atom32']}")


class TransformationTracker:
    """Tracks and applies geometric transformations to protein structure residues."""

    def __init__(self):
        self.transformations: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
        self.initial_coords: Dict[str, Dict[str, npt.NDArray[np.float64]]] = {}
        self.source_mapping: Dict[str, str] = {}

    def copy(self) -> "TransformationTracker":
        """Create a deep copy of the tracker."""
        new_tracker = TransformationTracker()
        new_tracker.transformations = defaultdict(list)
        for key, transforms in self.transformations.items():
            new_tracker.transformations[key] = [dict(t) for t in transforms]

        new_tracker.initial_coords = {
            key: {atom: coord.copy() for atom, coord in coords.items()}
            for key, coords in self.initial_coords.items()
        }

        new_tracker.source_mapping = dict(self.source_mapping)
        return new_tracker

    def update_from(self, other_tracker: "TransformationTracker") -> None:
        """Update this tracker with transformations from another tracker."""
        for key, transforms in other_tracker.transformations.items():
            self.transformations[key].extend([dict(t) for t in transforms])

        for key, coords in other_tracker.initial_coords.items():
            if key not in self.initial_coords:
                self.initial_coords[key] = {
                    atom: coord.copy() for atom, coord in coords.items()
                }

        self.source_mapping.update(other_tracker.source_mapping)

    def add_transformation(
        self,
        structure_id: str,
        chain_id: str,
        residue_id: str,
        transform_type: str,
        params: Dict[str, Any],
    ) -> None:
        """
        Record a transformation with source tracking.

        Args:
            structure_id: ID of structure being transformed
            chain_id: Chain identifier
            residue_id: Residue number
            transform_type: Type of transformation
            params: Transformation parameters
        """
        key = f"{chain_id}_{residue_id}"
        self.source_mapping[key] = structure_id
        transform = {
            "type": transform_type,
            "params": params.copy(),
            "step": len(self.transformations[key]),
            "source": structure_id,
        }
        self.transformations[key].append(transform)

    def apply_transformations_to_residue(
        self,
        structure: Structure,
        chain_id: str,
        residue_id: str,
        source_structure_id: Optional[str] = None,
    ):
        """Apply only transformations from the specified source structure."""
        key = f"{chain_id}_{residue_id}"
        if key not in self.transformations:
            return

        # Filter transformations by source if specified
        transforms = self.transformations[key]
        if source_structure_id:
            transforms = [t for t in transforms if t["source"] == source_structure_id]

        sorted_transforms = sorted(transforms, key=lambda x: x["step"])
        for transform in sorted_transforms:
            params = transform["params"]

            if transform["type"] == "side_chain":
                if "chi_type" in params:
                    if params["chi_type"] == "chi1":
                        self._apply_chi1_rotation(
                            structure, chain_id, int(residue_id), params["angle"]
                        )
                    elif params["chi_type"] == "chi2":
                        self._apply_chi2_rotation(
                            structure, chain_id, int(residue_id), params["angle"]
                        )
                else:
                    self._apply_rotation_matrix(
                        structure, chain_id, int(residue_id), params["rotation_matrix"]
                    )

            elif transform["type"] == "backbone":
                rotate_backbone(
                    structure,
                    chain_id,
                    int(residue_id),
                    params["angle"],
                    params["angle_type"],
                )

            elif transform["type"] == "relative_backbone":
                self._apply_relative_backbone(
                    structure, chain_id, int(residue_id), params["angle"]
                )

    def _apply_chi1_rotation(
        self, structure: Structure, chain_id: str, residue_id: int, angle: float
    ) -> None:
        """Apply chi1 rotation to a residue."""
        chain = structure[0][chain_id]
        for residue in chain:
            if residue.id[1] == residue_id:
                backbone_atoms = set(["N", "CA", "C", "O"])
                side_chain_atoms = [
                    atom for atom in residue if atom.name not in backbone_atoms
                ]
                if not side_chain_atoms:
                    return

                axis = get_chi1_axis(residue)
                for atom in side_chain_atoms:
                    if atom.name != "CB":
                        atom.coord = rotate_around_axis(
                            atom.coord, axis, angle, residue["CB"].coord
                        )

    def _apply_chi2_rotation(
        self, structure: Structure, chain_id: str, residue_id: int, angle: float
    ) -> None:
        """Apply chi2 rotation to a residue."""
        chain = structure[0][chain_id]
        for residue in chain:
            if residue.id[1] == residue_id:
                if "CG" not in residue:
                    return

                backbone_atoms = set(["N", "CA", "C", "O"])
                side_chain_atoms = [
                    atom for atom in residue if atom.name not in backbone_atoms
                ]

                axis = get_chi2_axis(residue)
                for atom in side_chain_atoms:
                    if atom.name not in ["CB", "CG"]:
                        atom.coord = rotate_around_axis(
                            atom.coord, axis, angle, residue["CG"].coord
                        )

    def _apply_rotation_matrix(
        self,
        structure: Structure,
        chain_id: str,
        residue_id: int,
        rotation_matrix: npt.NDArray[np.float64],
    ) -> None:
        """Apply general rotation matrix to residue side chain."""
        chain = structure[0][chain_id]
        for residue in chain:
            if residue.id[1] == residue_id:
                backbone_atoms = set(["N", "CA", "C", "O"])
                side_chain_atoms = [
                    atom for atom in residue if atom.name not in backbone_atoms
                ]

                center = residue["CB"].coord if "CB" in residue else residue["CA"].coord
                for atom in side_chain_atoms:
                    atom.coord = np.dot(rotation_matrix, atom.coord - center) + center

    def _apply_relative_backbone(
        self, structure: Structure, chain_id: str, residue_id: int, angle: float
    ) -> None:
        """Apply rotation relative to backbone plane."""
        chain = structure[0][chain_id]
        for residue in chain:
            if residue.id[1] == residue_id:
                backbone_atoms = set(["N", "CA", "C", "O"])
                side_chain_atoms = [
                    atom for atom in residue if atom.name not in backbone_atoms
                ]
                if not side_chain_atoms:
                    return

                backbone_normal = get_backbone_plane(residue)
                ca_coord = residue["CA"].coord
                cb_coord = residue["CB"].coord if "CB" in residue else None

                if cb_coord is not None:
                    ca_cb = cb_coord - ca_coord
                    rotation_axis = np.cross(backbone_normal, ca_cb)
                    rotation_axis = rotation_axis / np.linalg.norm(rotation_axis)
                    rotation = Rotation.from_rotvec(angle * rotation_axis)

                    for atom in side_chain_atoms:
                        atom.coord = ca_coord + rotation.apply(atom.coord - ca_coord)


def rotate_backbone(
    structure: Structure,
    chain_id: str,
    residue_id: int,
    angle: float,
    angle_type: str = "phi",
    tracker: Optional[TransformationTracker] = None,
    structure_id: Optional[str] = None,
) -> None:
    """
    Rotate around phi or psi angle with optional transformation tracking.

    Args:
        structure: PDB structure to modify
        chain_id: Chain identifier
        residue_id: Residue number to rotate
        angle: Rotation angle in radians
        angle_type: Type of rotation ('phi' or 'psi')
        tracker: Optional transformation tracker
        structure_id: Optional structure identifier for tracking
    """
    chain = structure[0][chain_id]
    for residue in chain:
        if residue.id[1] == int(residue_id):
            if residue.id[1] <= 1 or residue.id[1] >= len(chain) - 1:
                return

            atoms = get_phi_psi_atoms(residue)
            if not atoms:
                return

            if tracker and structure_id:
                tracker.add_transformation(
                    structure_id,
                    chain_id,
                    residue_id,
                    "backbone",
                    {"angle": angle, "angle_type": angle_type},
                )

            if angle_type == "phi":
                rotation_axis = atoms["phi"]["CA"] - atoms["phi"]["N"]
                rotation_axis = rotation_axis / np.linalg.norm(rotation_axis)
                center = atoms["phi"]["N"]

                for atom in residue:
                    if atom.name != "N":
                        atom.coord = rotate_around_axis(
                            atom.coord, rotation_axis, angle, center
                        )

                start_id, end_id = atoms["range"]
                for res_idx in range(residue.id[1] + 1, end_id + 1):
                    if res_idx in chain:
                        for atom in chain[res_idx]:
                            atom.coord = rotate_around_axis(
                                atom.coord, rotation_axis, angle, center
                            )

            elif angle_type == "psi":
                rotation_axis = atoms["psi"]["C"] - atoms["psi"]["CA"]
                rotation_axis = rotation_axis / np.linalg.norm(rotation_axis)
                center = atoms["psi"]["CA"]

                for atom in residue:
                    if atom.name not in ["N", "CA"]:
                        atom.coord = rotate_around_axis(
                            atom.coord, rotation_axis, angle, center
                        )

                start_id, end_id = atoms["range"]
                for res_idx in range(residue.id[1] + 1, end_id + 1):
                    if res_idx in chain:
                        for atom in chain[res_idx]:
                            atom.coord = rotate_around_axis(
                                atom.coord, rotation_axis, angle, center
                            )


def rotate_relative_to_backbone(
    structure: Structure,
    chain_id: str,
    residue_id: int,
    angle: float,
    tracker: Optional[TransformationTracker] = None,
    structure_id: Optional[str] = None,
) -> None:
    """
    Rotate side chain relative to backbone plane with optional transformation tracking.

    Args:
        structure: PDB structure to modify
        chain_id: Chain identifier
        residue_id: Residue number to rotate
        angle: Rotation angle in radians
        tracker: Optional transformation tracker
        structure_id: Optional structure identifier for tracking
    """
    chain = structure[0][chain_id]
    for residue in chain:
        if residue.id[1] == int(residue_id):
            backbone_atoms = set(["N", "CA", "C", "O"])
            side_chain_atoms = [
                atom for atom in residue if atom.name not in backbone_atoms
            ]
            if not side_chain_atoms:
                return

            backbone_normal = get_backbone_plane(residue)
            ca_coord = residue["CA"].coord
            cb_coord = residue["CB"].coord if "CB" in residue else None

            if cb_coord is not None:
                if tracker and structure_id:
                    tracker.add_transformation(
                        structure_id,
                        chain_id,
                        residue_id,
                        "relative_backbone",
                        {"angle": angle},
                    )

                ca_cb = cb_coord - ca_coord
                rotation_axis = np.cross(backbone_normal, ca_cb)
                rotation_axis = rotation_axis / np.linalg.norm(rotation_axis)
                rotation = Rotation.from_rotvec(angle * rotation_axis)

                for atom in side_chain_atoms:
                    atom.coord = ca_coord + rotation.apply(atom.coord - ca_coord)


def rotate_side_chain(
    structure: Structure,
    chain_id: str,
    residue_id: int,
    rotation_matrix: npt.NDArray[np.float64],
    tracker: Optional[TransformationTracker] = None,
    structure_id: Optional[str] = None,
) -> None:
    """
    Enhanced side chain rotation with transformation tracking.

    Args:
        structure: PDB structure to modify
        chain_id: Chain identifier
        residue_id: Residue number to rotate
        rotation_matrix: 3x3 rotation matrix
        tracker: Optional transformation tracker
        structure_id: Optional structure identifier for tracking
    """
    chain = structure[0][chain_id]
    for residue in chain:
        if residue.id[1] == int(residue_id):
            backbone_atoms = set(["N", "CA", "C", "O"])
            side_chain_atoms = [
                atom for atom in residue if atom.name not in backbone_atoms
            ]
            if not side_chain_atoms:
                return

            rotation_type = np.random.choice(
                ["chi1", "chi2", "random"], p=[0.6, 0.3, 0.1]
            )

            if rotation_type == "chi1":
                axis = get_chi1_axis(residue)
                center = residue["CA"].coord
                angle = np.random.uniform(-np.pi, np.pi)

                if tracker and structure_id:
                    tracker.add_transformation(
                        structure_id,
                        chain_id,
                        residue_id,
                        "side_chain",
                        {"chi_type": "chi1", "angle": angle},
                    )

                for atom in side_chain_atoms:
                    if atom.name != "CB":
                        atom.coord = rotate_around_axis(
                            atom.coord, axis, angle, residue["CB"].coord
                        )

            elif rotation_type == "chi2" and "CG" in residue:
                axis = get_chi2_axis(residue)
                if axis is not None:
                    center = residue["CB"].coord
                    angle = np.random.normal(0, np.pi / 6)

                    if tracker and structure_id:
                        tracker.add_transformation(
                            structure_id,
                            chain_id,
                            residue_id,
                            "side_chain",
                            {"chi_type": "chi2", "angle": angle},
                        )

                    for atom in side_chain_atoms:
                        if atom.name not in ["CB", "CG"]:
                            atom.coord = rotate_around_axis(
                                atom.coord, axis, angle, residue["CG"].coord
                            )
            else:
                center = residue["CB"].coord if "CB" in residue else residue["CA"].coord
                angle = np.random.normal(0, 0.1)
                axis = np.random.rand(3)
                axis /= np.linalg.norm(axis)
                rot = Rotation.from_rotvec(angle * axis)

                if tracker and structure_id:
                    tracker.add_transformation(
                        structure_id,
                        chain_id,
                        residue_id,
                        "side_chain",
                        {"rotation_matrix": rot.as_matrix()},
                    )

                for atom in side_chain_atoms:
                    atom.coord = np.dot(rot.as_matrix(), atom.coord - center) + center


def rotate_around_axis(
    point: npt.NDArray[np.float64],
    axis: npt.NDArray[np.float64],
    angle: float,
    center: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]:
    """Apply Rodrigues rotation formula ref:
    https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula."""
    k = axis / np.linalg.norm(axis)
    cos_t = np.cos(angle)
    sin_t = np.sin(angle)

    p = point - center
    p_rot = p * cos_t + np.cross(k, p) * sin_t + k * np.dot(k, p) * (1 - cos_t)

    return p_rot + center


def store_residue_coords(
    structure: Structure, chain_id: str, residue_id: int
) -> Dict[str, npt.NDArray[np.float64]]:
    """Store coordinates for a single residue."""
    chain = structure[0][chain_id]
    for residue in chain:
        if residue.id[1] == int(residue_id):
            return {atom.name: atom.coord.copy() for atom in residue}
    return {}


def restore_residue_coords(
    structure: Structure,
    chain_id: str,
    residue_id: int,
    coords: Dict[str, npt.NDArray[np.float64]],
) -> None:
    """Restore coordinates for a single residue."""
    chain = structure[0][chain_id]
    for residue in chain:
        if residue.id[1] == int(residue_id):
            for atom in residue:
                atom.coord = coords[atom.name]


def find_potential_matches(
    structures: Dict[str, Structure],
    residue1_info: Dict[str, str],
    residue2_info: Dict[str, str],
    residue3_info: Dict[str, str],
    is_divalent: bool = False,
    swap_structures: bool = False,
) -> List[Tuple[Dict[str, Any], float]]:
    """
    Find all potential matching residue triplets between structures and calculate relevant
    distances.

    Args:
        structures: Dictionary of all structures (copy1, copy2)
        residue1_info: information for first residue
        residue2_info: information for second residue
        residue3_info: information for third residue
        is_divalent: Whether this is a divalent crosslink
        swap_structures: Whether to swap structure assignments between copy1/copy2
    """
    matches = []

    def find_matching_residues(
        structure: Structure, res_type: str, position: str
    ) -> List[Tuple[str, int]]:
        """Find matching residues and return (chain_id, residue_id)"""
        matching = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.resname == res_type and str(residue.id[1]) == position:
                        matching.append((chain.id, residue.id[1]))
        return matching

    struct1 = "copy2" if swap_structures else "copy1"
    struct2 = "copy1" if swap_structures else "copy2"

    if is_divalent:
        LOG.debug(f"Divalent crosslink: R1 in {struct1}, R2 in {struct2}")
        r1_matches = find_matching_residues(
            structures[struct1], residue1_info["type"], residue1_info["position"]
        )
        r2_matches = find_matching_residues(
            structures[struct2], residue2_info["type"], residue2_info["position"]
        )
        for r1 in r1_matches:
            chain1, pos1 = r1
            for r2 in r2_matches:
                chain2, pos2 = r2

                crosslink = {
                    "R1": {
                        "structure_id": struct1,
                        "chain": chain1,
                        "position": str(pos1),
                        "type": residue1_info["type"],
                        "atom": residue1_info["atom1"],
                    },
                    "R2": {
                        "structure_id": struct2,
                        "chain": chain2,
                        "position": str(pos2),
                        "type": residue2_info["type"],
                        "atom": residue2_info["atom2"],
                    },
                    "R3": {
                        "structure_id": struct2,
                        "chain": "A",
                        "position": "1",
                        "type": "NONE",
                        "atom31": "",
                        "atom32": "",
                    },
                }

                try:
                    dist1, _ = get_distances(structures, crosslink)
                    LOG.debug(
                        f"Distance for R1({chain1}:{pos1})-R2({chain2}:{pos2}): {dist1:.2f}"
                    )
                    matches.append((crosslink, dist1))
                except (KeyError, ValueError) as e:
                    LOG.debug(f"Skipping invalid combination: {e}")
                    continue
    else:
        LOG.debug(f"Trivalent crosslink: R1/R2 in {struct1}, R3 in {struct2}")
        r1_matches = find_matching_residues(
            structures[struct1], residue1_info["type"], residue1_info["position"]
        )
        r2_matches = find_matching_residues(
            structures[struct1],  # Same structure as R1
            residue2_info["type"],
            residue2_info["position"],
        )
        r3_matches = find_matching_residues(
            structures[struct2],  # Different structure
            residue3_info["type"],
            residue3_info["position"],
        )
        for r1 in r1_matches:
            chain1, pos1 = r1
            for r2 in r2_matches:
                chain2, pos2 = r2
                for r3 in r3_matches:
                    chain3, pos3 = r3

                    crosslink = {
                        "R1": {
                            "structure_id": struct1,
                            "chain": chain1,
                            "position": str(pos1),
                            "type": residue1_info["type"],
                            "atom": residue1_info["atom1"],
                        },
                        "R2": {
                            "structure_id": struct1,
                            "chain": chain2,
                            "position": str(pos2),
                            "type": residue2_info["type"],
                            "atom": residue2_info["atom2"],
                        },
                        "R3": {
                            "structure_id": struct2,
                            "chain": chain3,
                            "position": str(pos3),
                            "type": residue3_info["type"],
                            "atom31": residue3_info["atom31"],
                            "atom32": residue3_info["atom32"],
                        },
                    }

                    try:
                        dist1, dist2 = get_distances(structures, crosslink)
                        LOG.debug(
                            f"Distances for R1({chain1}:{pos1})-R3({chain3}:{pos3}): {dist1:.2f}, "
                            f"R2({chain2}:{pos2})-R3({chain3}:{pos3}): {dist2:.2f}"
                        )
                        matches.append((crosslink, dist1))
                    except (KeyError, ValueError) as e:
                        LOG.debug(f"Skipping invalid combination: {e}")
                        continue

    matches.sort(key=lambda x: x[1])
    return matches


def select_best_matching_crosslinks(
    structures: Dict[str, Structure], crosslink_info: List[Dict[str, Any]]
) -> List[Dict[str, Dict[str, Any]]]:
    """Select the best matching crosslink pairs based on initial distances."""
    selected_matches = []

    for i, info in enumerate(crosslink_info):
        if not info:
            continue

        LOG.debug(f"\nProcessing crosslink {i+1}")
        LOG.debug("=" * 50)

        is_divalent = info.get("residue3_type") == "NONE"
        all_matches = []

        matches1 = find_potential_matches(
            structures,
            {
                "type": info["residue1_type"],
                "position": info["residue1_position"],
                "atom1": info["atom1"],
            },
            {
                "type": info["residue2_type"],
                "position": info["residue2_position"],
                "atom2": info["atom2"],
            },
            {
                "type": info["residue3_type"],
                "position": info["residue3_position"],
                "atom31": info["atom31"],
                "atom32": info["atom32"],
            },
            is_divalent=is_divalent,
        )
        all_matches.extend(matches1)

        matches2 = find_potential_matches(
            structures,
            {
                "type": info["residue1_type"],
                "position": info["residue1_position"],
                "atom1": info["atom1"],
            },
            {
                "type": info["residue2_type"],
                "position": info["residue2_position"],
                "atom2": info["atom2"],
            },
            {
                "type": info["residue3_type"],
                "position": info["residue3_position"],
                "atom31": info["atom31"],
                "atom32": info["atom32"],
            },
            is_divalent=is_divalent,
            swap_structures=True,
        )
        all_matches.extend(matches2)

        if all_matches:
            all_matches.sort(key=lambda x: x[1])
            best_match = all_matches[0]
            LOG.debug("\nFinal selection:")
            LOG.debug(f"Found {len(all_matches)} total possible matches")
            LOG.debug(f"Selected best match with distance {best_match[1]:.2f}")
            LOG.debug(
                f"R1: Chain {best_match[0]['R1']['chain']}, Pos {best_match[0]['R1']['position']}"
            )
            LOG.debug(
                f"R2: Chain {best_match[0]['R2']['chain']}, Pos {best_match[0]['R2']['position']}"
            )
            if not is_divalent:
                LOG.debug(
                    f"R3: Chain {best_match[0]['R3']['chain']}, Pos {best_match[0]['R3']['position']}"
                )
            selected_matches.append(best_match[0])
        else:
            LOG.warning("No valid matches found for crosslink specification")

    return selected_matches


def optimize_crosslink(
    structures: Dict[str, Structure],
    crosslink: Dict[str, Dict[str, Any]],
    tracker: TransformationTracker,
    max_steps: int = 20000,
    target_distance: float = 1.5,
    previous_best_distance: float = float("inf"),
) -> Tuple[Dict[str, Structure], TransformationTracker]:
    """
    Optimize crosslink geometry using Monte Carlo optimization.

    Args:
        structures: Dictionary of PDB structures
        crosslink: Crosslink specification dictionary
        tracker: Transformation tracker for recording moves
        max_steps: Maximum optimization steps
        target_distance: Target distance for optimization

    Returns:
        Tuple of optimized structures and transformation tracker
    """
    best_distance = previous_best_distance
    best_structures = {k: v.copy() for k, v in structures.items()}
    best_tracker = TransformationTracker()
    current_tracker = TransformationTracker()

    is_divalent = crosslink["R3"]["type"] == "NONE"
    residue_types = ["R1", "R2"] if is_divalent else ["R1", "R2", "R3"]
    target_distance = 1.5 if is_divalent else 0.1

    # Phase 1: Backbone exploration
    LOG.debug("\nPhase 1: Backbone exploration")
    angle_steps = np.linspace(-np.pi / 2, np.pi / 2, 8)
    for residue_type in residue_types:
        restype_best_tracker = TransformationTracker()
        for angle in angle_steps:
            structures_copy = {k: v.copy() for k, v in best_structures.items()}
            temp_tracker = TransformationTracker()
            residue = crosslink[residue_type]

            rotate_relative_to_backbone(
                structures_copy[residue["structure_id"]],
                residue["chain"],
                residue["position"],
                angle,
                tracker=temp_tracker,
                structure_id=residue["structure_id"],
            )

            dist1, dist2 = get_distances(structures_copy, crosslink)
            current_distance = dist1 if is_divalent else (dist1 + dist2)

            if current_distance < best_distance:
                best_distance = current_distance
                best_structures = {k: v.copy() for k, v in structures_copy.items()}
                restype_best_tracker = temp_tracker.copy()
                current_tracker.update_from(temp_tracker)
                if is_divalent:
                    LOG.debug(f"Improved distance: {dist1:.2f}")
                else:
                    LOG.debug(f"Improved: {dist1:.2f}, {dist2:.2f}")
            best_tracker = current_tracker.copy()

    # Phase 2: Main optimization
    LOG.debug("\nPhase 2: Main optimization")
    structures = {k: v.copy() for k, v in best_structures.items()}
    temperature = 1.0
    cooling_rate = 0.999
    min_temp = 0.2
    no_improvement_count = 0
    best_max_distance = float("inf")

    for step in range(max_steps):
        dist1, dist2 = get_distances(structures, crosslink)
        current_max_distance = dist1 if is_divalent else max(dist1, dist2)

        if is_divalent and dist1 <= target_distance + 1.0:
            LOG.debug(f"Target reached at step {step}")
            return best_structures, best_tracker

        if (
            not is_divalent
            and dist1 <= target_distance + 1.0
            and dist2 <= target_distance + 1.0
        ):
            LOG.debug(f"Target reached at step {step}")
            return best_structures, best_tracker

        if no_improvement_count > 500:
            LOG.debug(f"Resetting at step {step}")
            temperature = 1.0
            if random.random() < 0.5:
                structures = {k: v.copy() for k, v in best_structures.items()}
                current_tracker = best_tracker.copy()
            no_improvement_count = 0

        residue_type = (
            random.choice(["R1", "R2"])
            if is_divalent
            else (
                random.choice(["R1", "R3", "R3"])
                if dist1 > dist2
                else random.choice(["R2", "R3", "R3"])
            )
        )

        residue = crosslink[residue_type]
        structure = structures[residue["structure_id"]]

        old_coords = store_residue_coords(
            structure, residue["chain"], residue["position"]
        )
        temp_tracker = TransformationTracker()

        # Apply transformations
        if random.random() < 0.3 and temperature > 0.3:
            angle_type = random.choice(["phi", "psi"])
            angle = np.random.normal(0, 0.15 * temperature)
            rotate_backbone(
                structure,
                residue["chain"],
                residue["position"],
                angle,
                angle_type,
                tracker=temp_tracker,
                structure_id=residue["structure_id"],
            )
        else:
            if random.random() < 0.7:
                angle = np.random.normal(0, 0.3 * temperature)
                rotate_relative_to_backbone(
                    structure,
                    residue["chain"],
                    residue["position"],
                    angle,
                    tracker=temp_tracker,
                    structure_id=residue["structure_id"],
                )
            else:
                rotation = Rotation.random()
                rotate_side_chain(
                    structure,
                    residue["chain"],
                    residue["position"],
                    rotation.as_matrix(),
                    tracker=temp_tracker,
                    structure_id=residue["structure_id"],
                )

        new_dist1, new_dist2 = get_distances(structures, crosslink)
        new_max_distance = new_dist1 if is_divalent else max(new_dist1, new_dist2)

        # Acceptance criteria
        accept = False
        if is_divalent:
            if new_dist1 <= max(target_distance + 1.0, dist1):
                accept = True
            elif random.random() < math.exp(-(new_dist1 - dist1) / temperature):
                accept = True
        else:
            if min(new_dist1, new_dist2) < 1.5:
                accept = False
            elif new_dist1 <= max(target_distance + 1.0, dist1) and new_dist2 <= max(
                target_distance + 1.0, dist2
            ):
                accept = True
            elif new_max_distance < current_max_distance - 0.2:
                accept = True
            elif random.random() < math.exp(
                -(new_max_distance - current_max_distance) / temperature
            ):
                accept = True

        if accept:
            current_tracker.update_from(temp_tracker)
            if new_max_distance < best_max_distance:
                best_max_distance = new_max_distance
                best_structures = {k: v.copy() for k, v in structures.items()}
                best_tracker = current_tracker.copy()
                no_improvement_count = 0
            else:
                no_improvement_count += 1
        else:
            restore_residue_coords(
                structure, residue["chain"], residue["position"], old_coords
            )
            no_improvement_count += 1

        temperature = max(temperature * cooling_rate, min_temp)

        if step % 100 == 0:
            if is_divalent:
                LOG.debug(
                    f"Step {step}: Distance = {new_dist1:.2f}, T = {temperature:.4f}"
                )
            else:
                LOG.debug(
                    f"Step {step}: {new_dist1:.2f}, {new_dist2:.2f}, T = {temperature:.4f}"
                )

    return best_structures, best_tracker


def optimize_structure(
    initial_pdb: str,
    copy1_pdb: str,
    copy2_pdb: str,
    crosslink_info: List[Dict[str, Any]],
    optimized_pdb: str,
    previous_best_distance: float = float("inf"),
) -> Tuple[float, TransformationTracker, Path]:
    """
    Optimize protein structure to satisfy crosslinking constraints.

    Args:
        initial_pdb: Path to initial PDB file
        copy1_pdb: Path to first copy PDB - translated copy
        copy2_pdb: Path to second copy PDB - original copy
        crosslink_info: List of crosslink specifications
        optimized_pdb: Path to save optimized structure

    Returns:
        Tuple of total crosslink distance and transformation tracker
    """
    structures = {
        "initial": load_pdb(initial_pdb),
        "copy1": load_pdb(copy1_pdb),
        "copy2": load_pdb(copy2_pdb),
    }
    save_pdb(structures["copy1"], "initial_copy1.pdb")
    save_pdb(structures["copy2"], "initial_copy2.pdb")
    crosslinks = select_best_matching_crosslinks(structures, crosslink_info)
    master_tracker = TransformationTracker()

    # Calculate initial total distance before optimization
    initial_total_distance = 0
    for crosslink in crosslinks:
        dist1, dist2 = get_distances(structures, crosslink)
        initial_total_distance += dist1

    LOG.debug(
        f"Initial total distance before optimization: {initial_total_distance:.2f}"
    )

    for i, crosslink in enumerate(crosslinks):
        LOG.debug(f"\nOptimizing crosslink {i+1}")
        log_crosslink_info(crosslink, i)
        tracker = TransformationTracker()

        # If this is the first crosslink and we have a previous best, pass it
        crosslink_previous_best = previous_best_distance if i == 0 else float("inf")

        structures, crosslink_tracker = optimize_crosslink(
            structures,
            crosslink,
            tracker,
            previous_best_distance=crosslink_previous_best,  # Pass the parameter
        )
        master_tracker.update_from(crosslink_tracker)

    optimized_initial = load_pdb(initial_pdb)

    # Apply transformations to initial structure
    for crosslink in crosslinks:
        for res_type, structure_id in [
            ("R1", crosslink["R1"]["structure_id"]),
            ("R2", crosslink["R2"]["structure_id"]),
            ("R3", crosslink["R3"]["structure_id"]),
        ]:
            residue = crosslink[res_type]
            chain_id = residue["chain"]
            res_id = residue["position"]

            # Only apply transformations that came from the correct source structure
            master_tracker.apply_transformations_to_residue(
                optimized_initial,
                chain_id,
                res_id,
                source_structure_id=structure_id,  # Filter by source
            )

    save_pdb(optimized_initial, str(optimized_pdb))
    save_pdb(structures["copy1"], "optimized_copy1.pdb")
    save_pdb(structures["copy2"], "optimized_copy2.pdb")

    total_distance = 0
    for crosslink in crosslinks:
        dist1, dist2 = get_distances(structures, crosslink)
        total_distance += dist1
        if crosslink["R3"]["type"] == "NONE":
            LOG.debug(
                f"Final distance: {crosslink['R1']['type']}-{crosslink['R2']['type']} = {dist1:.2f}"
            )
        else:
            LOG.debug(
                f"Final distances: {crosslink['R1']['type']}-{crosslink['R3']['type']} = {dist1:.2f}, "
                f"{crosslink['R2']['type']}-{crosslink['R3']['type']} = {dist2:.2f}"
            )

    return total_distance, master_tracker, Path(optimized_pdb)
