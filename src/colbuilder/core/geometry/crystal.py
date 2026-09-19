# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0

from __future__ import annotations
import numpy as np
from typing import Dict, List, Optional, Tuple
from pathlib import Path

from colbuilder.core.utils.logger import setup_logger

LOG = setup_logger(__name__)


class Crystal:
    """
    Reads crystal information to generate crystal-symmetry matrix.
    Allows conversion between unit-cell shift & transformation matrix.

    Attributes
    ----------
    pdb_file : Path
        Path to the PDB file.
    is_line : Tuple[str, ...]
        Tuple of valid line types in PDB files.
    crystal : Dict[str, Optional[float]]
        Dictionary to store crystal parameters.
    """

    def __init__(self, pdb=None):
        self.pdb_file = Path(pdb) if pdb else None
        self.is_line = ("ATOM  ", "HETATM", "ANISOU")
        self.crystal = {k: None for k in ["a", "b", "c", "alpha", "beta", "gamma"]}

    def read_crystal(self, pdb: Optional[Path] = None) -> Dict[str, str]:
        """
        Read crystallographic information from PDB file.

        Parameters
        ----------
        pdb : Optional[Path], default=None
            Path to the PDB file. If None, uses the instance's pdb_file.

        Returns
        -------
        Dict[str, str]
            Dictionary of crystal parameters.
        """
        pdb = pdb or self.pdb_file
        with open(pdb.with_suffix(".pdb"), "r") as f:
            for line in f:
                if line.startswith("CRYST1"):
                    return dict(zip(self.crystal.keys(), line.split()[1:]))
        raise ValueError(f"No CRYST1 record found in {pdb}")

    def read_spacegroup(self, pdb: Optional[Path] = None) -> int:
        """
        Read space-group from crystallographic information in PDB file.

        Parameters
        ----------
        pdb : Optional[Path], default=None
            Path to the PDB file. If None, uses the instance's pdb_file.

        Returns
        -------
        int
            Space group number.
        """
        pdb = pdb or self.pdb_file
        with open(pdb.with_suffix(".pdb"), "r") as f:
            for line in f:
                if line.startswith("CRYST1"):
                    return int(line.split()[-2])
        raise ValueError(f"No CRYST1 record found in {pdb}")

    def get_default_transformation(self) -> List[float]:
        """
        Get a default transformation matrix for the crystal.

        Returns
        -------
        List[float]
            Default transformation matrix (identity matrix).

        Raises
        ------
        ValueError
            If the crystal symmetry matrix cannot be determined.
        """
        try:
            cs_matrix = self.read_cs_matrix()
            # Default transformation is the origin (0, 0, 0) in the unit cell
            default_shift = [0, 0, 0]
            return self.get_t_matrix(cs_matrix=cs_matrix, s_matrix=default_shift)
        except Exception as e:
            LOG.error(f"Failed to generate default transformation matrix: {str(e)}")
            raise ValueError("Could not generate default transformation matrix.") from e

    def read_cs_matrix(
        self,
        pdb: Optional[Path] = None,
        spacegroup: Optional[int] = None,
        crystal: Optional[Dict[str, str]] = None,
    ) -> np.ndarray:
        """
        Determine crystal-symmetry matrix CS based on crystallographic information & space group.

        Parameters
        ----------
        pdb : Optional[Path], default=None
            Path to the PDB file.
        spacegroup : Optional[int], default=None
            Space group number.
        crystal : Optional[Dict[str, str]], default=None
            Dictionary of crystal parameters.

        Returns
        -------
        np.ndarray
            Crystal-symmetry matrix.

        Raises
        ------
        ValueError
            If space group is not recognized.
        """
        if spacegroup is None:
            spacegroup = self.read_spacegroup(pdb)
        if crystal is None:
            crystal = self.read_crystal(pdb)

        if spacegroup == 1:
            ax, ay, az = float(crystal["a"]), 0, 0
            bx = float(crystal["b"]) * np.cos(np.deg2rad(float(crystal["gamma"])))
            by = float(crystal["b"]) * np.sin(np.deg2rad(float(crystal["gamma"])))
            bz = 0
            cx = float(crystal["c"]) * np.cos(np.deg2rad(float(crystal["beta"])))
            cy = float(crystal["c"]) * (
                (
                    np.cos(np.deg2rad(float(crystal["alpha"])))
                    - np.cos(np.deg2rad(float(crystal["beta"])))
                    * np.cos(np.deg2rad(float(crystal["gamma"])))
                )
                / np.sin(np.deg2rad(float(crystal["gamma"])))
            )
            cz = np.sqrt(
                np.power(float(crystal["c"]), 2) - np.power(cx, 2) - np.power(cy, 2)
            )
            return np.array([[ax, bx, cx], [ay, by, cy], [az, bz, cz]]).round(
                decimals=4
            )
        else:
            raise ValueError(
                f"Space-group not recognized: Crystal-rotation-matrix only for space-group 1 available. Got {spacegroup}"
            )

    def get_s_matrix(
        self,
        pdb: Optional[Path] = None,
        cs_matrix: Optional[np.ndarray] = None,
        t_matrix: Optional[List[float]] = None,
    ) -> List[int]:
        """
        Get shift matrix S from transformation matrix T using crystal-symmetry matrix CS: S = T \ CS.

        Parameters
        ----------
        pdb : Optional[Path], default=None
            Path to the PDB file.
        cs_matrix : Optional[np.ndarray], default=None
            Crystal-symmetry matrix.
        t_matrix : Optional[List[float]], default=None
            Transformation matrix.

        Returns
        -------
        List[int]
            Shift matrix.

        Raises
        ------
        ValueError
            If no transform matrix is provided.
        """
        if cs_matrix is None:
            cs_matrix = self.read_cs_matrix(pdb)
        if t_matrix is None:
            raise ValueError(
                "No transform-matrix given, hence no unit-cell shift matrix can be calculated."
            )
        return list(np.linalg.solve(cs_matrix, t_matrix).round(decimals=0).astype(int))

    def get_t_matrix(
        self,
        pdb: Optional[Path] = None,
        cs_matrix: Optional[np.ndarray] = None,
        s_matrix: Optional[List[int]] = None,
    ) -> List[float]:
        """
        Get transformation matrix T from shift matrix S using crystal-symmetry matrix CS: T = CS x S.

        Parameters
        ----------
        pdb : Optional[Path], default=None
            Path to the PDB file.
        cs_matrix : Optional[np.ndarray], default=None
            Crystal-symmetry matrix.
        s_matrix : Optional[List[int]], default=None
            Shift matrix.

        Returns
        -------
        List[float]
            Transformation matrix.

        Raises
        ------
        ValueError
            If no unit-cell shift matrix is provided.
        """
        if cs_matrix is None:
            cs_matrix = self.read_cs_matrix(pdb)
        if s_matrix is None:
            raise ValueError(
                "No unit-cell shift-matrix given, hence no transform matrix can be calculated."
            )
        return list(np.dot(cs_matrix, s_matrix).round(decimals=3).astype(float))

    def translate_crystal(self, pdb=None, translate=None, bool_system=False):
        """
        Translate center of gravity of unit cell to position [0, 0, 400].

        Parameters
        ----------
        pdb : Union[str, Path], default=None
            Path to the PDB file.
        translate : Optional[List[float]], default=None
            Translation vector.
        bool_system : bool, default=False
            Whether to translate the entire system.
        """

        pdb = Path(pdb) if pdb else self.pdb_file
        if not bool_system:
            translate = [0, 0, translate[2] - self.get_cog(pdb)]
        with open(pdb.with_suffix(".pdb"), "r") as f:
            atoms = f.readlines()
        with open(pdb.with_suffix(".pdb"), "w") as f:
            for line in atoms:
                if line[:6] in self.is_line:
                    x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                    new_z = round(z + translate[2], 3)
                    f.write(f"{line[:30]}{x:8.3f}{y:8.3f}{new_z:8.3f}{line[54:]}")
                else:
                    f.write(line)

    def get_cog(self, pdb=None):
        """
        Get center of gravity of PDB file.

        Parameters
        ----------
        pdb : Union[str, Path], default=None
            Path to the PDB file.

        Returns
        -------
        float
            Z-coordinate of the center of gravity.
        """
        pdb = Path(pdb) if pdb else self.pdb_file
        with open(pdb.with_suffix(".pdb"), "r") as f:
            z_coords = [float(line[46:54]) for line in f if line[:6] in self.is_line]
        return np.nanmean(z_coords)
