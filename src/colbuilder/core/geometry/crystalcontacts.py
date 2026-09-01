# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0

from typing import Dict, List, Optional, Union, Any
from pathlib import Path


class CrystalContacts:
    """
    Reads contact information from crystal_contact file and writes updated contact information for Chimera.

    Attributes:
        crystalcontacts_file (Path): Path to the crystal contacts file.
        t_matrix (Dict[float, List[float]]): Dictionary to store transformation matrices.
        models (List[float]): List of model IDs.
    """

    def __init__(self, crystalcontacts_file: Optional[Union[str, Path]] = None):
        self.crystalcontacts_file: Path = (
            Path(crystalcontacts_file) if crystalcontacts_file else Path()
        )
        self.t_matrix: Dict[float, List[float]] = {}
        self.models: List[float] = []

    def read_crystalcontacts(
        self, crystalcontacts_file: Optional[Union[str, Path]] = None
    ) -> List[str]:
        """
        Read crystal contacts information from Chimera contact output file.

        Args:
            crystalcontacts_file (Optional[Union[str, Path]]): Path to the crystal contacts file.

        Returns:
            List[str]: Lines from the crystal contacts file.
        """
        file_path = (
            Path(crystalcontacts_file)
            if crystalcontacts_file
            else self.crystalcontacts_file
        )
        with open(file_path.with_suffix(".txt"), "r") as f:
            return f.readlines()

    def read_t_matrix(
        self,
        crystalcontacts_file: Optional[Union[str, Path]] = None,
        crystalcontacts: Optional[List[str]] = None,
    ) -> Dict[float, List[float]]:
        """
        Read transformation matrix T from contact file.

        Args:
            crystalcontacts_file (Optional[Union[str, Path]]): Path to the crystal contacts file.
            crystalcontacts (Optional[List[str]]): List of lines from the crystal contacts file.

        Returns:
            Dict[float, List[float]]: Dictionary of transformation matrices.
        """
        if crystalcontacts is None:
            crystalcontacts = self.read_crystalcontacts(crystalcontacts_file)

        self.t_matrix.clear()
        for idx in range(0, len(crystalcontacts), 4):
            model_id = float(crystalcontacts[idx].split(" ")[1])
            self.t_matrix[model_id] = [
                float(crystalcontacts[idx + 1].split(" ")[-1]),
                float(crystalcontacts[idx + 2].split(" ")[-1]),
                float(crystalcontacts[idx + 3].split(" ")[-1]),
            ]
        return self.t_matrix

    def write_crystalcontacts(
        self,
        system: Optional[Any] = None,
        crystalcontacts_file: Optional[Union[str, Path]] = None,
    ) -> None:
        """
        Writes crystal contacts to txt file for Chimera.

        Args:
            system (Optional[Any]): System object containing model information.
            crystalcontacts_file (Optional[Union[str, Path]]): Path to the output crystal contacts file.
        """
        file_path = (
            Path(crystalcontacts_file)
            if crystalcontacts_file
            else self.crystalcontacts_file
        )

        with open(file_path.with_suffix(".txt"), "w") as f:
            if system is None:
                contacts = self.read_crystalcontacts(file_path)
                for key_cc, val_cc in self.t_matrix.items():
                    f.write(f"Model {key_cc}\n")
                    for i, val in enumerate(val_cc):
                        f.write(
                            f"         {'1' if i == 0 else '0'} {'1' if i == 1 else '0'} {'1' if i == 2 else '0'} {val:.3f}\n"
                        )
            else:
                for model in system.get_models():
                    f.write(f"Model {model}\n")
                    for i, val in enumerate(
                        system.get_model(model_id=model).transformation
                    ):
                        f.write(
                            f"         {'1' if i == 0 else '0'} {'1' if i == 1 else '0'} {'1' if i == 2 else '0'} {val:.3f}\n"
                        )
