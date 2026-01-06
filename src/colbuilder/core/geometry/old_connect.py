# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0

from __future__ import annotations
import numpy as np
from itertools import product
from typing import Dict, List, Optional, Any, Union, Set, Tuple
from pathlib import Path

from colbuilder.core.geometry import model
from colbuilder.core.utils.logger import setup_logger

LOG = setup_logger(__name__)


class Connect:
    """
    Get connections between all models in system OR between potentially added model and system.

    Attributes:
        system (Any): The system object containing models.
        pairs (Dict[float, Optional[float]]): Dictionary of model pairs.
        connect (Dict[float, List[float]]): Dictionary of connections.
        connect_file (Optional[Path]): Path to the connection file.
        external_connect (List[float]): List of external connections.
        is_line (Tuple[str, ...]): Tuple of valid line types in PDB files.
    """

    def __init__(
        self, system: Optional[Any] = None, connect_file: Optional[Path] = None
    ):
        self.system = system
        self.pairs: Dict[float, Optional[float]] = {
            key: None for key in self.system.get_models()
        }
        self.connect: Dict[float, List[float]] = {}
        self.connect_file = Path(connect_file) if connect_file else None
        self.external_connect: List[float] = []
        self.is_line: Tuple[str, ...] = ("ATOM  ", "HETATM", "ANISOU", "TER   ")

    def get_model_connect(self, system: Any, unit_cell: List[float]) -> bool:
        """
        Get connection between added model and all already existing models in system.

        Args:
            system (Any): The system object.
            unit_cell (List[float]): Unit cell parameters.

        Returns:
            bool: True if a connection is found, False otherwise.
        """
        transformation = system.crystal.get_t_matrix(s_matrix=unit_cell)
        add_ = model.Model(
            id="add", transformation=transformation, pdb_file=system.crystal.pdb_file
        )

        for ref_model in self.system.get_models():
            if self.get_connect(
                ref_model=system.get_model(model_id=ref_model), model=add_
            ):
                return True
        return False

    def get_contact_connect(self, system: Any) -> Dict[float, List[float]]:
        """
        Get connection between all models/contacts in system.

        Args:
            system (Any): The system object.

        Returns:
            Dict[float, List[float]]: Dictionary of connections.
        """
        self.pairs = {key: [] for key in self.system.get_models()}
        for ref_model, model in product(self.system.get_models(), repeat=2):
            if ref_model != model:
                if self.get_connect(
                    ref_model=system.get_model(model_id=ref_model),
                    model=system.get_model(model_id=model),
                ):
                    self.pairs[ref_model].append(model)

        return self.merge_contacts(pairs=self.pairs)

    def merge_contacts(
        self, pairs: Dict[float, List[float]]
    ) -> Dict[float, List[float]]:
        """
        Merges contacts to generate groups of connections
        """
        self.connect = {key: set([key]) for key in pairs}
        for ref_key, connected_models in pairs.items():
            for model in connected_models:
                self.connect[ref_key].add(model)
                if model in self.connect:
                    self.connect[model].add(ref_key)

        self.connect = {k: sorted(v) for k, v in self.connect.items()}

        return self.clean_contacts(contactpairs=self.connect)

    def get_external_connect_file(
        self, system: Any, connect_file: Optional[Path] = None
    ) -> Any:
        """
        Read external connect file and update system accordingly.

        Args:
            system (Any): The system object.
            connect_file (Optional[Path]): Path to the connection file.

        Returns:
            Any: Updated system object.
        """
        if connect_file:
            self.external_connect = [
                float(l.split(" ")[0].replace(".caps.pdb", ""))
                for l in open(connect_file.with_suffix(".txt"), "r")
            ]
        if np.min(self.external_connect) > 0:
            self.external_connect = [i - 1 for i in self.external_connect]

        for model_id in system.get_connect().keys():
            if model_id not in self.external_connect:
                system.get_model(model_id=model_id).connect = None

        return system

    def clean_contacts(
        self, contactpairs: Dict[float, List[float]]
    ) -> Dict[float, List[float]]:
        """
        Clean merged contacts to prevent multi-counting in system.
        """
        return contactpairs

    def _get_unique_connections(self, system) -> List[str]:
        """
        Get unique connections for all models, including those with only self-connections.
        """
        unique_connections: Set[str] = set()
        for model_id in system.get_models():
            model = system.get_model(model_id=model_id)
            if model.connect:
                connections = sorted(
                    set(model.connect)
                )  # Remove duplicates within each model's connections
                connection_str = " ".join(
                    f"{int(connect)}.caps.pdb" for connect in connections
                )
                unique_connections.add(f"{connection_str} ; {model.type}")
            else:
                unique_connections.add(f"{int(model_id)}.caps.pdb ; {model.type}")
        return sorted(unique_connections)

    def write_connect(self, system=None, connect_file=None):
        """
        Writes system of model connections to file, removing duplicates and including all models.
        """
        connect_file_path = Path(connect_file).with_suffix(".txt")
        unique_connections = self._get_unique_connections(system)

        with open(connect_file_path, "w") as f:
            for connection in unique_connections:
                f.write(f"{connection}\n")

    def print_connection_summary(self, system):
        """
        Prints a summary of connections for all models in the system.
        """
        print("Models in the system:")
        for model_id in system.get_models():
            model = system.get_model(model_id=model_id)
            if model.connect:
                print(
                    f"Model ID: {model_id}, Type: {model.type}, Connect: {model.connect}"
                )
            else:
                print(
                    f"Model ID: {model_id}, Type: {model.type}, Connect: [{model_id}] (self only)"
                )

        print("\nConnections as they will appear in the file:")
        for connection in self._get_unique_connections(system):
            print(connection)

    def get_connect(self, ref_model: Any, model: Any, cut_off: float = 3.0) -> bool:
        """
        Calculates distance between models: distance below cut_off (2.0 A) keep model.

        Args:
            ref_model (Any): Reference model.
            model (Any): Model to compare.
            cut_off (float): Distance cut-off in Angstroms.

        Returns:
            bool: True if models are connected, False otherwise.
        """
        if not hasattr(ref_model, "crosslink") or not hasattr(model, "crosslink"):
            return False
        if not ref_model.crosslink or not model.crosslink:
            return False

        for ref_c, c in product(ref_model.crosslink, model.crosslink):
            if np.linalg.norm(ref_c.position - c.position) < cut_off:
                return True
        return False

    def run_connect(
        self, system: Any, unit_cell: Optional[List[float]] = None
    ) -> Dict[float, List[float]]:
        """
        Wrapper to determine connection between all contacts (1) or
        between an added model and the current system (2).

        Args:
            system (Any): The system object.
            unit_cell (Optional[List[float]]): Unit cell parameters.

        Returns:
            Dict[float, List[float]]: Dictionary of connections.
        """
        has_crosslinks = any(
            hasattr(system.get_model(model_id=model_id), "crosslink")
            and system.get_model(model_id=model_id).crosslink
            for model_id in system.get_models()
        )

        if not has_crosslinks:
            LOG.debug("No crosslinks found in system - skipping connection analysis")
            return {}

        if unit_cell is None:
            return self.get_contact_connect(system=system)
        else:
            return self.get_model_connect(system=system, unit_cell=unit_cell)
