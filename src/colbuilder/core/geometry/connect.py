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
        self.marker_resnames: Set[str] = set(
            [
                "AGS",
                "APD",
                "LGX",
                "LPS",
                "LZS",
                "LZD",
                "L5Y",
                "L4Y",
                "L5X",
                "LY5",
                "LX5",
                "LY4",
                "LX4",
                "LYX",
                "LY3",
                "LY2",
                "LX3",
                "LX2",
                "LXX",
                "L3Y",
                "L2Y",
                "LXY",
                "L3X",
                "L2X",
                "LYY",
            ]
        )

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

    def _model_has_marker(self, model_id: float, caps_dir: Optional[Path], system: Any) -> bool:
        """
        Check if a model has any marker residue. Prefer scanning the caps file
        if available; fallback to the in-memory crosslink list.
        """
        candidate_paths: List[Path] = []
        if caps_dir:
            try:
                model_type = getattr(system.get_model(model_id=model_id), "type", None)
            except Exception:
                model_type = None

            candidate_paths.append(caps_dir / f"{int(float(model_id))}.caps.pdb")
            if model_type:
                candidate_paths.append(
                    caps_dir / str(model_type) / f"{int(float(model_id))}.caps.pdb"
                )
            # As a fallback, check all subdirectories in caps_dir
            if caps_dir.exists():
                for sub in caps_dir.iterdir():
                    if sub.is_dir():
                        candidate_paths.append(
                            sub / f"{int(float(model_id))}.caps.pdb"
                        )

            for caps_path in candidate_paths:
                if caps_path.exists():
                    try:
                        with caps_path.open("r") as f:
                            for line in f:
                                if not line.startswith(("ATOM", "HETATM")):
                                    continue
                                resn = line[17:20].strip()
                                if resn in self.marker_resnames:
                                    return True
                        # If we found a caps file but no marker, treat as marker-less
                        return False
                    except Exception:
                        continue

        # Only fall back to crosslink attribute if we could not read any caps file
        try:
            model_obj = system.get_model(model_id=model_id)
            return bool(getattr(model_obj, "crosslink", []))
        except Exception:
            return False

    def _connections_from_model_graph(self, system: Any) -> List[str]:
        """
        Fallback: build connectivity lines from the in-memory model.connect graph.
        Mirrors the behaviour of the previous implementation, ensuring geometry_gen
        runs without marker information still get usable connectivity.
        """
        unique_connections: Set[str] = set()
        for model_id in system.get_models():
            model_obj = system.get_model(model_id=model_id)
            if model_obj.connect:
                connections = sorted(set(model_obj.connect))
                connection_str = " ".join(
                    f"{int(connect)}.caps.pdb" for connect in connections
                )
                unique_connections.add(f"{connection_str} ; {model_obj.type}")
            else:
                unique_connections.add(f"{int(model_id)}.caps.pdb ; {model_obj.type}")
        return sorted(unique_connections)

    def _get_unique_connections(
        self, system: Any, connect_file_path: Optional[Path]
    ) -> List[str]:
        """
        Get unique connections for all models, including those with only self-connections.

        Prefers marker-based detection from caps files when markers are present;
        falls back to the system's model.connect graph when no markers are found
        (e.g., pure geometry_gen runs without crosslink markers).
        """
        caps_dir = connect_file_path.parent if connect_file_path else None
        if not caps_dir or not caps_dir.exists():
            return self._connections_from_model_graph(system)

        system_ids = set(system.get_models()) if system is not None else set()

        caps_map: Dict[float, Path] = {}
        # Include caps in current dir and immediate subdirs (A/B/C)
        for path in caps_dir.glob("*.caps.pdb"):
            try:
                mid = float(path.stem.split(".")[0])
            except ValueError:
                continue
            if system_ids and mid not in system_ids:
                continue
            caps_map[mid] = path
        for sub in caps_dir.iterdir():
            if sub.is_dir():
                for path in sub.glob("*.caps.pdb"):
                    try:
                        mid = float(path.stem.split(".")[0])
                    except ValueError:
                        continue
                    if system_ids and mid not in system_ids:
                        continue
                    caps_map[mid] = path

        marker_positions: Dict[float, List[np.ndarray]] = {}
        markerless: Set[float] = set()
        for mid, pth in caps_map.items():
            try:
                with pth.open("r") as f:
                    for line in f:
                        if not line.startswith(("ATOM", "HETATM")):
                            continue
                        resn = line[17:20].strip()
                        if resn not in self.marker_resnames:
                            continue
                        try:
                            pos = np.array(
                                [
                                    float(line[30:38]),
                                    float(line[38:46]),
                                    float(line[46:54]),
                                ]
                            )
                            marker_positions.setdefault(mid, []).append(pos)
                        except Exception:
                            continue
            except Exception:
                continue
            if mid not in marker_positions:
                markerless.add(mid)

        # No markers found anywhere -> fallback to existing connect graph
        if not marker_positions and not markerless:
            return self._connections_from_model_graph(system)

        connect_pairs: Dict[float, Set[float]] = {
            mid: set() for mid in marker_positions
        }
        mids = sorted(marker_positions.keys())
        for i, m1 in enumerate(mids):
            for m2 in mids[i + 1 :]:
                pos1 = marker_positions.get(m1, [])
                pos2 = marker_positions.get(m2, [])
                if not pos1 or not pos2:
                    continue
                min_dist = min(
                    np.linalg.norm(a - b) for a in pos1 for b in pos2
                )
                if min_dist < 4.0:
                    connect_pairs[m1].add(m2)
                    connect_pairs[m2].add(m1)

        unique_connections: Set[str] = set()
        for mid in mids:
            if system_ids and mid not in system_ids:
                continue
            model_type = getattr(system.get_model(model_id=mid), "type", "NC")
            partners = sorted(connect_pairs.get(mid, set()))
            if partners:
                parts = [f"{int(mid)}.caps.pdb"] + [
                    f"{int(p)}.caps.pdb" for p in partners
                ]
                unique_connections.add(" ".join(parts) + f" ; {model_type}")
            else:
                unique_connections.add(f"{int(mid)}.caps.pdb ; {model_type}")

        for mid in sorted(markerless):
            model_type = getattr(system.get_model(model_id=mid), "type", "NC")
            unique_connections.add(f"{int(mid)}.caps.pdb ; {model_type}")

        return sorted(unique_connections)

    def write_connect(self, system=None, connect_file=None):
        """
        Writes system of model connections to file, removing duplicates and including all models.
        """
        connect_file_path = Path(connect_file).with_suffix(".txt")
        unique_connections = self._get_unique_connections(system, connect_file_path)

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
