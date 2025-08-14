# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0

from typing import List, Dict, Any, Optional, Tuple, Set, Union
from pathlib import Path
import os

from colbuilder.core.utils.logger import setup_logger
from colbuilder.core.topology.crosslink import Crosslink

LOG = setup_logger(__name__)


class Itp:
    """Class for managing molecular topology files in the Martini 3 force field format.

    This class handles reading, processing, and merging topology files that combine
    Martini 3 force field parameters with Go-like potentials. It supports handling of:
    - Multiple molecular components and their connections
    - Position restraints and various bonded interactions
    - Go-model interactions and exclusions
    - Virtual sites and crosslinking between components

    The class maintains separate data structures for initial component-wise storage
    and final merged topologies, ensuring proper index handling and connectivity.
    """

    def __init__(self, system: Any = None, model_id: Optional[int] = None) -> None:
        """Initialize a new topology processor instance.

        Args:
            system: The molecular system containing component and connectivity information
            model_id: Unique identifier for the molecular model being processed

        The initialization creates empty data structures for:
        - Component-wise storage (atoms, bonds, angles, etc.)
        - Final merged topology elements
        - Crosslinking information and virtual site mappings
        """
        self.system = system

        # Component-wise storage arrays
        self.molecule: List[List[Any]] = self.allocate(model_id=model_id)
        self.mol_ends: List[List[Any]] = self.allocate(model_id=model_id)
        self.atoms: List[List[Any]] = self.allocate(model_id=model_id)
        self.posres: List[List[Any]] = self.allocate(model_id=model_id)
        self.bonds: List[List[Any]] = self.allocate(model_id=model_id)
        self.constraints: List[List[Any]] = self.allocate(model_id=model_id)
        self.angles: List[List[Any]] = self.allocate(model_id=model_id)
        self.exclusions: List[List[Any]] = self.allocate(model_id=model_id)
        self.go_exclusions: List[List[Any]] = self.allocate(model_id=model_id)
        self.dihedrals: List[List[Any]] = self.allocate(model_id=model_id)
        self.virtual_sites: List[List[Any]] = self.allocate(model_id=model_id)
        self.go_table: List[List[Any]] = self.allocate(model_id=model_id)
        self.pairs: List[List[Any]] = self.allocate(model_id=model_id)

        # Final merged topology structures
        self.final_atoms: List[List[Any]] = []
        self.final_posres: List[List[Any]] = []
        self.final_bonds: List[List[Any]] = []
        self.final_flex_bonds: List[List[Any]] = []
        self.final_constraints: List[List[Any]] = []
        self.final_angles: List[List[Any]] = []
        self.final_dihedrals: List[List[Any]] = []
        self.final_exclusions: List[List[Any]] = []
        self.final_go_exclusions: List[List[Any]] = []
        self.final_virtual_sites: List[List[Any]] = []
        self.final_pairs: List[List[Any]] = []

        # Crosslinking (simplified - just store the bonded parameters directly)
        self.crosslink_bonded: Dict[str, List[List[Any]]] = {
            'bonds': [], 'angles': [], 'dihedrals': []
        }
        self.vs_to_col: Dict[str, str] = {}
        self.delta_merge: int = 0
        self.no_line: Tuple[str, ...] = (
            "[",
            "\n",
            "#endif\n",
            "#ifdef",
            "#ifndef",
            "#include",
            ";[",
            ";",
        )

    def allocate(self, model_id: Optional[int] = None) -> List[List[Any]]:
        """Create storage arrays based on the number of molecular connections.

        Args:
            model_id: Identifier for the molecular model

        Returns:
            A list of empty lists, with one list per molecular connection
            in the specified model
        """
        # Allocate for ALL connections (including cross-connections)
        # Even though cross-connections don't have ITP files, we need the slots
        size = len(self.system.get_model(model_id=model_id).connect)
        return [[] for _ in range(size)]

    def read_model(
        self, model_id: Optional[int] = None, system: Optional[Any] = None
    ) -> None:
        """
        Read and merge all connected ITP files for a single model.

        Reads ITP files for ALL connections (both self and cross-connections) to ensure
        that crosslinks referencing atoms from different connections work properly.

        Args:
            model_id: Unique identifier for the model to process
            system: Reserved for future system-level configuration (not currently used)
        """
        cnt_con = 0
        connect_ids = self.system.get_model(model_id=model_id).connect

        for connect_id in connect_ids:
            # Check if ITP files exist for this connection
            itp_path = f"col_{int(model_id)}.{int(connect_id)}.itp"
            excl_path = f"col_{int(model_id)}.{int(connect_id)}_go-excl.itp"
            table_path = f"col_{int(model_id)}.{int(connect_id)}_go-table.itp"
            
            if os.path.exists(itp_path):
                # Read the ITP files for this connection
                LOG.debug(f"Reading ITP files for model {model_id} connection {connect_id}")
                self.read_itp(model_id=model_id, connect_id=connect_id, cnt_con=cnt_con)
                self.read_excl(model_id=model_id, connect_id=connect_id, cnt_con=cnt_con)
                self.read_table(model_id=model_id, connect_id=connect_id, cnt_con=cnt_con)
                LOG.debug(f"Successfully read ITP files for connection {connect_id}")
            else:
                # ITP files don't exist for this connection - leave empty but still allocate space
                LOG.debug(f"No ITP files found for model {model_id} → connection {connect_id} (expected for cross-connections)")
                
            cnt_con += 1
            
        LOG.debug(f"Processed {cnt_con} connections for model {model_id}")
        
        # Log what we actually read
        atoms_read = sum(len(atoms) for atoms in self.atoms if atoms)
        LOG.info(f"Total atoms read from ITP files for model {model_id}: {atoms_read}")

    def read_itp(
        self,
        model_id: Optional[int] = None,
        connect_id: Optional[int] = None,
        cnt_con: Optional[int] = None,
    ) -> None:
        """Read and parse a single ITP (Interaction Parameter) file.

        Reads an ITP file for a specific connection and parses its contents into
        appropriate data structures. Handles various section types including atoms,
        bonds, angles, dihedrals, and other molecular topology parameters.

        Args:
            model_id: Identifier for the model being processed
            connect_id: Identifier for the specific connection
            cnt_con: Counter index for the current connection

        Raises:
            FileNotFoundError: If the specified ITP file cannot be found
        """
        if model_id is None or connect_id is None:
            raise ValueError("model_id and connect_id cannot be None when reading ITP file.")
        itp_path = f"col_{int(model_id)}.{int(connect_id)}.itp"
        bonded_type = ""

        try:
            with open(itp_path, "r") as f:
                for line in f:
                    if line[0] == ";":
                        continue

                    if cnt_con is None:
                        raise ValueError("cnt_con cannot be None when appending to molecule.")
                    self.molecule[cnt_con].append(line)

                    # Parse section headers
                    if line == "[ atoms ]\n":
                        bonded_type = "atoms"
                    elif line == "[ position_restraints ]\n":
                        self.mol_ends[cnt_con] = [
                            int(self.molecule[cnt_con][-3].split(" ")[0])
                        ]
                        bonded_type = "posres"
                    elif line == "[ bonds ]\n":
                        bonded_type = "bonds"
                    elif line == "[ constraints ]\n":
                        bonded_type = "constraints"
                    elif line == "[ virtual_sitesn ]\n":
                        bonded_type = "virtualsites"
                    elif line == "[ angles ]\n":
                        bonded_type = "angles"
                    elif line == "[ dihedrals ]\n":
                        bonded_type = "dihedrals"
                    elif line == "[ exclusions ]\n":
                        bonded_type = "exclusions"

                    # Parse section content
                    if line.split(" ")[0] not in self.no_line:
                        tokens = [k for k in line.split(" ") if k and k != "\n"]

                        # Store tokens in appropriate data structure based on section type
                        if cnt_con is None:
                            raise ValueError("cnt_con cannot be None when parsing ITP sections.")
                        if bonded_type == "atoms":
                            self.atoms[cnt_con].append(tokens)
                        elif bonded_type == "posres":
                            self.posres[cnt_con].append(tokens)
                        elif bonded_type == "bonds":
                            self.bonds[cnt_con].append(tokens)
                        elif bonded_type == "constraints":
                            self.constraints[cnt_con].append(tokens)
                        elif bonded_type == "virtualsites":
                            self.virtual_sites[cnt_con].append(tokens)
                        elif bonded_type == "angles":
                            tokens[-1] = tokens[-1].replace("\n", "")
                            self.angles[cnt_con].append(tokens)
                        elif bonded_type == "exclusions":
                            self.exclusions[cnt_con].append(tokens)
                        elif bonded_type == "dihedrals":
                            self.dihedrals[cnt_con].append(tokens)

        except FileNotFoundError:
            LOG.error(f"ITP file not found: {itp_path}")
            raise

    def read_excl(
        self,
        model_id: Optional[int] = None,
        connect_id: Optional[int] = None,
        cnt_con: Optional[int] = None,
    ) -> None:
        """Read and parse a Go-model exclusion file.

        Processes exclusion definitions that specify which atom pairs should be
        excluded from non-bonded interactions in the Go-model potential.

        Args:
            model_id: Identifier for the model being processed
            connect_id: Identifier for the specific molecular connection
            cnt_con: Connection counter used for array indexing

        Raises:
            FileNotFoundError: If the specified exclusion file cannot be found
        """
        if model_id is None or connect_id is None:
            raise ValueError("model_id and connect_id cannot be None when reading exclusion file.")
        excl_path = f"col_{int(model_id)}.{int(connect_id)}_go-excl.itp"

        try:
            with open(excl_path, "r") as f:
                for line in f:
                    if line.split(" ")[0] not in self.no_line:
                        tokens = [
                            k
                            for k in line.split(" ")
                            if k and k != "\n" and k.strip() != ";"
                        ]
                        if tokens:
                            if cnt_con is None:
                                raise ValueError("cnt_con cannot be None when appending to go_exclusions.")
                            self.go_exclusions[cnt_con].append(tokens)

        except FileNotFoundError:
            LOG.error(f"Exclusion file not found: {excl_path}")
            raise

    def read_table(
        self,
        model_id: Optional[int] = None,
        connect_id: Optional[int] = None,
        cnt_con: Optional[int] = None,
    ) -> None:
        """Read and parse a Go-model interaction table file.

        Processes tabulated potential parameters that define the Go-model interactions
        between specific atom pairs.

        Args:
            model_id: Identifier for the model being processed
            connect_id: Identifier for the specific molecular connection
            cnt_con: Connection counter used for array indexing

        Raises:
            FileNotFoundError: If the specified table file cannot be found
        """
        if model_id is None or connect_id is None:
            raise ValueError("model_id and connect_id cannot be None when reading Go-table file.")
        table_path = f"col_{int(model_id)}.{int(connect_id)}_go-table.itp"

        try:
            with open(table_path, "r") as f:
                for line in f:
                    if line.split(" ")[0] not in self.no_line:
                        tokens = [
                            k
                            for k in line.split(" ")
                            if k and k != "\n" and k.strip() != ";"
                        ]
                        if tokens:
                            if cnt_con is None:
                                raise ValueError("cnt_con cannot be None when appending to go_table.")
                            self.go_table[cnt_con].append(tokens)

        except FileNotFoundError:
            LOG.error(f"Go-table file not found: {table_path}")
            raise

    def go_to_pairs(self, model_id: Optional[int] = None) -> None:
        """Convert Go-model table entries to pair interactions.

        Processes all molecular connections in a model to:
        1. Map virtual sites to corresponding column atoms
        2. Generate pair interaction parameters from Go-model definitions

        Args:
            model_id: Identifier for the model being processed
        """
        num_connections = len(self.system.get_model(model_id=model_id).connect)

        for cnt_con in range(num_connections):
            # Skip empty connections (cross-connections)
            if not self.atoms[cnt_con]:
                continue
            self.match_vs_to_pairs(cnt_con=cnt_con)
            self.get_pairs(cnt_con=cnt_con)

    def match_vs_to_pairs(self, cnt_con: Optional[int] = None) -> None:
        """Map virtual sites to their corresponding column atoms.

        Creates a mapping dictionary that associates virtual site IDs with
        their corresponding column atom indices. This mapping is used for
        converting Go-model interactions into pair parameters.

        Args:
            cnt_con: Counter index for the current molecular connection
        """
        if cnt_con is None:
            raise ValueError("cnt_con cannot be None when mapping virtual sites to pairs.")
        for atom_entry in self.atoms[cnt_con]:
            if atom_entry[1].startswith("col"):
                self.vs_to_col[atom_entry[1]] = atom_entry[0]
                atom_entry[1] = "col"

    def get_pairs(self, cnt_con: Optional[int] = None) -> None:
        """Generate pair interactions from Go-model table entries.

        Processes each Go-model table entry to create pair interaction parameters.
        Each pair entry contains:
        - Mapped atom indices for both interaction partners
        - Interaction parameters from the Go-table
        - Original virtual site IDs for reference

        Args:
            cnt_con: Counter index for the current molecular connection

        Note:
            The pair entry format follows: [atom1, atom2, param1, ..., param6, vs1_id, vs2_id]
            where atom1/2 are mapped column indices and vs1/2_id are original virtual site IDs
        """
        if cnt_con is None:
            raise ValueError("cnt_con cannot be None when generating pairs.")
        for table_entry in self.go_table[cnt_con]:
            vs1, vs2 = table_entry[0:2]

            pair_entry = [
                self.vs_to_col[vs1],  # First atom index
                self.vs_to_col[vs2],  # Second atom index
                *table_entry[2:8],  # Interaction parameters
                vs1,  # Original virtual site ID 1
                vs2,  # Original virtual site ID 2
            ]

            self.pairs[cnt_con].append(pair_entry)

    def merge_topology(self, cnt_con: Optional[int] = None) -> None:
        """
        Merge a connection's topology with the accumulated topology.

        Processes a single molecular connection's topology data and merges it into
        the final topology structures. Handles index adjustments and special cases
        for different interaction types.

        Args:
            cnt_con: Counter index for the current molecular connection being merged

        Note:
            Now handles both self-connections and cross-connections that have ITP data
        """
        # Skip if this connection has no data 
        if cnt_con is not None and not self.atoms[cnt_con]:
            LOG.debug(f"Skipping empty connection {cnt_con} in merge_topology")
            return
        
        # Update index offset based on previous connections that had data
        if cnt_con is not None and cnt_con != 0:
            # Find the last non-empty connection
            for i in range(cnt_con - 1, -1, -1):
                if self.mol_ends[i]:
                    prev_end = self.mol_ends[i]
                    if isinstance(prev_end, list):
                        prev_end = int(prev_end[0])
                    else:
                        prev_end = int(prev_end)
                    self.delta_merge += prev_end
                    LOG.debug(f"Adding {prev_end} to delta_merge from connection {i}, new delta: {self.delta_merge}")
                    break

        # Ensure cnt_con is not None before using as index
        if cnt_con is None:
            raise ValueError("cnt_con cannot be None when merging atoms.")
        
        LOG.debug(f"Merging connection {cnt_con} with {len(self.atoms[cnt_con])} atoms, delta_merge: {self.delta_merge}")
        
        # Process atoms with index adjustment
        merged_atoms = [
            [
                int(a[0]) + self.delta_merge,
                a[1],
                a[2],
                a[3],
                a[4],
                int(a[5]) + self.delta_merge,
                a[6],
            ]
            for a in self.atoms[cnt_con]
        ]
        self.final_atoms.extend(merged_atoms)

        # Process position restraints
        merged_posres = [
            [int(p[0]) + self.delta_merge, p[1], p[2], p[3], p[4]]
            for p in self.posres[cnt_con]
        ]
        self.final_posres.extend(merged_posres)

        # Process bonds and separate flexible bonds
        merged_bonds = [
            [
                int(b[0]) + self.delta_merge,
                int(b[1]) + self.delta_merge,
                b[2],
                b[3],
                b[4],
            ]
            for b in self.bonds[cnt_con]
        ]
        for bond in merged_bonds:
            if int(bond[-1]) == 1000000:
                self.final_flex_bonds.append(bond)
            else:
                self.final_bonds.append(bond)

        # Process angles
        merged_angles = [
            [
                int(a[0]) + self.delta_merge,
                int(a[1]) + self.delta_merge,
                int(a[2]) + self.delta_merge,
                a[3],
                a[4],
                a[5] + "\n",
            ]
            for a in self.angles[cnt_con]
        ]
        self.final_angles.extend(merged_angles)

        # Process dihedrals based on entry length
        merged_dihedrals = []
        if self.dihedrals[cnt_con]:  # Check if not empty
            dihedral_length = len(self.dihedrals[cnt_con][0])

            if dihedral_length == 8:
                for dih in self.dihedrals[cnt_con]:
                    indices = [int(idx) + self.delta_merge for idx in dih[:4]]
                    merged_dih = [*indices, *dih[4:]]
                    merged_dihedrals.append(merged_dih)

            elif dihedral_length == 7:
                merged_dihedrals = [
                    [
                        int(d[0]) + self.delta_merge,
                        int(d[1]) + self.delta_merge,
                        int(d[2]) + self.delta_merge,
                        int(d[3]) + self.delta_merge,
                        d[4],
                        d[5],
                        d[6],
                    ]
                    for d in self.dihedrals[cnt_con]
                ]

            elif dihedral_length == 6:
                merged_dihedrals = [
                    [
                        int(d[0]) + self.delta_merge,
                        int(d[1]) + self.delta_merge,
                        int(d[2]) + self.delta_merge,
                        int(d[3]) + self.delta_merge,
                        d[4],
                        d[5],
                    ]
                    for d in self.dihedrals[cnt_con]
                ]

        self.final_dihedrals.extend(merged_dihedrals)

        # Process constraints
        merged_constraints = [
            [int(c[0]) + self.delta_merge, int(c[1]) + self.delta_merge, c[2], c[3]]
            for c in self.constraints[cnt_con]
        ]
        self.final_constraints.extend(merged_constraints)

        # Process virtual sites
        merged_vsites = [
            [
                str(int(v[0]) + self.delta_merge),
                str(int(1)),
                str(int(v[2]) + self.delta_merge) + "\n",
            ]
            for v in self.virtual_sites[cnt_con]
        ]
        self.final_virtual_sites.extend(merged_vsites)

        # Process Go exclusions
        merged_go_excl = [
            [
                int(e[0]) + self.delta_merge,
                int(e[1]) + self.delta_merge,
                e[2],
                int(e[3]) + self.delta_merge,
                int(e[4]) + self.delta_merge,
            ]
            for e in self.go_exclusions[cnt_con]
        ]
        self.final_go_exclusions.extend(merged_go_excl)

        # Process exclusions
        merged_excl = []
        for excl in self.exclusions[cnt_con]:
            merged_entry = [str(int(idx) + self.delta_merge) for idx in excl if idx]
            if merged_entry:
                merged_entry[-1] += "\n"
                merged_excl.append(merged_entry)
        self.final_exclusions.extend(merged_excl)

        # Process pairs
        merged_pairs = [
            [
                int(p[0]) + self.delta_merge,
                int(p[1]) + self.delta_merge,
                p[2],
                format(float(p[3]), ".10f"),
                format(float(p[4]), ".10f"),
                ";",
                p[-2],
                p[-1] + "\n",
            ]
            for p in self.pairs[cnt_con]
        ]
        self.final_pairs.extend(merged_pairs)
        
        LOG.debug(f"Merged connection {cnt_con}: {len(merged_atoms)} atoms, final total: {len(self.final_atoms)}")

    def make_topology(
        self, model_id: Optional[int] = None, cnt_model: Optional[int] = None
    ) -> None:
        """
        Create a complete topology by merging connections and adding crosslinks.

        Creates a complete molecular topology by:
        1. Merging all component topologies with proper index adjustments
        2. Setting up crosslink structures (simplified approach)
        3. Writing the final topology and exclusion files

        Args:
            model_id: Identifier for the molecular model being processed
            cnt_model: Counter index used for output file naming
        """
        # Merge all connection topologies first
        for cnt_con in range(len(self.system.get_model(model_id=model_id).connect)):
            self.merge_topology(cnt_con=cnt_con)
        
        # Set up crosslinks using simplified approach (following working version)
        if len(self.system.get_model(model_id=model_id).connect) == 1:
            self.crosslink_bonded = {k: [] for k in ["bonds", "angles", "dihedrals"]}
            LOG.debug(f"Single connection model {model_id}: no crosslinks")
        else:
            try:
                crosslinker = Crosslink(cnt_model=cnt_model)
                self.crosslink_bonded = crosslinker.set_crosslink_bonded(cnt_model=cnt_model)
                
                # Log crosslink information for debugging
                if any(self.crosslink_bonded[k] for k in ['bonds', 'angles', 'dihedrals']):
                    LOG.info(f"Found crosslinks for model {model_id}:")
                    LOG.info(f"  Bonds: {len(self.crosslink_bonded['bonds'])}")
                    LOG.info(f"  Angles: {len(self.crosslink_bonded['angles'])}")
                    LOG.info(f"  Dihedrals: {len(self.crosslink_bonded['dihedrals'])}")
                else:
                    LOG.debug(f"No crosslinks found for model {model_id}")
                    
            except Exception as e:
                LOG.warning(f"Could not process crosslinks for model {model_id}: {str(e)}")
                # Initialize empty crosslink structures if processing fails
                self.crosslink_bonded = {k: [] for k in ["bonds", "angles", "dihedrals"]}

        # Write final topology files
        self.write_topology(cnt_model=cnt_model)
        self.write_excl(cnt_model=cnt_model)

    def write_topology(self, cnt_model: Optional[int] = None) -> None:
        """Write the complete molecular topology to an ITP file.

        Creates a structured topology file containing all merged molecular components
        and their interactions. Uses the simplified crosslink format from the working version.

        Args:
            cnt_model: Model counter used for output file naming

        Raises:
            PermissionError: If writing to the output file is not permitted
            Exception: For other file operation errors
        """
        if cnt_model is None:
            raise ValueError("cnt_model cannot be None when writing topology file.")
        output_path = f"col_{int(cnt_model)}.itp"

        try:
            with open(output_path, "w") as f:
                f.write("; Merging of topologies for models due to system\n")
                f.write("[ moleculetype ]\n")
                f.write(f"col_{cnt_model} 1\n")

                f.write("\n\n[ atoms ]\n")
                for atom in self.final_atoms:
                    f.write(
                        "{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}\n".format(
                            *[atom[i] for i in range(7)]
                        )
                    )

                f.write("\n[ position_restraints ]\n")
                f.write("#ifdef POSRES\n")
                for posre in self.final_posres:
                    f.write(" ".join(str(i) for i in posre))
                f.write("#endif\n")

                f.write("\n[ bonds ]\n")
                for bond in self.final_bonds:
                    f.write(" ".join(str(i) for i in bond))

                f.write("#ifdef FLEXIBLE\n; side chain flexible\n")
                for flex_bond in self.final_flex_bonds:
                    f.write(" ".join(str(i) for i in flex_bond))
                f.write("#endif\n")

                # Crosslink bonds (simplified format from working version)
                f.write("; crosslink bonds \n")
                for cb in self.crosslink_bonded['bonds']:
                    f.write(" ".join(str(i) for i in cb))

                f.write("\n[ pairs ]\n")
                for pair in self.final_pairs:
                    f.write(" ".join(str(i) for i in pair))

                f.write("\n[ constraints ]\n")
                f.write("#ifndef FLEXIBLE\n")
                for constraint in self.final_constraints:
                    f.write(" ".join(str(i) for i in constraint))
                f.write("#endif\n")

                f.write("\n[ virtual_sitesn ]\n")
                for vsite in self.final_virtual_sites:
                    f.write(" ".join(str(i) for i in vsite))

                f.write("\n[ angles ]\n")
                for angle in self.final_angles:
                    f.write(" ".join(str(i) for i in angle))

                # Crosslink angles (simplified format from working version)
                f.write("; crosslink angles \n")
                for ca in self.crosslink_bonded['angles']:
                    f.write(" ".join(str(i) for i in ca))

                f.write("\n[ dihedrals ]\n")
                for dihedral in self.final_dihedrals:
                    f.write(" ".join(str(i) for i in dihedral))

                # Crosslink dihedrals (simplified format from working version)
                f.write("; crosslink dihedrals \n")
                for cd in self.crosslink_bonded['dihedrals']:
                    f.write(" ".join(str(i) for i in cd))

                f.write("\n[ exclusions ]\n")
                for exclusion in self.final_exclusions:
                    f.write(" ".join(str(i) for i in exclusion))

        except PermissionError:
            LOG.error(f"Permission denied when writing to topology file: {output_path}")
            raise
        except Exception as e:
            LOG.error(f"Error writing topology file: {str(e)}")
            raise

    def write_excl(self, cnt_model: Optional[int] = None) -> None:
        """Write the merged Go-exclusions to an ITP file.

        Creates a file containing exclusion definitions for Go-model interactions,
        specifying which atom pairs should be excluded from non-bonded interactions.
        Each exclusion entry is written in space-separated format.

        Args:
            cnt_model: Model counter used for output file naming

        Raises:
            PermissionError: If writing to the output file is not permitted
            Exception: For other file operation errors
        """
        output_path = f"col_{cnt_model}_go-excl.itp"

        try:
            with open(output_path, "w") as f:
                f.write(";[ exclusions ]\n")
                for exclusion in self.final_go_exclusions:
                    # Write space-separated items, ensuring no trailing space
                    exclusion_str = " ".join(str(item) for item in exclusion)
                    f.write(f"{exclusion_str}\n")

        except PermissionError:
            LOG.error(
                f"Permission denied when writing to exclusion file: {output_path}"
            )
            raise
        except Exception as e:
            LOG.error(f"Error writing exclusion file: {str(e)}")
            raise