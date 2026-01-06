"""
Configuration Utilities for the ColBuilder Pipeline

This module provides utilities for managing and validating configuration settings in the ColBuilder
pipeline. It defines the `ColbuilderConfig` class, which encapsulates all configuration options,
validates input parameters, and resolves paths. The module also includes helper functions for
loading, validating, and updating configurations.

Key Features:
--------------
1. **Configuration Management**:
   - `ColbuilderConfig`: Centralized configuration class for managing pipeline settings.
   - Supports multiple operation modes: SEQUENCE, GEOMETRY, TOPOLOGY, MIX, and REPLACE.
   - Automatically resolves relative paths to absolute paths based on the working directory.

2. **Validation**:
   - Validates input files (FASTA, PDB) and configuration parameters.
   - Ensures proper values for critical parameters like `fibril_length`, `contact_distance`, and `solution_space`.
   - Raises `ConfigurationError` for invalid inputs or missing required files.

3. **Path Handling**:
   - Resolves paths for project data, working directories, and output files.
   - Provides default paths for essential resources like templates, force fields, and Chimera scripts.

4. **Dynamic Updates**:
   - Allows runtime updates to configuration settings.
   - Automatically adjusts operation modes based on updated flags.

5. **Integration with External Tools**:
   - Detects and configures Martinize2 for topology generation.
   - Supports Conda environments for Martinize2 execution.

Usage:
------
This module is designed to be used throughout the ColBuilder pipeline to manage and validate
configuration settings. It ensures consistency and simplifies error handling.

Example:
--------
Example:
--------
```python
from colbuilder.core.utils.config import get_config

# Initialize configuration
config = get_config(
    working_directory="/path/to/working_dir",
    sequence_generator=True,
    species="homo_sapiens",
    fasta_file="input.fasta"
)

# Access configuration attributes
print(config.working_directory)  # Output: /path/to/working_dir
print(config.mode)  # Output: OperationMode.SEQUENCE

# Validate paths and input files
config.validate_paths()

# Update configuration dynamically
config.update({"debug": True, "fibril_length": 300})
print(config.debug)  # Output: True
print(config.fibril_length)  # Output: 300

# Get output path for a file
output_path = config.get_output_path("output_file", ".pdb")
print(output_path)  # Output: /path/to/working_dir/output_file.pdb
```
"""

# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0

from enum import Flag, auto, Enum
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Union, Literal, Any
from pydantic import (
    BaseModel,
    Field,
    field_validator,
    computed_field,
    model_validator,
    ValidationInfo,
)
import yaml
import re
import os
import subprocess
import json
from functools import lru_cache

from .validators import BioformatValidator
from .exceptions import (
    ConfigurationError,
    SequenceGenerationError,
    GeometryGenerationError,
    TopologyGenerationError,
)

from colbuilder.core.utils.logger import setup_logger

LOG = setup_logger(__name__)


class OperationMode(Flag):
    """Operation modes for the Colbuilder pipeline."""

    NONE = 0
    SEQUENCE = auto()
    GEOMETRY = auto()
    TOPOLOGY = auto()
    MIX = auto()
    REPLACE = auto()


def resolve_relative_paths(config: Dict[str, Any], base_dir: Path) -> Dict[str, Any]:
    """Resolve relative paths in config to be absolute based on the given base directory."""
    path_keys = [
        "working_directory",
        "fasta_file",
        "mutated_pdb",
        "pdb_file",
        "crystalcontacts_file",
        "connect_file",
        "replace_file",
    ]

    result = config.copy()
    for key in path_keys:
        if (
            key in result
            and result[key] is not None
            and not isinstance(result[key], bool)
        ):
            path_str = str(result[key])
            # Skip empty strings
            if not path_str.strip():
                continue

            path = Path(path_str)
            if not path.is_absolute():
                result[key] = str((base_dir / path).resolve())

    # Handle list of paths in files_mix
    if "files_mix" in result and result["files_mix"]:
        result["files_mix"] = [
            (
                str((base_dir / Path(file_path)).resolve())
                if not Path(file_path).is_absolute()
                else file_path
            )
            for file_path in result["files_mix"]
        ]

    return result


class ColbuilderConfig(BaseModel):
    """Main configuration class for the Colbuilder pipeline."""

    # Operation mode
    mode: Optional[OperationMode] = Field(None, description="Operation mode")
    debug: bool = Field(default=False, description="Enable debug logging")
    working_directory: Optional[Path] = Field(
        default=Path.cwd(), description="Working directory"
    )
    config_file: Optional[Path] = Field(None, description="YAML configuration file")

    # PDB file generation mode (from sequence)
    sequence_generator: bool = Field(
        default=False, description="Run sequence generation"
    )
    species: str = Field(..., description="Species name")
    fasta_file: Optional[str] = Field(None, description="FASTA file")
    crosslink: bool = Field(default=False, description="Apply crosslinks")
    n_term_type: Optional[str] = Field(None, description="N-terminal type")
    c_term_type: Optional[str] = Field(None, description="C-terminal type")
    n_term_combination: Optional[str] = Field(
        None, description="N-terminal combination"
    )
    c_term_combination: Optional[str] = Field(
        None, description="C-terminal combination"
    )
    
    # Crosslink optimization settings
    crosslink_copies: List[str] = Field(
        default=["D0", "D5"],
        description="Pair of unit cell translations for crosslink optimization. "
                    "Must be exactly 2 elements from: D0, D1, D2, D3, D4, D5"
    )
    # Pre-mutated PDB input (skips homology modeling)
    mutated_pdb: Optional[Path] = Field(
        None, 
        description="Path to pre-mutated PDB file. If provided, skips homology modeling and applies only additional crosslinks"
    )
    
    # Additional crosslinks (applied when using mutated_pdb)
    additional_1_type: Optional[str] = Field(
        None, 
        description="First additional crosslink type to apply to mutated_pdb"
    )
    additional_2_type: Optional[str] = Field(
        None, 
        description="Second dditional crosslink type to apply to mutated_pdb"
    )
    additional_1_combination: Optional[str] = Field(
        None, 
        description="First additional crosslink position"
    )
    additional_2_combination: Optional[str] = Field(
        None, 
        description="Second additional crosslink position"
    )

    # Fibril geometry generation mode
    geometry_generator: bool = Field(
        default=False, description="Run geometry generation"
    )
    pdb_file: Optional[Path] = Field(None, description="Input PDB file")
    contact_distance: Optional[float] = Field(
        None, description="Contact distance for microfibril"
    )
    fibril_length: float = Field(
        default=None,  # Remove default here
        description="Length of microfibril in nanometers",
    )
    crystalcontacts_file: Optional[Path] = Field(
        None, description="Crystal contacts file"
    )
    connect_file: Optional[Path] = Field(None, description="Connect file")
    crystalcontacts_optimize: bool = Field(
        default=False, description="Optimize crystal contacts"
    )
    solution_space: Union[List[float], Tuple[float, float, float]] = Field(
        default=(1, 1, 1), description="Solution space"
    )
    pdb_first_line: Optional[str] = Field(
        default="CRYST1   39.970   26.950  677.900  89.24  94.59 105.58 P 1           2",
        description="Crystal contacts information",
    )

    # Mix crosslinks mode
    mix_bool: bool = Field(
        default=False, description="Generate a mixed crosslinked microfibril"
    )
    ratio_mix: Optional[Union[str, Dict[str, int]]] = Field(
        default=None, description="Ratio mix for crosslinks"
    )
    files_mix: Optional[Union[List[Path], Tuple[Path, ...]]] = Field(
        default=None, description="PDB files with different crosslink types"
    )

    # Replace crosslinks mode
    replace_bool: bool = Field(
        default=False, description="Generate a microfibril with less crosslinks"
    )
    auto_fix_unpaired: bool = Field(
        default=False,
        description="Automatically detect unpaired enzymatic crosslinks and replace them",
    )
    ratio_replace_scope: Literal["enzymatic", "non_enzymatic", "all"] = Field(
        default="non_enzymatic",
        description="Scope of residues considered for ratio-based replacement",
    )
    ratio_replace: Optional[float] = Field(
        None, description="Ratio of crosslinks to be replaced"
    )
    replace_file: Optional[Path] = Field(
        None, description="File with crosslinks to be replaced"
    )
    manual_replacements: Optional[List[str]] = Field(
        default=None,
        description="Manual replacement directives (<pdb> <RES> <resid> <chain>)",
    )

    # Topology generation mode
    topology_generator: bool = Field(
        default=False, description="Generate topology files"
    )
    force_field: Optional[str] = Field(None, description="Force field to be used")
    topology_debug: Optional[bool] = Field(
        default=False,
        description="Save all intermediate files used in building topology",
    )
    martinize2_command: Optional[str] = Field(
        None, description="Detected Martinize2 command"
    )
    martinize2_env: Optional[str] = Field(
        None, description="Detected Martinize2 environment"
    )
    use_conda_run: bool = Field(
        default=False, description="Whether to use conda run for Martinize2"
    )
    go_epsilon: float = Field(
        default=9.414,
        description="GO epsilon value for Martini3-CG parametrization (default: 9.414)",
    )

    # Path Configuration
    PROJECT_ROOT: Path = Field(
        default_factory=lambda: Path(__file__).resolve().parent.parent.parent
    )
    DATA_DIR: Path = Field(
        default_factory=lambda: Path(__file__).resolve().parent.parent.parent / "data"
    )
    HOMOLOGY_LIB_DIR: Path = Field(
        default_factory=lambda: Path(__file__).resolve().parent.parent.parent
        / "data"
        / "sequence"
    )
    TEMPLATE_PDB_PATH: Path = Field(
        default_factory=lambda: Path(__file__).resolve().parent.parent.parent
        / "data"
        / "sequence"
        / "template.pdb"
    )
    TEMPLATE_FASTA_PATH: Path = Field(
        default_factory=lambda: Path(__file__).resolve().parent.parent.parent
        / "data"
        / "sequence"
        / "template.fasta"
    )
    RESTYP_LIB_PATH: Path = Field(
        default_factory=lambda: Path(__file__).resolve().parent.parent.parent
        / "data"
        / "sequence"
        / "modeller"
        / "restyp_mod.lib"
    )
    TOP_HEAV_LIB_PATH: Path = Field(
        default_factory=lambda: Path(__file__).resolve().parent.parent.parent
        / "data"
        / "sequence"
        / "modeller"
        / "top_heav_mod.lib"
    )
    PAR_MOD_LIB_PATH: Path = Field(
        default_factory=lambda: Path(__file__).resolve().parent.parent.parent
        / "data"
        / "sequence"
        / "modeller"
        / "par_mod.lib"
    )
    CROSSLINKS_FILE: Path = Field(
        default_factory=lambda: Path(__file__).resolve().parent.parent.parent
        / "data"
        / "sequence"
        / "crosslinks.csv"
    )
    CHIMERA_SCRIPTS_DIR: Path = Field(
        default_factory=lambda: Path(__file__).resolve().parent.parent.parent
        / "chimera_scripts"
    )
    FORCE_FIELD_DIR: Path = Field(
        default_factory=lambda: Path(__file__).resolve().parent.parent.parent
        / "data"
        / "topology"
    )

    _species_map: Set[str] = {
        "homo_sapiens",
        "ailuropoda_melanoleuca",
        "danio_rerio",
        "mus_musculus",
        "mustela_putorius",
        "myotis_lucifugus",
        "otolemur_garnettii",
        "pan_troglodytes",
        "pongo_abelii",
        "rattus_norvegicus",
        "bos_taurus",
        "callithrix_jacchus",
        "canis_lupus",
        "loxodonta_africana",
        "oreochromis_niloticus",
        "oryzias_latipes",
        "pelodiscus_sinensis",
        "tetraodon_nigroviridis",
        "xiphophorus_maculatus",
    }

    def __init__(self, **data):
        LOG.debug(f"Initializing ColbuilderConfig with data: {data}")
        if "working_directory" in data:
            data["working_directory"] = Path(data["working_directory"]).resolve()
        else:
            data["working_directory"] = Path.cwd().resolve()

        super().__init__(**data)
        if self.mix_bool and self.ratio_mix is not None:
            self.ratio_mix = self._convert_ratio_mix(self.ratio_mix)
        self.solution_space = self._convert_to_tuple(self.solution_space)
        self.files_mix = tuple(self.files_mix) if self.files_mix else None
        self.set_mode()

    def get_project_data_path(self, relative_path: str) -> Path:
        """
        Get the path to a project data file, checking multiple locations.

        Args:
            relative_path: Path relative to the data directory

        Returns:
            Resolved path to the file
        """
        working_dir_path = self.working_directory / os.path.basename(relative_path)
        if working_dir_path.exists():
            return working_dir_path

        data_path = self.DATA_DIR / relative_path
        if data_path.exists():
            return data_path

        return data_path

    def get_output_path(self, basename: str, suffix: Optional[str] = None) -> Path:
        """
        Get a path for an output file in the working directory.

        Args:
            basename: Base name for the file
            suffix: Optional suffix to add (e.g., file extension)

        Returns:
            Path object for the output file
        """
        full_name = f"{basename}{suffix if suffix else ''}"
        return self.working_directory / full_name

    def model_post_init(self, __context: Any) -> None:
        """Post initialization validation and processing."""
        # Only validate species and handle fasta_file if sequence_generator is True
        if self.sequence_generator:
            if self.fasta_file is None:
                if self.species in self._species_map:
                    fasta_name = f"{self.species.replace('_', '')}.fasta"
                    self.fasta_file = str(
                        self.HOMOLOGY_LIB_DIR / "fasta_sequences" / fasta_name
                    )
                else:
                    error_info = ConfigurationError.get_error_info("CFG_ERR_003")
                    custom_message = f"Must provide fasta_file when using custom species: {self.species}"
                    full_message = f"ERROR: {error_info.message}. {custom_message}"
                    raise ConfigurationError(
                        message=full_message, error_code="CFG_ERR_003"
                    )
        # If sequence_generator is False, validate that PDB file is provided
        elif not self.sequence_generator and self.species not in self._species_map:
            if self.pdb_file is None:
                raise ConfigurationError(
                    "When using a custom species with sequence_generator=False, a PDB file must be provided",
                    error_code="CFG_ERR_003",
                )

        self.set_mode()

    @model_validator(mode="after")
    def validate_geometry_requirements(self) -> "ColbuilderConfig":
        """Validate geometry generation requirements."""
        if self.geometry_generator:
            if self.contact_distance is None and self.crystalcontacts_file is None:
                raise ConfigurationError(
                    "Either contact_distance or crystalcontacts_file must be provided for geometry generation",
                    error_code="CFG_ERR_006",
                )
            if self.contact_distance == 0 and self.crystalcontacts_file is None:
                pass
        return self

    @field_validator("contact_distance")
    def validate_contact_distance(cls, value: Optional[float]) -> Optional[float]:
        """Validate contact distance value."""
        if value is not None:
            if not isinstance(value, (int, float)):
                raise ConfigurationError(
                    "Contact distance must be a numeric value", error_code="CFG_ERR_006"
                )
            if value <= 0:
                raise ConfigurationError(
                    f"Contact distance must be positive, got {value}",
                    error_code="CFG_ERR_006",
                )
        return value

    @field_validator("fibril_length")
    def validate_fibril_length(cls, v: float) -> float:
        """Validate fibril length."""
        if v <= 0:
            raise ValueError("fibril_length must be greater than 0")
        LOG.debug(f"Validating fibril_length: {v}")
        return float(v)

    @field_validator("force_field")
    def validate_force_field(cls, value: Optional[str]) -> Optional[str]:
        """Validate force field value."""
        valid_fields = {"amber99", "martini3"}
        if value is not None and value not in valid_fields:
            raise ConfigurationError(
                f"Force field must be one of {valid_fields}, got {value}",
                error_code="CFG_ERR_006",
            )
        return value

    @model_validator(mode="after")
    def validate_species_requirements(self) -> "ColbuilderConfig":
        """Validate species and related requirements."""
        # If sequence generator is True, species must be from predefined list or fasta_file must be provided
        if self.sequence_generator:
            if self.species not in self._species_map and self.fasta_file is None:
                raise ConfigurationError(
                    f"When sequence_generator is True, species must be one of {self._species_map} or fasta_file must be provided",
                    error_code="CFG_ERR_003",
                )
        # If sequence generator is False and species is not from predefined list, PDB file must be provided
        elif self.species not in self._species_map:
            if self.pdb_file is None:
                raise ConfigurationError(
                    "When using a custom species with sequence_generator=False, a PDB file must be provided",
                    error_code="CFG_ERR_003",
                )
        return self

    @field_validator("species", mode="before")
    def convert_species_to_lowercase(cls, value: str) -> str:
        """Convert species name to lowercase."""
        return value.lower() if value else value

    def _convert_ratio_mix(self, value: Union[str, Dict[str, int]]) -> Dict[str, int]:
        """Convert string ratio mix to dictionary."""
        if isinstance(value, str):
            try:
                return {
                    item.split(":")[0]: int(item.split(":")[1])
                    for item in value.split()
                }
            except (ValueError, IndexError):
                raise ConfigurationError(
                    "Invalid ratio_mix format. Expected 'Type:percentage Type:percentage'",
                    error_code="CFG_ERR_004",
                )
        elif isinstance(value, dict):
            return value
        else:
            raise ConfigurationError(
                f"Invalid ratio_mix type. Expected string or dictionary, got {type(value).__name__}",
                error_code="CFG_ERR_004",
            )

    def _convert_to_tuple(
        self, value: Union[List[float], Tuple[float, float, float]]
    ) -> Tuple[float, float, float]:
        """Convert list to tuple for solution space."""
        if isinstance(value, list):
            return tuple(value)
        return value

    @field_validator("solution_space", mode="before")
    def validate_solution_space(cls, value):
        """Validate solution space format."""
        if isinstance(value, list):
            return tuple(value)
        return value
    
    @field_validator("crosslink_copies")
    def validate_crosslink_copies(cls, value: List[str]) -> List[str]:
        """Validate crosslink copies - must be exactly 2 valid translation names."""
        valid_names = ["D0", "D1", "D2", "D3", "D4", "D5"]
        if len(value) != 2:
            raise ConfigurationError(
                f"crosslink_copies must contain exactly 2 elements, got {len(value)}",
                error_code="CFG_ERR_006"
            )
        for item in value:
            if item not in valid_names:
                raise ConfigurationError(
                    f"Invalid translation name: {item}. Valid options: {valid_names}",
                    error_code="CFG_ERR_006"
                )
        if value[0] == value[1]:
            raise ConfigurationError(
                f"crosslink_copies must contain 2 different translations, got duplicate: {value[0]}",
                error_code="CFG_ERR_006"
            )
        
        return value
    
    @model_validator(mode="after")
    def validate_mutated_pdb_config(self) -> "ColbuilderConfig":
        """Validate mutated PDB configuration."""
        if self.mutated_pdb:
            # If mutated_pdb is provided, we need at least one additional crosslink
            has_additional = any([
                self.additional_1_type,
                self.additional_2_type
            ])
            
            if not has_additional and not self.crosslink:
                LOG.warning(
                    "mutated_pdb provided but no additional crosslinks specified. "
                    "Only optimization will be performed."
                )
            
            # Validate that additional crosslink positions are provided if types are specified
            if self.additional_1_type and not self.additional_1_combination:
                raise ConfigurationError(
                    "additional_1_combination required when additional_1_type is specified",
                    error_code="CFG_ERR_006"
                )
                
            if self.additional_2_type and not self.additional_2_combination:
                raise ConfigurationError(
                    "additional_2_combination required when additional_2_type is specified",
                    error_code="CFG_ERR_006"
                )
        
        return self
    
    @field_validator("mutated_pdb", mode="before")
    def validate_mutated_pdb_path(cls, value):
        """Validate mutated PDB path."""
        if value is not None:
            path = Path(value)
            if not path.is_absolute():
                # Will be resolved relative to working directory later
                return value
            return value
        return None

    @field_validator("files_mix", mode="before")
    def validate_files_mix(cls, value):
        """Validate files mix format."""
        return tuple(value) if value else None

    @property
    def ratio_mix(self) -> Dict[str, int]:
        """Get ratio mix property."""
        return self._ratio_mix

    @ratio_mix.setter
    def ratio_mix(self, value: Union[str, Dict[str, int]]):
        """Set ratio mix property."""
        self._ratio_mix = self._convert_ratio_mix(value)

    @field_validator("manual_replacements", mode="before")
    def validate_manual_replacements(
        cls, value: Optional[Union[str, List[str], Tuple[str, ...]]]
    ) -> Optional[List[str]]:
        """Normalize manual replacement directives."""
        if value is None:
            return None
        if isinstance(value, str):
            cleaned = [line.strip() for line in value.splitlines() if line.strip()]
            return cleaned or None
        if isinstance(value, (list, tuple)):
            cleaned = [str(v).strip() for v in value if str(v).strip()]
            return cleaned or None
        raise ConfigurationError(
            f"Invalid manual_replacements type: {type(value).__name__}",
            error_code="CFG_ERR_006",
        )

    def set_mode(self):
        """Set operation mode based on configuration flags."""
        self.mode = OperationMode.NONE
        if self.sequence_generator:
            self.mode |= OperationMode.SEQUENCE
        if self.geometry_generator:
            self.mode |= OperationMode.GEOMETRY
        if self.topology_generator:
            self.mode |= OperationMode.TOPOLOGY
        if self.mix_bool:
            self.mode |= OperationMode.MIX
        if self.replace_bool:
            self.mode |= OperationMode.REPLACE
        if self.auto_fix_unpaired:
            self.mode |= OperationMode.REPLACE

    def validate_paths(self):
        """Validate existence of required input paths and files."""
        input_paths = [
            self.PROJECT_ROOT,
            self.DATA_DIR,
            self.HOMOLOGY_LIB_DIR,
            self.TEMPLATE_PDB_PATH,
            self.TEMPLATE_FASTA_PATH,
            self.RESTYP_LIB_PATH,
            self.TOP_HEAV_LIB_PATH,
            self.PAR_MOD_LIB_PATH,
            self.CROSSLINKS_FILE,
            self.FORCE_FIELD_DIR,
            self.CHIMERA_SCRIPTS_DIR,
        ]

        # Only check paths that are specified
        optional_paths = [
            self.mutated_pdb,
            self.pdb_file,
            self.crystalcontacts_file,
            self.connect_file,
            self.replace_file,
        ]

        for path in input_paths:
            if isinstance(path, Path) and path is not None and not path.exists():
                # Try to find the file in alternative locations
                basename = os.path.basename(path)
                alternate_paths = [
                    self.working_directory / basename,  # Working directory
                    Path.cwd() / basename,  # Current directory
                    Path.cwd() / "data" / basename,  # Current directory's data folder
                ]

                found = False
                for alt_path in alternate_paths:
                    if alt_path.exists():
                        LOG.warning(
                            f"Required input found in alternate location: {alt_path}"
                        )
                        found = True
                        break

                if not found:
                    raise ConfigurationError(
                        f"Required input file or directory not found: {path}",
                        error_code="CFG_ERR_002",
                    )

        for path in optional_paths:
            if isinstance(path, Path) and path is not None and not path.exists():
                LOG.warning(f"Optional input file not found: {path}")

    @field_validator("ratio_mix", mode="before")
    def validate_ratio_mix(
        cls, value: Union[str, Dict[str, int], Tuple]
    ) -> Dict[str, int]:
        """Validate ratio mix format and values."""
        LOG.debug(f"Validating ratio_mix: {value}")
        if isinstance(value, str):
            try:
                result = {
                    item.split(":")[0]: int(item.split(":")[1])
                    for item in value.split()
                }
                if sum(result.values()) != 100:
                    raise ConfigurationError(
                        "Mix ratios must sum to 100%", error_code="CFG_ERR_004"
                    )
                return result
            except (ValueError, IndexError):
                raise ConfigurationError(
                    "Invalid ratio_mix format. Expected 'Type:percentage Type:percentage'",
                    error_code="CFG_ERR_004",
                )
        elif isinstance(value, dict):
            if sum(value.values()) != 100:
                raise ConfigurationError(
                    "Mix ratios must sum to 100%", error_code="CFG_ERR_004"
                )
            return value
        elif isinstance(value, tuple):
            LOG.info(f"ratio_mix is a tuple: {value}")
            try:
                result = {item.split(":")[0]: int(item.split(":")[1]) for item in value}
                LOG.info(f"Converted tuple to dictionary: {result}")
                if sum(result.values()) != 100:
                    raise ConfigurationError(
                        "Mix ratios must sum to 100%", error_code="CFG_ERR_004"
                    )
                return result
            except (ValueError, IndexError):
                raise ConfigurationError(
                    "Invalid ratio_mix format in tuple. Expected 'Type:percentage Type:percentage'",
                    error_code="CFG_ERR_004",
                )
        else:
            raise ConfigurationError(
                f"Invalid ratio_mix type. Expected string, dictionary, or tuple, got {type(value).__name__}",
                error_code="CFG_ERR_004",
            )

    @model_validator(mode="after")
    def validate_mix_config(cls, values):
        """Validate mixing configuration including fibril_length."""
        if values.mix_bool:
            if values.ratio_mix is None:
                raise ConfigurationError(
                    "ratio_mix is required when mix_bool is True",
                    error_code="CFG_ERR_004",
                )
            if not values.files_mix:
                raise ConfigurationError(
                    "files_mix is required when mix_bool is True",
                    error_code="CFG_ERR_004",
                )
            if values.fibril_length is None:
                values.fibril_length = values.get(
                    "fibril_length", 40.0
                )  # Default to 40 for mixing
            LOG.debug(f"Validated fibril_length for mixing: {values.fibril_length}")
        return values

    @model_validator(mode="after")
    def validate_replace_config(self) -> "ColbuilderConfig":
        """Validate replacement configuration."""
        if self.replace_bool:
            has_manual_replacements = bool(self.manual_replacements)
            if self.geometry_generator:
                if self.ratio_replace is None and not has_manual_replacements:
                    raise ConfigurationError(
                        "Either ratio_replace or manual_replacements must be specified when replace_bool is True with geometry generation",
                        error_code="CFG_ERR_006",
                    )
            elif self.ratio_replace is None and self.replace_file is None and not has_manual_replacements:
                raise ConfigurationError(
                    "Either ratio_replace, replace_file, or manual_replacements must be specified when replace_bool is True",
                    error_code="CFG_ERR_006",
                )
            if self.ratio_replace_scope not in {"enzymatic", "non_enzymatic", "all"}:
                raise ConfigurationError(
                    "ratio_replace_scope must be one of: enzymatic, non_enzymatic, all",
                    error_code="CFG_ERR_006",
                )
            if self.ratio_replace is not None:
                if not isinstance(self.ratio_replace, (int, float)):
                    raise ConfigurationError(
                        "ratio_replace must be a numeric value",
                        error_code="CFG_ERR_006",
                    )
                if self.ratio_replace < 0 or self.ratio_replace > 100:
                    raise ConfigurationError(
                        f"ratio_replace must be between 0 and 100, got {self.ratio_replace}",
                        error_code="CFG_ERR_006",
                    )
            if self.replace_file is not None and not isinstance(
                self.replace_file, Path
            ):
                self.replace_file = Path(self.replace_file)

            # Ensure the CHIMERA_SCRIPTS_DIR exists and contains swapaa.py
            chimera_scripts_path = self.CHIMERA_SCRIPTS_DIR
            if not chimera_scripts_path.exists():
                raise ConfigurationError(
                    f"CHIMERA_SCRIPTS_DIR not found: {chimera_scripts_path}",
                    error_code="CFG_ERR_002",
                )

            swapaa_script = chimera_scripts_path / "swapaa.py"
            if not swapaa_script.exists():
                LOG.warning(
                    f"swapaa.py script not found in {chimera_scripts_path}. "
                    f"This script is required for crosslink replacement."
                )

        return self

    def get_conda_env_path(self) -> Optional[str]:
        """
        Get the full path to the conda environment for martinize2.

        Returns:
            Optional[str]: Full path to the conda environment or None if not found
        """
        if not self.martinize2_env:
            return None

        try:
            # Run conda info command to get environment paths
            result = subprocess.run(
                ["conda", "info", "--envs", "--json"],
                capture_output=True,
                text=True,
                check=True,
            )
            env_data = json.loads(result.stdout)

            # Look for the specified environment in the paths
            for env_path in env_data["envs"]:
                if os.path.basename(env_path) == self.martinize2_env:
                    LOG.info(f"Found conda environment path: {env_path}")
                    return env_path

            LOG.warning(
                f"Could not find conda environment path for: {self.martinize2_env}"
            )
            return None
        except Exception as e:
            LOG.error(f"Error determining conda environment path: {e}")
            return None

    def update(self, new_config: Dict[str, Any]):
        """Update configuration with new values."""
        for key, value in new_config.items():
            if hasattr(self, key):
                setattr(self, key, value)
        self.set_mode()

    def __str__(self):
        """String representation of the configuration."""
        return f"ColbuilderConfig(mode={self.mode}, pdb_file={self.pdb_file}, output={self.output})"

    @property
    def output(self) -> str:
        """Get the output filename."""
        return f"collagen_fibril_{self.species}"

    class Config:
        """Pydantic configuration."""

        use_enum_values = True


def validate_input_files(config: ColbuilderConfig) -> None:
    """Validate input files if they are provided and required."""
    validator = BioformatValidator()

    try:
        if config.sequence_generator and config.fasta_file is not None:
            fasta_path = Path(config.fasta_file)
            warnings = validator.validate_input_files(fasta_path, None)
            if warnings:
                for warning in warnings:
                    LOG.warning(warning)

        if config.geometry_generator and config.pdb_file is not None:
            pdb_path = config.pdb_file
            warnings = validator.validate_input_files(
                fasta_path=Path(config.fasta_file) if config.fasta_file else None,
                pdb_path=pdb_path,
            )
            if warnings:
                for warning in warnings:
                    LOG.warning(warning)

    except SequenceGenerationError as e:
        raise ConfigurationError(
            f"ERROR: {str(e)}", original_error=e, error_code="CFG_ERR_005"
        )
    except GeometryGenerationError as e:
        raise ConfigurationError(
            f"ERROR: {str(e)}", original_error=e, error_code="CFG_ERR_005"
        )
    except Exception as e:
        raise ConfigurationError(
            "File validation failed", original_error=e, error_code="CFG_ERR_005"
        )


_config_instance = None


def get_config(
    existing_config: Optional[ColbuilderConfig] = None, **kwargs
) -> ColbuilderConfig:
    """Get or create the configuration singleton."""
    global _config_instance
    if (_config_instance is None) or (existing_config is not None):
        if existing_config:
            LOG.debug("Reusing existing ColbuilderConfig object in get_config")
            _config_instance = existing_config
        else:
            LOG.debug("Creating ColbuilderConfig object in get_config")
            _config_instance = ColbuilderConfig(**kwargs)
    return _config_instance


def load_yaml_config(yaml_file: Path) -> Dict[str, Any]:
    """Load configuration from YAML file."""
    try:
        with open(yaml_file, "r") as file:
            config_data = yaml.safe_load(file)
            LOG.debug(f"Loaded YAML configuration: {config_data}")
            return config_data
    except Exception as e:
        raise ConfigurationError(
            f"ERROR: {yaml_file}", original_error=e, error_code="CFG_ERR_001"
        )


def validate_config(config: Dict[str, Any]) -> ColbuilderConfig:
    """Validate configuration dictionary and create config object."""
    LOG.debug("Calling validate_config")
    try:
        return ColbuilderConfig(**config)
    except Exception as e:
        LOG.error(f"Configuration validation failed: {e}")
        raise ConfigurationError(
            "Configuration validation failed",
            original_error=e,
            error_code="CFG_ERR_001",
        )
