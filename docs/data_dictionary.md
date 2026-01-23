# ColBuilder Data Dictionary

## Purpose

This data dictionary provides a comprehensive reference for parameters, variables, and data structures used in ColBuilder. It serves several purposes:

- **Configuration Reference**: Helps users set up correct configuration files with appropriate parameter values
- **Error Troubleshooting**: Assists in understanding error messages related to specific parameters
- **Code Understanding**: Supports developers working with the ColBuilder codebase
- **Documentation**: Provides detailed explanations of data structures and their relationships

## Table of Contents

- [ColBuilder Data Dictionary](#colbuilder-data-dictionary)
  - [Purpose](#purpose)
  - [Table of Contents](#table-of-contents)
  - [Quick Reference](#quick-reference)
  - [Configuration Parameters](#configuration-parameters)
    - [General Options](#general-options)
    - [Sequence Generation Parameters](#sequence-generation-parameters)
      - [Crosslink Configuration](#crosslink-configuration)
    - [Geometry Generation Parameters](#geometry-generation-parameters)
    - [Mixing and Replacement Parameters](#mixing-and-replacement-parameters)
    - [Topology Generation Parameters](#topology-generation-parameters)
  - [Internal Data Structures](#internal-data-structures)
    - [Sequence Generation Structures](#sequence-generation-structures)
    - [Geometry Generation Structures](#geometry-generation-structures)
    - [Topology Generation Structures](#topology-generation-structures)
  - [Glossary](#glossary)

## Quick Reference

The most commonly used parameters for ColBuilder configuration:

| Parameter | Description | Example Value | Required? |
|-----------|-------------|---------------|-----------|
| species | Species name for collagen sequence | "homo_sapiens" | Yes |
| sequence_generator | Enable sequence generation | true | - |
| geometry_generator | Enable geometry generation | true | - |
| topology_generator | Enable topology generation | false | - |
| crosslink | Enable crosslinking | true | - |
| n_term_type | N-terminal crosslink type | "HLKNL" | If crosslink=true |
| c_term_type | C-terminal crosslink type | "HLKNL" | If crosslink=true |
| n_term_combination | string | N-terminal residue combination | "9.C - 947.A" | If crosslink=true |
| c_term_combination | string | C-terminal residue combination | "1047.C - 104.C" | If crosslink=true |
| mutated_pdb | Pre-mutated PDB for additional crosslinks | "rattusnorvegicus_N_PYD_C_PYD.pdb" | Optional |
| additional_1_type | First non-enzymatic crosslink | "Glucosepane" | Optional |
| additional_1_combination | First non-enzymatic combination | "523.A - 286.C" | If additional_1_type |
| additional_2_type | Second non-enzymatic crosslink | "Pentosidine" | Optional |
| additional_2_combination | Second non-enzymatic combination | "523.A - 286.C" | If additional_2_type |
| crosslink_copies | Translation pair for optimization | ["D0", "D5"] | Optional |
| fibril_length | Length of microfibril (nm) | 60.0 | For geometry |
| contact_distance | Contact distance (Å) | 20 | For geometry |
| force_field | Force field for topology | "amber99" | For topology |

## Configuration Parameters

### General Options

| Parameter | Type | Description | Valid Values | Default |
|-----------|------|-------------|--------------|---------|
| debug | boolean | Enable detailed debug logging | true/false | false |
| working_directory | Path | Working directory for input/output files | Any valid path | Current directory |
| config_file | Path | Path to configuration YAML file | Any valid path | None |
| mode | OperationMode | Computed operation mode flags | Combination of SEQUENCE, GEOMETRY, TOPOLOGY, MIX, REPLACE | Computed from other flags |

### Sequence Generation Parameters

| Parameter | Type | Description | Valid Values | Default |
|-----------|------|-------------|--------------|---------|
| sequence_generator | boolean | Enable sequence generation | true/false | false |
| species | string | Species name for sequence and crosslinks | One of supported species* or custom name | Required |
| fasta_file | Path | Path to input FASTA file | Valid FASTA file path or null | Auto-generated based on species |
| crosslink | boolean | Enable crosslinking in the model | true/false | false |
| mutated_pdb | Path | Input PDB for additional crosslinks | Valid PDB file path or null | None |
| additional_1_type | string | First non-enzymatic crosslink type | e.g., "Glucosepane" | None |
| additional_1_combination | string | Residue combination for additional_1 | Valid combination | None |
| additional_2_type | string | Second non-enzymatic crosslink type | e.g., "Pentosidine" | None |
| additional_2_combination | string | Residue combination for additional_2 | Valid combination | None |
| crosslink_copies | list | Translation pair for crosslink optimization | ["D0", "D5"] | ["D0", "D5"] |

*Supported species: homo_sapiens, pan_troglodytes, pongo_abelii, callithrix_jacchus, otolemur_garnettii, mus_musculus, rattus_norvegicus, bos_taurus, canis_lupus, ailuropoda_melanoleuca, mustela_putorius, myotis_lucifugus, loxodonta_africana, danio_rerio, oreochromis_niloticus, oryzias_latipes, tetraodon_nigroviridis, xiphophorus_maculatus, pelodiscus_sinensis

#### Crosslink Configuration

| Parameter | Type | Description | Valid Values | Example |
|-----------|------|-------------|--------------|---------|
| n_term_type | string | N-terminal crosslink type | "DPD", "DPL", "HLKNL", "LKNL", "NOCROSS", "PYD", "PYL", "deHHLNL", "deHLNL", "MOLD" | "HLKNL" |
| c_term_type | string | C-terminal crosslink type | Same as n_term_type | "HLKNL" |
| n_term_combination | string | N-terminal residue combination | Format: "ResNum.Chain - ResNum.Chain" | "9.C - 947.A" |
| c_term_combination | string | C-terminal residue combination | Format: "ResNum.Chain - ResNum.Chain" | "1047.C - 104.C" |

*Check available crosslinks and respective combinations at [src/colbuilder/data/sequence/crosslinks.csv](https://github.com/graeter-group/colbuilder/blob/main/src/colbuilder/data/sequence/crosslinks.csv)

**Shift Note**:
- For the additional non-enzymatic step (`mutated_pdb`), set `crosslink_copies`
  to the shift listed in `crosslinks.csv` for the selected crosslink/position.
  Enzymatic crosslinks default to D0-D5 unless specified otherwise.

**Validation Rules**:
- For human (homo_sapiens) HLKNL crosslinks, valid N-terminal combinations include: "5.B - 944.B", "9.C - 944.B", "9.C - 947.A", "947.A - 5.B"
- For human (homo_sapiens) HLKNL crosslinks, valid C-terminal combinations include: "104.C - 1047.A", "1047.A - 98.B", "1047.C - 104.C", "1047.C - 98.B"
- Combinations vary by species and crosslink type
- Format validation: must match the pattern "^\d+\.[A-C]\s*-\s*\d+\.[A-C]$"

### Geometry Generation Parameters

| Parameter | Type | Description | Valid Values | Default |
|-----------|------|-------------|--------------|---------|
| geometry_generator | boolean | Enable geometry generation | true/false | false |
| pdb_file | Path | Input PDB file (if sequence_generator=false) | Valid PDB file path | None |
| contact_distance | float | Contact distance for microfibril (Å) | Positive number (typically 15-40) | None (required) |
| fibril_length | float | Length of microfibril (nm) | 0 < value ≤ 334 | None (required)  |
| crystalcontacts_file | Path | File with crystal contacts | Valid file path | None |
| connect_file | Path | File with connection information | Valid file path | None |
| crystalcontacts_optimize | boolean | Optimize crystal contacts | true/false | true |
| solution_space | List/Tuple | Solution space dimensions [dx, dy, dz] | Three positive numbers | [1, 1, 1] |
| pdb_first_line | string | Crystal contacts information | Valid PDB CRYST1 line | Default crystal parameters |

**Note**: Either contact_distance or crystalcontacts_file must be provided when geometry_generator is true.

### Mixing and Replacement Parameters

| Parameter | Type | Description | Valid Values | Default |
|-----------|------|-------------|--------------|---------|
| mix_bool | boolean | Generate mixed crosslinked microfibril | true/false | false |
| ratio_mix | Dict or string | Ratio for mix-crosslink setup | Format: "Type:percentage Type:percentage" or Dict[str, int] | {} |
| files_mix | List of Paths | PDB files with different crosslink types | Valid PDB file paths | [] |
| replace_bool | boolean | Replace crosslinks with lysines | true/false | false |
| ratio_replace | float | Percentage of crosslinks to replace | 0-100 | None |
| ratio_replace_scope | string | Scope for replacement | "enzymatic", "non_enzymatic", "all" | "non_enzymatic" |
| replace_file | Path | Input PDB file of fibril with crosslinks to be replaced | Valid file path | None |
| manual_replacements | List[str] | Manual replacement directives | Valid directives | None |
| auto_fix_unpaired | boolean | Auto-detect unpaired enzymatic crosslinks | true/false | false |

**Validation Rules**:
- When mix_bool=true, ratio_mix and files_mix must be provided
- Percentages in ratio_mix must sum to 100
- When replace_bool=true, either geometry_generation=true or replace_file must be provided
- ratio_replace must be between 0 and 100

### Topology Generation Parameters

| Parameter | Type | Description | Valid Values | Default |
|-----------|------|-------------|--------------|---------|
| topology_generator | boolean | Generate topology files | true/false | false |
| force_field | string | Force field for simulations | "amber99", "martini3" | None |
| martinize2_command | string | Detected Martinize2 command path | Valid executable path | Auto-detected |
| martinize2_env | string | Detected Martinize2 conda environment | Valid conda environment name | Auto-detected |
| use_conda_run | boolean | Use conda run for Martinize2 | true/false | false |
| go_epsilon | float | GO epsilon for Martini3-CG parametrization | Positive number | 9.414 |

## Internal Data Structures

### Sequence Generation Structures

| Name | Type | Description |
|------|------|-------------|
| CrosslinkPair | class | Represents a pair of residues involved in a crosslink |
| ResiduePosition | class | Represents a residue position with residue type, atom name, and position string |
| SequenceGenerator | class | Manages the generation of collagen structure from sequences |
| OptimizationState | class | Tracks state during crosslink optimization |

**Key Methods in SequenceGenerator**:
- `generate()`: Main method to generate collagen structure
- `_run_alignment()`: Performs sequence alignment
- `_run_modelling()`: Runs MODELLER for structure generation
- `_load_crosslinks()`: Loads crosslink information
- `_apply_crosslinks()`: Applies crosslinks to structure
- `_finalize_output()`: Optimizes and finalizes output

**Files and Paths**:
- `TEMPLATE_PDB_PATH`: Path to template PDB file
- `TEMPLATE_FASTA_PATH`: Path to template FASTA file
- `RESTYP_LIB_PATH`: Path to MODELLER residue type library
- `TOP_HEAV_LIB_PATH`: Path to MODELLER topology library
- `PAR_MOD_LIB_PATH`: Path to MODELLER parameter library
- `CROSSLINKS_FILE`: Path to crosslinks CSV file
- `CHIMERA_SCRIPTS_DIR`: Path to UCSF Chimera scripts directory

### Geometry Generation Structures

| Name | Type | Description |
|------|------|-------------|
| System | class | Represents a molecular system with models, coordinates, and topology |
| Crystal | class | Represents the crystal structure of a collagen microfibril |
| GeometryService | class | Coordinates geometry operations |
| CrystalBuilder | class | Builds the initial crystal structure |
| CrosslinkMixer | class | Mixes different crosslink types |
| CrosslinkReplacer | class | Replaces crosslinks with standard amino acids |

**Key Methods in GeometryService**:
- `build_geometry()`: Main entry point for geometry generation
- `_handle_direct_replacement()`: Handles replacement without geometry generation
- `_handle_mixing_only()`: Handles mixing without geometry generation
- `_handle_full_generation()`: Handles complete geometry generation

**System Properties**:
- `crystal`: Crystal structure information
- `coordinates`: Atomic coordinates
- `models`: Dictionary of model objects
- `connect`: Connections between models
- `total_atoms`: Total number of atoms in system

**Temporary Resources**:
- Standard temp files: "replace.txt"
- Standard temp directories: "NC", "NCP", "D", "T"

### Topology Generation Structures

| Name | Type | Description |
|------|------|-------------|
| Amber | class | Handles amber99 force field topology generation |
| Martini | class | Handles martini3 force field topology generation |

**Key Functions**:
- `build_topology()`: Main entry point for topology generation
- `build_amber99()`: Builds amber99 topology
- `build_martini3()`: Builds martini3 topology

**Output Files**:
- `{species}_topology_files/`: Directory containing topology files
- `collagen_fibril_{species}.top`: Topology file
- `collagen_fibril_{species}.gro`: GROMACS coordinate file
- Various .itp files: Include topology files

**Temporary Files**:
- Force field files (amber99sb-star-ildnp.ff/ or similar)
- Various intermediate files (*.itp, *.CG.pdb, *.merge.pdb, etc.)

## Glossary

- **Collagen**: Main structural protein in connective tissues
- **Triple helix**: Three polypeptide chains wound around each other, forming the basic collagen structure
- **Microfibril**: Assembly of multiple collagen molecules in a quasi-hexagonal packing
- **Crosslink**: Covalent bond connecting collagen chains, providing structural stability
- **Homology modeling**: Method to predict protein structure based on related proteins
- **FASTA**: Text-based format for representing peptide sequences
- **MSA**: Multiple Sequence Alignment
- **MODELLER**: Software for homology modeling of protein structures
- **PDB**: Protein Data Bank file format
- **GROMACS**: Molecular dynamics simulation software
- **Force field**: Parameters and functions describing potential energy in simulations
- **Crystal contacts**: Contact points between molecules in a crystal structure
- **Bravais lattice**: Infinite array of discrete points with translational symmetry
- **Simulated annealing**: Optimization technique for crosslink positioning
