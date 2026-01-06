# ColBuilder Configuration Reference

## Table of Contents
- [ColBuilder Configuration Reference](#colbuilder-configuration-reference)
  - [Table of Contents](#table-of-contents)
  - [Introduction](#introduction)
  - [Configuration File Format](#configuration-file-format)
  - [Operation Mode Parameters](#operation-mode-parameters)
  - [Input Configuration Parameters](#input-configuration-parameters)
  - [Sequence Generation Parameters](#sequence-generation-parameters)
  - [Geometry Generation Parameters](#geometry-generation-parameters)
  - [Mixing and Replacement Parameters](#mixing-and-replacement-parameters)
  - [Topology Generation Parameters](#topology-generation-parameters)
  - [Common Configuration Examples](#common-configuration-examples)
    - [Basic Human Collagen Microfibril](#basic-human-collagen-microfibril)
    - [Bovine Collagen with Trivalent Crosslinks](#bovine-collagen-with-trivalent-crosslinks)
    - [Mixed Crosslinked Microfibril (80% Divalent + 20% Trivalent crosslink)](#mixed-crosslinked-microfibril-80-divalent--20-trivalent-crosslink)
    - [Microfibril and Coarse-Grained Topology Generation](#microfibril-and-coarse-grained-topology-generation)
  - [Parameter Interactions and Dependencies](#parameter-interactions-and-dependencies)
  - [Configuration Validation](#configuration-validation)

## Introduction

ColBuilder uses a YAML configuration file to control all aspects of the pipeline. This document provides a comprehensive reference for all available configuration parameters, their meanings, and how they interact with each other.

## Configuration File Format

Configuration files for ColBuilder are written in YAML format. A basic configuration file might look like this:

```yaml
# Basic configuration
mode: null
config_file: null
sequence_generator: true
geometry_generator: true
topology_generator: false
debug: false
working_directory: "./"

# Input Configuration
species: "homo_sapiens"

# Sequence Settings
fasta_file: null
mutated_pdb: null
crosslink: true
n_term_type: "HLKNL"
c_term_type: "HLKNL"
n_term_combination: "9.C - 947.A"
c_term_combination: "1047.C - 104.C"
additional_1_type: null
additional_1_combination: null
additional_2_type: null
additional_2_combination: null
crosslink_copies: ["D0", "D5"]

# Geometry Parameters
pdb_file: null
contact_distance: 20
fibril_length: 60.0
crystalcontacts_file: null
connect_file: null
crystalcontacts_optimize: false

# Mixing Options
mix_bool: false
ratio_mix: "D:70 T:30"
files_mix:
 - "human-D.pdb"
 - "human-T.pdb"

# Replacement Options
replace_bool: false
ratio_replace: 30
replace_file: null

# Topology Options
force_field: "amber99"
```

## Operation Mode Parameters

These parameters control which stages of the pipeline are executed.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `mode` | string, null | null | Specific operation mode if needed (deprecated, use individual flags instead) |
| `config_file` | string, null | null | Path to another config file (for nested configs) |
| `sequence_generator` | boolean | true | Enable/disable sequence generation stage |
| `geometry_generator` | boolean | true | Enable/disable geometry generation stage |
| `topology_generator` | boolean | false | Enable/disable topology generation stage |
| `debug` | boolean | false | Enable detailed debugging output |
| `working_directory` | string | "./" | Directory for input and output files |

**Notes**:
- The `sequence_generator`, `geometry_generator`, and `topology_generator` flags determine which pipeline stages are executed.
- If all three are set to `true`, the full pipeline will run in sequence.
- The `working_directory` parameter sets the base directory for all input and output files.

## Input Configuration Parameters

These parameters control basic input settings, particularly species selection.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `species` | string | "rattus_norvegicus" | Species for collagen sequence and structure |

**Available Species Options**:
- **Mammals (Primates)**: homo_sapiens, pan_troglodytes, pongo_abelii, callithrix_jacchus, otolemur_garnettii
- **Mammals (Rodents)**: mus_musculus, rattus_norvegicus
- **Mammals (Other)**: bos_taurus, canis_lupus, ailuropoda_melanoleuca, mustela_putorius, myotis_lucifugus, loxodonta_africana
- **Fish**: danio_rerio, oreochromis_niloticus, oryzias_latipes, tetraodon_nigroviridis, xiphophorus_maculatus
- **Reptiles**: pelodiscus_sinensis

**Notes**:
- The `species` parameter is used to select the appropriate sequence and crosslink information.
- Rat collagen (rattus_norvegicus) is used as the default structural template for all species.

## Sequence Generation Parameters

These parameters control the sequence generation stage (homology modeling).

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `fasta_file` | string, null | null | Path to custom FASTA file (if null, generated based on species) |
| `mutated_pdb` | string, null | null | Use a pre-mutated PDB as input for additional crosslinks |
| `crosslink` | boolean | true | Enable/disable crosslinking in the model |
| `n_term_type` | string | "HLKNL" | N-terminal crosslink type |
| `c_term_type` | string | "HLKNL" | C-terminal crosslink type |
| `n_term_combination` | string | depends on species | N-terminal residue combination for crosslinks |
| `c_term_combination` | string | depends on species | C-terminal residue combination for crosslinks |
| `additional_1_type` | string, null | null | Optional non-enzymatic crosslink (AGE) |
| `additional_1_combination` | string, null | null | Residue combination for the first additional crosslink |
| `additional_2_type` | string, null | null | Optional second non-enzymatic crosslink |
| `additional_2_combination` | string, null | null | Residue combination for the second additional crosslink |
| `crosslink_copies` | list | ["D0", "D5"] | Translation pair for crosslink optimization |

**Available Crosslink Types**:
- **DPD**: Deoxypyridinoline (trivalent crosslink)
- **DPL**: Deoxypyrrololine (trivalent)
- **HLKNL**: Hydroxylysino-5-ketonorleucine (divalent)
- **LKNL**: Lysino-5-ketonorleucine (divalent)
- **NOCROSS**: No crosslinking
- **PYD**: Pyridinoline (trivalent)
- **PYL**: Pyrrole (trivalent)
- **deHHLNL**: Dehydro-hydroxylysino-norleucine (divalent)
- **deHLNL**: Dehydro-lysino-norleucine (divalent)
- **Glucosepane**: Advanced glycation end-product crosslink (specific species)
- **Pentosidine**: Advanced glycation end-product crosslink (specific species)
- **MOLD**: Non-enzymatic, LYS-LYS derived crosslink

**Residue Combination Examples for Homo Sapiens with HLKNL**:
- **N-terminal combinations**: "5.B - 944.B", "9.C - 944.B", "9.C - 947.A", "947.A - 5.B"
- **C-terminal combinations**: "104.C - 1047.A", "1047.A - 98.B", "1047.C - 104.C", "1047.C - 98.B"

**Notes**:
- If `fasta_file` is null, ColBuilder will use built-in sequence data for the specified species.
- The `crosslink` parameter must be set to `true` for crosslinking to be applied.
- The residue combinations specify which residues are involved in crosslinking and must match valid combinations for the chosen species. A complete list of the species and combinations currently available can be found at [src/colbuilder/data/sequence/crosslinks.csv](https://github.com/graeter-group/colbuilder/blob/main/src/colbuilder/data/sequence/crosslinks.csv).
- Residue combinations follow the format: `[Residue Number].[Chain] - [Residue Number].[Chain]` or `[Residue Number].[Chain] - [Residue Number].[Chain] - [Residue Number].[Chain]` for trivalent crosslinks.
- When adding non-enzymatic crosslinks via `mutated_pdb`, use the crosslink-specific
  shift listed in `crosslinks.csv` for `crosslink_copies`. Enzymatic crosslinks
  use D0-D5 by default unless specified otherwise.

## Geometry Generation Parameters

These parameters control the generation of the microfibril structure.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `pdb_file` | string, null | null | Input PDB file (set to null if sequence_generator is true) |
| `contact_distance` | integer | 20 | Distance threshold for contacts (Å) |
| `fibril_length` | float | 60.0 | Length of the generated fibril (nm) |
| `crystalcontacts_file` | string, null | null | File with crystal contacts information |
| `connect_file` | string, null | null | File with connection information |
| `crystalcontacts_optimize` | boolean | false | Optimize crystal contacts during generation |

**Notes**:
- The `contact_distance` parameter controls the radial size of the microfibril. Larger values result in thicker fibrils but require more computation time and memory.
- The `fibril_length` parameter sets the length of the fibril in nanometers. Larger values create longer fibrils.
- If `crystalcontacts_optimize` is set to true, ColBuilder will attempt to optimize the packing of collagen molecules in a Bravais lattice.
- If `pdb_file` is provided, the sequence generation stage will be skipped and this PDB will be used as input for geometry generation.

## Mixing and Replacement Parameters

These parameters control advanced features for creating mixed crosslinked microfibrils or replacing crosslinks by standard residues.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `mix_bool` | boolean | false | Enable mixing of different crosslink types |
| `ratio_mix` | string | "D:70 T:30" | Format: "Type:percentage Type:percentage" |
| `files_mix` | list of strings | | Required if mix_bool is true, paths to PDB files with different crosslink types |
| `replace_bool` | boolean | false | Enable crosslink replacement (with lysines) |
| `ratio_replace` | integer | 30 | Percentage of crosslinks to replace |
| `ratio_replace_scope` | string | "non_enzymatic" | Scope for replacement: enzymatic, non_enzymatic, all |
| `replace_file` | string, null | null | File with crosslinks to be replaced |
| `manual_replacements` | list of strings | null | Manual replacement directives |
| `auto_fix_unpaired` | boolean | false | Auto-detect unpaired enzymatic crosslinks and replace |

**Notes**:
- The `mix_bool` feature allows creation of heterogeneous crosslinked microfibrils, which more closely resemble natural collagen.
- The `ratio_mix` parameter specifies the proportion of each crosslink type in the mixed microfibril.
- The `files_mix` parameter specifies the path to the PDB files of two collagen molecules, each with a different crosslink type.
- The `replace_bool` feature simulates partial crosslinking by replacing some crosslinks with unmodified lysine residues.
- The `ratio_replace` parameter controls what percentage of crosslinks should be replaced.
- The `replace_file` parameter specifies the path to the PDB file of a previously generated collagen microfibril.
- Topology generation can be enabled in the same run after mixing.

## Topology Generation Parameters

These parameters control the generation of topology files for molecular dynamics simulations.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `force_field` | string | "amber99" | Force field for topology generation |

**Available Force Field Options**:
- **amber99**: Standard Amber99 force field (add publication) 
- **martini3**: Martini 3 coarse-grained force field (add publication)

**Notes**:
- The `force_field` parameter selects which force field to use for generating topology files.
- The amber99 force field is recommended for most atomistic simulations of collagen.
- Custom force field parameters for collagen and crosslinks are included in ColBuilder.

## Common Configuration Examples

### Basic Human Collagen Microfibril

```yaml
species: "homo_sapiens"
sequence_generator: true
geometry_generator: true
crosslink: true
fibril_length: 60.0
contact_distance: 20
n_term_type: "HLKNL"
c_term_type: "HLKNL"
n_term_combination: "9.C - 947.A"
c_term_combination: "1047.C - 104.C"
```

### Bovine Collagen with Trivalent Crosslinks

```yaml
species: "bos_taurus"
sequence_generator: true
geometry_generator: true
crosslink: true
n_term_type: "PYD"
c_term_type: "PYD"
n_term_combination: "9.C	-	5.B	-	942.B"
c_term_combination: "1046.C	- 1046.A	- 103.C" 
fibril_length: 80.0
contact_distance: 20
```

### Mixed Crosslinked Microfibril (80% Divalent + 20% Trivalent crosslink)

```yaml
species: "homo_sapiens"
sequence_generator: false
geometry_generator: false
mix_bool: true
ratio_mix: "D:80 T:20"
files_mix:
 - "human-D.pdb"
 - "human-T.pdb"
```

### Microfibril and Coarse-Grained Topology Generation

```yaml
species: "homo_sapiens"
sequence_generator: false
geometry_generator: true
topology_generator: true
pdb_file: "path/to/collagen_molecule.pdb"
fibril_length: 70.0
contact_distance: 40
force_field: "martini3"
```

## Parameter Interactions and Dependencies

Understanding how parameters interact is important for successful use of ColBuilder:

1. **Pipeline Stage Dependencies**:
   - If `sequence_generator` is true, a new PDB with the structure of a collagen molecule will be generated, and any provided `pdb_file` will be ignored.
   - If `sequence_generator` is false but `geometry_generator` is true, you must provide a valid `pdb_file`.

2. **Crosslinking Dependencies**:
   - If `crosslink` is true, you must specify valid `n_term_type`, `c_term_type`, `n_term_combination`, and `c_term_combination` values.
   - The crosslink type and residue combinations must be compatible with the chosen species.

3. **Mixing and Replacement Dependencies**:
   - If `mix_bool` is true, you must provide valid files in `files_mix` and a proper ratio in `ratio_mix`.
   - If `replace_bool` is true, you must specify a valid percentage in `ratio_replace`. If you are not also generating a geometry (i.e., `generate_geometry` is false), you must also provide the path to a PDB file of a previously generated collagen fibril in `replace_file`.

4. **Geometry Generation Dependencies**:
   - If `crystalcontacts_optimize` is true, the geometry generation will take longer but may produce better-packed microfibrils.
   - The `contact_distance` parameter becomes irrelevant if `crystalcontacts_file` is provided.

## Configuration Validation

ColBuilder performs validation of the configuration file before running the pipeline:

1. **Required Parameters**: Checks that all required parameters are present
2. **Parameter Types**: Validates that parameters have the correct data types
3. **Valid Values**: Ensures parameters have valid values (e.g., valid species name, crosslink types)
4. **Consistency**: Verifies that parameter combinations are consistent and compatible
5. **File Existence**: Checks that any specified input files exist and are readable

If validation fails, ColBuilder will provide an error message indicating the specific issue with the configuration file.

For any issues with configuration, refer to the error message or enable debug mode (`debug: true`) for more detailed information about the validation process.
