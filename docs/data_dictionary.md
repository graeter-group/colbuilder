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
      - [Additional Crosslinks (Mutated PDB Workflow)](#additional-crosslinks-mutated-pdb-workflow)
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
| mutated_pdb | Pre-mutated PDB for additional crosslinks | "structure.pdb" or null | For mutated PDB workflow |
| crosslink | Enable crosslinking | true | - |
| n_term_type | N-terminal crosslink type | "HLKNL" | If crosslink=true |
| c_term_type | C-terminal crosslink type | "HLKNL" | If crosslink=true |
| n_term_combination | N-terminal residue combination | "9.C - 947.A" | If crosslink=true |
| c_term_combination | C-terminal residue combination | "1047.C - 104.C" | If crosslink=true |
| additional_1_type | First additional crosslink type | "Glucosepane" or null | For mutated PDB workflow |
| additional_1_combination | Position for first additional crosslink | "1008.A - 767.B" | If additional_1_type set |
| crosslink_copies | Periodic copies for optimization | ["D0", "D5"] | Optional |
| fibril_length | Length of microfibril (nm) | 60.0 | For geometry |
| contact_distance | Contact distance (Å) | 20 | For geometry |
| force_field | Force field for topology | "amber99" | For topology |

## Configuration Parameters

### General Options

| Parameter | Type | Description | Valid Values | Default |
|-----------|------|-------------|--------------|---------|
| debug | boolean | Keep intermediate files for debugging | true/false | false |
| working_directory | string/Path | Working directory for input/output files | Any valid path | current directory |
| config_file | Path | Path to configuration YAML file | Any valid path | None |

**Notes**:
- `working_directory` sets the base directory for all input and output files
- `debug: true` will preserve all intermediate files for troubleshooting

### Sequence Generation Parameters

| Parameter | Type | Description | Valid Values | Default |
|-----------|------|-------------|--------------|---------|
| sequence_generator | boolean | Enable sequence generation | true/false | false |
| species | string | Species name for sequence and crosslinks | One of supported species* | Required |
| mutated_pdb | string/Path/null | Pre-mutated PDB for adding additional crosslinks | Valid PDB file path or null | null |
| fasta_file | Path/null | Path to input FASTA file | Valid FASTA file path or null | Auto-generated based on species |
| crosslink | boolean | Enable crosslinking in the model | true/false | false |

*Supported (built-in) species: homo_sapiens, pan_troglodytes, pongo_abelii, callithrix_jacchus, otolemur_garnettii, mus_musculus, rattus_norvegicus, bos_taurus, canis_lupus, ailuropoda_melanoleuca, mustela_putorius, myotis_lucifugus, loxodonta_africana, danio_rerio, oreochromis_niloticus, oryzias_latipes, tetraodon_nigroviridis, xiphophorus_maculatus, pelodiscus_sinensis

**Notes**:
- **Mutated PDB workflow**: When `mutated_pdb` is provided, ColBuilder adds crosslinks to the existing structure instead of generating from sequence
- This is particularly useful for adding AGE crosslinks on top of existing enzymatic crosslinks
- The mutated PDB must have compatible terminal crosslinks specified in `n_term_*` and `c_term_*`
- **Custom species**: `species` is not restricted to the built-in list above — any name works as long as `fasta_file` is also provided (ColBuilder raises `CFG_ERR_003` otherwise). Crosslinking a genuinely novel species additionally requires adding matching rows to `crosslinks.csv` yourself, since crosslink positions are looked up by exact species name. See [`docs/examples/example8`](examples/example8/) for the full workflow, including generating a FASTA from an existing PDB with the bundled `pdb2fasta` script.

#### Crosslink Configuration

| Parameter | Type | Description | Valid Values | Example |
|-----------|------|-------------|--------------|---------|
| n_term_type | string | N-terminal crosslink type | See below* | "HLKNL" |
| c_term_type | string | C-terminal crosslink type | See below* | "HLKNL" |
| n_term_combination | string | N-terminal residue combination | Format varies by crosslink type** | "9.C - 947.A" |
| c_term_combination | string | C-terminal residue combination | Format varies by crosslink type** | "1047.C - 104.C" |

*Available crosslink types:

**Enzymatic Divalent (2 residues):**
- "HLKNL": Hydroxy-lysinonorleucine (mature divalent, most common)
- "LKNL": Lysinonorleucine (immature divalent)
- "deHLNL": Dehydro-hydroxylysinonorleucine (immature precursor)
- "deHHLNL": Dehydro-dihydroxylysinonorleucine (immature precursor)

**Enzymatic Trivalent (3 residues):**
- "PYD": Pyridinoline (hydroxylysine-derived, most common trivalent)
- "DPD": Deoxypyridinoline (lysine-derived trivalent)
- "PYL": Pyrrole (alternative trivalent pathway)
- "DPL": Deoxypyrrole (alternative trivalent pathway)

**Non-Enzymatic Divalent, LYS-LYS derived (2 residues):**
- "MOLD": Methylglyoxal-lysine dimer

**Non-Enzymatic (AGE), LYS-ARG derived (2 residues):**
- "Pentosidine": Well-characterized AGE crosslink
- "Glucosepane": Most abundant AGE in human tissue

Non-enzymatic types can be used as terminal crosslinks (`n_term_type`/`c_term_type`)
just like the enzymatic ones, or added on top of an existing structure via the
mutated PDB workflow (`additional_1_type`/`additional_2_type`) — most commonly
the latter, since AGE crosslinks are typically combined with an existing
enzymatic crosslink rather than used alone.

**Format for residue combinations:
- **Divalent**: "ResNum.Chain - ResNum.Chain" (e.g., "9.C - 947.A")
- **Trivalent**: "ResNum.Chain - ResNum.Chain - ResNum.Chain" (e.g., "6.B - 9.C - 946.A")

**Validation Rules**:
- For human (homo_sapiens) HLKNL crosslinks:
  - Valid N-terminal combinations: "5.B - 944.B", "9.C - 944.B", "9.C - 947.A", "947.A - 5.B"
  - Valid C-terminal combinations: "104.C - 1047.A", "1047.A - 98.B", "1047.C - 104.C", "1047.C - 98.B"
- For human (homo_sapiens) PYD crosslinks:
  - Valid N-terminal combination: "6.B - 9.C - 946.A"
  - Valid C-terminal combination: "1046.C - 1046.A - 103.C"
- Combinations vary by species and crosslink type
- Divalent format validation: must match pattern "^\d+\.[A-C]\s*-\s*\d+\.[A-C]$"
- Trivalent format validation: must match pattern "^\d+\.[A-C]\s*-\s*\d+\.[A-C]\s*-\s*\d+\.[A-C]$"

**Crosslink type validation against an input PDB**:
- When a `pdb_file` is provided together with `n_term_type`/`c_term_type` (and `crosslink` or `replace_bool` is enabled), ColBuilder checks that the crosslink residues actually present in the PDB match the declared types.
- This check considers `additional_1_type`/`additional_2_type` as well, not just the terminal types — so a mixed structure from the AGE workflow (e.g. terminal PYD plus an additional Glucosepane) validates correctly as long as every type actually present is declared somewhere across those four fields.
- A mismatch (e.g. a trivalent structure with only a divalent type declared) raises `GEO_ERR_008`, listing the detected crosslink residues.

#### Additional Crosslinks (Mutated PDB Workflow)

| Parameter | Type | Description | Valid Values | Example |
|-----------|------|-------------|--------------|---------|
| additional_1_type | string/null | First additional crosslink type | AGE crosslink types or null | "Glucosepane" |
| additional_1_combination | string/null | Position for first additional crosslink | "ResNum.Chain - ResNum.Chain" | "1008.A - 767.B" |
| additional_2_type | string/null | Second additional crosslink type (optional) | AGE crosslink types or null | "Pentosidine" |
| additional_2_combination | string/null | Position for second additional crosslink | "ResNum.Chain - ResNum.Chain" | "950.B - 710.C" |
| crosslink_copies | list of strings | Periodic copies for distance optimization | Exactly 2 different elements from D0-D5 | ["D0", "D5"] |

**Notes**:
- Additional crosslinks are used in the **mutated PDB workflow** to add AGE or other crosslinks on top of existing enzymatic crosslinks
- Requires `mutated_pdb` to be set to an existing PDB file
- **Workflow**: Run sequence generation with `mutated_pdb` and `additional_*` parameters separately, then use output in geometry generation, declaring the additional type(s) there too so crosslink-type validation matches the structure
- `crosslink_copies` must be exactly 2 different elements from D0, D1, D2, D3, D4, D5 (duplicates and any other value are rejected — there is no range syntax)
- Default `crosslink_copies` is ["D0", "D5"] if not specified

*Check available crosslinks and respective combinations at [src/colbuilder/data/sequence/crosslinks.csv](https://github.com/graeter-group/colbuilder/blob/main/src/colbuilder/data/sequence/crosslinks.csv)

### Geometry Generation Parameters

| Parameter | Type | Description | Valid Values | Default |
|-----------|------|-------------|--------------|---------|
| geometry_generator | boolean | Enable geometry generation | true/false | false |
| pdb_file | Path/null | Input PDB file (if sequence_generator=false) | Valid PDB file path or null | null |
| contact_distance | float | Contact distance for microfibril (Å) | Positive number (typically 15-40) | None (required*) |
| fibril_length | float | Length of microfibril (nm) | Positive number | None (required*) |
| crystalcontacts_file | Path/null | File with crystal contacts | Valid file path or null | null |
| connect_file | Path/null | File with connection information | Valid file path or null | null |
| crystalcontacts_optimize | boolean | Optimize crystal contacts | true/false | false |
| solution_space | List/Tuple | Solution space dimensions [dx, dy, dz] | Three positive numbers | [1, 1, 1] |
| pdb_first_line | string | Crystal contacts information | Valid PDB CRYST1 line | Default crystal parameters |

\* `fibril_length` is required for geometry generation, with one exception: in
mixing mode (`mix_bool: true`) it defaults to 40.0 nm if not set.
`contact_distance` is required unless `crystalcontacts_file` is provided instead.

**Notes**:
- Either `contact_distance` or `crystalcontacts_file` must be provided when geometry_generator is true
- When using mutated PDB workflow with additional crosslinks, use the output PDB from sequence generation as `pdb_file`
- **IMPORTANT**: Run sequence generation separately first when adding additional crosslinks, then use its output here

### Mixing and Replacement Parameters

| Parameter | Type | Description | Valid Values | Default |
|-----------|------|-------------|--------------|---------|
| mix_bool | boolean | Generate mixed crosslinked microfibril | true/false | false |
| ratio_mix | Dict or string | Ratio for mix-crosslink setup | Format: "Type:percentage Type:percentage" or Dict[str, int] | None |
| files_mix | List of Paths | PDB files with different crosslink types | Valid PDB file paths (≥2 files) | None |
| replace_bool | boolean | Replace crosslinks with standard residues | true/false | false |
| auto_fix_unpaired | boolean | Automatically detect crosslink markers left without a partner and replace them | true/false | false |
| ratio_replace | float | Percentage of eligible crosslinks to REMOVE (replace with standard residues) | 0-100 | None |
| ratio_replace_scope | string | Which crosslinks are eligible for ratio-based replacement | "enzymatic", "non_enzymatic", "all" | "enzymatic" |
| replace_file | Path/null | Input PDB file of fibril with crosslinks | Valid file path or null | null |
| manual_replacements | list of strings/null | Explicit, deterministic replacement directives | One `"<caps_file> <RES> <resid> <chain>"` entry per residue | null |

**Validation Rules**:
- When `mix_bool=true`:
  - `ratio_mix` must be provided with at least 2 types
  - `files_mix` must contain at least 2 valid PDB file paths
  - Percentages in `ratio_mix` must sum to 100
  - Number of files in `files_mix` must match number of types in `ratio_mix`
- When `replace_bool=true`:
  - Choose exactly one replacement mechanism: `ratio_replace` (0-100, with `ratio_replace_scope`) or `manual_replacements`
  - Either `geometry_generator=true` OR `replace_file` must be provided
  - If `geometry_generator=false`, must provide `replace_file`
  - `ratio_replace_scope` must be one of `enzymatic`, `non_enzymatic`, or `all`

**Notes**:
- **Mixing** creates heterogeneous microfibrils with different crosslink types (e.g., 80% divalent + 20% trivalent)
- **Replacement** simulates partial crosslinking or aged collagen by replacing some crosslinks with standard residues. There are two mechanisms:
  - **Ratio-based** (`ratio_replace` + `ratio_replace_scope`): a random percentage of eligible crosslinks, drawn from the chosen scope. `ratio_replace` is the fraction **removed**, not the remaining density — e.g. `ratio_replace: 70` removes 70% of eligible crosslinks, leaving 30%.
  - **Manual** (`manual_replacements`): exact residues named explicitly, for reproducible, targeted edits. Each entry targets a residue inside a per-model `{id}.caps.pdb` file, only written to disk under `debug: true`.
- `ratio_replace_scope` selects which crosslinks may be replaced by ratio-based replacement. The default `enzymatic` targets enzymatic crosslinks (HLKNL/PYD-derived residues); `non_enzymatic` targets AGE crosslinks (Glucosepane, Pentosidine) and MOLD; `all` considers both
- `auto_fix_unpaired` automatically finds crosslink markers that would otherwise be left without a geometric partner and converts them to standard residues. If you've already configured your own replacement (`ratio_replace`, `replace_file`, or `manual_replacements`), ColBuilder defers to it instead of silently overriding it with its own auto-fix list — it logs a warning rather than changing your config.
- Set `replace_file: null` to use geometry generation output for replacement

### Topology Generation Parameters

| Parameter | Type | Description | Valid Values | Default |
|-----------|------|-------------|--------------|---------|
| topology_generator | boolean | Generate topology files | true/false | false |
| force_field | string | Force field for simulations | "amber99", "martini3" | None |
| topology_debug | boolean | Save intermediate topology files | true/false | false |
| martinize2_command | string | Detected Martinize2 command path | Valid executable path | Auto-detected |
| go_epsilon | float | GO epsilon for Martini3-CG parametrization | Positive number | 9.414 |

**Notes**:
- **Topology-only mode**: When both `sequence_generator=false` and `geometry_generator=false` but `topology_generator=true`, ColBuilder generates topology from an existing fibril PDB
- Martini3 crosslink parametrization currently only covers PYD and HLKNL; other crosslink types (DPD, PYL, DPL, MOLD, and the non-enzymatic AGE types) are not yet parametrized for Martini3 and will not produce correct coarse-grained bonded terms — use `force_field: "amber99"` for those
- Set `topology_debug: true` to preserve intermediate files for troubleshooting

## Internal Data Structures

### Sequence Generation Structures

| Name | Type | Description |
|------|------|-------------|
| CrosslinkPair | class | Represents a pair of residues involved in a crosslink |
| CrosslinkPosition | class | Represents a residue position with residue type, atom name, and position string |
| SequenceGenerator | class | Manages the generation of collagen structure from sequences |
| OptimizationState | class | Tracks state during crosslink optimization |

**Key Methods in SequenceGenerator**:
- `generate()`: Main method to generate collagen structure
- `_run_alignment()`: Performs sequence alignment
- `_run_modelling()`: Runs MODELLER for structure generation
- `_load_crosslinks()`: Loads crosslink information
- `_apply_crosslinks()`: Applies crosslinks to structure
- `_apply_additional_crosslinks()`: Applies additional crosslinks in mutated PDB workflow
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
- Standard temp files: "replace.txt", "manual_replacements.txt"
- Standard temp directories: "NC", "D", "T", "M", "DT", "TD", "replace_manual", "replace_crosslinks"

### Topology Generation Structures

| Name | Type | Description |
|------|------|-------------|
| Amber | class | Handles amber99 force field topology generation |
| Martini | class | Handles martini3 force field topology generation |
| Itp | class | Manages ITP file generation and crosslink bonded terms |

**Key Functions**:
- `build_topology()`: Main entry point for topology generation
- `build_amber99()`: Builds amber99 topology
- `build_martini3()`: Builds martini3 topology
- `make_topology()`: Generates topology for individual models

**Output Files**:
- `{species}_topology_files/`: Directory containing topology files
- `collagen_fibril_{species}.top`: Topology file
- `collagen_fibril_{species}.gro`: GROMACS coordinate file
- Various .itp files: Include topology files (e.g., `col_0.itp`, `col_1.itp`, etc.)

**Temporary Files**:
- Force field files (amber99sb-star-ildnp.ff/ or similar)
- Various intermediate files (*.itp, *.CG.pdb, *.merge.pdb, etc.)

## Glossary

- **Collagen**: Main structural protein in connective tissues
- **Triple helix**: Three polypeptide chains wound around each other, forming the basic collagen structure
- **Microfibril**: Assembly of multiple collagen molecules in a quasi-hexagonal packing
- **Crosslink**: Covalent bond connecting collagen chains, providing structural stability
- **Enzymatic crosslink**: Crosslink formed through lysyl oxidase (LOX) enzyme activity
- **Non-enzymatic crosslink (AGE)**: Advanced Glycation End-product formed spontaneously, typically accumulating with aging
- **Divalent crosslink**: Crosslink connecting 2 residues (e.g., HLKNL, LKNL, MOLD)
- **Trivalent crosslink**: Crosslink connecting 3 residues (e.g., PYD, DPD)
- **Mutated PDB workflow**: Method to add additional crosslinks to existing structures
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
- **Topology-only mode**: Operating mode that generates topology files from existing fibril PDB without running geometry generation
- **ITP file**: GROMACS Include Topology file containing molecular topology information