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
    - [Adding AGE Crosslinks (Mutated PDB Workflow)](#adding-age-crosslinks-mutated-pdb-workflow)
    - [Topology-Only Generation](#topology-only-generation)
  - [Parameter Interactions and Dependencies](#parameter-interactions-and-dependencies)
  - [Configuration Validation](#configuration-validation)

## Introduction

ColBuilder uses a YAML configuration file to control all aspects of the pipeline. This document provides a comprehensive reference for all available configuration parameters, their meanings, and how they interact with each other.

Only set the options a given run actually needs — every parameter below has a
default, and most examples only declare a handful of them. See
[`docs/examples/`](examples/) for complete, runnable configs covering each
workflow described here.

## Configuration File Format

Configuration files for ColBuilder are written in YAML format. A basic configuration file might look like this:

```yaml
# Basic configuration
working_directory: "./"
debug: false

# Operation modes
sequence_generator: true
geometry_generator: true
topology_generator: false

# Input Configuration
species: "homo_sapiens"
mutated_pdb: null
fasta_file: null

# Sequence Settings
crosslink: true
n_term_type: "HLKNL"
c_term_type: "HLKNL"
n_term_combination: "9.C - 947.A"
c_term_combination: "1047.C - 104.C"

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
ratio_replace_scope: "enzymatic"   # enzymatic | non_enzymatic | all
replace_file: null
manual_replacements: null

# Topology Options
force_field: "amber99"
topology_debug: false
```

## Operation Mode Parameters

These parameters control which stages of the pipeline are executed.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `config_file` | string, null | null | Path to another config file (for nested configs) |
| `sequence_generator` | boolean | false | Enable/disable sequence generation stage |
| `geometry_generator` | boolean | false | Enable/disable geometry generation stage |
| `topology_generator` | boolean | false | Enable/disable topology generation stage |
| `debug` | boolean | false | Keep intermediate files for debugging |
| `working_directory` | string | current directory | Directory for input and output files |

**Notes**:
- All three stage flags (`sequence_generator`, `geometry_generator`, `topology_generator`) default to `false`; a config must opt into whichever stages it needs.
- If all three are set to `true`, the full pipeline runs in sequence in a single invocation.
- The `working_directory` parameter sets the base directory for all input and output files.
- **Topology-only mode**: When `topology_generator: true` and both `sequence_generator: false` and `geometry_generator: false`, ColBuilder operates in topology-only mode, generating topology files from an existing fibril PDB without running geometry generation.

## Input Configuration Parameters

These parameters control basic input settings, particularly species selection.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `species` | string | (required) | Species for collagen sequence and structure |
| `mutated_pdb` | string, null | null | Path to pre-mutated PDB for adding additional crosslinks |
| `fasta_file` | string, null | null | Path to custom FASTA file |

**Available Species Options**:
- **Mammals (Primates)**: homo_sapiens, pan_troglodytes, pongo_abelii, callithrix_jacchus, otolemur_garnettii
- **Mammals (Rodents)**: mus_musculus, rattus_norvegicus
- **Mammals (Other)**: bos_taurus, canis_lupus, ailuropoda_melanoleuca, mustela_putorius, myotis_lucifugus, loxodonta_africana
- **Fish**: danio_rerio, oreochromis_niloticus, oryzias_latipes, tetraodon_nigroviridis, xiphophorus_maculatus
- **Reptiles**: pelodiscus_sinensis

**Notes**:
- The `species` parameter is used to select the appropriate sequence and crosslink information; it is always required.
- Rat collagen (rattus_norvegicus) is used as the default structural template for all species.
- **Mutated PDB workflow**: Set `mutated_pdb` to the path of an existing PDB (e.g., output from a previous sequence generation run) to add additional crosslinks on top of existing ones. This is particularly useful for adding AGE crosslinks.
- **Custom species**: If `species` is not one of the built-in species listed above, `fasta_file` is required — ColBuilder raises a configuration error otherwise. This lets you build a fibril for a species with no bundled sequence, provided you supply its FASTA yourself (see [`docs/examples/example8`](examples/example8/) for the full "bring your own sequence" workflow, including generating that FASTA from an existing PDB with the bundled `pdb2fasta` script). Crosslinking a genuinely novel species additionally requires adding matching rows to `crosslinks.csv` yourself, since crosslink residue positions are looked up by exact species name.
- If `fasta_file` is provided for a built-in species, it overrides the bundled sequence.

## Sequence Generation Parameters

These parameters control the sequence generation stage (homology modeling).

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `crosslink` | boolean | false | Enable/disable crosslinking in the model |
| `n_term_type` | string, null | null | N-terminal crosslink type |
| `c_term_type` | string, null | null | C-terminal crosslink type |
| `n_term_combination` | string, null | null | N-terminal residue combination for crosslinks |
| `c_term_combination` | string, null | null | C-terminal residue combination for crosslinks |
| `additional_1_type` | string, null | null | First additional crosslink type (for mutated PDB workflow) |
| `additional_1_combination` | string, null | null | Position for first additional crosslink (required if `additional_1_type` is set) |
| `additional_2_type` | string, null | null | Second additional crosslink type (optional) |
| `additional_2_combination` | string, null | null | Position for second additional crosslink (required if `additional_2_type` is set) |
| `crosslink_copies` | list of strings | ["D0", "D5"] | Pair of unit cell translations for crosslink optimization (optional) |

**Available Crosslink Types**:

*Enzymatic Divalent (2 residues):*
- **HLKNL**: Hydroxylysino-5-ketonorleucine (mature divalent, most common)
- **LKNL**: Lysino-5-ketonorleucine (immature divalent)
- **deHLNL**: Dehydro-hydroxylysinonorleucine (immature precursor)
- **deHHLNL**: Dehydro-dihydroxylysinonorleucine (immature precursor)

*Non-Enzymatic Divalent, LYS-LYS derived (2 residues):*
- **MOLD**: Methylglyoxal-lysine dimer

*Enzymatic Trivalent (3 residues):*
- **PYD**: Pyridinoline (hydroxylysine-derived, most common trivalent)
- **DPD**: Deoxypyridinoline (lysine-derived trivalent)
- **PYL**: Pyrrole (alternative trivalent pathway)
- **DPL**: Deoxypyrrole (alternative trivalent pathway)

*Non-Enzymatic (AGE), LYS-ARG derived:*
- **Pentosidine**: Well-characterized AGE crosslink
- **Glucosepane**: Most abundant AGE in human tissue

**Residue Combination Format**:
- **Divalent crosslinks**: `"resid.chain - resid.chain"` (e.g., `"9.C - 947.A"`)
- **Trivalent crosslinks**: `"resid.chain - resid.chain - resid.chain"` (e.g., `"6.B - 9.C - 946.A"`)

**Residue Combination Examples for Homo Sapiens**:

*HLKNL (divalent):*
- **N-terminal**: "5.B - 944.B", "9.C - 944.B", "9.C - 947.A", "947.A - 5.B"
- **C-terminal**: "104.C - 1047.A", "1047.A - 98.B", "1047.C - 104.C", "1047.C - 98.B"

*PYD (trivalent):*
- **N-terminal**: "6.B - 9.C - 946.A"
- **C-terminal**: "1046.C - 1046.A - 103.C"

**Crosslink Optimization**:
- `crosslink_copies` specifies which periodic copies to use for distance optimization.
- Must be exactly 2 **different** elements from: D0, D1, D2, D3, D4, D5 — duplicates and any other value are rejected.
- Default is `["D0", "D5"]` if not specified.
- Custom values like `["D2", "D3"]` can improve optimization for specific crosslink positions (e.g. additional AGE crosslinks placed elsewhere on the structure than the terminal ones).

**Notes**:
- The `crosslink` parameter must be set to `true` for crosslinking to be applied.
- Terminal crosslinks (`n_term_*` and `c_term_*`) are typically applied in standard sequence generation.
- Additional crosslinks (`additional_*`) are used in the mutated PDB workflow to add AGE or other crosslinks on top of existing structures. If `additional_1_type`/`additional_2_type` is set, the matching `additional_1_combination`/`additional_2_combination` is required.
- A complete list of available species and combinations can be found at [src/colbuilder/data/sequence/crosslinks.csv](https://github.com/graeter-group/colbuilder/blob/main/src/colbuilder/data/sequence/crosslinks.csv).
- When using the mutated PDB workflow, run sequence generation separately first, then use the output in geometry generation.

**Crosslink type validation**:
- When you provide an input `pdb_file` together with `n_term_type`/`c_term_type` (and `crosslink` or `replace_bool` is enabled), ColBuilder checks that the crosslink residues present in the PDB are consistent with the crosslink types you've specified. Divalent types (HLKNL, LKNL, deHLNL, deHHLNL, MOLD) expect a divalent structure, and trivalent types (PYD, DPD, PYL, DPL) expect a trivalent structure.
- This check also considers `additional_1_type`/`additional_2_type`, not just the terminal types — so an AGE workflow's mixed structure (e.g. terminal PYD plus an additional Glucosepane) validates correctly as long as **all** the types actually present are declared somewhere across `n_term_type`/`c_term_type`/`additional_1_type`/`additional_2_type`.
- If the PDB contains crosslinks of a category that was not requested anywhere in those fields (for example, a trivalent PDB while only `n_term_type: "HLKNL"` is set), generation stops with error `GEO_ERR_008` and a message listing the detected crosslink residues. Update the type fields to match the structure, or supply a PDB whose crosslinks match the specified types.

## Geometry Generation Parameters

These parameters control the generation of the microfibril structure.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `pdb_file` | string, null | null | Input PDB file (set to null if sequence_generator is true) |
| `contact_distance` | float, null | (required*) | Distance threshold for contacts (Å) |
| `fibril_length` | float, null | (required*) | Length of the generated fibril (nm) |
| `crystalcontacts_file` | string, null | null | File with crystal contacts information |
| `connect_file` | string, null | null | File with connection information |
| `crystalcontacts_optimize` | boolean | false | Optimize crystal contacts during generation |

\* `fibril_length` is required for geometry generation, with one exception: in
mixing mode (`mix_bool: true`) it defaults to 40.0 nm if not set.
`contact_distance` is required unless `crystalcontacts_file` is provided instead.

**Notes**:
- The `contact_distance` parameter controls the radial size of the microfibril. Larger values result in thicker fibrils but require more computation time and memory.
- The `fibril_length` parameter sets the length of the fibril in nanometers. Larger values create longer fibrils.
- If `crystalcontacts_optimize` is set to true, ColBuilder will attempt to optimize the packing of collagen molecules in a Bravais lattice.
- If `pdb_file` is provided, the sequence generation stage will be skipped and this PDB will be used as input for geometry generation.
- **IMPORTANT**: When using additional crosslinks (mutated PDB workflow), run sequence generation separately first, then use its output as `pdb_file` in geometry generation.

## Mixing and Replacement Parameters

These parameters control advanced features for creating mixed crosslinked microfibrils or replacing crosslinks with standard residues.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `mix_bool` | boolean | false | Enable mixing of different crosslink types |
| `ratio_mix` | string, null | null | Format: "Type:percentage Type:percentage" (required if `mix_bool` is true) |
| `files_mix` | list of strings, null | null | Required if `mix_bool` is true, paths to PDB files with different crosslink types |
| `replace_bool` | boolean | false | Enable crosslink replacement (with standard residues) |
| `auto_fix_unpaired` | boolean | false | Automatically detect crosslink markers left without a partner and replace them |
| `ratio_replace` | float, null | null | Percentage of eligible crosslinks to REMOVE, i.e. replace with standard residues (0-100) |
| `ratio_replace_scope` | string | "enzymatic" | Which crosslinks are eligible for ratio-based replacement: `enzymatic`, `non_enzymatic`, or `all` |
| `replace_file` | string, null | null | File with crosslinks to be replaced (defaults to the geometry generation output if unset) |
| `manual_replacements` | list of strings, null | null | Explicit, deterministic replacement directives: one `"<caps_file> <RES> <resid> <chain>"` entry per residue |

**Notes**:
- The `mix_bool` feature allows creation of heterogeneous crosslinked microfibrils, which more closely resemble natural collagen.
- The `ratio_mix` parameter specifies the proportion of each crosslink type in the mixed microfibril. Percentages must sum to 100.
- The `files_mix` parameter specifies paths to PDB files of collagen molecules, each with a different crosslink type.
- The `replace_bool` feature simulates partial crosslinking or aged collagen by replacing some crosslinks with unmodified standard residues. There are two ways to choose which crosslinks get replaced — pick one:
  - **Ratio-based** (`ratio_replace` + `ratio_replace_scope`): a random percentage of eligible crosslinks, drawn from the chosen scope.
  - **Manual** (`manual_replacements`): exact residues named explicitly, for reproducible, targeted edits (e.g. removing one specific crosslink to test its mechanical contribution). Each entry targets a residue inside a per-model `{id}.caps.pdb` file, which is only written to disk under `debug: true`.
- `ratio_replace` is the fraction **removed**, not the remaining density: e.g. `ratio_replace: 70` removes 70% of eligible crosslinks, leaving 30%.
- The `ratio_replace_scope` parameter selects which crosslinks may be replaced by ratio-based replacement. **The default is `enzymatic`**, so ratio-based replacement targets enzymatic crosslinks (e.g. HLKNL/PYD-derived residues) unless you choose otherwise. Use `non_enzymatic` to target advanced glycation end-product (AGE) crosslinks (Glucosepane, Pentosidine) and MOLD, or `all` to consider both.
- `auto_fix_unpaired` automatically finds crosslink markers that would otherwise be left without a geometric partner (which would produce an incomplete, non-physical crosslink) and converts them to standard residues. If you've already requested your own replacement (`ratio_replace`, `replace_file`, or `manual_replacements`), ColBuilder defers to your request instead of silently overriding it with its own auto-fix list — it logs a warning rather than changing your config.
- The `replace_file` parameter specifies the path to a PDB file of a previously generated collagen microfibril. Set to null to use the geometry generation output.

## Topology Generation Parameters

These parameters control the generation of topology files for molecular dynamics simulations.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `force_field` | string, null | null | Force field for topology generation |
| `topology_debug` | boolean | false | Save intermediate topology files for debugging |
| `go_epsilon` | float | 9.414 | Go-model epsilon value for Martini3 coarse-grained parametrization |

**Available Force Field Options**:
- **amber99**: Standard Amber99 force field for atomistic simulations
- **martini3**: Martini 3 coarse-grained force field

**Notes**:
- The `force_field` parameter selects which force field to use for generating topology files; it is required whenever `topology_generator: true`.
- The amber99 force field is recommended for most atomistic simulations of collagen.
- Custom force field parameters for collagen and crosslinks are included in ColBuilder.
- Martini3 crosslink parametrization currently only covers PYD and HLKNL; other crosslink types (DPD, PYL, DPL, MOLD, and the non-enzymatic AGE types) are not yet parametrized for Martini3 and will not produce correct coarse-grained bonded terms — use `force_field: "amber99"` for those.
- Set `topology_debug: true` to keep intermediate files for troubleshooting topology generation issues.
- **Topology-only mode**: When both `sequence_generator` and `geometry_generator` are false, topology generation can process existing fibril PDB files directly.

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
n_term_combination: "9.C - 5.B - 942.B"
c_term_combination: "1046.C - 1046.A - 103.C" 
fibril_length: 80.0
contact_distance: 20
```

### Mixed Crosslinked Microfibril (80% Divalent + 20% Trivalent crosslink)

```yaml
species: "homo_sapiens"
sequence_generator: false
geometry_generator: true
mix_bool: true
ratio_mix: "D:80 T:20"
files_mix:
 - "homo_sapiens_N_HLKNL_C_HLKNL.pdb"
 - "homo_sapiens_N_PYD_C_PYD.pdb"
contact_distance: 20
fibril_length: 40.0
```

### Microfibril and Coarse-Grained Topology Generation

```yaml
species: "homo_sapiens"
sequence_generator: false
geometry_generator: true
topology_generator: true
pdb_file: "path/to/collagen_molecule.pdb"
fibril_length: 70.0
contact_distance: 30
force_field: "martini3"
```

### Adding AGE Crosslinks (Mutated PDB Workflow)

**Step 1**: Generate base structure with terminal crosslinks
```yaml
# config_step1.yaml
species: "rattus_norvegicus"
sequence_generator: true
geometry_generator: false
crosslink: true
n_term_type: "PYD"
c_term_type: "PYD"
n_term_combination: "6.B - 9.C - 946.A"
c_term_combination: "1046.C - 1046.A - 103.C"
```

**Step 2**: Add AGE crosslink to the output from step 1
```yaml
# config_step2.yaml
species: "rattus_norvegicus"
sequence_generator: true
geometry_generator: false
mutated_pdb: "rattusnorvegicus_N_PYD_C_PYD.pdb"
crosslink: true
n_term_type: "PYD"
c_term_type: "PYD"
n_term_combination: "6.B - 9.C - 946.A"
c_term_combination: "1046.C - 1046.A - 103.C"
additional_1_type: "Glucosepane"
additional_1_combination: "1055.C - 822.A"
crosslink_copies: ["D2", "D3"]
```

**Step 3**: Build fibril from the mutated structure
```yaml
# config_step3.yaml
species: "rattus_norvegicus"
sequence_generator: false
geometry_generator: true
topology_generator: true
pdb_file: "rattusnorvegicus_N_PYD_C_PYD+ADD1_Glucosepane.pdb"
crosslink: true
n_term_type: "PYD"
c_term_type: "PYD"
n_term_combination: "6.B - 9.C - 946.A"
c_term_combination: "1046.C - 1046.A - 103.C"
# Must also be declared: pdb_file above already contains Glucosepane (added
# in step 2), and crosslink-type validation checks against everything
# actually present in the structure, not just the terminal types.
additional_1_type: "Glucosepane"
fibril_length: 40.0
contact_distance: 20
force_field: "amber99"
```

### Topology-Only Generation

Generate topology files from an existing fibril PDB:

```yaml
species: "homo_sapiens"
sequence_generator: false
geometry_generator: false
topology_generator: true
pdb_file: "path/to/existing_fibril.pdb"
crosslink: true
n_term_type: "PYD"
c_term_type: "PYD"
n_term_combination: "6.B - 9.C - 946.A"
c_term_combination: "1046.C - 1046.A - 103.C"
force_field: "amber99"
```

**Note**: Crosslink parameters must match what's present in the input PDB file (including any `additional_*_type` present — see the crosslink type validation note above). If the input PDB has more than one crosslink type mixed together (e.g. output from a mixing or replacement step) such that no single set of type fields describes it, set `crosslink: false` instead — validation is skipped, since a single declared type wouldn't be meaningful for a genuinely mixed structure.

## Parameter Interactions and Dependencies

Understanding how parameters interact is important for successful use of ColBuilder:

1. **Pipeline Stage Dependencies**:
   - If `sequence_generator` is true, a new PDB with the structure of a collagen molecule will be generated, and any provided `pdb_file` will be ignored.
   - If `sequence_generator` is false but `geometry_generator` is true, you must provide a valid `pdb_file`.
   - **Topology-only mode**: When both `sequence_generator` and `geometry_generator` are false but `topology_generator` is true, ColBuilder generates topology from an existing fibril PDB.

2. **Crosslinking Dependencies**:
   - If `crosslink` is true, you must specify valid `n_term_type`, `c_term_type`, `n_term_combination`, and `c_term_combination` values.
   - The crosslink type and residue combinations must be compatible with the chosen species.
   - **Trivalent crosslinks (e.g., PYD, DPD, PYL, DPL)** require 3 residue positions per terminus.
   - **Divalent crosslinks (e.g., HLKNL, LKNL, MOLD)** require 2 residue positions per terminus.

3. **Mutated PDB Workflow Dependencies**:
   - When `mutated_pdb` is provided, sequence generation adds crosslinks to the existing structure.
   - Terminal crosslinks (`n_term_*`, `c_term_*`) must match those already in the mutated PDB.
   - Additional crosslinks (`additional_1_*`, `additional_2_*`) are added on top of existing ones; each `additional_N_type` requires the matching `additional_N_combination`.
   - **IMPORTANT**: Run sequence generation with additional crosslinks separately, then use the output in geometry generation, declaring the additional type(s) there too so crosslink-type validation matches the structure.

4. **Mixing and Replacement Dependencies**:
   - If `mix_bool` is true, you must provide valid files in `files_mix` and a proper ratio in `ratio_mix`.
   - Ratios in `ratio_mix` must sum to 100.
   - If `replace_bool` is true, choose exactly one replacement mechanism: `ratio_replace` (0-100, with `ratio_replace_scope`) or `manual_replacements`.
   - If `replace_bool` is true and `geometry_generator` is false, you must provide a `replace_file`.
   - `ratio_replace_scope` must be one of `enzymatic` (default), `non_enzymatic`, or `all`.
   - `auto_fix_unpaired` will not override a replacement mechanism you've already configured (`ratio_replace`, `replace_file`, or `manual_replacements`) — it only acts when none of those are set.

5. **Geometry Generation Dependencies**:
   - If `crystalcontacts_optimize` is true, the geometry generation will take longer but may produce better-packed microfibrils.
   - The `contact_distance` parameter becomes irrelevant if `crystalcontacts_file` is provided; one of the two is required.
   - `fibril_length` is required for geometry generation, except in mixing mode where it defaults to 40.0 nm.

6. **Crosslink Optimization Dependencies**:
   - `crosslink_copies` must contain exactly 2 elements.
   - Valid elements: D0, D1, D2, D3, D4, D5 — the two elements must be different from each other.
   - Default is `["D0", "D5"]` if not specified.

## Configuration Validation

ColBuilder performs validation of the configuration file before running the pipeline:

1. **Required Parameters**: Checks that all required parameters are present
2. **Parameter Types**: Validates that parameters have the correct data types
3. **Valid Values**: Ensures parameters have valid values (e.g., valid species name, crosslink types)
4. **Consistency**: Verifies that parameter combinations are consistent and compatible
5. **File Existence**: Checks that any specified input files exist and are readable
6. **Crosslink Validation**:
   - Verifies that crosslink types are valid
   - Checks that residue combinations match the required format (2 for divalent, 3 for trivalent)
   - Ensures combinations are compatible with the chosen species
   - Cross-checks the crosslink residues actually present in an input `pdb_file` against the declared type fields (see [Crosslink type validation](#sequence-generation-parameters) above)
7. **Ratio Validation**:
   - For mixing: Verifies ratios sum to 100
   - For replacement: Checks ratio is between 0-100

If validation fails, ColBuilder will provide an error message indicating the specific issue with the configuration file.

For any issues with configuration, refer to the error message or enable debug mode (`debug: true`) for more detailed information about the validation process.
