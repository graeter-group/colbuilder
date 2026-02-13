<div align="center">
    <h1>ColBuilder</h1>
    <p>Generate atomistic and coarse-grained models of collagen microfibrils with customizable crosslinking</p>
    <img src="https://github.com/user-attachments/assets/e09bda5a-04e4-46ad-b03f-3bcb3346a52f" alt="colbuilder-schematic-orange-compressed" width="80%">
</div>

## 📋 Table of Contents
- [📋 Table of Contents](#-table-of-contents)
- [📚 About](#-about)
  - [Key Features](#key-features)
- [🚀 Installation](#-installation)
  - [Prerequisites](#prerequisites)
  - [Step-by-Step Installation](#step-by-step-installation)
  - [Dependencies](#dependencies)
    - [PyMOL](#pymol)
    - [muscle (Multiple Sequence Alignment)](#muscle-multiple-sequence-alignment)
    - [UCSF Chimera](#ucsf-chimera)
    - [Modeller](#modeller)
- [🚀 Quick Start](#-quick-start)
- [⚙️ Operation Modes \& Workflow](#️-operation-modes--workflow)
  - [🧠 Understanding PDB Types](#-understanding-pdb-types)
  - [📊 Mode Summary Table](#-mode-summary-table)
  - [🔁 Valid Mode Combinations](#-valid-mode-combinations)
  - [✅ Valid Workflows](#-valid-workflows)
  - [🔧 Mode 4 (Mixing Crosslinks): Run Separately via Script](#-mode-4-mixing-crosslinks-run-separately-via-script)
- [🔗 Collagen Crosslinks](#-collagen-crosslinks)
  - [Enzymatic Crosslinks](#enzymatic-crosslinks)
  - [Non-Enzymatic Crosslinks (AGEs)](#non-enzymatic-crosslinks-ages)
  - [Crosslink Combinations](#crosslink-combinations)
- [📖 Usage Guide](#-usage-guide)
  - [Basic Usage](#basic-usage)
  - [Configuration Options](#configuration-options)
  - [Example Workflows](#example-workflows)
    - [Creating a Basic Human Collagen Microfibril](#creating-a-basic-human-collagen-microfibril)
    - [Generating a Crosslinked Bovine Microfibril](#generating-a-crosslinked-bovine-microfibril)
    - [Creating a Mixed Crosslinked (80% Divalent + 20% Trivalent) Human Collagen Microfibril from Collagen Molecules](#creating-a-mixed-crosslinked-80-divalent--20-trivalent-human-collagen-microfibril-from-collagen-molecules)
    - [Generating a Coarse-Grained Topology File for MD Simulation](#generating-a-coarse-grained-topology-file-for-md-simulation)
    - [Generating Topology Files from an Existing Fibril (Topology-Only Mode)](#generating-topology-files-from-an-existing-fibril-topology-only-mode)
    - [Creating a Fibril with PYD Crosslinks](#creating-a-fibril-with-pyd-crosslinks)
    - [Replacing Crosslinks in an Existing Fibril](#replacing-crosslinks-in-an-existing-fibril)
    - [Adding AGE Crosslinks to an Existing Structure (Mutated PDB Workflow)](#adding-age-crosslinks-to-an-existing-structure-mutated-pdb-workflow)
- [📚 Documentation](#-documentation)
- [🤝 Contributing](#-contributing)
- [📚 Publications \& Citation](#-publications--citation)
- [🙏 Acknowledgements](#-acknowledgements)

## 📚 About

**ColBuilder** is a specialized tool for generating atomistic and coarse-grained models of collagen microfibrils from single collagen molecules. Developed by the Gräter group at the Max Planck Institute for Polymer Research, it provides researchers with a flexible framework to create biologically relevant collagen structures for molecular dynamics simulations and structural studies.

### Key Features

- **Custom microfibril generation**: Create collagen microfibrils from individual molecules or amino acid sequences with precise control over structural parameters
- **Comprehensive crosslinking support**: Model enzymatic (divalent and trivalent) and non-enzymatic (e.g., AGE) crosslinks with flexible positioning
- **Highly configurable**: Adjust collagen sequence, fibril geometry, crosslink types and density to match your custom conditions
- **Multiple force fields**: Generate topology files for both atomistic (Amber99) and coarse-grained (Martini3) simulations
- **Flexible workflows**: Run complete pipeline or individual steps
- **Crosslink manipulation**: Mix different crosslink types or reduce crosslink density in existing fibrils
- **Simulation-ready output**: Generate complete topology files compatible with GROMACS and other major MD packages
- **Reproducible research**: Standardized approach to collagen modeling ensures consistency across studies

## 🚀 Installation

### Prerequisites

- Python 3.9 or later
- Git
- Conda package manager (we recommend [miniforge](https://github.com/conda-forge/miniforge))

### Step-by-Step Installation

1. **Create and activate a conda environment**:
   ```bash
   conda create -n colbuilder python=3.9
   conda activate colbuilder
   ```

2. **Clone the repository**:
   ```bash
   git clone git@github.com:graeter-group/colbuilder.git
   cd colbuilder
   ```

3. **Install ColBuilder**:
   ```bash
   pip install .
   ```

### Dependencies

ColBuilder requires several external tools to function properly:

#### PyMOL
```bash
conda install conda-forge::pymol-open-source
```

**Note**: If PyMOL fails due to missing `libnetcdf.so`, install:
```bash
conda install -c conda-forge libnetcdf==4.7.3
```

#### muscle (Multiple Sequence Alignment)
```bash
conda install bioconda::muscle
```

#### UCSF Chimera
1. Download the latest version of [UCSF Chimera](https://www.cgl.ucsf.edu/chimera/download.html) (64-bit recommended)
2. Make the binary executable and run the installer:
   ```bash
   cd ~/Downloads  # or wherever you downloaded the file
   chmod +x chimera*.bin
   ./chimera*.bin
   ```
3. Follow the installation prompts, preferably creating a symlink in a directory in your `$PATH`

**Note**: ColBuilder specifically requires UCSF Chimera, not the newer ChimeraX.

#### Modeller
1. Download [Modeller version 10.5](https://salilab.org/modeller/download_installation.html)
2. Follow the installation instructions provided
3. Add the following environment variables to your `.bashrc` or `.bash_profile`:
   ```bash
   export PYTHONPATH="/home/user/bin/modeller10.5/lib/x86_64-intel8/python3.3:$PYTHONPATH"
   export PYTHONPATH="/home/user/bin/modeller10.5/modlib:$PYTHONPATH"
   export LD_LIBRARY_PATH="/home/user/bin/modeller10.5/lib/x86_64-intel8:$LD_LIBRARY_PATH"
   ```
   (Adjust paths according to your installation location)

## 🚀 Quick Start

To verify your installation and run a basic example:

1. **Verify installation**:
   ```bash
   colbuilder --help
   ```

2. **Create a basic configuration file** (save as `config.yaml`):
   ```yaml
   # Basic human collagen microfibril configuration   
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

3. **Run ColBuilder**:
   ```bash
   colbuilder --config_file config.yaml
   ```

## ⚙️ Operation Modes & Workflow

ColBuilder operates through modular **modes**, each responsible for a different part of the collagen model-building pipeline. These modes can be combined in various ways or run separately using different configuration files.

### 🧠 Understanding PDB Types

ColBuilder produces or requires **two kinds of PDB files**:

- **Collagen triple helix molecule PDB**: A single ~300 nm-long collagen molecule (usually with specified crosslink residues). Output of **Mode 1**, input to **Modes 2** and **4**.
- **Collagen fibril PDB**: A full microfibril model composed of multiple triple helices arranged based on crystal geometry, length, and crosslinking. Output of **Modes 2, 4, or 5**, input to **Modes 3** and **5**.

Understanding this distinction is crucial for organizing your workflow correctly.

---

### 📊 Mode Summary Table

| # | Mode                   | Purpose                                                                 | Input(s)                                                       | Output                             | Can Run With Other Modes?   |
|---|------------------------|-------------------------------------------------------------------------|----------------------------------------------------------------|------------------------------------|------------------------------|
| 1 | `sequence_generator` | Generate a collagen triple helix molecule via homology modeling | `species` or custom FASTA | Triple helix PDB | Yes: with 2, 3, 5 |
| 2 | `geometry_generator` | Assemble a collagen fibril from a single triple helix | PDB from Mode 1 or custom PDB | Fibril PDB | Yes: with 1, 3, 5 |
| 3 | `topology_generator` | Generate topology files for GROMACS simulations | Fibril PDB (from Mode 2, 4, or 5) | `.top`, `.itp`, `.gro` | Yes: with 2, 4, 5 |
| 3* | `topology_generator` (topology-only) | Generate topology files from an existing fibril PDB without geometry generation | Pre-existing fibril PDB with crosslinks | `.top`, `.itp`, `.gro` | No, standalone mode |
| 4 | `mix_bool` | Generate a fibril by mixing multiple crosslink types | Two (or more) triple helix PDBs from Mode 1 | Mixed fibril PDB | No, requires separate script |
| 5 | `replace_bool` | Replace crosslinks in an existing fibril | Fibril PDB from Mode 2 or 4 | Modified fibril PDB | Yes: with 2, 3 |

---

### 🔁 Valid Mode Combinations

These combinations can be run **in a single config file**:

```yaml
# Example combination
sequence_generator: true
geometry_generator: true         # (optional)
topology_generator: true         # (optional)
replace_bool: true               # (optional)
```

### ✅ Valid Workflows

These mode combinations can be run **in a single configuration file**:

- ✅ **`1 + 2`** - Generate molecule and build fibril
- ✅ **`1 + 2 + 3`** - Complete pipeline: molecule → fibril → topology - [example1](docs/examples/example1)
- ✅ **`2 + 3`** - Build fibril and topology from custom triple helix PDB
- ✅ **`1 + 2 + 5`** - Generate, build, then replace crosslinks
- ✅ **`1 + 2 + 5 + 3`** - Generate, build, replace crosslinks by their standard amino acid, then create topology
- ✅ **`2 + 5`** - Build fibril then replace crosslinks - [example3](docs/examples/example3)
- ✅ **`2 + 5 + 3`** - Build fibril, replace crosslinks, then create topology
- ✅ **`3`** - Topology-only mode: generate topology from existing fibril PDB - [example4](docs/examples/example4)

---

### 🔧 Mode 4 (Mixing Crosslinks): Run Separately via Script

Mixing crosslinks (**Mode 4**) currently requires a separate workflow using (at least) two config files for triple helix generation and one for fibril construction (and optional topology generation):

[example2](docs/examples/example2)

```bash
# Example bash script for mixing crosslinks
colbuilder --config_file triple_helix_A.yaml
colbuilder --config_file triple_helix_B.yaml
colbuilder --config_file mix_geometry.yaml   # sets mix_bool: true and includes both PDBs
```

You can also chain this with `replace_bool` (Mode 5) or `topology_generator` (Mode 3) in the third config.

---

## 🔗 Collagen Crosslinks

ColBuilder supports a comprehensive range of collagen crosslinks, including both enzymatic and non-enzymatic types. Crosslinks are critical for collagen stability and mechanical properties.

### Enzymatic Crosslinks

Enzymatic crosslinks are formed through lysyl oxidase (LOX) activity and mature through various pathways:

**Divalent Crosslinks (2 residues):**
- **HLKNL** (Hydroxy-lysinonorleucine) - Mature divalent, most common
- **LKNL** (Lysinonorleucine) - Immature divalent
- **deHLNL** (Dehydro-hydroxylysinonorleucine) - Immature precursor
- **deHHLNL** (Dehydro-dihydroxylysinonorleucine) - Immature precursor

**Trivalent Crosslinks (3 residues):**
- **PYD** (Pyridinoline) - Hydroxylysine-derived, most common trivalent
- **DPD** (Deoxypyridinoline) - Lysine-derived trivalent
- **PYL** (Pyrrole) - Alternative trivalent pathway
- **DPL** (Deoxypyrrole) - Alternative trivalent pathway

### Non-Enzymatic Crosslinks (AGEs)

Advanced Glycation End-products form spontaneously and accumulate with aging:

- **Pentosidine** 
- **Glucosepane** 
- **MOLD** (Methylglyoxal-lysine dimer) - Glycation-derived crosslink

### Crosslink Combinations

Crosslinks form at specific positions in the collagen molecule:

**Terminal Crosslinks:**
- **N-terminal**: Typically at positions 9 and 947 (or species-equivalent)
- **C-terminal**: Typically at positions 1047 and 104 (or species-equivalent)

**Divalent combinations** (2 residues):
```yaml
n_term_combination: "9.C - 947.A"
c_term_combination: "1047.C - 104.C"
```

**Trivalent combinations** (3 residues):
```yaml
# PYD crosslink example
n_term_combination: "6.B - 9.C - 946.A"
c_term_combination: "1046.C - 1046.A - 103.C"
```

**Note**: Residue positions may vary slightly between species. See [crosslinks.csv](https://github.com/graeter-group/colbuilder/blob/main/src/colbuilder/data/sequence/crosslinks.csv) for species-specific combinations.

---

## 📖 Usage Guide

### Basic Usage

The general syntax for running ColBuilder is:

```bash
colbuilder --config_file config.yaml [OPTIONS]
```

### Configuration Options

ColBuilder uses YAML configuration files to define parameters. Here's a complete template with all available options:

```yaml
# ================================================================================
# ColBuilder Configuration File
# ================================================================================
# The pipeline consists of three main modules that can be run 
# independently or in combination:
#
# 1. SEQUENCE GENERATOR: Creates collagen triple helix structures from sequences
# 2. GEOMETRY GENERATOR: Builds collagen microfibrils from triple helices
# 3. TOPOLOGY GENERATOR: Generates force field topologies for simulations
#
# IMPORTANT: When using additional crosslinks with a pre-mutated PDB, 
# e.g. to add AGE crosslinks:
# - Run sequence generation first to apply and optimize additional crosslinks
# - Then run geometry generation in a separate execution using the output PDB
# ================================================================================

# Operation Mode
mode: null
config_file: null            # Path to another config file (for nested configs)  
sequence_generator: true     # Generate triple helix structure from sequence
geometry_generator: true     # Build microfibril from triple helix
topology_generator: false    # Generate simulation topology files
debug: false                 # Keep intermediate files for debugging
working_directory: "./"      # Directory for output files

# ================================================================================
# SEQUENCE GENERATION SETTINGS
# ================================================================================
# These settings control the generation of collagen triple helix structures.
# Two main workflows are supported:
# 1. Standard: Generate from FASTA sequence (leave mutated_pdb as null)
# 2. Mutated PDB: Add crosslinks to existing structure (must provide mutated_pdb)

# Species selection
species: "homo_sapiens"
# Available species options:
# Mammals (Primates): homo_sapiens, pan_troglodytes, pongo_abelii, callithrix_jacchus, otolemur_garnettii
# Mammals (Rodents): mus_musculus, rattus_norvegicus
# Mammals (Other): bos_taurus, canis_lupus, ailuropoda_melanoleuca, mustela_putorius, myotis_lucifugus, loxodonta_africana
# Fish: danio_rerio, oreochromis_niloticus, oryzias_latipes, tetraodon_nigroviridis, xiphophorus_maculatus
# Reptiles: pelodiscus_sinensis

# Mutated PDB workflow - provide pre-optimized structure to add additional crosslinks
mutated_pdb: null  # Path to pre-mutated PDB, or null for standard workflow
# Example: "rattusnorvegicus_N_PYD_C_PYD+ADD1_Glucosepane.pdb"

# Custom sequence input - overrides species selection if provided
fasta_file: null  # Path to custom FASTA file, or null to use available species sequence
# Example: "path/to/custom_sequence.fasta"

# ================================================================================
# CROSSLINK CONFIGURATION
# ================================================================================
# Terminal crosslinks are applied at N and C termini during standard generation.
# Additional crosslinks can be added anywhere (e.g., AGE products).

crosslink: true  # Enable crosslink application and optimization

# Terminal crosslinks (standard sequence generation)
n_term_type: "HLKNL"               # N-terminal crosslink type
c_term_type: "HLKNL"               # C-terminal crosslink type
n_term_combination: "9.C - 947.A"  # Format: "resid.chain - resid.chain"
c_term_combination: "1047.C - 104.C"

# Available crosslink types:
# Enzymatic Divalent: "HLKNL", "LKNL", "deHLNL", "deHHLNL"
# Enzymatic Trivalent: "PYD", "DPD", "PYL", "DPL" (3 residues required)
# Non-Enzymatic (AGE): "Pentosidine", "Glucosepane", "MOLD"

# Trivalent crosslink example (requires 3 residues):
# n_term_combination: "6.B - 9.C - 946.A"
# c_term_combination: "1046.C - 1046.A - 103.C"

# Additional crosslinks (for mutated PDB workflow)
# Uncomment to add crosslinks on top of existing structure:
# additional_1_type: "Pentosidine"
# additional_1_combination: "1008.A - 767.B"
# additional_2_type: null  # Optional second additional crosslink
# additional_2_combination: null

# Crosslink optimization settings
# crosslink_copies: ["D0", "D5"]  # Periodic copies for optimization, must be exactly 2 elements
# Options: D0-D5, D0-D1, D1-D2, D2-D3, D3-D4. If not specified, default is ["D0", "D5"]

# ================================================================================
# GEOMETRY GENERATION SETTINGS
# ================================================================================
# These settings control microfibril assembly from triple helices.
# The geometry generator arranges multiple copies of the input structure
# into a collagen microfibril based on crystal contacts.

# Input structure for geometry generation
pdb_file: null  # Path to input PDB, or null if sequence_generator is true
# IMPORTANT: When using additional crosslinks, run sequence generation 
# separately first, then use its output here

# Microfibril parameters
contact_distance: 20      # Maximum distance (Å) for crystal contacts
fibril_length: 70.0       # Target length of microfibril (nm)

# Advanced geometry options (usually left as null/false)
crystalcontacts_file: null      # Custom crystal contacts file
connect_file: null              # Custom connectivity information
crystalcontacts_optimize: false # Optimize crystal contacts
# Note: crystalcontacts_optimize is automatically true when using contact_distance

# ================================================================================
# MIXING OPTIONS
# ================================================================================
# Create heterogeneous microfibrils with different crosslink types.

mix_bool: false           # Enable mixing of different crosslink types
ratio_mix: "D:80 T:20"    # Ratio of each type (must sum to 100)
# Format: "TypeLabel:percentage TypeLabel:percentage"
# Example: "D:80 T:20" = 80% divalent, 20% trivalent
files_mix:  # PDB files for each crosslink type (required if mix_bool=true)
  - "collagen-molecule-crosslinkA.pdb"  # File for first type
  - "collagen-molecule-crosslinkB.pdb"  # File for second type
# Note: Labels (D, T, etc.) are arbitrary but must match ratio_mix

# ================================================================================
# REPLACEMENT OPTIONS
# ================================================================================
# Reduce crosslink density by replacing some with standard amino acids.
# Useful for modeling aged or diseased tissue with fewer crosslinks.

replace_bool: false      # Enable crosslink replacement
ratio_replace: 30        # Percentage of crosslinks to replace (0-100)
replace_file: null       # Input file, or null to use geometry output

# ================================================================================
# TOPOLOGY GENERATION SETTINGS
# ================================================================================
# Generate molecular dynamics topology files for the final structure.

force_field: "amber99"  # Force field to use
# Options: "amber99" (all-atom), "martini3" (coarse-grained)

topology_debug: false   # Save intermediate topology files

# Topology-Only Mode
# When topology_generator: true AND both sequence_generator: false and 
# geometry_generator: false, ColBuilder operates in topology-only mode,
# reading an existing fibril PDB and generating topology files without 
# running geometry generation. Crosslink parameters must match the input PDB.
```

For a complete list of configuration options, see the [detailed documentation](https://github.com/graeter-group/colbuilder/tree/main/docs).

### Example Workflows

#### Creating a Basic Human Collagen Microfibril

```yaml
# config_human_basic.yaml
species: "homo_sapiens"
sequence_generator: true
geometry_generator: true
crosslink: false
fibril_length: 40.0
contact_distance: 25
```

```bash
colbuilder --config_file config_human_basic.yaml
```

#### Generating a Crosslinked Bovine Microfibril

```yaml
# config_bovine_crosslinked.yaml
species: "bos_taurus"
sequence_generator: true
geometry_generator: true
crosslink: true
n_term_type: "HLKNL"
c_term_type: "HLKNL"
n_term_combination: "9.C - 946.A"
c_term_combination: "1046.C - 103.C"
fibril_length: 80.0
contact_distance: 15
```

```bash
colbuilder --config_file config_bovine_crosslinked.yaml
```

#### Creating a Mixed Crosslinked (80% Divalent + 20% Trivalent) Human Collagen Microfibril from Collagen Molecules

First, generate the two crosslink types:

```yaml
# config_divalent_molecule.yaml
species: "homo_sapiens"
sequence_generator: true
geometry_generator: false
crosslink: true
n_term_type: "HLKNL"
c_term_type: "HLKNL"
n_term_combination: "9.C - 947.A"
c_term_combination: "1047.C - 104.C"
```

```yaml
# config_trivalent_molecule.yaml
species: "homo_sapiens"
sequence_generator: true
geometry_generator: false
crosslink: true
n_term_type: "PYD"
c_term_type: "PYD"
n_term_combination: "6.B - 9.C - 946.A"
c_term_combination: "1046.C - 1046.A - 103.C"
```

Then mix them:

```yaml
# config_mixed_crosslinks.yaml
sequence_generator: false
geometry_generator: true
contact_distance: 25
fibril_length: 40.0
mix_bool: true
ratio_mix: "D:80 T:20"
files_mix:
 - "homo_sapiens_N_HLKNL_C_HLKNL.pdb"
 - "homo_sapiens_N_PYD_C_PYD.pdb"
```

```bash
colbuilder --config_file config_divalent_molecule.yaml
colbuilder --config_file config_trivalent_molecule.yaml
colbuilder --config_file config_mixed_crosslinks.yaml
```

#### Generating a Coarse-Grained Topology File for MD Simulation

```yaml
# config_topology.yaml
species: "homo_sapiens"
sequence_generator: false
geometry_generator: true
topology_generator: true
pdb_file: "path/to/template_collagen_molecule.pdb"
contact_distance: 30
fibril_length: 40.0
force_field: "martini3"
```

```bash
colbuilder --config_file config_topology.yaml
```

#### Generating Topology Files from an Existing Fibril (Topology-Only Mode)

If you already have a complete fibril PDB file (e.g., from a previous ColBuilder run or from another source) and only need to generate topology files for simulation:

```yaml
# config_topology_only.yaml
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

```bash
colbuilder --config_file config_topology_only.yaml
```

**Important notes for topology-only mode:**
- The input PDB must be a complete fibril (not a single triple helix)
- Crosslink information must match what's present in the PDB file
- The tool will automatically detect model connectivity from crosslink positions
- Works with fibrils containing divalent (D), trivalent (T), or mixed (M) crosslinks
- Supports both Amber99 and Martini3 force fields

#### Creating a Fibril with PYD Crosslinks

PYD (Pyridinoline) is a trivalent crosslink requiring three residues:

```yaml
# config_pyd_fibril.yaml
species: "rattus_norvegicus"
sequence_generator: true
geometry_generator: true
crosslink: true
n_term_type: "PYD"
c_term_type: "PYD"
n_term_combination: "6.B - 9.C - 946.A"
c_term_combination: "1046.C - 1046.A - 103.C"
fibril_length: 40.0
contact_distance: 20
```

```bash
colbuilder --config_file config_pyd_fibril.yaml
```

#### Replacing Crosslinks in an Existing Fibril

Reduce crosslink density by replacing a percentage with lysine:

```yaml
# config_replace_crosslinks.yaml
species: "homo_sapiens"
sequence_generator: false
geometry_generator: true
pdb_file: "homo_sapiens_N_HLKNL_C_HLKNL.pdb"
crosslink: true
n_term_type: "HLKNL"
c_term_type: "HLKNL"
n_term_combination: "9.C - 947.A"
c_term_combination: "1047.C - 104.C"
replace_bool: true
ratio_replace: 30  # Replace 30% of crosslinks
fibril_length: 40.0
contact_distance: 20
```

```bash
colbuilder --config_file config_replace_crosslinks.yaml
```

#### Adding AGE Crosslinks to an Existing Structure (Mutated PDB Workflow)

Add non-enzymatic crosslinks (AGEs) on top of existing enzymatic crosslinks. Examples for this workflow can be found in [docs/examples/example5](https://github.com/graeter-group/colbuilder/tree/main/docs/examples/example5):

- Sequence + Geometry + Topology generation for a collagen with a mixture of PYD and Glucosepane crosslinks [example5.1](https://github.com/graeter-group/colbuilder/tree/main/docs/examples/example5/example5.1)
- Sequence + Geometry + Replacement generation for a collagen with a mixture PYD and Glucosepane crosslinks [example5.2](https://github.com/graeter-group/colbuilder/tree/main/docs/examples/example5/example5.2)
- Sequence + Geometry + Mixing + Topology generation for different crosslink combinations [example5.3](https://github.com/graeter-group/colbuilder/tree/main/docs/examples/example5/example5.3)

In general, the most basic workflow happens in three steps:

**Step 1:** Generate base structure with enzymatic crosslinks

```yaml
# config_add_age_step1.yaml
species: "rattus_norvegicus"
sequence_generator: true
geometry_generator: false
topology_generator: false
crosslink: true

# Base terminal crosslinks
n_term_type: "PYD"
c_term_type: "PYD"
n_term_combination: "6.B - 9.C - 946.A"
c_term_combination: "1046.C - 1046.A - 103.C"
```

**Step 2**: Add AGE crosslinks to the output of step 1:

```yaml
# config_add_age_step2.yaml
species: "rattus_norvegicus"
sequence_generator: true
geometry_generator: false
topology_generator: false

mutated_pdb: "rattusnorvegicus_N_PYD_C_PYD.pdb"

crosslink: true
# Terminal crosslinks must match step 1
n_term_type: "PYD"
c_term_type: "PYD"
n_term_combination: "6.B - 9.C - 946.A"
c_term_combination: "1046.C - 1046.A - 103.C"

# Add Glucosepane AGE crosslink
additional_1_type: "Glucosepane"
additional_1_combination: "1055.C - 822.A"

# Customize crosslink optimization
crosslink_copies: ["D2", "D3"]
```

**Step 3:** Build fibril using the output from step 2:

```yaml
# config_add_age_step3.yaml
species: "rattus_norvegicus"
sequence_generator: false
geometry_generator: true
topology_generator: true

crosslink: true
n_term_type: "PYD"
c_term_type: "PYD"
n_term_combination: "6.B - 9.C - 946.A"
c_term_combination: "1046.C - 1046.A - 103.C"

# Use output from step 2
pdb_file: "rattusnorvegicus_N_PYD_C_PYD+ADD1_Glucosepane.pdb"

# Fibril parameters
fibril_length: 40.0
contact_distance: 20
force_field: "amber99"
```

```bash
# Run in three steps
colbuilder --config_file config_add_age_step1.yaml
colbuilder --config_file config_add_age_step2.yaml
colbuilder --config_file config_add_age_step3.yaml
```

**Why three steps?** The mutated PDB workflow requires separate execution to properly optimize the additional crosslinks before building the fibril geometry.

## 📚 Documentation

For detailed API documentation, advanced usage examples, and theoretical background:

- [User Guide](https://github.com/graeter-group/colbuilder/tree/main/docs/user_guide.md)
- [Configuration Reference](https://github.com/graeter-group/colbuilder/tree/main/docs/configuration.md)
- [Data Reference](https://github.com/graeter-group/colbuilder/tree/main/docs/data.md)
- [Data Dictionary](https://github.com/graeter-group/colbuilder/tree/main/docs/data_dictionary.md)
- [Example Gallery](https://github.com/graeter-group/colbuilder/tree/main/docs/examples)

## 🤝 Contributing

We welcome contributions to ColBuilder! Please see our [contributing guidelines](https://github.com/graeter-group/colbuilder/tree/main/CONTRIBUTING.md) for details on how to submit issues, pull requests, and code reviews.

## 📚 Publications & Citation

If you use ColBuilder in your research, please cite our paper:

https://academic.oup.com/bioinformatics/article/41/6/btaf278/8125020

A BibTeX entry is provided in the [CITATION.bibtex](https://github.com/graeter-group/colbuilder/tree/main/CITATION.bibtex) file.

If you perform coarse-grained simulations, please also cite the Martini 3 paper for collagen fibrils:

https://www.cell.com/biophysj/fulltext/S0006-3495(25)00663-0

## 🙏 Acknowledgements

ColBuilder is developed and maintained by the Gräter group at the Max Planck Institute for Polymer Research. We thank all contributors that have supported this work.

---

For questions, feedback, or support, please [open an issue](https://github.com/graeter-group/colbuilder/issues) on our GitHub repository.