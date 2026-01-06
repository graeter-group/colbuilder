# ColBuilder User Guide

## Table of Contents

1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Getting Started](#getting-started)
4. [Basic Workflows](#basic-workflows)
5. [Configuration Options](#configuration-options)
6. [Working with Different Species](#working-with-different-species)
7. [Customizing Collagen Structures](#customizing-collagen-structures)
8. [Advanced Features](#advanced-features)
9. [Troubleshooting](#troubleshooting)
10. [Best Practices](#best-practices)
11. [Next Steps](#next-steps)

## Introduction

ColBuilder is a specialized tool for generating atomistic models of collagen microfibrils from single collagen molecules. It offers a comprehensive pipeline for creating detailed, biologically relevant collagen structures for use in molecular dynamics simulations and structural studies.

### What is ColBuilder?

ColBuilder enables researchers to:
- Generate atomistic models of collagen microfibrils with precise control over their geometric features
- Customize collagen structures with different crosslink types and configurations
- Create models for various species (human, bovine, rodent, and more)
- Prepare structures for molecular dynamics simulations

### Pipeline Overview

ColBuilder implements a three-stage pipeline:

1. **Sequence Generation (Homology Modeling)**: Creates a collagen triple helix based on sequence information
2. **Geometry Generation**: Builds a microfibril structure with specified dimensions and properties
3. **Topology Generation**: Prepares files for molecular dynamics simulations

Each stage can be run independently or as a complete pipeline, allowing flexibility in your workflow.
In addition, ColBuilder supports adding non-enzymatic crosslinks on top of a
pre-mutated triple helix and mixing multiple PDB variants before topology
generation.

## Installation

### Prerequisites

Before installing ColBuilder, ensure you have the following external tools installed:

1. **Python 3.9 or later**
2. **Conda package manager** (we recommend [miniforge](https://github.com/conda-forge/miniforge))
3. **External dependencies:**
   - PyMOL - For visualization and some structural operations
   - MUSCLE - For sequence alignment
   - UCSF Chimera - For structural optimization
   - Modeller (version 10.5) - For homology modeling
   - GROMACS (optional) - For topology generation and molecular dynamics

### Installation Steps

1. **Create a virtual environment:**
   ```bash
   conda create -n colbuilder python=3.9
   conda activate colbuilder
   ```

2. **Clone the repository:**
   ```bash
   git clone git@github.com:graeter-group/colbuilder.git
   cd colbuilder
   ```

3. **Install ColBuilder:**
   ```bash
   pip install .
   ```

4. **Install PyMOL:**
   ```bash
   conda install conda-forge::pymol-open-source
   ```
   
   If you encounter issues with missing libraries:
   ```bash
   conda install -c conda-forge libnetcdf==4.7.3
   ```

5. **Install MUSCLE:**
   ```bash
   conda install muscle
   ```

6. **Install UCSF Chimera:**
   
   Download from [UCSF Chimera website](https://www.cgl.ucsf.edu/chimera/download.html) and follow installation instructions.
   
   Make the binary executable and run the installer:
   ```bash
   chmod +x chimera*.bin
   ./chimera*.bin
   ```
   Ensure Chimera is added to your PATH.

7. **Install Modeller 10.5:**
   
   Download from [Modeller website](https://salilab.org/modeller/download_installation.html) and follow installation instructions.
   
   Add required environment variables to your .bashrc file:
   ```bash
   export PYTHONPATH="/path/to/modeller10.5/lib/x86_64-intel8/python3.3:$PYTHONPATH"
   export PYTHONPATH="/path/to/modeller10.5/modlib:$PYTHONPATH"
   export LD_LIBRARY_PATH="/path/to/modeller10.5/lib/x86_64-intel8:$LD_LIBRARY_PATH"
   ```

8. **Verify installation:**
   ```bash
   colbuilder --help
   ```

## Getting Started

### Quick Start

Create a basic configuration file named `config.yaml`:

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

Run ColBuilder:
```bash
colbuilder --config_file config.yaml
```

This command will generate a collagen microfibril for human (homo_sapiens) collagen with HLKNL crosslinks.

### Understanding the Output

After successful execution, ColBuilder will generate several output files:

- **Sequence Generation Output:**
  - `{species}_alignment.fasta`: Multiple sequence alignment file
  - `{species}_N_{n_term_type}_C_{c_term_type}.pdb`: PDB file of the collagen triple helix

- **Geometry Generation Output:**
  - `collagen_fibril_{species}.pdb`: PDB file of the collagen microfibril

- **Topology Generation Output (if enabled):**
  - `{species}_topology_files/`: Directory containing topology files
  - `collagen_fibril_{species}.top`: Topology file
  - `collagen_fibril_{species}.gro`: GROMACS coordinate file

## Basic Workflows

### 1. Full Pipeline: Sequence to Topology

This workflow runs the complete pipeline from sequence to topology:

```yaml
# config_full.yaml
species: "homo_sapiens"
sequence_generator: true
geometry_generator: true
topology_generator: true
crosslink: true
n_term_type: "HLKNL"
c_term_type: "HLKNL"
n_term_combination: "9.C - 947.A"
c_term_combination: "1047.C - 104.C"
fibril_length: 60.0
contact_distance: 20
force_field: "amber99"
```

```bash
colbuilder --config_file config_full.yaml
```

### 2. Only Sequence Generation

This workflow only generates the collagen triple helix:

```yaml
# config_sequence.yaml
species: "homo_sapiens"
sequence_generator: true
geometry_generator: false
topology_generator: false
crosslink: true
n_term_type: "HLKNL"
c_term_type: "HLKNL"
n_term_combination: "9.C - 947.A"
c_term_combination: "1047.C - 104.C"
```

```bash
colbuilder --config_file config_sequence.yaml
```

### 3. Add Non-Enzymatic Crosslinks (Mutated-PDB Workflow)

This workflow starts from a pre-mutated triple-helix PDB (typically produced
by a standard enzymatic sequence generation step) and adds one or more
non-enzymatic crosslinks (AGEs).

```yaml
# config_non_enzymatic.yaml
sequence_generator: true
geometry_generator: false
topology_generator: false
species: "rattus_norvegicus"
mutated_pdb: "rattusnorvegicus_N_PYD_C_PYD.pdb"
crosslink: true
n_term_type: "PYD"
c_term_type: "PYD"
n_term_combination: "9.C - 5.B - 944.B"
c_term_combination: "1046.C - 1046.A - 103.C"
additional_1_type: "Glucosepane"
additional_1_combination: "523.A - 286.C"
crosslink_copies: ["D2", "D1"]
```

```bash
colbuilder --config_file config_non_enzymatic.yaml
```

This produces a new triple-helix PDB with the additional non-enzymatic
crosslink(s), which can then be used for geometry generation. In this second
sequence-generation step, set `crosslink_copies` to the required shift listed in
`src/colbuilder/data/sequence/crosslinks.csv`. Enzymatic crosslinks use D0-D5 by default.

### 4. Only Geometry Generation

This workflow generates the microfibril geometry from an existing PDB:

```yaml
# config_geometry.yaml
species: "homo_sapiens"
sequence_generator: false
geometry_generator: true
topology_generator: false
pdb_file: "path/to/triple_helix.pdb"
fibril_length: 60.0
contact_distance: 20
```

```bash
colbuilder --config_file config_geometry.yaml
```

### 5. Mix PDB Variants and Generate Topology

Mixing can be used to combine multiple crosslink configurations into a single
microfibril. Topology generation can follow directly after mixing.

```yaml
# config_mix_topology.yaml
sequence_generator: false
geometry_generator: false
mix_bool: true
topology_generator: true
ratio_mix: "A:80 B:20"
files_mix:
 - "variant_A.pdb"
 - "variant_B.pdb"
force_field: "amber99"
```

```bash
colbuilder --config_file config_mix_topology.yaml
```

## Configuration Options

ColBuilder is highly configurable through its YAML configuration file. The following sections detail the available options.

### General Options

```yaml
debug: false                # Enable detailed debug logging
working_directory: "./"     # Working directory for input and output files
```

### Sequence Generation Options

```yaml
sequence_generator: true    # Enable sequence generation
species: "homo_sapiens"     # Species name
fasta_file: null            # Custom FASTA file (if null, uses built-in species data)
mutated_pdb: null           # Use a pre-mutated PDB as input for additional crosslinks
crosslink: true             # Enable crosslinking
n_term_type: "HLKNL"        # N-terminal crosslink type
c_term_type: "HLKNL"        # C-terminal crosslink type
n_term_combination: "9.C - 947.A"  # N-terminal residue combination
c_term_combination: "1047.C - 104.C" # C-terminal residue combination
additional_1_type: null     # Optional non-enzymatic crosslink (e.g., Glucosepane)
additional_1_combination: null
additional_2_type: null     # Optional second non-enzymatic crosslink
additional_2_combination: null
crosslink_copies: ["D0", "D5"]  # Translation pair for crosslink optimization
# Non-enzymatic additions should use the shift listed in crosslinks.csv for that crosslink.
# Enzymatic crosslinks default to D0-D5 unless specified otherwise.
```

### Geometry Generation Options

```yaml
geometry_generator: true    # Enable geometry generation
pdb_file: null              # Input PDB file (if sequence_generator is false)
contact_distance: 20        # Distance threshold for contacts (Å)
fibril_length: 60.0         # Length of the generated fibril (nm)
crystalcontacts_file: null  # File with crystal contacts info
connect_file: null          # File with connections between contacts
crystalcontacts_optimize: false # Optimize crystal contacts
```

### Mixing Options

```yaml
mix_bool: false             # Enable mixing of different crosslink types
ratio_mix: "D:70 T:30"      # Ratio for different crosslink types
files_mix:                  # Required if mix_bool is true
 - "human-D.pdb"            # PDB file with type D crosslinks
 - "human-T.pdb"            # PDB file with type T crosslinks
```

When mixing is enabled, topology generation can also be enabled in the same
run by setting `topology_generator: true`.

### Replacement Options

```yaml
replace_bool: false         # Enable replacement of crosslinks with lysines
ratio_replace: 30           # Percentage of crosslinks to replace
replace_file: null          # File with crosslinks to be replaced
```

### Topology Options

```yaml
topology_generator: false   # Enable topology generation
force_field: "amber99"      # Force field to use (amber99 or martini3)
```

For a complete reference of configuration options, see the [Configuration Reference](configuration_reference.md).

## Working with Different Species

ColBuilder supports generating collagen structures for a variety of species:

### Supported Species

- **Mammals (Primates)**: homo_sapiens, pan_troglodytes, pongo_abelii, callithrix_jacchus, otolemur_garnettii
- **Mammals (Rodents)**: mus_musculus, rattus_norvegicus
- **Mammals (Other)**: bos_taurus, canis_lupus, ailuropoda_melanoleuca, mustela_putorius, myotis_lucifugus, loxodonta_africana
- **Fish**: danio_rerio, oreochromis_niloticus, oryzias_latipes, tetraodon_nigroviridis, xiphophorus_maculatus
- **Reptiles**: pelodiscus_sinensis

### Species-Specific Configuration

Different species have different residue combinations available for crosslinking. When configuring crosslinks, ensure you use valid combinations for your chosen species. Here are some examples:

#### Human (homo_sapiens)

```yaml
species: "homo_sapiens"
n_term_type: "HLKNL"
c_term_type: "HLKNL"
n_term_combination: "9.C - 947.A"
c_term_combination: "1047.C - 104.C"
```

#### Bovine (bos_taurus)

```yaml
species: "bos_taurus"
n_term_type: "HLKNL"
c_term_type: "HLKNL"
n_term_combination: "9.C - 946.A"
c_term_combination: "1046.C - 103.C"
```

#### Using Custom Sequences

For species not in the predefined list, provide a custom FASTA file:

```yaml
species: "custom_species"
fasta_file: "path/to/custom_sequence.fasta"
sequence_generator: true
crosslink: false
# Other parameters as needed
```

Crosslinks definition for custom_species must be added to the file in src/colbuilder/data/sequence/crosslinks.csv

## Customizing Collagen Structures

### Crosslink Types

ColBuilder supports various crosslink types found in collagen. The complete,
species-specific list is maintained in
`src/colbuilder/data/sequence/crosslinks.csv`; common examples include:

- **HLKNL**: Hydroxylysino-5-ketonorleucine (divalent crosslink)
- **LKNL**: Lysino-5-ketonorleucine (divalent)
- **PYD**: Pyridinoline (trivalent)
- **DPD**: Deoxypyridinoline (trivalent)
- **PYL**: Pyrrole crosslink (trivalent)
- **DPL**: (trivalent)
- **deHHLNL**: Dehydro-hydroxylysino-norleucine (divalent)
- **deHLNL**: Dehydro-lysino-norleucine (divalent)
- **NOCROSS**: No crosslinking
- **Glucosepane**: Advanced glycation end-product crosslink (specific species)
- **Pentosidine**: Advanced glycation end-product crosslink (specific species)
- **MOLD**: Non-enzymatic, LYS-LYS derived crosslink

Each crosslink type can be positioned in a few combination of selected Lysine residues, depending on the species being modeled. All available crosslinks and respective combinations for each species are listed at [src/colbuilder/data/sequence/crosslinks.csv](https://github.com/graeter-group/colbuilder/blob/main/src/colbuilder/data/sequence/crosslinks.csv)

### Fibril Dimensions

Control the size of your microfibril using these parameters:

```yaml
contact_distance: 20        # Controls radial size (larger values produce thicker fibrils)
fibril_length: 60.0         # Controls fibril length in nanometers
```

### Optimizing Crystal Contacts

For better-packed microfibrils, enable crystal contact optimization:

```yaml
crystalcontacts_optimize: true
solution_space: [1, 1, 1]   # Solution space for optimization [d_x, d_y, d_z]
```

## Advanced Features

### Mixing Crosslink Types

To create a heterogeneous microfibril with different crosslink types:

```yaml
mix_bool: true
ratio_mix: "D:70 T:30"      # 70% type D, 30% type T
files_mix:
 - "human-D.pdb"            # triple helix PDB with type D crosslinks
 - "human-T.pdb"            # triple helix PDB with type T crosslinks
```

This feature allows modeling of collagen structures with a mixture of crosslink
types and can be chained directly into topology generation.

### Replacing Crosslinks

To model partially crosslinked collagen (replacing some crosslinks with standard lysines):

```yaml
replace_bool: true
ratio_replace: 30                   # Replace 30% of crosslinks with lysines
replace_file: "original_fibril.pdb" # PDB with original collagen microfibril structure
```

### Using Different Force Fields

For molecular dynamics simulations, ColBuilder supports multiple force fields:

```yaml
topology_generator: true
force_field: "amber99"      # For all-atom simulations
```

```yaml
topology_generator: true
force_field: "martini3"     # For coarse-grained simulations
```

## Troubleshooting

### Common Issues

#### 1. PyMOL Issues

**Error:** "PyMOL not found" or library-related errors

**Solution:**
```bash
conda install conda-forge::pymol-open-source
conda install -c conda-forge libnetcdf==4.7.3
```

#### 2. UCSF Chimera Issues

**Error:** "Chimera not found" or scripts not running

**Solutions:**
- Ensure UCSF Chimera (not ChimeraX) is installed
- Add Chimera to your PATH
- Check that Chimera scripts are in the expected directory

#### 3. Modeller Issues

**Error:** Import errors with Modeller

**Solutions:**
- Verify environment variables are set correctly
- Ensure you're using Modeller 10.5 specifically
- Check the Modeller license key is properly installed

#### 4. Sequence Generation Errors

**Error:** "Sequence generation failed" or "Crosslink application failed"

**Solutions:**
- Verify species name is correctly spelled
- Ensure crosslink combinations are valid for the species
- Check that fasta_file path is correct if using a custom sequence

#### 5. Geometry Generation Errors

**Error:** "Geometry generation failed" or "Invalid contact distance"

**Solutions:**
- Ensure contact_distance is positive and reasonable (typically 15-25)
- Verify the input PDB file exists and is valid
- Check fibril_length is within reasonable limits (typically 30-334)

#### 6. Topology Generation Errors

**Error:** "Topology generation failed" or "Invalid force field"

**Solutions:**
- Verify force_field is either "amber99" or "martini3"
- Ensure GROMACS is installed if using amber99
- Check that Martinize2 is installed if using martini3

### Debugging

To enable detailed logging:

```bash
colbuilder --config_file config.yaml --debug
```

This will provide more information about what's happening during execution and save all intermediate files.

## Best Practices

### Performance Considerations

- **Memory Usage**: Larger - especially in diameter - fibrils require significant memory. Start with smaller fibril_length values (contact_distance: 20-60 Ang) before attempting larger structures.
- **Processing Time**: Sequence generation with crosslink optimization is one of the most time-consuming steps. Consider running this step separately and reusing the output for multiple geometry variations.

### Workflow Recommendations

1. **Iterative Approach**: Start with sequence generation only, validate the output, then proceed to geometry and topology generation.
   
2. **Validation at Each Step**: Visualize the output of each stage with molecular visualization software like PyMOL or VMD to ensure correctness.

3. **Parameter Exploration**: Test different contact_distance values to find the optimal radial packing for your system.

## Next Steps

After generating your collagen microfibril:

### Visualization

Visualize your structures using molecular visualization software:

```bash
vmd collagen_fibril_homo_sapiens.pdb
```

### Molecular Dynamics Simulations

If you've generated topology files, you can proceed to molecular dynamics simulations using GROMACS:

1. **Energy Minimization**
2. **Equilibration**
3. **Production Runs**

Refer to the GROMACS documentation for detailed instructions.

### Analysis

Analyze your collagen structure:
- Crosslink distribution and geometry
- Fibril packing and density
- Mechanical properties through simulations

### Extensions

ColBuilder can be extended to:
- Study the effects of mutations on collagen structure
- Model disease-relevant modifications
- Investigate mechanical properties through simulations
- Study interactions with other extracellular matrix components

For more information, see the [Configuration Reference](configuration.md) and the [Data Dictionary](data_dictionary.md).
