<div align="center">
    <h1>ColBuilder</h1>
    <p>Generate atomistic models of collagen microfibrils from single collagen molecules</p>
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
- [📖 Usage Guide](#-usage-guide)
  - [Basic Usage](#basic-usage)
  - [Configuration Options](#configuration-options)
  - [Example Workflows](#example-workflows)
    - [Creating a Basic Human Collagen Microfibril](#creating-a-basic-human-collagen-microfibril)
    - [Generating a Crosslinked Bovine Microfibril](#generating-a-crosslinked-bovine-microfibril)
    - [Creating a Mixed Crosslinked (80% Divalent + 20% Trivalent) Human Collagen Microfibril](#creating-a-mixed-crosslinked-80-divalent--20-trivalent-human-collagen-microfibril)
    - [Generating a Coarse-Grained Topology File for MD Simulation](#generating-a-coarse-grained-topology-file-for-md-simulation)
- [📚 Documentation](#-documentation)
- [🤝 Contributing](#-contributing)
- [📚 Publications \& Citation](#-publications--citation)
- [🙏 Acknowledgements](#-acknowledgements)

## 📚 About

**ColBuilder** is a specialized tool for generating atomistic models of collagen microfibrils from single collagen molecules. Developed by the Gräter group at the Max Planck Institute for Polymer Research, it provides researchers with a flexible framework to create biologically relevant collagen structures for molecular dynamics simulations and structural studies.

### Key Features

- **Custom microfibril generation**: Create collagen microfibrils from individual molecules or amino acid sequences with precise control over structural parameters
- **Highly configurable**: Adjust collagen sequence, fibril geometry, crosslink types and density to match your custom conditions
- **Simulation-ready output**: Generate atomistic and coarse-grained topology files compatible with major molecular dynamics packages
- **Reproducible research**: Standardized approach to collagen modeling to ensure consistency across studies

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
conda install muscle
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
   working_directory: "./"     
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

## 📖 Usage Guide

### Basic Usage

The general syntax for running ColBuilder is:

```bash
colbuilder --config_file config.yaml [OPTIONS]
```

### Configuration Options

ColBuilder uses YAML configuration files to define parameters. Here's a complete template with all available options:

```yaml
# Operation Mode
mode: null                   # Specific operation mode if needed
config_file: null            # Path to another config file (for nested configs)
sequence_generator: true     # Generate sequence from species
geometry_generator: true     # Generate fibril geometry
topology_generator: false    # Generate topology files
debug: false                 # Enable debug mode
working_directory: "./"      # Working directory for inputs/outputs

# Input Configuration
species: "homo_sapiens"      # Species for collagen sequence
# Available species options:
# Mammals (Primates): homo_sapiens, pan_troglodytes, pongo_abelii, callithrix_jacchus, otolemur_garnettii
# Mammals (Rodents): mus_musculus, rattus_norvegicus
# Mammals (Other): bos_taurus, canis_lupus, ailuropoda_melanoleuca, mustela_putorius, myotis_lucifugus, loxodonta_africana
# Fish: danio_rerio, oreochromis_niloticus, oryzias_latipes, tetraodon_nigroviridis, xiphophorus_maculatus
# Reptiles: pelodiscus_sinensis

# Sequence Settings
fasta_file: null             # Custom FASTA file path (if null, auto-generated based on species)
crosslink: true              # Enable crosslinking in the model
n_term_type: "HLKNL"         # N-terminal crosslink type (Options: "DPD", "DPL", "HLKNL", "LKNL", "PYD", "PYL", "deHHLNL", "deHLNL", "NONE") 
c_term_type: "HLKNL"         # C-terminal crosslink type (Options: "DPD", "DPL", "HLKNL", "LKNL", "PYD", "PYL", "deHHLNL", "deHLNL", "NONE")
n_term_combination: "9.C - 947.A"    # N-terminal residue combination
c_term_combination: "1047.C - 104.C" # C-terminal residue combination

# Geometry Parameters
pdb_file: null               # Input PDB file (set to null if sequence_generator is true)
contact_distance: 20         # Distance threshold for contacts (Å)
fibril_length: 70.0          # Length of the generated fibril (nm)
crystalcontacts_file: null   # File with crystal contacts 
connect_file: null           # File with connection information
crystalcontacts_optimize: false  # Optimize crystal contacts during generation

# Mixing Options (for mixed crosslinked microfibril)
mix_bool: false              # Enable mixing of different crosslink types
ratio_mix: "A:70 B:30"       # Format: "Type:percentage Type:percentage"
files_mix:                   # Required if mix_bool is true
 - "collagen-molecule-crosslinkA.pdb" # PDB file of collagen molecule with type A crosslinks (created by only setting squence and crosslinks = true (please look at the examples))
 - "collagen-molecule-crosslinkB.pdb" # PDB file of collagen molecule with type B crosslinks

# Replacement Options (for fewer crosslinks)
replace_bool: false          # Enable crosslink replacement
ratio_replace: 30            # Percentage of crosslinks to replace
replace_file: null           # File with crosslinks to be replaced (set to null if geometry_generation is true)

# Topology Options
force_field: "amber99"       # Force field for topology generation (Options: "amber99", "martini3")
```

For a complete list of configuration options, see the [detailed documentation](https://github.com/graeter-group/colbuilder/tree/main/docs).

### Example Workflows

#### Creating a Basic Human Collagen Microfibril

```yaml
# config_human_basic.yaml
working_directory: "./"      
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
working_directory: "./"  
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

```yaml
# config_mixed_crosslinks.yaml
species: "homo_sapiens"
sequence_generator: false
geometry_generator: false
mix_bool: true
ratio_mix: "D:80 T:20"
files_mix:
 - "human-D.pdb"
 - "human-T.pdb"
```

```bash
colbuilder --config_file config_mixed_crosslinks.yaml
```

#### Generating a Coarse-Grained Topology File for MD Simulation

```yaml
# config_topology.yaml
species: "homo_sapiens"
sequence_generator: false`
geometry_generator: true
topology_generator: true
pdb_file: "path/to/template_collagen_molecule.pdb"
force_field: "martini3"
```

```bash
colbuilder --config_file config_topology.yaml
```

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

https://www.biorxiv.org/content/10.1101/2024.12.10.627782v1

A BibTeX entry is provided in the [CITATION.cff](https://github.com/graeter-group/colbuilder/tree/main/CITATION.cff) file.

## 🙏 Acknowledgements

ColBuilder is developed and maintained by the Gräter group at the Max Planck Institute for Polymer Research. We thank all contributors that have supported this work.

---

For questions, feedback, or support, please [open an issue](https://github.com/graeter-group/colbuilder/issues) on our GitHub repository.
