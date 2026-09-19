<div align="center">
    <h1>ColBuilder</h1>
    <p>Generate atomistic and coarse-grained models of collagen microfibrils with customizable crosslinking</p>
    <a href="LICENSE.md"><img src="https://img.shields.io/badge/License-Apache%202.0-blue.svg" alt="License: Apache 2.0"></a>
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
  - [✅ Valid Workflows](#-valid-workflows)
- [🔗 Collagen Crosslinks](#-collagen-crosslinks)
  - [Enzymatic Crosslinks](#enzymatic-crosslinks)
  - [Non-Enzymatic Crosslinks](#non-enzymatic-crosslinks)
  - [Crosslink Combinations](#crosslink-combinations)
- [📖 Usage Guide](#-usage-guide)
  - [Basic Usage](#basic-usage)
  - [Configuration Options](#configuration-options)
  - [Example Workflows](#example-workflows)
- [📚 Documentation](#-documentation)
- [🤝 Contributing](#-contributing)
- [📚 Publications \& Citation](#-publications--citation)
- [🙏 Contributors](#-contributors)

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
| 4 | `mix_bool` | Generate a fibril by mixing multiple crosslink types | Two (or more) triple helix PDBs from Mode 1 | Mixed fibril PDB | Yes: with 3 (topology), in the same config; requires the input PDBs from separate Mode-1 runs first |
| 5 | `replace_bool` | Replace crosslinks in an existing fibril | Fibril PDB from Mode 2 or 4 | Modified fibril PDB | Yes: with 2, 3 |

---

### ✅ Valid Workflows

These mode combinations can be run **in a single configuration file**:

- ✅ **`1 + 2`** - Generate molecule and build fibril
- ✅ **`1 + 2 + 3`** - Complete pipeline: molecule → fibril → topology - [example1](docs/examples/example1)
- ✅ **`2 + 3`** - Build fibril and topology from custom triple helix PDB
- ✅ **`1 + 2 + 5`** - Generate, build, then replace crosslinks
- ✅ **`1 + 2 + 5 + 3`** - Generate, build, replace crosslinks with standard residues, then create topology
- ✅ **`2 + 5`** - Build fibril then replace crosslinks - [example3](docs/examples/example3)
- ✅ **`2 + 5 + 3`** - Build fibril, replace crosslinks, then create topology
- ✅ **`3`** - Topology-only mode: generate topology from existing fibril PDB - [example4](docs/examples/example4)
- ✅ **`4 + 3`** - Mix crosslink types and build topology in one config, once the input PDBs exist - [example5.3](docs/examples/example5/example5.3)

Mixing (**Mode 4**) itself needs at least two triple-helix PDBs to mix, so the
full workflow spans multiple config files: generate each crosslink variant
separately (Mode 1), then mix them (optionally with topology generation in
the same config):

```bash
colbuilder --config_file triple_helix_A.yaml
colbuilder --config_file triple_helix_B.yaml
colbuilder --config_file mix_geometry.yaml   # sets mix_bool: true and includes both PDBs
```

See [example2](docs/examples/example2) for a complete, runnable version of this.

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

### Non-Enzymatic Crosslinks

These form spontaneously rather than through enzyme activity, and can be used as terminal crosslinks or added on top of an existing structure (see the [mutated PDB workflow](docs/examples/example5)):

- **MOLD** (Methylglyoxal-lysine dimer) - Divalent, LYS-LYS derived
- **Glucosepane** - Divalent, LYS-ARG derived; most abundant advanced glycation end-product (AGE) in human tissue
- **Pentosidine** - Divalent, LYS-ARG derived; well-characterized AGE crosslink

**Note**: Martini3 crosslink parametrization currently only covers PYD and HLKNL. Other types require `force_field: "amber99"`.

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

ColBuilder uses YAML configuration files to define parameters. Only set what a
given run actually needs — every option has a default, and most configs only
declare a handful of them:

```yaml
species: "bos_taurus"
sequence_generator: true
geometry_generator: true
topology_generator: true

crosslink: true
n_term_type: "PYD"                        # trivalent: 3 residues per terminus
c_term_type: "PYD"
n_term_combination: "5.B - 9.C - 946.A"
c_term_combination: "1046.C - 1046.A - 103.C"

contact_distance: 20      # Å, radial size of the microfibril
fibril_length: 60.0       # nm, length of the microfibril

force_field: "martini3"   # or "amber99"
```

For the **complete reference** of every configuration option (mixing,
replacement, AGE/mutated-PDB workflow, custom species, all defaults and
validation rules), see the [Configuration Reference](docs/configuration.md).

### Example Workflows

[`docs/examples/`](docs/examples/) contains complete, runnable configs for
every workflow below — each with its own README explaining what it does and
what output to expect:

| Example | Workflow |
|---|---|
| [example1](docs/examples/example1) | Complete pipeline in one config: sequence → geometry → amber99 topology, divalent (HLKNL) crosslinks |
| [example2](docs/examples/example2) | Mixing two crosslink types (HLKNL + PYD) into one microfibril by ratio |
| [example3](docs/examples/example3) | Reducing crosslink density on an existing fibril via ratio-based replacement |
| [example4](docs/examples/example4) | Topology-only mode: amber99 and martini3 topology from an existing mixed-crosslink fibril |
| [example5](docs/examples/example5) | Adding non-enzymatic crosslinks (Glucosepane, Pentosidine, MOLD) via the mutated-PDB workflow — combined with replacement (5.2) and with mixing (5.3) |
| [example6](docs/examples/example6) | Complete pipeline with the martini3 (coarse-grained) force field |
| [example7](docs/examples/example7) | Deterministic, named-residue crosslink replacement (`manual_replacements`) |
| [example8](docs/examples/example8) | Bringing your own sequence: custom species via `fasta_file` and the bundled `pdb2fasta` script |

A quick illustration of the AGE crosslink workflow (example5), since it's the
least obvious one — it needs sequence generation run twice, then geometry and
topology from the result:

```yaml
# Step 1: base structure with terminal crosslinks
species: "rattus_norvegicus"
sequence_generator: true
crosslink: true
n_term_type: "PYD"
c_term_type: "PYD"
n_term_combination: "6.B - 9.C - 946.A"
c_term_combination: "1046.C - 1046.A - 103.C"
```

```yaml
# Step 2: add Glucosepane on top of step 1's output
species: "rattus_norvegicus"
sequence_generator: true
mutated_pdb: "rattusnorvegicus_N_PYD_C_PYD.pdb"
crosslink: true
n_term_type: "PYD"           # must match what's already in mutated_pdb
c_term_type: "PYD"
n_term_combination: "6.B - 9.C - 946.A"
c_term_combination: "1046.C - 1046.A - 103.C"
additional_1_type: "Glucosepane"
additional_1_combination: "1055.C - 822.A"
```

```yaml
# Step 3: build fibril + topology from step 2's output
species: "rattus_norvegicus"
geometry_generator: true
topology_generator: true
pdb_file: "rattusnorvegicus_N_PYD_C_PYD+ADD1_Glucosepane.pdb"
crosslink: true
n_term_type: "PYD"
c_term_type: "PYD"
n_term_combination: "6.B - 9.C - 946.A"
c_term_combination: "1046.C - 1046.A - 103.C"
additional_1_type: "Glucosepane"  # must also be declared here: pdb_file above
                                   # already contains it, and crosslink-type
                                   # validation checks against everything
                                   # actually present in the structure
fibril_length: 40.0
contact_distance: 20
force_field: "amber99"
```

```bash
colbuilder --config_file config_step1.yaml
colbuilder --config_file config_step2.yaml
colbuilder --config_file config_step3.yaml
```

## 📚 Documentation

For detailed API documentation, advanced usage examples, and theoretical background:

- [User Guide](docs/user_guide.md)
- [Configuration Reference](docs/configuration.md)
- [Data Reference](docs/data.md)
- [Data Dictionary](docs/data_dictionary.md)
- [Example Gallery](docs/examples)

## 🤝 Contributing

We welcome contributions to ColBuilder! Please see our [contributing guidelines](CONTRIBUTING.md) for details on how to submit issues, pull requests, and code reviews.

## 📚 Publications & Citation

* If you use ColBuilder in your research, please cite our paper: https://academic.oup.com/bioinformatics/article/41/6/btaf278/8125020. A BibTeX entry is provided in the [CITATION.bibtex](CITATION.bibtex) file.

* If you perform coarse-grained simulations, please also cite the Martini 3 paper for collagen fibrils: https://www.cell.com/biophysj/fulltext/S0006-3495(25)00663-0
* If you use AGE crosslinks in colbuilder, then please also cite "Introducing non-enzymatic crosslinks into atomistic simulations of collagen fibrils" https://www.biorxiv.org/content/10.64898/2026.03.13.711566v1

## 🙏 Contributors

ColBuilder is developed and maintained by the Gräter group at the Max Planck Institute for Polymer Research. 
* Debora Monego
* Johanna Buck
* Matthias Brosz
* Guido Giannetti (University of Vienna)
* Justin Pils (University of Vienna)

---

For questions, feedback, or support, please [open an issue](https://github.com/graeter-group/colbuilder/issues) on our GitHub repository.

ColBuilder is licensed under the [Apache License 2.0](LICENSE.md).
