<img src="https://github.com/user-attachments/assets/ea9609dc-9b11-4335-8d9d-a4d1bc9298cb" width="1500">

<div align="center">
    <h1>ColBuilder 2.0</h1>
    <p>Generate microfibrils from single collagen molecules</p>
</div>

---

[![conda-forge](https://img.shields.io/badge/conda-forge-orange.svg)]()
[![Documentation](https://img.shields.io/badge/docs-brightgreen.svg)](https://github.com/graeter-group/colbuilder/tree/main/docs)
[![license](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](http://www.apache.org/licenses/LICENSE-2.0)

# 📚 About

**ColBuilder** is a tool for generating atomistic models of collagen microfibrils from single collagen molecules. It offers flexibility in input parameters, structural optimization, and supports various crosslink types to mimic biological systems.

![ColBuilder Workflow](https://github.com/user-attachments/assets/ccbdfbdf-a78e-407c-899d-4ecc67f7cd07)

## Key Features

- Generate microfibrils from single collagen molecules
- Flexible input parameters (collagen sequence, fibril geometry, crosslink type and density)
- Structural optimization using UCSF Chimera
- Support for different crosslink types
- Generation of atomistic and coarse-grained topology files

The `colbuilder` tool was created by the Gräter group at the Max Planck Institute for Polymer Research, and will be regularly maintained.
If you find ColBuilder useful, please see the citation file for details on how to cite.

# 🚀 Getting Started

## Install ColBuilder

**Use a virtual environment for the ColBuilder project.** We recommend the [miniforge](https://github.com/conda-forge/miniforge) environment manager.

    conda create -n colbuilder python=3.9
    conda activate colbuilder

**Clone this repository and install ColBuilder:**

    git clone git@github.com:graeter-group/colbuilder.git
    cd colbuilder
    pip install .

**Install additional dependencies:**

A few external tools are required to run ColBuilder. You can download them according to the following instructions.

### PyMOL

    conda install conda-forge::pymol-open-source

**NOTE**: pymol cannot be called from python if the libnetcdf.so library is missing. If running `colbuilder --help` results in an error, we recommend installing the libnetcdf.so package:

    conda install -c conda-forge libnetcdf==4.7.3

### muscle

    conda install muscle

### UCSF Chimera

Download the latest version of [UCSF Chimera](https://www.cgl.ucsf.edu/chimera/download.html) for your OS. We recommend using the 64-bit version. After downloading the binary, navigate to the download directory, make the binary executable and install it. 

    cd ~/Downloads # or where you downloaded the chimera binary
    chmod +x chimera*.bin`
    ./chimera*.bin

The installer will ask you where you wish to create the symlink for UCSF Chimera. We recommend choosing a symlink that is already present in your `$PATH`, but if you choose not to do so, remember to make sure ColBuilder can find it.

**NOTE**: ColBuilder uses UCSF Chimera instead of the most recent version of the software, UCSF ChimeraX. 

### Modeller

Download [Modeller version 10.5](https://salilab.org/modeller/download_installation.html) for your OS and follow the INSTALLATION instructions. For Ubuntu users, we recommend downloading the Generic Unix Tarball version. During installation you will be asked where to install modeller. Please make sure you add this location and the following PYTHONPATH environment variable values to your bashrc file:

    export PYTHONPATH="/home/user/bin/modeller10.5/lib/x86_64-intel8/python3.3:$PYTHONPATH"
    export PYTHONPATH="/home/user/bin/modeller10.5/modlib:$PYTHONPATH"
    export LD_LIBRARY_PATH="/home/user/bin/modeller10.5/lib/x86_64-intel8:$LD_LIBRARY_PATH"

# 📖 Tutorial: The Collagen Microfibril

After compiling the `colbuilder` Command Line Tool, edit the configuration file and run ColBuilder with:

    colbuilder --config_file config.yaml

# 📊 Examples & Documentation

For more detailed examples and API documentation, please visit our [documentation folder](https://github.com/graeter-group/colbuilder/tree/main/docs).

# 📚 Publications

If you use ColBuilder 2.0 in your research, please cite our paper:

# 🙏 Acknowledgements

We would like to thank...
