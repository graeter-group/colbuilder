# ColBuilder Data Directory

This directory contains essential data files and templates used by ColBuilder throughout its pipeline. These files provide the foundation for sequence alignment, structural modeling, crosslinking, and molecular dynamics simulations.

## Directory Structure

```
data/
├── sequence/
│   ├── __init__.py
│   ├── crosslinks.csv
│   ├── template.fasta
│   ├── template.pdb
│   ├── fasta_sequences/
│   │   ├── ailuropodamelanoleuca.fasta
│   │   ├── bostaurus.fasta
|   │   ├── ...
│   │   └── xiphophorusmaculatus.fasta
│   └── modeller/
│       ├── par_mod.lib
│       ├── restyp_mod.lib
│       ├── top_heav_mod.lib
│       └── top_mod.lib
│
└── topology/
    ├── __init__.py
    ├── create_goVirt.py
    ├── modifications.mapping
    ├── selectors.py
    ├── amber99/
    │   ├── aminoacids.rtp
    │   └── modifications.ff
    ├── amber99sb-star-ildnp.ff/
    │   ├── aminoacids.arn
    │   ├── aminoacids.rtp
    │   ├── ...
    │   └── watermodels.dat
    ├── contactmap/
    │   ├── chemical_map.c
    │   ├── chemical_map.o
    │   ├── contact_map
    │   ├── contact_map.c
    │   ├── contact_map.h
    │   ├── contact_map.o
    │   ├── LICENSE
    │   ├── makefile
    │   ├── pdb_map.c
    │   ├── pdb_map.o
    │   ├── protein_map.c
    │   └── protein_map.o
    ├── martini300C-ff/
    │   ├── aminoacids.ff
    │   ├── citations.bib
    │   ├── general.ff
    │   ├── modifications.ff
    │   └── small_molecule_martini3.ff
    └── martini300C-mapping/
        ├── ace.amber99.map
        ├── ...
        ├── hyp.amber99.map
        ├── l4y.amber99.map
        ├── l5y.amber99.map
        ├── ly2.amber99.map
        ├── ly3.amber99.map
        ├── lyx.amber99.map
        ├── nme.amber99.map
        └── modifications.amber99.mapping
```

## Data File Descriptions

### Topology Files

The `topology/` directory contains force field files and utilities essential for the topology generation stage of ColBuilder:

#### Utility Files

- `__init__.py`: Python package initialization file
- `create_goVirt.py`: Script for generating Go-like potentials for Martini models
- `modifications.mapping`: Mapping information for residue modifications in Martini models
- `selectors.py`: Python utilities for atom and residue selection in Martini models

#### `amber99/`

Custom AMBER99 modifications for collagen simulations:

- `aminoacids.rtp`: Modified residue topology parameters that include custom collagen crosslinks
- `modifications.ff`: Force field modifications specific to collagen crosslinks

#### `amber99sb-star-ildnp.ff/`

Standard AMBER99SB-STAR-ILDNP force field files used by GROMACS:

- `aminoacids.rtp`, `aminoacids.arn`, `aminoacids.r2b`: Residue topology parameters for standard amino acids
- `aminoacids.c.tdb`, `aminoacids.n.tdb`: Terminal database for C and N termini
- `aminoacids.hdb`: Hydrogen database for amino acids
- `atomtypes.atp`: Atom type definitions
- `ffbonded.itp`: Bonded interaction parameters
- `ffnonbonded.itp`: Non-bonded interaction parameters
- `forcefield.itp`, `forcefield.doc`: Force field description and documentation
- Various solvent models: `spc.itp`, `spce.itp`, `tip3p.itp`, `tip4p.itp`, etc.
- `residuetypes.dat`, `specbond.dat`: Residue type definitions and special bonds

#### `contactmap/`

Tools for generating contact maps between residues:

- `contact_map`: Executable for creating contact maps
- Source files (`.c`, `.h`) and compiled objects (`.o`) for the contact map utility
- `makefile`: Build instructions for the contact map utility

#### `martini300C-ff/`

Martini 3.0 coarse-grained force field files:

- `aminoacids.ff`: Coarse-grained parameters for amino acids
- `general.ff`: General force field parameters
- `modifications.ff`: Parameters for modified residues including crosslinks
- `small_molecule_martini3.ff`: Parameters for small molecules
- `citations.bib`: Bibliography of references for the Martini force field

#### `martini300C-mapping/`

Mapping files for converting atomistic to coarse-grained representations:

- Mapping files for standard amino acids (e.g., `ala.amber99.map`, `gly.amber99.map`)
- Mapping files for modified residues including crosslinks (e.g., `l4y.amber99.map`, `ly3.amber99.map`)
- Mapping files for terminal caps (`ace.amber99.map`, `nme.amber99.map`)
- `modifications.amber99.mapping`: Definition of mapping for all modifications

### Sequence Files

The `sequence/` directory contains files used in the sequence generation (homology modeling) stage:

#### `modeller/` - Custom MODELLER Libraries

These libraries enable MODELLER to handle collagen-specific residues and crosslinks:

- `par_mod.lib`: Special parameters for collagen modeling
- `restyp_mod.lib`: Custom residue type definitions including hydroxyproline and crosslinked residues
- `top_heav_mod.lib`: Topology library for heavy atoms
- `top_mod.lib`: General topology library for MODELLER

#### `fasta_sequences/` - Species-Specific Sequences

Pre-defined collagen sequences for various species: homo_sapiens, pan_troglodytes, pongo_abelii, callithrix_jacchus, otolemur_garnettii, mus_musculus, rattus_norvegicus, bos_taurus, canis_lupus, ailuropoda_melanoleuca, mustela_putorius, myotis_lucifugus, loxodonta_africana, danio_rerio, oreochromis_niloticus, oryzias_latipes, tetraodon_nigroviridis, xiphophorus_maculatus, pelodiscus_sinensis

#### Other Essential Files

- `crosslinks.csv`: Comprehensive database of lysine-derived crosslinks with the following information:
  - Species
  - Terminal type (N or C)
  - Residue combinations
  - Crosslink types (HLKNL, LKNL, PYD, DPD, etc.)
  - Atom positions for forming crosslinks
  - Translation/shift information used for crosslink optimization (e.g., D0-D5, D2-D1)

- `template.fasta`: Template sequence file for *Rattus norvegicus* (rat) collagen type I. Used as a reference in the homology modeling process.

- `template.pdb`: Template structure file for rat collagen type I. Provides the structural basis for homology modeling of other species.

## How These Files Are Used

### 1. Sequence Generation Stage

During sequence generation, ColBuilder uses:

- `template.fasta` and `template.pdb` as the structural and sequence reference
- The species-specific FASTA file from `fasta_sequences/` as the target sequence
- `crosslinks.csv` to determine which residues to crosslink and how
- MODELLER libraries in `modeller/` to handle special residues and crosslinks

This process creates a collagen triple helix with appropriate crosslinks for the selected species.
When adding non-enzymatic crosslinks via a mutated PDB workflow, the shift listed in
`crosslinks.csv` should be used for `crosslink_copies`.

### 2. Geometry Generation Stage

The geometry generation stage primarily uses the output from sequence generation rather than files in this directory.

### 3. Topology Generation Stage

During topology generation, ColBuilder uses:

- Force field files from `topology/amber99sb-star-ildnp.ff/` for standard all-atom simulations
- Modified parameters from `topology/amber99/` for handling crosslinked residues
- For coarse-grained simulations:
  - `martini300C-ff/` directory for Martini 3.0 force field parameters
  - `martini300C-mapping/` directory for mapping atomistic structures to coarse-grained models
  - `create_goVirt.py` for generating Go-like potentials in coarse-grained models
  - utilities from `selectors.py` for atom selection and manipulation

## Customizing Data Files

Advanced users may modify these files to extend ColBuilder's capabilities:

### Adding New Crosslink Types

1. Add new entries to `crosslinks.csv` with the new crosslink information
2. Update `restyp_mod.lib` to define the new residue types
3. Update `top_heav_mod.lib` and `top_mod.lib` with the topology
4. Update `amber99/aminoacids.rtp` with parameters for the new crosslinks
5. For coarse-grained simulations, also add:
   - Mapping definitions in `martini300C-mapping/` (create files like `newcrosslink.amber99.map`)
   - Force field parameters in `martini300C-ff/modifications.ff`

### Using Different Reference Structures

1. Replace `template.fasta` and `template.pdb` with alternative reference structures
2. Ensure the new reference has the same chain labeling conventions
3. Update the sequence alignment in the homology modeling process

### Adding New Species

1. Create a new FASTA file in `fasta_sequences/` directory following the naming convention (e.g., `newspecies.fasta`)
2. Add appropriate entries to `crosslinks.csv` for the new species

### Customizing Force Field Parameters

1. Modify files in `amber99/` to adjust force field parameters for all-atom simulations
2. For coarse-grained simulations, modify files in `martini300C-ff/`
3. Create a new force field directory following the GROMACS naming convention

**Warning**: Modifications require expertise in structural biology, force field parameterization, and understanding of the ColBuilder pipeline. Always back up original files before making changes.

## File Formats

- `.csv`: Comma-separated values containing crosslink definitions
- `.fasta`: Text-based format for protein sequences
- `.pdb`: Protein Data Bank format for 3D structural data
- `.rtp`, `.atp`, `.itp`: GROMACS parameter and topology files
- `.lib`: MODELLER library files for residue definitions and parameters

## Additional Resources

- [MODELLER Documentation](https://salilab.org/modeller/manual/)
- [GROMACS Force Field Documentation](https://manual.gromacs.org/current/reference-manual/topologies/force-field-organization.html)
- [PDB File Format Specification](https://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html)
- [Martini 3 Force Field Documentation](http://cgmartini.nl)
