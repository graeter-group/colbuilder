<div align="center">
    <h1>Example 1</h1>
    <p>
        Enzymatic + non-enzymatic (Glucosepane) workflow: sequence generation,
        additional crosslinking, geometry generation, and topology generation.
    </p>
</div>

This example runs ColBuilder in three sequential stages to (1) place enzymatic
crosslinks in a collagen triple helix, (2) add a non-enzymatic crosslink
(Glucosepane) on top of the enzymatic structure, and (3) build the microfibril
geometry plus GROMACS topology.

## Run Order

```bash
# 1) Enzymatic crosslinks on a fresh sequence
colbuilder --config_file config_ex5_step1_enzymatic.yaml --debug

# 2) Add non-enzymatic crosslink(s) on the pre-mutated PDB
colbuilder --config_file config_ex5_step2_nonenzymatic.yaml --debug

# 3) Build geometry and topology using the final PDB
colbuilder --config_file config_ex5_step3.yaml --debug
```

## What Each Step Does

### 1) `config_ex5_step1_enzymatic.yaml`
Standard sequence generation. It builds a triple helix from the species sequence
and applies terminal enzymatic crosslinks (N- and C-terminal (i.e. PYD)).

Key outputs:
- `rattusnorvegicus_alignment.fasta`  
  Multiple sequence alignment used by the modeller step.
- `rattusnorvegicus_N_PYD_C_PYD.pdb`  
  Triple-helix model with enzymatic crosslinks applied and optimized.

### 2) `config_ex5_step2_nonenzymatic.yaml`
Mutated-PDB workflow. It takes the enzymatic PDB from step 1 and applies an
additional non-enzymatic crosslink (i.e. Glucosepane) at the specified residues.
Use the crosslink-specific shift for this second sequence step; it is listed in
`src/colbuilder/data/sequence/crosslinks.csv` alongside each crosslink entry.

Key outputs:
- `rattusnorvegicus_N_PYD_C_PYD+ADD1_Glucosepane.pdb`  
  The same triple helix, now carrying the additional non-enzymatic crosslink.

### 3) `config_ex5_step3.yaml`
Geometry generation builds the microfibril from the final PDB and then generates
the force-field topology files for simulation.

Key outputs:
- `collagen_fibril_rattus_norvegicus.pdb`  
  Final microfibril geometry 
- `rattus_norvegicus_topology_files/`  
  GROMACS topology folder containing:
  - `collagen_fibril_rattus_norvegicus.top` (topology entry point)
  - `collagen_fibril_rattus_norvegicus.gro` (coordinate file)
  - `col_*.itp` (per-chain/topology includes)
  - `posre_*.itp` (position restraints per chain)
  - `amber99sb-star-ildnp.ff/` (force-field directory)

## Notes
- The non-enzymatic step expects the enzymatic PDB from step 1 to exist in this
  directory.
- The geometry step uses the PDB from step 2 as its input.