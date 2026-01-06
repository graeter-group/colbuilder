<div align="center">
    <h1>Example 3B</h1>
    <p>
        Mixed-crosslink microfibril: HLKNL + PYD + MOLD combined with
        Glucosepane + Pentosidine and mixed by ratio.
    </p>
</div>

This example extends the mixed-crosslink workflow by adding MOLD. MOLD is
non-enzymatic, but it is LYS-LYS derived like enzymatic crosslinks, so it must
share the same terminal positions as HLKNL/PYD for mixing. Glucosepane and
Pentosidine are LYS-ARG derived and must share positions with each other (at
non-enzymatic sites), not with the enzymatic positions.

## Run Order

```bash
# 1) Enzymatic-position sequences (HLKNL, PYD, MOLD)
colbuilder --config_file HLKNL_seq.yaml --debug
colbuilder --config_file PYD_seq.yaml --debug
colbuilder --config_file MOLD_seq.yaml --debug

# 2) Add non-enzymatic crosslinks
colbuilder --config_file HLKNL+Glucosepane.yaml --debug
colbuilder --config_file PYD+Pentosidine.yaml --debug
colbuilder --config_file MOLD+Pentosidine.yaml --debug

# 3) Mix the three inputs and generate topology
colbuilder --config_file mix+topology.yaml --debug
```

## What Each Step Does

### 1) `HLKNL_seq.yaml`, `PYD_seq.yaml`, `MOLD_seq.yaml`
Standard sequence generation with terminal crosslinks at both N and C termini.
HLKNL and PYD are enzymatic; MOLD is LYS-LYS derived and must use the same
terminal positions for mix-compatibility.

Key outputs:
- `rattusnorvegicus_N_HLKNL_C_HLKNL.pdb`
- `rattusnorvegicus_N_PYD_C_PYD.pdb`
- `rattusnorvegicus_N_MOLD_C_MOLD.pdb`
- `rattusnorvegicus_alignment.fasta`

### 2) `HLKNL+Glucosepane.yaml`, `PYD+Pentosidine.yaml`, `MOLD+Pentosidine.yaml`
Mutated-PDB workflow that adds non-enzymatic crosslinks (AGEs) on top of each
enzymatic-position input. In this example, one AGE is added per input, but the
workflow supports multiple non-enzymatic crosslinks.
Use the crosslink-specific shift for this second sequence step; it is listed in
`src/colbuilder/data/sequence/crosslinks.csv` alongside each crosslink entry.

Key outputs:
- `rattusnorvegicus_N_HLKNL_C_HLKNL+ADD1_Glucosepane.pdb`
- `rattusnorvegicus_N_PYD_C_PYD+ADD1_Pentosidine.pdb`
- `rattusnorvegicus_N_MOLD_C_MOLD+ADD1_Pentosidine.pdb`

### 3) `mix+topology.yaml`
Mixes the three PDBs into a single microfibril using the requested ratio and
then generates GROMACS topologies.

Key settings in the YAML:
- `files_mix` lists the three input PDBs to mix.
- `ratio_mix: "A:60 B:20 C:20"` sets the composition.

Key outputs:
- `collagen_fibril_rattus_norvegicus.pdb`
- `rattus_norvegicus_topology_files/`
  - `collagen_fibril_rattus_norvegicus.top`
  - `collagen_fibril_rattus_norvegicus.gro`
  - `col_*.itp`
  - `posre_*.itp`
  - `amber99sb-star-ildnp.ff/`

## Notes
- HLKNL, PYD, and MOLD must share identical terminal position definitions.
- Glucosepane and Pentosidine must share the same non-enzymatic positions,
  distinct from the enzymatic positions.
