<div align="center">
    <h1>Example 3A</h1>
    <p>
        Mixed-crosslink microfibril: HLKNL + PYD (enzymatic) combined with
        Glucosepane + Pentosidine (non-enzymatic) and mixed by ratio.
    </p>
</div>

This example builds two compatible triple-helix inputs and then mixes them
into a single microfibril with a defined ratio. Compatibility is enforced by
using identical position definitions within each crosslink class before mixing.

All commands below are run from this directory:
`docs/examples/non-enzymatic/example3-mix/HLKNL+PYD+Glucosepane+Pentosidine`

## Run Order

```bash
# 1) Enzymatic-only sequences
colbuilder --config_file config_ex5_step1_HLKNL_seq.yaml --debug
colbuilder --config_file config_ex5_step3_PYD_seq.yaml --debug

# 2) Add non-enzymatic crosslinks
colbuilder --config_file config_ex5_step2_HLKNL+Glucosepane.yaml --debug
colbuilder --config_file config_ex5_step4_PYD+Pentosidine.yaml --debug

# 3) Mix the two variants and generate topology
colbuilder --config_file config_ex5_step5.yaml --debug
```

## What Each Step Does

### 1) `config_ex5_step1_HLKNL_seq.yaml` and `config_ex5_step3_PYD_seq.yaml`
Standard sequence generation with enzymatic crosslinks at both N and C termini
(HLKNL vs PYD).  
Important: HLKNL and PYD must use the same position definitions at the termini,
otherwise the two inputs are not mix-compatible.

Key outputs:
- `rattusnorvegicus_N_HLKNL_C_HLKNL.pdb`
- `rattusnorvegicus_N_PYD_C_PYD.pdb`
- `rattusnorvegicus_alignment.fasta`  
  Alignment used by the modelling step.

### 2) `config_ex5_step2_HLKNL+Glucosepane.yaml` and `config_ex5_step4_PYD+Pentosidine.yaml`
Mutated-PDB workflow that adds non-enzymatic crosslinks (AGEs) on top of the
enzymatic inputs. In this example we add one (Glucosepane or Pentosidine), but
the workflow supports multiple non-enzymatic crosslinks.  
Important: Glucosepane and Pentosidine must share the same position definitions
(distinct from the enzymatic positions) so the two inputs can be mixed later.
Use the crosslink-specific shift for this second sequence step; it is listed in
`src/colbuilder/data/sequence/crosslinks.csv` alongside each crosslink entry.

Key outputs:
- `rattusnorvegicus_N_HLKNL_C_HLKNL+ADD1_Glucosepane.pdb`
- `rattusnorvegicus_N_PYD_C_PYD+ADD1_Pentosidine.pdb`

### 3) `config_ex5_step5.yaml`
Mixes the two final PDBs into a single microfibril with the requested ratio and
then generates GROMACS topologies.

Key settings in the YAML:
- `files_mix` lists the two input PDBs to mix.
- `ratio_mix: "A:80 B:20"` sets the composition (80% variant A, 20% variant B).

Key outputs:
- `collagen_fibril_rattus_norvegicus.pdb`  
  Final mixed microfibril structure.
- `rattus_norvegicus_topology_files/`  
  GROMACS topology folder containing:
  - `collagen_fibril_rattus_norvegicus.top` (topology entry point)
  - `collagen_fibril_rattus_norvegicus.gro` (coordinate file)
  - `col_*.itp` (per-chain/topology includes)
  - `posre_*.itp` (position restraints per chain)
  - `amber99sb-star-ildnp.ff/` (force-field directory)

## Notes
- HLKNL and PYD are enzymatic crosslinks at both termini and must use identical
  position definitions to be mix-compatible.
- Glucosepane and Pentosidine are non-enzymatic (ARG/LYS-derived) and must share
  the same position definitions (distinct from the enzymatic ones).
- You can add multiple non-enzymatic crosslinks; this example shows one.
- Mixing only happens in the final stage; all earlier steps only write PDBs.