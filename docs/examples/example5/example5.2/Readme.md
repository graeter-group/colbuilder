<div align="center">
    <h1>Example 2</h1>
    <p>
        Enzymatic + non-enzymatic (Glucosepane) workflow with crosslink replacement
        using Chimera swapaa before topology generation.
    </p>
</div>

This example mirrors Example 5.1, but adds a replacement stage that converts a
user-defined fraction of crosslink markers into standard amino acids using the
`swapaa` Chimera script. This is useful when you want to reduce crosslink
density while keeping the rest of the model intact.

## Run Order

```bash
# 1) Enzymatic crosslinks on a fresh sequence
colbuilder --config_file config_ex5_step1_enzymatic.yaml --debug

# 2) Add non-enzymatic crosslink(s) on the pre-mutated PDB
colbuilder --config_file config_ex5_step2_nonenzymatic.yaml --debug

# 3) Build geometry, replace selected crosslinks, and generate topology
colbuilder --config_file config_ex5_step3.yaml --debug
```

## What Each Step Does

### 1) `config_ex5_step1_enzymatic.yaml`
Standard sequence generation. It builds a triple helix from the species sequence
and applies terminal enzymatic crosslinks (N- and C-terminal (i.e PYD)).

Key outputs:
- `rattusnorvegicus_alignment.fasta`  
  Multiple sequence alignment used by the modeller step.
- `rattusnorvegicus_N_PYD_C_PYD.pdb`  
  Triple-helix model with enzymatic crosslinks applied and optimized.

### 2) `config_ex5_step2_nonenzymatic.yaml`
Mutated-PDB workflow. It takes the enzymatic PDB from step 1 and applies an
additional non-enzymatic crosslink (Glucosepane) at the specified residues.
Use the crosslink-specific shift for this second sequence step; it is listed in
`src/colbuilder/data/sequence/crosslinks.csv` alongside each crosslink entry.

Key outputs:
- `rattusnorvegicus_N_PYD_C_PYD+ADD1_Glucosepane.pdb`  
  The same triple helix, now carrying the additional non-enzymatic crosslink.

### 3) `config_ex5_step3.yaml`
Builds the microfibril geometry, then replaces a fraction of crosslinks using
Chimera `swapaa`.

Key replacement settings in the YAML:
- `replace_bool: true` enables the replacement stage.
- `ratio_replace: 30` replaces 30% of eligible markers.
- `ratio_replace_scope: "all"` selects from enzymatic + non-enzymatic markers.

Key outputs:
- `collagen_fibril_rattus_norvegicus.pdb`  
  Final microfibril geometry after replacement.
- `manual_replacements.txt`  
  The actual replacement instructions applied by Chimera.

## Notes
- The non-enzymatic step expects the enzymatic PDB from step 1 to exist in this
  directory.
- The geometry step uses the PDB from step 2 as its input.
- Replacement uses Chimera `swapaa` to mutate the selected crosslink markers to
  standard residues while preserving the rest of the structure.