<div align="center">
    <h1>Example 5.3B</h1>
    <p>Mixed-crosslink microfibril: HLKNL+Glucosepane, PYD+Pentosidine, and MOLD+Pentosidine combined by ratio.</p>
</div>

Extends Example 5.3A with a third variant: MOLD (divalent, LYS-LYS derived,
non-enzymatic) + Pentosidine. All three enzymatic/LYS-LYS-derived positions
(HLKNL, PYD, MOLD) share the same terminal definitions, and both AGEs
(Glucosepane, Pentosidine) share the same non-enzymatic positions, so the three
variants can be mixed together into one microfibril.

## Usage

```bash
colbuilder --config_file config_step1_hlknl_sequence.yaml --debug
colbuilder --config_file config_step2_add_glucosepane.yaml --debug
colbuilder --config_file config_step3_pyd_sequence.yaml --debug
colbuilder --config_file config_step4_add_pentosidine.yaml --debug
colbuilder --config_file config_step5_mold_sequence.yaml --debug
colbuilder --config_file config_step6_add_pentosidine.yaml --debug
colbuilder --config_file config_step7_mix_topology.yaml --debug
```

## Files

- `config_step1_hlknl_sequence.yaml` — sequence generation, terminal HLKNL crosslinks (variant A).
- `config_step2_add_glucosepane.yaml` — mutated-PDB workflow, adds Glucosepane to variant A.
- `config_step3_pyd_sequence.yaml` — sequence generation, terminal PYD crosslinks (variant B).
- `config_step4_add_pentosidine.yaml` — mutated-PDB workflow, adds Pentosidine to variant B.
- `config_step5_mold_sequence.yaml` — sequence generation, terminal MOLD crosslinks (variant C).
- `config_step6_add_pentosidine.yaml` — mutated-PDB workflow, adds Pentosidine to variant C.
- `config_step7_mix_topology.yaml` — mixing-only mode: combines variants A, B, C at a 60:20:20 ratio
  (`ratio_mix: "A:60 B:20 C:20"`), then builds amber99 topology.

## Output Files

- `rattusnorvegicus_alignment.fasta` — MSA; only the last sequence-generation step's
  alignment survives here, since all three write to the same filename.
- `rattusnorvegicus_N_HLKNL_C_HLKNL.pdb`, `rattusnorvegicus_N_HLKNL_C_HLKNL+ADD1_Glucosepane.pdb` — variant A.
- `rattusnorvegicus_N_PYD_C_PYD.pdb`, `rattusnorvegicus_N_PYD_C_PYD+ADD1_Pentosidine.pdb` — variant B.
- `rattusnorvegicus_N_MOLD_C_MOLD.pdb`, `rattusnorvegicus_N_MOLD_C_MOLD+ADD1_Pentosidine.pdb` — variant C.
- `collagen_fibril_rattus_norvegicus.pdb` — final mixed microfibril.
- `rattus_norvegicus_topology_files/` — amber99 topology: `.top`, `.gro`, per-group
  `col_*.itp`/`posre_*.itp`, and the force-field directory.

## Notes

- The achieved 60:20:20 split can land off by one model either way between runs — which
  specific models get selected for each type isn't seeded.
