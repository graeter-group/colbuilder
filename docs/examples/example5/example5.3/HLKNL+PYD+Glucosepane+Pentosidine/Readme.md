<div align="center">
    <h1>Example 5.3A</h1>
    <p>Mixed-crosslink microfibril: HLKNL+Glucosepane combined with PYD+Pentosidine by ratio.</p>
</div>

Builds two independent enzymatic + non-enzymatic *rattus_norvegicus* structures
(HLKNL+Glucosepane and PYD+Pentosidine), then mixes them into one microfibril
at an 80:20 ratio and generates its amber99 topology.

## Usage

```bash
colbuilder --config_file config_step1_hlknl_sequence.yaml --debug
colbuilder --config_file config_step2_add_glucosepane.yaml --debug
colbuilder --config_file config_step3_pyd_sequence.yaml --debug
colbuilder --config_file config_step4_add_pentosidine.yaml --debug
colbuilder --config_file config_step5_mix_topology.yaml --debug
```

## Files

- `config_step1_hlknl_sequence.yaml` — sequence generation, terminal HLKNL crosslinks (variant A).
- `config_step2_add_glucosepane.yaml` — mutated-PDB workflow, adds Glucosepane to variant A.
- `config_step3_pyd_sequence.yaml` — sequence generation, terminal PYD crosslinks (variant B).
- `config_step4_add_pentosidine.yaml` — mutated-PDB workflow, adds Pentosidine to variant B.
- `config_step5_mix_topology.yaml` — mixing-only mode: combines variant A and B at an 80:20 ratio
  (`ratio_mix: "A:80 B:20"`), then builds amber99 topology.

## Output Files

- `rattusnorvegicus_alignment.fasta` — MSA; only the last of the two sequence-generation
  steps' alignment survives here, since both write to the same filename.
- `rattusnorvegicus_N_HLKNL_C_HLKNL.pdb`, `rattusnorvegicus_N_HLKNL_C_HLKNL+ADD1_Glucosepane.pdb` — variant A.
- `rattusnorvegicus_N_PYD_C_PYD.pdb`, `rattusnorvegicus_N_PYD_C_PYD+ADD1_Pentosidine.pdb` — variant B.
- `collagen_fibril_rattus_norvegicus.pdb` — final mixed microfibril.
- `rattus_norvegicus_topology_files/` — amber99 topology: `.top`, `.gro`, per-group
  `col_*.itp`/`posre_*.itp`, and the force-field directory.

## Notes

- The two variants must use identical terminal-crosslink position definitions
  (`n_term_combination`/`c_term_combination`) to be mix-compatible; likewise their
  additional (non-enzymatic) crosslink positions must match each other.
- The achieved 80:20 split can land off by one model either way between runs — which
  specific models get selected for each type isn't seeded.
