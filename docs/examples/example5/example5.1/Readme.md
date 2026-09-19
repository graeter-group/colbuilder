<div align="center">
    <h1>Example 5.1</h1>
    <p>Enzymatic + non-enzymatic (AGE) workflow: sequence generation, additional crosslinking, geometry generation, and topology generation.</p>
</div>

Builds a *rattus_norvegicus* triple helix with terminal PYD (enzymatic)
crosslinks, adds a non-enzymatic Glucosepane crosslink on top via the
mutated-PDB workflow, then builds the microfibril and its amber99 topology.

## Usage

```bash
colbuilder --config_file config_step1_pyd_sequence.yaml --debug
colbuilder --config_file config_step2_add_glucosepane.yaml --debug
colbuilder --config_file config_step3_geometry_topology.yaml --debug
```

## Files

- `config_step1_pyd_sequence.yaml` — sequence generation, terminal PYD crosslinks.
- `config_step2_add_glucosepane.yaml` — mutated-PDB workflow, adds Glucosepane on top of step 1.
- `config_step3_geometry_topology.yaml` — geometry + amber99 topology generation from step 2's output.

## Output Files

- `rattusnorvegicus_alignment.fasta` — MSA from step 1.
- `rattusnorvegicus_N_PYD_C_PYD.pdb` — step 1's output; step 2's `mutated_pdb` input.
- `rattusnorvegicus_N_PYD_C_PYD+ADD1_Glucosepane.pdb` — step 2's output; step 3's `pdb_file` input.
- `collagen_fibril_rattus_norvegicus.pdb` — final microfibril structure.
- `rattus_norvegicus_topology_files/` — amber99 topology: `.top`, `.gro`, per-group
  `col_*.itp`/`posre_*.itp`, and the force-field directory.

## Notes

- Step 3 repeats `crosslink`/`n_term_type`/`c_term_type`/`additional_1_type` even though no
  sequence generation happens there: these declare what's actually in `pdb_file` (PYD +
  Glucosepane), so crosslink-type validation passes instead of flagging a mismatch.
- `col_*.itp` group pairings (which models are crosslink-bonded) can shift model IDs
  slightly between otherwise-identical runs — the crystal-building step's numbering for
  later-assigned models isn't fully deterministic; the pairing structure itself is stable.
