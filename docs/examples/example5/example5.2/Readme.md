<div align="center">
    <h1>Example 5.2</h1>
    <p>Enzymatic + non-enzymatic workflow, then reducing crosslink density via ratio-based replacement.</p>
</div>

Builds the same PYD + Glucosepane *rattus_norvegicus* structure as Example
5.1, but replaces 30% of its crosslink markers with standard residues during
geometry generation, using Chimera `swapaa`.

## Usage

```bash
colbuilder --config_file config_step1_pyd_sequence.yaml --debug
colbuilder --config_file config_step2_add_glucosepane.yaml --debug
colbuilder --config_file config_step3_geometry_replace.yaml --debug
```

## Files

- `config_step1_pyd_sequence.yaml` — sequence generation, terminal PYD crosslinks.
- `config_step2_add_glucosepane.yaml` — mutated-PDB workflow, adds Glucosepane on top of step 1.
- `config_step3_geometry_replace.yaml` — geometry generation, then `replace_bool: true` with
  `ratio_replace: 30` and `ratio_replace_scope: "all"` to replace 30% of crosslinks
  (drawn from both PYD and Glucosepane markers) with standard residues.

No topology generation runs in this example (`topology_generator` isn't set).

## Output Files

- `rattusnorvegicus_alignment.fasta` — MSA from step 1.
- `rattusnorvegicus_N_PYD_C_PYD.pdb` — step 1's output; step 2's `mutated_pdb` input.
- `rattusnorvegicus_N_PYD_C_PYD+ADD1_Glucosepane.pdb` — step 2's output; step 3's `pdb_file` input.
- `collagen_fibril_rattus_norvegicus.pdb` — final microfibril, after replacement.
- `manual_replacements.txt` — unpaired crosslink markers colbuilder detected automatically
  before replacement (each would leave a dangling half-crosslink if left alone).
- `manual_replacements_applied.txt` — the full list of residues actually converted to
  standard amino acids by the ratio-based replacement, one `<caps_file> <resname> <resid>
  <chain>` entry per residue.

## Notes

- Because `ratio_replace` is set explicitly, colbuilder's auto-fix-unpaired step defers to
  it instead of silently forcing its own replacement list. `manual_replacements.txt` is
  informational (what it found), not what actually got replaced.
- Step 3 repeats `crosslink`/`n_term_type`/`c_term_type`/`additional_1_type` so that
  crosslink-type validation matches what's actually in `pdb_file`.
