<div align="center">
    <h1>Example 2</h1>
    <p>Mixing two crosslink types (divalent HLKNL and trivalent PYD) into a single microfibril by ratio.</p>
</div>

Builds two separate crosslinked sequences (HLKNL and PYD) for *homo_sapiens*,
then combines them into one mixed-crosslink microfibril at an 80:20 ratio
using mixing-only mode (no geometry/topology generation on the mix itself).

## Usage

```bash
colbuilder --config_file config_step1_hlknl_sequence.yaml   # 1) divalent (HLKNL) sequence
colbuilder --config_file config_step2_pyd_sequence.yaml     # 2) trivalent (PYD) sequence
colbuilder --config_file config_step3_mix_ratio.yaml        # 3) mix the two by ratio
```

## Files

- `config_step1_hlknl_sequence.yaml` — sequence generation, HLKNL (divalent) crosslinks.
- `config_step2_pyd_sequence.yaml` — sequence generation, PYD (trivalent) crosslinks.
- `config_step3_mix_ratio.yaml` — mixing-only mode, combines the two sequence PDBs at an 80:20 ratio.

## Output Files

- `homosapiens_alignment.fasta` — MSA output; only the last of the two sequence-generation
  steps' alignment survives here, since both write to the same filename.
- `homosapiens_N_HLKNL_C_HLKNL.pdb`, `homosapiens_N_PYD_C_PYD.pdb` — the two crosslinked
  sequence PDBs, consumed by `config_step3_mix_ratio.yaml`'s `files_mix`.
- `collagen_fibril_homo_sapiens.pdb` — final mixed-crosslink microfibril.

## Notes

- The achieved type distribution (target 80:20) can land off by one model either way
  between runs; which specific models get selected for each type isn't seeded.
