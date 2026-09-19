<div align="center">
    <h1>Example 3</h1>
    <p>Modifying an existing ColBuilder-generated fibril to decrease crosslink density.</p>
</div>

Removes 30% of the enzymatic (HLKNL) crosslinks from an existing *homo_sapiens*
fibril, chosen at random, via direct replacement mode (no geometry generation).

## Usage

```bash
colbuilder --config_file config_ratio_replace.yaml
```

## Files

- `config_ratio_replace.yaml` — direct replacement mode: `replace_bool: true`, `ratio_replace: 30`.

**Important:** `replace_file` points at `collagen_fibril_homo_sapiens.pdb` — the same name
colbuilder writes its output to. Running this example for real **overwrites the input
fibril in place**. Back it up first if you want to keep the original.

## Output Files

- `collagen_fibril_homo_sapiens.pdb` — the input fibril, overwritten in place by the
  replaced structure when run for real. Kept here as the pre-replacement input.
- `collagen_fibril_homo_sapiens_after_replace.pdb` — what that file becomes after running
  `config_ratio_replace.yaml` (kept separately here so the example stays reproducible).

