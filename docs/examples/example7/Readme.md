<div align="center">
    <h1>Example 7</h1>
    <p>Deterministic, named-residue crosslink replacement using <code>manual_replacements</code>.</p>
</div>

Removes exactly one specific PYD crosslink trio from a *rattus_norvegicus*
fibril by naming its residues directly, rather than a random ratio-based
selection (contrast with Example 3's `ratio_replace`) — useful for
reproducible, targeted edits, e.g. testing one crosslink's specific
mechanical contribution.

## Usage

`manual_replacements` targets residues inside per-model `{id}.caps.pdb`
files, which only exist once geometry has been built and are only kept on
disk under `debug: true` (in `.tmp/geometry_gen/`). Finding real targets is a
two-step process:

```bash
# 1) Build geometry with debug: true, so per-model caps files are kept
colbuilder --config_file config_step1_geometry_debug.yaml

# 2) Inspect a caps file for real crosslink residues, e.g.:
grep -E "LYX|LY2|LY3" .tmp/geometry_gen/T/58.caps.pdb | cut -c18-26 | sort -u
#   LY2 B   5
#   LY3 C   9
#   LYX C 103
# These three residues are one physical PYD crosslink (a trivalent trio —
# always replace all three together, never one alone).

# 3) Rerun with replace_bool + manual_replacements naming those residues
colbuilder --config_file config_step2_manual_replace.yaml
```

## Files

- `config_step1_geometry_debug.yaml` — plain geometry generation, `debug: true`, used only to
  discover real target residues as shown above.
- `config_step2_manual_replace.yaml` — same geometry config, plus `replace_bool: true` and
  `manual_replacements` naming all three residues of the one PYD trio to remove.

Both configs keep `debug: true`: step 1 needs it to keep the caps files at all, and step 2
keeps it for consistency and to inspect `.tmp/replace_crosslinks/` afterward.

## Output Files

- `rattusnorvegicus_N_PYD_C_PYD.pdb` — shared input PDB for both steps.
- `collagen_fibril_rattus_norvegicus.pdb` — final microfibril, with exactly that one PYD
  crosslink removed; every other crosslink is untouched.

Confirmed directly in the mutated caps file (`.tmp/replace_crosslinks/T/58.caps.pdb`): all
three targeted residues (chain B/5, chain C/9, chain C/103) are `LYS` — PYD's parent
residue — and no `LY2`/`LY3`/`LYX` markers remain in that model.

## Notes

- Which model gets which numeric ID depends on the crystal-building step and can shift
  slightly between otherwise-identical runs. To make this reproducible across runs, always rediscover real targets via step 1
  rather than assuming a fixed ID from this Readme applies to your own run.
- As with Example 3, `manual_replacements` (like `ratio_replace` and `replace_file`) is not
  silently overridden by colbuilder's auto-fix-unpaired mechanism (see Example 5.2's Readme).
