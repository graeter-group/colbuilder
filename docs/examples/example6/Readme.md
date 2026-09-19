<div align="center">
    <h1>Example 6</h1>
    <p>A complete workflow: sequence generation, geometry generation, and Martini3 (coarse-grained) topology generation in a single run.</p>
</div>

Builds a trivalent (PYD) crosslinked collagen microfibril for
*rattus_norvegicus* from scratch, and generates its Martini3 coarse-grained
topology, all from one config file — the same shape as Example 1, but
coarse-grained and with PYD instead of HLKNL.

## Usage

```bash
colbuilder --config_file config_full_pipeline_martini3.yaml
```

## Files

- `config_full_pipeline_martini3.yaml` — runs all three stages (sequence, geometry, topology) in one pass,
  with `force_field: "martini3"`.

## Output Files

- `rattusnorvegicus_alignment.fasta` — MSA from sequence generation.
- `rattusnorvegicus_N_PYD_C_PYD.pdb` — crosslinked sequence PDB, input to geometry generation.
- `collagen_fibril_rattus_norvegicus.pdb` — final microfibril structure.
- `rattus_norvegicus_martini3_topology_files/` — martini3 output (note the force-field-suffixed
  directory name — martini3 output always lands in `{species}_martini3_topology_files/`,
  unlike amber99's `{species}_topology_files/`): `.top`, a coarse-grained
  `collagen_fibril_CG_rattus_norvegicus.pdb` (not a `.gro` — convert with `gmx editconf`
  before `grompp` if you need one), `go-sites.itp`, `martini_v3.0.0*.itp`, and per-group
  `col_*.itp` files, including martinize2/vermouth's own per-connection intermediates
  (harmless — see Example 4's Readme).

## Notes

- Martini3 crosslink parametrization currently only covers PYD and HLKNL; other crosslink
  types (DPD, PYL, DPL, and the non-enzymatic AGE types) aren't yet parametrized for
  Martini3 and won't produce correct coarse-grained bonded terms — use `force_field:
  "amber99"` for those.
