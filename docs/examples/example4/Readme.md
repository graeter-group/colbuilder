<div align="center">
    <h1>Example 4</h1>
    <p>Topology-only workflow: generate topology files from an existing mixed crosslinked fibril PDB.</p>
</div>

Generates GROMACS topology for a pre-built, mixed-crosslink *rattus_norvegicus*
fibril, with no sequence or geometry generation. Run twice here to produce
both an atomistic (amber99) and a coarse-grained (martini3) topology from the
same input, for comparison.

## Usage

```bash
colbuilder --config_file config_topology_only_amber99.yaml    # -> rattus_norvegicus_topology_files/
colbuilder --config_file config_topology_only_martini3.yaml   # -> rattus_norvegicus_martini3_topology_files/
```

## Files

- `collagen_fibril_rattus_norvegicus_MIX.pdb` — pre-existing fibril with mixed
  (divalent + trivalent) crosslinks, provided as input.
- `config_topology_only_amber99.yaml` — topology-only mode, amber99.
- `config_topology_only_martini3.yaml` — topology-only mode, martini3.

Both configs set `crosslink: false`: `pdb_file` is a genuinely mixed structure, so
crosslink-type validation (which expects one declared type) isn't meaningful here.

## Output Files

- `rattus_norvegicus_topology_files/` — amber99 output: `collagen_fibril_rattus_norvegicus.top`,
  `.gro`, per-group `col_*.itp`/`posre_*.itp`, and the force-field directory.
- `rattus_norvegicus_martini3_topology_files/` — martini3 output: `.top`, a coarse-grained
  `collagen_fibril_CG_rattus_norvegicus.pdb` (not a `.gro` — convert with `gmx editconf`
  before `grompp` if you need one), `go-sites.itp`, `martini_v3.0.0*.itp`, and per-group
  `col_*.itp` files. Also includes martinize2/vermouth's own per-connection intermediate
  files (`col_<id>.<connect_id>*.itp`) — colbuilder copies these into the output directory
  alongside the final merged ones; they're harmless and not something to clean up by hand.

## Notes

- Switching force field only requires changing `force_field` (and rerunning), everything
  else about topology-only mode stays the same.
