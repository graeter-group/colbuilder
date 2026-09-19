<div align="center">
    <h1>Example 1</h1>
    <p>A complete workflow: sequence generation, geometry generation, and topology generation in a single run.</p>
</div>

Builds a divalent (HLKNL) crosslinked collagen microfibril for *homo_sapiens*
from scratch, and generates its amber99 topology, all from one config file.

## Usage

```bash
colbuilder --config_file config_full_pipeline.yaml
```

## Files

- `config_full_pipeline.yaml` — runs all three stages (sequence, geometry, topology) sequentially.

## Output Files

- `homosapiens_alignment.fasta` — MSA from sequence generation.
- `homosapiens_N_HLKNL_C_HLKNL.pdb` — crosslinked sequence PDB, input to geometry generation.
- `collagen_fibril_homo_sapiens.pdb` — final microfibril structure.
- `homo_sapiens_topology_files/` — amber99 topology: `collagen_fibril_homo_sapiens.top`,
  `collagen_fibril_homo_sapiens.gro`, per-group `col_*.itp`/`posre_*.itp`, and the force-field directory.

## Notes

- Intermediate files (`connect_from_colbuilder.txt`, `crystalcontacts_from_colbuilder*.txt`,
  `*_original.pdb`) are only written under `debug: true`, and even then land in
  `.tmp/geometry_gen/`, not the working directory.
