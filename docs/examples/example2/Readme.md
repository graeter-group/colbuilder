<div align="center">
    <h1>Example 2</h1>
    <p>  Creating two sequences with different crosslink types and mixing them to form a heterogeneous fibril
</p>
</div>

## Files and Usage

- `example2.sh`  
  → `uses colbuilder 3x`  
  → Usage:  
  ```bash
  ./example2.sh
- `config_human-D.yaml` and `config_human-T.yaml`
  → config files for colbuilder that produce a triplehelix with the respective crosslinks.
- `config_mix.yaml`
  → config file for mixing the two triplehelices accordning to the specification in the file.
## Output Files

- `collagen_fibril_homo_sapiens.pdb`  
  → Main structure output

- `homosapiens_N_PYD_C_PYD.pdb`  
  → triplehelix with PYD crosslink

- `homosapiens_N_HLKNL_C_HLKNL.pdb`  
  → triplehelix with HLKNL crosslink

- `connect_from_colbuilder.txt`
- `crystalcontacts_from_colbuilder_opt_id.txt`
- `crystalcontacts_from_colbuilder_opt.txt`
- `crystalcontacts_from_colbuilder.txt`
  → Additional outputs (gemoetry files)

- `output_terminal.txt`  
  → Captured example terminal output

