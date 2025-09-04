<div align="center">
    <h1>Example 3</h1>
    <p>  Modifying an existing ColBuilder-generated fibril to decrease crosslink density
</p>
</div>

## Files and Usage

- `config_replace.yaml`  
  → Input config file for `colbuilder`  
  → Usage:  
  ```bash
  colbuilder -f config_file config_replace.yaml
## Output Files

- `collagen_fibril_homo_sapiens.pdb`  
  → Main structure output

- `homosapiens_N_HLKNL_C_HLKNL_original.pdb`  
  → From sequence/geometry/crosslink settings

- `homosapiens_N_HLKNL_C_HLKNL.pdb`  
  → Modified/processed version

- `NC/`  
  → Directory containing additional outputs (e.g., chains or fragments)

- `output_terminal.txt`  
  → Captured example terminal output

