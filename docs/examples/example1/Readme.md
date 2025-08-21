<div align="center">
    <h1>Example 1</h1>
    <p> A complete workflow demonstrating sequence generation, geometry generation, and topology generation
</p>
</div>

## Files and Usage

- `config_replace.yaml`  
  → Input config file for `colbuilder`  
  → Usage:  
  ```bash
  colbuilder --config_file config_ex1.yaml

## Output Files:

- `collagen_fibril_homo_sapiens.pdb` 
- `gen_fibril_homo_sapiens.gro`
  → Main structure output 

# further output files (from sequence/gemetry/crosslink settings):
- `connect_from_colbuilder.txt`
- `crystalcontacts_from_colbuilder_opt_id.txt`
- `crystalcontacts_from_colbuilder_opt.txt`
- `crystalcontacts_from_colbuilder.txt`
-  homosapiens_alignment.fasta`
- `homosapiens_N_HLKNL_C_HLKNL_original.pdb`
- `homosapiens_N_HLKNL_C_HLKNL.pdb`
- `homo_sapiens_topology_files/`
   → folder with the topology files

- `output_terminal.txt`  
  → Captured example terminal output
