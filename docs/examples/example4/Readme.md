<div align="center">
    <h1>Example 4</h1>
    <p>Topology-only workflow: Generate topology files from an existing mixed crosslinked fibril PDB</p>
</div>

## Overview

This example demonstrates **topology-only mode**, where ColBuilder generates GROMACS topology files from an existing fibril PDB file without running sequence generation or geometry generation steps. This is useful when you:

- Already have a complete fibril structure
- Want to regenerate topology with different force field parameters
- Need to create topology for a fibril obtained from another source
- Want to switch between atomistic (amber99) and coarse-grained (martini3) representations

## Workflow Summary

```
Input: Existing Mixed Fibril PDB
         ↓
   [Topology Generation Only]
         ↓
Output: Topology Files (.top, .itp, .gro)
```

## Files and Usage

### Input Files

- **`config_ex4.yaml`**  
  → Configuration file for topology-only mode  
  → Usage:  
  ```bash
  colbuilder --config_file config_ex4.yaml
  ```

- **`collagen_fibril_rattus_norvegicus_MIX.pdb`** (provided by user)  
  → Pre-existing fibril PDB file with mixed crosslinks  
  → Can be from a previous ColBuilder run or any compatible source

## Output Files

### Main Topology Files

- **`rattus_norvegicus_topology_files/`**  
  → Directory containing all topology files for GROMACS

- **`collagen_fibril_rattus_norvegicus.top`**  
  → Main GROMACS topology file

- **`collagen_fibril_rattus_norvegicus.gro`**  
  → GROMACS structure file with coordinates

### Individual Model Topology Files

- **`col_0.itp`, `col_1.itp`, ..., `col_N.itp`**  
  → Individual topology files for each collagen model in the fibril  
  → Contains atom definitions, bonds, angles, dihedrals  
  → Includes crosslink bonded terms (bonds, angles, dihedrals)

## Important Notes

### Automatic Connectivity Detection

In topology-only mode, ColBuilder:
- Automatically detects model connectivity from crosslink positions
- Identifies which models are connected through crosslinks
- Works with divalent (D), trivalent (T), or mixed (M) crosslink patterns
- Generates appropriate bonded terms for inter-molecular crosslinks

### Force Field Options

You can generate topology files with either:

- **`force_field: "amber99"`** - All-atom representation  
  → Best for detailed atomistic simulations  
  → Larger system size, higher computational cost

- **`force_field: "martini3"`** - Coarse-grained representation  
  → Suitable for longer timescales and larger systems  
  → Reduced system size, lower computational cost  
  → Requires Martinize2 installation

### Switching Between Force Fields

To generate topology with a different force field:

1. Keep the same input PDB file
2. Change only the `force_field` parameter
3. Run ColBuilder again

Example for atomistic topology:
```yaml
force_field: "amber99"
```

Example for coarse-grained topology:
```yaml
force_field: "martini3"
```

## Common Use Cases

### 1. Regenerate Topology After Manual Edits

If you manually edited a fibril PDB file and need new topology:
```bash
colbuilder --config_file config_ex4.yaml
```

### 2. Create Both Atomistic and Coarse-Grained Topologies

Generate amber99 topology:
```yaml
# config_amber.yaml
force_field: "amber99"
# ... rest of config
```

Then generate martini3 topology:
```yaml
# config_martini.yaml
force_field: "martini3"
# ... rest of config
```

```bash
colbuilder --config_file config_amber.yaml
colbuilder --config_file config_martini.yaml
```

### 3. Topology from External Source

If you obtained a collagen fibril PDB from another source:
1. Identify the crosslink types present
2. Make sure they are available in colbuilder database and set appropriate `n_term_type`, `c_term_type` and combinations
3. Run topology-only mode

## Debug Mode

For troubleshooting, enable debug mode:

```yaml
debug: true
topology_debug: true
```

This preserves:
- Intermediate topology files
- Temporary PDB files
- Detailed processing logs

## References

- For crosslink types and combinations: [crosslinks.csv](https://github.com/graeter-group/colbuilder/blob/main/src/colbuilder/data/sequence/crosslinks.csv)
- For configuration details: [Configuration Reference](https://github.com/graeter-group/colbuilder/tree/main/docs/configuration.md)
- For topology generation: [User Guide](https://github.com/graeter-group/colbuilder/tree/main/docs/user_guide.md)