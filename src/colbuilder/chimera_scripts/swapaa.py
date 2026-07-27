#!/usr/bin/env python
"""
Chimera script to replace specified residues with lysine.

This script is called by the GeometryReplacer to swap amino acids in a PDB file.
Python 2.7 compatible version for use with Chimera.
"""

import sys
import os
from chimera import runCommand as rc

def main():
    """
    Process replacement instructions and execute Chimera commands.
    
    Expected arguments:
    1. Base filename of the replacement instructions file (without extension)
    2. System type directory
    """
    if len(sys.argv) < 3:
        print("Error: Missing required arguments.")
        print("Usage: chimera --nogui --script swapaa.py replacement_file system_type")
        sys.exit(1)
    
    # Get command line arguments
    file = str(sys.argv[1])
    system_type = str(sys.argv[2])
    
    # Add extension if not present
    if not file.endswith('.txt'):
        file = file + '.txt'
    
    # Check if file exists
    if not os.path.exists(file):
        print("Error: Replacement file not found: {0}".format(file))
        sys.exit(1)
    
    # Process each line of the replacement file
    with open(file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            # Parse replacement instruction
            replace = [m for m in line.split(' ') if m]
            
            if len(replace) < 4:
                print("Warning: Invalid replacement instruction: {0}".format(line))
                continue
                
            pdb_file = replace[0]
            residue_type = replace[1].upper() if len(replace) > 1 else "LYS"
            if residue_type not in {"LYS", "ARG"}:
                residue_type = "LYS"
            residue_id = replace[2]
            # Do NOT force lower-case: PDB chain IDs are upper-case (A/B/C). The old
            # .lower() built specs like '#0:770.b' which fail to match chain 'B' if
            # Chimera treats chain IDs case-sensitively, silently skipping the swap.
            chain_id = replace[3].replace('\n', '').strip()

            # Construct path to PDB file
            pdb_path = os.path.join(system_type, pdb_file)

            # Open PDB file
            print("Opening {0}".format(pdb_path))
            rc("open {0}".format(pdb_path))

            # Fail-proof chain matching: resolve chain_id against the chain IDs that
            # actually exist in the opened model, case-insensitively, and use the
            # real case. Works whether Chimera matches chains case-sensitively or
            # not; falls back to chain_id as given if the model can't be inspected.
            try:
                from chimera import openModels, Molecule
                mols = openModels.list(modelTypes=[Molecule])
                if mols:
                    present = [r.id.chainId for r in mols[0].residues if r.id.chainId]
                    exact = [c for c in present if c == chain_id]
                    ci = [c for c in present if c.lower() == chain_id.lower()]
                    if exact:
                        chain_id = exact[0]
                    elif ci:
                        chain_id = ci[0]
            except Exception as exc:
                print("Chain-id introspection failed ({0}); using '{1}' as-is".format(exc, chain_id))

            # Replace residue with requested type (default LYS)
            print("Replacing residue {0}.{1} with {2}".format(residue_id, chain_id, residue_type))
            rc("swapaa {0} #0:{1}.{2}".format(residue_type, residue_id, chain_id))
            
            # Write modified structure back to file
            print("Writing updated structure to {0}".format(pdb_path))
            rc("write #0 {0}".format(pdb_path))
            
            # Delete model to prepare for next iteration
            rc("del #0")
    
    print("Replacement operations completed successfully")

if __name__ == "__main__":
    main()