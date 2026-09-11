"""
This module provides functionality to convert PDB files into FASTA format sequences.

It extracts amino acid sequences from the alpha-carbon (CA) atoms in a PDB file and maps
three-letter amino acid codes to one-letter codes. The resulting sequences are formatted
in FASTA format, with each chain represented as a separate entry.

Key Features:
--------------
1. **PDB to FASTA Conversion**:
   - Extracts amino acid sequences from PDB files based on CA atom records.
   - Supports standard amino acids, modified residues (e.g., MSE), and collagen-specific residues
     (e.g., HYP, L4Y, L5Y).

2. **Chain Handling**:
   - Processes multiple chains in a PDB file.
   - Outputs each chain as a separate FASTA entry.

3. **Custom Residue Mapping**:
   - Maps three-letter amino acid codes to one-letter codes.
   - Handles non-standard residues and maps unsupported residues to gaps ('-').

Usage:
------
This module can be used as a standalone script or imported into other Python programs.
When run as a script, it takes a PDB file as input and prints the corresponding FASTA
sequences to the console.

Example (Standalone Script):
----------------------------
```bash
python pdb2fasta.py input_structure.pdb
```
"""

import re
import sys


def pdb_to_fasta(pdb_file):
    aa3to1 = {
        "ALA": "A",
        "VAL": "V",
        "PHE": "F",
        "PRO": "P",
        "MET": "M",
        "ILE": "I",
        "LEU": "L",
        "ASP": "D",
        "GLU": "E",
        "LYS": "K",
        "ARG": "R",
        "SER": "S",
        "THR": "T",
        "TYR": "Y",
        "HIS": "H",
        "CYS": "C",
        "ASN": "N",
        "GLN": "Q",
        "TRP": "W",
        "GLY": "G",
        "MSE": "M",
        "HYP": "O",
        "L4Y": "-",
        "L5Y": "-",
        "L4X": "-",
        "L5X": "-",
        "LY4": "-",
        "LY5": "-",
        "LX4": "-",
        "LX5": "-",
        "LY2": "-",
        "LY3": "-",
        "LYX": "-",
        "LX2": "-",
        "LX3": "-",
        "LXX": "-",
        "L2Y": "-",
        "L3Y": "-",
        "LXY": "-",
        "L2X": "-",
        "L3X": "-",
        "LYY": "-",
    }

    ca_pattern = re.compile(
        "^ATOM\s{2,6}\d{1,5}\s{2}CA\s[\sA]([A-Z]{3})\s([\s\w])|"
        "^ATOM\s{2,6}\d{1,5}\s{2}CA\s[\sA](HYP|L4Y|L5Y|L4X|L5X|LY4|LY5|LX4|LX5)\s([\s\w])|"
        "^ATOM\s{2,6}\d{1,5}\s{2}CA\s[\sA](LY2|LY3|LYX|LX2|LX3|LXX|L2Y|L3Y|LXY|L2X|L3X|LYY)\s([\s\w])|"
        "^HETATM\s{0,4}\d{1,5}\s{2}CA\s[\sA](MSE)\s([\s\w])"
    )

    chain_dict = {}
    chain_list = []

    with open(pdb_file, "r") as fp:
        for line in fp:
            if line.startswith("ENDMDL"):
                break
            match_list = ca_pattern.findall(line)
            if match_list:
                resn = match_list[0][0] + match_list[0][2] + match_list[0][4]
                chain = match_list[0][1] + match_list[0][3] + match_list[0][5]
                if chain in chain_dict:
                    chain_dict[chain] += aa3to1[resn]
                else:
                    chain_dict[chain] = aa3to1[resn]
                    chain_list.append(chain)

    fasta_content = ""
    for chain in chain_list:
        fasta_content += f">{pdb_file}:{chain}\n{chain_dict[chain]}\n"

    return fasta_content


def main():
    if len(sys.argv) != 2:
        print("Usage: pdb2fasta <pdb_file>")
        sys.exit(1)

    pdb_file = sys.argv[1]
    fasta_content = pdb_to_fasta(pdb_file)
    print(fasta_content)


if __name__ == "__main__":
    main()
