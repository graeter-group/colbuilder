#!/usr/bin/env python
"""
Chimera script to generate unit cell copies for crosslink optimization.
This script reads the CROSSLINK_COPIES environment variable to determine
which pair of translations to apply.
"""

import chimera
from chimera import runCommand
import os
import math
import json

def apply_unit_cell_translation(model, translation, cell_params):
    """Apply unit cell translation using cell parameters."""
    a, b, c = cell_params["a"], cell_params["b"], cell_params["c"]
    alpha, beta, gamma = map(
        math.radians, [cell_params["alpha"], cell_params["beta"], cell_params["gamma"]]
    )
    
    # Calculate unit cell vectors
    ax = a
    ay = 0
    az = 0
    
    bx = b * math.cos(gamma)
    by = b * math.sin(gamma)
    bz = 0
    
    cx = c * math.cos(beta)
    cy = c * (math.cos(alpha) - math.cos(beta) * math.cos(gamma)) / math.sin(gamma)
    cz = math.sqrt(c * c - cx * cx - cy * cy)
    
    # Apply translation
    na, nb, nc = translation
    for a in model.atoms:
        curr_coord = a.coord()
        new_x = curr_coord.x + na * ax + nb * bx + nc * cx
        new_y = curr_coord.y + na * ay + nb * by + nc * cy
        new_z = curr_coord.z + na * az + nb * bz + nc * cz
        a.setCoord(chimera.Point(new_x, new_y, new_z))

# Cell parameters from CRYST1 record for collagen I
CELL_PARAMS = {
    "a": 39.970,
    "b": 26.950,
    "c": 677.900,
    "alpha": 89.24,
    "beta": 94.59,
    "gamma": 105.58,
}

# Unit cell translations mapping
TRANSLATION_MAP = {
    "D4": (-4, 0, -4),  # For strand 1 (blue) - top (1-2 interactions)
    "D3": (-3, 0, -3),  # For strand 2 (green)
    "D2": (-2, 0, -2),  # For strand 3 (pink)
    "D1": (-1, 0, -1),  # For strand 4 (purple)
    "D0": (0, 0, 0),    # For strand 5 (orange) - reference position
    "D5": (-5, 0, -4),  # For strand 1 (blue) - bottom (5-1 interactions)
}

# Get configuration from environment
input_pdb = os.environ.get("INPUT_PDB")
crosslink_copies_str = os.environ.get("CROSSLINK_COPIES", '["D0", "D5"]')

try:
    # Parse the crosslink copies configuration
    crosslink_copies = json.loads(crosslink_copies_str)
    print("Requested copies: {}".format(crosslink_copies))
    
    # Validate we have exactly 2 copies
    if len(crosslink_copies) != 2:
        raise ValueError("Expected exactly 2 crosslink copies, got {}".format(len(crosslink_copies)))
    
    # Validate the copy names
    for copy_name in crosslink_copies:
        if copy_name not in TRANSLATION_MAP:
            raise ValueError("Invalid translation name: {}. Valid options: {}".format(
                copy_name, ", ".join(TRANSLATION_MAP.keys())))
    
except Exception as e:
    print("Error parsing CROSSLINK_COPIES: {}".format(e))
    print("Using default D0 and D5")
    crosslink_copies = ["D0", "D5"]

# Generate the two copies
print("\nGenerating 2 PDB copies:")
print("Copy | Translation | Vector")
print("-----|-------------|-------")

generated_pdbs = []
for i, copy_name in enumerate(crosslink_copies):
    translation = TRANSLATION_MAP[copy_name]
    
    # Open model
    model = chimera.openModels.open(input_pdb)[0]
    
    # Apply translation
    apply_unit_cell_translation(model, translation, CELL_PARAMS)
    
    # Save copy
    filename = "copy_{}.pdb".format(i)
    runCommand("write format pdb #{} {}".format(model.id, filename))
    generated_pdbs.append(filename)
    
    print("  {}  |      {}      | {}".format(i, copy_name, translation))
    
    # Close model
    runCommand("close #{}".format(model.id))

# Write the list of generated PDBs
with open("generated_pdbs.txt", "w") as f:
    for pdb in generated_pdbs:
        f.write(pdb + "\n")

print("\nGenerated 2 PDB files for crosslink optimization.")