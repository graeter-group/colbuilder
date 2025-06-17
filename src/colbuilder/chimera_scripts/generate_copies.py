import chimera
from chimera import runCommand
import os
import math


def apply_unit_cell_translation(model, translation, cell_params):
    """Apply unit cell translation using cell parameters."""
    a, b, c = cell_params["a"], cell_params["b"], cell_params["c"]
    alpha, beta, gamma = map(
        math.radians, [cell_params["alpha"], cell_params["beta"], cell_params["gamma"]]
    )

    ax = a
    ay = 0
    az = 0

    bx = b * math.cos(gamma)
    by = b * math.sin(gamma)
    bz = 0

    cx = c * math.cos(beta)
    cy = c * (math.cos(alpha) - math.cos(beta) * math.cos(gamma)) / math.sin(gamma)
    cz = math.sqrt(c * c - cx * cx - cy * cy)

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

# Required unit cell translations D0 & D5 (adjust commenting as needed)
TRANSLATIONS = [
    # (-4, 0, -4),  # D4 = For strand 1 (blue) - top (1-2 interactions)
    # (-3, 0, -3),  # D3 =For strand 2 (green)
    (-2, 0, -2),  # D2 =For strand 3 (pink)
    (-1, 0, -1),  # D1 = For strand 4 (purple)
    # (0, 0, 0),  # D0 = For strand 5 (orange) - reference position
    # (-5, 0, -4),  # D5 = For strand 1 (blue) - bottom (5-1 interactions)
]

input_pdb = os.environ.get("INPUT_PDB")

generated_pdbs = []
for i, trans in enumerate(TRANSLATIONS):
    model = chimera.openModels.open(input_pdb)[0]

    apply_unit_cell_translation(model, trans, CELL_PARAMS)

    filename = "copy_{}.pdb".format(i)
    runCommand("write format pdb #{} {}".format(model.id, filename))
    generated_pdbs.append(filename)
    print("  {}     |  Generated copy with translation {}".format(i, str(trans)))

    runCommand("close #{}".format(model.id))

with open("generated_pdbs.txt", "w") as f:
    for pdb in generated_pdbs:
        f.write(pdb + "\n")

print("\nGenerated {} PDB files.".format(len(generated_pdbs)))
