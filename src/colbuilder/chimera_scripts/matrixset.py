#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
import os
import sys
from chimera import runCommand as rc
from chimera import openModels, selection, Point

def print_usage():
    print("Usage: Set environment variables PDB_FILE, CRYSTALCONTACTS_FILE, SYSTEM_SIZE, and FIBRIL_LENGTH before running this script with Chimera")
    sys.exit(1)

pdb_file = os.environ.get('PDB_FILE')
crystalcontacts_file = os.environ.get('CRYSTALCONTACTS_FILE')
system_size = os.environ.get('SYSTEM_SIZE')
fibril_length = os.environ.get('FIBRIL_LENGTH')

if not all([pdb_file, crystalcontacts_file, system_size, fibril_length]):
    print_usage()

try:
    fibril_length = float(fibril_length) * 10 
    system_size = int(system_size)

    for i in range(system_size):
        rc("open " + pdb_file)
    
    start_pos = openModels.list()[0].atoms[0].coord()
    end_pos = openModels.list()[0].atoms[-1].coord()
    center_pos = (start_pos[2] + end_pos[2]) / 2.0
    start_pos[2] = center_pos - 0.5 * fibril_length
    end_pos[2] = center_pos + 0.5 * fibril_length

    rc("matrixset " + crystalcontacts_file + ".txt")

    for i, model in enumerate(openModels.list()):
        rc("write #{0} {1}.pdb".format(model.id, model.id))
        rc("del #{0}".format(model.id))
        rc("open {0}.pdb".format(i))
        os.remove("{0}.pdb".format(i))

    select_atoms = []
    for model in openModels.list():
        for residue in model.residues:
            mass = [atom.element.mass for atom in residue.atoms]
            com = Point([atom.coord() for atom in residue.atoms], mass)
            if start_pos[2] <= com[2] <= end_pos[2]:
                select_atoms.extend(residue.atoms)

    select_fibril = selection.ItemizedSelection()
    select_fibril.add(select_atoms)
    select_fibril.addImplied()
    selection.mergeCurrent(selection.REPLACE, select_fibril)
    rc("sel invert")
    rc("del sel")

    output_file = str(crystalcontacts_file) + "_id.txt"
    with open(output_file, 'w') as f:
        for model in openModels.list():
            rc("write #{0} {1}.pdb".format(model.id, model.id))
            f.write('Model {0}\n'.format(model.id))

    if os.path.exists(output_file):
        with open(output_file, 'r') as f:
            print(f.read())
    else:
        print("Error: Output file was not created: %s" % output_file)
        print("Directory contents: %s" % os.listdir(os.getcwd()))

except Exception as e:
    print("An error occurred: %s" % str(e))
    print("Current working directory: %s" % os.getcwd())
    print("Directory contents: %s" % os.listdir(os.getcwd()))
    sys.exit(1)