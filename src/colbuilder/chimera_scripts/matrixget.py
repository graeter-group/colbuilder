import os
import sys
from chimera import runCommand as rc

def print_usage():
    print("Usage: Set environment variables PDB_FILE, CONTACT_DISTANCE, and CRYSTALCONTACTS_FILE before running this script with Chimera")
    sys.exit(1)

pdb_file = os.environ.get('PDB_FILE')
d_contact = os.environ.get('CONTACT_DISTANCE')
crystalcontacts_file = os.environ.get('CRYSTALCONTACTS_FILE')

if not all([pdb_file, d_contact, crystalcontacts_file]):
    print_usage()

try:
    open_command = "open {}".format(pdb_file)
    rc(open_command)
    
    contact_command = "crystalcontacts #0 {} copies true schematic false".format(d_contact)
    rc(contact_command)
    
    matrix_command = "matrixget {}.txt".format(crystalcontacts_file)
    rc(matrix_command)
    
    if not os.path.exists("{}.txt".format(crystalcontacts_file)):
        print("Error: Crystal contacts file was not created: {}.txt".format(crystalcontacts_file))
except Exception as e:
    print("An error occurred: {}".format(str(e)))
    sys.exit(1)