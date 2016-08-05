import os
import simplejson as json

# WHERE THIS CONFIG FILE LIVES
CONFIG_PATH = os.path.dirname(os.path.realpath(__file__))

# LOCATION OF PDB RESIDUE TYPE JSON FILES
# GENERATED USING `https://github.com/harryjubb/pdb_residue_types`
# USE THE SCRIPTS IN THAT REPO TO UPDATE
PDB_RESIDUE_TYPES_BY_RESIDUE_JSON_FILE = os.path.join(
    CONFIG_PATH, 'data', 'pdb_residue_types_by_residue.json')

PDB_RESIDUE_TYPES_BY_RESIDUE = None

if os.path.exists(PDB_RESIDUE_TYPES_BY_RESIDUE_JSON_FILE):
    with open(PDB_RESIDUE_TYPES_BY_RESIDUE_JSON_FILE, 'rb') as fo:
        PDB_RESIDUE_TYPES_BY_RESIDUE = json.load(fo)
