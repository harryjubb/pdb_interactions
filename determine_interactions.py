import argparse
import logging
import os
import re

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import NeighborSearch

# MATCHES PDB ACCESSIONS
# E.G. 1XKK
PDB_REGEX = r'([0-9]{1}[a-zA-Z]{3})'
# MATCHES PDB-BIOMOLID
# E.G. 1XKK-1
PDB_BIOMOL_REGEX = PDB_REGEX + r'-([0-9]+)'
# MATCHES RSCB STYLE ASSEMBLIES
# E.G. 1XKK.pdb1
PDB_RCSB_ASM_REGEX = PDB_REGEX + r'\.pdb([0-9]+)'

PDB_REGEX = re.compile(PDB_REGEX)
PDB_BIOMOL_REGEX = re.compile(PDB_BIOMOL_REGEX)
PDB_RCSB_ASM_REGEX = re.compile(PDB_RCSB_ASM_REGEX)

if __name__ == '__main__':

    # ARGUMENT PARSING
    parser = argparse.ArgumentParser(description='''
# Determine Interactions

Calculate all atom-atom, residue-residue, and chain-chain interactions
in PDB files.
''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(
        'inputfile', type=str, help='Path to the PDB file to be analysed.')
    parser.add_argument('-i', '--interacting', type=float,
                        default=5.0, help='Distance cutoff for interactions.')
    parser.add_argument('-v', '--verbose',
                        action='store_true', help='Be chatty.')

    args = parser.parse_args()

    # SET ARGS TO CONSTANTS
    INPUT_FILE = args.inputfile
    INPUT_FILENAME = os.path.split(INPUT_FILE)[1]
    INTERACTION_THRESHOLD = args.interacting

    # LOGGING
    if args.verbose:
        logging.basicConfig(
            level=logging.INFO,
            format='%(levelname)s//%(asctime)s.%(msecs).03d//%(message)s',
            datefmt='%H:%M:%S')
    else:
        logging.basicConfig(
            level=logging.WARN,
            format='%(levelname)s//%(asctime)s.%(msecs).03d//%(message)s',
            datefmt='%H:%M:%S')

    logging.info('Program begin.')

    # DETECT PDB ACCESSION FROM FILENAME IF POSSIBLE
    PDB = ''
    BIOMOL = ''

    pdb_match = re.search(PDB_REGEX, INPUT_FILENAME)
    if pdb_match:
        PDB = pdb_match.group(1)

    pdb_biomol_match = re.search(PDB_BIOMOL_REGEX, INPUT_FILENAME)
    if pdb_biomol_match:
        PDB = pdb_biomol_match.group(1)
        BIOMOL_ID = pdb_biomol_match.group(2)

    pdb_rcsb_asm_match = re.search(PDB_RCSB_ASM_REGEX, INPUT_FILENAME)
    if pdb_rcsb_asm_match:
        PDB = pdb_rcsb_asm_match.group(1)
        BIOMOL_ID = pdb_rcsb_asm_match.group(2)

    # LOAD STRUCTURE
    structure = PDBParser().get_structure('structure', INPUT_FILE)
    structure_atoms = list(structure.get_atoms())

    logging.info('Loaded PDB structure (BioPython).')

    # CONSTRUCT KDTREE
    neighborsearch = NeighborSearch(structure_atoms)

    logging.info('Constructured NeighborSearch.')

    # GET INTERACTING ATOMS
    atom_pairs = neighborsearch.search_all

    # FINISH UP
    logging.info('Program end.')
