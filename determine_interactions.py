import argparse
import logging
import os
import re

import simplejson as json

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import NeighborSearch

LEVEL_MAP = {
    'A': 'atom',
    'R': 'residue',
    'C': 'chain'
}

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

This program assumes that the PDB file has already been cleaned and only takes
into account the first model.
''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(
        'inputfile', type=str, help='Path to the PDB file to be analysed.')
    parser.add_argument('-i', '--interacting', type=float,
                        default=5.0, help='Distance cutoff for interactions.')
    parser.add_argument('-o', '--outputs', type=str,
                        default='cr', help='Granularity of output: string'
                        'including letters "c" for chain level, "r" for '
                        'residue level, "a" for atom level.')
    parser.add_argument('-v', '--verbose',
                        action='store_true', help='Be chatty.')

    args = parser.parse_args()

    # SET ARGS TO CONSTANTS
    INPUT_FILE = args.inputfile
    INPUT_FILE_SPLITEXT = os.path.splitext(INPUT_FILE)[0]
    INPUT_FILENAME = os.path.split(INPUT_FILE)[1]
    INTERACTION_THRESHOLD = args.interacting
    OUTPUTS = args.outputs.upper()

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
    BIOMOL_ID = ''

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

    # GET INTERACTIONS
    logging.info('Calculating interactions...')
    for interaction_level in 'ARC':

        if interaction_level in OUTPUTS:

            logging.info('Calculating interactions for {}s...'.format(LEVEL_MAP[interaction_level]))

            pairs = neighborsearch.search_all(INTERACTION_THRESHOLD,
                                              level=interaction_level)

            logging.info('Search complete for {}s.'.format(LEVEL_MAP[interaction_level]))

            logging.info('Organising interactions for {}s...'.format(LEVEL_MAP[interaction_level]))
            interactions = {}

            for entities in pairs:

                entity1, entity2 = entities

                id1 = entity1.get_full_id()
                id2 = entity2.get_full_id()

                if interaction_level == 'A':
                    res1 = entity1.get_parent()
                    res2 = entity2.get_parent()

                    res1 = res1.resname
                    res2 = res2.resname

                    entity1 = ','.join([id1[2], str(id1[3][1]) + id1[3][2].strip() + '`' +  res1, entity1.name])
                    entity2 = ','.join([id2[2], str(id2[3][1]) + id2[3][2].strip() + '`' +  res2, entity2.name])

                elif interaction_level == 'R':
                    entity1 = ','.join([id1[2], str(id1[3][1]) + id1[3][2].strip() + '`' +  entity1.resname])
                    entity2 = ','.join([id2[2], str(id2[3][1]) + id2[3][2].strip() + '`' +  entity2.resname])


                elif interaction_level == 'C':
                    entity1 = entity1.id
                    entity2 = entity2.id

                if entity1 not in interactions:
                    interactions[entity1] = []

                if entity2 not in interactions:
                    interactions[entity2] = []

                if entity2 not in interactions[entity1]:
                    interactions[entity1].append(entity2)

                if entity1 not in interactions[entity2]:
                    interactions[entity2].append(entity1)

            for entity in interactions:
                interactions[entity] = sorted(interactions)

            logging.info('Organisation complete for {}s.'.format(LEVEL_MAP[interaction_level]))

            logging.info('Writing JSON for {}s...'.format(LEVEL_MAP[interaction_level]))

            pairs_json = {
                'input': INPUT_FILE,
                'pdb': PDB,
                'biomol_id': BIOMOL_ID,
                'level': LEVEL_MAP[interaction_level],
                'residue_interactions': interactions
            }

            with open('.'.join([INPUT_FILE_SPLITEXT, LEVEL_MAP[interaction_level], 'interactions', 'json']), 'wb') as fo:
                json.dump(pairs_json, fo)

            logging.info('JSON output written for {}s.'.format(LEVEL_MAP[interaction_level]))

    # FINISH UP
    logging.info('Program end.')
