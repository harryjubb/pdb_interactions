# PDB Interactions

[CLI](https://en.wikipedia.org/wiki/Command-line_interface) for calculating all atom-atom, residue-residue, and chain-chain interactions in [PDB](http://www.rcsb.org/pdb/home/home.do) files.

## Dependencies

Written for Python 2.7.x. To install dependencies with [`pip`](https://pip.pypa.io/en/stable/installing/):

    pip install -r requirements.txt

## Usage

    python determine_interactions.py -h

## Residue Typing

Outputting the entity type (e.g. protein, DNA, RNA, etc) is optional, and depends on
data files generated from `https://github.com/harryjubb/pdb_residue_types`.
