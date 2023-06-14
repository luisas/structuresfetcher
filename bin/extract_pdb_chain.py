#!/usr/bin/env python3

from Bio.PDB import PDBParser, Select,PDBIO
import sys

# Specify the chain ID to extract
pdb_file = sys.argv[1]
chain_id = sys.argv[2]
outfile = sys.argv[3]

# Custom class to select the desired chain
class ChainSelector(Select):
    def accept_chain(self, chain):
        return chain.get_id() == chain_id

# Create a PDB parser object
parser = PDBParser()
structure = parser.get_structure("input", pdb_file)
io = PDBIO()
io.set_structure(structure)
io.save(outfile, ChainSelector())