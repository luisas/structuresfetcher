#!/usr/bin/env python3
# coding: utf-8

# In[41]:


import Bio.PDB as bpdb
from Bio import PDB
import sys

start_res = int(sys.argv[1])
end_res = int(sys.argv[2])
structure = sys.argv[3]
outname = sys.argv[4]
chain_id = sys.argv[5]
shift= int(sys.argv[6])

# Prepare real start and end residues
start_res = start_res - shift + 1
end_res = end_res - shift + 1


class ResSelect(bpdb.Select):
    def accept_residue(self, res):
        # i checked, both inclusive!
        if res.id[1] >= start_res and res.id[1] <= end_res and res.parent.id == chain_id:
            return True
        else:
            return False

# Select the desired residues and save to a PDB file (temp.pdb)
s = bpdb.PDBParser().get_structure('input', structure)
io = bpdb.PDBIO()
io.set_structure(s)
io.save("temp.pdb", ResSelect())



pdb_io = PDB.PDBIO()
pdb_parser = PDB.PDBParser()
final_structure = pdb_parser.get_structure("temp", "temp.pdb")

for model in final_structure:
    for chain in model:
        for i, residue in enumerate(chain.get_residues()):
            res_id = list(residue.id)
            res_id[1] = i+1
            residue.id = tuple(res_id)

pdb_io.set_structure(final_structure)
pdb_io.save(outname + ".pdb")



d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}


# iterate each model, chain, and residue
for model in final_structure:
   for chain in model:
       seq = []
       for residue in chain:
           seq.append(d3to1[residue.resname])
sequence = ''.join(seq)
seq_id = outname
f = open(outname+".fa", "w")
f.write(">"+seq_id + '\n')
f.write(sequence+ '\n')
f.close()
