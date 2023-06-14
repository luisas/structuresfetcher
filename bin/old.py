#!/usr/bin/env python3
import Bio.PDB as bpdb
from Bio import PDB
import sys

structure = sys.argv[3]
start_res = int(sys.argv[1])
end_res = int(sys.argv[2])
outname = sys.argv[4]



try:
    class ResSelect(bpdb.Select):
        def accept_residue(self, res):
            # i checked, both inclusive!
            if res.id[1] >= start_res and res.id[1] <= end_res:
                return True
            else:
                return False

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

except:
    print("Error: ", structure)
    pass

