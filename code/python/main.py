from rdkit import Chem
from MolFileManipulation import Mol2ToMol
import subprocess

print("starting")

fileMol2Name = "ACPNEU.mol2"

obj1 = Mol2ToMol(fileMol2Name)

obj1.writeCppInput()

#subprocess.call("StereoisomerIdentifierRmsd.exe " + obj1.fileMol2Name + "-cpp.inp", shell=True)

#print("finished")

#rank = list(Chem.CanonicalRankAtoms(mol, breakTies=False)) -simetria
