from rdkit import Chem
from MolFileManipulation import Mol2ToMol
import subprocess

print("starting")

#subprocess.call("dir", shell=True)

fileMol2Name = "ACPNEU.mol2"

obj1 = Mol2ToMol(fileMol2Name)

obj1.writeCppInput()

#obj1.writeMolFile()

#mol = Chem.MolFromMolFile(obj1.fileMol2Name + '.mol')

#rank = list(Chem.CanonicalRankAtoms(mol, breakTies=True))

#encontrar quem e quem

#print('[',end='');print(*rank, sep=', ', end='');print(']')




print("finished")

#rank = list(Chem.CanonicalRankAtoms(mol, breakTies=False)) -simetria
