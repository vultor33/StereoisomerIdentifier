import subprocess
import glob
from rdkit import Chem
from MolFileManipulation import Mol2ToMol

print("starting")

try:
	obj1 = Mol2ToMol("CAHJAX.mol2")
	obj1.writeCppInput()
	obj1.runStereoisomerIdentifierRmsd()

except Exception as e:
	if str(e) == "Chem.MolFromMol2File failed":
		print(" rdkit coldn't read mol2 file\n")
	elif str(e) == "Metal number error":
		print(" only one metal is supported\n")
	else:
		print(str(e))

print("finished")
exit()


allMol2Files = glob.glob("*.mol2")

calculating = open("calculating.txt", "w")

for mol2 in allMol2Files:
	calculating.write("Calculating:  {:>20}".format(mol2))
	try:
		obj1 = Mol2ToMol(mol2)
		obj1.writeCppInput()
		obj1.runStereoisomerIdentifierRmsd()
			
	except Exception as e:
		if str(e) == "Chem.MolFromMol2File failed":
			calculating.write(" rdkit coldn't read mol2 file\n")
		elif str(e) == "Too many metals":
			calculating.write(" only one metal is supported\n")
		else:
			calculating.write( str(e) + "\n")

calculating.close()
