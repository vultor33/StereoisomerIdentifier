import subprocess
import glob
from rdkit import Chem
from MolFileManipulation import Mol2ToMol

print("starting")

allMol2Files = glob.glob("*.mol2")

calculating = open("calculating.txt", "w")
for mol2 in allMol2Files:
	calculating.write("Calculating:  {:>20}".format(mol2))
	
	try:
		obj1 = Mol2ToMol(mol2)
		obj1.writeCppInput()
		if obj1.sucess:
			subprocess.call("StereoisomerIdentifierRmsd.exe " + obj1.fileMol2Name + "-cpp.inp", shell=True)
			calculating.write(" sucess\n")
		else:
			calculating.write(" failed\n")
			
	except:
		calculating.write(" failed\n")

calculating.close()
print("finished")
