import subprocess
import glob
import ntpath
from rdkit import Chem
from MolFileManipulation import Mol2ToMol


#CALCULATE WITH LIST
calcFilesTemp = ["TALWOS"]
calcFiles = []
for file in calcFilesTemp:
	calcFiles.append("G:\\!CSD-database\\"+file+".search1.mol2")
	#desktop - calcFiles.append("C:\\Users\\frederico\\Desktop\\"+file+".search1.mol2")    

#allMol2Files = glob.glob("G:\\!CSD-database\\*.mol2")
#calcFiles = allMol2Files[0:43783]
#del allMol2Files

calculating = open("calculating.csv", "w")
calculating.write("CSD;Info;Metal;Formula;ID;RMSD")

for mol2 in calcFiles:
#	if ntpath.basename(mol2) != "EGALUV.search1.mol2":
#		continue
	calculating.write("\n{:<9};".format(mol2.partition(".")[0]))
	try:
		obj1 = Mol2ToMol(mol2)
		#obj1.keepIdentifierFiles()
		obj1.runStereoisomerIdentifierRmsd(calculating)
		
	except Exception as e:
		if str(e) == "Chem.MolFromMol2File failed":
			calculating.write("E.rdkit")
		elif str(e) == "Metal number error - 0 metals":
			calculating.write("E.0Metals")
		else:
			calculating.write("E.{}".format(str(e)))

calculating.close()


#REMOVER:
#	def _writeMolFile(self):
#	justLetters = lambda enterString : ''.join([i for i in enterString if i.isalpha()])