import subprocess
import glob
from rdkit import Chem
from MolFileManipulation import Mol2ToMol

flagOnlyOne = False

if flagOnlyOne:
	calculating = open("calone.txt", "w")
	print("starting")
	try:
		obj1 = Mol2ToMol("ABACRE.search1.mol2")
		obj1.runStereoisomerIdentifierRmsd(calculating)

	except Exception as e:
		if str(e) == "Chem.MolFromMol2File failed":
			print(" rdkit coldn't read mol2 file\n")
		elif str(e) == "Metal number error - 0 metals":
			print(" Couldn't find any metal at .mol2 file\n")
		elif str(e) == "Metal number error - more than one":
			print(" More than one metal at .mol2 file\n")
		elif str(e) == "Number of ligands error":
			print(" Number of ligands need to be between 1 and 8\n")
		elif str(e) == "chelations not well defined":
			print(" Error on defining formula")
		else:
			print(" Error:  ",str(e))
	print("finished")
	exit()

else:

	allMol2Files = glob.glob("*.mol2")
	calculating = open("calculating.txt", "w")

	for mol2 in allMol2Files:
		calculating.write("\nCalculating:  {:>20}  :;".format(mol2))
		try:
			obj1 = Mol2ToMol(mol2)
			obj1.runStereoisomerIdentifierRmsd(calculating)
			
		except Exception as e:
			if str(e) == "Chem.MolFromMol2File failed":
				calculating.write(";rdkit coldn't read mol2 file")
			elif str(e) == "Metal number error - 0 metals":
				calculating.write(";Couldn't find any metal at .mol2 file")
			elif str(e) == "Metal number error - more than one":
				calculating.write(";More than one metal at .mol2 file")
			elif str(e) == "Number of ligands error":
				calculating.write(";Number of ligands need to be between 1 and 8")
			elif str(e) == "chelations not well defined":
				calculating.write(";Error on defining formula")
			else:
				calculating.write(";Error:  {}\n".format(str(e)))
	
	calculating.close()


#REMOVER:

#	def _writeMolFile(self):
#	justLetters = lambda enterString : ''.join([i for i in enterString if i.isalpha()])