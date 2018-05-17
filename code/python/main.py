import subprocess
import glob
from rdkit import Chem
from MolFileManipulation import Mol2ToMol

flagOnlyOne = False

if flagOnlyOne:

	print("starting")
	try:
		obj1 = Mol2ToMol("PEBDOP.mol2")
		obj1.runStereoisomerIdentifierRmsd()

	except Exception as e:
		if str(e) == "Chem.MolFromMol2File failed":
			print(" rdkit coldn't read mol2 file\n")
		elif str(e) == "Metal number error":
			print(" Couldn't find any metal at .mol2 file\n")
		elif str(e) == "Number of ligands error":
			print(" Number of ligands need to be between 4 and 9\n")
		else:
			print(" Error:  ",str(e))
	finally:
		print("\n")
	print("finished")
	exit()

else:

	allMol2Files = glob.glob("*.mol2")
	calculating = open("calculating.txt", "w")

	for mol2 in allMol2Files:
		calculating.write("Calculating:  {:>20}".format(mol2))
		try:
			obj1 = Mol2ToMol(mol2)
			obj1.runStereoisomerIdentifierRmsd()
				
		except Exception as e:
			if str(e) == "Chem.MolFromMol2File failed":
				calculating.write(" rdkit coldn't read mol2 file")
			elif str(e) == "Too many metals":
				calculating.write(" Couldn't find any metal at .mol2 file")
			elif str(e) == "Number of ligands error":
				calculating.write(" Number of ligands need to be between 4 and 9")
			else:
				calculating.write(" Error:  " + str(e))
		finally:
			calculating.write("\n")
	
	calculating.close()


#REMOVER:

#	def _writeMolFile(self):
#	justLetters = lambda enterString : ''.join([i for i in enterString if i.isalpha()])