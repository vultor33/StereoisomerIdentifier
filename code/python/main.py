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
	calculating = open("calculating.csv", "w")
	calculating.write("CSD;Info;Metal;Formula;ID;RMSD")

	for mol2 in allMol2Files:
		if mol2 != "EGOCOT.search1.mol2":
			continue
		calculating.write("\n{:<9};".format(mol2.partition(".")[0]))
		try:
			obj1 = Mol2ToMol(mol2)
			obj1.runStereoisomerIdentifierRmsd(calculating)
			
		except Exception as e:
			if str(e) == "Chem.MolFromMol2File failed":
				calculating.write("E.rdkit")
			elif str(e) == "Metal number error - 0 metals":
				calculating.write("E.0Metals")
			else:
				calculating.write("E.{}".format(str(e)))
	
	calculating.close()



#READING AND ANALYZING RESULTS
#import csv
#with open('calculating.csv', 'r') as csvfile:
#	spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')
#	for row in spamreader:
#		if len(row) > 3:
#			print(row[4].split('-'))
#exit()



#REMOVER:

#	def _writeMolFile(self):
#	justLetters = lambda enterString : ''.join([i for i in enterString if i.isalpha()])