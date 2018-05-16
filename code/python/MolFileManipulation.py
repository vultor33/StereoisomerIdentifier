import os
import subprocess
from rdkit import Chem
from collections import Counter

#transfer to a utility functions

untilPoint = lambda enterString : enterString.partition(".")[0]

justLetters = lambda enterString : ''.join([i for i in enterString if i.isalpha()])

class Mol2ToMol:
	"""Class to handle mol2 and mol files"""

	def __init__(self, fileMol2Name):
		#Defining class variables
		self.__fileStream = ""
		self.__molName = ""
		self.__listAtoms = []
		self.__listBonds = []
		self.__metalsInMol2File = []
		self.__fileMol2Name = ""
		self.__equivalenceRank = []

		self.__fileMol2Name, fileExtension = os.path.splitext(fileMol2Name)
		if not fileExtension == '.mol2':
			return

		#Read mol2 file
		mol2Input = open(fileMol2Name,"r")
		self.__fileStream = mol2Input.read().splitlines()
		self._extractMol2Info()
		mol2Input.close()

		#Canonicalizing
		mol = Chem.MolFromMol2File(self.__fileMol2Name + '.mol2', removeHs = False)
		if mol == None:
			raise Exception('Chem.MolFromMol2File failed')
		self.__equivalenceRank = list(Chem.CanonicalRankAtoms(mol, breakTies=False))

	def printInfo(self):
		print("molName: ",self.__molName)
		print("\n\nAtoms:\n\n",*self.__listAtoms,sep="\n")
		print("\n\nBonds:\n\n",*self.__listBonds,sep="\n")
		print("metal bond lines: ",self.__metalBondsLines)

	def runStereoisomerIdentifierRmsd(self):
		for iMetal in self.__metalsInMol2File:
			self._writeCppInput(int(iMetal) + 1)
			subprocess.call("StereoisomerIdentifierRmsd.exe " + self.__fileMol2Name + "-cpp.inp", shell=True)
	


	def _extractMol2Info(self):
		i = 0
		while i < len(self.__fileStream):
			if self.__fileStream[i].find(self.__flagMolecule) > -1:
				i+=1
				self.__molName = self.__fileStream[i]
				i+=1		
				auxNumbers = self.__fileStream[i].split()
				nAtoms = int(auxNumbers[0])
				nBonds = int(auxNumbers[1])

			if self.__fileStream[i].find(self.__flagAtom) > -1:
				self.__listAtoms = self.__fileStream[i+1:i+nAtoms+1]
				i += nAtoms + 1

			if self.__fileStream[i].find(self.__flagBond) > -1:
				self.__listBonds = self.__fileStream[i+1:i+nBonds+1]
				break
			i+=1

		i = 0
		while i < len(self.__listAtoms):
			listAtomsColumns = self.__listAtoms[i].split()
			if untilPoint(listAtomsColumns[5]) in self.__allMetals:
				self.__metalsInMol2File.append(i)
			i+=1
			
		if len(self.__metalsInMol2File) < 1:
			raise Exception("Metal number error")

	def _generateMolecularFormula(self, rankL, ligandsBondedToMetal):	
		for i in ligandsBondedToMetal:
			atomsColumns = self.__listAtoms[i-1].split()
			rankL.append(self.__equivalenceRank[i-1])
	
		rankLCounting = Counter(rankL)
		rankLCountingElems = []
		rankLCountingKeys = []
		for comp in rankLCounting:
			rankLCountingElems.append(rankLCounting[comp])
			rankLCountingKeys.append(comp)

		zipped = list(zip(rankLCountingElems, rankLCountingKeys))
		zipped.sort()
		rankLCountingElems, rankLCountingKeys = zip(*zipped)
		ligandsAmount = list(rankLCountingElems)
		ligandsAmount.reverse()
		ligandTypes = self.__alphabet[:len(ligandsAmount)]
		molecularFormula = "M"
		i = 0		
		while i < len(ligandTypes):
			if ligandsAmount[i] == 1:
				molecularFormula += ligandTypes[i]
			else:
				molecularFormula += ligandTypes[i] + str(ligandsAmount[i])
			i+=1
		newRankL = []
		for iRank in rankL:
			newRankL.append(rankLCountingKeys.index(iRank))
		for i in range(len(newRankL)):
			rankL[i] = newRankL[i]
		return molecularFormula


	def _writeCppInput(self, iMetal):
		ligandsBondedToMetal = []
		i = 0
		while i < len(self.__listBonds):
			auxB1 = int(self.__listBonds[i].split()[1])
			auxB2 = int(self.__listBonds[i].split()[2])
			if auxB1 == iMetal:
				ligandsBondedToMetal.append(auxB2)		
			if auxB2 == iMetal:
				ligandsBondedToMetal.append(auxB1)
			i+=1

		if len(ligandsBondedToMetal) > 9 or len(ligandsBondedToMetal) < 4:
			raise Exception("Number of ligands error")

		rankL = []
		formula = self._generateMolecularFormula(rankL, ligandsBondedToMetal)
		cppInput = open(self.__fileMol2Name + "-cpp.inp", "w")
		cppInput.write(formula + "\n")
		atomsColumns = self.__listAtoms[iMetal - 1].split()
		cppInput.write("{:>10}{:>10}{:>10}{:>5}\n".format(
		atomsColumns[2],
		atomsColumns[3],
		atomsColumns[4],
		"-1"))
		listRankL = iter(rankL)
		for i in ligandsBondedToMetal:
			atomsColumns = self.__listAtoms[i-1].split()
			cppInput.write("{:>10}{:>10}{:>10}{:>5}\n".format(
			atomsColumns[2],
			atomsColumns[3],
			atomsColumns[4],
			next(listRankL)))
		cppInput.write("end\n")
		cppInput.close()








	

	# DEACTIVATED - loading mol2 file
	def _writeMolFile(self):
		raise Exception("_writeMolFile deactivated")

		# This function need the following rules
		#self.__iMetal = metalsInMol2File[0] + 1
		#i = 0
		#while i < len(self.__listBonds):
		#	auxB1 = int(self.__listBonds[i].split()[1])
		#	auxB2 = int(self.__listBonds[i].split()[2])
		#	if auxB1 == self.__iMetal:
		#		self.__ligandBondedToMetal.append(auxB2)		
		#		self.__metalBondsLines.append(i)
		#	if auxB2 == self.__iMetal:
		#		self.__ligandBondedToMetal.append(auxB1)
		#		self.__metalBondsLines.append(i)
		#	i+=1
	
	
		molFile = open(self.__fileMol2Name + ".mol", "w")
		molFile.write(self.__molName)
		molFile.write("\n")
		molFile.write("StereoisomerIdentifier\n\n")
		molFile.write("{:>3}{:>3}".format(
		str(len(self.__listAtoms)),
		str(len(self.__listBonds) - len(self.__metalBondsLines))))
		molFile.write("  0     0  0  0  0  0  0999 V2000")
		molFile.write("\n")
		i = 0
		
		for atomLine in self.__listAtoms:
			listAtomsColumns = atomLine.split()
			molFile.write("{:>10}{:>10}{:>10}".format(
			listAtomsColumns[2],
			listAtomsColumns[3],
			listAtomsColumns[4]))
			molFile.write(" ")
			molFile.write(untilPoint(listAtomsColumns[5]))
			molFile.write("  0  0  0  0  0  0  0  0  0  0  0  0")
			molFile.write("\n")
		
		i = 0
		while i < len(self.__listBonds):
			if i in self.__metalBondsLines:
				i+=1
				continue
			listBondsColumns = self.__listBonds[i].split()
			molFile.write("{:>3}{:>3}".format(
			listBondsColumns[1],
			listBondsColumns[2]))
			molFile.write("  1  0  0  0")
			molFile.write("\n")
			i+=1
			
		molFile.write("M  END\n")
		molFile.close()

	



	#constant
	__flagMolecule = "@<TRIPOS>MOLECULE"
	__flagAtom = "@<TRIPOS>ATOM"
	__flagBond = "@<TRIPOS>BOND"
	__alphabet = ['a','b','c','d','e','f','g','h','i','j','k','l']
	__allMetals = ["Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
	"Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
	"La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
	"Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg",
	"Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr",
	"Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Uub"]

if __name__ ==  "__main__":
	print("Module to handle Mol files")
