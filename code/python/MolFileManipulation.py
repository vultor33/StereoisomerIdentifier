import os
import subprocess
from rdkit import Chem
from collections import Counter
from FormulaHandling import FormulaHandling

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
		if mol is None:
			raise Exception('Chem.MolFromMol2File failed')
		self.__equivalenceRank = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
		

	def printInfo(self):
		print("molName: ",self.__molName)
		print("\n\nAtoms:\n\n",*self.__listAtoms,sep="\n")
		print("\n\nBonds:\n\n",*self.__listBonds,sep="\n")
		print("metal bond lines: ",self.__metalBondsLines)

	def runStereoisomerIdentifierRmsd(self, outputFile_):
		for iMetal in self.__metalsInMol2File:
			try:
				outputFile_.write(self.__listAtoms[iMetal].split()[1] + " ; ")
			
				status = self._writeCppInput(int(iMetal) + 1)
				if status == 1:
					outputFile_.write("L-1;0;")
				else:
					subprocess.call("StereoisomerIdentifierRmsd.exe " + self.__fileMol2Name + "-cpp.inp", shell=True)
					cppOutput_ = open(self.__fileMol2Name + "-cpp.inp.log","r")
					cppStream_ = cppOutput_.read().splitlines()
					outputFile_.write(cppStream_[1] + ";" + cppStream_[2] + ";")
			
			except Exception as e:
				if str(e) == "Chem.MolFromMol2File failed":
					outputFile_.write(" rdkit coldn't read mol2 file;;")
				elif str(e) == "Metal number error - 0 metals":
					outputFile_.write(" Couldn't find any metal at .mol2 file;;")
				elif str(e) == "Metal number error - more than one":
					outputFile_.write(" More than one metal at .mol2 file;;")
				elif str(e) == "Number of ligands error":
					outputFile_.write(" Number of ligands need to be between 1 and 8;;")
				elif str(e) == "chelations not well defined":
					outputFile_.write(" Error on defining formula;;")
				else:
					outputFile_.write(" Error:  {}\n".format(str(e)))

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
			raise Exception("Metal number error - 0 metals")


	def _writeCppInput(self, iMetal):
		#Finding ligands and chelations
		iChelates = []
		ligandsBondedToMetal  = []
		self._findChelations(iMetal, iChelates, ligandsBondedToMetal)

		if len(ligandsBondedToMetal) > 8:
			raise Exception("Number of ligands error")

		if len(ligandsBondedToMetal) == 1:
			return 1

		rankL = []
		for i in ligandsBondedToMetal:
			atomsColumns = self.__listAtoms[i-1].split()
			rankL.append(self.__equivalenceRank[i-1])
		objF = FormulaHandling()
		objF.generateMolecularFormula(rankL,iChelates)
		objFenum = FormulaHandling()
		objFenum.generateEnumerationFormula(rankL,iChelates)
		rankToTypesMap = objF.calculateNewMapBetweenAtomsAndTypes(rankL,iChelates)

		cppInput = open(self.__fileMol2Name + "-cpp.inp", "w")

		if objF.getFormula() == objFenum.getFormula():
			cppInput.write(objF.getFormula() + "\n")
			cppInput.write("chelates:  {}".format(len(iChelates)))
			for chel in iChelates:
				cppInput.write("  cI-length:  {}  cI:".format(len(chel)))
				for chelI in chel:
					cppInput.write(" {} ".format(chelI))
			cppInput.write("\n")
		else:
			cppInput.write(objFenum.getFormula() + "\n\n")

		atomsColumns = self.__listAtoms[self.__metalsInMol2File[0]].split()
		cppInput.write("{:>5}{:>10}{:>10}{:>10}{:>5}\n".format(
		atomsColumns[1],
		atomsColumns[2],
		atomsColumns[3],
		atomsColumns[4],
		"-1"))
		rankIndex = iter(rankL)
		for i in ligandsBondedToMetal:
			atomsColumns = self.__listAtoms[i-1].split()
			cppInput.write("{:>5}{:>10}{:>10}{:>10}{:>5}\n".format(
			atomsColumns[1],
			atomsColumns[2],
			atomsColumns[3],
			atomsColumns[4],
			rankToTypesMap[next(rankIndex)]))
		cppInput.write("end\n")
		cppInput.close()

		return 0
		


	def _findChelations(self, iMetal, iChelates, ligandsBondedToMetal):
		#building graph		
		graph = {}
		for i in range(len(self.__listAtoms)):
			graph[i+1] = self._getAtomBonds(i+1)
	
		auxLigandsBondedToMetal = graph[iMetal]
		for iLigand in auxLigandsBondedToMetal:
			ligandsBondedToMetal.append(iLigand)
		
		del graph[iMetal]
		
		chelates = []
		for i in range(len(ligandsBondedToMetal) - 1):
			alreadyFound = False
			for chel in chelates:
				if i in chel:
					alreadyFound = True
					break
			if alreadyFound:
				continue

			j = i + 1
			auxChelates = [i]
			while j < len(ligandsBondedToMetal):
				pathToJ = self._findAnyGraphPath(graph, ligandsBondedToMetal[i], ligandsBondedToMetal[j])
				if not pathToJ is None:
					auxChelates.append(j)
				j+=1
			chelates.append(auxChelates)		

		auxChelates = [s for s in chelates if len(s) != 1]
		for iChel in auxChelates:
			iChelates.append(iChel)
		
		
	
	def _findAnyGraphPath(self, graph, start, end, path=[]):
		path = path + [start]
		if start == end:
			return path
		if not start in graph:
			return None
		for node in graph[start]:
			if node not in path:
				newpath = self._findAnyGraphPath(graph, node, end, path)
				if newpath: return newpath
		return None	

	def _getAtomBonds(self, atomI):
		atomsBonded = []
		for i in range(len(self.__listBonds)):
			auxB1 = int(self.__listBonds[i].split()[1])
			auxB2 = int(self.__listBonds[i].split()[2])
			if auxB1 == atomI:
				atomsBonded.append(auxB2)		
			if auxB2 == atomI:
				atomsBonded.append(auxB1)
		return atomsBonded





	

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
	"Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Uub",
	"Al","Ga","Ge","In","Sn","Sb","Ti","Pb","Bi","Po"]

if __name__ ==  "__main__":
	print("Module to handle Mol files")
