import os
import ntpath
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
		self.__graphFindTries = 0


		self.__fileMol2Name, fileExtension = os.path.splitext(fileMol2Name)
		self.__fileMol2Name = ntpath.basename(self.__fileMol2Name)
		if not fileExtension == '.mol2':
			return

		#Read mol2 file
		mol2Input = open(fileMol2Name,"r")
		self.__fileStream = mol2Input.read().splitlines()
		self._extractMol2Info()
		mol2Input.close()

		self._writeMolFile()

		#Canonicalizing
		#mol = Chem.MolFromMol2File(fileMol2Name, removeHs = False)
		mol = Chem.MolFromMolFile(self.__fileMol2Name + '.mol', removeHs = False)
		if os.path.isfile(self.__fileMol2Name + '.mol'):
			os.remove(self.__fileMol2Name + '.mol')
		if mol is None:
			raise Exception('Chem.MolFromMol2File failed')
		self.__equivalenceRank = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
		#Chem.AssignStereochemistry(mol, flagPossibleStereoCenters=True)
		#self.__equivalenceRank = [int(a.GetProp('_CIPRank')) for a in mol.GetAtoms()]
		

	def printInfo(self):
		print("molName: ",self.__molName)
		print("\n\nAtoms:\n\n",*self.__listAtoms,sep="\n")
		print("\n\nBonds:\n\n",*self.__listBonds,sep="\n")
		print("metal bond lines: ",self.__metalBondsLines)

	def runStereoisomerIdentifierRmsd(self, outputFile_):
		if len(self.__metalsInMol2File) == 2:
			metal1 = untilPoint(self.__listAtoms[self.__metalsInMol2File[0]].split()[5])
			metal2 = untilPoint(self.__listAtoms[self.__metalsInMol2File[1]].split()[5])
			metalConnected = False
			try:
				metalConnected = self._checkIfIsConnected(self.__metalsInMol2File[0],self.__metalsInMol2File[1])
			except Exception as e:
				metalConnected = False
			
			if metal1 == metal2 and metalConnected:
				outputFile_.write("Dimetal;")
			else:
				outputFile_.write("Mix;")
		else:
			outputFile_.write("Mix;")
			
		for iMetal in self.__metalsInMol2File:
			try:
				outputFile_.write(self.__listAtoms[iMetal].split()[1] + ";")
		
				info = self._writeCppInput(int(float(iMetal)) + 1)
				if info[0] == 1:
					outputFile_.write(info[1] + ";" + info[2] + "-L-1;0;")
				else:
					subprocess.call("StereoisomerIdentifierRmsd.exe " + self.__fileMol2Name + "-cpp.inp", shell=True)
					cppOutput_ = open(self.__fileMol2Name + "-cpp.inp.log","r")
					cppStream_ = cppOutput_.read().splitlines()
					cppOutput_.close()
					if cppStream_[1] == "failed":
						outputFile_.write(info[1] + ";" + info[2] + "-E.Polyedron;" + cppStream_[2] + ";")
					elif cppStream_[1] == "rmsdfailed":
						outputFile_.write(info[1] + ";E.RMSD;" + cppStream_[2] + ";")
					else:
						outputFile_.write(info[1] + ";" + info[2] + "-" + cppStream_[1] + ";" + cppStream_[2] + ";")
		
		
			except Exception as e:
				if str(e) == "Metal number error - 0 metals":
					outputFile_.write("E.NoMetal;;;")
				elif str(e) == "No ligands":
					outputFile_.write("E.NoLigands;;;")
				elif str(e) == "Metal number error - more than one":
					outputFile_.write("E.ManyMetals;;;")
				elif str(e) == "Number of ligands error":
					outputFile_.write("E.NL>8;;;")
				elif str(e) == "chelations not well defined":
					outputFile_.write("E.Formula;;;")
				elif str(e) == "path error":
					outputFile_.write("E.Graph;;;")
				else:
					outputFile_.write("E.{};;;".format(str(e)))
		
		fileLog = self.__fileMol2Name + "-cpp.inp.log"
		fileInp = self.__fileMol2Name + "-cpp.inp"
		if os.path.isfile(fileLog):
			os.remove(fileLog)
		if os.path.isfile(fileInp):
			os.remove(fileInp)



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

		returnInfo = []
		if len(ligandsBondedToMetal) == 0:
			raise Exception("No ligands")
		if len(ligandsBondedToMetal) == 1:
			returnInfo.append(1)
			returnInfo.append('a')
			returnInfo.append('a')
			return returnInfo

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

		atomsColumns = self.__listAtoms[iMetal-1].split()
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

		returnInfo.append(0)
		returnInfo.append(objF.getFormula())
		returnInfo.append(objFenum.getFormula())
		return returnInfo

	def _checkIfIsConnected(self, atom1, atom2):		
		graph = {}
		for i in range(len(self.__listAtoms)):
			graph[i+1] = self._getAtomBonds(i+1)
		
		self.__graphFindTries = 0
		pathToJ = self._findAnyGraphPath(graph, int(float(atom1)) + 1, int(float(atom2)) + 1)
		return not pathToJ is None


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
				self.__graphFindTries = 0
				pathToJ = self._findAnyGraphPath(graph, ligandsBondedToMetal[i], ligandsBondedToMetal[j])
				if not pathToJ is None:
					auxChelates.append(j)
				j+=1
			chelates.append(auxChelates)		

		auxChelates = [s for s in chelates if len(s) != 1]
		for iChel in auxChelates:
			iChelates.append(iChel)
		
		
	
	def _findAnyGraphPath(self, graph, start, end, path=[]):
		self.__graphFindTries +=1
		if self.__graphFindTries > 100000:
			raise Exception("path error")

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
		linesWithMetalBonds = []
		for iMetal in self.__metalsInMol2File:
			for i in range(len(self.__listBonds)):
				auxB1 = int(self.__listBonds[i].split()[1])
				auxB2 = int(self.__listBonds[i].split()[2])
				if auxB1 == iMetal+1:
					if not i in linesWithMetalBonds:
						linesWithMetalBonds.append(i)
				if auxB2 == iMetal+1:
					if not i in linesWithMetalBonds:
						linesWithMetalBonds.append(i)

		molFile = open(self.__fileMol2Name + ".mol", "w")
		molFile.write(self.__molName)
		molFile.write("\n")
		molFile.write("StereoisomerIdentifier\n\n")
		molFile.write("{:>3}{:>3}".format(
		str(len(self.__listAtoms)),
		str(len(self.__listBonds) - len(linesWithMetalBonds)))) # remove all metal bonds
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
			if i in linesWithMetalBonds:
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
	"Al","Ga","Ge","In","Sn","Sb","Tl","Pb","Bi","Po"]

if __name__ ==  "__main__":
	print("Module to handle Mol files")
