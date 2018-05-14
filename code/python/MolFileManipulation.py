from rdkit import Chem
from collections import Counter
import os

#transfer to a utility functions
justLetters = lambda enterString : ''.join([i for i in enterString if not i.isdigit()])

class Mol2ToMol:
	"""Class to handle mol2 and mol files"""

	#constant
	__flagMolecule = "@<TRIPOS>MOLECULE"
	__flagAtom = "@<TRIPOS>ATOM"
	__flagBond = "@<TRIPOS>BOND"
	__allMetals = ["Eu", "Fe", "V"]

	def __init__(self, fileMol2Name):
		#Defining class variables
		self.__fileStream = ""
		self.__molName = ""
		self.__listAtoms = []
		self.__listBonds = []
		self.__iMetal = 0
		self.__metalBondsLines = [] # indexed to __listBonds
		self.__ligandBondedToMetal = [] # start with 1 - subtract one to get corresponding atoms
		self.fileMol2Name = ""

		fileName, fileExtension = os.path.splitext(fileMol2Name)
		if not fileExtension == '.mol2':
			return
		self.fileMol2Name = fileName
	
		#Read mol2 file
		mol2Input = open(fileMol2Name,"r")
		self.__fileStream = mol2Input.read().splitlines()
		self._extractMol2Info()
		mol2Input.close()

	def printInfo(self):
		print("molName: ",self.__molName)
		print("\n\nAtoms:\n\n",*self.__listAtoms,sep="\n")
		print("\n\nBonds:\n\n",*self.__listBonds,sep="\n")
		print("metal bond lines: ",self.__metalBondsLines)
		print("ligand bond: ",self.__ligandBondedToMetal)
		for iMetalBond in self.__metalBondsLines:
			print("bond: ",self.__listBonds[iMetalBond])
		for iLigandBond in self.__ligandBondedToMetal:
			print("ligands: ",self.__listAtoms[iLigandBond-1])


	def writeCppInput(self):
		self.writeMolFile()
		mol = Chem.MolFromMolFile(self.fileMol2Name + '.mol')
		rank = list(Chem.CanonicalRankAtoms(mol, breakTies=False))

		xL = [] ; yL = [] ; zL = [] ; rankL = []
		for i in self.__ligandBondedToMetal:
			atomsColumns = self.__listAtoms[i-1].split()
			xL.append(atomsColumns[2])
			yL.append(atomsColumns[3])
			zL.append(atomsColumns[4])
			rankL.append(rank[i-1])
		zipped = list(zip(rankL, xL, yL, zL))
		zipped.sort()
		rankL, xL, yL, zL = zip(*zipped)
		
		
		
		#function que me de a composicao dessa coisa
		#menor quantidade menor letra.
		#atribuir novos valores para o dicionario
		#aplicar no rankL - substituindo os valores antigos pelos novos
		#
		
		print(Counter(rankL))
		dict = Counter(rankL)
		dictElements = []
		dictKeys = []
		for comp in dict:
			dictElements.append(dict[comp])
			dictKeys.append(comp)
		zipped2 = list(zip(dictElements, dictKeys))
		zipped2.sort()
		dictElements, dictKeys = zip(*zipped2)
		print("keys:  ",dictKeys)
		print("elements:  ", dictElements)

		cppInput = open(self.fileMol2Name + "-cpp.inp", "w")
		atomsColumns = self.__listAtoms[self.__iMetal - 1].split()
		cppInput.write("{:>10}{:>10}{:>10}{:>5}\n".format(
		atomsColumns[2],
		atomsColumns[3],
		atomsColumns[4],
		"-1"))
	
		for i in range(len(xL)):
			cppInput.write("{:>10}{:>10}{:>10}{:>5}\n".format(
			xL[i],
			yL[i],
			zL[i],
			rankL[i]))

		cppInput.write("end\n")
		cppInput.close()
	



	

	def writeMolFile(self):
		molFile = open(self.fileMol2Name + ".mol", "w")
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
			molFile.write(justLetters(listAtomsColumns[1]))
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

	def _extractMol2Info(self):
		i = 0
		# !!!!!!ADICIONAR PROTECAO CONTRA USUARIOS AQUI!!!!!!
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
		metalsInMol2File = []
		while i < len(self.__listAtoms):
			listAtomsColumns = self.__listAtoms[i].split()
			for atom in self.__allMetals:
				if listAtomsColumns[1].find(atom) > -1:
					metalsInMol2File.append(i)
			i+=1

		if len(metalsInMol2File) > 1:
			print("Too many metals")
			return

		self.__iMetal = metalsInMol2File[0] + 1
		i = 0
		while i < len(self.__listBonds):
			auxB1 = int(self.__listBonds[i].split()[1])
			auxB2 = int(self.__listBonds[i].split()[2])
			if auxB1 == self.__iMetal:
				self.__ligandBondedToMetal.append(auxB2)		
				self.__metalBondsLines.append(i)
			if auxB2 == self.__iMetal:
				self.__ligandBondedToMetal.append(auxB1)
				self.__metalBondsLines.append(i)
		
			i+=1

	


if __name__ ==  "__main__":
	print("Module to handle Mol files")
