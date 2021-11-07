import os
import io
import ntpath
from polycip.ErrorMessages import ErrorMessages
from polycip.Utilities import Utilities


# Um objeto para lidar com o formato mol
# Um objeto para lidar com o executavel cpp
# Um struct com list atoms, bond, equivalence rank e etc.

class MolFileHandling:
	"""Class to handle mol files"""

	def __init__(self, fileName):
		#objects
		self.__errorMessages_ = ErrorMessages()
		self.__util_ = Utilities()

		#constants
		self.__flagMolecule = "@<TRIPOS>MOLECULE"
		self.	__flagAtom = "@<TRIPOS>ATOM"
		self.__flagBond = "@<TRIPOS>BOND"
		self.__temporaryFilePreffix = "temporaryMol-"

		#variables
		self.__molFile = None
		self.__molFormat = ""
		self.__fileMol2Name = ""
		self.__fileStream = ""
		self.__molName = ""
		self.__listAtoms = []
		self.__listBonds = []
		self.__listOfBondsInPairs = []
		self.__metalsInMol2File = []

		#initialization
		self._initializeObject(fileName)


	def __del__(self):
		if self.__molFile != None:
			if not self.__molFile.closed:
				self.__molFile.close()
			if os.path.isfile(self.getTemporaryMolName() + self.__molFormat):
				os.remove(self.getTemporaryMolName() + self.__molFormat)
				
	def getTemporaryMolName(self):
		return self.__temporaryFilePreffix + self.getBaseFileName()

	def getBaseFileName(self):
		return ntpath.basename(self.__fileMol2Name)

	def getMetalsInMol2FileList(self):
		return self.__metalsInMol2File

	def getListAtoms(self):
		return self.__listAtoms

	def getListBonds(self):
		return self.__listBonds
	
	def getBondsInPairs(self):
		return self.__listOfBondsInPairs

	def writeMolFile(self):
		self.__molFormat = ".mol"
		self._openMolFile()
		self._writeMolFileFormat()
		self.__molFile.close()

	def writemolstream(self):
		self.__molFile = io.StringIO()
		self._writeMolFileFormat()		
		return self.__molFile.getvalue()

	def writeMol2File(self):
		self.__molFormat = ".mol2"
		self._openMolFile()
		self._writeMol2FileFormat()
		self.__molFile.close()


# PRIVATE
	def _initializeObject(self,fileName):
		self._checkExtension(fileName)
		self._readMol2File()

	def _writeMolFileFormat(self):
		linesWithMetalBonds = self._getListOfLinesWithMetalBonds()
		self.__molFile.write(self.__molName)
		self.__molFile.write("\n")
		self.__molFile.write("StereoisomerIdentifier\n\n")
		self.__molFile.write("{:>3}{:>3}".format(
		str(len(self.__listAtoms)),
		str(len(self.__listBonds) - len(linesWithMetalBonds)))) # remove all metal bonds
		self.__molFile.write("  0     0  0  0  0  0  0999 V2000")
		self.__molFile.write("\n")
		i = 0
		
		for atomLine in self.__listAtoms:		
			listAtomsColumns = atomLine.split()
		
			atomName = self.__util_.untilPoint(listAtomsColumns[5])
			if atomName == 'Du':
				atomName = '*'
			self.__molFile.write("{:>10}{:>10}{:>10}".format(
			listAtomsColumns[2],
			listAtomsColumns[3],
			listAtomsColumns[4]))
			self.__molFile.write(" ")
			self.__molFile.write(atomName)
			self.__molFile.write("  0  0  0  0  0  0  0  0  0  0  0  0")
			self.__molFile.write("\n")
		
		i = 0
		while i < len(self.__listBonds):
			if i in linesWithMetalBonds:
				i+=1
				continue
			listBondsColumns = self.__listBonds[i].split()
			self.__molFile.write("{:>3}{:>3}".format(
			listBondsColumns[1],
			listBondsColumns[2]))
			self.__molFile.write("  1  0  0  0")
			self.__molFile.write("\n")
			i+=1
			
		self.__molFile.write("M  END\n")
		
	def _writeMol2FileFormat(self):
		linesWithMetalBonds = self._getListOfLinesWithMetalBonds()
		self.__molFile.write('# MOL2 file generated by StereoisomerIdentifier\n\n')
		self.__molFile.write('@<TRIPOS>MOLECULE\n')
		self.__molFile.write(self.__molName + '\n')

		newNumberOfBonds = len(self.__listBonds) - len(linesWithMetalBonds)
		self.__molFile.write("{:>6}{:>6}{:>6}\n".format(
		str(len(self.__listAtoms)),
		str(newNumberOfBonds),# remove all metal bonds
		"0")) 
		self.__molFile.write("SMALL\nNO_CHARGE\n****\nGenerated from StereoisomerIdentifier\n\n")
		self.__molFile.write("@<TRIPOS>ATOM\n")

		printAtoms = False
		printBonds = False
		kBonds = 0
		newKbonds = 1
		
		for line in self.__fileStream:

			if "@<TRIPOS>ATOM" in line:
				printAtoms = True
				continue

			if "@<TRIPOS>BOND" in line:
				self.__molFile.write("@<TRIPOS>BOND\n")
				printBonds = True
				printAtoms = False
				continue
			
			if kBonds == len(self.__listBonds):
				break
			
			if printAtoms:
				self.__molFile.write(line + '\n')
				
			if printBonds:
				if kBonds in linesWithMetalBonds:
					kBonds += 1
					continue
				else:
					kBonds += 1
					lineList = line.split()
					self.__molFile.write("{:>7}{:>5}{:>5} {:>2}\n".format(
					str(newKbonds),
					lineList[1],
					lineList[2],
					lineList[3]))
					newKbonds += 1


	def _openMolFile(self):
		temporaryMolName = self.getTemporaryMolName()
		self.__molFile = open(temporaryMolName + self.__molFormat, "w")
		
	def _readMol2File(self):
		mol2Input = open(self.__fileMol2Name + ".mol2","r")
		self.__fileStream = mol2Input.read().splitlines()
		self._extractMol2Info()
		self._groupBondsInPairs()
		mol2Input.close()

	def _extractMol2Info(self):
		for i in range(len(self.__fileStream)):
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

		for i in range(len(self.__listAtoms)):
			listAtomsColumns = self.__listAtoms[i].split()
			if self.__util_.isMetal(listAtomsColumns[5]):
				self.__metalsInMol2File.append(i)
			
		if len(self.__metalsInMol2File) < 1:
			raise Exception(self.__errorMessages_.getNoMetalError())

	def _checkExtension(self, fileName):
		self.__fileMol2Name = ""
		self.__fileMol2Name, fileExtension = os.path.splitext(fileName)
		if fileExtension != '.mol2':
			raise Exception(self.__errorMessages_.getMolFileError())

	def _getListOfLinesWithMetalBonds(self):
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
		
		return linesWithMetalBonds

	def _groupBondsInPairs(self):
		for i in range(len(self.__listBonds)):
			auxB1 = int(self.__listBonds[i].split()[1])
			auxB2 = int(self.__listBonds[i].split()[2])
			self.__listOfBondsInPairs.append([auxB1,auxB2])