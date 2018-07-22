from FormulaHandling import FormulaHandling
from Utilities import Utilities
from ErrorMessages import ErrorMessages
import itertools
import math


class AllMolecularFormulasGenerator:
	"""Class to enumerate all possible formulas with all possible chelations
	
    Parameters
    ----------
	Coordination number\n
	File stream

    Returns
    -------
	All molecular formulas for a given coordination number

    Example
    -------
	objG = AllMolecularFormulasGenerator()\n
	objG.generateAllFormulas(coordinationNumber)\n
	objG.printAllFormulasToStream(fileStream)

    Note
    -------
	This file comes with two additional functions:\n
	useAllMolecularFormulasGenerator(nCoordinationMax)\n
	testAllMolecularFormulasGenerator(nCoordinationMax, fileReference)
	"""
	def __init__(self):
		self.__errorMessages_ = ErrorMessages()
		self.__OUTPUT_FILE_SUGGESTION = "DataAllFormulas.py" # need to be open outside the object
		self.__coordinationNumber = 0
		self.__allFormulas = []
		self.__allEnumerationFormulas = []
		self.__allReferences = []
		self.__allChelations = []

	def getOutputFileName(self):
		return self.__OUTPUT_FILE_SUGGESTION

	def generateAllFormulas(self, coordinationNumber):
		self.__coordinationNumber = coordinationNumber
		self._generateAllMonodentatesFormulas()
		self._generateAllChelatesFormulas()

	def printAllFormulasToStream(self, allFormFile):
		tempListName ="aux" + str(self.__coordinationNumber) +"FormList"
		allFormFile.write(tempListName +" = []\n")
		for i in range(len(self.__allFormulas)):
			allFormFile.write("{}.append([{},\'{}\',\'{}\',{},{}])\n".format(
			tempListName,
			self.__coordinationNumber,
			self.__allFormulas[i],
			self.__allEnumerationFormulas[i],
			self.__allReferences[i],
			self.__allChelations[i]))
		allFormFile.write("allFormList.append(" + tempListName + ")\n")

	def _generateAllMonodentatesFormulas(self):
		rangeToDefineCombinatios = list(range(0,self.__coordinationNumber))
		allPossibleCombinations = list(itertools.combinations_with_replacement(rangeToDefineCombinatios, self.__coordinationNumber))
		for combination in allPossibleCombinations:
			self._addFormula(combination,[])

	def _generateAllChelatesFormulas(self):
		allChelationsCombinations = self._generateAllChelateCombinations()
		allMonodentatePriorities = list(self.__allReferences)
		for monodentatePriorities in allMonodentatePriorities:  #use monodentates as base to chelations definitions
			for chelCombination in allChelationsCombinations:
				allChelations = self._defineAllChelations(chelCombination)
				for chelations in allChelations:
					self._addSafeChelateFormula(monodentatePriorities,chelations)

	def _addFormula(self, priorities, chelations):
			objP = FormulaHandling()
			objP.generateMolecularFormula(priorities, chelations)
			if not objP.getFormula() in self.__allFormulas:
				self.__allFormulas.append(objP.getFormula())
				self.__allReferences.append(objP.getReferenceLine())
				self.__allChelations.append(objP.getCanonChelation())
				self._addEnumerationFormula(priorities,chelations)

	def _addSafeChelateFormula(self,priorities,chelations):
		try:
			self._addFormula(priorities, chelations)

		except Exception as e:
			if str(e) == self.__errorMessages_.getChelationDefinitionError():
				pass
			else:
				print(self.__errorMessages_.getUnexpecterError(str(e)))

	def _addEnumerationFormula(self, priorities,chelations):
			objE = FormulaHandling()
			objE.generateEnumerationFormula(priorities, chelations)
			self.__allEnumerationFormulas.append(objE.getFormula())

	def _generateAllChelateCombinations(self):
		allChelationsCombinations = []
		maxPossibleNumberOfChelates = int(math.floor(self.__coordinationNumber/2))
		for numberOfChelates in range(1,maxPossibleNumberOfChelates+1):
			genericCombinations = list(itertools.combinations_with_replacement(list(range(2,self.__coordinationNumber + 1)), numberOfChelates))
			for chelI in genericCombinations:
				sumOfChelateCoordination = sum(chelI)
				if sumOfChelateCoordination <= self.__coordinationNumber:
					allChelationsCombinations.append(chelI)

		return allChelationsCombinations

	def _defineAllChelations(self, chelCombination):
		chelList = []
		for chelCombinationI in chelCombination:
			chelReach = list(range(0,self.__coordinationNumber))
			chelList.append(list(itertools.combinations(chelReach, chelCombinationI)))

		allChelations = list(itertools.product(*chelList))
		return allChelations


################################################################################################################################################################

def useAllMolecularFormulasGenerator(nCoordinationMax):
	allFormulas_ = AllMolecularFormulasGenerator()
	fileName = allFormulas_.getOutputFileName()
	allFormFile = open(fileName, "w")
	allFormFile.write("allFormList = [[0]]\n")
	for i in range(1,nCoordinationMax + 1):
		objG = AllMolecularFormulasGenerator()
		objG.generateAllFormulas(i)
		objG.printAllFormulasToStream(allFormFile)
	allFormFile.close()

def testAllMolecularFormulasGenerator():
	print('TESTING AllMolecularFormulasGenerator')

	useAllMolecularFormulasGenerator(6)

	util_ = Utilities()
	allFormulas_ = AllMolecularFormulasGenerator()
	outputFileName = allFormulas_.getOutputFileName()
	testFilePath = "TestFiles\\" + outputFileName + '-reference'
	
	if util_.isOutOfOrderFilesEqual(outputFileName,testFilePath):
		print('PASSED')
	else:
		print('NOT PASSED')


if __name__ ==  "__main__":
	testAllMolecularFormulasGenerator()






