

class ErrorMessages:
	def __init__(self):
		self.__IOErrorMessage = "Could not open:  "
		self.__cipApplicationError = "Chem.MolFromMol2File failed"
		self.__molFileError = "Couldnt open mol file"
		self.__noMetalError = "Metal number error - 0 metals"
		self.__noLigandsError = "No ligands"
		self.__maxLigandsError = "Number of ligands error"
		self.__chelationDefinitionError = "chelations not well defined"
		self.__molecularFormulaError = 'Error on molecular formula conflicts'
		self.__molecularFormulasConflictTrade = 'Error on formulasConflictTrade - check input for errors'
		self.__redefineChelateTypes = 'Error on _redefineChelatesTypes: new ranks couldnt be set'
		self.__cppFileMissing = "couldnt find cpp log file"
		self.__graphError = "path error"
		
	def getGraphError(self):
		return self.__graphError
		
	def getIOError(self):
		return self.__IOErrorMessage
	
	def getCipApplicationError(self):
		return self.__cipApplicationError
	
	def getNoMetalError(self):
		return self.__noMetalError
	
	def getNoLigandsError(self):
		return self.__noLigandsError
	
	def getMaxLigandsError(self):
		return self.__maxLigandsError
	
	def getChelationDefinitionError(self):
		return self.__chelationDefinitionError
	
	def getMolecularFormulaError(self):
		return self.__molecularFormulaError

	def getFirstFormulaHasPreferenceOverTheSecondError(self):
		return self.__molecularFormulasConflictTrade

	def getRedefineChelateTypesError(self):
		return self.__redefineChelateTypes

	def getMolFileError(self):
		return self.__molFileError
		
	def getCppFileMissing(self):
		return self.__cppFileMissing

	def getUnexpecterError(self, error):
		unexpected = 'Unexpected error: ' + error + '\n Please, contact developers'
		return unexpected
		
	