import glob
from ErrorMessages import ErrorMessages
from StereoisomerIdentifier import StereoisomerIdentifier
from Utilities import Utilities


class RunStereoisomerIdentifier:	
	"""Receives a list of files and identify the stereisomer of all of them
	
    Notes
    ----------
	See testRunStereoisomerIdentifier below for uses
	"""
	
	def __init__(self, pathInputFiles, pathOutputFiles, extension):
		self.__errorMessages_ = ErrorMessages()
		self.__HEADER = "CSD;Info;Metal;Formula;ID;RMSD"
		self.__pathInputFiles = pathInputFiles
		self.__OUTPUT_FILE_NAME = pathOutputFiles + "calculating.csv"
		self.__extension = extension
		self.__calcFiles = []
		self._openCalculatingFile()

	def __del__(self):
		if not self.__calculatingFile.closed:
			self.__calculatingFile.close()

	def getOutputFileName(self):
		return self.__OUTPUT_FILE_NAME
		
	def runFromList(self, csdRefcodesList):
		self._buildInputFileNames(csdRefcodesList, self.__extension)
		self._calcAllFiles()

	def runAllFilesFromInputDirectory(self):
		self.__calcFiles = glob.glob(self.__pathInputFiles + '*' + self.__extension)
		self._calcAllFiles()

	def runAllFilesFromInputDirectoryWithLimits(self,indexInit,indexFinal):
		tempFileNames = glob.glob(self.__pathInputFiles + '*' + self.__extension)
		self.__calcFiles = tempFileNames[indexInit:indexFinal]
		del tempFileNames
		self._calcAllFiles()
		

	def _calcAllFiles(self):
		for molFile in self.__calcFiles:
			self.__calculatingFile.write("\n{:<9};".format(molFile.partition(".")[0]))
			self._calcFile(molFile)

		self.__calculatingFile.close()

	def _calcFile(self, molFile):
			try:
				obj1 = StereoisomerIdentifier(molFile)
				#obj1.keepIdentifierFiles() fornecer opcoes para ativar ou desativar isso aqui
				obj1.runStereoisomerIdentifierRmsd(self.__calculatingFile)

			except Exception as e:
				if str(e) == self.__errorMessages_.getCipApplicationError():
					self.__calculatingFile.write("E.rdkit")  # UNIFICAR TODAS AS MENSAGENS DE ERRO
				elif str(e) == self.__errorMessages_.getNoMetalError():
					self.__calculatingFile.write("E.0Metals")
				else:
					self.__calculatingFile.write("E.{}".format(str(e)))
	
	def _buildInputFileNames(self, csdRefcodesList, extension):
		for csdRefcode in csdRefcodesList:
			self.__calcFiles.append(self.__pathInputFiles + csdRefcode + extension)

	def _openCalculatingFile(self):
		self.__calculatingFile = open(self.__OUTPUT_FILE_NAME, "w")
		self.__calculatingFile.write(self.__HEADER)
		
############################################################################################################
	

def testRunStereoisomerIdentifier():
	print('TESTING RunStereoisomerIdentifier')
	pathInput = "TestFiles\\"
	pathOutput = "TestFiles\\"
	extension = ".search1.mol2"
	calcFilesTemp = [
			"AVIYAH",
			"TALWOS",
			"BOGNIV01",
			"FAHWAM",
			"ITADAK",
			"LILHIZ",
			"NAGWIB",
			"OTETOY",
			"PIDYOR",
			"QEZCUU",
			"VACGUC01"]
	run_ = RunStereoisomerIdentifier(pathInput,pathOutput, extension)
	run_.runFromList(calcFilesTemp)
	
	util_ = Utilities()
	if util_.isOrderedFilesEqual(run_.getOutputFileName(),run_.getOutputFileName() + '-reference'):
		print('PASSED')
	else:
		print('NOT PASSED')

if __name__ ==  "__main__":
	testRunStereoisomerIdentifier()
	

	
	
	
	
	