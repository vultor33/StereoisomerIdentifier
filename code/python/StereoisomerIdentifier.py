from ErrorMessages import ErrorMessages
from Utilities import Utilities
from ExternalPrioritiesObtainer import ExternalPrioritiesObtainer
from MolFileHandling import MolFileHandling
from CppExecutableHandling import CppExecutableHandling
from BimetallicAnalysis import BimetallicAnalysis

class StereoisomerIdentifier:
	"""Class to identy the stereoisomer of a given organometallic complex"""

	def __init__(self, fileMol2Name):
		self.__errorMessages_ = ErrorMessages()
		self.__util_ = Utilities()
		self.__molFileHandling_ = MolFileHandling(fileMol2Name)
		self.__molFileHandling_.writeMolFile()
		prioritiesObtainer_ = ExternalPrioritiesObtainer(self.__molFileHandling_.getTemporaryMolName())
		self.__equivalenceRank = prioritiesObtainer_.getPriorities()
		self.__cppExec_ = CppExecutableHandling(self.__molFileHandling_, self.__equivalenceRank)
		self.__cppExec_.setKeepIdentierFiles(False)
		self.__bimetallicAnalysis_ = BimetallicAnalysis(self.__molFileHandling_, self.__equivalenceRank)

	def activateKeepIdentifier(self):
		self.__cppExec_.setKeepIdentierFiles(True)

	def runStereoisomerIdentifierRmsd(self, outputFile_):
		try:
			bimetallicWriteSummary = self.__bimetallicAnalysis_.checkIfItIsDimmetallic()
			outputFile_.write(bimetallicWriteSummary)
			
		except:
			outputFile_.write("Mix;")

		for iMetal in self.__molFileHandling_.getMetalsInMol2FileList(): #Evaluate the stereoisomer of each metal
			try:
				metalName = self.__molFileHandling_.getListAtoms()[iMetal].split()[1] + ";"
				outputFile_.write(metalName)
				cppResultSummary = self.__cppExec_.runCppCode(iMetal)
				outputFile_.write(cppResultSummary)
					
			except Exception as e:
				exceptionWriteMessage = self._exceptionWriteMessage(str(e))
				outputFile_.write(exceptionWriteMessage) 


	def _exceptionWriteMessage(self, exceptionString):
		if exceptionString == self.__errorMessages_.getNoMetalError():
			return "E.NoMetal;;;"
		elif exceptionString == self.__errorMessages_.getNoLigandsError():
			return "E.NoLigands;;;"
		elif exceptionString == self.__errorMessages_.getMaxLigandsError():
			return "E.NL>8;;;"
		elif exceptionString == self.__errorMessages_.getChelationDefinitionError():
			return "E.Formula;;;"
		elif exceptionString == self.__errorMessages_.getGraphError():
			return "E.Graph;;;"
		else:
			return "E.{};;;".format(exceptionString)

if __name__ ==  "__main__":
	print("Module to handle Mol files")
