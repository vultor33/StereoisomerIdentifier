import os
from ErrorMessages import ErrorMessages
from PrioritiesObtainer import PrioritiesObtainer
from shutil import copyfile
import subprocess

class CppExecutableHandling:
	"""Class to handle cpp executable and its files\n
		This class call PrioritiesObtainer.py for priorites and formulas\n
		Build cpp file\n
		Execute\n
		Read the results\n
		And return an appropriate summary
	
    Options
    ----------
	Activate setKeepIdentierFiles to see the files that go to executable
	"""
	
	def __init__(self, molFileHandlingObject, equivalenceRank):
		self.__errorMessages_ = ErrorMessages()
		self.__molFileHandling_ = molFileHandlingObject
		self.__equivalenceRank = equivalenceRank
		self.__keepIdentifierFiles = False
		
		self.__cppRmsdFailureFlag = "rmsdfailed"
		self.__cppPolyhedronFailureFlag = "failed"

	def __del__(self):
		if os.path.isfile(self.__molFileHandling_.getBaseFileName() + "-cpp.inp"):
			os.remove(self.__molFileHandling_.getBaseFileName() + "-cpp.inp")
		if os.path.isfile(self.__molFileHandling_.getBaseFileName() + "-cpp.inp.log"):
			os.remove(self.__molFileHandling_.getBaseFileName() + "-cpp.inp.log")

	def setKeepIdentierFiles(self, haveToKeep):
		self.__keepIdentifierFiles = haveToKeep

	def runCppCode(self, iMetal):
		prior_ = self._generatePriorities(iMetal)
		
		if prior_.isMetalOfCoordinationOne():
			return "a;a-L-1;0;"
		
		self._writeCppInput(iMetal, prior_)
		execCommand = "./StereoisomerIdentifierRmsd.exe " + self.__molFileHandling_.getBaseFileName() + "-cpp.inp"
		subprocess.call(execCommand, shell=True)
		cppResultSummary = self._readCppOutput(prior_)

		return cppResultSummary
		

	def _readCppOutput(self, prior_):
		cppOutput_ = open(self.__molFileHandling_.getBaseFileName() + "-cpp.inp.log","r")
		cppStream_ = cppOutput_.read().splitlines()
		cppOutput_.close()
		cppResultSummary = prior_.getDirectFormula() + ";"
		if len(cppStream_) == 0:
			raise Exception(self.__errorMessages_.getCppFileMissing())
		if cppStream_[1] == self.__cppRmsdFailureFlag:
			cppResultSummary += "E.RMSD;" + cppStream_[2] + ";"
		elif cppStream_[1] == self.__cppPolyhedronFailureFlag:
			cppResultSummary += prior_.getEnumerationFormula() + "-E.Polyedron;" + cppStream_[2] + ";"
		elif len(cppStream_) > 0:
			cppResultSummary += prior_.getEnumerationFormula() + "-" + cppStream_[1] + ";" + cppStream_[2] + ";"

		return cppResultSummary

	def _writeCppInput(self, iMetal, prior_):
		cppInput = open(self.__molFileHandling_.getBaseFileName() + "-cpp.inp", "w")

		self._writeTitle(cppInput, prior_)
		metalCoordinates = self.__molFileHandling_.getListAtoms()[iMetal].split()
		self._writeAtomFormat(cppInput,metalCoordinates,'-1',self.__equivalenceRank[iMetal])

		kIndex = 0
		for i in prior_.getLigandsBondedToMetalI():
			donorAtomsCoordinates = self.__molFileHandling_.getListAtoms()[i-1].split()
			localPriority = prior_.getPrioritesOfMetalI()[kIndex]
			globalPriority = self.__equivalenceRank[i-1]
			self._writeAtomFormat(cppInput,donorAtomsCoordinates,localPriority,globalPriority)
			kIndex += 1
			
		cppInput.write("end\n")
		cppInput.close()
		self._copyFileIfItIsToKeep(iMetal)


	def _generatePriorities(self, iMetal):
		prior_ = PrioritiesObtainer(self.__molFileHandling_, self.__equivalenceRank)
		prior_.calculateLigandsPriorities(iMetal)
		self.__directFormula = prior_.getDirectFormula()
		self.__enumerationFormula = prior_.getEnumerationFormula()
		return prior_
		

	def _copyFileIfItIsToKeep(self, iMetal):
		if self.__keepIdentifierFiles:
			metalStr = self.__molFileHandling_.getListAtoms()[iMetal].split()[1]
			fileCppNameKeep = self.__molFileHandling_.getBaseFileName() + "-" + metalStr + "-cpp.inp"
			copyfile(self.__molFileHandling_.getBaseFileName() + "-cpp.inp", fileCppNameKeep)

	def _writeTitle(self, cppInput, prior_):
		if prior_.isEnumerationEqualDirectFormula():
			cppInput.write(prior_.getDirectFormula() + "\n")
			cppInput.write("chelates:  {}".format(len(prior_.getChelatesOfMetalI())))
			for chel in prior_.getChelatesOfMetalI():
				cppInput.write("  cI-length:  {}  cI:".format(len(chel)))
				for chelI in chel:
					cppInput.write(" {} ".format(chelI))
			cppInput.write("\n")
		else:
			cppInput.write(prior_.getEnumerationFormula() + "\n\n")

	def _writeAtomFormat(self, cppInput, atomsList, localPriority, globalPriority):
			cppInput.write("{:>5}{:>10}{:>10}{:>10}{:>5}{:>5}\n".format(
					atomsList[1],
					atomsList[2],
					atomsList[3],
					atomsList[4],
					localPriority,
					globalPriority))
	
