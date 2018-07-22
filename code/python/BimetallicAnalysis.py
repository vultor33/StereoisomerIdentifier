from Utilities import Utilities
from Graph import Graph
from PrioritiesObtainer import PrioritiesObtainer


class BimetallicAnalysis:
	""" Evaluate if the unitary cell have a bimetallic complex """

	def __init__(self, molFileHandling_, equivalenceRank):
		self.__molFileHandling_ = molFileHandling_
		self.__equivalenceRank = equivalenceRank
		self.__util_ = Utilities()
		pairsOfBonds = self.__molFileHandling_.getBondsInPairs()
		self.__graph_ = Graph(pairsOfBonds)

	def checkIfItIsDimmetallic(self):
		metalsInMol2File = self.__molFileHandling_.getMetalsInMol2FileList()
		listAtoms = self.__molFileHandling_.getListAtoms()
		if len(metalsInMol2File) == 2:
			metal1 = self.__util_.untilPoint(listAtoms[metalsInMol2File[0]].split()[5])
			metal2 = self.__util_.untilPoint(listAtoms[metalsInMol2File[1]].split()[5])
			metalIndex1 = int(float(metalsInMol2File[0])) + 1
			metalIndex2 = int(float(metalsInMol2File[1])) + 1
			metalConnected = False
			try:
				metalConnected = self.__graph_.isConnected(metalIndex1,metalIndex2)
			except Exception:
				return "Mix;"
			
			if metal1 == metal2 and metalConnected:
				if self._isBimetallicComplexesHaveSameExternalPriorities(metalsInMol2File[0],metalsInMol2File[1]):
					return "Dimetal;"
				else:
					return "DimetalNF;"
			else:
				return "Mix;"
		else:
			return "Mix;"
		
	def _isBimetallicComplexesHaveSameExternalPriorities(self, iMetal1, iMetal2):
		priorities1 = self._getOriginalPriorities(iMetal1)
		priorities2 = self._getOriginalPriorities(iMetal2)
		priorities1.sort()
		priorities2.sort()
		return priorities1 == priorities2

	def _getOriginalPriorities(self, iMetal):
		prior_ = PrioritiesObtainer(self.__molFileHandling_, self.__equivalenceRank)
		prior_.calculateOnlyExternalPriorities(iMetal)
		return prior_.getExternalPrioritiesOfMetalI()
		
