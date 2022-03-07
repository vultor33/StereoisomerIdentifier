from identifier.FormulaHandling import FormulaHandling
from identifier.ErrorMessages import ErrorMessages
from identifier.Graph import Graph

MAXCOORDINATION = 10

class PrioritiesObtainer:
	
	def __init__(self, molFileHandlingObject, externalPriorites):
		self.__errorMessages_ = ErrorMessages()
		self.__molFileHandling_ = molFileHandlingObject
		self.__externalPriorites = externalPriorites

		self.__prioritiesOfMetalI = []
		self.__directFormula_ = FormulaHandling()
		self.__enumerationFormula_ = FormulaHandling()
		self.__chelatesOfMetalI = []
		self.__ligandsBondedToMetalI = []

	def calculateLigandsPriorities(self, iMetal):
		self._findChelations(iMetal)
		self._identifierLimitationsExceptions()

		externalPrioritiesI = self.getExternalPrioritiesOfMetalI()

		self._generateMolecularFormulas(externalPrioritiesI)
		self._defineFinalPriorities(externalPrioritiesI)	

	def calculateOnlyExternalPriorities(self, iMetal):
		self._findChelations(iMetal)
		self._identifierLimitationsExceptions()

	def getPrioritesOfMetalI(self):
		return self.__prioritiesOfMetalI

	def getDirectFormula(self):
		return self.__directFormula_.getFormula()
	
	def getEnumerationFormula(self):
		return self.__enumerationFormula_.getFormula()
	
	def getChelatesOfMetalI(self):
		return self.__chelatesOfMetalI
		
	def getLigandsBondedToMetalI(self):
		return self.__ligandsBondedToMetalI
	
	def getAllPriorities(self):
		return self.__externalPriorites
	
	def getExternalPrioritiesOfMetalI(self):
		priorities = []
		for i in self.__ligandsBondedToMetalI:
			priorities.append(self.__externalPriorites[i-1])
		return priorities
				
	
	def isEnumerationEqualDirectFormula(self):
		return self.getDirectFormula() == self.getEnumerationFormula()
	
	def isMetalOfCoordinationOne(self):
		return len(self.__ligandsBondedToMetalI) == 1

	def _defineFinalPriorities(self, externalPriorities):
		InputPrioritiesToFinalPrioritiesMap = self.__enumerationFormula_.getInputPrioritiesToFinalPrioritiesMap()
		self.__prioritiesOfMetalI = [InputPrioritiesToFinalPrioritiesMap[x] for x in externalPriorities]

	def _generateMolecularFormulas(self, externalPriorities):
		self.__directFormula_.generateMolecularFormula(externalPriorities,self.__chelatesOfMetalI)
		self.__enumerationFormula_.generateEnumerationFormula(externalPriorities,self.__chelatesOfMetalI)


	def _identifierLimitationsExceptions(self):
		if len(self.__ligandsBondedToMetalI) > MAXCOORDINATION:
			raise Exception(self.__errorMessages_.getMaxLigandsError())
		elif len(self.__ligandsBondedToMetalI) == 0:
			raise Exception(self.__errorMessages_.getNoLigandsError())


	def _findChelations(self, iMetal):
		graph_ = self._defineGraphAndLigandsBondedToMetal(iMetal)
		chelates = self._findAllChelatesPositions(graph_)
		self._defineChelatesOfMetalI(chelates)


	def _findAllChelatesPositions(self, graph_):
		chelates = []
		for i in range(len(self.__ligandsBondedToMetalI) - 1):
			if self._isInChelates(i, chelates):
				continue
			auxChelates = self._findAllDonorAtomsBondedToI(i, graph_)
			chelates.append(auxChelates)		
			
		return chelates


	def _defineChelatesOfMetalI(self, chelates):
		auxChelates = [s for s in chelates if len(s) != 1]
		for iChel in auxChelates:
			self.__chelatesOfMetalI.append(iChel)


	def _findAllDonorAtomsBondedToI(self, i, graph_):
			j = i + 1
			auxChelates = [i]
			while j < len(self.__ligandsBondedToMetalI):
				isConnected = graph_.isConnected(self.__ligandsBondedToMetalI[i], self.__ligandsBondedToMetalI[j])
				if isConnected:
					auxChelates.append(j)
				j+=1
				
			return auxChelates


	def _isInChelates(self, i, chelates):
			for chel in chelates:
				if i in chel:
					return True
			
			return False


	def _defineGraphAndLigandsBondedToMetal(self, iMetal):		
		pairsOfBonds = self.__molFileHandling_.getBondsInPairs()
		graph_ = Graph(pairsOfBonds)
		auxLigandsBondedToMetal = graph_.getNodeConnections(iMetal+1) #iMetal starts with 1 in mol file
		for iLigand in auxLigandsBondedToMetal:
			self.__ligandsBondedToMetalI.append(iLigand)
		graph_.remove(iMetal+1)
		return graph_


