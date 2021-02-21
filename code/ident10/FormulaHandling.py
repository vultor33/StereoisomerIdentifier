from collections import Counter
import operator
from ErrorMessages import ErrorMessages
from Utilities import Utilities

def applyNewCIP(priorities, chelations):
    fhandling_ = FormulaHandling()
    fhandling_.generateMolecularFormula(priorities, chelations)
    priorMap = fhandling_.getInputPrioritiesToFinalPrioritiesMap()
    newPriorities = [priorMap[x] for x in priorities]
    return newPriorities

class FormulaHandling:
	"""Class to generate and transform molecular formulas

    Parameters
    -------
	Set of priorities\n
	Set of chelations

    Returns
    -------
	Molecular formula\n
	Final priorities\n
	Final chelations\n
	Relation of input priorities to final priorities
    Notes
    ----------
	Use testAllMolecularFormulasGenerator() and testRunStereoisomerIdentifier() for tests
	"""

	def __init__(self):
		self.__errorMessages_ = ErrorMessages()
		self.__util_ = Utilities()

		self.__inputPriorities = []
		self.__inputChelations = []

		#Ex: received: rank: [10,9,14,45,45]. chelations: [[1,2],[3,4]]
		self.__finalFormula = ''                            # Molecular formula. Ex: Ma(A2)(AB)
		self.__finalPriorites = []                     # Line of priority types. Ex: [0,1,1,2,3]
		self.__finalChelations = []                          # Chelation related to the referenceLineVector. Ex: [[1,2],[3,4]]

		# These are relations of canon types and received types
		self.__inputPrioritiesToFinalPrioritiesMap = {}                # ranks were transformed. this dict relates old types to new types.
		self.__priorityOfEachLigand = {}                    # Rank of each type and chelate. Ex: {10:[0],9:[1,2]}
		self.__formulaOfEachLigand = {}                     # Ex: {9:'ab',10:'a'}
		self.__chelationOfEachLigand = {}                   # Keys: first type of the chelating group. values: all other types ordered.  Ex: {1:[1,2]}
		self.__amountOfEachLigand = {}                      # Polidentates are counted as one and it's first type is used. Ex:{0:1,1:1}

		self.__alphabet = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','e','s','t','u','v','x','w','y','z']
		
	def getFormula(self):
		return self.__finalFormula
		
	def getReferenceLine(self):
		return self.__finalPriorites
	
	def getCanonChelation(self):
		return self.__finalChelations
		
	def getInputPrioritiesToFinalPrioritiesMap(self):
		return self.__inputPrioritiesToFinalPrioritiesMap	

	def generateEnumerationFormula(self, priorities, chelations):
		self.__inputPriorities = priorities
		self.__inputChelations = chelations

		isDuplicate = self._checkIfThereAreEqualFormulaChelates()

		if isDuplicate:
			self.generateMolecularFormula(priorities, chelations)
		else:
			self.generateMolecularFormula(priorities, [])


	def generateMolecularFormula(self, rank, chelations):
		self.__inputPriorities = rank
		self.__inputChelations = chelations

		self._countLigandsAndDefineDictChelations()

		self._defineFormulaOfEachLigand()

		dictPriorities = self._rankFormulasByItsCanonicalPriorities() # Apply the new 3 CIP rules on each formula to rank them

		self._defineFinalFormula(dictPriorities)

		self._defineAllFinalPrioritiesAndChelations(dictPriorities)

		self._redefineChelatesTypes() # Until here we have an hierarchy on priorities - mono, bidentate, tridentate and etc. This function changes to: numbers come first. All previous hierarchy is keeped.


##################################################################################################


	def _definePriorityOfEachLigand(self, currentPriority, canonicalPriority):
		for iCan in range(len(canonicalPriority)):
			canonicalPriority[iCan] += currentPriority
		return canonicalPriority

	def _defineAllFinalPrioritiesAndChelations(self, dictPriorities):
		chelationBitePosition = 0
		currentPriorityIndex = 0
		for i in range(len(dictPriorities.keys())):
			keyI = self.__util_.dictFindKeyByValue(dictPriorities,i)
			currentPriority = self._obtainCanonicalPriorityByItsFormula(self.__formulaOfEachLigand[keyI])
			lastPriority = currentPriority[-1]
			
			self.__priorityOfEachLigand[keyI] = self._definePriorityOfEachLigand(currentPriorityIndex, currentPriority)
			self._defineFinalPrioritiesAndChelations(keyI, chelationBitePosition)

			chelationBitePosition += len(self.__amountOfEachLigand[keyI]*self.__priorityOfEachLigand[keyI])
			currentPriorityIndex += lastPriority + 1
	
	def _defineFinalPrioritiesAndChelations(self, keyI, chelationBitePosition):
		for i in  range(self.__amountOfEachLigand[keyI]):
			chelTemp = []
			for iCan in self.__priorityOfEachLigand[keyI]:
				self.__finalPriorites += [iCan]
				chelTemp.append(chelationBitePosition)
				chelationBitePosition+=1
			if len(chelTemp) > 1:
				self.__finalChelations.append(chelTemp)
	

##################################################################################################

	def _countLigandsAndDefineDictChelations(self):
		self._checkIfChelationsPointToSameDonorAtom()
		self._countChelatedLigandsAndAdvanceDictChelations()
		self._countMonodentateLigandsAndAdvanceDictChelations()
		self._checkIfThereAreDifferentChelatedLigandsWithSameTypes()

	def _checkIfChelationsPointToSameDonorAtom(self):
		i = 0
		j = 0
		while i < len(self.__inputChelations) - 1:
			j = i + 1
			while j < len(self.__inputChelations):
				for chelI in self.__inputChelations[i]:
					for chelJ in self.__inputChelations[j]:
						if chelI == chelJ:
							raise Exception(self.__errorMessages_.getChelationDefinitionError())
				j+=1
			i+=1

	def _countChelatedLigandsAndAdvanceDictChelations(self):
		for chel in self.__inputChelations:
			lisRankChel = []
			for iChel in chel:
				lisRankChel.append(self.__inputPriorities[iChel])
			lisRankChel.sort()
			if lisRankChel[0] in self.__chelationOfEachLigand:
				if self.__chelationOfEachLigand[lisRankChel[0]] != lisRankChel: #different chelations and same types arent allowed
					raise Exception(self.__errorMessages_.getChelationDefinitionError())
				self.__amountOfEachLigand[lisRankChel[0]] += 1
			else:
				self.__chelationOfEachLigand[lisRankChel[0]] = lisRankChel
				self.__amountOfEachLigand[lisRankChel[0]] = 1

	def _countMonodentateLigandsAndAdvanceDictChelations(self):
		allChelPositios = []
		for chelList in self.__inputChelations:
			for chelI in chelList:
				allChelPositios += [chelI]

		for iAux in range(len(self.__inputPriorities)):
			if iAux in allChelPositios:
				continue
			if self.__inputPriorities[iAux] in self.__chelationOfEachLigand:
				if self.__chelationOfEachLigand[self.__inputPriorities[iAux]] != [self.__inputPriorities[iAux]]:
					raise Exception(self.__errorMessages_.getChelationDefinitionError())
				self.__amountOfEachLigand[self.__inputPriorities[iAux]] += 1
			else:
				self.__chelationOfEachLigand[self.__inputPriorities[iAux]] = [self.__inputPriorities[iAux]]
				self.__amountOfEachLigand[self.__inputPriorities[iAux]] = 1

	def _checkIfThereAreDifferentChelatedLigandsWithSameTypes(self):
		for keyChel in self.__chelationOfEachLigand:
			for chelRank in self.__chelationOfEachLigand[keyChel]:
				for keyChel2 in self.__chelationOfEachLigand:
					if keyChel == keyChel2:
						continue
					for chelRank2 in self.__chelationOfEachLigand[keyChel2]:
						if chelRank == chelRank2:
							raise Exception(self.__errorMessages_.getChelationDefinitionError())

		for chel in self.__inputChelations:
			lisRankChel = []
			for iChel in chel:
				lisRankChel.append(self.__inputPriorities[iChel])
			lisRankChel.sort()
			if lisRankChel[0] in self.__chelationOfEachLigand:
				if self.__chelationOfEachLigand[lisRankChel[0]] != lisRankChel:
					raise Exception(self.__errorMessages_.getChelationDefinitionError())

##################################################################################################

	def _checkIfThereAreEqualFormulaChelates(self):
		priorityOfDonorAtoms = self._calculatePrioritiesOfEachDonorAtomAtChelation(self.__inputPriorities, self.__inputChelations)
		duplicate = False
		for i in range(len(self.__inputChelations)-1):
			for j in range(len(self.__inputChelations)):
				if i == j:
					continue
				for chelI in self.__inputChelations[i]:
					for chelJ in self.__inputChelations[j]:
						if priorityOfDonorAtoms[chelI] == priorityOfDonorAtoms[chelJ]:
							duplicate = True

		return duplicate

	def _calculatePrioritiesOfEachDonorAtomAtChelation(self, priorities, chelations):
		prioritiesOfDonorAtoms = {}
		for chel in chelations:
			for chelI in chel:
				prioritiesOfDonorAtoms[chelI] = priorities[chelI]
		
		return prioritiesOfDonorAtoms
	
##################################################################################################

	def _defineFormulaOfEachLigand(self):
		for key in self.__chelationOfEachLigand:
			self.__formulaOfEachLigand[key] = self._calculateFormula(self.__chelationOfEachLigand[key])

		listSorted = False  #bubble method
		while not listSorted: 
			listSorted = True
			for key1 in self.__formulaOfEachLigand:
				for key2 in self.__formulaOfEachLigand:
					if key1 == key2:
						continue
					elif self.__formulaOfEachLigand[key1] != self.__formulaOfEachLigand[key2]:
						continue
					else:
						listSorted = self._solveFormulaPriorityConflict(key1,key2)

	def _solveFormulaPriorityConflict(self, key1, key2):
		if self.__amountOfEachLigand[key1] > self.__amountOfEachLigand[key2]:
			self.__formulaOfEachLigand[key2] = self._advanceFormula(self.__formulaOfEachLigand[key2])
			return False
		elif self.__amountOfEachLigand[key1] < self.__amountOfEachLigand[key2]:
			self.__formulaOfEachLigand[key1] = self._advanceFormula(self.__formulaOfEachLigand[key1])
			return False
		else:
			if key1 < key2:
				self.__formulaOfEachLigand[key2] = self._advanceFormula(self.__formulaOfEachLigand[key2])
				return False
			elif key1 > key2:
				self.__formulaOfEachLigand[key1] = self._advanceFormula(self.__formulaOfEachLigand[key1])
				return False
			else:
				raise Exception(self.__errorMessages_.getMolecularFormulaError())
		
		return True

	#from formula: last term must be the largest and first the lowest
	def _advanceFormula(self,formula):
		types = self.__util_.justLetters(formula)
		deltaIndex = self.__alphabet.index(types[len(types)-1]) - self.__alphabet.index(types[0]) + 1
		newFormula = ''
		for iFor in formula:
			if iFor.isdigit():
				newFormula += iFor
				continue
			newFormula += self.__alphabet[deltaIndex + self.__alphabet.index(iFor)]
	
		return newFormula


	def _calculateFormula(self, rank):
		rankCounting = Counter(rank)
		rankCountingElems = []
		for comp in rankCounting:
			rankCountingElems.append(rankCounting[comp])

		rankCountingElems.sort()
		ligandsAmount = list(rankCountingElems)
		ligandsAmount.reverse()
		ligandTypes = self.__alphabet[:len(ligandsAmount)]
		molecularFormula = ""
		i = 0		
		while i < len(ligandTypes):
			if ligandsAmount[i] == 1:
				molecularFormula += ligandTypes[i]
			else:
				molecularFormula += ligandTypes[i] + str(ligandsAmount[i])
			i+=1
		return molecularFormula

##################################################################################################


	def _rankFormulasByItsCanonicalPriorities(self):
		dictPriorities = {}  #Priorities of each key, add 1 if it is bidentate and etc. This are used to define self.__priorityOfEachLigand
		for key in self.__formulaOfEachLigand:
			dictPriorities[key] = 0

		# bubble method to define priorities. All starts at zero, changes are decided by self._formulasConflictTrade
		trade  = True
		while trade:
			trade = False
			for key1 in self.__formulaOfEachLigand:
				for key2 in self.__formulaOfEachLigand:
					if key1 == key2:
						continue
					if dictPriorities[key1] == dictPriorities[key2]:
						trade = True
						if self._firstFormulaHasPreferenceOverTheSecond(self.__formulaOfEachLigand[key1],self.__formulaOfEachLigand[key2]):
							dictPriorities[key1] += 1
						else:
							dictPriorities[key2] += 1
							
		return dictPriorities

	def _firstFormulaHasPreferenceOverTheSecond(self, formula1, formula2):
		canon1 = self._obtainCanonicalPriorityByItsFormula(formula1)
		canon2 = self._obtainCanonicalPriorityByItsFormula(formula2)
	
		if len(canon1) != len(canon2):
			return len(canon1) > len(canon2)
		else:
			i = 0
			while i < len(canon1):
				if canon1[i] != canon2[i]:
					return canon1[i] > canon2[i]
				i+=1
			
			if self.__alphabet.index(formula1[0]) != self.__alphabet.index(formula2[0]):
				return self.__alphabet.index(formula1[0]) > self.__alphabet.index(formula2[0])
			else:
				raise Exception(self.__errorMessages_.getFirstFormulaHasPreferenceOverTheSecondError())


	def _obtainCanonicalPriorityByItsFormula(self, formula):
		k = 0
		rankCanon = []
		for i in range(len(formula)):
			if formula[i].isdigit():
				for iR in range(int(float(formula[i]))-1):
					rankCanon.append(k-1)
			else:
				rankCanon.append(k)
				k+=1

		return rankCanon
	
##################################################################################################

	def _defineFinalFormula(self, dictPriorities):
		for i in range(len(dictPriorities.keys())):
			keyI = self.__util_.dictFindKeyByValue(dictPriorities,i)
			if len(self.__formulaOfEachLigand[keyI]) > 1:
				self.__finalFormula += '(' + str(self.__formulaOfEachLigand[keyI]).upper() + ')'
			else:
				self.__finalFormula += str(self.__formulaOfEachLigand[keyI])
			if self.__amountOfEachLigand[keyI] > 1:
				self.__finalFormula += str(self.__amountOfEachLigand[keyI])

##################################################################################################


	def _redefineChelatesTypes(self):
		prioritiesOfDonorAtoms = self._calculatePrioritiesOfEachDonorAtomAtChelation(self.__finalPriorites,self.__finalChelations)
		
		prioritiesTransformMap = self._definePrioritiesTransformMap()
		
		self._redefinePriorityOfEachLigand(prioritiesTransformMap)

		self._defineInputPrioritiesToFinalPrioritiesMap()

		self._redefineFinalPrioritiesAndChelations(prioritiesTransformMap, prioritiesOfDonorAtoms)


	def _redefineFinalPrioritiesAndChelations(self, prioritiesTransformMap, prioritiesOfDonorAtoms):
		amountOfEachPriority = Counter(self.__finalPriorites)
		
		auxTypes = []
		for key in amountOfEachPriority:
			for i in range(amountOfEachPriority[key]):
				auxTypes.append(prioritiesTransformMap[key])

		auxTypes.sort()
		for i in range(len(auxTypes)):
			self.__finalPriorites[i] = auxTypes[i]
		for i in range(len(self.__finalChelations)):
			for j in range(len(self.__finalChelations[i])):
				self.__finalChelations[i][j] = auxTypes.index(prioritiesTransformMap[prioritiesOfDonorAtoms[self.__finalChelations[i][j]]])
				auxTypes[self.__finalChelations[i][j]] = -1
		

	def _defineInputPrioritiesToFinalPrioritiesMap(self):
		#Initial type: self.__chelationOfEachLigand[keyI] (but I have to order here
		#Final type:  self.__priorityOfEachLigand[keyI][i]. 
		for keyI in self.__chelationOfEachLigand:
			#dictChelations comes like this: [3, 9, 9, 10, 11]. I need repeating priorities to come first.
			mapBetweenChelatesAndTypes = self._calculateNewMapBetweenAtomsAndTypes(self.__chelationOfEachLigand[keyI])
			rankInTheInteriorOfTheChelate = list(self.__chelationOfEachLigand[keyI])
			for iRank in range(len(rankInTheInteriorOfTheChelate)):
				rankInTheInteriorOfTheChelate[iRank] = mapBetweenChelatesAndTypes[rankInTheInteriorOfTheChelate[iRank]]
			rankInTheInteriorOfTheChelate.sort()
			for i in range(len(self.__chelationOfEachLigand[keyI])):
				finalValue = self.__priorityOfEachLigand[keyI][i]
				internal = rankInTheInteriorOfTheChelate[i]
				internalToCipRank = self.__util_.dictFindKeyByValue(mapBetweenChelatesAndTypes,internal)

				if internalToCipRank in self.__inputPrioritiesToFinalPrioritiesMap:
					if self.__inputPrioritiesToFinalPrioritiesMap[internalToCipRank] != finalValue:
						raise Exception(self.__errorMessages_.getRedefineChelateTypesError())
	
				self.__inputPrioritiesToFinalPrioritiesMap[internalToCipRank] = finalValue
		

	def _redefinePriorityOfEachLigand(self, rankTypesTransformMap):
		auxDictCanonRanks = {}
		for keyRanks in self.__priorityOfEachLigand:
			newRankRelation = []
			for elem in self.__priorityOfEachLigand[keyRanks]:
				newRankRelation.append(rankTypesTransformMap[elem])
			auxDictCanonRanks[keyRanks] = newRankRelation
		self.__priorityOfEachLigand = dict(auxDictCanonRanks)


	def _definePrioritiesTransformMap(self):
		amountOfEachPriority = Counter(self.__finalPriorites)
		prioritiesTransformMap = {} # keys:old rank - value:new rank | always: [0 0 0 1 2 3 4]
		oldRankCounting = dict(amountOfEachPriority)
		for i in range(len(oldRankCounting)):
			maxKey = max(oldRankCounting.items(), key=operator.itemgetter(1))[0]
			listMax = [maxKey]
			for keysRank in oldRankCounting:
				if keysRank == maxKey:
					continue
				if oldRankCounting[keysRank] == oldRankCounting[maxKey]:
					listMax.append(keysRank)
			maxKey = min(listMax) #taking min of this list of equal values. This gives all hierarchy information
			prioritiesTransformMap[maxKey] = i
			del oldRankCounting[maxKey]
		
		return prioritiesTransformMap

	def _calculateNewMapBetweenAtomsAndTypes(self, rank):
		newMap = {}
		sortedRank = list(rank)
		sortedRank.sort()
		rankCounting = Counter(sortedRank)
		oldRankCounting = dict(rankCounting)
		for i in range(len(oldRankCounting)):
			maxKey = max(oldRankCounting.items(), key=operator.itemgetter(1))[0]
			newMap[maxKey] = i
			del oldRankCounting[maxKey]
				
		return newMap


if __name__ ==  "__main__":
	print("Module to handle molecular formulas")
	
	
	
