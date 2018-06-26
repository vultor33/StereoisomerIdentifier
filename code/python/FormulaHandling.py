from collections import Counter
import itertools
import math
import operator

alphabet = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','e','s','t','u','v','x','w','y','z'] # gerar alfabeto com o python

justLetters = lambda enterString : ''.join([i for i in enterString if i.isalpha()])

class FormulaHandling:
	"""Class to generate and transform molecular formulas"""
	
	def __init__(self):

		#Ex: received: rank: [10,9,14,45,45]. chelations: [[1,2],[3,4]]
		self.__finalFormula = ''                            # Molecular formula. Ex: Ma(A2)(AB)
		self.__referenceLineVector = []                     # Line of priority types. Ex: [0,1,1,2,3]
		self.__canonChelation = []                          # Chelation related to the referenceLineVector. Ex: [[1,2],[3,4]]
		

		# These are relations of canon types and received types
		self.__newRankTypesTransformMap = {}                # ranks were transformed. this dict relates old types to new types.
		self.__dictCanonRanks = {}                          # Rank of each type and chelate. Ex: {10:[0],9:[1,2]}
		self.__dictFormulas = {}                            # Formula of each type. Ex: {9:'ab',10:'a'}
		self.__dictChelations = {}                          # Keys: first type of the chelating group. values: all other types ordered.  Ex: {1:[1,2]}
		self.__dictNumbers = {}                             # Number of each type of ligands. Polidentates are counted as one and it's first type is used. Ex:{0:1,1:1}
		
		
	def printInfo(self):
		print('Canon chelation:  ',self.__canonChelation)
		print('formula:  ',self.__finalFormula)
		print('referencelinevector:  ',self.__referenceLineVector)
		print('Canon relation:  ',self.__dictCanonRanks)		
		print('Formulas:  ',self.__dictFormulas)
		print('Formula numbers:  ',self.__dictNumbers)
		print('Chelations from types:  ',self.__dictChelations)

	def getFormula(self):
		return self.__finalFormula
		
	def getReferenceLine(self):
		return self.__referenceLineVector
	
	def getCanonChelation(self):
		return self.__canonChelation
		
	def getRankTransformMap(self):
		return self.__newRankTypesTransformMap	

	def calculateAllChelateCombinations(self,size):
		allChelComb = []
		chelPossible = list(range(2,size + 1))
		chelNumber = 1 
		while chelNumber <= int(math.floor(size/2)):
			allChel = list(itertools.combinations_with_replacement(chelPossible, chelNumber))
			for chelI in allChel:
				sum = 0
				for iSum in chelI:
					sum += iSum
				if sum <= size:
					allChelComb.append(chelI)
			chelNumber += 1
		return allChelComb

	def calculateAllChelatesReaches(self, size, chelSize):
		chelReach = list(range(0,size))
		return list(itertools.combinations(chelReach, chelSize))

	#use chelations only if there is repeating ligands
	def generateEnumerationFormula(self, rank, chelations):
		oldMap = {}
		for chel in chelations:
			for chelI in chel:
				oldMap[chelI] = rank[chelI]

		duplicate = False
		for i in range(len(chelations)-1):
			for j in range(len(chelations)):
				if i == j:
					continue
				for chelI in chelations[i]:
					for chelJ in chelations[j]:
						if oldMap[chelI] == oldMap[chelJ]:
							duplicate = True
		if duplicate:
			self.generateMolecularFormula(rank, chelations)
		else:
			self.generateMolecularFormula(rank, [])

	
	def generateMolecularFormula(self, rank, chelations):
		if self._chelationsConflict(chelations):
			raise Exception('chelations not well defined')

		#count chelates and define self.__dictChelations (see above)	
		for chel in chelations:
			lisRankChel = []
			for iChel in chel:
				lisRankChel.append(rank[iChel])
			lisRankChel.sort()
			if lisRankChel[0] in self.__dictChelations:
				if self.__dictChelations[lisRankChel[0]] != lisRankChel: #different chelations and same types arent allowed
					raise Exception('chelations not well defined')
				self.__dictNumbers[lisRankChel[0]] += 1
			else:
				self.__dictChelations[lisRankChel[0]] = lisRankChel
				self.__dictNumbers[lisRankChel[0]] = 1

		allChelPositios = []
		for chelList in chelations:
			for chelI in chelList:
				allChelPositios += [chelI]

		# count monodentates
		for iAux in range(len(rank)):
			if iAux in allChelPositios:
				continue
			if rank[iAux] in self.__dictChelations:
				if self.__dictChelations[rank[iAux]] != [rank[iAux]]:
					raise Exception('chelations not well defined')
				self.__dictNumbers[rank[iAux]] += 1
			else:
				self.__dictChelations[rank[iAux]] = [rank[iAux]]
				self.__dictNumbers[rank[iAux]] = 1

		#different chelations and same types arent allowed
		for keyChel in self.__dictChelations:
			for chelRank in self.__dictChelations[keyChel]:
				for keyChel2 in self.__dictChelations:
					if keyChel == keyChel2:
						continue
					for chelRank2 in self.__dictChelations[keyChel2]:
						if chelRank == chelRank2:
							raise Exception('chelations not well defined')
								
						
			
				


		for chel in chelations:
			lisRankChel = []
			for iChel in chel:
				lisRankChel.append(rank[iChel])
			lisRankChel.sort()
			if lisRankChel[0] in self.__dictChelations:
				if self.__dictChelations[lisRankChel[0]] != lisRankChel: #different chelations and same type arent allowed
					raise Exception('chelations not well defined')



		for key in self.__dictChelations:
			self.__dictFormulas[key] = self._calculateFormula(self.__dictChelations[key])


		trade = True
		while trade:
			trade = False
			for key1 in self.__dictFormulas:
				for key2 in self.__dictFormulas:
					if key1 == key2:
						continue
					elif self.__dictFormulas[key1] == self.__dictFormulas[key2]:
						#solving conflict
						if self.__dictNumbers[key1] > self.__dictNumbers[key2]:
							self.__dictFormulas[key2] = self._advanceFormula(self.__dictFormulas[key2])
							trade = True
						elif self.__dictNumbers[key1] < self.__dictNumbers[key2]:
							self.__dictFormulas[key1] = self._advanceFormula(self.__dictFormulas[key1])
							trade = True
						else:
							if key1 < key2:
								self.__dictFormulas[key2] = self._advanceFormula(self.__dictFormulas[key2])
								trade = True
							elif key1 > key2:
								self.__dictFormulas[key1] = self._advanceFormula(self.__dictFormulas[key1])
								trade = True
							else:
								raise Exception('Error on molecular formula conflicts')
	

		dictPriorities = {}  #Priorities of each key, add 1 if it is bidentate and etc. This are used to define self.__dictCanonRanks
		for key in self.__dictFormulas:
			dictPriorities[key] = 0
	
		# bubble method to define priorities. All starts at zero, changes are decided by self._formulasConflictTrade
		trade  = True
		while trade:
			trade = False
			for key1 in self.__dictFormulas:
				for key2 in self.__dictFormulas:
					if key1 == key2:
						continue
					if dictPriorities[key1] == dictPriorities[key2]:
						trade = True
						if self._formulasConflictTrade(self.__dictFormulas[key1],self.__dictFormulas[key2]):
							dictPriorities[key1] += 1
						else:
							dictPriorities[key2] += 1

		i = 0
		j = 0
		k = 0
		while i < len(self.__dictFormulas.keys()):
			keyI = self._dictFindKeyByValue(dictPriorities,i)
			canon = self._formulasCanonRank(self.__dictFormulas[keyI])
			tempK = canon[len(canon) - 1]
			for iCan in range(len(canon)):
				canon[iCan] += k
			self.__dictCanonRanks[keyI] = canon
			for iTemp in range(self.__dictNumbers[keyI]):
				chelTemp = []
				for iCan in canon:
					self.__referenceLineVector += [iCan]
					chelTemp.append(j)
					j+=1
				if len(chelTemp) > 1:
					self.__canonChelation.append(chelTemp)
			
			
			if len(self.__dictFormulas[keyI]) > 1:
				self.__finalFormula += '(' + str(self.__dictFormulas[keyI]).upper() + ')'
			else:
				self.__finalFormula += str(self.__dictFormulas[keyI])
			if self.__dictNumbers[keyI] > 1:
				self.__finalFormula += str(self.__dictNumbers[keyI])
	
			k += tempK + 1
			i+=1


		# Until here we have an hierarchy on priorities - mono, bidentate, tridentate and etc. This function changes to: numbers come first. All previous hierarchy is keeped.
		self._redefineChelatesTypes(self.__referenceLineVector,self.__canonChelation)


	def _redefineChelatesTypes(self, rank, chelatesList):
		# Types are redefined. Bidentates stick to the types

		#if chelatesList == []:
		#	return
		
		#rank type each received chelate is pointing
		oldMap = {}
		for chel in chelatesList:
			for chelI in chel:
				oldMap[chelI] = rank[chelI]
				

		rankTypesTransformMap = {} # keys:old rank - value:new rank | always: [0 0 0 1 2 3 4]
		rankCounting = Counter(rank)
		oldRankCounting = dict(rankCounting)
		for i in range(len(oldRankCounting)):
			maxKey = max(oldRankCounting.items(), key=operator.itemgetter(1))[0]
			listMax = [maxKey]
			for keysRank in oldRankCounting:
				if keysRank == maxKey:
					continue
				if oldRankCounting[keysRank] == oldRankCounting[maxKey]:
					listMax.append(keysRank)
			maxKey = min(listMax) #taking min of this list of equal values. This gives all hierarchy information
			rankTypesTransformMap[maxKey] = i
			del oldRankCounting[maxKey]

		auxDictCanonRanks = {}
		for keyRanks in self.__dictCanonRanks:
			newRankRelation = []
			for elem in self.__dictCanonRanks[keyRanks]:
				newRankRelation.append(rankTypesTransformMap[elem])
			auxDictCanonRanks[keyRanks] = newRankRelation
			
		self.__dictCanonRanks = dict(auxDictCanonRanks)

		#Initial type: self.__dictChelations[keyI] (but I have to order here
		#Final type:  self.__dictCanonRanks[keyI][i]. 
		for keyI in self.__dictChelations:
			#dictChelations comes like this: [3, 9, 9, 10, 11]. I need repeating priorities to come first.
			mapBetweenChelatesAndTypes = self._calculateNewMapBetweenAtomsAndTypes(self.__dictChelations[keyI])
			rankInTheInteriorOfTheChelate = list(self.__dictChelations[keyI])
			for iRank in range(len(rankInTheInteriorOfTheChelate)):
				rankInTheInteriorOfTheChelate[iRank] = mapBetweenChelatesAndTypes[rankInTheInteriorOfTheChelate[iRank]]
			rankInTheInteriorOfTheChelate.sort()
			for i in range(len(self.__dictChelations[keyI])):
				finalValue = self.__dictCanonRanks[keyI][i]
				internal = rankInTheInteriorOfTheChelate[i]
				internalToCipRank = self._dictFindKeyByValue(mapBetweenChelatesAndTypes,internal)

				if internalToCipRank in self.__newRankTypesTransformMap:
					if self.__newRankTypesTransformMap[internalToCipRank] != finalValue:
						raise Exception('Error on _redefineChelatesTypes: new ranks couldnt be set')
	
				self.__newRankTypesTransformMap[internalToCipRank] = finalValue

		auxTypes = []
		for key in rankCounting:
			for i in range(rankCounting[key]):
				auxTypes.append(rankTypesTransformMap[key])

		#rever isso aqui tambem
		auxTypes.sort()
		for i in range(len(auxTypes)):
			rank[i] = auxTypes[i]
		for i in range(len(chelatesList)):
			for j in range(len(chelatesList[i])):
				chelatesList[i][j] = auxTypes.index(rankTypesTransformMap[oldMap[chelatesList[i][j]]])
				auxTypes[chelatesList[i][j]] = -1


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

			
	def _chelationsConflict(self, chelations):
		# two chelates pointing to the same ligand
		i = 0
		j = 0
		while i < len(chelations) - 1:
			j = i + 1
			while j < len(chelations):
				for chelI in chelations[i]:
					for chelJ in chelations[j]:
						if chelI == chelJ:
							return True
				j+=1
			i+=1
		
		return False


	def _dictFindKeyByValue(self, dict, search_value):
		for key, value in dict.items():
		    if value == search_value:
        		return key

	def _calculateFormula(self, rank):
		rankCounting = Counter(rank)
		rankCountingElems = []
		for comp in rankCounting:
			rankCountingElems.append(rankCounting[comp])

		rankCountingElems.sort()
		ligandsAmount = list(rankCountingElems)
		ligandsAmount.reverse()
		ligandTypes = alphabet[:len(ligandsAmount)]
		molecularFormula = ""
		i = 0		
		while i < len(ligandTypes):
			if ligandsAmount[i] == 1:
				molecularFormula += ligandTypes[i]
			else:
				molecularFormula += ligandTypes[i] + str(ligandsAmount[i])
			i+=1
		return molecularFormula

	#from formula: last term must be the largest and first the lowest
	def _advanceFormula(self,formula):
		types = justLetters(formula)
		deltaIndex = alphabet.index(types[len(types)-1]) - alphabet.index(types[0]) + 1
		newFormula = ''
		for iFor in formula:
			if iFor.isdigit():
				newFormula += iFor
				continue
			newFormula += alphabet[deltaIndex + alphabet.index(iFor)]
	
		return newFormula
	
	def _formulasCanonRank(self, formula):
		i = 0
		k = 0
		rankCanon = []
		while i < len(formula):
			if formula[i].isdigit():
				for iR in range(int(float(formula[i]))-1):
					rankCanon.append(k-1)
			else:
				rankCanon.append(k)
				k+=1
			i+=1
		return rankCanon
	
	# Don't trade (False) if formula1 should come first
	def _formulasConflictTrade(self, formula1, formula2):
		canon1 = self._formulasCanonRank(formula1)
		canon2 = self._formulasCanonRank(formula2)
	
		if len(canon1) != len(canon2):
			return len(canon1) > len(canon2)
		else:
			i = 0
			while i < len(canon1):
				if canon1[i] != canon2[i]:
					return canon1[i] > canon2[i]
				i+=1
			
			if alphabet.index(formula1[0]) != alphabet.index(formula2[0]):
				return alphabet.index(formula1[0]) > alphabet.index(formula2[0])
			else:
				raise Exception('Error on formulasConflictTrade - check input for errors')

		
	

		
		



























##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################




class GenerateAllFormulas:
	"""Class to enumerate all formulas with all possible chelations"""
	# How to use
	#objG = GenerateAllFormulas()
	#objG.generateAllFormulas(8)
	#objG.printFile()

	def __init__(self):
		self.__allFormulas = []
		self.__allReferences = []
		self.__allChelations = []

	def generateAllFormulas(self, size):
		allComb = list(itertools.combinations_with_replacement(list(range(0,size)), size))
		for perm in allComb: # monodentates
			objP = FormulaHandling()
			objP.generateMolecularFormula(perm, [])
			if not objP.getFormula() in self.__allFormulas:
				self.__allFormulas.append(objP.getFormula())
				self.__allReferences.append(objP.getReferenceLine())
				self.__allChelations.append(['None'])
				
		allChelComb = self._calculateAllChelateCombinations(size)
		monodentateReferences = list(self.__allReferences)
		for ref in monodentateReferences:  #use monodentates as base to chelations definitions
			for chelCom in allChelComb:
				chelList = []
				for chelComI in chelCom:
					chelReach = list(range(0,size))
					chelList.append(list(itertools.combinations(chelReach, chelComI)))

				allChelationsToFormulas = list(itertools.product(*chelList))
				for chelations in allChelationsToFormulas:
					try:
						objC = FormulaHandling()
						objC.generateMolecularFormula(ref, chelations)
						if not objC.getFormula() in self.__allFormulas:
							self.__allFormulas.append(objC.getFormula())
							self.__allReferences.append(objC.getReferenceLine())
							self.__allChelations.append(objC.getCanonChelation())
		
					except Exception as e:
						if str(e) == "chelations not well defined":
							pass
						else:
							print(" Error:  ",str(e))

	def redefineFormulas(self):
		#Modify the function: self.printFile to its final form
		allFormFile = open("allFormulasRedefined.py", "w")
		allFormFile.write("allFormList = []\n")
	
		ncoord = 1
		allFormFile.write("aux1FormList = []\n")
		for i in range(len(allFormList1)):
			objFormula = FormulaHandling()
			if allFormList1[i][2] == ['None']:
				objFormula.generateMolecularFormula(allFormList1[i][1],[])
			else:
				objFormula.generateMolecularFormula(allFormList1[i][1],allFormList1[i][2])
			allFormList1[i][1] = objFormula.getReferenceLine()
			allFormList1[i][2] = objFormula.getCanonChelation()
			objFormula2 = FormulaHandling()
			if allFormList1[i][2] == ['None']:
				objFormula2.generateEnumerationFormula(allFormList1[i][1],[])
			else:
				objFormula2.generateEnumerationFormula(allFormList1[i][1],allFormList1[i][2])
		
			allFormFile.write("aux1FormList.append([{},\'{}\',\'{}\',{},{}])\n".format(
			ncoord,
			allFormList1[i][0],
			objFormula2.getFormula(),
			objFormula.getReferenceLine(),
			objFormula.getCanonChelation()))
		allFormFile.write("allFormList.append(aux1FormList)\n")
	
		ncoord = 2
		allFormFile.write("aux2FormList = []\n")
		for i in range(len(allFormList2)):
			objFormula = FormulaHandling()
			if allFormList2[i][2] == ['None']:
				objFormula.generateMolecularFormula(allFormList2[i][1],[])
			else:
				objFormula.generateMolecularFormula(allFormList2[i][1],allFormList2[i][2])
			allFormList2[i][1] = objFormula.getReferenceLine()
			allFormList2[i][2] = objFormula.getCanonChelation()
			objFormula2 = FormulaHandling()
			if allFormList2[i][2] == ['None']:
				objFormula2.generateEnumerationFormula(allFormList2[i][1],[])
			else:
				objFormula2.generateEnumerationFormula(allFormList2[i][1],allFormList2[i][2])
		
			allFormFile.write("aux2FormList.append([{},\'{}\',\'{}\',{},{}])\n".format(
			ncoord,
			allFormList2[i][0],
			objFormula2.getFormula(),
			objFormula.getReferenceLine(),
			objFormula.getCanonChelation()))
		allFormFile.write("allFormList.append(aux2FormList)\n")
	
		ncoord = 3
		allFormFile.write("aux3FormList = []\n")
		for i in range(len(allFormList3)):
			objFormula = FormulaHandling()
			if allFormList3[i][2] == ['None']:
				objFormula.generateMolecularFormula(allFormList3[i][1],[])
			else:
				objFormula.generateMolecularFormula(allFormList3[i][1],allFormList3[i][2])
			allFormList3[i][1] = objFormula.getReferenceLine()
			allFormList3[i][2] = objFormula.getCanonChelation()
			objFormula2 = FormulaHandling()
			if allFormList3[i][2] == ['None']:
				objFormula2.generateEnumerationFormula(allFormList3[i][1],[])
			else:
				objFormula2.generateEnumerationFormula(allFormList3[i][1],allFormList3[i][2])
	
			allFormFile.write("aux3FormList.append([{},\'{}\',\'{}\',{},{}])\n".format(
			ncoord,
			allFormList3[i][0],
			objFormula2.getFormula(),
			objFormula.getReferenceLine(),
			objFormula.getCanonChelation()))
		allFormFile.write("allFormList.append(aux3FormList)\n")
	
		ncoord = 4
		allFormFile.write("aux4FormList = []\n")
		for i in range(len(allFormList4)):
			objFormula = FormulaHandling()
			if allFormList4[i][2] == ['None']:
				objFormula.generateMolecularFormula(allFormList4[i][1],[])
			else:
				objFormula.generateMolecularFormula(allFormList4[i][1],allFormList4[i][2])
			allFormList4[i][1] = objFormula.getReferenceLine()
			allFormList4[i][2] = objFormula.getCanonChelation()
			objFormula2 = FormulaHandling()
			if allFormList4[i][2] == ['None']:
				objFormula2.generateEnumerationFormula(allFormList4[i][1],[])
			else:
				objFormula2.generateEnumerationFormula(allFormList4[i][1],allFormList4[i][2])
		
			allFormFile.write("aux4FormList.append([{},\'{}\',\'{}\',{},{}])\n".format(
			ncoord,
			allFormList4[i][0],
			objFormula2.getFormula(),
			objFormula.getReferenceLine(),
			objFormula.getCanonChelation()))
		allFormFile.write("allFormList.append(aux4FormList)\n")
	
		ncoord = 5
		allFormFile.write("aux5FormList = []\n")
		for i in range(len(allFormList5)):
			objFormula = FormulaHandling()
			if allFormList5[i][2] == ['None']:
				objFormula.generateMolecularFormula(allFormList5[i][1],[])
			else:
				objFormula.generateMolecularFormula(allFormList5[i][1],allFormList5[i][2])
			allFormList5[i][1] = objFormula.getReferenceLine()
			allFormList5[i][2] = objFormula.getCanonChelation()
			objFormula2 = FormulaHandling()
			if allFormList5[i][2] == ['None']:
				objFormula2.generateEnumerationFormula(allFormList5[i][1],[])
			else:
				objFormula2.generateEnumerationFormula(allFormList5[i][1],allFormList5[i][2])
		
			allFormFile.write("aux5FormList.append([{},\'{}\',\'{}\',{},{}])\n".format(
			ncoord,
			allFormList5[i][0],
			objFormula2.getFormula(),
			objFormula.getReferenceLine(),
			objFormula.getCanonChelation()))
		allFormFile.write("allFormList.append(aux5FormList)\n")
	
		ncoord = 6
		allFormFile.write("aux6FormList = []\n")
		for i in range(len(allFormList6)):
			objFormula = FormulaHandling()
			if allFormList6[i][2] == ['None']:
				objFormula.generateMolecularFormula(allFormList6[i][1],[])
			else:
				objFormula.generateMolecularFormula(allFormList6[i][1],allFormList6[i][2])
			allFormList6[i][1] = objFormula.getReferenceLine()
			allFormList6[i][2] = objFormula.getCanonChelation()
			objFormula2 = FormulaHandling()
			if allFormList6[i][2] == ['None']:
				objFormula2.generateEnumerationFormula(allFormList6[i][1],[])
			else:
				objFormula2.generateEnumerationFormula(allFormList6[i][1],allFormList6[i][2])
		
			allFormFile.write("aux6FormList.append([{},\'{}\',\'{}\',{},{}])\n".format(
			ncoord,
			allFormList6[i][0],
			objFormula2.getFormula(),
			objFormula.getReferenceLine(),
			objFormula.getCanonChelation()))
		allFormFile.write("allFormList.append(aux6FormList)\n")
	
		ncoord = 7
		allFormFile.write("aux7FormList = []\n")
		for i in range(len(allFormList7)):
			objFormula = FormulaHandling()
			if allFormList7[i][2] == ['None']:
				objFormula.generateMolecularFormula(allFormList7[i][1],[])
			else:
				objFormula.generateMolecularFormula(allFormList7[i][1],allFormList7[i][2])
			allFormList7[i][1] = objFormula.getReferenceLine()
			allFormList7[i][2] = objFormula.getCanonChelation()
			objFormula2 = FormulaHandling()
			if allFormList7[i][2] == ['None']:
				objFormula2.generateEnumerationFormula(allFormList7[i][1],[])
			else:
				objFormula2.generateEnumerationFormula(allFormList7[i][1],allFormList7[i][2])
		
			allFormFile.write("aux7FormList.append([{},\'{}\',\'{}\',{},{}])\n".format(
			ncoord,
			allFormList7[i][0],
			objFormula2.getFormula(),
			objFormula.getReferenceLine(),
			objFormula.getCanonChelation()))
		allFormFile.write("allFormList.append(aux7FormList)\n")
	
		ncoord = 8
		allFormFile.write("aux8FormList = []\n")
		for i in range(len(allFormList8)):
			objFormula = FormulaHandling()
			if allFormList8[i][2] == ['None']:
				objFormula.generateMolecularFormula(allFormList8[i][1],[])
			else:
				objFormula.generateMolecularFormula(allFormList8[i][1],allFormList8[i][2])
			allFormList8[i][1] = objFormula.getReferenceLine()
			allFormList8[i][2] = objFormula.getCanonChelation()
			objFormula2 = FormulaHandling()
			if allFormList8[i][2] == ['None']:
				objFormula2.generateEnumerationFormula(allFormList8[i][1],[])
			else:
				objFormula2.generateEnumerationFormula(allFormList8[i][1],allFormList8[i][2])
		
			allFormFile.write("aux8FormList.append([{},\'{}\',\'{}\',{},{}])\n".format(
			ncoord,
			allFormList8[i][0],
			objFormula2.getFormula(),
			objFormula.getReferenceLine(),
			objFormula.getCanonChelation()))
		allFormFile.write("allFormList.append(aux8FormList)\n")
	
		allFormFile.close()
		
	

	def printFile(self):
		allFormFile = open("allFormulas.py", "w")
		allFormFile.write("allFormList = []\n")
		for i in range(len(self.__allFormulas)):
			allFormFile.write("allFormList.append([\'{}\',{},{}])\n".format(
			self.__allFormulas[i],
			self.__allReferences[i],
			self.__allChelations[i]))
			
		allFormFile.close()			

	def _calculateAllChelateCombinations(self,size):
		allChelComb = []
		chelNumber = 1 
		while chelNumber <= int(math.floor(size/2)):
			allChel = list(itertools.combinations_with_replacement(list(range(2,size + 1)), chelNumber))
			for chelI in allChel:
				sum = 0
				for iSum in chelI:
					sum += iSum
				if sum <= size:
					allChelComb.append(chelI)
			chelNumber += 1
		return allChelComb



#########################################################################
#########################################################################
# TESTE ALL FORM LIST
#for x in allFormList:
#	len(x)
#	for y in x:
#		if y == 0:
#			continue
#		if len(y) > 3:
#			obj1 = FormulaHandling()
#			obj1.generateMolecularFormula(y[3], y[4])
#			print(y)
#			if y[3] != obj1.getReferenceLine():
#				print("Problem on monodentates definitions")
#				exit()
#			if y[4] != obj1.getCanonChelation():
#				print("Problem on chelate definitions")
#				exit()
#########################################################################
#########################################################################


if __name__ ==  "__main__":
	print("Module to handle molecular formulas")
	
	
	
