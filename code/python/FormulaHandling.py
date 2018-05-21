from collections import Counter
import itertools
import math


alphabet = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','e','s','t','u','v','x','w','y','z'] # gerar alfabeto com o python

justLetters = lambda enterString : ''.join([i for i in enterString if i.isalpha()])

class FormulaHandling:
	"""Class to generate and transform molecular formulas"""
	
	def __init__(self):
		self.__canonChelation = []
		self.__finalFormula = ''
		self.__referenceLineVector = []
		self.__dictCanonRanks = {}
		self.__dictPriorities = {} 
		self.__dictFormulas = {}
		self.__dictChelations = {}
		self.__dictNumbers = {}
		
	def printInfo(self):
		print('Canon chelation:  ',self.__canonChelation)
		print('formula:  ',self.__finalFormula)
		print('referencelinevector:  ',self.__referenceLineVector)
		print('Canon relation:  ',self.__dictCanonRanks)		
		print('Priorities:  ',self.__dictPriorities)
		print('Formulas:  ',self.__dictFormulas)
		print('Formula numbers:  ',self.__dictNumbers)
		print('Chelations from types:  ',self.__dictChelations)

	def getFormula(self):
		return self.__finalFormula
		
	def getReferenceLine(self):
		return self.__referenceLineVector
	
	def getCanonChelation(self):
		return self.__canonChelation

	def generateMolecularFormula(self, rank, chelations):	
		if self._chelationsConflict(chelations):
			raise Exception('chelations not well defined')
	
		for chel in chelations:
			lisRankChel = []
			for iChel in chel:
				lisRankChel.append(rank[iChel])
			lisRankChel.sort()
			if lisRankChel[0] in self.__dictChelations:
				if self.__dictChelations[lisRankChel[0]] != lisRankChel:
					raise Exception('chelations not well defined')
				self.__dictNumbers[lisRankChel[0]] += 1
			else:
				self.__dictChelations[lisRankChel[0]] = lisRankChel
				self.__dictNumbers[lisRankChel[0]] = 1

		allChelPositios = []
		for chelList in chelations:
			for chelI in chelList:
				allChelPositios += [chelI]
	
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
	
		for key in self.__dictFormulas:
			self.__dictPriorities[key] = 0
	
		trade  = True
		while trade:
			trade = False
			for key1 in self.__dictFormulas:
				for key2 in self.__dictFormulas:
					if key1 == key2:
						continue
					if self.__dictPriorities[key1] == self.__dictPriorities[key2]:
						trade = True
						if self._formulasConflictTrade(self.__dictFormulas[key1],self.__dictFormulas[key2]):
							self.__dictPriorities[key1] += 1
						else:
							self.__dictPriorities[key2] += 1

		i = 0
		j = 0
		k = 0
		while i < len(self.__dictFormulas.keys()):
			keyI = self._dictFindKeyByValue(self.__dictPriorities,i)
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


	def calculateAllChelateCombinations(self,size):
		allChelComb = []
		chelPossible = []
		i = 0
		while i < size - 1:
			chelPossible.append(i + 2)
			i+=1
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

	def calculateAllChelatesReaches(self,size, chelSize):
		chelReach = []
		i = 0
		while i < size:
			chelReach.append(i)
			i+=1
		return list(itertools.combinations(chelReach, chelSize))
		
	def _chelationsConflict(self, chelations):
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

		
	

		
		


if __name__ ==  "__main__":
	print("Module to handle molecular formulas")