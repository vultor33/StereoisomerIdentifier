from Enumeration import Enumeration
from collections import Counter

alphabet = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','e','s','t','u','v','x','w','y','z'] # gerar alfabeto com o python

justLetters = lambda enterString : ''.join([i for i in enterString if i.isalpha()])


def calculateFormula(rank):
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

#last term must be the largest and first the lowest
def advanceFormula(formula):
	types = justLetters(formula)
	deltaIndex = alphabet.index(types[len(types)-1]) - alphabet.index(types[0]) + 1
	newFormula = ''
	for iFor in formula:
		if iFor.isdigit():
			newFormula += iFor
			continue
		newFormula += alphabet[deltaIndex + alphabet.index(iFor)]
	
	return newFormula
	

def formulasCanonRank(formula):
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
	
# if formula1 comes first: False
def formulasConflictTrade(formula1, formula2):
	canon1 = formulasCanonRank(formula1)
	canon2 = formulasCanonRank(formula2)

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

def generateMolecularFormula(rank, chelations):	
	dictChelations = {}
	dictNumbers = {}
	for chel in chelations:
		lisRankChel = []
		for iChel in chel:
			lisRankChel.append(rank[iChel])
		lisRankChel.sort()
		if lisRankChel[0] in dictChelations:
			if dictChelations[lisRankChel[0]] != lisRankChel:
				raise Exception('chelations not well defined')
			dictNumbers[lisRankChel[0]] += 1
		else:
			dictChelations[lisRankChel[0]] = lisRankChel
			dictNumbers[lisRankChel[0]] = 1

	allChelValues = []
	for chelList in dictChelations.values():
		allChelValues += chelList

	for iRank in rank:
		if iRank in allChelValues:
			continue
		if iRank in dictChelations:
			dictNumbers[iRank] += 1
		else:
			dictChelations[iRank] = [iRank]
			dictNumbers[iRank] = 1

	dictFormulas = {}
	for key in dictChelations:
		dictFormulas[key] = calculateFormula(dictChelations[key])


	trade = True
	while trade:
		trade = False
		for key1 in dictFormulas:
			for key2 in dictFormulas:
				if key1 == key2:
					continue
				elif dictFormulas[key1] == dictFormulas[key2]:
					print('conflit: ',key1, ' and ', key2)

					#solving conflict
					if dictNumbers[key1] > dictNumbers[key2]:
						dictFormulas[key2] = advanceFormula(dictFormulas[key2])
						trade = True
					elif dictNumbers[key1] < dictNumbers[key2]:
						dictFormulas[key1] = advanceFormula(dictFormulas[key1])
						trade = True
					else:
						if key1 < key2:
							dictFormulas[key2] = advanceFormula(dictFormulas[key2])
							trade = True
						elif key1 > key2:
							dictFormulas[key1] = advanceFormula(dictFormulas[key1])
							trade = True
						else:
							raise Exception('Error on molecular formula conflicts')
	
	dictPriorities = {}
	for key in dictFormulas:
		dictPriorities[key] = 0
	
	trade  = True
	while trade:
		trade = False
		for key1 in dictFormulas:
			for key2 in dictFormulas:
				if key1 == key2:
					continue
				if dictPriorities[key1] == dictPriorities[key2]:
					trade = True
					if formulasConflictTrade(dictFormulas[key1],dictFormulas[key2]):
						dictPriorities[key1] += 1
					else:
						dictPriorities[key2] += 1
	
	
					
						
				
	# USAR AS PRIORIDADES PARA CONSTRUIR A FORMULA E DEFINIR UMA NUMERACAO UNICA PARA TODOS
		
		
	print('PRIORITIES:   ',dictPriorities)
	print(dictFormulas)
	print(dictChelations)
	print(dictNumbers)
	return
		
	

	# contagem das cores
	rankCounting = Counter(rank)
	rankCountingElems = []
	rankCountingKeys = []
	for comp in rankCounting:
		rankCountingElems.append(rankCounting[comp])
		rankCountingKeys.append(comp)

	zipped = list(zip(rankCountingElems, rankCountingKeys))
	zipped.sort()
	rankCountingElems, rankCountingKeys = zip(*zipped)
	print('elems:  ', rankCountingElems)
	print('keys:  ',rankCountingKeys)
	
	ligandsAmount = list(rankCountingElems)
	ligandsAmount.reverse()
	ligandTypes = alphabet[:len(ligandsAmount)]
	molecularFormula = "M"
	i = 0		
	while i < len(ligandTypes):
		if ligandsAmount[i] == 1:
			molecularFormula += ligandTypes[i]
		else:
			molecularFormula += ligandTypes[i] + str(ligandsAmount[i])
		i+=1
		
	# For C++
	newRank = []
	for iRank in rank:
		newRank.append(rankCountingKeys.index(iRank))
	for i in range(len(newRank)):
		rank[i] = newRank[i]
	print(rank)
	return molecularFormula


		
		
		
		
		
# 	

print('conf:  ',formulasConflictTrade('a3b', 'a3b2'))

print('f: ',calculateFormula([7,5,7,9,9,7,8,11,12,13,7,3,3,2,2,15]))

print('advance formula:   ',advanceFormula('ac2d2fg5'))

print('formula:  ',generateMolecularFormula([7,5,7,9,9,7,8,11,12,13,7,3,3,2,2,15,2,2,15], [[0,2,3],[4,5,10],[1,8,9],[11,12],[13,14,15],[16,17,18]]))

exit()


substratalScaffoldTemp = [[1,2,3,4],[1,3,2,4],[1,2,4,3]]
allRotations = [[1,2,3,0],[2,3,0,1],[3,0,1,2],[0,3,2,1],[2,1,0,3],[1,0,3,2],[3,2,1,0]]
colorFirst = [1,1,2,2]
chelateFirst = []
stereoRemoved = []
rcw = []

obj1 = Enumeration()
obj1.readRotations(41)
obj1.readScaffold(41)
exit()
obj1.setSubstratal(substratalScaffoldTemp)
del allRotations
del substratalScaffoldTemp

obj1.makeEnumeration(colorFirst, chelateFirst, rcw, stereoRemoved)

print('rcw:  ', rcw)
print('removed: ', stereoRemoved)
