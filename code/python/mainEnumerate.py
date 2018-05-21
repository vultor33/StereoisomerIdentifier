from Enumeration import Enumeration
from FormulaHandling import FormulaHandling
from collections import Counter
import itertools
import math

objF = FormulaHandling()
		
objF.generateMolecularFormula([7,5,7,9,9,7,8,11,12,13,7,3,3,2,2,15,2,2,15], [[0,2,3],[4,5,10],[1,8,9],[11,12],[13,14,15],[16,17,18]])

allPerm = list(itertools.product([1, 2, 3, 4], repeat = 4))

#allPerm = list(itertools.combinations('0123', 2))

#print(allPerm)

allFormulas = []
allReference = []

for perm in allPerm:
	objP = FormulaHandling()
	objP.generateMolecularFormula(perm, [])
	if not objP.getFormula() in allFormulas:
		allFormulas.append(objP.getFormula())
		allReference.append(objP.getReferenceLine())

# Construir uma forma de soltar todos os quelados
# nao permitir que os quelados fiquem encavalados dois atomos apontando para o mesmo cara.
# 3 com problema - resolver

#print(allFormulas)
#print(allReference)
#print(allBi)
#allBi = list(itertools.combinations([2,3,4], 2))
#allBi = list(itertools.combinations_with_replacement([2,3,4,5,6,7], 3))
#print(allBi)
		
#print(objP.calculateAllChelateCombinations(8))
#print(objP.calculateAllChelatesReaches(4, 4))
#exit()

objA = FormulaHandling()
n = 4
allChelComb = objA.calculateAllChelateCombinations(n)
allFormulas = []

for ref in allReference:
	for chelCom in allChelComb:
		chelList = []
		for chelComI in chelCom:
			chelList.append(objA.calculateAllChelatesReaches(n, chelComI))
		allChelationsToFormulas = list(itertools.product(*chelList))
		for chelations in allChelationsToFormulas:
			try:
				objC = FormulaHandling()
				objC.generateMolecularFormula(ref, chelations)
				if not objC.getFormula() in allFormulas:
					allFormulas.append(objC.getFormula())

			except Exception as e:
				if str(e) == "chelations not well defined":
					print(end='')
				else:
					print(" Error:  ",str(e))

print(allFormulas)

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
