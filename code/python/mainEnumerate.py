from Enumeration import Enumeration
from FormulaHandling import FormulaHandling
from collections import Counter
import itertools

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

print(allFormulas)
print(allReference)
allBi = list(itertools.combinations([0,1,2,3], 3))
print(allBi)
print('testing ref')
for ref in allReference:
	for bi in allBi:
		objP = FormulaHandling()
		objP.generateMolecularFormula(ref, [bi])
		objP.printInfo()
		print(objP.getFormula())


#depois aplica os bidentados nas formulas daqui


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
