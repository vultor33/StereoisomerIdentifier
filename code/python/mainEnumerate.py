from Enumeration import Enumeration
from DataAllFormulas import allFormList8
from FormulaHandling import FormulaHandling
from collections import Counter


allFormList = allFormList8

objFormula = FormulaHandling()

i = 288
print(allFormList8[i][0])
print(allFormList8[i][1])
print(allFormList8[i][2])

#print(objFormula.generateMolecularFormula(allFormList8[i][1],allFormList8[i][2]))
print(objFormula.generateMolecularFormula(allFormList8[i][1],[]))

# na verdade os bidentados só precisam da quantidade de atomos que eles tinham antes
# ja que tudo pode variar livremente.

print(objFormula.getFormula())
print(objFormula.getReferenceLine())
print(objFormula.getCanonChelation())


rank = allFormList8[i][1]
rankCounting = Counter(rank)
rankCountingElems = []
rankCountingKeys = []
for comp in rankCounting:
	rankCountingElems.append(rankCounting[comp])
	rankCountingKeys.append(comp)

print('keys:   ',rankCountingKeys)
print('elems:  ',rankCountingElems)
print('list r: ',list(range(len(rankCountingKeys))))


exit()



geoCode = 60

obj1 = Enumeration(geoCode)

for i in range(len(allFormList)):
#	if i != 9:
#		continue
	if allFormList[i][2] == ['None']:
		obj1.makeEnumeration(allFormList[i][0], allFormList[i][1], [])
	else:
		obj1.makeEnumeration(allFormList[i][0], allFormList[i][1], allFormList[i][2])
