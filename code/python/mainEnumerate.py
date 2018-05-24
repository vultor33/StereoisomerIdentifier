from Enumeration import Enumeration
from DataAllFormulas import allFormList
from FormulaHandling import FormulaHandling
from collections import Counter
import operator

geoCode = 61
ncoord = 6

obj1 = Enumeration(geoCode)
enumeratedFormulas = []

countingFile = open("{}-counting.csv".format(geoCode),"w")
for i in range(len(allFormList[ncoord])):
	if i != 59:
		continue
	print('i:   ',i)
	if allFormList[ncoord][i][2] in enumeratedFormulas:
		continue

	counting = []
	obj1.makeEnumeration(allFormList[ncoord][i][2], allFormList[ncoord][i][3], allFormList[ncoord][i][4], counting)
	enumeratedFormulas.append(allFormList[ncoord][i][2])
	
	countingFile.write("{};".format(allFormList[ncoord][i][2]))
	for countI in counting:
		countingFile.write("{};".format(countI))
	countingFile.write("\n")