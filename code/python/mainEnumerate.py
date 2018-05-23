from Enumeration import Enumeration
from allFormulasRedefined import allFormList
from FormulaHandling import FormulaHandling
from collections import Counter
import operator

geoCode = 60
ncoord = 6 
ncoord -= 1 #list starts at 0

obj1 = Enumeration(geoCode)
enumeratedFormulas = []

for i in range(len(allFormList[ncoord])):
#	if i != 9:
#		continue

	if allFormList[ncoord][i][2] in enumeratedFormulas:
		continue

	obj1.makeEnumeration(allFormList[ncoord][i][2], allFormList[ncoord][i][3], allFormList[ncoord][i][4])
	enumeratedFormulas.append(allFormList[ncoord][i][2])