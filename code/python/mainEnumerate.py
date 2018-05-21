from Enumeration import Enumeration
from FormulaHandling import FormulaHandling
from FormulaHandling import GenerateAllFormulas
from collections import Counter
import itertools
import math
import csv
import json
from DataAllFormulas import allFormList4

#for form in allFormList8:
#	print(form)

#substratalScaffoldTemp = [[1,2,3,4],[1,3,2,4],[1,2,4,3]]
#allRotations = [[1,2,3,0],[2,3,0,1],[3,0,1,2],[0,3,2,1],[2,1,0,3],[1,0,3,2],[3,2,1,0]]
#colorFirst = [1,1,2,2]
#chelateFirst = []

obj1 = Enumeration()
obj1.readRotations(41)
obj1.readScaffold(41)
#exit()
#obj1.setSubstratal(substratalScaffoldTemp)
#del allRotations
#del substratalScaffoldTemp

for i in range(len(allFormList4)):
#	if i != 8:
#		continue
	stereoRemoved = []
	rcw = []
	if allFormList4[i][2] == ['None']:
		obj1.makeEnumeration(allFormList4[i][1], [], rcw, stereoRemoved)
	else:
		obj1.makeEnumeration(allFormList4[i][1], allFormList4[i][2], rcw, stereoRemoved)
	
	print(allFormList4[i][0],'  Total:  ',3 - len(stereoRemoved))
#	exit()

	
#print('rcw:  ', rcw)
#print('removed: ', stereoRemoved)
