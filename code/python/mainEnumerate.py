from Enumeration import Enumeration
from FormulaHandling import FormulaHandling
from FormulaHandling import GenerateAllFormulas
from collections import Counter
import itertools
import math
import csv
import json
from DataAllFormulas import allFormList4

for form in allFormList8:
	print(form)

exit()

with open('allFormulas.csv', 'r') as csvfile:
	spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')
	for row in spamreader:
		#print(row)
		saveList = []
		for elem in row[1]:
			if elem.isdigit():
				saveList.append(int(float(elem)))
		#print(saveList)
		lista2 = json.load(row[2])
		print(lista2)
		#print('lista:  ',row[1])
		#print (', '.join(row))

exit()


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
