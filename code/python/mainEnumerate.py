from Enumeration import Enumeration
from DataAllFormulas import allFormList6

allFormList = allFormList6

geoCode = 60

obj1 = Enumeration()
obj1.readRotations(geoCode)
obj1.readScaffold(geoCode)

for i in range(len(allFormList)):
	if i != 6:
		continue
	if allFormList[i][2] == ['None']:
		obj1.makeEnumeration(allFormList[i][1], [])
	else:
		obj1.makeEnumeration(allFormList[i][1], allFormList[i][2])
