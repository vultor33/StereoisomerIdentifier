from Enumeration import Enumeration
from DataAllFormulas import allFormList6

allFormList = allFormList6

geoCode = 60

obj1 = Enumeration(geoCode)
#obj1.readRotations(geoCode)
#obj1.readScaffold(geoCode)

for i in range(len(allFormList)):
	if i != 9:
		continue
	if allFormList[i][2] == ['None']:
		obj1.makeEnumeration(allFormList[i][0], allFormList[i][1], [])
	else:
		obj1.makeEnumeration(allFormList[i][0], allFormList[i][1], allFormList[i][2])
