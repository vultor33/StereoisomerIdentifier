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
	stereoRemoved = []
	rcw = []
	if allFormList[i][2] == ['None']:
		obj1.makeEnumeration(allFormList[i][1], [], rcw, stereoRemoved)
	else:
		obj1.makeEnumeration(allFormList[i][1], allFormList[i][2], rcw, stereoRemoved)
	
	print(allFormList[i][0], '   remov:  ',len(stereoRemoved))
	print(rcw)
