from Enumeration import Enumeration
from FormulaHandling import FormulaHandling
from FormulaHandling import GenerateAllFormulas
from collections import Counter
import itertools
import math


objG = GenerateAllFormulas()

objG.generateAllFormulas(4)

objG.printFile()

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
