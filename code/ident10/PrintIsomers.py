import ReadStereoisomerList
from Isomers import Isomers
import PrintMol2

def print_isomers(shape, composition):
    isomerList = ReadStereoisomerList.defineIsomerList(shape, composition, 'cpp')
    allIsomers = Isomers(isomerList)
    for name, colors, chelates, geo in allIsomers:
        PrintMol2.printMol2(name + '.mol2', colors, chelates, geo)