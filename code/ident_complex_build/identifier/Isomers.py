from identifier import BuildTargetComplex
from identifier import Rmsd

def findisomer(isomerList, csd_xyz, csd_colors):
    allIsomers = Isomers(isomerList)
    allRmsd = {}
    for isoName, isoColors, isoChelates, isoGeo in allIsomers:
        isoCoords, isoColorsDent = BuildTargetComplex.addChelation(isoGeo, isoColors, isoChelates)
        rmsd = Rmsd.MarquesMethod(csd_xyz, csd_colors, isoCoords, isoColorsDent)
        allRmsd[rmsd] = isoName
    minRmsd = min(list(allRmsd.keys()))
    return allRmsd[minRmsd]

class Isomers:
    def __init__(self, isomerList):
        self.shape = isomerList['shape']
        self.coordinationNumber = isomerList['ncoord']
        self.composition = isomerList['composition']
        self.colors = isomerList['colors']
        self.chelates = isomerList['chelates']
        self.isomers = isomerList['isomers']
        self.isomerNames = iter(isomerList['isomers'].keys())
        self.refGeometry = isomerList['refGeometry']
            
    def __iter__(self):
        self.isomerNames = iter(self.isomers.keys())
        return self

    def __next__(self):
        name = next(self.isomerNames)
        permutation = self.isomers[name]
        colors = self.colors.copy()
        chelates = self.chelates.copy()
        refGeometry = self.refGeometry.copy()
        isomerColors = BuildTargetComplex.applySymmetry(colors, permutation)
        isomerChelates = BuildTargetComplex.applySymmetryChelates(chelates, permutation)
        return name, isomerColors, isomerChelates, refGeometry
