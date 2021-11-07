from MolFileHandling import MolFileHandling
from ExternalPrioritiesObtainer import ExternalPrioritiesObtainer
from PrioritiesObtainer import PrioritiesObtainer
import numpy as np
import math
import glob

def refineCoordinationPolyhedron(fileMol2Name):
    prior_, molFileHandling_ = extractCoordinationPolyhedron(fileMol2Name)
    name = molFileHandling_.getBaseFileName()
    colors = prior_.getPrioritesOfMetalI()
    composition = prior_.getDirectFormula()
    chelates = prior_.getChelatesOfMetalI()

    lcoords = []
    metalCoordinates = molFileHandling_.getListAtoms()[molFileHandling_.getMetalsInMol2FileList()[0]]
    metalCoordinates = metalCoordinates.split()
    lcoordsi = [float(x) for x in metalCoordinates[2:5]]
    lcoords.append(lcoordsi)
    
    for i in prior_.getLigandsBondedToMetalI():
        donorAtomsCoordinates = molFileHandling_.getListAtoms()[i-1].split()
        lcoordsi = [float(x) for x in donorAtomsCoordinates[2:5]]
        lcoords.append(lcoordsi)
    lcoords = centerAndNormalize(lcoords)
    return name, lcoords, colors, composition, chelates

def extractCoordinationPolyhedron(fileMol2Name):
    molFileHandling_ = MolFileHandling(fileMol2Name)
    molstream = molFileHandling_.writemolstream()
    extPrior_ = ExternalPrioritiesObtainer(molFileHandling_.getTemporaryMolName(), molstream)
    prior_ = PrioritiesObtainer(molFileHandling_, extPrior_.getPriorities())
    metalsList = molFileHandling_.getMetalsInMol2FileList()
    prior_ = generatePriorities(self, molFileHandling_, extPrior_, metalsList[0])
    return prior_, molFileHandling_

def generatePriorities(self, molFileHandling_, extPrior_, iMetal):
    prior_ = PrioritiesObtainer(molFileHandling_, extPrior_.getPriorities())
    prior_.calculateLigandsPriorities(iMetal)
    return prior_

def centerAndNormalize(lcoords):
    X = []
    for line in lcoords[1:]:
        line = np.array(line)
        line -= np.array(lcoords[0])
        X.append(line)
    X = np.array(X)
    for i in range(len(X)):
        norm = math.sqrt(X[i].dot(X[i]))
        X[i] = X[i]/norm        
    return X
