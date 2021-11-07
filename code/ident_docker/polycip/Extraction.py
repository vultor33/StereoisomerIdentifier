from polycip.MolFileHandling import MolFileHandling
from polycip.ExternalPrioritiesObtainer import ExternalPrioritiesObtainer
from polycip.PrioritiesObtainer import PrioritiesObtainer
from polycip import BuildTargetComplex
import numpy as np
import math
import glob


def refine(fileMol2Name):
    molFileHandling_ = MolFileHandling(fileMol2Name)
    metalsList = molFileHandling_.getMetalsInMol2FileList()
    extPrior_ = externalPriorities(molFileHandling_)    
    externalPrior = extPrior_.getPriorities()
    response = {}
    response['complexes'] = {}
    
    for imetal in metalsList:
        response['complexes'][imetal] = {}
        metalCoordinates = molFileHandling_.getListAtoms()[molFileHandling_.getMetalsInMol2FileList()[imetal]]
        metalCoordinates = metalCoordinates.split()
        response['complexes'][imetal]['MetalSymbol'] = metalCoordinates[1]
        response['complexes'][imetal]['MetalPos'] = molFileHandling_.getMetalsInMol2FileList()[imetal] + 1

        prior_ = generatePriorities(molFileHandling_, externalPrior, imetal)
        ligands_file_index = prior_.getLigandsBondedToMetalI()
        response['complexes'][imetal]['LigandsOnMol2File'] = ligands_file_index

        lcoords = []
        lcoordsi = [float(x) for x in metalCoordinates[2:5]]
        lcoords.append(lcoordsi)
        globalPriorities = []
        for i in prior_.getLigandsBondedToMetalI():
            globalPriorities.append(externalPrior[i-1])
            donorAtomsCoordinates = molFileHandling_.getListAtoms()[i-1].split()            
            lcoordsi = [float(x) for x in donorAtomsCoordinates[2:5]]            
            lcoords.append(lcoordsi)
        lcoords = centerAndNormalize(lcoords)
        response['complexes'][imetal]['localPriorities'] = prior_.getPrioritesOfMetalI()
        response['complexes'][imetal]['globalPriorities'] = globalPriorities
        response['complexes'][imetal]['composition'] = prior_.getDirectFormula()
        response['complexes'][imetal]['chelates'] = prior_.getChelatesOfMetalI()
        response['complexes'][imetal]['LigandsCoords'] = lcoords.tolist()
        isoCoords, isoColorsDent = BuildTargetComplex.addChelation(
            response['complexes'][imetal]['LigandsCoords'], 
            response['complexes'][imetal]['localPriorities'],
            response['complexes'][imetal]['chelates'])
        xyzCoords = BuildTargetComplex.polyToXyz(isoCoords, isoColorsDent)
        response['complexes'][imetal]['xyz'] = xyzCoords
    
    return response

def externalPriorities(molFileHandling_):
    molstream = molFileHandling_.writemolstream()
    extPrior_ = ExternalPrioritiesObtainer(molFileHandling_.getTemporaryMolName(), molstream)
    return extPrior_ 

def generatePriorities(molFileHandling_, externalPrior, iMetal):
    prior_ = PrioritiesObtainer(molFileHandling_, externalPrior)
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
    
    
