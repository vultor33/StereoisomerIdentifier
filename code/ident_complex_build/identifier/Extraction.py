from identifier.MolFileHandling import MolFileHandling
from identifier.ExternalPrioritiesObtainer import ExternalPrioritiesObtainer
from identifier.PrioritiesObtainer import PrioritiesObtainer
from identifier.FormulaHandling import applyNewCIP
from identifier import BuildTargetComplex
import numpy as np
import math
import glob


def refine(fileMol2Name, fileStream = ""):
    molFileHandling_ = MolFileHandling(fileMol2Name, fileStream)
    metalsList = molFileHandling_.getMetalsInMol2FileList()
    try:
        extPrior_ = externalPriorities(molFileHandling_)
        externalPrior = extPrior_.getPriorities()
    except:
        return {'complexes': [], 'error': 'Rdkit error, try another complex'}    
    response = {}
    response['complexes'] = {}
    
    for imetal in metalsList:
        response['complexes'][imetal] = {}
        metalCoordinates = molFileHandling_.getListAtoms()[imetal]
        metalCoordinates = metalCoordinates.split()
        response['complexes'][imetal]['MetalSymbol'] = metalCoordinates[1]
        response['complexes'][imetal]['MetalPos'] = imetal + 1

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
    
    return response


def prioritiesOnLigands(fileStream, connectedToMetal, chelations):
    fileMol2Name = ' .mol2' # dummy - info is passed with a stream
    fileStream = fileStream.splitlines()
    molFileHandling_ = MolFileHandling(fileMol2Name, fileStream)
    try:
        extPrior_ = externalPriorities(molFileHandling_)
        externalPrior = extPrior_.getPriorities()
    except:
        return {'error': 'Rdkit error, try another complex'}    

    oldCip = []
    for iconnect in connectedToMetal:
        oldCip.append(externalPrior[iconnect])
    newCip = applyNewCIP(oldCip, chelations)
    return {'originalCip': oldCip, 'newCip': newCip}


    response = {}
    response['complexes'] = {}
    print('Priorities: ', externalPrior)



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



#####################################################################################################
#####################################################################################################
# COMPLEX BUILD SPECIFIC FUNCTIONS ##################################################################
#####################################################################################################
#####################################################################################################


def sortLigands(cip, connectedToMetal, chelations):
    taken = []
    ligands = {}
    lk = 0
    priorSort = []
    for chel in chelations:
        ligands[lk] = {}
        ligands[lk]['connectPositions'] = chel
        ligands[lk]['priorities'] =  [cip[c] for c in chel]
        ltype = 'AB'
        priorSortTerm = 10000 + min(ligands[lk]['priorities'])
        if ligands[lk]['priorities'][0] == ligands[lk]['priorities'][1]: # BIDENTATES ONLY
            ltype = 'AA'
            priorSortTerm = 100 + min(ligands[lk]['priorities'])
        ligands[lk]['ligandType'] = ltype
        priorSort.append(priorSortTerm)
        taken += chel
        lk +=1

    for c in range(len(connectedToMetal)):
        if c not in taken:
            ligands[lk] = {}
            ligands[lk]['connectPositions'] = [c]
            ligands[lk]['priorities'] = [cip[c]]
            ligands[lk]['ligandType'] = 'a'
            priorSort.append(cip[c])
            lk += 1
    
    sortedLigands = {}
    sortedList = np.argsort(priorSort)
    for i in range(len(sortedList)):
        sortedLigands[i] = dict(ligands[sortedList[i]])

    return sortedLigands
    
    
        











    
