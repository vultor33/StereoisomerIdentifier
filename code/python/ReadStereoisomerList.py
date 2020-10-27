import re
import os
import numpy as np


###################################################################################
# FORMAT SPECIFIC READ ############################################################
###################################################################################

def defineIsomerList(shape, composition, geometries):
    ncoord = len(geometries[shape]['coordinates'])
    lines = readFile(shape, ncoord, composition)
    isomerList = {}
    isomerList['shape'] = shape
    isomerList['ncoord'] = ncoord
    isomerList['composition'] = composition
    isomerList['colors'] = readColors(lines[0])
    isomerList['chelates'] = readChelates(lines[0])
    isomerList['isomers'] = definePermutations(lines)
    return isomerList

def readFile(shape, ncoord, composition):
    stereo_path = findStereoisomerlistFile(shape, ncoord, composition)
    f = open(stereo_path, "r")
    return f.readlines()

def definePermutations(lines):
    isomerPermutation = {}
    for i in range(1, len(lines)):
        line = lines[i]
        if line[0] == 'R' or line[0] == 'S' or line[0] == 'G':
            isomerType = line[0]
            isomerCount = 0
            continue
        else:
            isomerCount += 1
            permutation = line.split(';')[0].rstrip().split(' ')
            permutation = [int(x) for x in permutation]
            isomerName = isomerType + '-' + str(isomerCount) + ' '
            isomerPermutation[isomerName] = permutation
    return isomerPermutation

def findStereoisomerlistFile(shape, ncoord, composition):
    folder = 'Stereoisomerlist'
    folder = os.path.join(folder, 'CN' + str(ncoord))
    folder = os.path.join(folder, shape)
    folder = os.path.join(folder,shape + '-' + composition + '.csv')
    return folder

def readColors(line0):
    comp = line0.rstrip().split('  ')
    comp = comp[3].split(':')[1]
    comp = comp.lstrip().split(' ')
    comp = [int(x) for x in comp]
    return np.array(comp)

def readChelates(line):
    line = line.rstrip()
    line = re.sub('\s+',' ',line)
    line = line.split(' ')
    chel = []
    i = 0
    while i < len(line):
        if line[i] == 'cI:':
            chelI = []
            for j in range(i + 1, len(line)):
                if line[j] == 'cI-length:': break
                chelI.append(int(line[j]))
            chel.append(chelI)
            i = j
        i += 1
        
    return chel




