import numpy as np
import math

def build(geometry, permutation, colors, chelates):
    colors = applySymmetry(colors, permutation)
    chelates = applySymmetryChelates(chelates, permutation)
    targetGeo, targetColors = addChelation(geometry, colors, chelates)
    return targetGeo, targetColors
        
def applySymmetry(geometry, symm):
    geoSymm = geometry.copy()
    for i in range(len(symm)):
        geoSymm[i] = geometry[symm[i]]
    return geoSymm

def applySymmetryChelates(chelates, permutation):
    chelatesNew = []
    for chel in chelates:
        chel = [permutation.index(x) for x in chel]
        chelatesNew.append(chel)
    return chelatesNew

def addChelation(X, colors, chelates = []):
    targetCp = np.array(X).copy()
    ncoord = targetCp.shape[0]
    targetCp = addChelationPoints(targetCp, chelates)
    nchelatePoints = targetCp.shape[0] - ncoord
    targetColors = colors.copy()
    targetColors = np.append(targetColors, -1 * np.ones(nchelatePoints, dtype=int))
    return targetCp, targetColors

def addChelationPoints(lcoords, chelates):
    for chel in chelates:
        chel_size = len(chel)
        for i in range(chel_size):
            chelI = chel[i]
            for j in range(i + 1, chel_size):
                chelJ = chel[j]
                p1, p2 = create_dummy_points(lcoords[chelI], lcoords[chelJ])
                lcoords = np.vstack([lcoords,p1,p2])
    return lcoords
                
def create_dummy_points(atom1, atom2):
    auxPoint = atom2 - atom1
    aux1 = atom1 + (auxPoint * 0.25)
    aux2 = atom1 + (auxPoint * 0.75)
    return aux1, aux2                


###################################################################################
# USER DEFINED ####################################################################
###################################################################################

import re

colors = ['C', 'O', 'N', 'P', 'I', 'B', 'Se', 'Ca', 'Na', 'Se']
colorDent = 'H'
colorMetal = 'Au'

def polyToXyz(isoCoords, isoColorsDent):
    xyzString = ''
    xyzString += str(len(isoColorsDent) + 1) + '\n\n'
    xyzString += 'Au  0.0  0.0  0.0\n'
    for i in range(len(isoCoords)):
        coords = str(isoCoords[i][0]) + ' ' + str(isoCoords[i][1]) + str(isoCoords[i][2])
        if isoColorsDent[i] == -1:
            xyzString += colorDent + '  ' + coords + '\n'
        else:
            xyzString += colors[isoColorsDent[i]] + '  ' + coords + '\n'
    return xyzString



def readinput(filename):
    f = open(filename, "r")
    lines = f.readlines()
    lines = [x.rstrip().lstrip() for x in lines]
    composition = lines[0]
    chelates = readchelates(lines[1])
    metalName, metalCoord, _ = linetocoordinates(lines[2])
    atomNames = []
    ligandsCoords = []
    priorities = []
    for i in range(3, len(lines) - 1):
        atomName, coord, priority = linetocoordinates(lines[i])
        atomNames.append(atomName)
        ligandsCoords.append(coord)
        priorities.append(priority)
    metalCoord = np.array(metalCoord)
    ligandsCoords = np.array(ligandsCoords)
    ligandsCoords -= metalCoord # CENTERED ON METAL
    ligandsCoords = normalize(ligandsCoords) 
    return ligandsCoords, priorities, composition, chelates

def readchelates(line):
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

def linetocoordinates(line):
    line = line.split('   ')
    atomName = line[0]
    coord = [float(line[1]), float(line[2]), float(line[3])]
    priority = int(line[4])
    return atomName, coord, priority

def normalize(X):
    NormCons = 1./np.linalg.norm(X,axis=1)
    NormCons.resize(X.shape[0],1)
    X = NormCons*X
    return X


