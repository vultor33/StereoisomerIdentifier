import json
import matplotlib.pyplot as plt
import numpy as np
from collections import OrderedDict

with open('geometries2.json', 'w') as file:
    json.dump(json_data, file)
    
with open('geometries.json', 'r') as j:
    json_data = json.load(j, object_pairs_hook=OrderedDict)
    
# CODE TO ADD ROTATIONS TO JSON GEOMETRIES
geoCode = 90
geoFlag = 'TCTPR-9'
filename = 'rotations/rotations-' + str(geoCode) + '.txt'
rot_flags = list(json_data[geoFlag]['rotations'])
f = open(filename, "r")
rot_lines = f.readlines()
for i in range(len(rot_flags)):
    line = rot_lines[i].rstrip().split('  ')
    line = [int(x) for x in line]
    json_data[geoFlag]['rotations'][rot_flags[i]] = line
        
# CODE TO ADD COORDINATION NUMBER        
for cpcode in json_data:
    ncoord = int(cpcode.split('-')[1])
    json_data[cpcode]['ncoord'] = ncoord
    
# DEFINICAO DE UM JSON QUE MAPEIA AS COMPOSICOES DO PYTHON COM AS DO CPP
folder = os.path.join('StereoisomerListCPP', 'CN10')
folder = os.path.join(folder, 'TD-10')
_, _, filenames = next(os.walk(folder))
cpp_python_map = {}
for file in filenames:
    if 'counting' in file:
        continue
    fName = os.path.join(folder, file)
    f = open(fName, "r")
    lines = f.readlines()
    colors = ReadStereoisomerList.readColorsCPP(lines[0])
    chelates = ReadStereoisomerList.readChelatesCPP(lines[0])
    fHandling_ = FormulaHandling.FormulaHandling()
    fHandling_.generateMolecularFormula(colors, chelates)
    composition = fHandling_.getFormula()
    cpp_python_map[composition] = file.split('-')[2]
with open('cpp_composition_map.json', 'w') as file:
    json.dump(cpp_python_map, file)  # TRANSFERIR ESSE ARQUIVO PARA A PASTA 'StereosiomerListCpp'
    
    
    
# COORDENACAO 10
import json
import numpy as np
import os
import re
import math
from IPython.display import Image

import Rmsd
import PrintMol2
import CompareStructures
import BuildTargetComplex
import ReadStereoisomerList
import FormulaHandling
from Isomers import Isomers
import CsdExtraction

def isbidentatemax(chelates):
    for chel in chelates:
        if len(chel) > 2:
            return False
    return True
    
_, _, filenames = next(os.walk('G:\\!CSD-database'))

coord10 = {}

for file in filenames:
    fileMol2Name = os.path.join('G:\\!CSD-database', file)
    try:
        name, lcoords, colors, composition, chelates = CsdExtraction.refineCoordinationPolyhedron(fileMol2Name)
    except:
        continue
        
    if len(lcoords) != 10:
        continue
    if not isbidentatemax(chelates):
        continue
        
    try:
        shape, _ = CompareStructures.findShape(lcoords)
    except:
        pass
    
    coord10[fileMol2Name] = {}
    coord10[fileMol2Name]['isomer'] = [name, lcoords, colors, composition, chelates]
    coord10[fileMol2Name]['shape'] = shape

for filename in coord10:
    name, lcoords, colors, composition, chelates = coord10[filename]['isomer']
    try:
        lcoords.tolist()
    except:
        continue
    coord10[filename]['isomer'] = name, lcoords.tolist(), colors, composition, chelates
    
with open('coord10.json', 'w') as file:
    json.dump(coord10, file)



