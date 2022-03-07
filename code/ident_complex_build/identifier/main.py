import numpy as np

from identifier import Extraction
from identifier import CompareStructures
from identifier import BuildTargetComplex
from identifier.Isomers import findisomer
from identifier import ReadStereoisomerList

def identify(mol2Stream, permutations = True):
    try:
        payload = run_identify(mol2Stream, permutations)
        return payload
    except Exception as e:
        return {'complexes': [], 'error': str(e)}

def run_identify(mol2Stream, permutations):

    payload = Extraction.refine('dummy.mol2', mol2Stream)
    if 'error' in payload:
        return payload
   
    # add shape
    for i in payload['complexes'].keys():
        lcoords = payload['complexes'][i]['LigandsCoords']
        csd_shape, rmsd_value = CompareStructures.findShape(np.array(lcoords))
        payload['complexes'][i]['shape'] = csd_shape
        payload['complexes'][i]['shape_rmsd'] = rmsd_value
    
    if not permutations:
        return payload
    
    # add xyz_pchelates
    for i in payload['complexes'].keys():
        lcoords = payload['complexes'][i]['LigandsCoords']
        colors = payload['complexes'][i]['localPriorities']
        chelates = payload['complexes'][i]['chelates']  
        xyz_pchealtes, colors_pchelates = BuildTargetComplex.addChelation(lcoords, colors, chelates)
        payload['complexes'][i]['xyz_pchealtes'] = xyz_pchealtes.tolist()
        payload['complexes'][i]['colors_pchelates'] = colors_pchelates.tolist()

    # add final isomer
    for i in payload['complexes'].keys():
        xyz_pchealtes = np.array(payload['complexes'][i]['xyz_pchealtes'])
        colors_pchelates = payload['complexes'][i]['colors_pchelates']
        shape = payload['complexes'][i]['shape']
        composition = payload['complexes'][i]['composition']
        isomerList = ReadStereoisomerList.defineIsomerList(shape, composition,'cpp')
        finalIsomer = findisomer(isomerList, xyz_pchealtes, colors_pchelates)
        payload['complexes'][i]['FinalIsomer'] = finalIsomer
    
    
    return payload
    


