import numpy as np

import Rmsd


def findShape(X, geometries):
    ncoord = X.shape[0]
    colors = np.zeros(ncoord)
    minRmsd = 1e10
    shape = ''
    for geo in geometries:
        if geometries[geo]['ncoord'] == ncoord:
            coordinates = geometries[geo]['coordinates']
            rmsd = Rmsd.MarquesMethod(X, colors, coordinates, colors)
            if rmsd < minRmsd:
                shape = geo
                minRmsd = rmsd
    
    if len(shape) == 0:
        raise Exception('Shape not found, check input for errors')
    return shape, minRmsd
