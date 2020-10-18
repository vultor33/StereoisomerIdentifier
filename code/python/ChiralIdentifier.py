import numpy as np
import Rmsd

ZERO_THRESHOULD = 1e-2

def isChiral(X, colors):
    Y = X.copy()
    Y[:,1] *= -1  # XZ plane reflection
    rmsd = Rmsd.MarquesMethod(X, colors, Y, colors)
    return rmsd > ZERO_THRESHOULD

def test_isChiral():
    X = np.array([[-0.57735027, -0.57735027, -0.57735027],
        [-0.57735027,  0.57735027,  0.57735027],
        [ 0.57735027,  0.57735027, -0.57735027],
        [ 0.57735027, -0.57735027,  0.57735027]]) 
    assert(isChiral(X, [1,1,1,1]) == False)
    assert(isChiral(X, [1,2,3,4]))


