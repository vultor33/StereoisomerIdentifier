import scipy
from scipy import optimize
import numpy as np

from Rotations import QuaternionAlign, ApplyRandomRotation

MAX_ITER = 300
BIG_NUMBER = 1e10
COLOR_LIST = ['red','blue','green','magenta']  # TODO - auto-generate this colors

def MarquesMethod(X, colorsX, Y, colorsY):
    '''
    This function implements the method shown in:
    o How Different Are Two Chemical Structures? J. Chem. Inf. Model. 2010, 50, 2129–2140 
	  DOI: 10.1021/ci100219f
	
    INPUTS
        o X[n - atoms rows, 3 - xyz coordinates] - first molecule
        o colorsX - vector of atoms chemical composition for molecule X
        o Y[n - atoms rows, 3 - xyz coordinates] - second molecule
        o colorsX - vector of atoms chemical composition for molecule Y
    OUTPUTS
        o minRmsd - minimum rmsd of X and Y molecules.	
    '''
	
    minRMSD = BIG_NUMBER
    for i in range(MAX_ITER):
        Yrotated = ApplyRandomRotation(Y)
        Yhungarian = ApplyHungarianPermutation(X, colorsX, Yrotated, colorsY)
        RMSDMeasure, Yalligned = QuaternionAlign(X, Yhungarian)
    
        if np.isnan(RMSDMeasure): continue 
        if RMSDMeasure > minRMSD:
            continue
        else:
            minRMSD = RMSDMeasure
    
    return minRMSD


def ApplyHungarianPermutation(X, colorsX, Y, colorsY):
    CostMatrix = GenerateCostMatrix(X, colorsX, Y, colorsY)    
    rows, cols= scipy.optimize.linear_sum_assignment(CostMatrix)
    Yassigned = Y.copy()
    Yassigned[rows] = Y[cols]
    return Yassigned

# GENERATE COST MATRIX USING A MODIFIED DISTANCE MATRIX 
# no caso do gabriel e a mesma molecula entao colorsX == colorsY
def GenerateCostMatrix(X, colorsX, Y, colorsY):
    N = X.shape[0]
    M = Y.shape[0]        
    CostMatrix = np.zeros((N,M))
    for n in range(N):
        for m in range(M):
            p = X[n]
            q = Y[m]
            cp = colorsX[n]
            cq = colorsY[m]
            r = np.linalg.norm(p-q)
                
            if cp==cq: #Do not pair vertices that are not of the same type
                CostMatrix[n,m]=r
            else:
                # TODO - check this r multiplication
                CostMatrix[n,m]=1e10*r #hence why pairs of different vertex types get a strong cost
    return CostMatrix
    
    
def print_chimera(X, Y, colors):
    '''OUTPUTTING BILD FILE FOR CHIMERA'''
    referenceBildFile = open("refpoly.bild",'w')
    rotatedBildFile = open("rotpoly.bild",'w')
            
    line = ".color white\n"
    line += ".dot 0 0 0\n"
    referenceBildFile.write(line)
    rotatedBildFile.write(line)
    for n in range(poly.shape[0]):
        color = colors[n]
        color = COLOR_LIST[color]
                
        line = ""
        x,y,z = X[n]
        line += ".color {}\n".format(color)
        line += ".dot {} {} {}\n".format(x,y,z)
        line += ".v {} {} {} 0 0 0\n".format(x,y,z)
        referenceBildFile.write(line)
                
        line = ""               
        x,y,z = Y[n]
        line += ".color {}\n".format(color)
        line += ".dot {} {} {}\n".format(x,y,z)
        line += ".v {} {} {} 0 0 0\n".format(x,y,z)
        rotatedBildFile.write(line)
    referenceBildFile.close()
    rotatedBildFile.close()    
    

def test_MarquesRmsd():
    ZERO_TOL = 1e-7
    X = np.array([[-0.95665764, -0.18814179,  0.        ],
       [ 0.69913244, -0.57902081,  0.        ],
       [ 0.97260633, -0.06803016,  0.        ],
       [-0.71508113,  0.83519276,  0.        ]])
    Y = np.array([[-0.95665764,  0.18814179,  0.        ],
       [ 0.69913244,  0.57902081,  0.        ],
       [ 0.97260633,  0.06803016,  0.        ],
       [-0.71508113, -0.83519276,  0.        ]])
    
    rmsd0 = MarquesMethod(X, [0,0,0,0], Y, [0,0,0,0])
    assert(abs(rmsd0) <= ZERO_TOL)
    rmsd1 = MarquesMethod(X, [1,0,0,0], Y, [0,1,0,0])
    assert(abs(rmsd1 - 0.23601106351121656) <= ZERO_TOL)
       
       

