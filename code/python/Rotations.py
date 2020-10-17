import numpy as np
import math


#############################################################################################################
# EULER #####################################################################################################
#############################################################################################################


def ApplyEulerAngles(
        X, 
        EulerAngles, 
        center = np.zeros(3)):
    '''
    Rotates geometry X around the Euler Angles
    Inputs
    o X: n,3 matrix. where rows are atoms and columns are (x,y,z) coordinates 
    o EulerAngles: euler angles in radians
    o center: reference translation to apply rotation matrix. Use X center to make specific rotations.
    Outputs
    o X: rotated X matrix
    '''

    theta,psi,phi = EulerAngles
    costheta = np.cos(theta)
    sintheta = np.sin(theta)
    cospsi = np.cos(psi)
    sinpsi = np.sin(psi)
    cosphi = np.cos(phi)
    sinphi = np.sin(phi)

    Rx = np.zeros((3,3))
    Ry = np.zeros((3,3))
    Rz = np.zeros((3,3))

    Rx[0,0] =  1.
    Rx[1,1] =  costheta
    Rx[1,2] =  sintheta
    Rx[2,1] = -sintheta
    Rx[2,2] =  costheta
    
    Ry[0,0]=  cospsi
    Ry[0,2]= -sinpsi
    Ry[1,1]=  1.
    Ry[2,0]=  sinpsi
    Ry[2,2]=  cospsi

    Rz[0,0]=  cosphi
    Rz[0,1]=  sinphi
    Rz[1,0]= -sinphi
    Rz[1,1]=  cosphi
    Rz[2,2]=  1. 
    
    R = np.linalg.multi_dot([Rz,Ry,Rx])

    X += -center
    X  = np.dot(X,R)
    X += +center
    return X


def ApplyRandomRotation(X):
    theta = np.random.uniform(-1,1)
    psi = np.random.uniform(-1,1)
    phi = np.random.uniform(-1,1)
    EulerAngles = np.array([theta,psi,phi])*2*math.pi
    return ApplyEulerAngles(X.copy(), EulerAngles)

    
        
def test_ApplyEulerAngles():   
    euler = (0, 0, 0.4*math.pi)
    X = np.zeros((4,3))
    X[1,0] = 1
    X[2,0] = 2
    X[3,1] = 1
    X -= X.mean(axis=0)
    NormCons = 1./np.linalg.norm(X,axis=1)
    NVertices = X.shape[0]
    NormCons.resize(NVertices,1)
    X = NormCons*X
    Y = ApplyEulerAngles(X.copy(),euler)
    X_ref = np.array([[-0.9486833 , -0.31622777,  0.        ],
       [ 0.70710678, -0.70710678,  0.        ],
       [ 0.98058068, -0.19611614,  0.        ],
       [-0.70710678,  0.70710678,  0.        ]])
    Y_ref = np.array([[ 0.00759122, -0.99997119,  0.        ],
       [ 0.89100652,  0.4539905 ,  0.        ],
       [ 0.48953362,  0.87198442,  0.        ],
       [-0.89100652, -0.4539905 ,  0.        ]])
    assert(np.allclose(X,X_ref))
    assert(np.allclose(Y,Y_ref))    


#############################################################################################################
# QUATERNIONS ###############################################################################################
#############################################################################################################

def QuaternionAlign(X,Y):
    ''' Given two sets of points find a rotation that minimizes 
        the RMSD of the two sets.
    INPUTS
    o X - reference geometry
    o Y - geometry to rotate
    OUTPUTS
    o rmsd - final rmsd between structures
    o Yprime - rotated Y structure'''
        
    rmsd,R = QuaternionRotation(X,Y)
    Yprime = np.dot(Y,R) 
    return [rmsd,Yprime]


def QuaternionRotation(X,Y):
    N = X.shape[0]
    # X and Y must be centered at the origin of the coordinate system.

    ## CONSTRUCT KEARSLEY ALIGNMENT MATRIX ##
    # KEARSLEY, S. On the orthogonal transformation used for structural comparison. Acta Cryst. (1989). A45, 208-210
    Xm = Y[:,0]-X[:,0]
    Xp = Y[:,0]+X[:,0]
    
    Ym = Y[:,1]-X[:,1]
    Yp = Y[:,1]+X[:,1]
    
    Zm = Y[:,2]-X[:,2]
    Zp = Y[:,2]+X[:,2]
    
    Q = np.zeros((4,4))

    # MAIN DIAGONAL
    Q[0,0] = (Xm*Xm + Ym*Ym + Zm*Zm).sum()
    Q[1,1] = (Xm*Xm + Yp*Yp + Zp*Zp).sum()
    Q[2,2] = (Xp*Xp + Ym*Ym + Zp*Zp).sum()
    Q[3,3] = (Xp*Xp + Yp*Yp + Zm*Zm).sum()
    
    
    ## UPPER TRIANGLE ##

    # COLUMN 0
    
    # COLUMN 1
    Q[0,1] = (Yp*Zm - Ym*Zp).sum()

    # COLUMN 2
    Q[0,2] = (Xm*Zp - Xp*Zm).sum()
    Q[1,2] = (Xm*Ym - Xp*Yp).sum()
    
    # COLUMN 3
    Q[0,3] = (Xp*Ym - Xm*Yp).sum()
    Q[1,3] = (Xm*Zm - Xp*Zp).sum()
    Q[2,3] = (Ym*Zm - Yp*Zp).sum()



    ## LOWER TRIANGLE ##

    # COLUMN 0
    Q[1,0] = Q[0,1]
    Q[2,0] = Q[0,2]
    Q[3,0] = Q[0,3]
    
    # COLUMN 1
    Q[2,1] = Q[1,2]
    Q[3,1] = Q[1,3]

    # COLUMN 2
    Q[3,2] = Q[2,3]
    
    # COLUMN 3
    
    # OBTAIN QUATERNION EIGENVECTORS
    eigenvals, eigenvecs = np.linalg.eigh(Q)


    ## CONSTRUCT ROTATION MATRIX ##
    R = np.zeros((3,3))
    a = 0
    l = eigenvals[a]
    if abs(l) <= 1e-7:    
        l = 0
    elif l < 0: # TODO - raise warning
        print(Q)
        print('eigenvals: ', eigenvals,l,N)
    rmsd = np.sqrt(l/N)
    q1,q2,q3,q4 = eigenvecs[:,a] # Eigenvector with lower eigenvalue
    #q4,q3,q2,q1 = eigenvecs[:,a] # Eigenvector with lower eigenvalue
    #r = Rotation.from_quat([q4,q3,q2,q1])
    #R = r.as_matrix()
    #quit()

    q11 = q1*q1
    q22 = q2*q2
    q33 = q3*q3
    q44 = q4*q4
    
    q12 = q1*q2
    q13 = q1*q3
    q14 = q1*q4
    
    q23 = q2*q3
    q24 = q2*q4
    
    q34 = q3*q4
    
    
    # MAIN DIAGONAL
    R[0,0] = q11 + q22 - q33 - q44
    R[1,1] = q11 + q33 - q22 - q44
    R[2,2] = q11 + q44 - q22 - q33

    # UPPER TRIANGLE
    R[0,1] = 2*(q23+q14)
    
    R[0,2] = 2*(q24-q13)
    R[1,2] = 2*(q34+q12)
    
    # LOWER TRIANGLE
    R[1,0] = 2*(q23-q14)
    R[2,0] = 2*(q24+q13)

    R[2,1] = 2*(q34-q12)
 
    return [rmsd,R]

def test_QuaternionAlign():
    X = np.array([[-0.9486833 , -0.31622777,  0.        ],
        [ 0.70710678, -0.70710678,  0.        ],
        [ 0.98058068, -0.19611614,  0.        ],
        [-0.70710678,  0.70710678,  0.        ]])
    Y = np.array([[ 0.00759122, -0.99997119,  0.        ],
        [ 0.89100652,  0.4539905 ,  0.        ],
        [ 0.48953362,  0.87198442,  0.        ],
        [-0.89100652, -0.4539905 ,  0.        ]])
    rmsd, Yprime = QuaternionAlign(X, Y)
    Yprime_ref = np.array([[-0.9486833 , -0.31622777,  0.        ],
        [ 0.70710678, -0.70710678,  0.        ],
        [ 0.98058068, -0.19611614,  0.        ],
        [-0.70710678,  0.70710678,  0.        ]])
    assert(np.allclose(Yprime, Yprime_ref))    











