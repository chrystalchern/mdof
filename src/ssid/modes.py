import numpy as np
import scipy.linalg as sl
from numpy import pi

def condeig(a): # TODO: make this match matlab source code for condeig
    """
    vals, vecs, conds = condeig(A) Computes condition numbers for the
    eigenvalues of a matrix. The condition numbers are the reciprocals
    of the cosines of the angles between the left and right eigenvectors.
    Inspired by Arno Onken's Octave code for condeig.

    https://github.com/macd/rogues/blob/master/rogues/utils/condeig.py
    """
    m, n = a.shape
    # eigenvalues, left and right eigenvectors
    lamr, vl, vr = sl.eig(a, left=True, right=True)
    vl = vl.T.conj()
    # Normalize vectors
    for i in range(n):
        vl[i, :] = vl[i, :] / np.sqrt(abs(vl[i, :] ** 2).sum())
    # Condition numbers are reciprocal of the cosines (dot products) of the
    # left eignevectors with the right eigenvectors.
    c = abs(1 / np.diag(np.dot(vl, vr))) 
    return vr, lamr, c

def modes(realization, dt):

    A,_,C,_ = realization
    # eigendecomp A
    Psi,Gam,cnd = condeig(A)  # eigenvectors (Psi) & eigenvalues (Gam) of the matrix A

    # get damping and frequencies from eigendecomp of A
    Lam = (np.log(Gam))/dt
    Omega = np.real((Lam*np.conj(Lam))**0.5)  # radians per second. taking the real part because np keeps a +0j.
    freq = Omega/(2*pi) # cycles per second (Hz)
    damp = -np.real(Lam)/Omega

    # get modeshapes from C and eigendecomp of A
    modeshape = C @ Psi

    # weed out unique roots: get indices of roots that only show up once, and
    # the index of the first of each pair.
    _, notroots = np.unique(freq.round(decimals=5), return_index=True)
    
    # print(notroots)
    modes = {str(i):
                {'cnd': cnd[i],   # condition number of the eigenvalue
                'freq': freq[i],  # identified frequency
                'damp': damp[i],  # identified damping ratio
                'modeshape': modeshape[:,i]
                }
            for i in range(len(freq)) if i not in notroots
            }

    return modes
