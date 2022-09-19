from cmath import isclose
from numpy import pi, log, sign
from numpy.linalg import eig
import numpy as np
import scipy.linalg as sl

def condeig(a):
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
    vl = vl.T
    # Normalize vectors
    for i in range(n):
        vl[i, :] = vl[i, :] / np.sqrt(abs(vl[i, :] ** 2).sum())
    # Condition numbers are reciprocal of the cosines (dot products) of the
    # left eignevectors with the right eigenvectors.
    c = abs(1 / np.diag(np.dot(vl, vr)))
    return vr, lamr, c

def ComposeModes(dt, A, B, C, D, debug=False, **kwds)->dict:
    
    n = A.shape[0]
    m = C.shape[0]

    v, d, cnd = condeig(A)  # eigenvectors (d) & eiegenvalues (v) of the matrix A
    # kit = np.log(d)         # logarithm of the eigenvalues

    # a) Determination of modal frequencies (Eqs. 3.46 & 3.39)
    sj1 = np.log(d)/dt            # dt is the time step
    freq1 = (np.real(sj1*np.conj(sj1))**0.5)/(2*pi)

    roots = []

    # selection of proper roots
    if np.isclose(freq1[0], freq1[1]):
        roots.append(0)

    if np.isclose(freq1[n-1], freq1[n-2]):
        roots.append(n-1)

    for i in range(2, n-2):
        if np.isclose(freq1[i], freq1[i+1]) or np.isclose(freq1[i], freq1[i-1]):
            roots.append(i)
 

    # b) Determination of damping ratios (Eqs. 3.46 & 3.39)
    damp1 = -((np.real(sj1))/(2*pi*freq1))

    # Represent the identified frequency & damping information
    # of the proper roots in a matrix

    # NOTE: These values are stored in one array like this for legacy reasons,
    # but this should be cleaned and a proper data structure should be used!
    freqdmp = np.array([
               [freq1[i],   # first column: identified frequency
                damp1[i],   # second column: identified damping ratio
                cnd[i]]     # condition number of the eigenvalue
            for i in roots
    ])

    # c) Determination of mode shapes
    modes_raw = C@v   # mode shapes (Eq. 3.40)

    modeshape = np.zeros((m,len(roots)), dtype=complex)

    # extract mode shapes from mod corresponding to a frequency
    for q,root in enumerate(roots):
        modeshape[:m,q] = modes_raw[:m, root]

    for q,root in enumerate(roots):
        om  = np.argmax(abs(np.real(modeshape[:,q])))
        mx  = abs(np.real(modeshape[om,q]))
        modeshape[:,q] = np.real(modeshape[:,q])/mx*sign(np.real(modeshape[om,q]))

    if debug:
        return locals()
    return freqdmp, modeshape, None, v, d


