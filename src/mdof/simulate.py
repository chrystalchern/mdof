import numpy as np
import scipy

def _simulate(system, u, w=None):
    (A, B, C, D) = system
    p,n = C.shape
    T, _ = u.shape
    # Initialize
    x = np.empty((T,n))

    x[0] = 0
    for t in range(T-1):
        x[t+1] = A.dot(x[t]) + B.dot(u[t])

    y = C.dot(x.T).T + D.dot(w.T).T

    return x, y

def simulate(system, u, w=None, backend=None):
    (A, B, C, D) = system
    if backend is None:
        return scipy.signal.dlsim(scipy.signal.dlti(A, B, C, D), u[:, j])
    
    else:
        return _simulate((A, B, C, D), u, w=w)[1]

