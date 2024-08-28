import numpy as np


def form_observability(A,C,no):
    """
    Construct observability matrix of order `no` from `A` and `C` system matrices.
    """
    p,n = C.shape
    Observability = np.empty((no,p,n))
    Observability[0,:,:] = C
    A_pwr = A
    for pwr in range(1,no):
        Observability[pwr,:,:] =  C@A_pwr
        A_pwr = A@A_pwr
    Observability = Observability.reshape((no*p,n))
    return Observability


def form_controllability(A,B,nc):
    """
    Construct controllability matrix of order `nc` from `A` and `B` system matrices.
    """
    n,q = B.shape
    Controllability = np.empty((n,q,nc))
    Controllability[:,:,0] = B
    A_pwr = A
    for pwr in range(1,nc):
        Controllability[:,:,pwr] =  A_pwr@B
        A_pwr = A@A_pwr
    Controllability = Controllability.reshape((n,q*nc))
    return Controllability


def block_hankel(series):
    pass
