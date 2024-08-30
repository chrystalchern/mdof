import numpy as np


def stack_powers(A,n_pwr, axis=None):
    """
    Construct array storing the powers of `A`.
    """
    n,_ = A.shape
    powers_A = np.empty((n_pwr,n,n))
    A_p = np.eye(n)
    for pwr in range(0,n_pwr):
        powers_A[pwr] = A_p
        A_p = A@A_p
    return powers_A


def form_observability(A,C,no):
    """
    Construct observability matrix of order `no` from `A` and `C` system matrices.
    """
    p,n = C.shape
    Observability = C@stack_powers(A, no)
    Observability = Observability.reshape((no*p,n))
    return Observability


def form_controllability(A,B,nc):
    """
    Construct controllability matrix of order `nc` from `A` and `B` system matrices.
    """
    n,q = B.shape
    # Controllability = np.empty((n,q,nc))
    # Controllability[:,:,0] = B
    # A_pwr = A
    # for pwr in range(1,nc):
    #     Controllability[:,:,pwr] =  A_pwr@B
    #     A_pwr = A@A_pwr
    # Controllability = Controllability.reshape((n,q*nc))

    # TODO: check form_controllability
    import warnings 
    warnings.warn("form_controllability is not yet tested")

    Controllability = np.array([Ap@B for Ap in stack_powers(A,nc)])
    Controllability = Controllability.reshape((n,q*nc))
    return Controllability


def block_hankel(series):
    pass
