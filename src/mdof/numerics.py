import numpy as np
import numpy.linalg
import scipy


def lsq_solver(options):
    return lambda *args: numpy.linalg.lstsq(*args, rcond=options.get('rcond',None))[0]


def lin_solver():
    pass


def svd_routine(lapack_driver="gesvd", **kwds):
    # default lapack_driver="gesdd"
    def _svd(*args):
        U,s,Vh = scipy.linalg.svd(*args, lapack_driver=lapack_driver, **kwds)
        return U,s,Vh.T.conj()
    return _svd


def _condeig(A):
    """
    vals, vecs, conds = condeig(A) Computes condition numbers for the
    eigenvalues of a matrix. The condition numbers are the reciprocals
    of the cosines of the angles between the left and right eigenvectors.
    Inspired by Arno Onken's Octave code for condeig.

    https://github.com/macd/rogues/blob/master/rogues/utils/condeig.py
    """
    m,n = A.shape
    # eigenvalues, left and right eigenvectors
    lamr, vl, vr = scipy.linalg.eig(A, left=True, right=True)
    vl = vl.T.conj()
    # Normalize vectors
    for i in range(n):
        vl[i, :] = vl[i, :] / np.sqrt(abs(vl[i, :] ** 2).sum())
    # Condition numbers are reciprocal of the cosines (dot products) of the
    # left eignevectors with the right eigenvectors.
    c = abs(1 / np.diag(np.dot(vl, vr))) 
    return (vr,lamr,c)

from functools import lru_cache, wraps

def np_cache(function):
    @lru_cache()
    def cached_wrapper(hashable_array):
        array = np.array(hashable_array)
        return function(array)

    @wraps(function)
    def wrapper(array):
        return cached_wrapper(tuple(array))

    # copy lru_cache attributes over too
    wrapper.cache_info = cached_wrapper.cache_info
    wrapper.cache_clear = cached_wrapper.cache_clear

    return wrapper

_cached_condeig = np_cache(_condeig)


def linear_interpolate(x, y, target_x):
    sorted_indices = sorted(range(len(x)), key = lambda i: x[i])
    x = x[sorted_indices]
    y = y[sorted_indices]
    i1 = max(np.where(x<=target_x)[0])
    i2 = min(np.where(x>=target_x)[0])
    x1 = x[i1]
    x2 = x[i2]
    y1 = y[i1]
    y2 = y[i2]
    target_y = y1 + (target_x-x1)*(y2-y1)/(x2-x1)    
    return target_y