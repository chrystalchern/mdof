import numpy as np
lin_solve = np.linalg.solve
import multiprocessing
from functools import partial
import warnings
try:
    from tqdm import tqdm as progress_bar

except:
    def progress_bar(arg, **kwds): return arg

from . import numerics

import scipy

import numpy as np


def _block_tr(nrow, nblock, ncol, a, flag: int):
    """
    Routine to perform block transposition

    Args:
        nrow (int): Number of rows in each block
        nblock (int): Number of blocks
        ncol (int): Number of columns in each block
        a (numpy.ndarray): Input matrix
        flag (int): Set to 1 for wide matrix, 0 for tall matrix

    Returns:
        numpy.ndarray: Block transposition (not the matrix transposed)
    """
    at = []
    if flag == 1:
        for i in range(nblock):
            nstart = i*ncol
            nend = i * ncol
            at.append(a[:, i*ncol:(i+1)*ncol])
        at = np.concatenate(at, axis=0)
    else:
        for i in range(nblock):
            nstart = i*nrow
            nend = i * nrow
            at.append(a[i*nrow:(i+1)*nrow, :])
        at = np.concatenate(at, axis=1)
    return at


def _minoe_cmpcc(inputs, outputs, A, C):
    pass

def _minoe_socit(inputs, outputs, A, C):
    pass

def ac2bd(inputs, outputs, A, C, option=None):
    pass
