from __future__ import division
import sys, time, functools, inspect, types

import numpy as np
from numpy.random import randn
from numpy.lib.stride_tricks import as_strided as ast
import matplotlib.pyplot as plt

##############################
#  numerical linear algebra  #
##############################

def AR_striding(data,nlags):
    # I had some trouble with views and as_strided, so copy if not contiguous
    data = np.asarray(data)
    if not data.flags.c_contiguous:
        data = data.copy(order='C')
    if data.ndim == 1:
        data = np.reshape(data,(-1,1))
    sz = data.dtype.itemsize
    return ast(
            data,
            shape=(data.shape[0]-nlags,data.shape[1]*(nlags+1)),
            strides=(data.shape[1]*sz,sz))

def project_rowspace(A,B):
    return A.dot(B.T.dot(np.linalg.solve(B.dot(B.T),B)))

def project_rowspace_slow(A,B):
    return A.dot(B.T.dot(solve_psd(B.dot(B.T),B)))

def solve_psd(A,b,overwrite_b=False,return_chol=False):
    from scipy.linalg.lapack import dposv
    if return_chol:
        L, x, _ = dposv(A,b,lower=1,overwrite_b=overwrite_b)
        return np.tril(L), x
    else:
        return dposv(A,b,overwrite_b=overwrite_b)[1]

# http://web.stanford.edu/group/mmds/slides2010/Martinsson.pdf
def thin_svd_randomized(A,k):
    n = A.shape[1]
    Omega = randn(n,k)
    Y = A.dot(Omega)
    Q, R = np.linalg.qr(Y)
    B = Q.T.dot(A)
    Uhat, Sigma, V = np.linalg.svd(B)
    U = Q.dot(Uhat)
    return U, Sigma, V.T

def thin_svd(A,k):
    U, s, VT = np.linalg.svd(A)
    return U[:,:k], s[:k], VT[:k]


def normalize(a):
    return a / a.sum()

#########
#  viz  #
#########

def plot_eigvals(mat,*args,**kwargs):
    evals = np.linalg.eigvals(mat)
    print(sorted(evals,key=np.real))
    plt.plot(np.real(evals),np.imag(evals),*args,**kwargs)
    plt.axis([-1,1,-1,1])

def plot_singularvals(mat,*args,**kwargs):
    svals = normalize(np.linalg.svd(mat)[1])
    print(svals)
    plt.plot(svals)

##########
#  text  #
##########

def printnow(s):
    sys.stdout.write(s + '\n')
    sys.stdout.flush()

def print_enter_exit(func):
    global indentlevel
    indentlevel = 0

    @functools.wraps(func)
    def wrapped(*args,**kwargs):
        global indentlevel
        printnow('   '*indentlevel + 'entering %s...' % func.__name__)
        indentlevel += 1

        tic = time.time()
        out = func(*args,**kwargs)
        elapsed = time.time() - tic

        indentlevel -= 1
        printnow('   '*indentlevel + '...done in %0.3f seconds!' % elapsed)

        return out

    return wrapped

# python abuse! see these:
# http://stackoverflow.com/a/1095621
# http://stackoverflow.com/a/8951874
def attach_print_enter_exit():
    frame = inspect.stack()[1]
    module = inspect.getmodule(frame[0])
    for k,v in vars(module).items():
        if isinstance(v, types.FunctionType) and inspect.getmodule(v) == module:
            vars(module)[k] = print_enter_exit(v)
