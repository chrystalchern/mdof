# Adapted from:
#   https://github.com/mattjj/py4sid/raw/master/estimation.py
#
import numpy as np
import scipy
import matplotlib.pyplot as plt

import sys, time, functools, inspect, types

from numpy.random import randn
from numpy.lib.stride_tricks import as_strided

from mdof.numerics import thin_svd, solve_psd

##############################
#  numerical linear algebra  #
##############################

def _AR_striding(data,nlags):
    # I had some trouble with views and as_strided, so copy if not contiguous
    data = np.asarray(data)
    if not data.flags.c_contiguous:
        data = data.copy(order='C')
    if data.ndim == 1:
        data = np.reshape(data,(-1,1))

    sz = data.dtype.itemsize
    return as_strided(
            data,
            shape=(data.shape[0]-nlags,data.shape[1]*(nlags+1)),
            strides=(data.shape[1]*sz,sz))



def _project_rowspace(A,B, slow=False):
    if slow:
        return A.dot(B.T.dot(solve_psd(B.dot(B.T),B)))
    else:
        return A.dot(B.T.dot(np.linalg.solve(B.dot(B.T),B)))



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


####################
#  Alg 3 in Ch. 3  #
####################

def _compute_Ys(y,i):
    p = y.shape[1]
    Y = _AR_striding(y,2*i-1).T
    Yp, Yf = Y[:i*p].copy(), Y[i*p:].copy()
    Ypp, Yfm = Y[:(i+1)*p].copy(), Y[(i+1)*p:].copy()
    return (Yp, Yf), (Ypp, Yfm)

def _compute_Os(Y,Ym):
    (Yp,Yf) = Y
    (Ypp,Yfm) = Ym
    Oi = _project_rowspace(Yf,Yp)
    Oim1 = _project_rowspace(Yfm,Ypp)
    return Oi, Oim1


def _compute_gammas(Oi,nhat,p):
    U,s,VT = thin_svd(Oi,nhat, random=True)
    Gamma_i = U.dot(np.diag(np.sqrt(s)))
    Gamma_im1 = Gamma_i[:-p]
    return Gamma_i, Gamma_im1

def _compute_Xs(Gamma_i,Gamma_im1,Oi,Oim1):
    Xihat = solve_psd(Gamma_i.T.dot(Gamma_i), Gamma_i.T.dot(Oi))
    Xip1hat = solve_psd(Gamma_im1.T.dot(Gamma_im1), Gamma_im1.T.dot(Oim1))
    return Xihat, Xip1hat


def system_parameters_from_states(Xihat, Xip1hat, Yii, nhat):
    ACT = np.linalg.lstsq(Xihat.T, np.vstack((Xip1hat, Yii)).T)[0]
    Ahat, Chat = ACT.T[:nhat], ACT.T[nhat:]
    rho = np.vstack((Xip1hat,Yii)) - ACT.T.dot(Xihat)

    QSR = rho.dot(rho.T) / rho.shape[1]
    Q, S, R = QSR[:nhat,:nhat], QSR[:nhat,nhat:], QSR[nhat:, nhat:]
    Bhat = scipy.linalg.sqrtm(Q)
    Dhat = scipy.linalg.sqrtm(R)

    return Ahat, Bhat, Chat, Dhat


def estimate_parameters_4sid(y,i,nhat=None,return_xs=False):
    p = y.shape[1]

    (Yp, Yf), (Ypp, Yfm) = _compute_Ys(y,i)
    Oi, Oim1 = _compute_Os((Yp,Yf),(Ypp,Yfm))

    if nhat is None:
        viz_rank2(y,i,k=10)
        plt.show()
        nhat = int(input('n = '))

    Gamma_i, Gamma_im1 = _compute_gammas(Oi,nhat,p)
    Xihat, Xip1hat     = _compute_Xs(Gamma_i,Gamma_im1,Oi,Oim1)

    Yii = Yf[:p]
    Ahat, Bhat, Chat, Dhat = system_parameters_from_states(Xihat, Xip1hat, Yii, nhat)

    if return_xs:
        return (Ahat, Bhat, Chat, Dhat), Xihat.T
    else:
        return Ahat, Bhat, Chat, Dhat


def estimate_parameters_moments(outputs, i, nhat):
    ####################
    #  Alg 2 in Ch. 3  #
    ####################

    y = outputs
    p = y.shape[1]

    ## my version, fast but does it work as well?
    Y = _AR_striding(y,2*i-1).T
    l = Y.shape[0]//2
    U, s, VT = thin_svd(np.cov(Y)[:l,l:],nhat)

    ## slow version in book
    # (Yp, Yf), _ = _compute_Ys(y,i)
    # Oi = _project_rowspace(Yf,Yp)
    # U,s,VT = thin_svd(Oi,nhat)

    Gamma_i = U.dot(np.diag(np.sqrt(s)))

    ## estimate C
    C = Gamma_i[:p]

    ## estimate stable A (other methods on p. 54)
    A = np.linalg.lstsq(Gamma_i,np.vstack((Gamma_i[p:],np.zeros((p,Gamma_i.shape[1])))))[0]

    # TODO compute B, D
    # TODO why is C's first singular value so big? maybe need bigger i here

    return A, C


#
#
#
def viz_rank(Oi,k=None):
    if k is None:
        U, s, VT = np.linalg.svd(Oi)
    else:
        U, s, VT = thin_svd(Oi,k, random=True)

    plt.figure()
    plt.stem(np.arange(s.shape[0]),s)

def viz_rank2(y,i,k=None):
    Y = _AR_striding(y,i-1).T
    l = Y.shape[0]//2

    if k is None:
        _, s, _ = np.linalg.svd(np.cov(Y)[:l,l:])
    else:
        _, s, _ = thin_svd(np.cov(Y)[:l,l:],k, random=True)

    plt.figure()
    plt.stem(np.arange(s.shape[0]),s)


