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


def _form_IU(dimI, u):
    r"""
    Broadcast the scalar elements of the vector `u` as :math:`[Iu_1, Iu_2, \cdots, Iu_q]`.
    Same as:
    q = u.shape[0]
    Idd = np.eye(dimI)
    IU = np.empty((dimI,dimI*q))
    for i in range(q):
        IU[:,i*dimI:(i+1)*dimI] = u[i]*Idd
    """
    assert u.ndim == 1
    Idd = np.eye(dimI)
    IU = np.kron(u,Idd)
    return IU


def _blk_3(i, CA, U):
    return i, np.einsum('kil,klj->ij', CA[:i,:,:], U[-i:,:,:])

from .construct import stack_powers
from .simulate import simulate
def ac2bd(inputs, outputs, A, C, **options):
    """
    Compute the ``B`` and ``D`` system matrices from the ``inputs`` and
    ``outputs`` arrays and the ``A`` and ``C`` system matrices.
    """
    n, _ = A.shape

    p,N = outputs.shape
    q,_ = inputs.shape

    
    lsq_solve = numerics.lsq_solver(options.get("lsq", {}))
    
    #
    # Setting up the Phi matrix
    #
    Phi  = np.zeros((p*N, n+p*q+n*q), dtype=complex)

    # First block column of Phi 
    CA_powers = C@stack_powers(A, n_pwr=N)
    # for i in range(N):
    #     Phi[i*p:(i+1)*p, :n] = CA_powers[i]
    Phi[:,:n] = CA_powers.reshape(N*p,n)

    # Second block column of Phi
    for i in range(N):
        Phi[i*p:(i+1)*p, n:n+p*q] = _form_IU(p,inputs[:,i])

    # Third block column of Phi
    Un = np.array([_form_IU(n,inputs[:,i]) for i in range(N)])
    assert Un.shape == (N,n,n*q)
    if False:
        threads = options.get("threads",6)
        chunk = options.get("chunk", 200)
        # Execute a loop in parallel that looks something like:
        #
        # for i in range(1,N):
        #     Phi[i*p:(i+1)*p, r+p*q:] = sum(  CA_powers[j]@Un[i-j-1] for j in range(i)   )
        #
        with multiprocessing.Pool(threads) as pool:
            for i,res in progress_bar(
                    pool.imap_unordered(
                        partial(_blk_3, CA=CA_powers,U=np.flip(Un,0)),
                        range(1,N),
                        chunk
                    ),
                    total = N
                ):
                Phi[i*p:(i+1)*p, n+p*q:] = res
    else:
        zero_q = np.zeros((p, 1))
        ei = np.zeros((n,1))
        start = n+p*q
        offset = 0
        for j in range(q):
            for i in range(n):
                ei[:] = 0
                ei[i] = 1
                # _, z, _ = scipy.signal.dlsim(scipy.signal.dlti(A, ei, C, zero_q), inputs[j,:])
                z = simulate((A,ei,C,zero_q),inputs[j,:][:,None].T).T
                assert z.shape == (N,p)
                z = z.flatten()
                Phi[:,start+offset] = z
                offset += 1
    
    # Solve the least squares problem :
    #    Theta = Phi \ y
    y = outputs.T.flatten()
    Theta = lsq_solve(Phi,y)

    dcol = Theta[n:n+p*q]
    bcol = Theta[n+p*q:n+p*q+n*q]
    D = dcol.reshape(p,q)
    B = bcol.reshape(q,n).T
    if np.max(D.imag) > 1e-4:
        warnings.warn(f"matrix D imaginary part as large as {np.max(D.imag)}")
    if np.max(B.imag) > 1e-4:
        warnings.warn(f"matrix B imaginary part as large as {np.max(B.imag)}")
    D = D.real
    B = B.real

    return (B,D)