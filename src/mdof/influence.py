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


def _minoe_cmpcc(inputs, outputs, A, C):
    pass


def _minoe_socit(inputs, outputs, A, C):
    pass


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


def _form_IU(dimI, u):
    assert len(u.shape) == 1
    q = u.shape[0]
    Idd = np.eye(dimI)
    IU = np.empty((dimI,dimI*q))
    for i in range(q):
        IU[:,i*dimI:(i+1)*dimI] = u[i]*Idd
    return IU


def _blk_3(i, CA, U):
    return i, np.einsum('kil,klj->ij', CA[:i,:,:], U[-i:,:,:])


def _form_powers(C,A,ns,p,r,pow_impl,p_max=10):
    CA_powers = np.zeros((ns, p, r), dtype=complex)
    if pow_impl==1:
        CA_powers[0] = C
        A_p = A
        for pwr in range(1,ns):
            CA_powers[pwr] = C@A_p
            A_p = A@A_p
        return CA_powers
    if pow_impl==2:
        Lam, Vl = scipy.linalg.eig(A)
        Vr = scipy.linalg.inv(Vl)
        CA_powers[0] = C
        Lam_p = Lam
        for pwr in range(1,ns):
            CA_powers[pwr] =  C@Vl@np.diag(Lam_p)@Vr
            Lam_p = Lam_p*Lam
        return CA_powers
    if pow_impl==3:
        Lam, Vl = scipy.linalg.eig(A)
        Vr = scipy.linalg.inv(Vl)
        CA_powers[0] = C
        Lam_p = Lam
        for pwr in range(1,ns):
            CA_powers[pwr] =  C@Vl@np.diag(Lam_p)@Vr
            Lam_p = [w**min(pwr, p_max) for w in Lam]
        return CA_powers
    

from .simulate import simulate
def _ac2bd(u, y, A, C):
    """
    Compute the system matrices ``B`` and ``D`` from the system matrices ``A`` and ``C``.
    """
    N, p = y.shape
    _, q = u.shape
    n, _ = A.shape

    # Form observability matrix up to the data length for computing x(0)
    Phi = C
    temp = C
    for k in range(2, N + 1):
        temp = np.dot(temp, A)
        Phi = np.vstack((Phi, temp))

    # Add the input matrix associated with D matrix
    Phi = np.hstack((Phi, np.zeros((N * p, p * q))))
    for j in range(q):
        for i in range(p):
            Phi[i:N*p+1:p, n+j*p+i] = u[:,j]

    # Add the imaginary output matrix associated with B matrix
    null_feedthrough = np.zeros((p, 1))
    for j in range(q):
        for i in range(n):
            ei = np.zeros((n, 1))
            ei[i,0] = 1
            _,z,_ = scipy.signal.dlsim(scipy.signal.dlti(A,ei,C,null_feedthrough),
                                         u[:,j])
            # z = simulate((A,ei,C,null_feedthrough),u[:,j][:,None].T).T
            z = _block_tr(1,N,p,z,0)
            Phi = np.hstack((Phi,z.T))
    ycol = _block_tr(1,N,p,y,0)
    Theta = np.linalg.pinv(Phi) @ ycol.T
    x0 = Theta[:n, 0]
    D = _block_tr(p,q,1,Theta[n:n+p*q],0)
    B = _block_tr(n,q,1,Theta[n+p*q:n+p*q+n*q],0)

    return B,D,x0


def ac2bd(inputs, outputs, A, C, **options):
    """
    Compute the ``B`` and ``D`` system matrices from the ``inputs`` and
    ``outputs`` arrays and the ``A`` and ``C`` system matrices.
    """
    r = options.get('r')

    p,N = outputs.shape
    q,_ = inputs.shape
    
    threads = options.get("threads",6)
    chunk = options.get("chunk", 200)
    
    lsq_solve = numerics.lsq_solver(options.get("lsq", {}))
    
    #
    # Setting up the Phi matrix
    #
    Phi  = np.zeros((p*N, r+p*q+r*q), dtype=complex)

    # First block column of Phi
    CA_powers = _form_powers(C,A,N,p,r,pow_impl=1)

    for df in range(N):
        Phi[df*p:(df+1)*p, :r] = CA_powers[df]

    # Second block column of Phi
    for i in range(N):
        Phi[i*p:(i+1)*p, r:r+p*q] = _form_IU(p,inputs[:,i])

    # Third block column of Phi
    Un = np.array([_form_IU(r,inputs[:,i]) for i in range(N)])
    assert Un.shape == (N,r,r*q)

    # Execute a loop in parallel that looks something like:
    #    for i in  range(1,N):
    #        Phi[] = _blk_3(i, CA_Powers, np.flip(...))

    # for i in range(1,N):
    #     Phi[i*p:(i+1)*p, r+p*q:] = np.sum([CA_powers[j]@Un[i-j-1] for j in range(i)],0)

    # for i in progress_bar(range(1,N), total=N-1):
    #     Phi[i*p:(i+1)*p, r+p*q:] = sum(  CA_powers[j]@Un[i-j-1] for j in range(i)   )

    with multiprocessing.Pool(threads) as pool:
        for i,res in progress_bar(
                pool.imap_unordered(
                    partial(_blk_3, CA=CA_powers,U=np.flip(Un,0)),
                    range(1,N),
                    chunk
                ),
                total = N
            ):
            Phi[i*p:(i+1)*p, r+p*q:] = res

    # zero_q = np.zeros((q, 1))
    # for j in range(q):
    #     for i in range(r):
    #         ei_r = np.zeros((r,1))
    #         ei_r[i] = 1
    #         _, z, _ = scipy.signal.dlsim(scipy.signal.dlti(A, ei_r, C, zero_q), inputs[j,:])
    #         assert z.shape == (N,p)
    #         z = z.flatten()
    #         Phi[:,r+p*q+j*i] = z
            
    y = outputs.T.flatten()
    # y = outputs.flatten()
    Theta = lsq_solve(Phi,y)
    # teta = np.linalg.pinv(Phi) @ y

    x0 = Theta[:r]
    dcol = Theta[r:r+p*q]
    bcol = Theta[r+p*q:r+p*q+r*q]

    D = np.zeros((p,q), dtype=complex)
    B = np.zeros((r,q), dtype=complex)
    for i in range(q):
        D[:,i] = dcol[i*p:(i+1)*p]
    if np.max(D.imag) > 1e-4:
        warnings.warn(f"matrix D imaginary part as large as {np.max(D.imag)}")
    D = D.real
    for i in range(q):
        B[:,i] = bcol[i*r:(i+1)*r]
    if np.max(B.imag) > 1e-4:
        warnings.warn(f"matrix B imaginary part as large as {np.max(B.imag)}")
    B = B.real

    return B,D