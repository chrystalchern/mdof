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


def form_IU(dimI, u):
    pass


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
    

def _ac2bd(u, y, a, c):
    """
    Compute the system matrices ``B`` and ``D`` from the system matrices ``A`` and ``C``.
    """
    nd, m = y.shape
    _, r = u.shape
    n, _ = a.shape

    # Form observability matrix up to the data length for computing x(0)
    phi = c
    temp = c
    for k in range(2, nd + 1):
        temp = np.dot(temp, a)
        phi = np.vstack((phi, temp))

    # Add the input matrix associated with D matrix
    phi = np.hstack((phi, np.zeros((nd * m, m * r))))
    for j in range(r):
        for i in range(m):
            phi[i:nd*m+1:m, n+j*m+i] = u[:,j]

    # Add the imaginary output matrix associated with B matrix
    d = np.zeros((m, 1))
    for j in range(r):
        for i in range(n):
            b = np.zeros((n, 1))
            b[i, 0] = 1
            _, z, _ = scipy.signal.dlsim(scipy.signal.dlti(a, b, c, d), u[:, j])
            z = _block_tr(1, nd, m, z, 0)
            phi = np.hstack((phi, z.T))
    z = _block_tr(1, nd, m, y, 0)
    b = np.linalg.pinv(phi) @ z.T
    x0 = b[:n, 0]
    d = _block_tr(m, r, 1, b[n:n+m*r], 0)
    b = _block_tr(n, r, 1, b[n+m*r:n+m*r+n*r], 0)

    return b, d, x0


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
    
    # Setting up the Phi matrix
    Phi  = np.zeros((p*N, r+p*q+r*q), dtype=complex)

    # First block column of Phi
    CA_powers = _form_powers(C,A,N,p,r,pow_impl=1)

    for df in range(N):
        Phi[df*p:(df+1)*p, :r] = CA_powers[df]

    # Second block column of Phi
    Ipp = np.eye(p)
    for i in range(N):
        Phi[i*p:(i+1)*p, r:r+p*q] = np.kron(Ipp,inputs[:,i])

    # Third block column of Phi
    Irr = np.eye(r)
    Un = np.array([np.kron(Irr,inputs[:,i]) for i in range(N)])
    assert Un.shape == (N,r,r*q)

    # Execute a loop in parallel that looks something like:
    #    for i in  range(1,N):
    #        Phi[] = _blk_3(i, CA_Powers, np.flip(...))

    # for i in range(1,N):
    #     Phi[i*p:(i+1)*p, r+p*q:] = np.sum([CA_powers[j]@Un[i-j-1] for j in range(i)],0)

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

    y = outputs[:,:N].flatten()
    teta = lsq_solve(Phi,y)
    # teta = np.linalg.pinv(Phi) @ y

    x0 = teta[:r]
    dcol = teta[r:r+p*q]
    bcol = teta[r+p*q:r+p*q+r*q]

    D = np.zeros((p,q), dtype=complex)
    B = np.zeros((r,q), dtype=complex)
    for wq in range(q):
        D[:,wq] = dcol[wq*p:(wq+1)*p]
    if np.max(D.imag) > 1e-4:
        warnings.warn(f"matrix D imaginary part as large as {np.max(D.imag)}")
    D = D.real
    for ww in range(q):
        B[:,ww] = bcol[ww*r:(ww+1)*r]
    if np.max(B.imag) > 1e-4:
        warnings.warn(f"matrix B imaginary part as large as {np.max(B.imag)}")
    B = B.real

    return B,D