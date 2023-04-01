import numpy as np
import scipy
linsolve = np.linalg.solve
lsqminnorm = lambda *args: np.linalg.lstsq(*args, rcond=None)[0]
import multiprocessing
from functools import partial

# Y = output data: response to unit impulse, or "impulse response," or "Markov parameters".
# dimensions of Y: p x q x nt, where nt = number of timesteps = number of Markov parameters = number of blocks
# nc = number of block rows in Hankel matrix = order of controllability matrix
# no = number of block columns in Hankel matrix = order of observability matrix
# r = reduced model order = dimension of reduced A = newly assumed dimension of state variable
def era(Y,no,nc,r=None):
    if not r:
        r = int(Y.shape[2]/2)
    p,q = Y.shape[:2] # p = number of outputs, q = number of inputs
    
    # get D from first p x q block of impulse response
    Dr = Y[:,:,0]  # first block of output data
    
    assert Y.shape[2] >= no+nc   # make sure there are enough timesteps to assemble this size of Hankel matrix
    
    # make impulse response into Hankel matrix and shifted Hankel matrix
    H = np.zeros((p*(no), q*(nc+1)))
    for i in range(no):
        for j in range(nc+1):
            H[p*i:p*(i+1), q*j:q*(j+1)] = Y[:,:,i+j+1]
    H0 = H[:,:-q]
    H1 = H[:,q:]
    assert H0.shape == H1.shape == (p*(no), q*(nc))

    # reduced SVD of Hankel matrix
    def _svd(*args):
        U,S,V = scipy.linalg.svd(*args, lapack_driver="gesvd")
        return U,S,V.T.conj()
    U,S,V = _svd(H0)
    SigmaInvSqrt = np.diag(S[:r]**-0.5)
    SigmaSqrt = np.diag(S[:r]**0.5)
    Ur = U[:,:r]
    Vr = V[:,:r]

    # get A from SVD and shifted Hankel matrix
    Ar = SigmaInvSqrt @ Ur.T.conj() @ H1 @ Vr @ SigmaInvSqrt

    # get B and C
    Br = (SigmaSqrt @ Vr.T.conj())[:,:q]
    Cr = (Ur @ SigmaSqrt)[:p,:]

    return (Ar,Br,Cr,Dr)


# a = (alpha) number of block rows in Hankel of correlation matrix
# b = (beta) number of block columns in Hankel of correlation matrix
# l = initial lag for data correlations
# g = lags (gap) between correlation matrices
def era_dc(Y,no,nc,a=0,b=0,l=0,g=1,r=None):
    if not r:
        r = int(Y.shape[2]/2)
    p,q = Y.shape[:2] # p = number of outputs, q = number of inputs

    # get D from first p x q block of impulse response
    Dr = Y[:,:,0]  # first block of output data
    
    assert Y.shape[2] >= l+(a+1+b+1)*g+no+nc   # make sure there are enough timesteps to assemble the Hankel matrices
    
    # Hankel matrix of impulse response (Markov parameters)
    H = np.zeros((p*(no), q*(nc+l+(a+1+b+1)*g)))
    for i in range(no):
        for j in range(nc+l+(a+1+b+1)*g):
            H[p*i:p*(i+1), q*j:q*(j+1)] = Y[:,:,i+j+1]
    H0 = H[:,:q*nc]
    assert H0.shape == (p*(no), q*(nc))

    dimR = p*no # Dimension of square correlation matrices
    dimHRl = (dimR*(a+1), dimR*(b+1)) # Dimension of Hankel matrix of correlation matrices
    HRl = np.zeros(dimHRl) # Hankel matrix of correlation matrices
    HRl1 = np.zeros(dimHRl) # Shifted Hankel matrix of correlation matrices
    for i in range(a+1):
        for j in range(b+1):
            Hl = H[:, q*(l+1+(i+j)*g):q*(l+1+(i+j)*g+nc)]
            Hl1 = H[:, q*(l+1+(i+j)*g+1):q*(l+1+(i+j)*g+nc+1)]
            assert Hl.shape == Hl1.shape == (p*(no), q*(nc))
            R = Hl@H0       # correlation matrix
            R1 = Hl1@H0     # shifted correlation matrix
            assert R.shape == R1.shape == (dimR,dimR)
            HRl[dimR*i:dimR*(i+1), dimR*j:dimR*(j+1)] = R
            HRl1[dimR*i:dimR*(i+1), dimR*j:dimR*(j+1)] = R1

    # reduced SVD of Hankel matrix of correlation matrices
    def _svd(*args):
        U,S,V = scipy.linalg.svd(*args, lapack_driver="gesvd")
        return U,S,V.T.conj()
    U,S,V = _svd(HRl)
    SigmaInvSqrt = np.diag(S[:r]**-0.5)
    SigmaSqrt = np.diag(S[:r]**0.5)
    Ur = U[:,:r]
    Vr = V[:,:r]

    # get A from SVD and shifted Hankel matrix of correlation matrices
    Ar = SigmaInvSqrt @ Ur.T.conj() @ HRl1 @ Vr @ SigmaInvSqrt

    # get B and C
    Br = ((SigmaInvSqrt @ Ur.T.conj())[:,:dimR] @ H0)[:,:q]
    Cr = (Ur @ SigmaSqrt)[:p,:]

    return (Ar,Br,Cr,Dr)


def _blk_3(i, CA, U):
    return i, np.einsum('kil,klj->ij', CA[:i,:,:], U[-i:,:,:])

def srim(
    dati,
    dato,
    debug = False,
    full     : bool = True,
    verbose  : bool = False,
    pool_size: int  = 6,
    **config
):
    """
    mro $(p)$ determines order of the observer Kalman ARX filter used in OKID-ERA-DC.
    orm $(n)$ determines size of the state-space model used for representing the system.

    Returns
    =======
     freqdampSRIM:
        variable is a matrix that includes the information of identified
        frequencies, damping ratios & validation of the modes with MPC & EMAC criteria.
        Each row of freqdamp corresponds to a mode. Columns are as follows:
        1) frequency, 2) damping ratio, 3) order index, 4) condition number, 5) MPC.
        If values in columns 5 is > 0.5, identified mode is valid.
     modeshapeSRIM:
         stores the mode shape information for identified modes.
     RMSEpredSRIM:
         root mean square error of the predicted output from
         identified parameters with respect to the actual output

    ## SRIM Methodology

    More information on SRIM algorithm can be found in Sections 3.4.4 & 3.4.5 of
    (Arici & Mosalam, 2006). Equations below refer to this report. SRIM is a MIMO
    SI method that is based on state space identification using least squares and
    consists of the following steps:


    2a. Determine output (y) & input (u) vectors [Eqs. 3.58 & 3.60].
    2b. Compute the correlation terms & the coefficient matrix (Eqs. 3.68 & 3.69).
    2c. Obtain observability matrix using full or partial decomposition (Eqs. 3.72 & 3.74).
    2d. Use the observability matrix to compute system matrices A, B & C, in which modal 
        information is embedded.
    2e. Obtain the modal information from matrices A, B & C.
    2f. Spatial & temporal validation of the identified modes.
    2g. Back calculate (estimate) the output accelerations with the state-space system &
        check against the actual output accelerations.

    For orm = 2, one mode is found, for orm = 4, two modes are found.
    For case 1, one mode is transverse & the other is torsion.
    For all other cases, the second mode is a higher mode.
    Sometimes higher orm still gives fewer modes, e.g. orm = 8 for case 1 gives
    three modes, but one of them is invalid according to the EMAC & MPC criteria.
    same orm in OKID-ERA-DC is used. It can be changed if needed.

    """
    
    #
    # Convenience argument handling
    #
    p = config.get("p", config.get("mro"))         # # steps used for the identification (ie, prediction horizon)
    n = n1 = config.get("n", config.get("orm", 4))  # Order of the model.

    if isinstance(dati, list):
        dati = np.array([i.data for i in dati]).T
    elif issubclass(dati.__class__, dict):
        dati = dati.data

    if isinstance(dato, list):
        dato = np.array([i.data for i in dato]).T
    elif issubclass(dato.__class__, dict):
        dato = dato.data

    if len(dati.shape) < 2:
        dati = np.atleast_2d(dati).T
    if len(dato.shape) < 2:
        dato = np.atleast_2d(dato).T

    if verbose:
        progress_bar = tqdm
    else:
        progress_bar = lambda arg, **kwds: (i for i in arg)

    dn = config.get("dn", None) or dati.shape[0]

    # 2a. Compute y (output) and u (input) vectors (Eqs. 3.58 & 3.60)

    # Note that main Step 2 develops Eq. 3.57.
    # Therefore, it is not part of the code.
    # Accordingly, the code continues with Step 2a to compute the output & input vectors.

    # Calculate the usable size of the data matrix
    # dn = size(dat,1)/div;       # total # time steps after decimating
    nsizS = dn-1-p+2

    l,m = dato.shape # m is the number of output channels
    _,r = dati.shape # r is the number of input channels


    assert p >= n/m + 1

    ypS = np.zeros((m*p,nsizS))
    upS = np.zeros((r*p,nsizS))

    # Compute Y (output) & U (input) vectors (Eqs. 3.58 & 3.60 Arici 2006)
    for b in range(p):
        ypS[b*m:(b+1)*m,:nsizS+1] = dato[b:nsizS+b, :].T
        upS[b*r:(b+1)*r,:nsizS+1] = dati[b:nsizS+b, :].T


    # 2b. Compute the correlation terms and the coefficient matrix 
    #     (Eqs. 3.68 & 3.69).

    # Compute the correlation terms (Eq. 3.68)
    Ryy = ypS@ypS.T/nsizS
    Ruu = upS@upS.T/nsizS
    Ruy = upS@ypS.T/nsizS

    assert Ryy.shape[0] == Ryy.shape[1] == p*m
    assert Ruy.shape[0] == p*r
    assert Ruy.shape[1] == p*m

    # Compute the correlation matrix (Eq. 3.69)
    Rhh = Ryy - Ruy.T@linsolve(Ruu,Ruy)

    # 2c. Obtain observability matrix using full or partial decomposition 
    #     (Eqs. 3.72 & 3.74).
    if full:
        # Full Decomposition Method
        un,*_ = np.linalg.svd(Rhh,0)           # Eq. 3.74
        Observability = un[:,:n]                          # Eq. 3.72
        A = lsqminnorm(Observability[:(p-1)*m,:], Observability[m:p*m,:])
        C = Observability[:m,:]
    else:
        # Partial Decomposition Method
        un,*_ = np.linalg.svd(Rhh[:,:(p-1)*m+1],0)
        Observability = un[:,:n]
        A = lsqminnorm(Observability[:(p-1)*m,:], Observability[m:p*m,:])
        C = un[:m,:]

    # Computation of system matrices B & D
    # Output Error Minimization

    # Setting up the Phi matrix
    Phi  = np.zeros((m*nsizS, n+m*r+n*r))
    CA_powers = np.zeros((nsizS, m, A.shape[1]))
    CA_powers[0, :, :] = C
    A_p = A
    for pwr in range(1,nsizS):
        CA_powers[pwr,:,:] =  C@A_p
        A_p = A@A_p

    # First block column of Phi
    for df in range(nsizS):
        Phi[df*m:(df+1)*m,:n] = CA_powers[df,:,:]

    # Second block column of Phi
    Imm = np.eye(m)
    for i in range(nsizS):
        Phi[i*m:(i+1)*m, n:n+m*r] = np.kron(dati[i,:],Imm)

    # Third block column of Phi
    In1n1 = np.eye(n)
    cc = n + m*r + 1
    dd = n + m*r + n*r

    krn = np.array([np.kron(dati[i,:],In1n1) for i in range(nsizS)])

    with multiprocessing.Pool(pool_size) as pool:
        for i,res in progress_bar(
                pool.imap_unordered(
                    partial(_blk_3,CA=CA_powers,U=np.flip(krn,0)),
                    range(1,nsizS),
                    200
                ),
                total = nsizS
            ):
            Phi[i*m:(i+1)*m,cc-1:dd] = res

    y = dato[:nsizS,:].flatten()

    teta = lsqminnorm(Phi,y)

    x0 = teta[:n1]
    dcol = teta[n1:n1+m*r]
    bcol = teta[n1+m*r:n1+m*r+n1*r]

    D = np.zeros((m,r))
    B = np.zeros((n,r))
    for wq in range(r):
        D[:,wq] = dcol[wq*m:(wq+1)*m]

    for ww in range(r):
        B[:,ww] = bcol[ww*n:(ww+1)*n]

    assert A.shape[0] == A.shape[1] == n
    if debug:
        return locals()
    return A,B,C,D
