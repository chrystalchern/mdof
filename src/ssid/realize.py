import numpy as np
import scipy
linsolve = np.linalg.solve
lsqminnorm = lambda *args: np.linalg.lstsq(*args, rcond=None)[0]
import multiprocessing
from functools import partial

# Y = output data: response to unit impulse, or "impulse response," or "Markov parameters".
# dimensions of Y: p x q x nt, where nt = number of timesteps = number of Markov parameters = number of blocks
# no = number of block rows in Hankel matrix = order of observability matrix
# nc = number of block columns in Hankel matrix = order of controllability matrix
# r = reduced model order = dimension of reduced A = newly assumed dimension of state variable
def era(Y,no=None,nc=None,r=None,**options):
    
    if r is None:
        r = min(20, int(Y.shape[2]/2))
    p,q = Y.shape[:2] # p = number of outputs, q = number of inputs
    
    # get D from first p x q block of impulse response
    Dr = Y[:,:,0]  # first block of output data
    
    # size of Hankel matrix
    if no is None:
        if nc is None:
            no = nc = min(300, (Y.shape[2]-1)/2)
        else:
            no = min(300, Y.shape[2]-1-nc)
    elif nc is None:
        nc = min(300, Y.shape[2]-1-no)
    else:
        # make sure there are enough timesteps to assemble this size of Hankel matrix
        assert Y.shape[2] >= no+nc
    
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
def era_dc(Y,no=None,nc=None,a=0,b=0,l=0,g=1,r=None,**options):
    if r is None:
        r = int(Y.shape[2]/2)
    p,q = Y.shape[:2] # p = number of outputs, q = number of inputs

    # get D from first p x q block of impulse response
    Dr = Y[:,:,0]  # first block of output data

    # size of Hankel matrix
    if no is None:
        if nc is None:
            no = nc = min(300, int((Y.shape[2]-1)/2))
        else:
            no = min(300, int(Y.shape[2]-1-nc))
    elif nc is None:
        nc = min(300, int(Y.shape[2]-1-no))
    else:
        # make sure there are enough timesteps to assemble the Hankel matrices
        assert Y.shape[2] >= l+(a+1+b+1)*g+no+nc
    
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
    input,
    output,
    no = None,
    r = None,
    debug = False,
    full     : bool = True,
    pool_size: int  = 6,
    **config
):
    """
    no: number of steps used for identification (prediction horizon),
        and the order of the observability matrix.
        analagous to order (number of autoregressors) of the observer Kalman
        ARX filter used in OKID-ERA-DC.
    r:  size of the state-space model used for representing the system.

    Returns
    =======
    A,B,C,D:
        State-space model matrices

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

    For r = 2, one mode is found, for r = 4, two modes are found.
    For case 1, one mode is transverse & the other is torsion.
    For all other cases, the second mode is a higher mode.
    Sometimes higher r still gives fewer modes, e.g. r = 8 for case 1 gives
    three modes, but one of them is invalid according to the EMAC & MPC criteria.
    same r in OKID-ERA-DC is used. It can be changed if needed.

    """
    
    #
    # Convenience argument handling
    #

    if isinstance(input, list):
        input = np.array([i.data for i in input]).T
    elif issubclass(input.__class__, dict):
        input = input.data

    if isinstance(output, list):
        output = np.array([i.data for i in output]).T
    elif issubclass(output.__class__, dict):
        output = output.data

    if len(input.shape) < 2:
        input = np.atleast_2d(input).T
    if len(output.shape) < 2:
        output = np.atleast_2d(output).T

    progress_bar = lambda arg, **kwds: (i for i in arg)

    assert input.shape[0] == output.shape[0]
    nt = config.get("nt", None) or input.shape[0]

    if no is None:
        no = min(300, nt)

    if r is None:
        r = min(10, int(no/2))
    # 2a. Compute y (output) and u (input) vectors (Eqs. 3.58 & 3.60)

    # Note that main Step 2 develops Eq. 3.57.
    # Therefore, it is not part of the code.
    # Accordingly, the code continues with Step 2a to compute the output & input vectors.

    # Calculate the usable size of the data matrix
    # nt = size(dat,1)/div;       # total # time steps after decimating
    ns = nt-1-no+2

    _,p = output.shape # p is the number of output channels
    _,q = input.shape # q is the number of input channels

    assert no >= r/p + 1    # make sure prediction horizon is large enough that
                            # observability matrix is full rank

    ypS = np.zeros((p*no,ns))
    upS = np.zeros((q*no,ns))

    # Compute Y (output) & U (input) vectors (Eqs. 3.58 & 3.60 Arici 2006)
    for b in range(no):
        ypS[b*p:(b+1)*p,:ns+1] = output[b:ns+b, :].T
        upS[b*q:(b+1)*q,:ns+1] = input[b:ns+b, :].T


    # 2b. Compute the correlation terms and the coefficient matrix 
    #     (Eqs. 3.68 & 3.69).

    # Compute the correlation terms (Eq. 3.68)
    Ryy = ypS@ypS.T/ns
    Ruu = upS@upS.T/ns
    Ruy = upS@ypS.T/ns

    assert Ryy.shape[0] == Ryy.shape[1] == no*p
    assert Ruy.shape[0] == no*q
    assert Ruy.shape[1] == no*p

    # Compute the correlation matrix (Eq. 3.69)
    Rhh = Ryy - Ruy.T@linsolve(Ruu,Ruy)

    # 2c. Obtain observability matrix using full or partial decomposition 
    #     (Eqs. 3.72 & 3.74).
    if full:
        # Full Decomposition Method
        un,*_ = np.linalg.svd(Rhh,0)           # Eq. 3.74
        Observability = un[:,:r]                          # Eq. 3.72
        A = lsqminnorm(Observability[:(no-1)*p,:], Observability[p:no*p,:])
        C = Observability[:p,:]
    else:
        # Partial Decomposition Method
        un,*_ = np.linalg.svd(Rhh[:,:(no-1)*p+1],0)
        Observability = un[:,:r]
        A = lsqminnorm(Observability[:(no-1)*p,:], Observability[p:no*p,:])
        C = un[:p,:]

    # Computation of system matrices B & D
    # Output Error Minimization

    # Setting up the Phi matrix
    Phi  = np.zeros((p*ns, r+p*q+r*q))
    CA_powers = np.zeros((ns, p, A.shape[1]))
    CA_powers[0, :, :] = C
    A_p = A
    for pwr in range(1,ns):
        CA_powers[pwr,:,:] =  C@A_p
        A_p = A@A_p

    # First block column of Phi
    for df in range(ns):
        Phi[df*p:(df+1)*p,:r] = CA_powers[df,:,:]

    # Second block column of Phi
    Ipp = np.eye(p)
    for i in range(ns):
        Phi[i*p:(i+1)*p, r:r+p*q] = np.kron(input[i,:],Ipp)

    # Third block column of Phi
    In1n1 = np.eye(r)
    cc = r + p*q + 1
    dd = r + p*q + r*q

    krn = np.array([np.kron(input[i,:],In1n1) for i in range(ns)])

    with multiprocessing.Pool(pool_size) as pool:
        for i,res in progress_bar(
                pool.imap_unordered(
                    partial(_blk_3,CA=CA_powers,U=np.flip(krn,0)),
                    range(1,ns),
                    200
                ),
                total = ns
            ):
            Phi[i*p:(i+1)*p,cc-1:dd] = res

    y = output[:ns,:].flatten()

    teta = lsqminnorm(Phi,y)

    x0 = teta[:r]
    dcol = teta[r:r+p*q]
    bcol = teta[r+p*q:r+p*q+r*q]

    D = np.zeros((p,q))
    B = np.zeros((r,q))
    for wq in range(q):
        D[:,wq] = dcol[wq*p:(wq+1)*p]

    for ww in range(q):
        B[:,ww] = bcol[ww*r:(ww+1)*r]

    assert A.shape[0] == A.shape[1] == r
    if debug:
        return locals()
    return A,B,C,D

# input = input data. dimensions of input: q x nt, where nt = number of timesteps.
# output = response data due to input data. dimensions of output: p x nt.
# no = number of steps used for identification (prediction horizon),
    # and the order of the observability matrix.
    # analagous to order (number of autoregressors) of the observer Kalman
    # ARX filter used in OKID-ERA-DC.
# r = size of the state-space model used for representing the system.
# Juang 1997, "System Realization Using Information Matrix," Journal of Guidance, Control, and Dynamics
def srim2(input,output,no=None,r=None,full=True,**options):
    input = np.atleast_2d(input)
    output = np.atleast_2d(output)

    p = output.shape[0]
    q = input.shape[0]
    nt = output.shape[1]
    assert output.shape[1] == input.shape[1]

    if no is None:
        no = min(300, nt)

    if r is None:
        r = min(10, int(no/2))

    # maximum possible number of columns in the Y and U data matrices
    ns = nt-1-no+2

    assert no >= r/p + 1    # make sure prediction horizon is large enough that
                            # observability matrix is full rank (Juang Eq. 8)

    Yno = np.zeros((p*no,ns))
    Uno = np.zeros((q*no,ns))

    # Construct Y (output) & U (input) data matrices (Eqs. 3.58 & 3.60 Arici 2006)
    for i in range(no):
        Yno[i*p:(i+1)*p,:] = output[:,i:ns+i]
        Uno[i*q:(i+1)*q,:] = input[:,i:ns+i]


    # 2b. Compute the correlation terms and the coefficient matrix 
    #     (Eqs. 3.68 & 3.69).

    # Compute the correlation terms (Eq. 3.68)
    Ryy = Yno@Yno.T/ns
    Ruu = Uno@Uno.T/ns
    Ruy = Uno@Yno.T/ns

    assert Ryy.shape[0] == Ryy.shape[1] == no*p
    assert Ruy.shape[0] == no*q
    assert Ruy.shape[1] == no*p

    # Compute the correlation matrix (Eq. 3.69)
    Rhh = Ryy - Ruy.T@linsolve(Ruu,Ruy)

    # 2c. Obtain observability matrix using full or partial decomposition 
    #     (Eqs. 3.72 & 3.74).
    if full:
        # Full Decomposition Method
        un,*_ = np.linalg.svd(Rhh,0)           # Eq. 3.74
        Observability = un[:,:r]               # Eq. 3.72
        A = lsqminnorm(Observability[:(no-1)*p,:], Observability[p:no*p,:])
        C = Observability[:p,:]
    else:
        # Partial Decomposition Method
        un,*_ = np.linalg.svd(Rhh[:,:(no-1)*p+1],0)
        Observability = un[:,:r]
        A = lsqminnorm(Observability[:(no-1)*p,:], Observability[p:no*p,:])
        C = un[:p,:]

    # Computation of system matrices B & D
    # Output Error Minimization

    # Setting up the Phi matrix
    Phi  = np.zeros((p*ns, r+p*q+r*q))
    CA_powers = np.zeros((ns, p, A.shape[1]))
    CA_powers[0, :, :] = C
    A_p = A
    for pwr in range(1,ns):
        CA_powers[pwr,:,:] =  C@A_p
        A_p = A@A_p

    # First block column of Phi
    for df in range(ns):
        Phi[df*p:(df+1)*p,:r] = CA_powers[df,:,:]

    # Second block column of Phi
    Ipp = np.eye(p)
    for i in range(ns):
        Phi[i*p:(i+1)*p, r:r+p*q] = np.kron(input[i,:],Ipp)

    # Third block column of Phi
    In1n1 = np.eye(r)
    cc = r + p*q + 1
    dd = r + p*q + r*q

    krn = np.array([np.kron(input[i,:],In1n1) for i in range(ns)])

    with multiprocessing.Pool(pool_size) as pool:
        for i,res in progress_bar(
                pool.imap_unordered(
                    partial(_blk_3,CA=CA_powers,U=np.flip(krn,0)),
                    range(1,ns),
                    200
                ),
                total = ns
            ):
            Phi[i*p:(i+1)*p,cc-1:dd] = res

    y = output[:ns,:].flatten()

    teta = lsqminnorm(Phi,y)

    x0 = teta[:r]
    dcol = teta[r:r+p*q]
    bcol = teta[r+p*q:r+p*q+r*q]

    D = np.zeros((p,q))
    B = np.zeros((r,q))
    for wq in range(q):
        D[:,wq] = dcol[wq*p:(wq+1)*p]

    for ww in range(q):
        B[:,ww] = bcol[ww*r:(ww+1)*r]

    assert A.shape[0] == A.shape[1] == r

    return A,B,C,D
