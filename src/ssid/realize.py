
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

## ERA ##
# Y = output data: response to unit impulse, or "impulse response," or "Markov parameters".
# dimensions of Y: p x q x nt, where nt = number of timesteps = number of Markov parameters = number of blocks
# no = number of block rows in Hankel matrix = order of observability matrix
# nc = number of block columns in Hankel matrix = order of controllability matrix
# r = reduced model order = dimension of reduced A = newly assumed dimension of state variable
def era(Y,no=None,nc=None,r=None,**options):
    p,q,nt = Y.shape # p = number of outputs, q = number of inputs, nt = number of timesteps
    if r is None:
        r = min(20, int(nt/2))

    # get D from first p x q block of impulse response
    Dr = Y[:,:,0]  # first block of output data

    # size of Hankel matrix
    if no is None:
        if nc is None:
            no = nc = min(300, int((nt-1)/2))
        else:
            no = min(300, int(nt-1-nc))
    elif nc is None:
        nc = min(300, int(nt-1-no))
    else:
        # make sure there are enough timesteps to assemble this size of Hankel matrix
        assert nt >= no+nc

    # make impulse response into Hankel matrix and shifted Hankel matrix
    H = np.zeros((p*(no), q*(nc+1)))
    for i in range(no):
        for j in range(nc+1):
            H[p*i:p*(i+1), q*j:q*(j+1)] = Y[:,:,i+j+1]
    H0 = H[:,:-q]
    H1 = H[:,q:]
    assert H0.shape == H1.shape == (p*(no), q*(nc))

    # reduced SVD of Hankel matrix
    _svd = numerics.svd_routine(**options.get("svd", {}))

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

## ERA-DC ##
# a = (alpha) number of block rows in Hankel of correlation matrix
# b = (beta) number of block columns in Hankel of correlation matrix
# l = initial lag for data correlations
# g = lags (gap) between correlation matrices
def era_dc(Y,no=None,nc=None,a=0,b=0,l=0,g=1,r=None,**options):
    p,q,nt = Y.shape # p = number of outputs, q = number of inputs, nt = number of timesteps
    if r is None:
        r = int(nt/2)

    # get D from first p x q block of impulse response
    Dr = Y[:,:,0]  # first block of output data

    # size of Hankel matrix
    if no is None:
        if nc is None:
            no = nc = min(300, int(nt/2-1))
        else:
            no = min(300, int(nt/2-1))
    elif nc is None:
        nc = min(300, int(nt-no-2))
    # make sure there are enough timesteps to assemble the Hankel matrices
    assert nt >= l+(a+1+b+1)*g+no+nc

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
    _svd = numerics.svd_routine(**options.get("svd", {}))

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


## SRIM ##
# inputs = input data. dimensions of inputs: q x nt, where nt = number of timesteps.
# outputs = response data due to input data. dimensions of outputs: p x nt.
# no = number of steps used for identification (prediction horizon),
    # and the order of the observability matrix.
    # analagous to order (number of autoregressors) of the observer Kalman
    # ARX filter used in OKID-ERA-DC.
# r = size of the state-space model used for representing the system.
# Juang 1997, "System Realization Using Information Matrix," Journal of Guidance, Control, and Dynamics
def srim(inputs,outputs,no=None,r=None,full=True,pool_size=6,**options):
    lsq_solve = numerics.lsq_solver(options.get("lsq", {}))

    if len(inputs.shape) == 1:
        inputs = inputs[None,:]
    if len(outputs.shape) == 1:
        outputs = outputs[None,:]

    if inputs.shape[0] > inputs.shape[1]:
        warnings.warn("input data has more channels (dim 1) than timesteps (dim 2)")
    if outputs.shape[0] > outputs.shape[1]:
        warnings.warn("output data has more channels (dim 1) than timesteps (dim 2)")

    q,nt = inputs.shape
    p = outputs.shape[0]
    assert nt == outputs.shape[1]

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

    # progress_bar = lambda arg, **kwds: (i for i in arg)

    # Construct Y (output) & U (input) data matrices (Eqs. 3.58 & 3.60 Arici 2006)
    for i in range(no):
        Yno[i*p:(i+1)*p,:] = outputs[:,i:ns+i]
        Uno[i*q:(i+1)*q,:] = inputs[:,i:ns+i]


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
    Rhh = Ryy - Ruy.T@lin_solve(Ruu,Ruy)

    # 2c. Obtain observability matrix using full or partial decomposition 
    #     (Eqs. 3.72 & 3.74).
    if full:
        # Full Decomposition Method
        un,*_ = np.linalg.svd(Rhh,0)           # Eq. 3.74
        Observability = un[:,:r]               # Eq. 3.72
        A = lsq_solve(Observability[:(no-1)*p,:], Observability[p:no*p,:])
        C = Observability[:p,:]
    else:
        # Partial Decomposition Method
        un,*_ = np.linalg.svd(Rhh[:,:(no-1)*p+1],0)
        Observability = un[:,:r]
        A = lsq_solve(Observability[:(no-1)*p,:], Observability[p:no*p,:])
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
        Phi[i*p:(i+1)*p, r:r+p*q] = np.kron(inputs[:,i],Ipp)

    # Third block column of Phi
    In1n1 = np.eye(r)
    cc = r + p*q + 1
    dd = r + p*q + r*q

    krn = np.array([np.kron(inputs[:,i],In1n1) for i in range(ns)])

    # Execute a loop in parallel that looks something like:
    #    for i in  range(1,ns):
    #        Phi[] = _blk_3(i, CA_Powers, np.flip(...))

    with multiprocessing.Pool(pool_size) as pool:
        for i,res in progress_bar(
                pool.imap_unordered(
                    partial(_blk_3, CA=CA_powers,U=np.flip(krn,0)),
                    range(1,ns),
                    200
                ),
                total = ns
            ):
            Phi[i*p:(i+1)*p,cc-1:dd] = res

    y = outputs[:,:ns].flatten()

    teta = lsq_solve(Phi,y)

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

def subspace(outputs,no=None,r=None,cov_driven=True,**options):
    outputs = np.atleast_2d(outputs)
    p,nt = outputs.shape
    assert nt > p

    if no == None:
        no = min(300, int(nt/2))
    if r == None:
        r = min(10,no-1)

    if cov_driven:
        Toep = np.zeros((p*no,p*no))
        for i in range(2*no-1):
            covs = np.zeros((p,p,no-1))
            for k in range(no-1):
                covs[:,:,k] = outputs[:,k+i]@outputs[k]
            Ri = np.sum(covs,axis=3)/no
            for j in range(i):
                for l in range(j):
                    Toep[j,-l-1] = Ri

    else:
        pass

    return
