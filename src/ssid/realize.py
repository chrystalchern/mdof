import numpy as np
import scipy

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