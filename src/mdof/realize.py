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


def _ac2bd(u, y, a, c):
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
    from scipy.io import savemat
    savemat("out/phiz_ours.mat", {'phi_ours':phi,'z_ours':z})
    b = np.linalg.pinv(phi) @ z.T
    x0 = b[:n, 0]
    d = _block_tr(m, r, 1, b[n:n+m*r], 0)
    b = _block_tr(n, r, 1, b[n+m*r:n+m*r+n*r], 0)

    return b, d, x0


def ac2bd(inputs, outputs, r, A, C, **options):

    p,N = outputs.shape
    q,_ = inputs.shape
    
    threads = options.get("threads",6)
    chunk = options.get("chunk", 200)
    
    lsq_solve = numerics.lsq_solver(options.get("lsq", {}))
    
    # Setting up the Phi matrix
    Phi  = np.zeros((p*N, r+p*q+r*q), dtype=complex)

    # First block column of Phi
    CA_powers = _form_powers(C,A,N,p,r,options.get("pow_impl",1),options.get("pmax",10))
    # CA_powers = np.zeros((N, p, r))
    # CA_powers[0] = C
    # A_p = A
    # for pwr in range(1,N):
    #     CA_powers[pwr] =  C@A_p
    #     A_p = A@A_p
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


def _blk_3(i, CA, U):
    return i, np.einsum('kil,klj->ij', CA[:i,:,:], U[-i:,:,:])


def _form_powers(C,A,ns,p,r,pow_impl,p_max):
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


def srim(inputs,outputs,**options):
    """
    System realization from input and output data, with output error minimization method.
    System Realization Using Information Matrix (SRIM) [4]_.
    
    :param inputs:  input time history. dimensions: :math:`(q,nt)`, where
                    :math:`q` = number of inputs, and :math:`nt` = number of timesteps
    :type inputs:   array
    :param outputs: output response history.
                    dimensions: :math:`(p,nt)`, where :math:`p` = number of outputs, and
                    :math:`nt` = number of timesteps
    :type outputs:  array
    :param horizon: number of steps used for identification (prediction horizon).
                    default: :math:`min(300, nt)`
    :type horizon:  int, optional
    :param order:   model order. default: :math:`min(20,` ``horizon``:math:`/2)`
    :type order:    int, optional
    :param full:    if True, full SVD. default: True
    :type full:     bool, optional
    :param find:    "ABCD" or "AC". default: "ABCD"
    :type find:     string, optional
    :param threads: number of threads used during the output error minimization method.
                    default: 6
    :type threads:  int, optional
    :param chunk:   chunk size in output error minimization method. default: 200
    :type chunk:    int, optional

    :return:        realization in the form of state space coefficients ``(A,B,C,D)``
    :rtype:         tuple of arrays

    References
    ----------
    .. [4]  Juang, J. N. (1997). System realization using information matrix. Journal
            of Guidance, Control, and Dynamics, 20(3), 492-500.
            (https://doi.org/10.2514/2.4068)
    """
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

    no = options.get("no",
         options.get("horizon",
                     min(300, nt)))

    r = options.get("r", 
        options.get("order",
                    min(20, int(no/2))))

    full = options.get("full", True)
    find = options.get("find","ABCD")

    # maximum possible number of columns in the Y and U data matrices
    ns = options.get("ns",nt-1-no+2)

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

    assert A.shape[0] == A.shape[1] == r

    if "b" not in find.lower() and "d" not in find.lower():
        return (A,None,C,None)
    
    # Computation of system matrices B & D
    # Output Error Minimization
    use_juang_ac2bd = options.get("use_juang_ac2bd", False)
    if use_juang_ac2bd:
        B, D, x0 = _ac2bd(inputs.T, outputs.T, A, C)
    else:
        B, D = ac2bd(inputs, outputs, r, A, C, **options)
    return A,B,C,D

    
def era(Y,**options):
    """
    System realization from Markov parameters (discrete impulse response data).
    Ho-Kalman / Eigensystem Realization Algorithm (ERA) [1]_ [2]_.

    :param Y:       Markov parameters. dimensions: :math:`(p,q,nt)`, where :math:`p` = number of outputs,
                    :math:`q` = number of inputs, and :math:`nt` = number of Markov parameters.
    :type Y:        array
    :param horizon: number of block rows in Hankel matrix = order of observability matrix.
                    default: :math:`min(150, (nt-1)/2)`
    :type horizon:  int, optional
    :param nc:      number of block columns in Hankel matrix = order of controllability matrix.
                    default: :math:`min(150, max(nt-1-` ``horizon``:math:`, (nt-1)/2))`
    :type nc:       int, optional
    :param order:   model order. default: :math:`min(20,` ``horizon``:math:`/2)`
    :type order:    int, optional

    :return:        realization in the form of state space coefficients ``(A,B,C,D)``
    :rtype:         tuple of arrays

    References
    ----------
    .. [1]  Ho, Β. L., & Kálmán, R. E. (1966). Effective construction of linear state-variable models
            from input/output functions: Die Konstruktion von linearen Modeilen in der Darstellung
            durch Zustandsvariable aus den Beziehungen für Ein-und Ausgangsgrößen. at-Automatisierungstechnik,
            14(1-12), 545-548. (https://doi.org/10.1524/auto.1966.14.112.545)
    .. [2]  Juang, J. N., & Pappa, R. S. (1985). An eigensystem realization algorithm for modal parameter
            identification and model reduction. Journal of guidance, control, and dynamics, 8(5), 620-627.
            (https://doi.org/10.2514/3.20031)
    """
    p,q,nt = Y.shape # p = number of outputs, q = number of inputs, nt = number of timesteps

    no = options.get("no",
         options.get("horizon",
                     None))
    nc = options.get("nc",
                     None)

    # get D from first p x q block of impulse response
    Dr = Y[:,:,0]  # first block of output data

    # size of Hankel matrix
    if no is None:
        no = min(150, int((nt-1)/2))
    if nc is None:
        nc = min(150, max(nt-1-no, int((nt-1)/2)))
    # make sure there are enough timesteps to assemble this size of Hankel matrix
    assert nt >= no+nc

    r = options.get("r", 
        options.get("order",
                    min(20, int(no/2))))

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


def era_dc(Y,**options):
    """
    System realization from Markov parameters (discrete impulse response data).
    Eigensystem Realization Algorithm with Data Correlations (ERA/DC) [3]_.


    :param Y:       Markov parameters. dimensions: :math:`(p,q,nt)`, where :math:`p` = number of outputs,
                    :math:`q` = number of inputs, and :math:`nt` = number of Markov parameters.
    :type Y:        array
    :param horizon: number of block rows in Hankel matrix = order of observability matrix.
                    default: :math:`min(150, (nt-1)/2)`
    :type horizon:  int, optional
    :param nc:      number of block columns in Hankel matrix = order of controllability matrix.
                    default: :math:`min(150, max(nt-1-` ``horizon``:math:`, (nt-1)/2))`
    :type nc:       int, optional
    :param order:   model order. default: :math:`min(20,` ``horizon``:math:`/2)`
    :type order:    int, optional
    :param a:       :math:`(\\alpha)` number of block rows in Hankel of correlation matrix. default: 0
    :type a:        int, optional
    :param b:       :math:`(\\beta)` number of block columns in Hankel of correlation matrix. default: 0
    :type b:        int, optional
    :param l:       initial lag for data correlations. default: 0
    :type l:        int, optional
    :param g:       lags (gap) between correlation matrices. default: 1
    :type g:        int, optional

    :return:        realization in the form of state space coefficients ``(A,B,C,D)``
    :rtype:         tuple of arrays

    References
    ----------
    .. [3]  Juang, J. N., Cooper, J. E., & Wright, J. R. (1987). An eigensystem realization algorithm
            using data correlations (ERA/DC) for modal parameter identification.
            (https://ntrs.nasa.gov/citations/19870035963)    
    """
    p,q,nt = Y.shape # p = number of outputs, q = number of inputs, nt = number of timesteps

    no = options.get("no",
         options.get("horizon",
                     None))
    nc = options.get("nc",
                     None)
    a = options.get("a", 0)
    b = options.get("b", 0)
    l = options.get("l", 0)
    g = options.get("g", 1)

    # get D from first p x q block of impulse response
    Dr = Y[:,:,0]  # first block of output data

    # size of Hankel matrix
    if no is None:
        no = min(150, int((nt-1)/2))
    if nc is None:
        nc = min(150, max(nt-1-no, int((nt-1)/2)))
    # make sure there are enough timesteps to assemble the Hankel matrices
    assert nt >= l+(a+1+b+1)*g+no+nc

    r = options.get("r", 
        options.get("order",
                    min(20, int(no/2))))

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
            R = Hl@H0.T       # correlation matrix
            R1 = Hl1@H0.T     # shifted correlation matrix
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
    