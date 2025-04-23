import numpy as np
lin_solve = np.linalg.solve
import multiprocessing
from functools import partial
import warnings
try:
    from tqdm import tqdm as progress_bar
except:
    def progress_bar(arg, **kwds): return arg
from .import numerics

from .evolve import obsv2ac
from .influence import ac2bd


def srim(inputs,outputs,**options):
    r"""
    System realization from input and output data, with output error minimization method.
    System Realization Using Information Matrix (SRIM) [1]_.
    
    :param inputs:  input time history. dimensions: :math:`(q,nt)`, where
                    :math:`q` = number of inputs, and :math:`nt` = number of time samples
    :type inputs:   array
    :param outputs: output response history.
                    dimensions: :math:`(p,nt)`, where :math:`p` = number of outputs, and
                    :math:`nt` = number of time samples
    :type outputs:  array
    :param horizon: (optional) number of samples used for identification (prediction horizon).
                    default: :math:`\min(300, nt)`
    :type horizon:  int
    :param order:   (optional) model order. default: :math:`\min(20,` ``horizon``:math:`/2)`
    :type order:    int
    :param full:    (optional) if True, full SVD. default: True
    :type full:     bool
    :param find:    (optional) "ABCD" or "AC". default: "ABCD"
    :type find:     string
    :param threads: (optional) number of threads used during the output error minimization method.
                    default: 6
    :type threads:  int
    :param chunk:   (optional) chunk size in output error minimization method. default: 200
    :type chunk:    int

    :return:        realization in the form of state space coefficients ``(A,B,C,D)``
    :rtype:         tuple of arrays

    References
    ----------
    .. [1]  Juang, J. N. (1997). System realization using information matrix. Journal
            of Guidance, Control, and Dynamics, 20(3), 492-500.
            (https://doi.org/10.2514/2.4068)
    """

    if len(inputs.shape) == 1:
        inputs = inputs[None,:]
    if len(outputs.shape) == 1:
        outputs = outputs[None,:]

    if inputs.shape[0] > inputs.shape[1]:
        warnings.warn("input data has more channels (dim 1) than time samples (dim 2)")
    if outputs.shape[0] > outputs.shape[1]:
        warnings.warn("output data has more channels (dim 1) than time samples (dim 2)")

    q,nt = inputs.shape
    p = outputs.shape[0]
    assert nt == outputs.shape[1]

    no = options.get("no",
         options.get("horizon",
                     min(300, nt)))

    n = options.get("n", 
        options.get("r",
        options.get("order",
                    min(20, int(no/2)))))

    full = options.get("full", True)
    find = options.get("find","ABCD")

    # maximum possible number of columns in the Y and U data matrices
    ns = options.get("ns",nt-1-no+2)

    assert no >= n/p + 1    # make sure prediction horizon is large enough that
                            # observability matrix is full rank (Juang Eq. 8)

    Yno = np.zeros((p*no,ns))
    Uno = np.zeros((q*no,ns))

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
    if full:
        un,*_ = np.linalg.svd(Rhh,full_matrices=False)                # Eq. 3.74 Arici 2006; Eq. 25 Juang
    else:
        un,*_ = np.linalg.svd(Rhh[:,:(no-1)*p+1],full_matrices=False) # Eq. 3.74 Arici 2006; Eq. 30 Juang
    
    Observability = un[:,:n]                        # Eq. 3.72 Arici 2006; Eq. 27 Juang

    # Computation of system matrices A & C
    A,C = obsv2ac(Observability, no, p, **options)

    assert A.shape[0] == A.shape[1] == n
    assert C.shape[0] == p
    assert C.shape[1] == n

    if "b" not in find.lower() and "d" not in find.lower():
        return (A,None,C,None)
    
    # Computation of system matrices B & D
    # Output Error Minimization
    B, D = ac2bd(inputs, outputs, A, C, **options)

    return (A,B,C,D)


def era(Y,**options):
    r"""
    System realization from Markov parameters (discrete impulse response data).
    Ho-Kalman / Eigensystem Realization Algorithm (ERA) [2]_ [3]_.

    :param Y:       Markov parameters. dimensions: :math:`(p,q,nt)`, where :math:`p` = number of outputs,
                    :math:`q` = number of inputs, and :math:`nt` = number of Markov parameters.
    :type Y:        array
    :param horizon: (optional) number of block rows in Hankel matrix = order of observability matrix.
                    default: :math:`\min(150, (nt-1)/2)`
    :type horizon:  int
    :param nc:      (optional) number of block columns in Hankel matrix = order of controllability matrix.
                    default: :math:`\min(150, max(nt-1-` ``horizon``:math:`, (nt-1)/2))`
    :type nc:       int
    :param order:   (optional) model order. default: :math:`\min(20,` ``horizon``:math:`/2)`
    :type order:    int

    :return:        realization in the form of state space coefficients ``(A,B,C,D)``
    :rtype:         tuple of arrays

    References
    ----------
    .. [2]  Ho, Β. L., & Kálmán, R. E. (1966). Effective construction of linear state-variable models
            from input/output functions: Die Konstruktion von linearen Modeilen in der Darstellung
            durch Zustandsvariable aus den Beziehungen für Ein-und Ausgangsgrößen. at-Automatisierungstechnik,
            14(1-12), 545-548. (https://doi.org/10.1524/auto.1966.14.112.545)
    .. [3]  Juang, J. N., & Pappa, R. S. (1985). An eigensystem realization algorithm for modal parameter
            identification and model reduction. Journal of guidance, control, and dynamics, 8(5), 620-627.
            (https://doi.org/10.2514/3.20031)
    """
    p,q,nt = Y.shape # p = number of outputs, q = number of inputs, nt = number of time samples

    no = options.get("no",
         options.get("horizon",
                     None))
    nc = options.get("nc",
                     None)

    # get D from first p x q block of impulse response
    D = Y[:,:,0]  # first block of output data

    # size of Hankel matrix
    if no is None:
        no = min(150, int((nt-1)/2))
    if nc is None:
        nc = min(150, max(nt-1-no, int((nt-1)/2)))
    # make sure there are enough time samples to assemble this size of Hankel matrix
    assert nt >= no+nc

    n = options.get("n",
        options.get("r", 
        options.get("order",
                    min(20, int(no/2)))))

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
    SigmaInvSqrt = np.diag(S[:n]**-0.5)
    SigmaSqrt = np.diag(S[:n]**0.5)
    Ur = U[:,:n]
    Vr = V[:,:n]

    # get A from SVD and shifted Hankel matrix
    A = SigmaInvSqrt @ Ur.T.conj() @ H1 @ Vr @ SigmaInvSqrt

    # get B and C
    B = (SigmaSqrt @ Vr.T.conj())[:,:q]
    C = (Ur @ SigmaSqrt)[:p,:]

    return (A,B,C,D)


def era_dc(Y,**options):
    r"""
    System realization from Markov parameters (discrete impulse response data).
    Eigensystem Realization Algorithm with Data Correlations (ERA/DC) [4]_.


    :param Y:       Markov parameters. dimensions: :math:`(p,q,nt)`, where :math:`p` = number of outputs,
                    :math:`q` = number of inputs, and :math:`nt` = number of Markov parameters.
    :type Y:        array
    :param horizon: (optional) number of block rows in Hankel matrix = order of observability matrix.
                    default: :math:`\min(150, (nt-1)/2)`
    :type horizon:  int
    :param nc:      (optional) number of block columns in Hankel matrix = order of controllability matrix.
                    default: :math:`\min(150, max(nt-1-` ``horizon``:math:`, (nt-1)/2))`
    :type nc:       int
    :param order:   (optional) model order. default: :math:`\min(20,` ``horizon``:math:`/2)`
    :type order:    int
    :param a:       (optional) :math:`(\\alpha)` number of block rows in Hankel of correlation matrix. default: 0
    :type a:        int
    :param b:       (optional) :math:`(\\beta)` number of block columns in Hankel of correlation matrix. default: 0
    :type b:        int
    :param l:       (optional) initial lag for data correlations. default: 0
    :type l:        int
    :param g:       (optional) lags (gap) between correlation matrices. default: 1
    :type g:        int

    :return:        realization in the form of state space coefficients ``(A,B,C,D)``
    :rtype:         tuple of arrays

    References
    ----------
    .. [4]  Juang, J. N., Cooper, J. E., & Wright, J. R. (1987). An eigensystem realization algorithm
            using data correlations (ERA/DC) for modal parameter identification.
            (https://ntrs.nasa.gov/citations/19870035963)    
    """
    p,q,nt = Y.shape # p = number of outputs, q = number of inputs, nt = number of time samples

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
    D = Y[:,:,0]  # first block of output data

    # size of Hankel matrix
    if no is None:
        no = min(150, int((nt-1)/2))
    if nc is None:
        nc = min(150, max(nt-1-no, int((nt-1)/2)))
    # make sure there are enough time samples to assemble the Hankel matrices
    assert nt >= l+(a+1+b+1)*g+no+nc

    n = options.get("n",
        options.get("r", 
        options.get("order",
                    min(20, int(no/2)))))

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
    SigmaInvSqrt = np.diag(S[:n]**-0.5)
    SigmaSqrt = np.diag(S[:n]**0.5)
    Ur = U[:,:n]
    Vr = V[:,:n]

    # get A from SVD and shifted Hankel matrix of correlation matrices
    A = SigmaInvSqrt @ Ur.T.conj() @ HRl1 @ Vr @ SigmaInvSqrt

    # get B and C
    B = ((SigmaInvSqrt @ Ur.T.conj())[:,:dimR] @ H0)[:,:q]
    C = (Ur @ SigmaSqrt)[:p,:]

    return (A,B,C,D)


from mdof.utilities import n4sid_utils
def n4sid(inputs, outputs, **options):
    """
    Numerical Algorithms for Subspace State Space System Identification
    """
   
    i = options.get("i", None)
    j = options.get("j", None) 
    m = inputs.shape[0]  
    l = outputs.shape[0]

    if options.get('simple', False):
        return # TODO: replace this with the Deterministic Algorithm 1


    stacked_hankel = n4sid_utils.stacked_hankel(inputs, outputs, j, 0, 2*i-1)

    _, R = np.linalg.qr(stacked_hankel.T)
    R_inv = np.linalg.inv(R.T)  
    RT = R.T                  
    R_inv_blocks, RT_blocks = n4sid_utils.partition_R_matrices(R_inv, RT, i, j, m, l)
    new_matrix_R5614, new_matrix_RT1414, new_matrix_R6615, new_matrix_RT1515 = n4sid_utils.compute_projection_matrices(RT_blocks, R_inv_blocks)
    Li1, _, Li3 = n4sid_utils.compute_Li_matrices(new_matrix_R5614, new_matrix_RT1414, i, l, m)
    Li_11, _, Li_13 = n4sid_utils.compute_Li_1_matrices(new_matrix_R6615, new_matrix_RT1515, i, l, m)
    _, U1, Sigma1, k = n4sid_utils.compute_gamma_and_svd(Li1, Li3, i, l, m, RT_blocks, threshold=1e-3)
    verbose = options.get("verbose", True)
    if verbose:
        print("Selected system order (k):", k)
    A, B, C, D, Qs, Ss, Rs = n4sid_utils.compute_state_space_matrices(U1, Sigma1, RT_blocks, i, l, m, Li1, Li3, Li_11, Li_13, k, j)

    
    
    return A,B,C,D, Qs, Ss, Rs



def deterministic(inputs, outputs, i, j):
    m = inputs.shape[0]  
    l = outputs.shape[0]
    #Step 1 
    def construct_Hankel(data, j, start, finish):
        dim_data, nps = data.shape
        Hankel = np.empty(((finish - start + 1) * dim_data, j))
        for k in range(start, finish + 1):
            Hankel[k * dim_data:(k + 1) * dim_data, :] = data[:, start + k:start + j + k]
        return Hankel
    # Build the Hankel matrices for both input and output data.
    U_big = construct_Hankel(inputs, j, 0, 2 * i - 1) 
    Y_big = construct_Hankel(outputs, j, 0, 2 * i - 1) 

    
    # Partition the Hankel matrices into "past" (p) and "future" (f) blocks.
    U_p = U_big[0 : i*m, :]  
    U_f = U_big[i*m : 2*i*m, :] 
    Y_p = Y_big[0 : i*l, :]
    Y_f = Y_big[i*l : 2*i*l, :]   
    # Combine the past input and output data into one matrix.
    W_p = np.vstack([U_p, Y_p]) 
    # Construct additional shifted blocks needed for the oblique projection:
    U_p_plus = U_big[0 : (i+1)*m, :] # Past inputs with an extra block row at the bottom of the same Hankel matrix
    U_f_minus = U_big[(i-1)*m : 2*(i-1)*m, :] # Future inputs with the first block row removed of the same Hankel matrix
    U_p_minus = U_big[0 : (i-1)*m, :] # Past inputs excluding the last block row of the same Hankel matrix
    U_f_plus = U_big[(i+1)*m : 2*(i+1)*m, :] # Future inputs with an extra block row at the bottom of the same Hankel matrix

    Y_p_plus = Y_big[0 : (i+1)*l, :]
    Y_f_minus = Y_big[(i-1)*l : 2*(i-1)*l, :]
    Y_f_plus = Y_big[(i+1)*l : 2*(i+1)*l, :]
    Y_p_minus = Y_big[0 : (i-1)*l, :]

    W_p_plus = np.vstack([U_p_plus, Y_p_plus])
    W_p_minus = np.vstack([U_p_minus, Y_p_minus])


    #Step 2 Calculate the oblique projections
    def oblique_projection_rows(A, B, C):
        # Compute intermediate products.
        CCt = C @ C.T  # (r, r)
        CBt = C @ B.T  # (r, q)
        BCt = B @ C.T  # (q, r)
        BBt = B @ B.T  # (q, q)
        M_top = np.hstack([CCt, CBt])  # (r, r+q)
        M_bot = np.hstack([BCt, BBt])  # (q, r+q)
        M = np.vstack([M_top, M_bot])  # ((r+q), (r+q))
        # Compute the pseudoinverse of M.
        M_dag = np.linalg.pinv(M)     
        ACt = A @ C.T   # shape (p, r)
        ABt = A @ B.T   # shape (p, q)
        top_block = np.hstack([ACt, ABt])   # (p, r+q)
        prod_block = top_block @ M_dag
        prod_block_Cpart = prod_block[:, :C.shape[0]]  # (p, r)
        P = prod_block_Cpart @ C
        return P
    # Compute the oblique projection for the future output block.
    O_i = oblique_projection_rows(Y_f, U_f, W_p)
    # Compute a shifted oblique projection for the "minus" block.
    O_i_minus = oblique_projection_rows(Y_f_minus, U_f_minus, W_p_plus)



    # Calculate the SVD of the weighted oblique projection
    W1 = np.eye(O_i.shape[0])  # W1 Chosen as the identity matrix.
    W2 = W_p.T @ np.linalg.pinv(W_p @ W_p.T) @ W_p

    #W2 = np.eye(O_i.shape[1])  # W2 Chosen as the identity matrix.
    U, S, Vt = np.linalg.svd(W1 @ O_i @ W2, full_matrices=False)

    # 3. Determine the order by inspecting the singular values in S and partition the SVD accordingly to obtain Up and Si
    rank_n = np.sum(S > 1e-6*max(S))  
    S1 = np.diag(S[:rank_n]) 
    U1 = U[:, :rank_n] 
    V1 = Vt[:rank_n, :]

    # Compute the extended observability matrix.
    Gamma_i = np.linalg.inv(W1) @ U1 @ np.sqrt(S1)
    Gamma_i_pinv = np.linalg.pinv(Gamma_i)
    Gamma_i_minus = Gamma_i_minus = Gamma_i[:-l, :]
    Gamma_i_minus_pinv = np.linalg.pinv(Gamma_i_minus)
    X_f_d = Gamma_i_pinv @ O_i
    X_f_d_next = Gamma_i_minus_pinv @ O_i_minus

    # Solve the set of linear equations for A, B, C and D
    i_th_block_U = U_big[(i-1)*m : i*m, :]
    i_th_block_Y = Y_big[(i-1)*l : i*l, :]
    M_left = np.vstack((X_f_d_next, i_th_block_Y))
    M_right = np.vstack((X_f_d, i_th_block_U))
    sys_matrix = M_left @ np.linalg.pinv(M_right)
    A, B = sys_matrix[:rank_n, :rank_n], sys_matrix[:rank_n, rank_n:]
    C, D = sys_matrix[rank_n:, :rank_n], sys_matrix[rank_n:, rank_n:]
    
    return A, B, C, D

