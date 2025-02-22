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
                    :math:`q` = number of inputs, and :math:`nt` = number of timesteps
    :type inputs:   array
    :param outputs: output response history.
                    dimensions: :math:`(p,nt)`, where :math:`p` = number of outputs, and
                    :math:`nt` = number of timesteps
    :type outputs:  array
    :param horizon: (optional) number of steps used for identification (prediction horizon).
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
        warnings.warn("input data has more channels (dim 1) than timesteps (dim 2)")
    if outputs.shape[0] > outputs.shape[1]:
        warnings.warn("output data has more channels (dim 1) than timesteps (dim 2)")

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
    p,q,nt = Y.shape # p = number of outputs, q = number of inputs, nt = number of timesteps

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
    # make sure there are enough timesteps to assemble this size of Hankel matrix
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
    D = Y[:,:,0]  # first block of output data

    # size of Hankel matrix
    if no is None:
        no = min(150, int((nt-1)/2))
    if nc is None:
        nc = min(150, max(nt-1-no, int((nt-1)/2)))
    # make sure there are enough timesteps to assemble the Hankel matrices
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

# Can delete anything after this line
 

from mdof.utilities import n4sid_utils
'''
The code implements the N4SID (Subspace Identification) algorithm to identify a state-space model of a linear time-invariant system based on the input data (inputs) and output data (outputs). The main process is as follows:

1. Construct the Hankel Matrix  
   The input and output data are organized into a Hankel matrix suitable for subspace identification using the functions construct_Hankel and stacked_Hankel.
2. QR Decomposition and Partitioning  
   The Hankel matrix is decomposed using QR decomposition. The resulting matrix is then partitioned into several sub-blocks according to specified row and column sizes, which will be used for subsequent projection matrix and subspace calculations.
3. Projection Matrix Calculation
   By multiplying the partitioned blocks (using np.dot), a series of intermediate results (such as Li1, Li2, etc.) are obtained. These results are further combined using np.hstack and np.vstack to construct a new matrix used for extracting the system's observable subspace.
4. Singular Value Decomposition (SVD) and Order Selection
   An SVD is performed on a key matrix, result_matrix_2. The function choose_k_gap is used to determine the system order k based on the gap between singular values. Then, the matrices U and Sigma are partitioned into their principal and secondary parts, which are used to compute the subspace matrix Gamma_i, among others.
5. Least-Squares Solution for the System Matrices
   Further combinations of the truncated U and Sigma matrices are used to construct and solve a set of linear equations, ultimately obtaining the state-space matrices A, B, C, and D.
6. Estimation of Noise Covariance Matrices (Qs, Ss, Rs) 
   The residual matrix, ρ₂ (rho_2), is computed, and based on this residual, the noise covariance matrices Qs, Ss, and Rs are estimated.
7. Simulation and Visualization
   The identified system (A, B, C, D) is used with the input data in the simulate function to obtain the predicted output, output_prediction. Then, Matplotlib is used to plot both the actual output and the predicted output for comparison. To correct a sign inversion issue, the predicted output is multiplied by -1 (i.e., -output_prediction.T) before plotting.

Through these steps, the function ultimately returns the four system matrices (A, B, C, D).
'''
def n4sid(inputs, outputs, **options):
    verbose=options.get('verbose',False)
    i = options.get("i", None)
    j = options.get("j", None) 
    m = inputs.shape[0]  
    l = outputs.shape[0]
    
    U = n4sid_utils.construct_Hankel(inputs, j, 0, 2 * i - 1)
    Y = n4sid_utils.construct_Hankel(outputs, j, 0, 2 * i - 1)
    
    #STEP_1 computer hankel matrix
    hankel_matrix = n4sid_utils.stacked_Hankel(inputs, outputs, j, 0, 2*i-1)
    #STEP2_Performing QR decomposition and partitioning the R matrix
    Q, R = np.linalg.qr(hankel_matrix.T)
    RT, QT = R.T, Q.T
    RT_row_sizes = [3*i, 3, 3*(i-1), i, 1, i-1] 
    RT_col_sizes = [3*i, 3, 3*(i-1), i, 1, i-1] 
    QT_row_sizes = [3*i, 3, 3*(i-1), i, 1, i-1]  
    QT_col_sizes = [j]  
    RT_blocks = n4sid_utils.slice_matrix(RT, RT_row_sizes, RT_col_sizes)
    QT_blocks = n4sid_utils.slice_matrix(QT, QT_row_sizes, QT_col_sizes)
    R_row_sizes = [3*i, 3, 3*(i-1), i, 1, i-1] 
    R_col_sizes = [3*i, 3, 3*(i-1), i, 1, i-1] 
    R_blocks = n4sid_utils.slice_matrix(R, R_row_sizes, R_col_sizes)
    #STEP_3 Calculating the projection
    R5614 = n4sid_utils.extract_and_combine_blocks(RT_blocks, row_range=(4, 6), col_range=(0, 4))
    RT1414 = n4sid_utils.extract_and_combine_blocks(R_blocks, row_range=(0, 4), col_range=(0, 4))
    result_matrix = np.dot(R5614, RT1414)

    l_i = l * i 
    m_i = m * i 
    Li1 = result_matrix[:l_i, :m_i] 
    Li2 = result_matrix[:l_i, m_i:2*m_i]  
    Li3 = result_matrix[:l_i, 2*m_i:]
    R6615 = n4sid_utils.extract_and_combine_blocks(RT_blocks, row_range=(5, 6), col_range=(0, 5))
    RT1515 = n4sid_utils.extract_and_combine_blocks(R_blocks, row_range=(0, 5), col_range=(0, 5))
    result_matrix_1 = np.dot(R6615, RT1515)
    l_i = l * i 
    m_i = m * i 
    Li_11 = result_matrix_1[:l*(i-1), :m*(i+1)] 
    Li_12 = result_matrix_1[:l*(i-1), m*(i+1):2*m_i]  
    Li_13 = result_matrix_1[:l*(i-1), 2*m_i:]
    #STEP4_Determine gamma and n through the SVD
    Zero_Matrix_1 = np.zeros((l_i, m_i))
    Combined_Matrix = np.hstack((Li1, Zero_Matrix_1, Li3))
    R1414 = n4sid_utils.extract_and_combine_blocks(RT_blocks, row_range=(0, 4), col_range=(0, 4))
    R5514 = n4sid_utils.extract_and_combine_blocks(RT_blocks, row_range=(4, 5), col_range=(0, 4))
    result_matrix_2 = np.dot(Combined_Matrix, R1414)
    U, Sigma, VT = np.linalg.svd(result_matrix_2, full_matrices=False)
    k = n4sid_utils.choose_k_gap(Sigma, threshold_ratio=10)
    print(k)
    U1 = U[:, :k]
    U2 = U[:, k:]
    Sigma1 = np.diag(Sigma[:k])
    Sigma2 = np.diag(Sigma[k:])
    Sigma1_half = np.sqrt(Sigma1) 
    Gamma_i = np.dot(U1, Sigma1_half)
    Gamma_i_pinv = np.linalg.pinv(Gamma_i)
    Gamma_i_minus_1 = Gamma_i[:-l, :]
    Gamma_i_minus_1_pinv = np.linalg.pinv(Gamma_i_minus_1)
    #Step_5 Construct all the matrices required for the least squares method.
    U1_truncated = U1[:-l, :]
    U1_truncated_pinv = np.linalg.pinv(U1_truncated)
    U1_transposed = U1.T
    Sigma1_neg_half = np.diag(1 / np.sqrt(np.diag(Sigma1)))
    R1514 = n4sid_utils.extract_and_combine_blocks(RT_blocks, row_range=(0, 5), col_range=(0, 4))
    R5514 = n4sid_utils.extract_and_combine_blocks(RT_blocks, row_range=(4, 5), col_range=(0, 4))
    R2214 = n4sid_utils.extract_and_combine_blocks(RT_blocks, row_range=(1, 2), col_range=(0, 4))
    Zero_Matrix_2 = np.zeros((l*(i-1), m*(i-1)))
    Combined_Matrix_2 = np.hstack((Li_11, Zero_Matrix_2, Li_13))
    intermediate_result_1 = np.dot(np.dot(np.dot(Sigma1_neg_half, U1_truncated_pinv), Combined_Matrix_2), R1514)
    final_matrix_3 = np.vstack((intermediate_result_1, R5514))
    intermediate_result_2 = np.dot(np.dot(np.dot(Sigma1_neg_half, U1_transposed), Combined_Matrix), R1414)
    final_matrix_4 = np.vstack((intermediate_result_2, R2214))
    #Step_6 Solve for the ABCD matrix.
    AAA = final_matrix_4
    BBB = final_matrix_3
    AAA_pinv = np.linalg.pinv(AAA)
    script_L = np.dot(BBB, AAA_pinv) 
    script_L11 = script_L[:k, :k]    
    script_L12 = script_L[:k, k:k+m]   
    script_L21 = script_L[k:k+l, :k]     
    script_L22 = script_L[k:k+l, k:k+m] 
    A = script_L11
    B = script_L12
    C = script_L21
    D = script_L22
    #Step_7 Solve for the QsSsRs matrix.
    rho_2 = BBB - np.dot(script_L, AAA)
    rho1_2 = rho_2[:k, :]
    rho2_2 = rho_2[k:, :]
    rho1_2_transposed = rho1_2.T
    rho2_2_transposed = rho2_2.T

    Qs = (1/j) * np.dot(rho1_2, rho1_2_transposed)
    Ss = (1/j) * np.dot(rho1_2, rho2_2_transposed)
    Rs = (1/j) * np.dot(rho2_2, rho2_2_transposed)

    system = (A, B, C, D)
    
    from mdof.simulate import simulate
    output_prediction = simulate(system, inputs)
    import matplotlib.pyplot as plt
    plt.plot(outputs.T)
    print(outputs)
    plt.show()
    plt.plot(-output_prediction.T)
    plt.show()

    if verbose:
        print("Hankel matrix:", hankel_matrix)

    return A, B, C, D
