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


import numpy as np
from scipy.linalg import svd
from numpy.linalg import pinv

class StateSpaceModel:
    def __init__(self, A, B, C, D, K=None):
        self.A = A
        self.B = B
        self.C = C
        self.D = D
        

    def __repr__(self):
        return f"StateSpaceModel(A={self.A}, B={self.B}, C={self.C}, D={self.D})"

def n4sid(data, nx, Ts=1, method='auto', **kwargs):
    """
    Estimate state-space model using subspace methods.
    Arguments:
    - data: (numpy array) The input-output data matrix.
    - nx: (int or str) Model order or 'best' for automatic model order determination.
    - Ts: (float) Sample time of the model, set to 0 for continuous models.
    - method: (str) 'auto' or other methods defined by user.
    - kwargs: Additional keyword arguments for future extensions or specific settings.
    Returns:
    - StateSpaceModel: The estimated state-space model.
    """
    # Data preprocessing
    U, s, Vh = svd(data, full_matrices=False)
    n = min(len(s), nx) if isinstance(nx, int) else np.argmax(s < 1e-10) + 1

    # Construct the state-space model matrices
    A = np.random.randn(n, n)  # Placeholder for actual computation
    B = np.random.randn(n, 1)  # Placeholder for actual computation
    C = np.random.randn(1, n)  # Placeholder for actual computation
    D = np.zeros((1, 1))       # Assuming no direct feedthrough



    model = StateSpaceModel(A, B, C, D)
    return model




from scipy.linalg import qr


def construct_hankel_matrix(u, y, i, j):
   
    # Ensure the sequences are long enough to avoid index errors.
    if len(u) < 2*i-1 or len(y) < 2*i-1:
        raise ValueError("The length of input or output sequences is insufficient for matrix construction.")
    
    # Construct U_0|2i-1.
    U_matrix = np.array([u[k:k+j] for k in range(2*i-1)])  # Rows are 2i-1, columns are j.
    
    # Construct Y_0|2i-1.
    Y_matrix = np.array([y[k:k+j] for k in range(2*i-1)])  # Rows are 2i-1, columns are j.

    # Normalize the matrix by dividing by sqrt(j).
    H_matrix = np.vstack((U_matrix, Y_matrix)) / np.sqrt(j)

    return H_matrix

# Example input, replace these with actual data.
u = np.random.rand(10)  # Example input data
y = np.random.rand(10)  # Example output data
i = 3  # Example value for i
j = 2  # Example value for j

H_matrix = construct_hankel_matrix(u, y, i, j)

# Perform QR decomposition.
Q, R = np.linalg.qr(H_matrix)

# Define blocks of the R matrix.
R_blocks = {
    "R11": R[0:i, 0:i],
    "R21": R[i:i+1, 0:i],
    "R22": R[i:i+1, i:i+1],
    "R31": R[i+1:2*i, 0:i],
    "R32": R[i+1:2*i, i:i+1],
    "R33": R[i+1:2*i, i+1:2*i],
    "R41": R[2*i:3*i, 0:i],
    "R42": R[2*i:3*i, i:i+1],
    "R43": R[2*i:3*i, i+1:2*i],
    "R44": R[2*i:3*i, 2*i:3*i],
    "R51": R[3*i:3*i+1, 0:i],
    "R52": R[3*i:3*i+1, i:i+1],
    "R53": R[3*i:3*i+1, i+1:2*i],
    "R54": R[3*i:3*i+1, 2*i:3*i],
    "R55": R[3*i:3*i+1, 3*i:3*i+1],
    "R61": R[3*i+1:4*i, 0:i],
    "R62": R[3*i+1:4*i, i:i+1],
    "R63": R[3*i+1:4*i, i+1:2*i],
    "R64": R[3*i+1:4*i, 2*i:3*i],
    "R65": R[3*i+1:4*i, 3*i:3*i+1],
    "R66": R[3*i+1:4*i, 3*i+1:4*i]
}

# Define blocks of the Q matrix.
Q_blocks = {
    "Q1": Q[0:i, 0:j],        
    "Q2": Q[i:i+1, 0:j],      
    "Q3": Q[i+1:2*i, 0:j],    
    "Q4": Q[2*i:3*i, 0:j],    
    "Q5": Q[3*i:3*i+1, 0:j],  
    "Q6": Q[3*i+1:4*i-1, 0:j] 
}

# Print the dimensions of each block to ensure they're correct.
for key, block in R_blocks.items():
    print(f"{key} Shape: {block.shape}")

for key, block in Q_blocks.items():
    print(f"{key} Shape: {block.shape}")

# Print QR decomposition results.
print("\nQ Matrix:")
print(Q)

print("\nR Matrix:")
print(R)


from mdof.utilities import n4sid_utils
def n4sid(inputs, outputs, **options):
    inputs, outputs = n4sid_utils.preprocess_data(u, y)
    i = options.get("i", 3)
    j = options.get("j", 26) 
    m = inputs.shape[0]  
    l = outputs.shape[0]

    U_hankel = n4sid_utils.construct_hankel(inputs, j, 0, 2*i-1)
    Y_hankel = n4sid_utils.construct_hankel(outputs, j, 0, 2*i-1)
    stacked_hankel = n4sid_utils.stacked_hankel(inputs, outputs, j, 0, 5)

    Q, R = np.linalg.qr(stacked_hankel.T)
    RT, QT = R.T, Q.T
    RT_blocks, QT_blocks = n4sid_utils.partition_R_matrix(RT, QT, i, j, m, l)
    new_matrix_R5614, new_matrix_RT1414 = n4sid_utils.compute_projection_matrices(RT_blocks, R_blocks, i, l, m)
    Li1, Li2, Li3, Li_11, Li_12, Li_13 = n4sid_utils.compute_Li_matrices(new_matrix_R5614, new_matrix_RT1414, i, l, m)
    Gamma_i, U1, Sigma1 = n4sid_utils.compute_gamma_and_svd(Li1, Li3, i, l, m)
    A, B, C, D = n4sid_utils.compute_state_space_matrices(U1, Sigma1, RT_blocks, R_blocks, i, l, m)

    n4sid_utils.get_dominant_coords()
    
    return A,B,C,D


