import numpy as np
from numpy import zeros, pi, size
import scipy.linalg
from scipy.linalg import fractional_matrix_power as matpow

linsolve = scipy.linalg.solve

lsqminnorm = lambda *args: np.linalg.lstsq(*args, rcond=None)[0]

def parse_okid(args, config):
    """
    p determines order of the observer Kalman ARX filter used in OKID-ERA-DC.
    n determines size of the state-space model used for representing the system.
    """
    return config

# y = output data: forced response. dimensions of y: p x nt, where nt = number of timesteps
# u = input data: input. dimensions of u: q x nt.
# m = number of Markov parameters (impulse response timesteps) to solve for
def get_impulse(y, u):
    p = y.shape[0]
    q = u.shape[0]
    nt = y.shape[1]
    assert y.shape[1] == u.shape[1]
    # form a blockwise upper triangular matrix of inputs
    U = np.zeros((q*nt, nt))
    for i in range(nt):
        for j in range(nt-i):
            U[q*j:q*(j+1),i] = u[:,i]
    Y_impulse = y @ np.linalg.pinv(U)
    assert Y_impulse.shape == (p, q*nt)
    return Y_impulse

def okid2(y, u, m):
    p = y.shape[0]
    q = u.shape[0]
    nt = y.shape[1]
    assert y.shape[1] == u.shape[1]
    # form a blockwise upper triangular matrix of inputs
    U = np.zeros((q*nt, nt))
    for i in range(nt):
        for j in range(nt-i):
            U[q*j:q*(j+1),i] = u[:,i]
    Y_impulse = y @ np.linalg.pinv(U)
    assert Y_impulse.shape == (p, q*nt)
    return Y_impulse

def OKID_brunton(y,u,m):
    
    p = y.shape[0]
    q = u.shape[0]
    nt = y.shape[1]
    assert y.shape[1] == u.shape[1]
    
    # Form data matrices y and V
    V = np.zeros((q+(q+p)*m,nt))
    for i in range(nt):
        V[:q,i] = u[:q,i]
        
    for i in range(1,m+1):
        for j in range(nt-i):
            vtemp = np.concatenate((u[:,j],y[:,j]))
            V[q+(i-1)*(q+p):q+i*(q+p),i+j] = vtemp
    
    # Solve for observer Markov parameters Ybar
    Ybar = y @ np.linalg.pinv(V,rcond=10**(-3))
    
    # Isolate system Markov parameters H, and observer gain M
    D = Ybar[:,:q] # feed-through term (or D matrix) is the first term
    
    Y = np.zeros((p,q,m))
    Ybar1 = np.zeros((p,q,m))
    Ybar2 = np.zeros((p,q,m))
    
    for i in range(m):
        Ybar1[:,:,i] = Ybar[:,q+(q+p)*i : q+(q+p)*i+q]
        Ybar2[:,:,i] = Ybar[:,q+(q+p)*i+q : q+(q+p)*(i+1)]
    
    Y[:,:,0] = Ybar1[:,:,0] + Ybar2[:,:,0] @ D
    for k in range(1,m):
        Y[:,:,k] = Ybar1[:,:,k] + Ybar2[:,:,k] @ D
        for i in range(k-1):
            Y[:,:,k] += Ybar2[:,:,i] @ Y[:,:,k-i-1]
            
    H = np.zeros((D.shape[0],D.shape[1],m+1))
    H[:,:,0] = D
    
    for k in range(1,m+1):
        H[:,:,k] = Y[:,:,k-1]
        
    return H

# Y = output data: response to unit impulse, or "impulse response," or "Markov parameters".
# dimensions of Y: p x q x nt, where nt = number of timesteps = number of Markov parameters = number of blocks
# mc = number of block rows in Hankel matrix = order of controllability matrix
# mo = number of block columns in Hankel matrix = order of observability matrix
# p = number of outputs
# q = number of inputs
# r = reduced model order = dimension of reduced A = newly assumed dimension of state variable
def era(Y,mo,mc,p,q,r):
    # get D from first p x q block of impulse response
    Dr = Y[:,:,0]  # first block of output data
    
    assert Y.shape[:2] == (p,q)  # sanity check that we're passing in the right output data
    assert Y.shape[2] >= mo+mc   # make sure there are enough timesteps to assemble this size of Hankel matrix
    
    # make impulse response into Hankel matrix and shifted Hankel matrix
    H = np.zeros((p*(mo), q*(mc+1)))
    for i in range(mo):
        for j in range(mc+1):
            H[p*i:p*(i+1), q*j:q*(j+1)] = Y[:,:,i+j+1]
    H0 = H[:,:-q]
    H1 = H[:,q:]
    assert H0.shape == H1.shape == (p*(mo), q*(mc))

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

    return (Ar,Br,Cr,Dr,S)

# l = initial lag for data correlations
# g = lags between correlation matrices
# a = number of block rows in Hankel of correlation matrix
# b = number of block columns in Hankel of correlation matrix
def era_dc(Y,mo,mc,p,q,r,l,g,a,b):
    # get D from first p x q block of impulse response
    Dr = Y[:,:,0]  # first block of output data
    
    assert Y.shape[:2] == (p,q)  # sanity check that we're passing in the right output data
    assert Y.shape[2] >= l+(a+1+b+1)*g+mo+mc   # make sure there are enough timesteps to assemble the Hankel matrices
    
    # make impulse response into Hankel matrix
    H = np.zeros((p*(mo), q*(mc+l+(a+1+b+1)*g))) # Hankel of Markov parameters
    for i in range(mo):
        for j in range(mc+l+(a+1+b+1)*g):
            H[p*i:p*(i+1), q*j:q*(j+1)] = Y[:,:,i+j+1]
    H0 = H[:,:q*mc]
    assert H0.shape == (p*(mo), q*(mc))

    dimR = p*mo # Dimension of square correlation matrices
    dimHRl = (dimR*(a+1), dimR*(b+1)) # Dimension of Hankels of correlation matrices
    HRl = np.zeros(dimHRl) # Hankel of correlation matrices
    HRl1 = np.zeros(dimHRl) # Shifted Hankel of correlation matrices
    for i in range(a+1):
        for j in range(b+1):
            Hl = H[:, q*(l+1+(i+j)*g):q*(l+1+(i+j)*g+mc)]
            Hl1 = H[:, q*(l+1+(i+j)*g+1):q*(l+1+(i+j)*g+mc+1)]
            assert Hl.shape == Hl1.shape == (p*(mo), q*(mc))
            R = Hl@H0
            R1 = Hl1@H0
            assert R.shape == R1.shape == (dimR,dimR)
            HRl[dimR*i:dimR*(i+1), dimR*j:dimR*(j+1)] = R
            HRl1[dimR*i:dimR*(i+1), dimR*j:dimR*(j+1)] = R1

    # reduced SVD of Hankel matrix
    def _svd(*args):
        U,S,V = scipy.linalg.svd(*args, lapack_driver="gesvd")
        return U,S,V.T.conj()
    U,S,V = _svd(HRl)
    SigmaInvSqrt = np.diag(S[:r]**-0.5)
    SigmaSqrt = np.diag(S[:r]**0.5)
    Ur = U[:,:r]
    Vr = V[:,:r]

    # get A from SVD and shifted Hankel matrix
    Ar = SigmaInvSqrt @ Ur.T.conj() @ HRl1 @ Vr @ SigmaInvSqrt

    # get B and C
    Br = ((SigmaInvSqrt @ Ur.T.conj())[:,:dimR] @ H0)[:,:q]
    Cr = (Ur @ SigmaSqrt)[:p,:]

    return (Ar,Br,Cr,Dr,S)


def _era_ya(kmax, m, l, r, Y, n, svd="gesvd"):
    # ERA/DC
    #
    # Form Hankel matrix (H) from the Markov parameters (Y) (Eq. 3.80)

    # Obtain Hankel Matrix of Zeroth Order & First Order
    H0,H1 = np.zeros((2, kmax*m, l*r))
    #H0 = hankel(Y,0)
    for hj in range(kmax):
        for jh in range(l):
            H0[hj*m:hj*m+m, jh*r:jh*r+r] = Y[jh+hj+1]
            H1[hj*m:hj*m+m, jh*r:jh*r+r] = Y[jh+hj+2]
    
    if svd in ["gesvd", "gesdd"]:
        import scipy.linalg
        def _svd(*args):
            U,S,V = scipy.linalg.svd(*args, lapack_driver=svd)
            return U,S,V.T.conj()

    ## 2c. Use H matrix to compute system matrices A, B & C
    R1,sis,S1 = _svd(H0) # singular value decomposition
    Sis = np.diag(sis[:n])

    Siv = np.diag(1/np.sqrt(sis[:n]))
    # A: state transition matrix (Eqs. 3.32 & 3.82)
    A  = Siv@R1[:,:n].T@H1@S1[:,:n]@Siv
    # B: input influence matrix (Eqs. 3.32 & 3.83)
    B = (np.sqrt(Sis)@S1[:,:n].T)[:,:r]
    Observability = R1[:,:n]@np.sqrt(Sis)
    # C: output influence matrix (Eqs. 3.34 & 3.84)
    C  = Observability[:m,:]

    return A,B,C


def okid(dati, dato, svd="gesvd", debug=False, **config):
    """
    PART 2: OKID-ERA-DC (Observer Kalman filter Identification -
            Eigen Realization with Direct Correlations)
    -----------------------------------------------------------------

    # Inputs
    Modelparameters
        div = 1;     # A parameter used for decimating data. 1 uses entire data without downsampling.
        mro = 10;    # Model reduction order
        orm = 4;     # Order of the model. # of computed and plotted modes dep on orm.
                     # For orm = 2, one mode is found, for orm = 4, two modes are found.
    svd:
        GESDD: a two-stage algorithm. It first reduces the input matrix to bidiagonal form via 
        Householder reflection, then divides the bidiagonal matrix into smaller matrices to calculate 
        singular values and the unitary matrices U and V. Finally, the singular values of the entire 
        input matrix is simply the those of the smaller sub-matrices.
        
        GESVD: also a two-stage algorithm where the first step is identical to GESDD. The second step 
        can be done using iterative QR decomposition to derive singular values. More mathematical details are available here.

    Sometimes higher orm still gives fewer modes, e.g. orm = 8 for case 1 gives
    three modes, but one of them is invalid according to the EMAC & MPC criteria.
    kmax = 100;  #Number of computed Markov parameters, indicated as 1000 on page 43 of
    (Arici & Mosalam, 2006). However, it was input as 100 in the code.
    kmax = 100 runs much faster & both kmax = 100 & 1000 give the same results.

    # Outputs
    1. freqdamp variable is a matrix that includes the information of identified
       frequencies, damping ratios & validation of the modes with MPC & EMAC criteria
       Each row of freqdamp corresponds to a mode. Columns are as follows:
       1) frequency, 2) damping ratio, 3) order index, 4) condition number, 5) input EMAC,
       6) output EMAC, 7) MPC. If values in columns 5-7 are > 0.5, identified mode is valid.
    2. modeshape stores the mode shape information for identified modes.
    3. RMSEpred: root mean square error of the predicted output from
       identified parameters with respect to the actual output.
    4. Markovparamerror: root mean square error used to validate accurate
       computation of Markov parameters

    # Description of the Methodology
    More information on OKID-ERA-DC method can be found in
    Section 3.4.6 of (Arici & Mosalam, 2006). Equations below refer to this report.
    OKID-ERA-DC is a Multiple Input - Multiple Output (MIMO) System Identification (SI) method
    that consists of the following steps:
    1. Data pre-processing (baseline correction, filtering & decimation)
    2. Identify observer Kalman filters using Observer Kalman filter Identification (OKID)
       methodology. This step basically consists of developing a simple observer model of
       a system from a MIMO ARX structure, Eq. (3.76), which is broken into these 6 steps:
   
    3a. Determine Markov parameters (M) in a least squares sense, Eq. (3.76).

    3b. Establish the Hankel matrix (H) from the Markov parameters, (Eq. 3.80).
    3c. Use H to compute system matrices A, B & C, in which modal information is embedded.
    // 3d. Obtain the modal information from matrices A, B & C.
    3e. Spatial & temporal validation of the identified modes.
    3f. Back calculate (estimate) the output accelerations with the state-space system &
        check against the actual output accelerations.
    """ 
    dt = to = config.get("dt", None) or dati["time_step"]
    p = config.get("p", config.get("mro")) # # steps used for the identification (ie, prediction horizon)
    n = n1 = config.get("n", config.get("orm", 4))  # Order of the model.
    
    if svd in ["gesvd", "gesdd"]:
        import scipy.linalg
        def _svd(*args):
            U,S,V = scipy.linalg.svd(*args, lapack_driver=svd)
            return U,S,V.T.conj()

    elif svd in ["xla", "jax"]:
        from jax.scipy.linalg import svd as _svd

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
        
    dati = dati[:-1,:]
    dato = dato[:-1,:]

    #verbose = True
    #dn = config.get("dn", None) or dati.shape[0]
    kmax    = config["kmax"]
    n = config["orm"]   # assign order of model input to variable n, consistent with Eqs. 3.81-3.84
    p = config["mro"]   # assign input model reduction order to variable p, consistent with Eq. 3.76


    l,m = dato.shape # m is the number of output channels
    _,r = dati.shape # r is the number of input channels

    # ASSERT SIZE(DATO,1) == SIZE(DATI,1)

    ## 2a. Obtain Observer Markov Parameters
    # The Markov parameters are computed in two steps:
    # i)  The Observer Markov matrices are computed from Eq. 3.76 using a linear regression approach
    # ii) Compute the system Markov parameters from the Markov matrices using recursive relations

    # Compute matrix U that represents ARX equation of current output on p time steps of past output
    # & input values (Eq. 3.76)
    U = np.zeros(((m+r)*p+r, l))
    U[:r,:] = dati.T
    for b in range(1,p+1):
        U[(b-1)*(r+m)+r:(b-1)*(r+m)+r+r+m, b:] = [*dati[:-b,:r].T, *dato[:-b, :m].T]


    # i) Compute the matrix of Observer Markov Parameter Matrix (M)
    #    in Eq 3.76 using Linear Regression
    uu,wr,V = _svd(U,0)
    pg = (r+m)*p+r
    for i in range((r+m)*p+r):
        if wr[i] <= 0.001:
            pg = i
            break

    s = np.diag(wr)

    pss = V[:,:pg]@linsolve(s[:pg,:pg], uu[:,:pg].T)
    M = dato.T@pss           # M: Observer Markov Parameter Matrix

    # Fit for multiple regression
    # RMSE between actual output & y on left hand side in Eq. 3.76. It should be 
    # quite small (e.g., 10^-3) for accurately computed Markow parameters.
    # ypreo = M@U
    # markovParamError = 1/m*sum((
    #     sum((dato[:,i] - ypreo[i,:].T)**2)/sum(dato[:,i]**2)
    #     for i in range(m)
    # ))  


    ## ii) Compute Markov parameters (Y) using recursive relations in Eqs. 3.78 & 3.79
    #      (Equation 6.21 in Juang 1994)

    # Matrix D is directly equal to the Observer Markov parameter matrix (Eq. 3.77)
    D = M[:, :r]  # D: Direct Transmission term

    Y = [D]
    # First p steps (Eq. 3.78)
    for i in range(p):
        Y.append(
            M[:, r+i*(r+m):r+i*(r+m)+r] + sum(
                M[:, r+j*(r+m)+r:r+j*(r+m)+r+m]@Y[i-j]
                for j in range(i+1)
            )
        )
    # Remainder (Eq. 3.79)
    for i in range(p, l+kmax):
        sumt = zeros((m,r))
        for j in range(p):
            sumt += M[:,r+j*(r+m)+r:r+j*(r+m)+r+m]@Y[i-j]

        Y.append(sumt)

    # The Markow parameters Y have been computed.

    A,B,C = _era_ya(kmax, m, l, r, Y, _svd, n)

    if debug:
        return locals()

    return A,B,C,D


def validate(freqdmp, v, system, **config):
    """2e. Validation Analysis

    Two criteria are used for selection of identified genuine modes
    (in the sense of spatial & temporal consistency).

    a) Modal Phase Collinearity (MPC) testing spatial consistency of identification results.
       Modes having MPC value above 0.5 (mpc parameter below) are considered as genuine modal quantities.

    b) Exted Modal Amplitude Coherence (EMAC), evaluates temporal consistency of the identification results.
       Both output EMAC & input EMAC can be computed. Input EMAC requires the controllability matrix.
       Because the controllability matrix is not estimated by all considered SI methods,
       this criterion is computed, but not used.
       Modes with output EMAC values < 0.5 are considered spurious & therefore not reported.

    """
    mpc = MPC(n, v, system)

    # Add the input EMAC, output EMAC, and MPC to the matrix freqdamp
    for lih in range(size(freqdmp)[0]):
        freqdmp[lih,4] = emacif(freqdmp[lih,2])
        freqdmp[lih,5] = emacof(freqdmp[lih,2])
        freqdmp[lih,6] = mpc[freqdmp(lih,3)]
        if freqdmp[lih,5]>0.5 and freqdmp[lih,7]>0.5:
            # valid
            return True
        else:
            # not valid
            return False

def MPC(n, v, system):
    """a) Modal Phase Collinearity (MPC) [Eqs. 3.85-3.87]"""
    _,__,C,___ = system
    modes_raw = C@v
    _, n = modes_raw.shape
    sxx, syy, sxy = np.zeros((3, *modes_raw.shape))
    nu, mpc = np.zeros((2, n))
    lam = np.zeros((2, n))
    for i in range(n):
        sxx[:,i] = np.real(modes_raw[:,i]).T@np.real(modes_raw[:,i])
        syy[:,i] = np.imag(modes_raw[:,i]).T@np.imag(modes_raw[:,i])
        sxy[:,i] = np.real(modes_raw[:,i]).T@np.imag(modes_raw[:,i])
        nu[i]    = (syy[:,i]-sxx[:,i])/(2*sxy[:,i])
        lam[1,i] = (sxx[:,i]+syy[:,i])/2 + sxy[:,i]*np.sqrt(nu[i]**2+1);
        lam[2,i] = (sxx[:,i]+syy[:,i])/2 - sxy[:,i]*np.sqrt(nu[i]**2+1);
        mpc[i]   = ((lam[0,i]-lam[1,i])/(lam[0,i]+lam[1,i]))**2;
    return mpc

def EMAC_Matrix(n, m, pto, ptop, debug=False):
    emac = np.zeros((n, m))
    for i in range(n):
        for j in range(m):
            Rij = min(
                (abs(pto[j,i])/abs(ptop[j,i])),
                (abs(ptop[j,i])/abs(pto[j,i]))
            )
            Pij = np.angle(pto[j,i]/ptop[j,i])
            Wij = max(
                1 - abs(Pij)/(pi/4),
                0
            )
            emac[i,j] = Rij*Wij
    return emac

def EMAC_Variation(n,m,pto,ptop):
    pass


def EnergyCondensedEMAC(n,m,emac,phi,debug=False):
    "Equation 3.93"
    return np.array([
        sum(emac[i,j]*abs(phi[i,j])**2 
            for j in range(m)
        )/np.real(phi[i,:].conjugate().transpose()@phi[i,:])
        for i in range(n) 
    ])

def MAC(shape,v, v_inv, A, B):
    n,m,r,l = shape
    lamb = v_inv@A@v
    bkh = v_inv@B
    for i in range(n):
        for j in range(l):
            qhat[i,j*r+1:j*r+r] = bkh[i,:]*(lamb[i,i])**j

    selsiz = min(size(qlin),size(qhat))

    for hnd in range(n):
        ql = qlin[hnd,:selsiz(2)]
        qh = qhat[hnd,:selsiz(2)]
        mac[hnd] = abs(ql*qh.T)/(abs(ql*ql.T)*abs(qh*qh.T))**0.5


def EMAC_Validate(shape, kmax, v, d, dt, system):
    "b) Exted Modal Amplitude Coherence (EMAC)"
    A,B,C,D = system
    n,m,r,l = shape

def OutputEMAC(n,m,kmax,A,C,Observability, debug=False,d=None,v=None, **_):
    """Output EMAC (Eqs. 3.88-3.89)"""
    if d is None:
        assert v is None
        from .ExtractModes import condeig 
        v, d, _ = condeig(A)
    pto = Observability[-m:,:]@v # the identified value at T0 ( last block row)
    ptop = C@v@np.diag(d**(kmax-1))
    emaco  = EMAC_Matrix(n,m,pto,ptop)
    if debug:
        return locals()
    return EnergyCondensedEMAC(n, m, emaco, (C@v).T)


def InputEMAC(Ctrl,A,B):
    #
    # Input EMAC
    #
    # # EMAC Input Variation
    # for i in range(l):
    #     qtovar = qlin[:,i*r:(i+1)*r]
    #     qtopvar = (np.diag(d)**i)*inm
    #     emaci = EMAC_Matrix(n,m,qtovar,qtopvar)
    #     emacivar[:,i] = EnergyCondensedEMAC(n,m,emaci,inm);

    # Pick the last block column
    v, d  = A.eigen()
    v_inv = np.linalg.inv(v)
    inm   = linsolve(v,B)   # Initial mode contribution
    qlin  = v_inv@Ctrl;     # Modal controllability Matrix (F' in Pappa 1993)
    qto   = qlin[:, (l-1)*r+1:l*r]
    qtop  = d**(l-1)@inm
    emaci = EMAC_Matrix(n, m, qto, qtop)
    return EnergyCondensedEMAC(n, m, emaci, inm)
    


if __name__ == "__main__":
    from pathlib import Path
    import quakeio
    # import ssid as si
    import okid
    import numpy as np
    channels = dict( # PAINTER RIO DELL
        inputs  = [17, 3, 20],
        outputs = [ 9, 7 , 4]
    )
    channels = dict( # HAYWARD
        inputs  = [17, 18, 24, 25],
        outputs = [19, 20 , 22, 23]
    )
    for file in Path("hwd_recs").glob("RioDell_P*.zip"):
        event = quakeio.read(file)
        try:
            inputs = [
                event.match("r", file_name=f".*{chan}.*").accel
                for chan in channels["inputs"]
            ]
            outpts = [
                event.match("r", file_name=f".*{chan}.*").accel
                for chan in channels["outputs"]
            ]
        except:
            print(f"failed {file.name}")    
    for file in Path("painter").glob("RioDell_P*.zip"):
        event = quakeio.read(file)
        try:
            inputs = [
                event.match("r", file_name=f".*{chan}.*").accel
                for chan in channels["inputs"]
            ]
            outpts = [
                event.match("r", file_name=f".*{chan}.*").accel
                for chan in channels["outputs"]
            ]
        except:
            print(f"failed {file.name}")

        else:
            dt = inputs[0]["time_step"]
            V = OKID.okid(inputs, outpts, dt=dt, kmax=500, mro=10, orm=4, verbose=True)
            break

            A,B,C,D = OKID.okid(inputs, outpts, dt=dt, kmax=500, mro=10, orm=4, verbose=True)
            freqdmpSRIM, modeshapeSRIM, *_ = si.ComposeModes(dt, A, B, C, D)
            # Add validation
            print(si.IdentifiedSystem(dt, A, B, C, D))
            # print(file, np.real(1/freqdmpSRIM[:,0]))
