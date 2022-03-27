import numpy as np
from numpy import zeros, pi, size
import scipy.linalg

linsolve = scipy.linalg.solve

lsqminnorm = lambda *args: np.linalg.lstsq(*args, rcond=None)[0]

def okid(dati, dato, svd="gesvd", **config):
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
        _svd = lambda *args: scipy.linalg.svd(*args, lapack_driver=svd)
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

    verbose = True
    dn = config.get("dn", None) or dati.shape[0]
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
    uu,wr,v = _svd(U,0)
    pg = (r+m)*p+r
    for lop in range((r+m)*p+r):
        if wr[lop] <= 0.001:
            pg = lop
            break

    s = np.diag(wr)

    pss = v.T[:,:pg]@linsolve(s[:pg,:pg], uu[:,:pg].T)
    M = dato.T@pss           # M: Observer Markov Parameter Matrix

    # Fit for multiple regression
    # RMSE between actual output & y on left hand side in Eq. 3.76. It should be 
    # quite small (e.g., 10^-3) for accurately computed Markow parameters.
    ypreo = M@U
    markovParamError = 1/m*sum((
        sum((dato[:,i] - ypreo[i,:].T)**2)/sum(dato[:,i]**2)
        for i in range(m)
    ))  


    ## ii) Compute Markov parameters (Y) using recursive relations in Eqs. 3.78 & 3.79
    #      (Equation 6.21 in Juang 1994)

    # Matrix D is directly equal to the Observer Markov parameter matrix (Eq. 3.77)
    M1 = D = M[:, :r]  # D: Direct Transmission term

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


    #
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
    
    ## 2c. Use H matrix to compute system matrices A, B & C
    R1,sis,S1 = _svd(H0) # singular value decomposition
    Sis = np.diag(sis[:n])

    Siv = np.diag(1/np.sqrt(sis[:n]))
    # A: state transition matrix (Eqs. 3.32 & 3.82)
    A  = Siv@R1[:,:n].T@H1@S1[:,:n]@Siv
    # B: input influence matrix (Eqs. 3.32 & 3.83)
    B = (np.sqrt(Sis)@S1[:,:n].T)[:,:r]
    Pb = R1[:,:n]@np.sqrt(Sis)
    # C: output influence matrix (Eqs. 3.34 & 3.84)
    C  = Pb[:m,:]
    return locals()
    return A,B,C,D


def validate(inputs, outputs, system, **config):
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
    A,B,C,D = system
    def mpc(modes_raw):
        """a) Modal Phase Collinearity (MPC) [Eqs. 3.85-3.87]"""
        for q in range(n):
            sxx[:,q] = np.real(modes_raw[:,q]).T@np.real(modes_raw[:,q])
            syy[:,q] = np.imag(modes_raw[:,q]).T@np.imag(modes_raw[:,q])
            sxy[:,q] = np.real(modes_raw[:,q]).T@np.imag(modes_raw[:,q])
            nu[q]    = (syy[:,q]-sxx[:,q])/(2*sxy[:,q])
            lam[1,q] = (sxx[:,q]+syy[:,q])/2 + sxy[:,q]*np.sqrt(nu[q]**2+1);
            lam[2,q] = (sxx[:,q]+syy[:,q])/2 - sxy[:,q]*np.sqrt(nu[q]**2+1);
            mpc[q]   = ((lam[0,q]-lam[1,q])/(lam[0,q]+lam[1,q]))**2;


    def emac():
        "b) Exted Modal Amplitude Coherence (EMAC)"
        v_inv = np.linalg.inv(v)
        qlin = v_inv@Qb;     # Controllability Matrix used for the input-EMAC
        plin = Pb@v;         # Observability Matrix used for the output-EMAC

        lamb = v_inv@A@v
        bkh = v_inv@B;
        for hn in range(n):
            for ll in range(l):
                qhat[hn,ll*r+1:ll*r+r] = bkh[hn,:]*(lamb[hn,hn])**ll

        selsiz = min(size(qlin),size(qhat))

        for hnd in range(n):
            ql = qlin[hnd,:selsiz(2)]
            qh = qhat[hnd,:selsiz(2)]
            mac[hnd] = abs(ql*qh.T)/(abs(ql*ql.T)*abs(qh*qh.T))**0.5

        #
        # Output EMAC (Eqs. 3.88-3.89)
        #
        # Pick the last block row
        pto = plin[(kmax-1)*m+1:m*kmax,:]  # the identified value at T0
        for ds in range(n):
            ptop[:,ds] = modes_raw[:,ds]*exp(sj1[ds]*dt*(kmax-1));

        # Computation of rij
        for qa in range(n):
            for qz in range(m):
                Rij[qa,qz] = min((abs(pto(qz,qa))/abs(ptop(qz,qa))),(abs(ptop(qz,qa))/abs(pto(qz,qa))));
                Pij = angle(pto(qz,qa)/ptop(qz,qa))
                Pijn[qa,qz] = Pij
                if abs(Pij) <= pi/4:
                    Wij[qa,qz] = 1 - abs(Pij)/(pi/4)
                else:
                    Wij[qa,qz] = 0

                emaco[qa,qz] = Rij[qa,qz]*Wij[qa,qz];      # emaco is the ouput emac

        #
        # Input EMAC
        #
        # Pick the last block column
        qto = qlin[:, (l-1)*r+1:l*r]
        qtop = d**(l-1)*inm;

        # EMAC Input Variation
        for er in range(l):
            qtovar = qlin[:,er*r:(er+1)*r];
            qtopvar = d**er*inm;
            # Computation of rik
            for qak in range(n):
                for qzk in range(r):
                    Rik[qak,qzk] = min(
                            (abs(qtovar[qak,qzk])/abs(qtopvar[qak,qzk])),
                            (abs(qtopvar[qak,qzk])/abs(qtovar[qak,qzk]))
                    )
                    Pik = angle(qtovar[qak,qzk]/qtopvar[qak,qzk])
                    if abs(Pik)<=pi/4:
                        Wik[qak,qzk] = 1 - abs(Pik)/(pi/4)
                    else:
                        Wik[qak,qzk] = 0.0
                    emaci[qak,qzk] = Rik[qak,qzk]*Wik[qak,qzk]

            # Weight for emaci
            emacif = zeros((n,1))
            for xc in range(n):
                sumi = 0;
                for lw in range(r):
                    sumi = emaci[xc,lw]*(inm[xc,lw]*inm[xc,lw].T) + sumi

                emacif[xc] = sumi/(inm[xc,:]*inm[xc,:].T)
            emacivar[:,er] = emacif;

        # Computation of rik
        for qak in range(n):
            for qzk in range(r):
                Rik[qak,qzk] = min((abs(qto[qak,qzk])/abs(qtop[qak,qzk])),(abs(qtop[qak,qzk])/abs(qto(qak,qzk))))
                Pik = angle(qto[qak,qzk]/qtop[qak,qzk]);
                if abs(Pik) <= pi/4:
                    Wik[qak,qzk] = 1-abs(Pik)/(pi/4)
                else:
                    Wik[qak,qzk] = 0.0

                emaci[qak,qzk] = Rik[qak,qzk]*Wik[qak,qzk]

        # Computation of final Input and Ouput EMAC
        for i in range(n):
            # Weight for emaco
            sumo = 0;
            for la in range(m):
                sumo = emaco[xc,la]*abs(modes_raw[la,xc])**2+sumo;

            #Weight for emaci
            sumi = 0;
            for lw in range(r):
                sumi = emaci[xc,lw]*abs(inm[xc,lw])**2+sumi 

            emacof[i] = sumo/((modes_raw[:,xc].T*modes_raw[:,xc]))  # emacof is the final output EMAC
            emacif[i] = sumi/(inm[xc,:]*inm[xc,:].T)                # emacif is the final input EMAC
            emac[i] = emacof[i]*emacif[i] 

    # Add the input EMAC, output EMAC, and MPC to the matrix freqdamp
    for lih in range(size(freqdmp)[0]):
        freqdmp[lih,4] = emacif(freqdmp[lih,2])
        freqdmp[lih,5] = emacof(freqdmp[lih,2])
        freqdmp[lih,6] = mpc[freqdmp(lih,3)]
        if freqdmp[lih,5]>0.5 and freqdmp[lih,7]>0.5:
            validationm=' valid'
        else:
            validationm=' not valid'

if __name__ == "__main__":
    from pathlib import Path
    import quakeio
    # import ssid as si
    import OKID
    import numpy as np
    channels = dict( # PAINTER RIO DELL
        inputs  = [17, 3, 20],
        outputs = [ 9, 7 , 4]
    )
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
