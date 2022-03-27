# [freqdmp, modeshape, RMSEpred, markovParamError] = 
def OKID_ERA_DC(dati, dato, dt, config):
#-----------------------------------------------------------------
# PART 2: OKID-ERA-DC (Observer Kalman filter Identification -
#         Eigen Realization with Direct Correlations)
#-----------------------------------------------------------------

# # Description of the Methodology:
#
# More information on OKID-ERA-DC method can be found in
# Section 3.4.6 of (Arici & Mosalam, 2006). Equations below refer to this report.
# OKID-ERA-DC is a Multiple Input - Multiple Output (MIMO) System Identification (SI) method
# that consists of the following steps:
# 1. Data pre-processing (baseline correction, filtering & decimation)
# 2. Identify observer Kalman filters using Observer Kalman filter Identification (OKID)
#    methodology. This step basically consists of developing a simple observer model of
#    a system from a MIMO ARX structure, Eq. (3.76), which is broken into these 6 steps:
#
# 2a. Determine Markov parameters (M) in a least squares sense, Eq. (3.76).

# 3b. Establish the Hankel matrix (H) from the Markov parameters, (Eq. 3.80).
# 3c. Use H to compute system matrices A, B & C, in which modal information is embedded.
# 3d. Obtain the modal information from matrices A, B & C.
# 3e. Spatial & temporal validation of the identified modes.
# 3f. Back calculate (estimate) the output accelerations with the state-space system &
#     check against the actual output accelerations.
#
# # Inputs
#
# In terms of input and output channels, five cases can be considered as
# described below.

# Note that the input and output channels are common to OKID-ERA-DC and
# SRIM.
#
# Modelparameters
#div = 1;     # A parameter used for decimating data. 1 uses entire data without downsampling.
#mro = 10;    # Model reduction order
#orm = 4;     # Order of the model. # of computed and plotted modes dep on orm.
#For orm = 2, one mode is found, for orm = 4, two modes are found.

# For case 1, one mode is transverse and the other is torsion.
# For all other cases, the second mode is a higher mode.
# Sometimes higher orm still gives fewer modes, e.g. orm = 8 for case 1 gives
# three modes, but one of them is invalid according to the EMAC & MPC criteria.
# kmax = 100;  #Number of computed Markov parameters, indicated as 1000 on page 43 of
# (Arici & Mosalam, 2006). However, it was input as 100 in the code.
# kmax = 100 runs much faster & both kmax = 100 & 1000 give the same results.

# Important output variables:
#  1. freqdamp variable is a matrix that includes the information of identified
#     frequencies, damping ratios & validation of the modes with MPC & EMAC criteria
#     Each row of freqdamp corresponds to a mode. Columns are as follows:
#     1)frequency, 2)damping ratio, 3)order index, 4)condition number, 5)input EMAC,
#     6)output EMAC, 7)MPC. If values in columns 5-7 are > 0.5, identified mode is valid.
#  2. modeshape stores the mode shape information for identified modes.
#  3. RMSEpred: root mean square error of the predicted output from
#     identified parameters with respect to the actual output.
#  4. Markovparamerror: root mean square error used to validate accurate
#     computation of Markov parameters

dn  = config.dn;
mro = config.mro;
orm = config.orm;
kmax = config.kmax;
verbose = true;

#d  = size(dat,1); # total number of time steps
#nc = size(dat,2);

## 2a. Obtain Observer Markov Parameters

# Note that main Step 2 develops Eq. 3.76. Therefore, it is not part of the code.
# Accordingly, the code continues with Step 2a to compute Observer Markov parameters in Eq. 3.76.
# The Markov parameters are computed in two steps:
# i) The Observer Markov matrices are computed from Eq. 3.76 using a linear regression approach
# ii) Compute the system Markov parameters from the Markov matrices using recursive relations

# defining several parameters
temsiz = size(dato);
temsizi = size(dati);
m = size(dato,2)       # number of output channels
l = size(dato,1)       ## of rows of dato = # of rows of dati = # of acceleration samples per channel
r = size(dati,2)       # number of input channels
# ASSERT SIZE(DATO,1) == SIZE(DATI,1)

n = orm                # assign order of model input to variable n, consistent with Eqs. 3.81-3.84
p = mro                # assign input model reduction order to variable p, consistent with Eq. 3.76

# Compute matrix U that represents ARX equation of current output on p time steps of past output
# & input values (Eq. 3.76)
U = zeros(((m+r)*p+r, l))
U[:r,:] = dati.T
for b in range(p+1):
    U[(b-2)*(r+m)+1+r:(b-2)*(r+m)+r+r+m, b:l] = \
            [[dati[:l-b+1,:r].T], [dato[:l-b+1, :m].T]]


# i) Compute the matrix of Observer Markov Parameter Matrix (M) 
#    in Eq 3.76 using Linear Regression
[uu,s,v] = svd(U,0);     # svd: Singular Value Decomposition
                         #      s is a diagonal matrix with the singular values
wr = diag(s)             # extract singular values from the diagonal matrix
pg = (r+m)*p+r
for lop in range((r+m)*p+r):
    if wr[lop] <= 0.001
        pg = lop;
        break



pss = v[:,:pg]@linsolve(s[1:pg,:pg], uu[:,:pg].T)
M = dato.T*pss           #M: Observer Markov Parameter Matrix

# Fit for multiple regression
ypreo = M@U

for i in range(m):
    temsump = sum((dato(:,i)-ypreo(i,:).T).^2)
    Jpre(i) = temsump/(sum(dato(:,i).^2))


markovParamError = sum(Jpre)/m;  # RMSE between actual output & y on left hand side in Eq. 3.76
# It should be quite small (e.g., 10^-3) for accurately computed Markow parameters


# ii) Compute Markov parameters (Y) using recursive relations in Eqs. 3.78 & 3.79

# Matrix D is directly equal to the Observer Markov parameter matrix (Eq. 3.77)
D = M[:, :r]  # D: Direct Transmission term, one of 4 system matrices of state space model

# Eqs. 3.31-3.34 define the four system matrices, A, B, C & D
Y[1] = D;
# First p steps (Eq. 3.78)
for ol in range(p):
    sumt = zeros((m,r))
    for lok in range(ol):
        sumt = M[:, r+(lok-1)*(r+m)+1+r:r+(lok-1)*(r+m)+r+m]*Y[ol-lok+1] + sumt

    Y[ol+1] = M[:, r+(ol-1)*(r+m)+1:r+(ol-1)*(r+m)+r] + sumt


# From p+1 to rest (Eq. 3.79)
for ol in range(p+1, dn+kmax):
    sumt = zeros((m,r))
    for lok in range(p):
        sumt += M[:,r+(lok-1)*(r+m)+1+r:r+(lok-1)*(r+m)+r+m]*Y[ol-lok+1]

    Y[ol+1] = +sumt

# Now, the Markow parameters Y have been computed.

## 2b. Establish Hankel matrix (H) from the Markov parameters (Y) (Eq. 3.80)

# psz = 1000;   psz is replaced by kmax and used as an imput parameter
# Obtain Hankel Matrix of Zeroth Order & First Order
for hj in range(kmax):
    for jh in range(l):
        H0[(hj-1)*m+1:hj*m, (jh-1)*r+1:jh*r] = Y[jh+hj]
        H1[(hj-1)*m+1:hj*m, (jh-1)*r+1:jh*r] = Y[jh+hj+1]



## 2c. Use H matrix to compute system matrices A, B & C

R1,Sis,S1 = svd(H0) # singular value decomposition


# A: state transition matrix (Eqs. 3.32 & 3.82)
A  = Sis[:n,:n]**(-0.5)*R1[:,:n].T*H1*S1[:,:n]*Sis[:n,:n]**(-0.5)
Qb = Sis[:n,:n]**0.5*S1[:,:n].T
# B: input influence matrix (Eqs. 3.32 & 3.83)
B  = Qb[:,:r]
Pb = R1[:,:n]*Sis[:n,:n]**0.5 
# C: output influence matrix (Eqs. 3.34 & 3.84)
C  = Pb[:m,:]















## 2d. Obtain the modal information from the system matrices A, B & C
# This includes determination of: a) modal frequencies, b) damping ratios & c) mode shapes
[freqdmp, modeshape, sj1, v, d] = ExtractModes(dt,A,B,C,D);
modes_raw = C*v;
inm = v\B;               #initial modal contribution

## 2e. Validation Analysis
# Two criteria are used for selection of identified genuine modes, in terms of spatial & temporal consistency.
# a) Modal Phase Collinearity (MPC) testing spatial consistency of identification results.
#    Modes having MPC value above 0.5 (mpc parameter below) are considered as genuine modal quantities.
# b) Exted Modal Amplitude Coherence (EMAC), evaluates temporal consistency of the identification results.
#    Both output EMAC & input EMAC can be computed. Input EMAC requires the controllability matrix.
#    Because the controllability matrix is not estimated by all considered SI methods,
#    this criterion is computed, but not used.
#    Modes with output EMAC values < 0.5 are considered spurious & therefore not reported.

# a) Modal Phase Collinearity (MPC) [Eqs. 3.85-3.87]
for q in range(n):
    sxx(:,q) = real(modes_raw(:,q)).T*real(modes_raw(:,q));
    syy(:,q) = imag(modes_raw(:,q)).T*imag(modes_raw(:,q));
    sxy(:,q) = real(modes_raw(:,q)).T*imag(modes_raw(:,q));
    nu(q) = (syy(:,q)-sxx(:,q))/(2*sxy(:,q));
    lam(1,q) = (sxx(:,q)+syy(:,q))/2+sxy(:,q)*(nu(q)^2+1)^0.5;
    lam(2,q) = (sxx(:,q)+syy(:,q))/2-sxy(:,q)*(nu(q)^2+1)^0.5;
    mpc(q) = ((lam(1,q)-lam(2,q))/(lam(1,q)+lam(2,q)))^2;


# b) Exted Modal Amplitude Coherence (EMAC)
v_inv = v^-1;
qlin = v_inv*Qb;     # Controllability Matrix used for the input-EMAC
plin = Pb*v;     # Observability Matrix used for the output-EMAC

lamb = v_inv*A*v
bkh = v_inv*B;
for hn in range(n):
    for ll = 0:l-1
        qhat(hn,ll*r+1:ll*r+r) = bkh(hn,:)*(lamb(hn,hn))^ll;


selsiz = min(size(qlin),size(qhat));

for hnd in range(n):
    ql = qlin[hnd,:selsiz(2)]
    qh = qhat[hnd,:selsiz(2)]
    mac(hnd) = abs(ql*qh.T)/(abs(ql*ql.T)*abs(qh*qh.T))**0.5


# Output EMAC (Eqs. 3.88-3.89)
# Pick the last block row
pto = plin[(kmax-1)*m+1:m*kmax,:]  # the identified value at T0
for ds in range(n):
    ptop[:,ds] = modes_raw(:,ds)*exp(sj1(ds)*dt*(kmax-1));


# Computation of rij
for qa in range(n):
    for qz in range(m):
        Rij[qa,qz] = min((abs(pto(qz,qa))/abs(ptop(qz,qa))),(abs(ptop(qz,qa))/abs(pto(qz,qa))));
        Pij = angle(pto(qz,qa)/ptop(qz,qa))
        Pijn[qa,qz] = Pij
        if abs(Pij) <= pi/4
            Wij[qa,qz] = 1 - abs(Pij)/(pi/4)
        else
            Wij[qa,qz] = 0

        emaco[qa,qz] = Rij[qa,qz]*Wij[qa,qz];      # emaco is the ouput emac



# Input EMAC
# Pick the last block column
qto = qlin[:, (l-1)*r+1:l*r]
qtop = d**(l-1)*inm;

# EMAC Input Variation
for er in range(l):
    qtovar = qlin(:,(er-1)*r+1:er*r);
    qtopvar = d^(er-1)*inm;
    # Computation of rik
    for qak in range(n):
        for qzk in range(r):
            Rik(qak,qzk) = min((abs(qtovar(qak,qzk))/abs(qtopvar(qak,qzk))),(abs(qtopvar(qak,qzk))/abs(qtovar(qak,qzk))));
            Pik = angle(qtovar(qak,qzk)/qtopvar(qak,qzk));
            if abs(Pik)<=pi/4
                Wik(qak,qzk) = 1-abs(Pik)/(pi/4);
            else
                Wik(qak,qzk) = 0;

            emaci(qak,qzk) = Rik(qak,qzk)*Wik(qak,qzk);


    # Weight for emaci
    emacif = zeros(n,1);
    for xc in range(n):
        sumi = 0;
        for lw in range(r):
            sumi = emaci(xc,lw)*(inm(xc,lw)*inm(xc,lw).T)+sumi;

        emacif(xc) = sumi/(inm(xc,:)*inm(xc,:).T);

    emacivar(:,er) = emacif;


# Computation of rik
for qak in range(n):
    for qzk in range(r):
        Rik(qak,qzk) = min((abs(qto(qak,qzk))/abs(qtop(qak,qzk))),(abs(qtop(qak,qzk))/abs(qto(qak,qzk))));
        Pik = angle(qto(qak,qzk)/qtop(qak,qzk));
        if abs(Pik)<=pi/4
          	Wik(qak,qzk) = 1-abs(Pik)/(pi/4);
       	else
            Wik(qak,qzk) = 0;

        emaci(qak,qzk) = Rik(qak,qzk)*Wik(qak,qzk);



# Computation of final Input and Ouput EMAC
for i in range(n):
    # Weight for emaco
    sumo = 0;
    for la in range(m):
        sumo = emaco[xc,la]*abs(modes_raw(la,xc))**2+sumo;

    #Weight for emaci
    sumi = 0;
    for lw in range(r):
        sumi = emaci(xc,lw)*abs(inm(xc,lw))**2+sumi 

    emacof[i] = sumo/((modes_raw(:,xc).T*modes_raw(:,xc)))            #emacof is the final output EMAC
    emacif[i] = sumi/(inm(xc,:)*inm(xc,:).T)              #emacif is the final input EMAC
    emac[i] = emacof[i]*emacif[i] 


# Add the input EMAC, output EMAC, and MPC to the matrix freqdamp
for lih in range(size(freqdmp):,1):
    freqdmp(lih,5) = emacif(freqdmp(lih,3))
    freqdmp(lih,6) = emacof(freqdmp(lih,3))
    freqdmp(lih,7) = mpc(freqdmp(lih,3))
    if freqdmp(lih,6)>0.5 and freqdmp(lih,7)>0.5:
        validationm=' valid'
    else:
        validationm=' not valid'


