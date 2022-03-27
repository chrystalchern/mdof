# Standard library
import multiprocessing
from functools import partial
# Dependencies
import numpy as np
from tqdm import tqdm

from ExtractModes import *


linsolve = np.linalg.solve
lsqminnorm = lambda *args: np.linalg.lstsq(*args,rcond=None)[0]

def blk_3(i, CA, U):
    return i, np.einsum('kil,klj->ij', CA[:i,:,:], U[-i:,:,:])

def _srim(dati, dato, config, full=True):
    """
    # Description of the Methodology:
   
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
   
    For orm = 2, one mode is found, for orm = 4, two modes are found.
    For case 1, one mode is transverse & the other is torsion.
    For all other cases, the second mode is a higher mode.
    Sometimes higher orm still gives fewer modes, e.g. orm = 8 for case 1 gives
    three modes, but one of them is invalid according to the EMAC & MPC criteria.
    same orm in OKID-ERA-DC is used. It can be changed if needed.
   
    Important output variables:
     1. freqdampSRIM variable is a matrix that includes the information of identified
        frequencies, damping ratios & validation of the modes with MPC & EMAC criteria.
        Each row of freqdamp corresponds to a mode. Columns are as follows:
        1)frequency, 2)damping ratio, 3)order index, 4)condition number, 5)MPC.
        If values in columns 5 is > 0.5, identified mode is valid.
     2. modeshapeSRIM stores the mode shape information for identified modes.
     3. RMSEpredSRIM: root mean square error of the predicted output from
        identified parameters with respect to the actual output (currently
        commented out).
    """

    dn = config.dn
    to = config.to
    dt = to
    p = config.p     # # steps used for the identification. Referred to as the prediction horizon in literature
    n1 = config.orm  # Order of the model.


    # 2a. Compute y (output) and u (input) vectors (Eqs. 3.58 & 3.60)

    # Note that main Step 2 develops Eq. 3.57.
    # Therefore, it is not part of the code.
    # Accordingly, the code continues with Step 2a to compute the output & input vectors.

    # Calculate the usable size of the data matrix
    #dn = size(dat,1)/div;       # total # time steps after decimating
    nsizS = dn-1-p+2
    n = n1

    l,m = dato.shape
    _,r = dati.shape # r is the number of input channels

    ypS = np.zeros((r*p,nsizS))     
    upS = np.zeros((r*p,nsizS))

    # Compute y (output) & u (input) vectors (Eqs. 3.58 & 3.60)
    for b in range(p):
        ypS[b*m:(b+1)*m,:nsizS+1] = dato[b:nsizS+b, :].T
        upS[b*r:(b+1)*r,:nsizS+1] = dati[b:nsizS+b, :].T


    # 2b. Compute the correlation terms and the coefficient matrix 
    #     (Eqs. 3.68 & 3.69).

    # Compute the correlation terms (Eq. 3.68)
    Ryy = ypS@ypS.T/nsizS
    Ruu = upS@upS.T/nsizS
    Ruy = upS@ypS.T/nsizS

    # Compute the correlation matrix (Eq. 3.69)
    Rhh = Ryy - Ruy.T@linsolve(Ruu,Ruy)

    #
    # 2c. Obtain observability matrix using full or partial decomposition 
    #     (Eqs. 3.72 & 3.74).
    #
    if full:
        # Full Decomposition Method
        un,*_ = np.linalg.svd(Rhh,0)           # Eq. 3.74
        Op = un[:,:n1]                         # Eq. 3.72
        A = lsqminnorm(Op[:(p-1)*m,:], Op[m:p*m,:])
        C = Op[:m,:]
    else:
        # Partial Decomposition Method
        un,*_ = np.linalg.svd(Rhh[:,:(p-1)*m+1],0)
        Op = un[:,:n1]
        A = lsqminnorm(Op[:(p-1)*m,:], Op[m:p*m,:])
        C = un[:m,:]


    # Computation of system matrices B & D
    # Output Error Minimization

    # Setting up the Phi matrix
    Phi  = np.zeros((m*nsizS, n1+m*r+n1*r))
    CA_powers = np.zeros((nsizS, m, A.shape[1]))
    CA_powers[0, :, :] = C
    A_p = A
    for pwr in range(1,nsizS):
        CA_powers[pwr,:,:] =  C@A_p
        A_p = A@A_p

    # First block column of Phi
    for df in range(nsizS):
        Phi[df*m:(df+1)*m,:n1] = CA_powers[df,:,:]

    # Second block column of Phi
    Imm = np.eye(m)
    for i in range(nsizS):
        Phi[i*m:(i+1)*m,n1:n1+m*r] = np.kron(dati[i,:],Imm)

    # Third block column of Phi
    In1n1 = np.eye(n1)
    cc = n1+m*r + 1
    dd = n1+m*r+n1*r

    krn = np.array([np.kron(dati[i,:],In1n1) for i in range(nsizS)])

    with multiprocessing.Pool(6) as pool:
        for i,res in tqdm(
                pool.imap_unordered(
                    partial(blk_3,CA=CA_powers,U=np.flip(krn,0)),
                    range(1,nsizS),
                    200
                ),
                total = nsizS
            ):
            Phi[i*m:(i+1)*m,cc-1:dd] = res

    y = dato[:nsizS,:].flatten()

    teta = lsqminnorm(Phi,y)

    x0 = teta[:n1]
    dcol = teta[n1:n1+m*r]
    bcol = teta[n1+m*r:n1+m*r+n1*r]

    D = np.zeros((m,r))
    B = np.zeros((n,r))
    for wq in range(r):
        D[:,wq] = dcol[wq*m:(wq+1)*m]

    for ww in range(r):
        B[:,ww] = bcol[ww*n:(ww+1)*n]

    return A,B,C,D

#PY
# #% 2e. Obtain the modal information from the system matrices A & C
# # This includes determination of: a) modal frequencies, b) damping ratios & c) mode shapes
# # c) Determination of mode shapes
#     mod = C1@vS                 # mode shapes (Eq. 3.40), v is the eigenvectors of matrix A
# 

#% 2f. Validation Analysis

# Two criteria are used for selection of identified genuine modes, in terms of spatial & temporal consistency.
# a) Modal Phase Collinearity (MPC) testing spatial consistency of identification results.
#    Modes having MPC value above 0.5 (mpc parameter below) are considered as genuine modal quantities.
# b) Extended Modal Amplitude Coherence (EMAC), evaluates temporal consistency of the identification results.
#    Both output EMAC & input EMAC can be computed. Input EMAC requires the controllability matrix.
#    Because the controllability matrix is not estimated by all considered SI methods,
#    this criterion is computed, but not used.
#    Modes with output EMAC values < 0.5 are considered spurious & therefore not reported.

# a) Modal Phase Collinearity (MPC) [Eqs. 3.85-3.87]
#Py
##    for q in range(n):
##        a = real(mod[:,q])
##        b = imag(mod[:,q])
##        sxx[:,q] = a.T*a
##        syy[:,q] = b.T*b
##        sxy[:,q] = a.T*b
##        nu[q] = (syy[:,q]-sxx[:,q])/(2*sxy[:,q])
##        lam[1,q] = (sxx[:,q]+syy[:,q])/2+sxy[:,q]*(nu(q)**2+1)**0.5
##        lam[2,q] = (sxx[:,q]+syy[:,q])/2-sxy[:,q]*(nu(q)**2+1)**0.5
##        mpc[q] = ((lam[0,q]-lam[1,q])/(lam[0,q]+lam[1,q]))**2


# b) Extended Modal Amplitude Coherence (EMAC)

# Only EMAC Output is computed as there is no Controllability Matrix

# Note that the computations are commented out as the matrix B is needed

#%KKKKK
##PY
##    plin = Op1@vS                # Observability Matrix used for the output-EMAC
##    lamb = linsolve(vS,A1)*vS
##    bkh = linsolve(vS,B)
### Pick the last block row
##    pto = plin((p-1)*m+1:m*p,:)  # the identified value at T0
##    for ds in range(n):
##        ptop[:,ds] = mod[:,ds]*exp(sj1S(ds)*to*(p-1))
##
### Computation of rij
##    for qa in range(n):
##        for qz in range(m):
##            Rij(qa,qz) = min((abs(pto(qz,qa))/abs(ptop(qz,qa))),(abs(ptop(qz,qa))/abs(pto(qz,qa))))
##            Pij = angle(pto(qz,qa)/ptop(qz,qa))
##            Pijn(qa,qz) = Pij
##            if abs(Pij) <= pi/4:
##                Wij[qa,qz] = 1-abs(Pij)/(pi/4)
##            else:
##                Wij[qa,qz] = 0
##
##            emaco[qa,qz] = Rij[qa,qz]*Wij[qa,qz]
##
##
### Computation of final emac
##    for xc in range(n):
##        # Weight for emaco
##        sumo = 0.0
##        for la in range(m):
##            sumo = emaco(xc,la)*abs(mod(la,xc))**2+sumo
##        emacof[xc] = sumo/((mod[:,xc].T*mod[:,xc]))
##        emac[xc] = emaco[xc]
#%KKKKK

# Add the MPC to the matrix freqdampSRIM
#    for lih = 1:size(freqdmpSRIM,1)
#        freqdmpSRIM[lih,5] = emacof(freqdmpSRIM(lih,3))
#        freqdmpSRIM[lih,6] = mpc(freqdmpSRIM(lih,3))
#        if freqdmpSRIM[lih,5]>0.5 and freqdmpSRIM[lih,6]>0.5:
#            validationm = ' valid'
#        else:
#            validationm = ' not valid'
#
#        scroutput = strcat('Mode',num2str(lih), ...
#            ': Output EMAC =  ',num2str(freqdmpSRIM(lih,5)),...
#            ', MPC =  ',num2str(freqdmpSRIM(lih,6)),...
#            ' -->',' SRIM Identified Mode ',...
#            num2str(lih), ' is',validationm)
#        sprintf(scroutput)


#% 2g. Back calculate (estimate) output accelerations with state-space system &
#%     check against actual output accelerations

# Note that the computations are commented out as the matrix B is needed

# Prediction using state space model
#%KKKKK
##PY
##    ms1 = modstruc(A1,B,C1,D,zeros(n,m),x0)
##    th1 = ms2th(ms1,'d')
##    e,r = resid([dato dati],th1)
##    simy = idsim([dati],th1);                # simy represents the estimated accelerations
##
##    for i in range(m):
##        temsum = sum((dato[:,i]-simy[:,i]).**2)
##        Jm[i] = temsum/(sum(dato[:,i].**2));     # Root mean square error of estimated accelerations
##
##    RMSEpredSRIM = sum(Jm)/m
###%KKKKK
##    return freqdmpSRIM,modeshapeSRIM,RMSEpredSRIM

if __name__ == "__main__":
    import sys
    import quakeio
    from pathlib import Path
    channels = [[17, 3, 20], [9, 7, 4]]

    if len(sys.argv) < 2:
        data_dir = Path("RioDell_Petrolia_Processed_Data")

        first_input = quakeio.read(data_dir/f"CHAN{channels[0][0]:03d}.v2")
        npoints = len(first_input.accel.data)
        inputs, outputs = np.zeros((2,npoints,len(channels[0])))

        # Inputs
        inputs[:,0] = first_input.accel.data
        for i,inp in enumerate(channels[0][1:]):
            inputs[:,i+1] = quakeio.read(data_dir/f"CHAN{inp:03d}.V2").accel.data

        # Outputs
        for i,inp in enumerate(channels[1]):
            outputs[:,i] = quakeio.read(data_dir/f"CHAN{inp:03d}.V2").accel.data

        dt = first_input.accel["time_step"]

    else:
        event = quakeio.read(sys.argv[-1])
        inputs = np.array([
            event.at(file_name=f"CHAN{i:03d}.V2").accel.data for i in channels[0]
        ]).T
        outputs = np.array([
            event.at(file_name=f"CHAN{i:03d}.V2").accel.data for i in channels[1]
        ]).T
        npoints = len(inputs[:,0])
        dt = event.at(file_name=f"CHAN{channels[0][0]:03d}.V2").accel["time_step"]




    class T: pass
    config_srim = T()
    config_srim.p   =  5
    config_srim.to  = dt
    config_srim.dn  = npoints - 1
    config_srim.orm =  4
    A,B,C,D = _srim(inputs, outputs, config_srim)
    freqdmpSRIM, modeshapeSRIM, *_ = ExtractModes(dt, A, B, C, D)
    print(1/freqdmpSRIM[:,0])


