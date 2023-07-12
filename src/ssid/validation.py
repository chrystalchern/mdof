import numpy as np
from numpy import pi

def EnergyCondensedEMAC(emac,phi):
    n,p = emac.shape
    assert emac.shape == phi.shape
    "Equation 3.93"
    return np.array([
        sum(emac[i,j]*abs(phi[i,j])**2 
            for j in range(p)
        )/np.real(phi[i,:].conjugate().transpose()@phi[i,:])
        for i in range(n) 
    ])

def EMAC_Matrix(Phi_final, Phi_final_hat):
    p,n = Phi_final.shape
    emac = np.zeros((n,p))
    for i in range(n):
        for j in range(p):
            Rij = min(
                (abs(Phi_final[j,i])/abs(Phi_final_hat[j,i])),
                (abs(Phi_final_hat[j,i])/abs(Phi_final[j,i]))
            )
            Pij = np.angle(Phi_final[j,i]/Phi_final_hat[j,i])
            Wij = max(
                1 - abs(Pij)/(pi/4),
                0
            )
            emac[i,j] = Rij*Wij
    return emac

# nt = number of timesteps (for OKID-ERA-DC and OKID-ERA, number of Markov parameters)
# Psi = eigenvectors of A, as columns of a matrix. can be reused from modal.system_modes().
# Gam = eigenvalues of A, as items of a vector. can be reused from modal.system_modes().
# Observability matrix can be reused from realize.srim().
# p = number of outputs
# n = model order (number of system variables)
def OutputEMAC(A,C,nt,Observability=None,Psi=None,Gam=None):
    p,n = C.shape
    assert A.shape == (n,n)
    """Output EMAC (Eqs. 3.88-3.89)"""
    if Gam is None:
        assert Psi is None
        from ssid.modal import condeig 
        Psi,Gam,_ = condeig(A)
    if Observability is None:
        Observability = np.empty((nt,p,n))
        Observability[0,:,:] = C
        A_pwr = A
        for pwr in range(1,nt):
            Observability[pwr,:,:] =  C@A_pwr
            A_pwr = A@A_pwr
        Observability = Observability.reshape((nt*p,n))
    Phi_final = Observability[-p:,:]@Psi                    # identified modal observability at the last timestep (last block row)
    Phi_final_hat = C@Psi@np.diag(Gam**(nt-1))              # expected modal observability at the last timestep
    assert Phi_final.shape == Phi_final_hat.shape == (p,n)
    emaco  = EMAC_Matrix(Phi_final,Phi_final_hat)
    return EnergyCondensedEMAC(emaco,(C@Psi).T)             # DIM: nxp, nxp

def MPC(A,C,Psi=None):
    """a) Modal Phase Collinearity (MPC) [Eqs. 3.85-3.87]"""
    if Psi is None:
        from ssid.modal import condeig 
        Psi,_,_ = condeig(A)
    _,n = C.shape
    assert Psi.shape == (n,n)
    modes_raw = C@Psi
    s11, s22, s12 = np.zeros((3,n))
    nu, mpc = np.zeros((2,n))
    lam = np.zeros((2,n))
    for i in range(n):
        mode_i = modes_raw[:,i]
        s11[i] = np.real(mode_i).conjugate().transpose()@np.real(mode_i)
        s22[i] = np.imag(mode_i).conjugate().transpose()@np.imag(mode_i)
        s12[i] = np.real(mode_i).conjugate().transpose()@np.imag(mode_i)
        nu[i]    = (s22[i]-s11[i])/(2*s12[i])
        lam[0,i] = (s11[i]+s22[i])/2 + s12[i]*np.sqrt(nu[i]**2+1)
        lam[1,i] = (s11[i]+s22[i])/2 - s12[i]*np.sqrt(nu[i]**2+1)
        mpc[i]   = ((lam[0,i]-lam[1,i])/(lam[0,i]+lam[1,i]))**2
    return mpc