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
            # print(f"{Phi_final[j,i]=}, {Phi_final_hat[j,i]=}")
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


def OutputEMAC(A,C,Psi=None,Gam=None,**options):
    """
    :param outlook:         number of timesteps for which to consider temporal consistency. default: 100
    :param Psi:             eigenvectors of A, as columns of a matrix. can be reused from modal.system_modes().
    :param Gam:             eigenvalues of A, as items of a vector. can be reused from modal.system_modes().
    :param Observability:   Observability matrix; can be reused from :func:`mdof.realize.srim`.
    """
    # """
    # :param p:               number of outputs
    # :param n:               model order (number of system variables)

    # Examples
    # --------
    # >>> a = np.random.randn(9, 6) + 1j*np.random.randn(9, 6)
    # >>> b = np.random.randn(2, 7, 8, 3) + 1j*np.random.randn(2, 7, 8, 3)

    # Reconstruction based on full SVD, 2D case:

    # >>> U, S, Vh = np.linalg.svd(a, full_matrices=True)
    # >>> U.shape, S.shape, Vh.shape
    # ((9, 9), (6,), (6, 6))
    # >>> np.allclose(a, np.dot(U[:, :6] * S, Vh))
    # True
    # >>> smat = np.zeros((9, 6), dtype=complex)
    # >>> smat[:6, :6] = np.diag(S)
    # >>> np.allclose(a, np.dot(U, np.dot(smat, Vh)))
    # True

    # Reconstruction based on reduced SVD, 2D case:

    # >>> U, S, Vh = np.linalg.svd(a, full_matrices=False)
    # >>> U.shape, S.shape, Vh.shape
    # ((9, 6), (6,), (6, 6))
    # >>> np.allclose(a, np.dot(U * S, Vh))
    # True
    # >>> smat = np.diag(S)
    # >>> np.allclose(a, np.dot(U, np.dot(smat, Vh)))
    # True

    # Reconstruction based on full SVD, 4D case:

    # >>> U, S, Vh = np.linalg.svd(b, full_matrices=True)
    # >>> U.shape, S.shape, Vh.shape
    # ((2, 7, 8, 8), (2, 7, 3), (2, 7, 3, 3))
    # >>> np.allclose(b, np.matmul(U[..., :3] * S[..., None, :], Vh))
    # True
    # >>> np.allclose(b, np.matmul(U[..., :3], S[..., None] * Vh))
    # True

    # Reconstruction based on reduced SVD, 4D case:

    # >>> U, S, Vh = np.linalg.svd(b, full_matrices=False)
    # >>> U.shape, S.shape, Vh.shape
    # ((2, 7, 8, 3), (2, 7, 3), (2, 7, 3, 3))
    # >>> np.allclose(b, np.matmul(U * S[..., None, :], Vh))
    # True
    # >>> np.allclose(b, np.matmul(U, S[..., None] * Vh))
    # True
    # """
    p,n = C.shape
    assert A.shape == (n,n)
    no = options.get("outlook",
         options.get("no",
         options.get("horizon",
                     100)))
    Observability = options.get("Observability",
                                None)
    """Output EMAC (Eqs. 3.88-3.89)"""
    if Gam is None:
        assert Psi is None
        from mdof.modal import condeig 
        Psi,Gam,_ = condeig(A)
    if Observability is None:
        Observability = np.empty((no,p,n))
        Observability[0,:,:] = C
        A_pwr = A
        for pwr in range(1,no):
            Observability[pwr,:,:] =  C@A_pwr
            A_pwr = A@A_pwr
        Observability = Observability.reshape((no*p,n))
    Phi_final = Observability[-p:,:]@Psi                    # identified modal observability at the last timestep (last block row)
    Phi_final_hat = C@Psi@np.diag(Gam**(no-1))              # expected modal observability at the last timestep
    assert Phi_final.shape == Phi_final_hat.shape == (p,n)
    emaco  = EMAC_Matrix(Phi_final,Phi_final_hat)
    return EnergyCondensedEMAC(emaco,(C@Psi).T)             # DIM: nxp, nxp

def MPC(A,C,Psi=None):
    """a) Modal Phase Collinearity (MPC) [Eqs. 3.85-3.87]"""
    if Psi is None:
        from mdof.modal import condeig 
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
        # print(f"{mode_i=}, {s11[i]=}, {s22[i]=}, {s12[i]=}")
        if s12[i] == 0:
            nu[i] = np.nan
        else:
            nu[i]= (s22[i]-s11[i])/(2*s12[i])
        lam[0,i] = (s11[i]+s22[i])/2 + s12[i]*np.sqrt(nu[i]**2+1)
        lam[1,i] = (s11[i]+s22[i])/2 - s12[i]*np.sqrt(nu[i]**2+1)
        mpc[i]   = ((lam[0,i]-lam[1,i])/(lam[0,i]+lam[1,i]))**2
    return mpc