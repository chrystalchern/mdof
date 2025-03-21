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


from .construct import form_observability
def OutputEMAC(A,C,**options):
    """
    :param outlook:         number of timesteps for which to consider temporal consistency.
                            default: 100
    """
    p,n = C.shape
    assert A.shape == (n,n)
    no = options.get("outlook",
         options.get("no",
         options.get("horizon",
                     100)))
    Observability = form_observability(A,C,no)
    """Output EMAC (Eqs. 3.88-3.89)"""
    vals, vecs = np.linalg.eig(A)
    Phi_final = Observability[-p:,:] @ vecs                  # identified modal observability at the last timestep (last block row)
    Phi_final_hat = C @ vecs @ np.diag(vals**(no-1))         # expected modal observability at the last timestep
    assert Phi_final.shape == Phi_final_hat.shape == (p,n)
    emaco  = EMAC_Matrix(Phi_final,Phi_final_hat)
    return EnergyCondensedEMAC(emaco,(C@vecs).T)             # DIM: nxp, nxp


def MPC(A,C):
    """a) Modal Phase Collinearity (MPC) [Eqs. 3.85-3.87]"""
    _, vecs = np.linalg.eig(A)
    _,n = C.shape
    assert vecs.shape == (n,n)
    modes_raw = C@vecs
    s11, s22, s12 = np.zeros((3,n))
    nu, mpc = np.zeros((2,n))
    lam = np.zeros((2,n))
    for i in range(n):
        mode_i = modes_raw[:,i]
        s11[i] = np.real(mode_i).conjugate().transpose()@np.real(mode_i)
        s22[i] = np.imag(mode_i).conjugate().transpose()@np.imag(mode_i)
        s12[i] = np.real(mode_i).conjugate().transpose()@np.imag(mode_i)
        if s12[i] == 0:
            nu[i] = np.nan
        else:
            nu[i]= (s22[i]-s11[i])/(2*s12[i])
        lam[0,i] = (s11[i]+s22[i])/2 + s12[i]*np.sqrt(nu[i]**2+1)
        lam[1,i] = (s11[i]+s22[i])/2 - s12[i]*np.sqrt(nu[i]**2+1)
        mpc[i]   = ((lam[0,i]-lam[1,i])/(lam[0,i]+lam[1,i]))**2
    return mpc


def eigenfilter(A, test, by_index=False, return_indices=False):
    """
    Filters away undesirable, e.g., unstable, modes of the matrix `A`.
    1. Finds the eigenvalues `val` of `A` for which `test` evaluates to True.
    2. Recomputes a "filtered" matrix, `filtered_A`, with those eigenvalues set to zero.
    Test is evaluated on the eigenvalue (`test(val)`) by default (`by_index=False').
    Test is evaluated on the index of the eigenvalue if `by_index=True`.
    """
    vals,vecs = np.linalg.eig(A)

    indices = []
    for i,val in enumerate(vals):
        if test(i if by_index else val):
            indices.append(i)
            vals[i] = 0

    filtered_A = vecs @ np.diag(vals) @ np.linalg.inv(vecs)

    if return_indices:
        return filtered_A, indices
    else:
        return filtered_A


import inspect
def test_string(test):
    test_source_string = inspect.getsource(test)
    if "lambda" in test_source_string:
        return "\n".join(test_source_string.split(":")[-1].split("\n")[:-1])
    else:
        return "\n".join(test_source_string.split("return")[-1].split("\n")[:-1])
    

def stabilize_discrete(A, verbose=False, list_filtered_modes=False):
    # if test is True, the mode is filtered out.
    test = lambda v: np.abs(v) > 1
    filtered_A, indices = eigenfilter(A, test, return_indices=True)
    real_mode_indices = np.unique([i//2 for i in indices]).tolist()
    if verbose:
        print(f"""
        removing mode indices {indices} (real mode indices {real_mode_indices})
        because magnitude of eigenvalue is greater than 1 ({test_string(test)})
        """)
    if list_filtered_modes:
        return filtered_A, real_mode_indices
    else:
        return filtered_A


def stabilize_continuous(A, verbose=False):
    # if test is True, the mode is filtered out.
    test = lambda v: np.real(v) > 0
    filtered_A, indices = eigenfilter(A, test, return_indices=True)
    if verbose:
        real_mode_indices = np.unique([i//2 for i in indices]).tolist()
        print(f"""
        removing mode indices {indices} (real mode indices {real_mode_indices})
        because real part of eigenvalue is positive ({test_string(test)})
        """)
    return filtered_A


def filter_consistent(A, C, verbose=False, **options):
    # energy condensed output EMAC (extended modal amplitude coherence)
    energy_condensed_emaco = OutputEMAC(A,C)
    threshold = options.get('threshold',0.5)
    # if test is True, the mode is filtered out.
    def test(i):
        return energy_condensed_emaco[i] < threshold
    filtered_A, indices = eigenfilter(A, test, by_index=True, return_indices=True)
    if verbose:
        real_mode_indices = np.unique([i//2 for i in indices]).tolist()
        print(f"""
        removing mode indices {indices} (real mode indices {real_mode_indices})
        because EMAC < {threshold} ({test_string(test)})
        """)
    return filtered_A
    

def filter_collinear(A, C, verbose=False, **options):
    # MPC (modal phase collinearity)
    mpc = MPC(A,C)
    threshold = options.get('threshold',0.5)
    # if test is True, the mode is filtered out.
    def test(i):
        return mpc[i] < threshold
    filtered_A, indices = eigenfilter(A, test, by_index=True, return_indices=True)
    if verbose:
        real_mode_indices = np.unique([i//2 for i in indices]).tolist()
        print(f"""
        removing mode indices {indices} (real mode indices {real_mode_indices})
        because MPC < {threshold} ({test_string(test)})
        """)
    return filtered_A

        