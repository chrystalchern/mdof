import numpy as np
from numpy import pi

def EnergyCondensedEMAC(n,m,emac,phi,debug=False):
    "Equation 3.93"
    return np.array([
        sum(emac[i,j]*abs(phi[i,j])**2 
            for j in range(m)
        )/np.real(phi[i,:].conjugate().transpose()@phi[i,:])
        for i in range(n) 
    ])

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



