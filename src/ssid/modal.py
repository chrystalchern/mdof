import numpy as np
import scipy.linalg as sl
from numpy import pi
from ssid.validation import OutputEMAC, MPC

def condeig(a): # TODO: make this match matlab source code for condeig
    """
    vals, vecs, conds = condeig(A) Computes condition numbers for the
    eigenvalues of a matrix. The condition numbers are the reciprocals
    of the cosines of the angles between the left and right eigenvectors.
    Inspired by Arno Onken's Octave code for condeig.

    https://github.com/macd/rogues/blob/master/rogues/utils/condeig.py
    """
    m, n = a.shape
    # eigenvalues, left and right eigenvectors
    lamr, vl, vr = sl.eig(a, left=True, right=True)
    vl = vl.T.conj()
    # Normalize vectors
    for i in range(n):
        vl[i, :] = vl[i, :] / np.sqrt(abs(vl[i, :] ** 2).sum())
    # Condition numbers are reciprocal of the cosines (dot products) of the
    # left eignevectors with the right eigenvectors.
    c = abs(1 / np.diag(np.dot(vl, vr))) 
    return vr, lamr, c

def system_modes(realization, dt, decimation=1, Observability=None, nt=None):

    dt = dt*decimation

    if nt is None:
        nt = 100

    A,_,C,_ = realization
    # eigendecomp A
    Psi,Gam,cnd = condeig(A)  # eigenvectors (Psi) & eigenvalues (Gam) of the matrix A

    # get damping and frequencies from eigendecomp of A
    Lam = (np.log(Gam))/dt
    # TODO: maybe clean
    Omega = np.real((Lam*np.conj(Lam))**0.5)  # radians per second. taking the real part because numpy keeps a +0j.
    # Omega = np.real_if_close((Lam*np.conj(Lam))**0.5) # use real_if_close
    freq = Omega/(2*pi) # cycles per second (Hz)
    Omega.__str__() # helps to avoid "divide by zero" numpy warning message
    # Lam_real = np.real(Lam)
    # Lam_real+=1
    # Lam_real-=1
    # Omega_real = np.real(Omega)
    # Omega_real+=1
    # Omega_real-=1
    # del Omega
    # del Lam
    # import gc
    # gc.collect()
    # print(Omega_real)
    # print(Lam_real)
    # damp = -Lam_real/Omega_real
    damp = -np.real(Lam)/Omega

    # get modeshapes from C and eigendecomp of A
    modeshape = C@Psi

    # energy condensed output EMAC (extended modal amplitude coherence)
    energy_condensed_emaco = OutputEMAC(A,C,nt,Observability=Observability,Psi=Psi,Gam=Gam)

    # MPC (modal phase collinearity)
    mpc = MPC(A,C,Psi=Psi)

    # weed out unique roots: get indices of (1) roots that only show up once, and
    # (2) the first of each pair.
    _, notroots = np.unique(freq.round(decimals=5), return_index=True)
    
    # print(notroots)
    modes = {str(i):
                {'freq': freq[i],  # identified frequency
                'damp': damp[i],   # identified damping ratio
                'modeshape': modeshape[:,i],  # identified modeshape
                'cnd': cnd[i],     # condition number of the eigenvalue
                'energy_condensed_emaco': energy_condensed_emaco[i],  # energy condensed output emac
                'mpc': mpc[i],     # MPC
                }
            for i in range(len(freq)) if i not in notroots # and energy_condensed_emaco[i] > 0.5 and mpc[i] > 0.5
            }

    return modes

def spectrum_modes(periods, amplitudes, nmodes=1):
    highest_amplitude_indices = np.argpartition(-amplitudes, range(nmodes))[:nmodes]  # TODO: peak picking algorithm
    fundamental_periods = periods[highest_amplitude_indices]
    fundamental_amplitudes = amplitudes[highest_amplitude_indices]
    return (fundamental_periods, fundamental_amplitudes)