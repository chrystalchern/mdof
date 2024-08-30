"""
This module implements functions that extract modal information
from a state space realization or transfer function.
"""
import numpy as np
import scipy.linalg as sl
from numpy import pi
from .validation import OutputEMAC, MPC
from .numerics import _condeig
    
def system_modes(realization, dt, n_peaks=None, sorted_by=None, sort_descending=True, filter_by=None, filter_lim=None, **options):
    """
    Modal identification from a state space system realization.

    :param realization:     realization in the form of state space coefficients ``(A,B,C,D)``
    :type realization:      tuple of arrays
    :param dt:              timestep.
    :type dt:               float
    :param n_peaks:         (optional) number of peaks to return. By default, all peaks are returned.
    :type n_peaks:          int
    :param sorted_by:       (optional) characteristic by which to sort peaks.
                            choose out of the following: 'period', 'MPC', 'EMAC', ...
                            sorting order is highest to lowest (descending), unless `sort_descending` is set to `False`.
    :type sorted_by:        string
    :param sort_descending: (optional) sort from highest to lowest
                            default: `True`
    :type sort_descending:  boolean
    :param filter_by:       (optional) list of characteristics by which to filter peaks.
                            example: ['MPC', 'EMAC'].
                            choose out of the following: 'period', 'MPC', 'EMAC', 'frequency', 'damping', 'condition_number', ...
    :type filter_by:        list of strings
    :param filter_lim:      (optional) list of lower and upper limits by which to filter.
                            example: [(0.5,1.0), (0.5,1.0)]
    :type filter_lim:       list of tuples of floats
    :param decimation:      (optional) decimation factor. default: 1
    :type decimation:       int
    :param Observability:   (optional) Observability matrix; can be reused from :func:`mdof.realize.srim`.
                            default: None
    :type Observability:    array

    :return:                system modes, including natural frequencies, damping ratios, mode shapes,
                            condition numbers, and modal validation metrics EMAC and MPC.
    :rtype:                 dictionary
    """
    decimation = options.get("decimation",
                             1)
    Observability = options.get("Observability",
                                None)
    
    dt = dt*decimation

    A,_,C,_ = realization
    # eigendecomp A
    Psi,Gam,cnd = _condeig(A)  # eigenvectors (Psi) & eigenvalues (Gam) of the matrix A

    # get damping and frequencies from eigendecomp of A
    Lam = (np.log(Gam))/dt
    Omega = np.abs(Lam) # radians per second.
    freq = Omega/(2*pi) # cycles per second (Hz)
    damp = -np.real(Lam)/Omega

    # get modeshapes from C and eigendecomp of A
    modeshape = C@Psi

    # energy condensed output EMAC (extended modal amplitude coherence)
    energy_condensed_emaco = OutputEMAC(A,C,**options)

    # MPC (modal phase collinearity)
    mpc = MPC(A,C)

    # for perfect data, the modes of the state space model come in pairs.
    # each pair's corresponding eigenvalues and eigenvectors are complex conjugates.
    # this means that we need to weed out unique roots: 
    # get indices of (1) roots that only show up once, and (2) the first of each pair.
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
            for i in range(len(freq)) if i not in notroots
            }

    if filter_by is not None:
        for i,filterkey in enumerate(filter_by):
            try:
                lim = filter_lim[i]
            except IndexError as e:
                print(e)
                print("Modes not filtered. filter_by and filter_lim must have matching lengths")
                return modes

            filterkey = filterkey.lower()
            if filterkey=='emac':
                filterkey = 'energy_condensed_emaco'
            elif filterkey=='period':
                filterkey = 'freq'
                lim = (1/lim[1], 1/lim[0])
            elif filterkey=='frequency':
                filterkey = 'freq'
            elif filterkey=='damping':
                filterkey = 'damp'
            elif filterkey=='condition_number':
                filterkey = 'cnd'

            try:
                modes = {k:v
                        for k,v in modes.items()
                        if v[filterkey]>=lim[0] and v[filterkey]<=lim[1]}
            except IndexError as e:
                print(e)
                print("Modes not filtered. filter_lim must be a list of tuples of two floats each.")
                return modes
            except TypeError as e:
                print(e)
                print("Modes not filtered. filter_lim must be a list of tuples of two floats each.")
                return modes
            except KeyError as e:
                print(e)
                print("Modes not filtered. Invalid filter_by.")
                return modes

    if sorted_by is not None:
        sorted_by = sorted_by.lower()
        if sorted_by=='period':
            sorted_by = 'freq'
            sort_descending = False
        sorted_modes = sorted(modes.items(), key=lambda mode: mode[1][sorted_by], reverse=sort_descending)
        modes = {k:v for k,v in sorted_modes}
    if n_peaks is not None:
        modes = list(modes.items())[:n_peaks]
        modes = {k:v for k,v in modes}

    return modes

    
def spectrum_modes(periods, amplitudes, n_peaks=None, sorted_by=None, **options):
    """
    Modal identification from a transfer function.

    :param periods:     transfer function periods
    :type periods:      array
    :param amplitudes:  transfer function amplitudes
    :type amplitudes:   array
    :param n_peaks:     (optional) number of peaks to return. By default, all peaks are returned.
    :type n_peaks:      int
    :param sorted_by:   (optional) characteristic by which to sort peaks. choose out of the following: 'height', ...
                        sorting order is highest to lowest (descending).
    :type sorted_by:    string
    :param prominence:  (optional) prominence of selected peaks
    :type prominence:   float

    :return:            (fundamental_periods, fundamental_amplitudes)
    :rtype:             tuple
    """
    from scipy.signal import find_peaks
    # height = options.get("height", 0.4)
    # width = options.get("width", 0.2)
    # rel_height = options.get("rel_height", 0.1)
    prominence = options.get("prominence", max(amplitudes)*0.3)
    
    # peak_indices, _ = find_peaks(amplitudes, height=height, width=width, rel_height=rel_height, prominence=prominence)
    peak_indices, _ = find_peaks(amplitudes, prominence=prominence)

    if sorted_by=='height':
        peak_indices = sorted(peak_indices, key=lambda peak: amplitudes[peak], reverse=True)
    fundamental_periods = periods[peak_indices]
    fundamental_amplitudes = amplitudes[peak_indices]
    
    if n_peaks is None:
        return (fundamental_periods,fundamental_amplitudes)
    else:
        return (fundamental_periods[:n_peaks],fundamental_amplitudes[:n_peaks])
    