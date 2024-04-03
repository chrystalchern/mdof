"""

.. currentmodule:: mdof

"""

# default imports

from .system import system
from .modal import system_modes, _condeig
from .transform import fdd

# from control.matlab import impulse as _impulse
# def impulse(system, t, **kwds):
#     from control import ss
#     dt = t[1] - t[0]
#     a,t = _impulse(ss(*system, dt), t, **kwds)
#     return a.squeeze()*dt,t


def eigid(inputs, outputs, **options):
    """
    System eigenvalues ``Gam`` and eigenvectors ``Psi`` identification from ``input`` and ``output`` arrays.

    :param inputs:      input time history. dimensions: :math:`(q,nt)`, where
                        :math:`q` = number of inputs, and :math:`nt` = number of timesteps
    :type inputs:       array
    :param outputs:     output response history.
                        dimensions: :math:`(p,nt)`, where :math:`p` = number of outputs, and
                        :math:`nt` = number of timesteps
    :type outputs:      array
    :param dt:          timestep.
    :type dt:           float
    :param decimation:  decimation factor. default: 1
    :type decimation:   int, optional

    :return:            (``Gam``, ``Psi``)
    :rtype:             tuple of 1D array, ND array
    """
    A,_,_,_ = sysid(inputs, outputs, **options)
    decimation = options.get("decimation", 1)
    # eigendecomp A
    Psi,Gam,_ = _condeig(A)  # eigenvectors (Psi) & eigenvalues (Gam) of the matrix A
    return (Gam,Psi)

    
def outid(outputs, dt, **options):
    """
    Fundamental frequencies ``F`` and modeshapes ``Phi`` identification from ``output`` array.

    :param outputs:     output response history.
                        dimensions: :math:`(p,nt)`, where :math:`p` = number of outputs, and
                        :math:`nt` = number of timesteps
    :type outputs:      array
    :param dt:          timestep.
    :type dt:           float

    :return:            (``F``, ``Phi``)
    :rtype:             tuple of 1D array, ND array
    """

    frequencies,U,S = transform.fdd(outputs=outputs, step=dt)

    import numpy as np
    p,nt = outputs.shape
    F = np.empty(p)

    from scipy.signal import find_peaks
    for i in range(p):
        amplitudes = S[i,:]
        prominence = options.get("prominence", max(amplitudes)*0.3)
        peaks, _ = find_peaks(amplitudes, prominence=prominence)
        F[i] = frequencies[peaks[0]]
        if i==0:
            Phi = U[:,:,peaks[0]]
    
    return (F,Phi)
    

def modes(inputs, outputs, dt, **options):
    """
    Modal identification using state space system realization from ``input`` and ``output`` arrays.

    :param inputs:      input time history. dimensions: :math:`(q,nt)`, where
                        :math:`q` = number of inputs, and :math:`nt` = number of timesteps
    :type inputs:       array
    :param outputs:     output response history.
                        dimensions: :math:`(p,nt)`, where :math:`p` = number of outputs, and
                        :math:`nt` = number of timesteps
    :type outputs:      array
    :param dt:          timestep.
    :type dt:           float
    :param method:      system identification method. default is "srim", other options are
                        "okid-era" and "okid-era-dc".
    :type method:       string, optional
    :param decimation:  decimation factor. default: 8
    :type decimation:   int, optional

    :return:            system modes, including natural frequencies, damping ratios, mode shapes,
                        condition numbers, and modal validation metrics EMAC and MPC.
    :rtype:             dictionary
    """
    realization = system(inputs, outputs, **options)

    return modal.system_modes(realization, dt)

    
def sysid(inputs, outputs, **options):
    """
    State space system realization (``A``, ``B``, ``C``, ``D``) from ``input`` and ``output`` arrays.

    :param inputs:      input time history. dimensions: :math:`(q,nt)`, where
                        :math:`q` = number of inputs, and :math:`nt` = number of timesteps
    :type inputs:       array
    :param outputs:     output response history.
                        dimensions: :math:`(p,nt)`, where :math:`p` = number of outputs, and
                        :math:`nt` = number of timesteps
    :type outputs:      array
    :param method:      system identification method. default is "srim", other options are "okid-era" and "okid-era-dc".
    :type method:       string, optional
    :param decimation:  decimation factor. default: 8
    :type decimation:   int, optional

    :return:            system realization in the form of state space coefficients (``A``, ``B``, ``C``, ``D``)
    :rtype:             tuple of arrays
    """
    return system(inputs, outputs, **options)