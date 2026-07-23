"""

.. currentmodule:: mdof

"""

# default imports
from .system import system
from .realization import Realization
from . import transform
from . import modal
import numpy as np

__all__ = [
    # core object
    "Realization",
    # identification
    "sysid",
    "system",
    "modes",
    "outid",
    "eigid",
    # prediction
    "predict",
    # submodules
    "transform",
    "modal",
]


def outid(outputs, dt, **options):
    """
    Fundamental periods ``P`` and modeshapes ``Phi`` identification
    from ``output`` array.

    :param outputs:     output response history.
                        dimensions: :math:`(p,nt)`, where
                        :math:`p` = number of outputs, and
                        :math:`nt` = number of timesteps
    :type outputs:      array
    :param dt:          timestep.
    :type dt:           float

    :return:            (``P``, ``Phi``)
    :rtype:             tuple of 1D array, ND array
    """

    frequencies,U,S = transform.fdd(outputs=outputs, step=dt)
        
    p,n_spectrum_pts = S.shape
    P = np.empty(p)

    from scipy.signal import find_peaks
    for i in range(p):
        amplitudes = S[i,:]
        prominence = options.get("prominence", max(amplitudes)*0.3)
        peaks, _ = find_peaks(amplitudes, prominence=prominence)
        P[i] = 1/(frequencies[peaks[0]])
        if i==0:
            Phi = U[:,:,peaks[0]]
    
    return P,Phi
    

def modes(inputs, outputs, dt, **options):
    """
    Fundamental periods ``P`` and modeshapes ``Phi`` identification
    from ``input`` and ``output`` arrays.

    :param inputs:      input time history.
                        dimensions: :math:`(q,nt)`, where
                        :math:`q` = number of inputs, and
                        :math:`nt` = number of timesteps.
    :type inputs:       array
    :param outputs:     output response history.
                        dimensions: :math:`(p,nt)`, where
                        :math:`p` = number of outputs, and
                        :math:`nt` = number of timesteps
    :type outputs:      array
    :param dt:          timestep, in seconds.
    :type dt:           float
    :param method:      system identification method. default is
                        "srim", other options are "okid-era" and
                        "okid-era-dc".
    :type method:       string, optional

    :return:            (``P``, ``Phi``)
    :rtype:             tuple of 1D array, ND array
    """

    realization = system(inputs, outputs, dt=dt, **options)

    if options.get("stabilize", False):
        realization = realization.stabilize()

    modes = realization.modes()
    P = [1/mode['freq'] for mode in modes.values()]
    Phi = np.array([mode['modeshape'] for mode in modes.values()])

    return P,Phi


def eigid(inputs, outputs, **options):
    r"""
    State space system eigenvalues ``vals`` and eigenvectors ``vecs``
    identification from `input`` and ``output`` arrays. This is
    the eigendecomposition of the discrete system state transition
    matrix, :math:`\mathbf{A}`.

    :param inputs:      input time history. dimensions:
                        :math:`(q,nt)`, where :math:`q` = number of
                        inputs, and :math:`nt` = number of timesteps
    :type inputs:       array
    :param outputs:     output response history.
                        dimensions: :math:`(p,nt)`, where
                        :math:`p` = number of outputs, and
                        :math:`nt` = number of timesteps
    :type outputs:      array
    :param decimation:  decimation factor. default: no decimation
    :type decimation:   int, optional

    :return:            (``vals``, ``vecs``)
    :rtype:             tuple of 1D array, ND array
    """
    
    A,_,_,_ = sysid(inputs, outputs, **options)
    vals,vecs = np.linalg.eig(A)
    return vals,vecs


def sysid(inputs, outputs, dt=None, **options):
    """
    State space system realization (``A``, ``B``, ``C``, ``D``) from
    ``input`` and ``output`` arrays.

    :param inputs:      input time history. dimensions:
                        :math:`(q,nt)`, where :math:`q` = number of
                        inputs, and :math:`nt` = number of timesteps
    :type inputs:       array
    :param outputs:     output response history.
                        dimensions: :math:`(p,nt)`, where
                        :math:`p` = number of outputs, and
                        :math:`nt` = number of timesteps
    :type outputs:      array
    :param dt:          timestep of the data, in seconds. If given,
                        it is stored on the returned realization
                        (adjusted for any decimation).
    :type dt:           float, optional
    :param method:      system identification method. default is
                        "srim", other options are "okid-era",
                        "okid-era-dc", and "deterministic".
    :type method:       string, optional
    :param decimation:  decimation factor. default: no decimation
    :type decimation:   int, optional

    :return:            state-space realization ``(A,B,C,D)`` as a
                        :class:`mdof.Realization`, which unpacks like
                        the tuple but also carries the effective ``dt``
                        and a provenance record.
    :rtype:             :class:`mdof.Realization`
    """

    return system(inputs, outputs, dt=dt, **options)


def predict(realization, inputs, **options):
    """
    Prediction of ``output`` from system realization
    (``A``, ``B``, ``C``, ``D``) and ``inputs``.

    :param realization: realization in the form of state space
                        coefficients ``(A,B,C,D)``
    :type realization:  tuple of arrays
    :param inputs:      input time history on which to predict the
                        output response.
                        dimensions: :math:`(q,nt)`, where
                        :math:`q` = number of inputs, and
                        :math:`nt` = number of timesteps
    :type inputs:       array

    :return:            output response history.
                        dimensions: :math:`(p,nt)`, where
                        :math:`p` = number of outputs, and
                        :math:`nt` = number of timesteps
    :rtype:             array
    """

    from mdof.simulate import simulate
    out_pred = simulate(realization, inputs)
    return out_pred