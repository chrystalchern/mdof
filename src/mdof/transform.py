import numpy as np
from scipy.fft import fft, fftfreq
from scipy import signal
from .numerics import decimate

def power_transfer(inputs, outputs, step, **options):
    """
    Power spectrum transfer function from input and output data.

    :param inputs:      input time history. dimensions: :math:`(q,nt)`, where
                        :math:`q` = number of inputs, and :math:`nt` = number of timesteps
    :type inputs:       array
    :param outputs:     output response history.
                        dimensions: :math:`(p,nt)`, where :math:`p` = number of outputs, and
                        :math:`nt` = number of timesteps
    :type outputs:      array
    :param step:        timestep.
    :type step:         float
    :param decimation:  decimation factor. default: 1
    :type decimation:   int, optional

    :return:            (periods, amplitudes)
    :rtype:             tuple of arrays
    """
    decimation = options.get("decimation", None)

    if decimation is not None:
        inputs = decimate(inputs, decimation=decimation)
        outputs = decimate(outputs, decimation=decimation)
        step = step*decimation

    assert len(inputs) == len(outputs)
    input_transform = power_spectrum(inputs, step, **options)
    output_transform = power_spectrum(outputs, step, **options)
    return (input_transform[0], output_transform[1]/input_transform[1])

def response_transfer(inputs, outputs, step, **options):
    """
    Response spectrum transfer function from input and output data.
    
    :param inputs:      input time history. dimensions: :math:`(q,nt)`, where
                        :math:`q` = number of inputs, and :math:`nt` = number of timesteps
    :type inputs:       array
    :param outputs:     output response history.
                        dimensions: :math:`(p,nt)`, where :math:`p` = number of outputs, and
                        :math:`nt` = number of timesteps
    :type outputs:      array
    :param step:        timestep.
    :type step:         float
    :param pseudo:      if True, uses pseudo accelerations. default: False
    :type pseudo:       bool, optional
    :param decimation:  decimation factor. default: 1
    :type decimation:   int, optional

    :return:            (periods, amplitudes)
    :rtype:             tuple of arrays
    """
    pseudo = options.get("pseduo", False)
    decimation = options.get("decimation", None)

    if decimation is not None:
        inputs = decimate(inputs, decimation=decimation)
        outputs = decimate(outputs, decimation=decimation)
        step = step*decimation

    from sdof import spectrum
    Din,  _,  Ain = spectrum(inputs,  step, **options)
    Dout, _, Aout = spectrum(outputs, step, **options)
    periods = Din[0]

    if pseudo:
        input_spectrum = Din[1,:]*(2*np.pi/periods)**2
    else:
        input_spectrum = Ain[1]

    if pseudo:
        output_spectrum = Dout[1,:]*(2*np.pi/periods)**2
    else:
        output_spectrum = Aout[1]

    return (periods, output_spectrum/input_spectrum)

def fourier_transfer(inputs, outputs, step, **options):
    """
    Fourier spectrum transfer function from input and output data.

    :param inputs:      input time history. dimensions: :math:`(q,nt)`, where
                        :math:`q` = number of inputs, and :math:`nt` = number of timesteps
    :type inputs:       array
    :param outputs:     output response history.
                        dimensions: :math:`(p,nt)`, where :math:`p` = number of outputs, and
                        :math:`nt` = number of timesteps
    :type outputs:      array
    :param step:        timestep.
    :type step:         float
    :param decimation:  decimation factor. default: 1
    :type decimation:   int, optional

    :return:            (periods, amplitudes)
    :rtype:             tuple of arrays
    """
    decimation = options.get("decimation", None)

    if decimation is not None:
        inputs = decimate(inputs, decimation=decimation)
        outputs = decimate(outputs, decimation=decimation)
        step = step*decimation

    assert len(inputs) == len(outputs)
    input_transform = fourier_spectrum(inputs, step, **options)
    output_transform = fourier_spectrum(outputs, step, **options)
    return (input_transform[0], output_transform[1]/input_transform[1])

def power_spectrum(series, step, period_band=None, **options):
    """
    Power spectral density of a signal, as a function of period.

    :param series:      time series.
    :type series:       1D array
    :param step:        timestep.
    :type step:         float
    :param period_band: minimum and maximum period of interest, in seconds.
    :type period_band:  tuple, optional

    :return:            (periods, amplitudes)
    :rtype:             tuple of arrays.
    """
    periods, amplitudes = fourier_spectrum(series, step, **options)
    if period_band is not None:
        period_indices = np.logical_and(periods>period_band[0], periods<period_band[1])
        periods = periods[period_indices]
        amplitudes = amplitudes[period_indices]
    return np.array([periods, (np.abs(amplitudes))**2])

# def _power_spectrum(series, step, **options):  # equivalent to power_spectrum()
#     frequencies, power_spectral_density = signal.periodogram(series, step)
#     return (1/frequencies, power_spectral_density)

def fourier_spectrum(series, step, period_band=None, **options):
    """
    Fourier amplitude spectrum of a signal, as a function of period.

    :param series:      time series.
    :type series:       1D array
    :param step:        timestep.
    :type step:         float
    :param period_band: minimum and maximum period of interest, in seconds.
    :type period_band:  tuple, optional

    :return:            (periods, amplitudes)
    :rtype:             tuple of arrays.
    """
    assert len(series.shape) == 1
    N = len(series)
    frequencies = fftfreq(N,step)[1:N//2+2]
    amplitudes = 2.0/N*np.abs(fft(series)[1:N//2+2])
    periods = 1/frequencies
    if period_band is not None:
        period_indices = np.logical_and(periods>period_band[0], periods<period_band[1])
        periods = periods[period_indices]
        amplitudes = amplitudes[period_indices]
    return np.array([periods, amplitudes])

def _newmark():
    pass
 