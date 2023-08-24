import numpy as np
from scipy.fft import fft, fftfreq
from .numerics import decimate

# def power_transfer(inputs, outputs, step): pass

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
    input_transform[0]=np.real(input_transform[0]) # prevents unwarranted "divide by zero" warning
    return (1/input_transform[0], output_transform[1]/input_transform[1])


def fourier_spectrum(series, step, period_band=None, **options):
    """
    Fourier amplitude spectrum from a signal.

    :param series:      time series.
    :type series:       1D array
    :param step:        timestep.
    :type step:         float
    :param period_band: minimum and maximum period of interest, in seconds.

    :return:            (frequencies, amplitudes)
    :rtype:             tuple of arrays.
    """
    assert len(series.shape) == 1
    N = len(series)
    frequencies = fftfreq(N,step)[1:N//2]
    amplitudes = 2.0/N*np.abs(fft(series)[1:N//2])
    if period_band is not None:
        frequency_band = (1/period_band[1], 1/period_band[0])
        frequency_indices = np.logical_and(frequencies>frequency_band[0], frequencies<frequency_band[1])
        frequencies = frequencies[frequency_indices]
        amplitudes = amplitudes[frequency_indices]
    return np.array([frequencies, amplitudes])

def _newmark():
    pass
