import numpy as np
from scipy.fft import fft, fftfreq
from scipy import signal

## Transfer functions    
def power_transfer(inputs, outputs, step, **options):
    """
    Power spectrum transfer function from input and output data.

    :param inputs:      input time history.
    :type inputs:       1D array
    :param outputs:     output response history.
    :type outputs:      1D array
    :param step:        timestep.
    :type step:         float
    :param period_band: (optional) minimum and maximum period of interest, in seconds.
    :type period_band:  tuple

    :return:            (periods, amplitudes)
    :rtype:             tuple of arrays
    """

    assert len(inputs) == len(outputs)
    input_transform = power_spectrum(inputs, step, **options)
    output_transform = power_spectrum(outputs, step, **options)
    return (input_transform[0], output_transform[1]/input_transform[1])


def fourier_transfer(inputs, outputs, step, **options):
    """
    Fourier spectrum transfer function from input and output data.

    :param inputs:      input time history.
    :type inputs:       1D array
    :param outputs:     output response history.
    :type outputs:      1D array
    :param step:        timestep.
    :type step:         float
    :param period_band: (optional) minimum and maximum period of interest, in seconds.
    :type period_band:  tuple

    :return:            (periods, amplitudes)
    :rtype:             tuple of arrays
    """

    assert len(inputs) == len(outputs)
    input_transform = fourier_spectrum(inputs, step, **options)
    output_transform = fourier_spectrum(outputs, step, **options)
    return (input_transform[0], output_transform[1]/input_transform[1])


def response_transfer(inputs, outputs, step, **options):
    """
    Response spectrum transfer function from input and output data.
    
    :param inputs:      input time history.
    :type inputs:       1D array
    :param outputs:     output response history.
    :type outputs:      1D array
    :param step:        timestep.
    :type step:         float
    :param pseudo:      (optional) if True, uses pseudo accelerations. default: False
    :type pseudo:       bool
    :param period_band: (optional) minimum and maximum period of interest, in seconds.
    :type period_band:  tuple

    :return:            (periods, amplitudes)
    :rtype:             tuple of arrays
    """
    pseudo = options.get("pseudo", False)
    if 'period_band' in options.keys():
        pmin, pmax = options['period_band']
        options['periods'] = np.linspace(pmin, pmax, 200)

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
    

### Spectra
def fdd_spectrum(outputs, step, **options):
    """
    Frequency Domain Decomposition spectrum [1]_ from output data.

    :param outputs:     output response history.
                        dimensions: :math:`(p,nt)`, where :math:`p` = number of outputs, and
                        :math:`nt` = number of timesteps
    :type outputs:      ND array.
    :param step:        timestep.
    :type step:         float
    :param period_band: (optional) minimum and maximum period of interest, in seconds.
    :type period_band:  tuple

    :return:            (periods, amplitudes)
    :rtype:             tuple of arrays
    
    References
    ----------
    .. [1]  Brincker, R., Zhang, L., & Andersen, P. (2001). Modal identification of
            output-only systems using frequency domain decomposition. Smart materials
            and structures, 10(3), 441. (https://doi.org/10.1088/0964-1726/10/3/303
    """
    frequencies, _, S = fdd(outputs, step)
    periods = 1/frequencies
    amplitudes = S[0,:]
    if 'period_band' in options.keys():
        pmin, pmax = options['period_band']
        period_indices = np.logical_and(periods>pmin, periods<pmax)
        periods = periods[period_indices]
        amplitudes = amplitudes[period_indices]
    return (periods, amplitudes)


def fdd(outputs, step):
    """
    Frequency Domain Decomposition [2]_ from output data.

    :param outputs:     output response history.
                        dimensions: :math:`(p,nt)`, where :math:`p` = number of outputs, and
                        :math:`nt` = number of timesteps
    :type outputs:      ND array.
    :param step:        timestep.
    :type step:         float

    :return:            (frequencies, ```U```, ```S```)
    :rtype:             tuple of arrays
    
    References
    ----------
    .. [2]  Brincker, R., Zhang, L., & Andersen, P. (2001). Modal identification of
            output-only systems using frequency domain decomposition. Smart materials
            and structures, 10(3), 441. (https://doi.org/10.1088/0964-1726/10/3/303
    """
    
    if len(outputs.shape) == 1:
        outputs = outputs[None,:]
    p,nt = outputs.shape
    transform_length = (nt//2)-1

    Gyy = np.zeros((p,p,transform_length))
    U = np.zeros((p,p,transform_length))
    S = np.zeros((p,transform_length))

    frequencies = fftfreq(nt,step)[1:nt//2]     
            
    for i,seriesi in enumerate(outputs):
        for j,seriesj in enumerate(outputs[i:]):
            amplitudesi = 2.0/nt*np.abs(fft(seriesi)[1:nt//2])
            amplitudesj = 2.0/nt*np.abs(fft(seriesj)[1:nt//2])
            Gyy[i,j] = Gyy[j,i] = amplitudesi * amplitudesj

    for i in range(transform_length):
        U[:,:,i],S[:,i],_ = np.linalg.svd(Gyy[:,:,i])

    return (frequencies, U, S)


# def _power_spectrum(series, step, **options):  # equivalent to power_spectrum()
#     frequencies, power_spectral_density = signal.periodogram(series, step)
#     return (1/frequencies, power_spectral_density)
def power_spectrum(series, step, period_band=None, **options):
    """
    Power spectrum of a signal, as a function of period (i.e., periodogram).

    :param series:      time series.
    :type series:       1D array
    :param step:        timestep.
    :type step:         float
    :param period_band: (optional) minimum and maximum period of interest, in seconds.
    :type period_band:  tuple

    :return:            (periods, amplitudes)
    :rtype:             tuple of arrays.
    """
    periods, amplitudes = fourier_spectrum(series, step, **options)
    if period_band is not None:
        period_indices = np.logical_and(periods>period_band[0], periods<period_band[1])
        periods = periods[period_indices]
        amplitudes = amplitudes[period_indices]
    return np.array([periods, (np.abs(amplitudes))**2])


def fourier_spectrum(series, step, period_band=None, **options):
    """
    Fourier amplitude spectrum of a signal, as a function of period.

    :param series:      time series.
    :type series:       1D array
    :param step:        timestep.
    :type step:         float
    :param period_band: (optional) minimum and maximum period of interest, in seconds.
    :type period_band:  tuple

    :return:            (periods, amplitudes)
    :rtype:             tuple of arrays.
    """
    assert len(series.shape) == 1
    N = len(series)
    frequencies = fftfreq(N,step)[1:N//2]
    amplitudes = 2.0/N*np.abs(fft(series)[1:N//2])
    periods = 1/frequencies
    if period_band is not None:
        period_indices = np.logical_and(periods>period_band[0], periods<period_band[1])
        periods = periods[period_indices]
        amplitudes = amplitudes[period_indices]
#     else:
#         period_band_warning = '''
# Warning: Recommend specifying a period band 
# to omit frequencies lower than Nyquist.
# For example, a common sampling rate is 0.01 seconds;
# periods should then be longer than 0.02 seconds by a good margin.
# '''
#         import warnings
#         warnings.warn(period_band_warning, category='UserWarning')

    return np.array([periods, amplitudes])

 