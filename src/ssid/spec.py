import numpy as np
from scipy.fft import fft, fftfreq
from .numerics import decimate

def power_transfer(inputs, outputs, step): pass

def response_transfer(inputs, outputs, step, pseudo=False, decimation=None, **kwds):

    if decimation is not None:
        inputs = decimate(inputs, decimation=decimation)
        outputs = decimate(outputs, decimation=decimation)
        step = step*decimation

    from sdof import spectrum
    Din,  _,  Ain = spectrum(inputs,  step, **kwds)
    Dout, _, Aout = spectrum(outputs, step, **kwds)
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

def fourier_transfer(inputs, outputs, step, decimation=None, **kwds):

    if decimation is not None:
        inputs = decimate(inputs, decimation=decimation)
        outputs = decimate(outputs, decimation=decimation)
        step = step*decimation

    assert len(inputs) == len(outputs)
    input_transform = fspec(inputs, step, **kwds)
    output_transform = fspec(outputs, step, **kwds)
    input_transform[0]=np.real(input_transform[0]) # prevents unwarranted "divide by zero" warning
    return (1/input_transform[0], output_transform[1]/input_transform[1])

def pspec(series, step):
    import scipy.signal

def rspec(series, step):
    pass

def fspec(series, step, period_band=None, **kwds):
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
