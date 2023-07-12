import numpy as np
from scipy.fft import fft, fftfreq

def power_transfer(inputs, outputs, step): pass

def response_transfer(inputs, outputs, step, **kwds):
    from sdof import spectrum
    input_spectrum, *_  = spectrum(inputs,  step, **kwds)
    output_spectrum, *_ = spectrum(outputs, step, **kwds)
    return (input_spectrum[0], output_spectrum[1]/input_spectrum[1])

def fourier_transfer(inputs, outputs, step, **kwds):
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
