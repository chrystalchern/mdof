import numpy as np
from scipy.fft import fft, fftfreq

def power_transfer(inputs, outputs, step): pass

def response_transfer(inputs, outputs, step): pass

def fourier_transfer(inputs, outputs, step):
    input_transform = fspec(inputs, step)
    output_transform = fspec(outputs, step)
    assert input_transform[0].all() == output_transform[0].all()
    return (1/input_transform[0], output_transform[1]/input_transform[1])

def pspec(series, step):
    import scipy.signal
def rspec(series, step):
    pass
def fspec(series, step):
    assert len(series.shape) == 1
    N = len(series)
    return np.array([fftfreq(N,step)[:N//2], 2.0/N*np.abs(fft(series)[:N//2])])

def _newmark():
    pass
