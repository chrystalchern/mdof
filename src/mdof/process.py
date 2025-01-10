import numpy as np

def decimate(series, decimation):
    if isinstance(series, np.ndarray):
        if len(series.shape) == 1:
            return series[np.arange(0,len(series),decimation)]
        else:
            return series[:,np.arange(0,series.shape[1],decimation)]
    if isinstance(series, list):
        return np.asarry(series)[np.arange(0,len(series),decimation)]
    

def power(series):
    """
    Compute the power of a signal as the sum of squared absolute
    values divided by duration.
    """
    series = np.asarray(series)
    if series.ndim != 1:
        raise ValueError("Series must be a 1D array.")
    N = len(series)
    return (1/N)*np.sum(np.abs(series)**2)


def snr(signal, noise):
    """
    Compute the signal to noise ratio given the pure signal
    and corresponding noise.
    """
    if signal.shape != noise.shape:
        raise ValueError("Signal and noise must have equal shape.")
    return np.array([power(s)/power(n) for s,n in zip(signal,noise)])