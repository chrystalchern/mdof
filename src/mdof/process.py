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


from scipy.signal import detrend
def baseline_correct(signal, mode='linear'):
    """
    Subtract a trend from a signal. 
    Modes include constant or linear.
    """
    return detrend(signal, type=mode)

def window(signal, time, interval, overlap):
    """
    Separate a `signal` and corresponding `time` array into 
    windowed intervals. Each interval has a length of `interval`
    seconds. The windows overlap by `overlap` seconds.
    Omits the first and last intervals.
    """
    if overlap >= interval:
        raise ValueError("Overlap must be smaller than interval.")
    
    signal = np.array(signal)
    time = np.array(time)
    if signal.ndim != 1:
        raise ValueError("Signal must be a list or 1D array.")
    if time.ndim != 1:
        raise ValueError("Time must be a list or 1D array.")
    if signal.shape != time.shape:
        raise ValueError(f"Array dimension mismatch. signal shape is {signal.shape} while time shape is {time.shape}.")
    
    npts = len(signal)
    dt = time[1]-time[0]
    interval_steps = int(interval//dt)
    overlap_steps = int(overlap//dt)
    step = interval_steps - overlap_steps
    indices = np.arange(0,npts,step)
    time_windows = [time[i:i+interval_steps] for i in indices][1:-1]
    signal_windows = [signal[i:i+interval_steps] for i in indices][1:-1]
    
    return time_windows, signal_windows

import matplotlib.pyplot as plt
import mdof.transform
def moving_window_transfer(time, signal, period_range=(0,1), method='fourier'):
    times, signals = window(signal, time, 10, 1)
    n_windows = len(times)
    n_transform_pts = len(times[0])//2-1
    time_grid = np.empty((n_windows,n_transform_pts))
    period_grid = np.empty((n_windows,n_transform_pts))
    amplitude_grid = np.empty((n_windows,n_transform_pts))
    for i,(t,w) in enumerate(zip(times, signals)):
        time_center = np.mean(t)
        time_grid[i] = time_center*np.ones(n_transform_pts)
        periods, amplitudes = mdof.transform.fourier_spectrum(w,t[1]-t[0])
        period_grid[i] = periods
        amplitude_grid[i] = amplitudes
    fig,ax = plt.subplots()
    ax.pcolormesh(time_grid,period_grid,amplitude_grid, shading='gouraud')
    ax.set_ylim(period_range)
    return fig

if __name__ == "__main__":
    time = np.arange(0,200,0.01)
    mysignal = np.sin((5/2*np.pi)*time) + np.random.random(len(time))
    moving_window_transfer(time, mysignal)