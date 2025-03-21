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
    windowed intervals. Omits the first and last intervals.
    Each interval has a length of `interval` seconds.
    The windows overlap by `overlap` seconds.
    The `signal` can be a 1D array (shape (`nt`,)) or 
    2D array (shape (`n`,`nt`)), where `nt` is the number of
    time samples and `n` is the number of signal channels.

    Returns 1D array `time_windows` (`nt`,)
    and 2D array `signal_windows` (`n`,`nt`)
    """
    if overlap >= interval:
        raise ValueError("Overlap must be smaller than interval.")
    
    signal = np.atleast_2d(signal)
    time = np.atleast_1d(time)
    if signal.ndim > 2:
        raise ValueError("Signal must be a list, 1D, or 2D array.")
    if time.ndim != 1:
        raise ValueError("Time must be a list or 1D array.")
    if signal.shape[-1] != time.shape[-1]:
        raise ValueError(f"Array dimension mismatch. signal shape is {signal.shape} while time shape is {time.shape}.")
    
    npts = signal.shape[1]
    dt = time[1]-time[0]
    interval_steps = int(interval//dt)
    overlap_steps = int(overlap//dt)
    step = interval_steps - overlap_steps
    indices = np.arange(0,step*((npts-overlap_steps)//step),step)
    time_windows = [time[i:i+interval_steps] for i in indices]
    signal_windows = [signal[:,i:i+interval_steps] for i in indices]
    
    return time_windows, signal_windows

import matplotlib.pyplot as plt
import plotly.graph_objects as go
def plot_moving_window(time_grid, period_grid, amplitude_grid, period_range=(0,1), plotly=False, **options):
    """
    Display a signal's moving window frequency response.
    For use with `moving_window_spectrum`, `moving_window_transfer`,
    or similar function that provides grid-shaped frequency
    response spectra at each time window of a time series.
    """
    x = time_grid.T[0]
    y = period_grid[0]
    z = amplitude_grid.T
    if plotly:
        fig = go.Figure(data=go.Heatmap(
            z = z,
            x = x,
            y = y,
            zsmooth = 'best',
            colorscale = 'viridis',
            coloraxis = 'coloraxis'
        ))
        fig.update_layout(
            xaxis_title = "time (s)",
            yaxis_title = "Period (s)",
            yaxis_range=period_range)
        fig.update_coloraxes(
            showscale = False)
    else:
        if 'figax' in options.keys():
            fig,ax = options['figax']
        else:
            fig,ax = plt.subplots()
        ax.pcolormesh(time_grid,period_grid,amplitude_grid, shading='gouraud')
        ax.set_xlabel("Time (s)")
        ax.set_ylim(period_range)
        ax.set_ylabel("Period (s)")
    return fig

import mdof.transform
def moving_window_spectrum(time, signal, interval=10, overlap=1, method='fourier', **options):
    """
    Obtain grid data for plotting a colormesh or heatmap of the
    moving window frequency response spectrum of a signal.
    """
    time_windows, signal_windows = window(signal, time, interval, overlap)
    # `window` gives `signal_windows` as a list of 2D arrays. To use with spectrum methods,
    # `signal` must be single output and passed into `window` as a 1D array or list.
    if signal_windows[0].shape[0] != 1:
        raise ValueError("signal must be a list or 1D array")
    
    # choose the spectrum method
    if method == 'fourier':
        spectrum = mdof.transform.fourier_spectrum
        n_transform_pts = len(time_windows[0])//2-1
    elif method == 'power':
        spectrum = mdof.transform.power_spectrum
        n_transform_pts = len(time_windows[0])//2-1
    elif method == 'response':
        spectrum = mdof.transform.response_spectrum
        n_transform_pts = options.get('n_transform_pts',200)

    n_windows = len(time_windows)
    dt = time_windows[0][1] - time_windows[0][0]

    time_grid = np.empty((n_windows,n_transform_pts))
    period_grid = np.empty((n_windows,n_transform_pts))
    amplitude_grid = np.empty((n_windows,n_transform_pts))

    normalize = options.get('normalize', False)
    for i,(t,w) in enumerate(zip(time_windows, signal_windows)):
        time_center = np.mean(t)
        time_grid[i] = time_center*np.ones(n_transform_pts)
        periods, amplitudes = spectrum(w[0],dt,**options)
        period_grid[i] = periods
        if normalize:
            amplitudes = amplitudes/np.max(amplitudes)
        amplitude_grid[i] = amplitudes
    return time_grid, period_grid, amplitude_grid

def signal_spectrum(time, signal, period_range=(0,1), method='fourier', **options):
    """
    Display the signal's time series response
    together with its moving window frequency response. 
    """
    time_grid, period_grid, amplitude_grid = moving_window_spectrum(time,
                                                                    signal,
                                                                    method=method,
                                                                    **options)
    fig,ax = plt.subplots(2,1, figsize=(12,5), sharex=True)
    ax[0].plot(time,signal)
    ax[0].set_ylabel("Signal ($m/s^2$)")
    ax[0].set_title("Time Series Response")
    plot_moving_window(time_grid,
                        period_grid,
                        amplitude_grid,
                        period_range,
                        figax = (fig,ax[1]),
                        **options)
    ax[1].set_title("Moving Window Frequency Response")
    return fig

def moving_window_transfer(time, inputs, outputs, interval=10, overlap=1, period_range=(0,1), method='fourier', **options):
    """
    Obtain grid data for plotting a colormesh or heatmap of the
    moving window frequency response spectrum of a transfer function.
    """
    time_windows, input_windows = window(inputs, time, interval, overlap)
    _, output_windows = window(outputs, time, interval, overlap)

    # choose the transfer function method
    if method=='fourier':
        transfer = mdof.transform.fourier_transfer
        n_transform_pts = len(time_windows[0])//2-1
        # fourier transfer works with SISO only
        if (input_windows[0].shape[0] != 1) or (output_windows[0].shape[0] != 1):
            raise ValueError("fourier transfer method works with single input, single output only")
    elif method=='response':
        transfer = mdof.transform.response_transfer
        n_transform_pts = options.get('n_transform_pts',200)
        # response transfer works with SISO only
        if (input_windows[0].shape[0] != 1) or (output_windows[0].shape[0] != 1):
            raise ValueError("response transfer method works with single input, single output only")
    elif method=='fdd':
        transfer = mdof.transform.fdd_spectrum
        n_transform_pts = len(time_windows[0])//2-1

    n_windows = len(time_windows)
    dt = time_windows[0][1] - time_windows[0][0]

    time_grid = np.empty((n_windows,n_transform_pts))
    period_grid = np.empty((n_windows,n_transform_pts))
    amplitude_grid = np.empty((n_windows,n_transform_pts))

    normalize = options.get('normalize', False)
    for i,(t,wi,wo) in enumerate(zip(time_windows,input_windows,output_windows)):
        time_center = np.mean(t)
        time_grid[i] = time_center*np.ones(n_transform_pts)
        if method in ['fourier', 'response']:
            periods, amplitudes = transfer(wi[0],wo[0],dt,**options)
        elif method == 'fdd':
            periods, amplitudes = transfer(wo,dt,**options)
        if method != 'response':
            # periods are in reverse order for fourier and fdd
            periods = periods[::-1]
            amplitudes = amplitudes[::-1]
        period_grid[i] = periods
        if normalize:
            amplitudes = amplitudes/np.max(amplitudes)
        amplitude_grid[i] = amplitudes
    return time_grid, period_grid, amplitude_grid


if __name__ == "__main__":
    tf = 100
    dt = 0.01
    npts = tf/dt
    time = np.arange(0,tf,dt)
    Tn1 = 1
    Tn2 = 2
    mysignal = np.concatenate((np.sin((2*np.pi/Tn1)*time[:int(npts/2)]),np.sin((2*np.pi/Tn2)*time[int(npts/2):]))) + np.random.random(len(time))
    fig = signal_spectrum(time, mysignal, period_range=(0,Tn2+1))
    plt.show()
