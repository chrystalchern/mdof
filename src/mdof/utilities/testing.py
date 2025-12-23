import mdof
from mdof import modal
from time import time
from control import ss, forced_response
import numpy as np

class Timer:
    def __init__(self):
        # Timer starts
        self.start = self.last = time()
        self.laps = []  

    def lap(self, message=None):
        if message is None:
            message = ""
        # The current lap-time
        laptime = round((time() - self.last), 2)
        self.laps.append(laptime)

        # Printing the lap number,
        # lap-time and total time
        print(f"{self.laps[-1]:<10} {message}")


        # Updating the previous total time
        # and lap number
        self.last = time()


def test_method(method, inputs, outputs, dt, t, **conf):
    time0 = time()
    A,B,C,D = mdof.system(method=method, inputs=inputs, outputs=outputs, **conf)
    time1 = time()
    modedict = modal.system_modes((A,B,C,D),dt,**conf)
    model = {
                "time":    time1-time0,
                "ypred":   forced_response(ss(A,B,C,D,dt), T=t, U=inputs, squeeze=False, return_x=False).outputs,
                "modes":   modedict,
                "period":  [1/value['freq'] for value in modedict.values()],
                "damping": [value['damp'] for value in modedict.values()]
            }
    return model

def mode_statistics(mode_results, key):
    import numpy as np
    values = [result[key] for results in mode_results for result in results]
    mean = np.mean(values)
    std = np.std(values)
    closest_values = [results[np.argmin(np.abs(mean-[result[key] for result in results]))] for results in mode_results]
    return [
        dict(**item, distance=(item[key]-mean)/std) for item in closest_values
    ]

def mode_set(mode_results, key):
    return 


def husid(accel_series, intensity_measure='arias'):
    if intensity_measure == 'arias':
        # Compute Arias Intensity (cumulative of squared acceleration)
        # intensity = np.tril(np.ones(len(accel_series)))@accel_series**2
        intensity = np.cumsum(accel_series**2)
    elif intensity_measure == 'isaacson':
        # Compute Isaacson's Intensity (cumulative of squared acceleration but might use different methodology)
        intensity = np.cumsum(accel_series**2)  # Example, may vary depending on specific implementation
    elif intensity_measure == 'cav':
        # Compute Cumulative Absolute Velocity (CAV) (cumulative of absolute acceleration)
        intensity = np.cumsum(np.abs(accel_series))
    elif intensity_measure == 'pga':
        # Compute Peak Ground Acceleration (cumulative max absolute value)
        intensity = np.array([np.max(np.abs(accel_series[:i]))
                              for i in range(len(accel_series))])
    elif intensity_measure == 'pgv':
        # Compute Peak Ground Velocity (cumulative max absolute value)
        veloc_series = np.cumsum(accel_series)  # Integrating to get velocity
        intensity = np.array([np.max(np.abs(veloc_series[:i]))
                              for i in range(len(veloc_series))])
    else:
        raise ValueError(f"Unknown intensity measure: {intensity_measure}")
    # Normalize intensity by the maximum value
    cumulative_husid = intensity / intensity[-1]
    return cumulative_husid

def intensity_bounds(accel_series, lb=0.005, ub=0.995, intensity_measure='arias'):
    cumulative_husid = husid(accel_series, intensity_measure)
    ilb = next(x for x, val in enumerate(cumulative_husid) if val > lb)
    iub = next(x for x, val in enumerate(cumulative_husid) if val > ub)
    return (ilb, iub)

def truncate_by_bounds(series, bounds):
    ilb, iub = bounds
    if series.ndim == 1:
        return series[ilb:iub]
    elif series.ndim == 2:
        return series[:,ilb:iub]

from scipy.signal import correlate, correlation_lags
import matplotlib.pyplot as plt

def align_signals(signal1, signal2, times=None, verbose=False, max_lag_allowed=None):
    """
    Align two 1D signals to maximize their cross-correlation.
    If max_lag_allowed is given, then the signals will only
    be aligned if it does not cause a lag > max_lag_allowed
    
    :param signal1:         first signal, 1D np.ndarray
    :param signal2:         second signal, 1D np.ndarray
    :param times:           optional, time array. must be
                            provided if using max_lag_allowed
    :param verbose:         if True or >=1, prints feedback
                            regarding alignment. If >=2, plots
                            lags vs correlations. Default False.
    :param max_lag_allowed: optional: maximum time, in seconds,
                            that the signals may be lagged for 
                            alignment
    """

    assert len(signal1) == len(signal2), "Please give two signals of equal length."
    assert signal1.ndim == 1 and signal2.ndim == 1, "Both signals must be 1D np.ndarray."
    npts = len(signal1)

    # Compute the cross-correlation between signal1 and signal2
    correlation = correlate(signal1, signal2, mode='full')
    
    # Get the lags corresponding to the cross-correlation
    lags = correlation_lags(len(signal1), len(signal2), mode='full')

    if verbose >= 2:
        _,ax=plt.subplots()
        ax.scatter(lags,correlation)
        ax.set_xlabel("lag")
        ax.set_ylabel("correlation")
        plt.show()
    
    # Sort lags from max to min correlations
    corr_idxs_sorted = np.argsort(-correlation)
    lags_sorted = lags[corr_idxs_sorted]

    # Find the lag corresponding to the maximum correlation
    # Only allow lags shorted than the max allowed, if specified
    if max_lag_allowed is not None:
        if times is None:
            raise ValueError("times must be given if max_lag_allowed applied.")
        dt = times[1]-times[0]
        lag_times = abs(lags_sorted)*dt
        lags_allowed_indices = np.where(lag_times < max_lag_allowed)[0]
        if dt > max_lag_allowed:
            if verbose:
                print(f"Time step ({dt} s) is larger than max allowed lag ({max_lag_allowed} s); no alignment applied.")
                lag = 0
        elif len(lags_allowed_indices) == 0: # TODO: confirm that this case is different than the previous.
            lag = 0
        else:
            lag = lags_sorted[lags_allowed_indices[0]]
    else:
        lag = lags_sorted[0]

    # Trim the signals if there is a lag
    if lag > 0:
        signal1_aligned = signal1[lag:]
        signal2_aligned = signal2[:len(signal2[lag:])]
    elif lag < 0:
        signal2_aligned = signal2[-lag:]
        signal1_aligned = signal1[:len(signal2[-lag:])]
    else:
        signal1_aligned = signal1
        signal2_aligned = signal2

    if verbose and time is not None:
        lag_time = lag*dt
        print(f"Signals aligned with lag time {lag_time} s.")

    # TODO: confirm that lag always trims signals to equal lengths
    assert len(signal1_aligned) == len(signal2_aligned), "trimmed signals aren't equal lengths"

    if times is None:
        times = np.arange(len(signal1_aligned))
    else:
        times = times[:len(signal1_aligned)]
    
    return lag, signal1_aligned, signal2_aligned, times