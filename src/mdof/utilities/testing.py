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

def align_signals(signal1, signal2, times=None, verbose=False):

    assert signal1.shape == signal2.shape, "Please give two signals of equal shape."
    signal1 = np.atleast_2d(signal1)
    signal2 = np.atleast_2d(signal2)
    npts = signal1.shape[1]

    max_lag = 0
    s1_aligned = []
    s2_aligned = []

    for s1,s2 in zip(signal1,signal2):
        # Compute the cross-correlation between signal1 and signal2
        correlation = correlate(s1, s2, mode='full')
        
        # Get the lags corresponding to the cross-correlation
        lags = correlation_lags(len(s1), len(s2), mode='full')

        if verbose:
            _,ax=plt.subplots()
            ax.scatter(lags,correlation)
            ax.set_xlabel("lag")
            ax.set_ylabel("correlation")
            plt.show()
        
        # Find the lag corresponding to the maximum correlation
        lag = lags[np.argmax(correlation)]
        
        # Trim the signals
        if lag > 0:
            s1_aligned.append(s1[lag:])
            s2_aligned.append(s2[:len(s1[lag:])])
        elif lag < 0:
            s2_aligned.append(s2[-lag:])
            s1_aligned.append(s1[:len(s2[-lag:])])
        else:
            s1_aligned.append(s1)
            s2_aligned.append(s2)

        max_lag = max(abs(lag), max_lag)

    new_npts = npts-max_lag
    signal1_aligned = np.array([s[:new_npts] for s in s1_aligned])
    signal2_aligned = np.array([s[:new_npts] for s in s2_aligned])

    if times is None:
        times = np.arange(0,new_npts,1)
    else:
        times = times[:new_npts]
    
    return max_lag, signal1_aligned, signal2_aligned, times