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
        intensity = np.tril(np.ones(len(accel_series)))@accel_series**2
    cumulative_husid = intensity/intensity[-1]
    return cumulative_husid

def intensity_bounds(accel_series, lb=0.005, ub=0.995, intensity_measure='arias'):
    cumulative_husid = husid(accel_series, intensity_measure)
    ilb = next(x for x, val in enumerate(cumulative_husid) if val > lb)
    iub = next(x for x, val in enumerate(cumulative_husid) if val > ub)
    return (ilb, iub)

def truncate_by_bounds(series, bounds):
    ilb, iub = bounds
    return series[ilb:iub]
