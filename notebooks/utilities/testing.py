import ssid
from ssid import modes
from time import time
from control.matlab import ss, lsim

def test_method(method, input, output, dt, t, **conf):
    time0 = time()
    A,B,C,D = ssid.system(method=method, input=input, output=output, **conf)
    modedict = modes.modes((A,B,C,D),dt)
    model = {
                "time":    time()-time0,
                "ypred":   lsim(ss(A,B,C,D,dt),input,t)[0],
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