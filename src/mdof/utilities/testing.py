import mdof
from mdof import modal
from time import time
from control import ss, forced_response

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