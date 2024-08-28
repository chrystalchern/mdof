import numpy as np

def decimate(series, decimation):
    if isinstance(series, np.ndarray):
        if len(series.shape) == 1:
            return series[np.arange(0,len(series),decimation)]
        else:
            return series[:,np.arange(0,series.shape[1],decimation)]
    if isinstance(series, list):
        return np.asarry(series)[np.arange(0,len(series),decimation)]
    