
def lsq_solver(options):
    import numpy.linalg
    return lambda *args: numpy.linalg.lstsq(*args, rcond=None)[0]

def lin_solver():
    pass

def svd_routine(**kwds):
    import scipy
    def _svd(*args):
        U,S,V = scipy.linalg.svd(*args, lapack_driver="gesvd")
        return U,S,V.T.conj()
    return _svd

def form_observability():
    pass

def decimate(series, decimation):
    import numpy as np
    if isinstance(series, np.ndarray):
        if len(series.shape) == 1:
            return series[np.arange(0,len(series),decimation)]
        else:
            return series[:,np.arange(0,series.shape[1],decimation)]
    if isinstance(series, list):
        return np.asarry(series)[np.arange(0,len(series),decimation)]