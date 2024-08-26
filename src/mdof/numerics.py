
def lsq_solver(options):
    import numpy.linalg
    return lambda *args: numpy.linalg.lstsq(*args, rcond=options.get('rcond',None))[0]

def lin_solver():
    pass

def svd_routine(lapack_driver="gesvd", **kwds):
    import scipy
    # default lapack_driver="gesdd"
    def _svd(*args):
        U,S,V = scipy.linalg.svd(*args, lapack_driver=lapack_driver, **kwds)
        return U,S,V.T.conj()
    return _svd

def form_observability():
    pass

def form_controllability():
    pass
    
def block_hankel(series):
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

def linear_interpolate(x, y, target_x):
    sorted_indices = sorted(range(len(x)), key = lambda i: x[i])
    x = x[sorted_indices]
    y = y[sorted_indices]
    i1 = max(np.where(x<=target_x)[0])
    i2 = min(np.where(x>=target_x)[0])
    x1 = x[i1]
    x2 = x[i2]
    y1 = y[i1]
    y2 = y[i2]
    target_y = y1 + (target_x-x1)*(y2-y1)/(x2-x1)    
    return target_y

cm2g = 0.0010197162129779