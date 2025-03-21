from sdof import _fsdof_integrate, _sdof_config
from ctypes import c_double, POINTER

def integrate(m,c,k,f,dt, u0=0.0, v0=0.0,
              out  =  None,
              alpha_m: float = 1.0,
              alpha_f: float = 1.0,
              beta   : float = 0.25,
              gamma  : float = 0.50
    ):

    import numpy as np
    if out is None:
        output = np.empty((3,len(f)))
    else:
        output = out
    output[:2,0] = u0, v0

    config = _sdof_config(
                alpha_m = alpha_m,
                alpha_f = alpha_f,
                beta    = beta,
                gamma   = gamma
    )

    _fsdof_integrate(config, m, c, k, 1.0, len(f), np.asarray(f).ctypes.data_as(POINTER(c_double)), dt, output)
    return output