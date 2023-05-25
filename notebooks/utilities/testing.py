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