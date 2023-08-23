# default imports

from .system import system

from control.matlab import impulse as _impulse


def impulse(system, t, **kwds):
    from control import ss
    dt = t[1] - t[0]
    a,t = _impulse(ss(*system, dt), t, **kwds)
    return a.squeeze()*dt,t
