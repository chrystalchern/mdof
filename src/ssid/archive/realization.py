
class TimeSeries:
    pass

class StateSpaceModel:
    pass

class Realization:
    def system(self, *args, **kwds):
        a,b,c,d = ssid.system(*args, **kwds)
        self.a = a
        self.b = b
        self.c = c
        return a,b,c,d
    
    def system_modes(self, realization, dt, **kwds):
        if self.decimation is not None:
            dt *= decimation
        return ssid.system_modes(realization, dt, **kwds)

realization = Realization()

