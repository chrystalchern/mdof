from typing import Any


class Config(dict):
    def __setattr__(self, name, value):
        self[name]  = value
        
    def __getattribute__(self, __name: str) -> Any:
        if __name in self:
            return self[__name]
        return super().__getattribute__(__name)


def _decimate(series, d):
    import numpy as np
    return series[np.arange(0,len(series),d)]


def extract_channels(event, channels, decimate=1, response='accel'):
    import numpy as np
    import sys

    def find(chan):
        return event.match("l", station_channel=f"{chan}")

    data = np.asarray([
        _decimate(getattr(find(chan),response).data, d=decimate)
        for chan in channels if find(chan) is not None
    ])

    if len(data) == 0:
        raise ValueError("No channels found")
    
    elif len(data) != len(channels):
        print(f"Only extracted {len(data)} channels, {len(channels)-len(data)} missing.", file=sys.stderr)

    dt = find(channels[0]).accel["time_step"]*decimate

    return data, dt


def list_files(directory, pattern):
    from pathlib import Path
    return Path(directory).glob(pattern)


def create_time_vector(nt,dt):
    '''
    nt: number of timesteps
    dt: timestep
    tf: final time
    '''
    import numpy as np
    tf = nt*dt                             # final time
    times = np.linspace(0, tf, nt)
    return nt, dt, tf, times


def gather_system_parameters(T=[3,2],m=10,c=0.1,Phi=None):
    '''
    T: natural periods
    m: mass coefficient
    c: damping coefficient
    k: stiffness coefficient
    Phi: mode shapes
    '''

    import numpy as np
    T = np.array(T)
    k = (2*np.pi/T)**2*m
    if Phi is None:
        Phi = np.eye(len(T))

    print(f"Periods: {T}")
    print(f"Modeshapes: \n{Phi}")
    return T, m, c, k, Phi


def generate_input(times, series=None, seed=None):
    import numpy as np
    if seed is None:
        seed = 23
    rng = np.random.default_rng(seed=seed)
    if series is None:
        input_motion = sum(np.sin(times*i/10/np.pi)*rng.normal(0,0.1) for i in range(200))
    else:
        input_motion = series
    return input_motion


import json
import numpy as np
class json_serialize(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, complex):
            return str(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)