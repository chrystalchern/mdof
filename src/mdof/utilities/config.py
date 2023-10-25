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

def extract_channels(event, channels, decimate=1):
    import numpy as np
    import sys

    def find(chan):
        #event.match("r", file_name=f".*{chan}\.[vV]2")
        return event.match("l", station_channel=f"{chan}")

    data = np.asarray([
        _decimate(find(chan).accel.data, d=decimate)
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