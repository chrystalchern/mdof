class Config(dict):
    def __setattr__(self, name, value):
        self[name]  = value


def _decimate(series, d):
    import numpy as np

    return series[np.arange(0,len(series),d)]

def extract_channels(event, channels, decimate=1):
    import numpy as np

    data = np.array([
        _decimate(event.match("r", file_name=f".*{chan}\.[vV]2").accel.data, d=decimate)
        for chan in channels
    ]).T

    dt = event.match("r", file_name=rf".*{channels[0]}\.[vV]2").accel["time_step"]*decimate 
    return data, dt


def list_files(directory, pattern):
    from pathlib import Path
    return Path(directory).glob(pattern)