import ast
import sys
import json

import ssid
import ssid.modal
import quakeio
import numpy as np
from .okid import parse_okid


HELP = """
ssid [-p|w <>...] <method> <event> <inputs> <outputs>
ssid [-p|w <>...] <method> -i <inputs> -o <outputs>

-i/--inputs  FILE...   Input data
-o/--outputs FILE...   Output data

Methods:
  <method> <outputs>  <plots>
    srim    {dftmcs}   {m}
    okid    {dftmcs}   {m}
    spec    {ft}       {a}
    four    {ft}       {a}
    test

Outputs options
-p/--plot
    a/--accel-spect
-w/--write
    ABCD   system matrices
    d      damping
    freq   frequency
    cycl   cyclic frequency
    t      period
    m      modes
    c      condition-number
"""

class JSON_Encoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        return json.JSONEncoder.default(self, obj)


def _decimate(series, d):
    import numpy as np
    return series[np.arange(0,len(series),d)]

def extract_channels(event, channels, decimate=1, permissive=True):
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
        raise ValueError(f"No channels found, requested '{channels}'")

    elif len(data) != len(channels):
        if permissive:
            print(f"Only extracted {len(data)} channels, {len(channels)-len(data)} missing.", file=sys.stderr)
        else:
            raise ValueError("Could not extract all desired channels")

    dt = find(channels[0]).accel["time_step"]*decimate 
    return data, dt


def parse_srim(argi, config):
    help="""
    SRIM -- System Identification with Information Matrix

    Parameters
    no          order of the observer Kalman ARX filter (formerly p).
    r           size of the state-space model used for 
                representing the system. (formerly orm/n)
    """
    config.update({"p"  :  5, "orm":  4})

    #argi = iter(args)
    channels = [[], []]
    for arg in argi:
        if arg == "--arx-order":
            config["no"] = int(next(argi))

        elif arg == "--dt":
            config["dt"] = float(next(argi))

        elif arg == "--ss-size":
            config["r"] = int(next(argi))

        elif arg in ["--help", "-h"]:
            print(help)
            sys.exit()

        elif arg == "--inputs":
            # inputs = next(argi)[1:-1].split(",")
            inputs = ast.literal_eval(next(argi))
            if isinstance(inputs, str):
                channels[0] = [int(inputs)]
            else:
                channels[0] = list(map(int, inputs))

        elif arg == "--outputs":
            # outputs = next(argi)[1:-1].split(",")
            outputs = ast.literal_eval(next(argi))
            if isinstance(outputs, str):
                channels[1] = [int(outputs)]
            else:
                channels[1] = list(map(int, outputs))

        elif arg == "--":
            continue

        else:
            config["event_file"] = arg

    event = quakeio.read(config["event_file"], exclusions=["filter*"])
    try:
        inputs,  dt = extract_channels(event, channels[0], decimate=8)
        outputs, dt = extract_channels(event, channels[1], decimate=8)
    except Exception as e:
        print(json.dumps({"error": str(e), "data": []}))
        return

    config["dt"] = dt

    try:
        A,B,C,D = ssid.system(inputs=inputs, outputs=outputs, full=True, **config)
        ss_modes = ssid.modal.system_modes((A,B,C,D),dt)

    except Exception as e:
        print(json.dumps({"error": str(e), "data": []}))
        return

    output = [
        {
            "period":  1/mode["freq"],
            "frequency": mode["freq"],
            "damping":   mode["damp"],
            "emac":      mode["energy_condensed_emaco"],
            "mpc":       mode["mpc"],
        }
        for mode in sorted(ss_modes.values(), key=lambda item: item["freq"])
    ]

    print(json.dumps({"data": output}, cls=JSON_Encoder, indent=4))
    return config


def parse_args(args):
    outputs = []
    sub_parsers = {
        "srim": parse_srim,
        "test": parse_srim,
        "okid": parse_okid,
        "okid": parse_okid
    }
    config = {
            "method": None,
            "operation": None,
            "protocol": None
    }
    config["channels"] = channels = [[], []]

    argi = iter(args[1:])
    for arg in argi:
        if arg == "-p":
            outputs.append(next(arg))

        elif arg in ["--help", "-h"]:
            print(HELP)
            sys.exit()

        elif arg == "--protocol":
            config["protocol"] = next(arg)

        elif arg == "--inputs":
            inputs = next(argi)[1:-1].split(",")
            if isinstance(inputs, str):
                channels[0] = [int(inputs)]
            else:
                channels[0] = list(map(int, inputs))

        elif arg == "--outputs":
            outputs = next(argi)[1:-1].split(",")
            if isinstance(outputs, str):
                channels[1] = [int(outputs)]
            else:
                channels[1] = list(map(int, outputs))

        # Positional args
        elif config["method"] is None:
            config["method"] = arg
            break

        # elif config["protocol"] and config["operation"] is None:
        #     config["operation"] = arg

        elif arg == "--":
            continue

        else:
            print(HELP)
            sys.exit()


    return sub_parsers[config.pop("method")](argi, config), outputs


def main():
    parse_args(sys.argv)

if __name__ == "__main__":
    main()
    sys.exit()
