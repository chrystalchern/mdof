import ast
import sys
import json

import ssid
import ssid.modal
import quakeio
import numpy as np
# from .okid import parse_okid

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


# def _decimate(series, decimation, **options):
#     import numpy as np
#     return series[np.arange(0,len(series),decimation)]


# def extract_channels(event, channels, decimate=1, permissive=True):
def extract_channels(event, channels, permissive=True):
    import numpy as np
    import sys

    def find(chan):
        #event.match("r", file_name=f".*{chan}\.[vV]2")
        return event.match("l", station_channel=f"{chan}")

    data = np.asarray([
        # _decimate(find(chan).accel.data, d=decimate)
        find(chan).accel.data
        for chan in channels if find(chan) is not None
    ])

    if len(data) == 0:
        raise ValueError(f"No channels found, requested '{channels}'")

    elif len(data) != len(channels):
        if permissive:
            print(f"Only extracted {len(data)} channels, {len(channels)-len(data)} missing.", file=sys.stderr)
        else:
            raise ValueError("Could not extract all desired channels")

    # dt = find(channels[0]).accel["time_step"]*decimate
    dt = find(channels[0]).accel["time_step"]
    return data, dt



def parse_time(argi, config, method=None):
    help = f"""\
    ssid {method} <event>


    --inputs  <int>
    --outputs <int>

    Options
    --threads <int>
    """
    damping  = None
    decimate = 1
    channels = [None, None]
    for arg in argi:
        if arg == "--threads" and method == "response":
            config["threads"] = int(next(argi))

        elif arg == "--inputs":
            channels[0] = next(argi)

        elif arg == "--outputs":
            channels[1] = next(argi)

        elif arg == "--damping" and method == "response":
            damp = next(argi)
            try:
                config["damping"] = [float(damp)]
            except:
                config["damping"] = ast.literal_eval(damp)
        
        elif arg == "--periodband":
            config["period_band"] = tuple(float(x) for x in (next(argi).split(" ")))

        elif arg == "-h" or arg == "--help":
            print(help)
            sys.exit()

        else:
            config["event_file"] = arg

    event = quakeio.read(config["event_file"], exclusions=["filter*"])

    try:
        # inputs,  dt = extract_channels(event, [channels[0]], decimate=decimate)
        inputs,  dt = extract_channels(event, [channels[0]])
        # outputs, dt = extract_channels(event, [channels[1]], decimate=decimate)
        outputs, dt = extract_channels(event, [channels[1]])
    except Exception as e:
        print(json.dumps({"error": str(e), "data": []}))
        return

    import ssid.spec
    from ssid.modal import spectrum_modes
    f = {
        "response": ssid.spec.response_transfer,
        "fourier":  ssid.spec.fourier_transfer,
    }[method]

    periods, amplitudes = spectrum_modes(
                          *f(inputs=inputs.flatten(), outputs=outputs.flatten(), step = dt, **config)
                        )
    output = [
            {"period": period, "amplitude": amplitude}
            for period, amplitude in zip(periods, amplitudes)
    ]
    print(json.dumps({"data": output}, cls=JSON_Encoder, indent=4))


def parse_srim(argi, config, method=None):
    help="""
    SRIM -- System Identification with Information Matrix

    Parameters
    no          order of the observer Kalman ARX filter (formerly p).
    r           size of the state-space model used for
                representing the system. (formerly orm/n)
    """
    config.update({"p"  :  5, "orm":  4})
    decimate = 8

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
            channels[0] = ast.literal_eval(next(argi))

        elif arg == "--outputs":
            channels[1] = ast.literal_eval(next(argi))

        elif arg == "--decimate":
            config["decimate"] = int(next(argi))

        elif arg == "--":
            continue

        else:
            config["event_file"] = arg

    event = quakeio.read(config["event_file"], exclusions=["filter*"])
    try:
        # inputs,  dt = extract_channels(event, [channels[0]], decimate=decimate)
        inputs,  dt = extract_channels(event, channels[0])
        # outputs, dt = extract_channels(event, [channels[1]], decimate=decimate)
        outputs, dt = extract_channels(event, channels[1])
    except Exception as e:
        print(json.dumps({"error": str(e), "data": []}))
        return

    config["dt"] = dt

    try:
        realization = ssid.system(inputs=inputs, outputs=outputs, full=True, method=method, **config)
        ss_modes = ssid.modal.system_modes(realization,dt,decimation=config.get("decimate",1))

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
        "srim":     parse_srim,
        "response": parse_time,
        "fourier":  parse_time,
        "test": parse_srim,
        "okid": parse_srim,
        "okid": parse_srim
    }
    method = None
    config = {
            "operation": None,
            "protocol": None
    }

    argi = iter(args[1:])
    for arg in argi:
        if arg == "-p":
            outputs.append(next(arg))

        elif arg in ["--help", "-h"]:
            print(HELP)
            sys.exit()

        elif arg == "--protocol":
            config["protocol"] = next(arg)

        # Positional args
        elif method is None:
            method = arg
            break

        elif arg == "--":
            continue

        else:
            print(HELP)
            sys.exit()


    return sub_parsers[method](argi, config, method=method), outputs


def main():
    parse_args(sys.argv)

if __name__ == "__main__":
    main()
    sys.exit()
