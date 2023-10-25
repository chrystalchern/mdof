# gets called when running mdof from the command line:
#    python -m mdof ...
#
import ast
import sys
import json

import mdof
import mdof.modal
import quakeio
import numpy as np
# from .okid import parse_okid

HELP = """
mdof [-p|w <>...] <method> <event> <inputs> <outputs>
mdof [-p|w <>...] <method> -i <inputs> -o <outputs>

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

    # dt = find(channels["inputs"]).accel["time_step"]*decimate
    dt = find(channels[0]).accel["time_step"]
    return data, dt



def parse_time(argi, config, channels, method=None):
    help = f"""\
    mdof {method} <event>


    --inputs  <int>
    --outputs <int>

    Options
    --threads <int>
    """
    damping  = None
    for arg in argi:
        if arg == "--threads" and method == "response":
            config["threads"] = int(next(argi))

        elif arg == "--config":
            conf_arg = json.loads(next(argi))
            for i in conf_arg.pop("channels", []):
                channels[i["type"]+"s"].append(i["id"])
            for i in {"inputs", "outputs"}:
                if i in conf_arg:
                    channels[i] = conf_arg.pop(i)
            config.update(conf_arg)

        elif arg == "--inputs":
            channels["intpus"] = next(argi)

        elif arg == "--outputs":
            channels["outputs"] = next(argi)

        elif arg == "--decimate":
            config["decimate"] = int(next(argi))

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
        # inputs,  dt = extract_channels(event, [channels["inputs"]], decimate=decimate)
        inputs,  dt = extract_channels(event, [channels["inputs"]])
        # outputs, dt = extract_channels(event, [channels["outputs"]], decimate=decimate)
        outputs, dt = extract_channels(event, [channels["outputs"]])
    except Exception as e:
        print(json.dumps({"error": str(e), "data": []}))
        return

    import mdof.transform
    from mdof.modal import spectrum_modes
    f = {
        "response": mdof.transform.response_transfer,
        "fourier":  mdof.transform.fourier_transfer,
    }[method]

    periods, amplitudes = spectrum_modes(
                          *f(inputs=inputs.flatten(), outputs=outputs.flatten(), step = dt, **config)
                        )
    output = [
            {"period": period, "amplitude": amplitude}
            for period, amplitude in zip(periods, amplitudes)
    ]
    print(json.dumps({"data": output}, cls=JSON_Encoder, indent=4))


def parse_srim(argi, config, channels, method=None):
    help="""
    SRIM -- System Identification with Information Matrix

    Parameters
    no          order of the observer Kalman ARX filter (formerly p).
    r           size of the state-space model used for
                representing the system. (formerly orm/n)
    """
#   config.update({"p"  :  5, "orm":  4})

    for arg in argi:
        if arg == "--arx-order":
            config["no"] = int(next(argi))

        elif arg == "--config":
            conf_arg = json.loads(next(argi))
            for i in conf_arg.pop("channels", []):
                channels[i["type"]+"s"].append(i["id"])
            for i in {"inputs", "outputs"}:
                if i in conf_arg:
                    channels[i] = conf_arg.pop(i)
            config.update(conf_arg)

        elif arg == "--dt":
            config["dt"] = float(next(argi))

        elif arg == "--ss-size":
            config["r"] = int(next(argi))

        elif arg in ["--help", "-h"]:
            print(help)
            sys.exit()

        elif arg == "--inputs":
            channels["inputs"] = json.loads(next(argi))

        elif arg == "--outputs":
            channels["outputs"] = json.loads(next(argi))

        elif arg == "--decimate":
            config["decimate"] = int(next(argi))

        elif arg == "--":
            continue

        else:
            config["event_file"] = arg

    event = quakeio.read(config["event_file"], exclusions=["filter*"])
    try:
        # inputs,  dt = extract_channels(event, [channels["inputs"]], decimate=decimate)
        inputs,  dt = extract_channels(event, channels["inputs"])
        # outputs, dt = extract_channels(event, [channels["outputs"]], decimate=decimate)
        outputs, dt = extract_channels(event, channels["outputs"])
    except Exception as e:
        print(json.dumps({"error": str(e), "data": []}))
        return

    config["dt"] = dt

    try:
        realization = mdof.system(inputs=inputs, outputs=outputs, full=True, method=method, **config)
        ss_modes = mdof.modal.system_modes(realization,dt,decimation=config.get("decimate",1))

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
    config = {}
    channels = {
#           "operation": None,
#           "protocol": None,
            "inputs":  [],
            "outputs": []
    }

    argi = iter(args[1:])
    for arg in argi:
        if arg == "-p":
            outputs.append(next(arg))

        elif arg == "--config":
            conf_arg = json.loads(next(argi))
            method = conf_arg.pop("method", method)
            for i in conf_arg.pop("channels", []):
                channels[i["type"]+"s"].append(i["id"])
            config.update(conf_arg)

        elif arg in ["--help", "-h"]:
            print(HELP)
            sys.exit()

#       elif arg == "--protocol":
#           config["protocol"] = next(arg)

        # Positional args
        elif method is None:
            method = arg
            break

        elif arg == "--":
            continue

        else:
            print(HELP)
            sys.exit()


    return sub_parsers[method](argi, config, channels, method=method), outputs


def main():
    parse_args(sys.argv)

if __name__ == "__main__":
    main()
    sys.exit()

