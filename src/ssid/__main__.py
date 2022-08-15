from .srim import parse_srim
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

def parse_brace2(args, config):
    pass

def parse_args(args):
    outputs = []
    parsers = {
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
        elif arg == "--setup":
            install_me()
            sys.exit()
        elif arg in ["--help", "-h"]:
            print(HELP)
            sys.exit()
        elif arg == "--protocol":
            config["protocol"] = next(arg)

        # Positional args
        elif config["method"] is None:
            config["method"] = arg
            break

        elif config["protocol"] and config["operation"] is None:
            config["operation"] = arg

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
        elif arg == "--":
            continue

        else:
            break

    return parsers[arg](argi, config), outputs

if __name__ == "__main__":
    import sys
    import quakeio
    import numpy as np
    from pathlib import Path
    method = None

    config, out_ops = parse_args(sys.argv)
    sys.exit()

    channels = config.get("channels", [[17, 3, 20], [9, 7, 4]])

    if config["method"] == "test":
        data_dir = Path("RioDell_Petrolia_Processed_Data")

        first_input = quakeio.read(data_dir/f"CHAN{channels[0][0]:03d}.V2")
        npoints = len(first_input.accel.data)
        inputs, outputs = np.zeros((2,npoints,len(channels[0])))

        # Inputs
        inputs[:,0] = first_input.accel.data
        for i,inp in enumerate(channels[0][1:]):
            inputs[:,i+1] = quakeio.read(data_dir/f"CHAN{inp:03d}.V2").accel.data

        # Outputs
        for i,inp in enumerate(channels[1]):
            outputs[:,i] = quakeio.read(data_dir/f"CHAN{inp:03d}.V2").accel.data

        dt = first_input.accel["time_step"]
        config["dt"] = dt

    elif "event_file" in config:
        event = quakeio.read(config["event_file"])
        try:
            inputs = [
                event.match("l", station_channel=f"{i}").accel.data for i in channels[0]
            ]
        except AttributeError:
            print(f"{channels[0]}", "\n"*3, file=sys.stderr)
            raise

        try:
            outputs = [
                event.match("l", station_channel=f"{i}").accel.data for i in channels[1]
            ]
        except AttributeError:
            print(f"{channels[1]}", "\n"*3, file=sys.stderr)
            raise

        #npoints = len(inputs[:,0])
        dt = event.match("l", station_channel=f"{channels[0][0]}").accel["time_step"]
        config["dt"] = dt


    A,B,C,D = srim(inputs, outputs, **config)

    freqdmpSRIM, modeshapeSRIM, *_ = ComposeModes(dt, A, B, C, D)

    if not out_ops:
        print(f"period: {np.real(1/freqdmpSRIM[:,0])}")
    elif "freq" in out_ops:
        print(f"frequency: {freqdmpSRIM[:,0]}")
    elif "cycl" in out_ops:
        print(f"cyclic_frequency: {2*np.pi*freqdmpSRIM[:,0]}")


