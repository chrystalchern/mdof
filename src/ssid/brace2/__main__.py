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


def parse_args(args):
    outputs = []
    parsers = {
        "srim": parse_srim,
        "test": parse_srim,
        "okid": parse_okid
    }
    config = {}
    argi = iter(args[1:])
    for arg in argi:
        if arg == "-p":
            outputs.append(next(arg))
        if arg == "--setup":
            install_me()
            sys.exit()
        elif arg in ["--help", "-h"]:
            print(HELP)
            sys.exit()
        else:
            config["method"] = arg
            return parsers[arg](argi, config), outputs

if __name__ == "__main__":
    import sys
    import quakeio
    import numpy as np
    from pathlib import Path
    channels = [[17, 3, 20], [9, 7, 4]]
    method = None

    config, out_ops = parse_args(sys.argv)

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
        inputs = np.array([
            event.at(file_name=f"CHAN{i:03d}.V2").accel.data for i in channels[0]
        ]).T
        outputs = np.array([
            event.at(file_name=f"CHAN{i:03d}.V2").accel.data for i in channels[1]
        ]).T
        npoints = len(inputs[:,0])
        dt = event.at(file_name=f"CHAN{channels[0][0]:03d}.V2").accel["time_step"]
        config["dt"] = dt

    #print(config)

    A,B,C,D = srim(inputs, outputs, **config)

    freqdmpSRIM, modeshapeSRIM, *_ = ComposeModes(dt, A, B, C, D)

    if not out_ops:
        print(f"period: {np.real(1/freqdmpSRIM[:,0])}")
    elif "freq" in out_ops:
        print(f"frequency: {freqdmpSRIM[:,0]}")
    elif "cycl" in out_ops:
        print(f"cyclic_frequency: {2*np.pi*freqdmpSRIM[:,0]}")


