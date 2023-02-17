import sys
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

    return sub_parsers[arg](argi, config), outputs


def main():
    config, out_ops = parse_args(sys.argv)

if __name__ == "__main__":
    main()
    sys.exit()
