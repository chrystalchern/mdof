from ssid import markov, realize

def system(input, output, **options):
    # if not method:
    #     method = "okid-era-dc"
    y = output
    f = input
    m = 300
    no = 140
    nc = 140
    a = 0
    b = 0
    l = 10
    g = 3
    r_okid_era_dc = 2
    Y = markov.okid(y, f, m)
    A,B,C,D = realize.era_dc(Y, no, nc, a, b, l, g, r=r_okid_era_dc)

    return A,B,C,D