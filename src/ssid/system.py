from ssid import markov, realize
from control.matlab import lsim
from matplotlib import pyplot as plt

def system(input, output, **options):
    # if not method:
    #     method = "okid-era-dc"
    y = output
    f = input
    m = 900
    no = 425
    nc = 425
    r = 20
    a = 0
    b = 0
    l = 10
    g = 3
    r_okid_era_dc = 10
    Y = markov.okid(y, f, m)
    A,B,C,D = realize.era_dc(Y, no, nc, a, b, l, g, r=r_okid_era_dc)

    return A,B,C,D

def plot_prediction(realization):
    A,B,C,D,dt = realization
    y_ssid = lsim(ss(A,B,C,D,dt),f,t)[0]
    fig, ax = plt.subplots(figsize=(8,4))
    ax.plot(t,y,label="original")
    ax.plot(t,y_ssid,"--",label=f"system id")
    ax.set_xlabel("time (s)", fontsize=12)
    ax.set_ylabel(r"response $y(t)$", fontsize=12)
    ax.set_title("Displacement Response", fontsize=14)
    ax.legend(fontsize=12);