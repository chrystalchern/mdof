import numpy as np
from matplotlib import pyplot as plt
import scienceplots

plt.style.use(["science", "no-latex"])# , "notebook"])

nln = "\n"

def print_modes(modes, Tn=None, zeta=None):

    header = "         T         \N{Greek Small Letter Zeta}"

    if Tn is not None:
        header += "          T % error"

    if zeta is not None:
        header += "    \N{Greek Small Letter Zeta} % error"

    print("Spectral quantities:")
    print(f"{header}")
    for mode in modes.values():
        f = mode["freq"]
        z = mode["damp"]
        row = f"      {1/f: <9.4}  {z: <9.4}"
        if Tn is not None:
            row += f"    {100*(1/f-Tn)/(Tn): <9.4}"
        if zeta is not None:
            row += f"    {100*(z-zeta)/zeta: <9.4}"
        print(row)

def plot_models(models, Tn, zeta):
    fig, ax = plt.subplots(2, 3, constrained_layout=True, sharex=True, figsize=(13,6.5))

    period = [models[method]["period"][0] for method in models]
    ax[0,0].bar(["true"]+list(models), [Tn]+period)
    ax[0,0].set_title("Periods")# , fontsize=14)
    ax[0,0].set_ylabel("Period (s)")# , fontsize=13)

    period_errors = [100*(p-Tn)/Tn for p in period]
    ax[1,0].bar(models.keys(), period_errors, color=None, edgecolor="k", linewidth=0.5)
    ax[1,0].set_title("Period Errors")# , fontsize=14)
    ax[1,0].set_ylabel("Percent Error (%)")# , fontsize=13)
    ax[1,0].set_xlabel("Method")# , fontsize=13)

    damp = [models[method]["damping"][0] for method in models]
    ax[0,1].bar(["true"]+list(models), [zeta]+damp)
    ax[0,1].set_title("Damping")# , fontsize=14)
    ax[0,1].set_ylabel("Damping Ratio")# , fontsize=13)

    damping_errors = [100*(d-zeta)/zeta for d in damp]
    ax[1,1].bar(models.keys(), damping_errors, color=None, edgecolor="k", linewidth=0.5)
    ax[1,1].set_title("Damping Errors")# , fontsize=14)
    ax[1,1].set_ylabel("Percent Error (%)")# , fontsize=13)
    ax[1,1].set_xlabel("Method")# , fontsize=13)

    ax[0,2].axis('off')

    times_list = [models[method]["time"] for method in models]
    ax[1,2].bar(models.keys(), times_list, color=None, edgecolor="k", linewidth=0.5)
    ax[1,2].set_title("Runtime")# , fontsize=14)
    ax[1,2].set_ylabel("time (s)")# , fontsize=13)
    ax[1,2].set_xlabel("Method")# , fontsize=13)

    # for axi in fig.axes:
    #     axi.tick_params(axis='x', labelsize=12)
    #     axi.tick_params(axis='y', labelsize=12)

    for i,error in zip([0,1,2],[period_errors,damping_errors,times_list]):
        rects = ax[1,i].patches
        for rect, label in zip(rects, error):
            label = f"{label: <5.3}"
            height = rect.get_height()
            ax[1,i].text(
                rect.get_x() + rect.get_width() / 2, height, label, ha="center", va="bottom"
            )
    
    fig.suptitle("Spectral Quantity Prediction with System Identification") #,fontsize=14)

def plot_io(input, output, t, title=None):
    fig, ax = plt.subplots(1,2,figsize=(12,4))
    # fig, ax = plt.subplots(1,2,figsize=(8,3))
    ax[0].plot(t,input)
    ax[0].set_xlabel("time (s)")# , fontsize=13)
    ax[0].set_ylabel("input")# , fontsize=13)
    ax[1].plot(t,output)
    ax[1].set_xlabel("time (s)")# , fontsize=13)
    ax[1].set_ylabel("output")# , fontsize=13)
    fig.suptitle(title, fontsize=14)

def plot_pred(ytrue, models, t, title=None):
    fig, ax = plt.subplots(figsize=(8,4))
    ax.plot(t,ytrue,label="true")
    if type(models) is np.ndarray:
        ax.plot(t,models,"--",label=f"prediction")
    else:
        for method in models:
            ax.plot(t,models[method]["ypred"],"--",label=method)
    ax.set_xlabel("time (s)")# , fontsize=13)
    ax.set_ylabel("output")# , fontsize=13)
    fig.legend(fontsize=12, frameon=True, framealpha=1)    
    fig.suptitle(title, fontsize=14)