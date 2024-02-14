import numpy as np
from matplotlib import pyplot as plt
import plotly.graph_objects as go

try:
    import scienceplots
    plt.style.use(["poster"])# ,"science", "notebook"])
except ImportError:
    pass
except OSError:
    plt.style.use("notebook")

nln = "\n"

DEFAULT_PLOTLY_COLORS=[
    '#1f77b4',  # muted blue
    '#ff7f0e',  # safety orange
    '#2ca02c',  # cooked asparagus green
    '#d62728',  # brick red
    '#9467bd',  # muted purple
    '#8c564b',  # chestnut brown
    '#e377c2',  # raspberry yogurt pink
    '#7f7f7f',  # middle gray
    '#bcbd22',  # curry yellow-green
    '#17becf'   # blue-teal
]

def print_modes(modes, Tn=None, zeta=None):

    header = "       T(s)        \N{Greek Small Letter Zeta}        EMACO      MPC       EMACO*MPC"

    if Tn is not None:
        header += "     T % error"

    if zeta is not None:
        header += "    \N{Greek Small Letter Zeta} % error"

    print("Spectral quantities:")
    print(f"{header}")
    for mode in sorted(modes.values(), key=lambda x: x["freq"]):
        f = mode["freq"]
        z = mode["damp"]
        emaco = mode["energy_condensed_emaco"]
        mpc = mode["mpc"]
        row = f"      {1/f: <9.4}  {z: <9.4}  {emaco: <9.4}  {mpc: <9.4}  {emaco*mpc: <9.4}"
        if Tn is not None:
            row += f"    {100*(1/f-Tn)/(Tn): <9.4}"
        if zeta is not None:
            row += f"    {100*(z-zeta)/zeta: <9.4}"
        print(row)
    print("Mean Period(s):", np.mean([1/v["freq"] for v in modes.values()]))
    print("Standard Dev(s):", np.std([1/v["freq"] for v in modes.values()]))

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

    for axi in fig.axes:
        axi.tick_params(axis='x', rotation=45)#, labelsize=12)
        # axi.tick_params(axis='y', labelsize=12)

    for i,error in zip([0,1,2],[period_errors,damping_errors,times_list]):
        rects = ax[1,i].patches
        for rect, label in zip(rects, error):
            label = f"{label: <5.3}"
            height = rect.get_height()
            ax[1,i].text(
                rect.get_x() + rect.get_width() / 2, height, label, ha="center", va="bottom"
            )
    
    fig.suptitle("Spectral Quantity Prediction with System Identification",fontsize=16)

def plot_io(inputs, outputs, t, title=None):
    fig, ax = plt.subplots(1,2,figsize=(8,3),constrained_layout=True)
    # fig, ax = plt.subplots(1,2,figsize=(8,3))
    color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
    if len(inputs.shape) > 1:
        for i in range(inputs.shape[0]):
            ax[0].plot(t,inputs[i,:])
    else:
        ax[0].plot(t,inputs)
    ax[0].set_xlabel("time (s)")# , fontsize=13)
    ax[0].set_ylabel("inputs")# , fontsize=13)
    if len(outputs.shape) > 1:
        for i in range(outputs.shape[0]):
            ax[1].plot(t,outputs[i,:], color=color_cycle[0], alpha=0.5)
    else:
        ax[1].plot(t,outputs)
    ax[1].set_xlabel("time (s)")# , fontsize=13)
    ax[1].set_ylabel("outputs")# , fontsize=13)
    fig.suptitle(title, fontsize=16)

def plot_pred(ytrue, models, t, title=None):
    color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
    fig, ax = plt.subplots(figsize=(8,4))
    if len(ytrue.shape) > 1:
        for i in range(ytrue.shape[0]):
            ax.plot(t,ytrue[i,:],label=["true" if i==0 else None][0],color=color_cycle[0],alpha=0.5)
    else:
        ax.plot(t,ytrue,label="true")
    if type(models) is np.ndarray:
        if len(models.shape) > 1:
            for i in range(models.shape[0]):
                ax.plot(t,models[i,:],"--",label=["prediction" if i==0 else None][0],color=color_cycle[1],alpha=0.5)
        else:
            ax.plot(t,models,"--",label=f"prediction")
    else:
        for method in models:
            if len(models[method]["ypred"].shape) > 1:
                for i in range(models[method]["ypred"].shape[0]):
                    ax.plot(t,models[method]["ypred"][i,:],"--",label=method)
            else:
                ax.plot(t,models[method]["ypred"],"--",label=method)
    ax.set_xlabel("time (s)")# , fontsize=13)
    ax.set_ylabel("outputs")# , fontsize=13)
    fig.legend(fontsize=12, frameon=True, framealpha=1)    
    ax.set_title(title) #, fontsize=14)

def plot_transfer(models, title=None, labels=None, plotly=False):
    if plotly:
        import plotly.graph_objects as go
        layout = go.Layout(
            title=title,
            xaxis=dict(
                title="Period (s)"
            ),
            yaxis=dict(
                title="Amplitude"
            ),
            width=600, height=300,
            margin=dict(l=70, r=20, t=20, b=20))
        fig = go.Figure(layout=layout)
        if type(models) is np.ndarray:
            if len(models.shape) > 2:
                for i in range(models.shape[0]):
                    fig.add_trace(go.Scatter(x=models[i,0],y=models[i,1]/max(models[i,1]),name=labels[i]))
            else:
                fig.add_trace(go.Scatter(x=models[0],y=models[1]/max(models[1]),name=labels))
        else:
            for method in models:
                fig.add_trace(go.Scatter(x=models[method][0],y=models[method][1]/max(models[method][1]),name=method))
        fig.show(renderer="notebook_connected")
    else:
        fig, ax = plt.subplots(figsize=(6,3))
        if type(models) is np.ndarray:
            if len(models.shape) > 2:
                for i in range(models.shape[0]):
                    ax.plot(models[i,0],models[i,1]/max(models[i,1]),label=labels[i])
            else:
                ax.plot(models[0],models[1]/max(models[1]),label=labels)
        else:
            for method in models:
                ax.plot(models[method][0],models[method][1]/max(models[method][1]),label=method)
        ax.set_xlabel("Period (s)")
        ax.set_ylabel("Amplitude")
        if labels is not None:
            ax.legend()#fontsize=12)
        ax.set_title(title)#, fontsize=14)

class FrequencyContent:
    def __init__(self, scale, period, xlabel, ylabel, xlimits) -> None:
        self.scale = scale
        self.period = period
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.xlimits = xlimits
        self.num_traces = 0
        layout = go.Layout(
            title=None,
            xaxis=dict(
                title=self.xlabel,
                range=xlimits
            ),
            yaxis=dict(
                title=self.ylabel
            ),
            width=500, height=300,
            margin=dict(l=70, r=20, t=20, b=20))
        self.fig = go.Figure(layout=layout)

    def add(self, periods, amplitudes=None, label=None):
        if self.period:
            x_data = periods
        else:
            x_data = 1/periods
        if amplitudes is not None:
            if self.scale:
                y_data = amplitudes/max(amplitudes)
            else:
                y_data = amplitudes
            self.fig.add_trace(go.Scatter(x=x_data,y=y_data,mode='lines',name=label,showlegend=True,line_color=DEFAULT_PLOTLY_COLORS[self.num_traces%len(DEFAULT_PLOTLY_COLORS)]))
        else:
            N = len(x_data)
            x_data_vlines = np.array([x_data,x_data,np.full(N,None)]).transpose().reshape(-1)
            y_data_vlines = np.array([np.zeros(N),np.ones(N),np.full(N,None)]).transpose().reshape(-1)
            self.fig.add_trace(go.Scatter(x=x_data_vlines,y=y_data_vlines,mode='lines',name=label,showlegend=True,line_dash="dash",line_color=DEFAULT_PLOTLY_COLORS[self.num_traces%len(DEFAULT_PLOTLY_COLORS)]))
        self.num_traces += 1        

    def get_figure(self):
        return self.fig

def make_hover_data(data, ln=None):
    import numpy as np
    if ln is None:
        items = np.array([d.values for d in data])
        keys = data[0].keys()
    else:
        items = np.array([list(data.values())]*ln)
        keys = data.keys()
    return {
        "hovertemplate": "<br>".join(f"{k}: %{{customdata[{v}]}}" for v,k in enumerate(keys)),
        "customdata": list(items),
    }

if __name__ == "__main__":

    import numpy as np
    periods = np.array([0.1,0.3,0.5,0.7])
    amplitudes = np.array([0.1,0.2,1,0.2])

    plot = FrequencyContent(scale=True, period=True, xlabel="Period (s)", ylabel="Amplitude", xlimits=[0,0.5])

    plot.add(periods, amplitudes, label="R1")
    plot.add(periods[:2], label="SRIM")
    plot.add(periods[2:], label="SRIM")
    plot.add(periods, amplitudes-0.05, label="R1")


    fig:go.Figure = plot.get_figure()
    print(fig.to_json())
