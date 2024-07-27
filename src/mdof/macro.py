from mdof import system, modal
from matplotlib import pyplot as plt
import numpy as np

def stabilization(inputs, outputs, dt, plotly=False, **options):
    """
    Stabilization diagram of modal properties from state space system identification, based on model order.

    :param inputs:      input time history. dimensions: :math:`(q,nt)`, where
                        :math:`q` = number of inputs, and :math:`nt` = number of timesteps
    :type inputs:       array
    :param outputs:     output response history.
                        dimensions: :math:`(p,nt)`, where :math:`p` = number of outputs, and
                        :math:`nt` = number of timesteps
    :type outputs:      array
    :param orders:      (optional) range of system orders. default: (2,100,2)
    :type orders:       tuple (`start`, `stop`, `step`)
    :param method:      (optional) system identification method. default is "srim", other options are "okid-era" and "okid-era-dc".
    :type method:       string

    :return:            system realization in the form of state space coefficients ``(A,B,C,D)``
    :rtype:             tuple of arrays
    """
    orders = options.get("orders", (2,100,2)) 
    mode_properties = np.array([[order, 1/mode["freq"], mode["damp"]] 
                                for order in range(*orders) 
                                for mode in modal.system_modes(system(inputs, outputs, order=order, **options), dt, **options).values()])
    if plotly:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
        fig = make_subplots(rows=1, cols=2, shared_yaxes=True, y_title="Model Order")
        fig.add_trace(go.Scatter(x=mode_properties[:,1],y=mode_properties[:,0], mode="markers", marker=dict(size=8)),
                      row=1, col=1)
        fig.add_trace(go.Scatter(x=mode_properties[:,2],y=mode_properties[:,0], mode="markers", marker=dict(size=8)),
                      row=1, col=2)
        fig.update_layout(height=300, margin=dict(l=70, r=20, t=20, b=20), showlegend=False)
        fig.layout.xaxis.title="Period (s)"
        fig.layout.xaxis2.title="Damping Ratio"
        fig.show(renderer="notebook_connected")
    else:
        fig, ax = plt.subplots(2,1,figsize=(6,6))
        ax[0].plot(mode_properties[:,1],mode_properties[:,0],'o',markersize=4)
        ax[0].set_title("Periods")
        ax[0].set_ylabel("Model Order")
        ax[1].plot(mode_properties[:,2],mode_properties[:,0],'o',markersize=4)
        ax[1].set_title("Damping Ratios")
        ax[1].set_ylabel("Model Order")

    return fig


DEFAULT_PLOTLY_COLORS=[
    '#1f77b4',  # muted blue
    '#ff7f0e',  # safety orange
    '#2ca02c',  # cooked asparagus green
    '#d62728',  # brick red
    '#9467bd',  # muted purple
    # '#8c564b',  # chestnut brown
    '#e377c2',  # raspberry yogurt pink
    '#7f7f7f',  # middle gray
    '#bcbd22',  # curry yellow-green
    '#17becf'   # blue-teal
]
from itertools import cycle
import plotly.graph_objects as go

class FrequencyContent:
    def __init__(self, scale, period, xlabel, ylabel, **options) -> None:
        self.scale = scale
        self.period = period
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.xlimits = options.get('xlimits', (0.1,3))
        self.colors = options.get('colors', cycle(DEFAULT_PLOTLY_COLORS))
        self.num_traces = 0
        layout = go.Layout(
            title=None,
            xaxis=dict(
                title=self.xlabel,
                range=self.xlimits
            ),
            yaxis=dict(
                title=self.ylabel
            ),
            width=700, height=300,
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
            self.fig.add_trace(go.Scatter(x=x_data,y=y_data,mode='lines',name=label,showlegend=True,line_color=next(self.colors)))
        else:
            N = len(x_data)
            x_data_vlines = np.array([x_data,x_data,np.full(N,None)]).transpose().reshape(-1)
            y_data_vlines = np.array([np.zeros(N),np.ones(N),np.full(N,None)]).transpose().reshape(-1)
            self.fig.add_trace(go.Scatter(x=x_data_vlines,y=y_data_vlines,mode='lines',name=label,showlegend=True,line_dash="dash",line_color=next(self.colors)))
        self.num_traces += 1        

    def get_figure(self):
        return self.fig
