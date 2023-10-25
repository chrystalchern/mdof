from mdof import modes
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
    :param orders:      range of system orders. default: (2,100,2)
    :type orders:       tuple (`start`, `stop`, `step`), optional
    :param method:      system identification method. default is "srim", other options are "okid-era" and "okid-era-dc".
    :type method:       string, optional
    :param decimation:  decimation factor. default: 1
    :type decimation:   int, optional

    :return:            system realization in the form of state space coefficients ``(A,B,C,D)``
    :rtype:             tuple of arrays
    """
    orders = options.get("orders", (2,100,2)) 
    mode_properties = np.array([[order, 1/mode["freq"], mode["damp"]] for order in range(*orders) for mode in modes(inputs, outputs, dt, order=order, **options).values()])
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


