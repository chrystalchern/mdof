import mdof
from mdof.utilities import extract_channels
from mdof.validation import stabilize_discrete
from mdof.utilities.testing import intensity_bounds, truncate_by_bounds, align_signals

import numpy as np


def _get_error(true, test, metric='l2_norm', normalized=True):
    """
    """
    if not isinstance(true, np.ndarray):
        true = np.array(true)
    if not isinstance(test, np.ndarray):
        test = np.array(test)
    if true.ndim != 1:
        Warning("The true series has incorrect dimensions for this operation (must be a 1D array or list).")
    if test.ndim != 1:
        Warning("The test series has incorrect dimensions for this operation (must be a 1D array or list).")
    assert true.shape == test.shape, f"Shapes are different for true series ({true.shape}) and test series ({test.shape})."
    if metric == 'l2_norm':
        error = np.linalg.norm(test-true)
        if normalized:
            return error/np.linalg.norm(true)
        return error
    if metric == 'l2_norm_max_normalized':
        error = np.linalg.norm(test-true)
        if normalized:
            return error/max(true)
        return error
    if metric == 'are':
        error = sum(np.abs(test)-np.abs(true))
        if normalized:
            return error/sum(np.abs(true))
        return error
    if metric == 'sym':
        error = sum(test-true)
        if normalized:
            return error/(sum(test+true)/2)
        return error
    else:
        return NotImplementedError("This metric is not recognized. The metrics available are 'l2_norm', 'are', and 'sym'.")


def _get_period_from_discrete_A_eigval(eigval, dt):
    return (2*np.pi) / np.abs( (np.log(eigval))/dt )


class Realization:
    def __init__(self, event, 
                 chan_conf, 
                 conf, 
                 stabilize: bool = True, 
                 windowed_intensity: bool = False,
                 verbose: bool = False,
                 response = 'accel',
                 cm: bool = True):
                 # Assumes units of cm/s/s. Converts to g.
                 # TODO: make unit non-specific.

        if verbose:
            units = 0.0010197162129779 if (cm and response=='accel') else 1
            peak_accel = np.abs(event['peak_accel']*units)
            event_date = event['event_date']
            print("peak acceleration (g):", peak_accel)
            print("event date/time:", event_date)

        inputs, dti = extract_channels(event, chan_conf['inputs'], response=response)
        outpts, dto = extract_channels(event, chan_conf['outputs'], response=response)
        assert abs(dti - dto) < 1e-5, "The input timestep is different from the output timestep."

        if windowed_intensity:
            bounds = intensity_bounds(outpts[0], lb=0.001, ub=0.999)
            inputs = truncate_by_bounds(inputs, bounds)
            outpts = truncate_by_bounds(outpts, bounds)

        realization = mdof.system(method='srim', inputs=inputs, outputs=outpts, **conf)

        if stabilize:
            A_stable, real_mode_indices = stabilize_discrete(A=realization[0], verbose=verbose, list_filtered_modes=True)
            modes_removed = np.array([real_mode_indices,[
                _get_period_from_discrete_A_eigval(np.linalg.eig(realization[0])[0][2*i], dti) for i in real_mode_indices
            ]])
            realization = (A_stable,*realization[1:])
            self._modes_removed = modes_removed
        else:
            self._modes_removed = None
        
        self._stabilized = stabilize
        self._inputs = inputs 
        self._outpts = outpts 
        self._realization = realization
        self._dt = dti
        self._chan_conf = chan_conf
        self._windowed_intensity = windowed_intensity
        self._response = response
        self._cm = cm
        self._verbose = verbose

    def predict(self, event):
        return Prediction(self._realization, event, self._chan_conf, self._windowed_intensity, self._response, verbose=self._verbose)



class Prediction:
    def __init__(self, realization, event, chan_conf, windowed_intensity, response='accel', align=True, verbose=False):
        self._windowed_intensity = windowed_intensity

        inputs, dti = extract_channels(event, chan_conf['inputs'], response=response)
        outpts, dto = extract_channels(event, chan_conf['outputs'], response=response)
        assert abs(dti - dto) < 1e-5, "The input timestep is different from the output timestep."

        self._dt = dti

        if windowed_intensity:
            bounds = intensity_bounds(outpts[0], lb=0.001, ub=0.999)
            inputs = truncate_by_bounds(inputs, bounds)
            outpts = truncate_by_bounds(outpts, bounds)

        self._inputs = inputs

        times = np.linspace(0,self._dt*outpts.shape[1],outpts.shape[1])

        out_pred = mdof.simulate.simulate(system=realization, inputs=inputs)
        if align:
            max_lag, outpts, out_pred, times = align_signals(signal1=outpts, signal2=out_pred, times=times)

        if False:
            print(f"Signals aligned by cross-correlation with a lag of up to {max_lag} timesteps.")

        self._times = times
        self._outpts = outpts
        self._out_pred = out_pred


    def _truncate(self, series):
        return (truncate_by_bounds(series, intensity_bounds(self._outpts[0])))

    def channels(self):
        return range(self._outpts.shape[0]) 

    def time(self, ic):
        """
        ic: index of channel 
        """
        if self._windowed_intensity:
            return self._times
        else:
            return self._truncate(self._times)


    def test(self, ic):
        """
        ic: index of channel 
        """
        if self._windowed_intensity:
            return self._out_pred[ic]
            # print(f"Ch {chan_conf['outputs'][j]}: {METRIC_NAMES[metric]} = {error}")
        else:
            return self._truncate(self._out_pred[ic])


    def true(self, ic):
        """
        ic: index of channel 
        """
        if self._windowed_intensity:
            return self._outpts[ic]
        else:
            return self._truncate(self._outpts[ic])


    def error(self, ic, metric):
        """
        ic: index of channel 
        """
        return _get_error(true=self.true(ic), test=self.test(ic), metric=metric)
    