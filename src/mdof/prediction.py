import mdof
from mdof.utilities import extract_channels
from mdof.validation import stabilize_discrete
from mdof.utilities.testing import intensity_bounds, truncate_by_bounds, align_signals
from typing import Optional, Literal

import numpy as np


def _parse_norm(ord_like):
    """
    Map user-friendly norm spec to numpy-compatible 'ord':
    accepts 1,2,'L1','L2','max','supremum','infinity' (-> 1,2,np.inf).
    """
    if ord_like in (1, '1', 'L1', 'l1'):
        return 1
    if ord_like in (2, '2', 'L2', 'l2'):
        return 2
    if ord_like in ('max', 'supremum', 'infinity', 'inf', np.inf):
        return np.inf
    raise ValueError(f"Unsupported norm: {ord_like}")

def _vector_norm(x: np.ndarray, ord_like=2, axis=-1) -> np.ndarray:
    """
    Vector norm along 'axis' without altering units:
    - L2: sqrt(sum(x^2)),  L1: sum(|x|),  L∞: max(|x|)
    """
    ord_val = _parse_norm(ord_like)
    if ord_val == 2:
        return np.sqrt(np.sum(x * x, axis=axis))
    elif ord_val == 1:
        return np.sum(np.abs(x), axis=axis)
    elif ord_val == np.inf:
        return np.max(np.abs(x), axis=axis)
    
def get_error_new(
    y_true,
    y_pred,
    *,
    normalized: bool = True,
    denominator_norm = 2,   # 1,2,'L1','L2','max','supremum','infinity'
    numerator_norm = 2,     # 1,2,'L1','L2','max','supremum','infinity'
    averaged: bool = False,
    axis: Optional[int] = -1,   # None => flatten and treat as one vector
    averaging_mode: Optional[Literal['mean','median','max','quantile','weighted','global']] = None,
    weights: Optional[np.ndarray] = None,  # for 'weighted' mode
    quantile: float = 0.95                 # for 'quantile' mode
):
    y_true = np.asarray(y_true)
    y_pred = np.asarray(y_pred)
    if y_true.shape != y_pred.shape:
        raise ValueError(f"Shapes differ: {y_true.shape} vs {y_pred.shape}")
    
    # decide reduction axis
    if axis is None:
        # treat the entire array as one vector (single sample)
        reduce_axis = tuple(range(y_true.ndim)) if y_true.ndim > 0 else None
    else:
        reduce_axis = axis

    # per-sample numerator norm
    diff = y_pred - y_true
    num = _vector_norm(diff, numerator_norm, axis=reduce_axis)

    # optional per-sample normalization (units cancel here)
    den = _vector_norm(y_true, denominator_norm, axis=reduce_axis)
    err = num / den

    # if not averaging, return per-sample array
    if not averaged and averaging_mode is None:
        return err

    # choose/resolve averaging mode
    if averaging_mode is None:
        # default scheme: L∞ -> max; otherwise -> mean
        if _parse_norm(numerator_norm) == np.inf:
            averaging_mode = 'max'
        else:
            averaging_mode = 'mean'

    # aggregate across the remaining sample dimensions
    # (after removing `reduce_axis`, 'err' holds one scalar per sample)
    if averaging_mode == 'mean':
        return float(np.mean(err))
    elif averaging_mode == 'median':
        return float(np.median(err))
    elif averaging_mode == 'max':
        return float(np.max(err))
    elif averaging_mode == 'quantile':
        return float(np.quantile(err, quantile))
    elif averaging_mode == 'weighted':
        if weights is None:
            raise ValueError("weights must be provided for averaging_mode='weighted'.")
        w = np.asarray(weights, dtype=float)
        if w.shape != err.shape:
            raise ValueError(f"weights shape {w.shape} must match err shape {err.shape}.")
        ws = w / (np.sum(w) if np.sum(w) != 0 else 1.0)
        return float(np.sum(ws * err))
    elif averaging_mode == 'global':
        # Global normalization (optional alternative): combine numerators/denominators before dividing.
        # Note: this is also dimensionless when normalized=True.
        if normalized:
            # recompute per-sample norms and sum them globally
            num_sum = np.sum(_vector_norm(diff, numerator_norm, axis=reduce_axis))
            den_sum = np.sum(_vector_norm(y_true, denominator_norm, axis=reduce_axis))
            return float(num_sum / max(den_sum))
        else:
            # global absolute error (sum of per-sample norms)
            return float(np.sum(_vector_norm(diff, numerator_norm, axis=reduce_axis)))
    else:
        raise ValueError(f"Unknown averaging_mode: {averaging_mode}")


def norm(y, ntime, order:int=1, averaging_mode='mean'):
    if averaging_mode is None:
        norm = np.linalg.norm(y, order)
        norm_manual = (np.sum(np.abs(y)**order))**(1/order)
        assert norm == norm_manual
    elif averaging_mode == 'mean':
        norm = (np.sum(np.abs(y)**order)/ntime)**(1/order)
    return norm

def _get_error(ytrue, ypred, numerator_norm, denominator_norm, numerator_averaged=True, denominator_averaged=True):
    ntime = len(ytrue)
    assert (len(ytrue) == len(ypred)), f"{len(ytrue)=}, {len(ypred)=}, but they must match"
    if numerator_norm == 'L1':
        numerator_norm = 1
    if numerator_norm == 'L2':
        numerator_norm = 2
    if denominator_norm == 'L1':
        denominator_norm = 1
    if denominator_norm == 'L2':
        denominator_norm = 2

    num_mode = 'mean' if numerator_averaged else None
    den_mode = 'mean' if denominator_averaged else None

    if isinstance(numerator_norm, int):
        numerator = norm(ytrue-ypred, ntime, order=numerator_norm, averaging_mode=num_mode)
    elif numerator_norm in ['max', 'supremum', 'infinity', 'uniform', 'inf', np.inf]:
        numerator = np.max(np.abs(ytrue - ypred))
    else:
        raise ValueError(f"Unsupported numerator_norm: {numerator_norm}")
    
    if isinstance(denominator_norm, int):
        denominator = norm(ytrue, ntime, order=denominator_norm, averaging_mode=den_mode)
    elif denominator_norm in ['max', 'supremum', 'infinity', 'uniform', 'inf', np.inf]:
        denominator = np.max(np.abs(ytrue))
    else:
        raise ValueError(f"Unsupported denominator_norm: {denominator_norm}")
    
    return numerator/denominator


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
            A_stable, indices = stabilize_discrete(A=realization[0], verbose=verbose, list_filtered_modes=True)
            eigvals = np.linalg.eigvals(realization[0])
            modes_removed = np.array([
                indices,
                [_get_period_from_discrete_A_eigval(eigvals[i], dti) for i in indices]
            ])
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
        #return _get_error(true=self.true(ic), test=self.test(ic), metric=metric)
        return _get_error(
            ytrue=self.true(ic),
            ypred=self.test(ic),
            numerator_norm=2,
            denominator_norm=2,
            numerator_averaged=True,
            denominator_averaged=True
        )
    