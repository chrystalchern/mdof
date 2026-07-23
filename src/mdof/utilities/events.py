"""
Seismic-event workflow built on top of the core :class:`mdof.Realization`.

An :class:`Event` identifies a state-space model from the recorded input and
output channels of a seismic event. It delegates identification,
stabilization, and mode-removal bookkeeping to :class:`mdof.Realization`,
and adds domain-specific conveniences: channel extraction, intensity windowing,
unit reporting, and prediction/validation against another event.

See :mod:`mdof.realization` and :mod:`mdof.modal` for the general
state-space machinery.
"""

import numpy as np

import mdof
from mdof.utilities import extract_channels
from mdof.utilities.testing import (
    intensity_bounds,
    truncate_by_bounds,
    align_signals,
)

# cm/s/s -> g
_CM_PER_S2_TO_G = 0.0010197162129779


def norm(y, ntime, order: int = 1, averaging_mode='mean'):
    if averaging_mode is None:
        n = np.linalg.norm(y, order)
        n_manual = (np.sum(np.abs(y)**order))**(1/order)
        assert n == n_manual
    elif averaging_mode == 'mean':
        n = (np.sum(np.abs(y)**order)/ntime)**(1/order)
    return n


def _get_error(ytrue, ypred, numerator_norm, denominator_norm,
               numerator_averaged=True, denominator_averaged=True):
    ntime = len(ytrue)
    assert len(ytrue) == len(ypred), \
        f"{len(ytrue)=}, {len(ypred)=}, but they must match"
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
        numerator = norm(ytrue-ypred, ntime, order=numerator_norm,
                         averaging_mode=num_mode)
    elif numerator_norm in ['max', 'supremum', 'infinity', 'uniform',
                            'inf', np.inf]:
        numerator = np.max(np.abs(ytrue - ypred))
    else:
        raise ValueError(f"Unsupported numerator_norm: {numerator_norm}")

    if isinstance(denominator_norm, int):
        denominator = norm(ytrue, ntime, order=denominator_norm,
                           averaging_mode=den_mode)
    elif denominator_norm in ['max', 'supremum', 'infinity', 'uniform',
                              'inf', np.inf]:
        denominator = np.max(np.abs(ytrue))
    else:
        raise ValueError(f"Unsupported denominator_norm: {denominator_norm}")

    return numerator/denominator


class Event:
    """
    Identify a state-space model from a recorded seismic event.

    :param event:       a parsed event record (e.g. from ``quakeio``) whose
                        channels are addressable by ``chan_conf``.
    :param chan_conf:   mapping with ``'inputs'`` and ``'outputs'`` channel
                        lists.
    :param conf:        identification options forwarded to
                        :func:`mdof.system` (e.g. ``order``, ``horizon``).
    :param stabilize:   filter unstable discrete modes. default: ``True``.
    :param windowed_intensity: truncate to the strong-motion window before
                        identification. default: ``False``.
    :param response:    channel response to extract. default: ``'accel'``.
    :param cm:          data are in cm/s/s (converted to g for reporting).

    The identified model is a core :class:`mdof.Realization`, available as
    :attr:`realization`; stabilization details are read from its provenance.
    """

    def __init__(self, event, chan_conf, conf,
                 stabilize: bool = True,
                 windowed_intensity: bool = False,
                 verbose: bool = False,
                 response='accel',
                 cm: bool = True):

        if verbose:
            units = _CM_PER_S2_TO_G if (cm and response == 'accel') else 1
            print("peak acceleration (g):",
                  np.abs(event['peak_accel']*units))
            print("event date/time:", event['event_date'])

        inputs, dti = extract_channels(event, chan_conf['inputs'],
                                       response=response)
        outpts, dto = extract_channels(event, chan_conf['outputs'],
                                       response=response)
        assert abs(dti - dto) < 1e-5, \
            "The input timestep is different from the output timestep."

        if windowed_intensity:
            bounds = intensity_bounds(outpts[0], lb=0.001, ub=0.999)
            inputs = truncate_by_bounds(inputs, bounds)
            outpts = truncate_by_bounds(outpts, bounds)

        # Identify system
        realization = mdof.system(method='srim', inputs=inputs,
                                  outputs=outpts, dt=dti, **conf)

        # Delegate stabilization + mode-removal bookkeeping to the core object.
        if stabilize:
            realization = realization.stabilize(verbose=verbose)

        self._realization = realization
        self._inputs = inputs
        self._outpts = outpts
        self._dt = dti
        self._chan_conf = chan_conf
        self._windowed_intensity = windowed_intensity
        self._response = response
        self._cm = cm
        self._verbose = verbose

    # -- read-through to the core realization --------------------------------

    @property
    def realization(self):
        """The underlying core :class:`mdof.Realization`."""
        return self._realization

    @property
    def dt(self):
        return self._dt

    @property
    def stabilized(self):
        """Whether unstable modes were filtered out."""
        return bool(self._realization.provenance.get("stabilized", False))

    @property
    def modes_removed(self):
        """
        Modes removed during stabilization as ``{"indices": [...],
        "periods": [...]}``, or ``None`` if the model was not stabilized.
        """
        return self._realization.provenance.get("modes_removed", None)

    def modes(self, **options):
        """Modal properties of the identified system (see
        :meth:`mdof.Realization.modes`)."""
        return self._realization.modes(**options)

    def predict(self, event, **options):
        """Validate the identified model against another ``event``."""
        return Prediction(
            self._realization, event, self._chan_conf,
            self._windowed_intensity, self._response,
            verbose=self._verbose, **options,
        )


class Prediction:
    def __init__(self, realization, event, chan_conf, windowed_intensity,
                 response='accel', align=True, verbose=False):
        self._windowed_intensity = windowed_intensity

        inputs, dti = extract_channels(event, chan_conf['inputs'],
                                       response=response)
        outpts, dto = extract_channels(event, chan_conf['outputs'],
                                       response=response)
        assert abs(dti - dto) < 1e-5, \
            "The input timestep is different from the output timestep."

        self._dt = dti

        if windowed_intensity:
            bounds = intensity_bounds(outpts[0], lb=0.001, ub=0.999)
            inputs = truncate_by_bounds(inputs, bounds)
            outpts = truncate_by_bounds(outpts, bounds)

        self._inputs = inputs

        times = np.linspace(0, self._dt*outpts.shape[1], outpts.shape[1])

        out_pred = mdof.simulate.simulate(system=realization, inputs=inputs)
        if align:
            _, outpts, out_pred, times = align_signals(
                signal1=outpts, signal2=out_pred, times=times)

        self._times = times
        self._outpts = outpts
        self._out_pred = out_pred

    def _truncate(self, series):
        return truncate_by_bounds(series, intensity_bounds(self._outpts[0]))

    def channels(self):
        return range(self._outpts.shape[0])

    def time(self, ic):
        """ic: index of channel"""
        if self._windowed_intensity:
            return self._times
        return self._truncate(self._times)

    def test(self, ic):
        """ic: index of channel"""
        if self._windowed_intensity:
            return self._out_pred[ic]
        return self._truncate(self._out_pred[ic])

    def true(self, ic):
        """ic: index of channel"""
        if self._windowed_intensity:
            return self._outpts[ic]
        return self._truncate(self._outpts[ic])

    def error(self, ic, metric):
        """ic: index of channel"""
        return _get_error(
            ytrue=self.true(ic),
            ypred=self.test(ic),
            numerator_norm=2,
            denominator_norm=2,
            numerator_averaged=True,
            denominator_averaged=True,
        )
