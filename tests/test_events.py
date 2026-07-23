"""
Tests for the seismic-event workflow (:mod:`mdof.utilities.events`).

`extract_channels` is mocked so these tests don't require real
recorded event data: synthetic input/output signals from a known
SISO system are fed in. The focus is that `Event` builds on the core
:class:`mdof.Realization` and reads stabilization bookkeeping off its
provenance rather than tracking it.
"""

import numpy as np
import pytest

import mdof
from mdof.realization import Realization
import mdof.utilities.events as events


CHAN_CONF = {"inputs": ["1"], "outputs": ["2"]}
DT = 0.02


def _known_siso(seed=0, n=800):
    rng = np.random.default_rng(seed)
    A = np.array([[0.6, 0.2], [0.0, 0.5]])
    B = np.array([[1.0], [0.5]])
    C = np.array([[1.0, 0.0]])
    D = np.zeros((1, 1))
    u = rng.standard_normal((1, n))
    y = mdof.simulate.simulate((A, B, C, D), u)
    return u, y


@pytest.fixture
def mock_channels(monkeypatch):
    u, y = _known_siso()

    def fake_extract(event, channels, response="accel"):
        if channels == CHAN_CONF["inputs"]:
            return u, DT
        return y, DT

    monkeypatch.setattr(events, "extract_channels", fake_extract)
    return u, y


def _event():
    return {"peak_accel": 1.0, "event_date": "2026-01-01"}


def test_event_wraps_core_realization(mock_channels):
    ev = events.Event(_event(), CHAN_CONF, conf={}, stabilize=False)
    assert isinstance(ev.realization, Realization)
    assert ev.dt == pytest.approx(DT)
    # dt threaded through identification -> carried on the core object
    assert ev.realization.dt == pytest.approx(DT)


def test_event_reads_stabilization_from_provenance(mock_channels):
    ev = events.Event(_event(), CHAN_CONF, conf={}, stabilize=True)
    assert ev.stabilized is True
    # modes_removed comes straight from the core
    # realization's provenance
    assert (ev.modes_removed
            is ev.realization.provenance["modes_removed"])
    assert "indices" in ev.modes_removed


def test_event_not_stabilized_has_no_modes_removed(mock_channels):
    ev = events.Event(_event(), CHAN_CONF, conf={}, stabilize=False)
    assert ev.stabilized is False
    assert ev.modes_removed is None


def test_event_modes_delegates_to_core(mock_channels):
    ev = events.Event(_event(), CHAN_CONF, conf={}, stabilize=True)
    m = ev.modes()
    assert isinstance(m, dict) and len(m) > 0


def test_backcompat_alias_from_prediction_module(mock_channels):
    # old import path still works and points at the renamed class
    from mdof.prediction import Realization as SeismicRealization
    assert SeismicRealization is events.Event
    ev = SeismicRealization(
        _event(), CHAN_CONF, conf={}, stabilize=True)
    assert isinstance(ev.realization, Realization)
