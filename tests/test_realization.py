"""
Tests for the core :class:`mdof.Realization` object and its wiring
into ``mdof.system`` / ``mdof.sysid``:

1. **Backward compatibility** -- a ``Realization`` unpacks, indexes,
   and slices exactly like the legacy ``(A, B, C, D)`` tuple.
2. **dt safety** -- the effective (post-decimation) timestep is
   carried on the object, so modal analysis cannot be silently off by
   the decimation factor.
"""

import numpy as np
import pytest

import mdof
from mdof.realization import Realization


def _toy_realization(n=2, q=1, p=1, dt=None, provenance=None):
    A = np.eye(n) * 0.9
    B = np.ones((n, q))
    C = np.ones((p, n))
    D = np.zeros((p, q))
    return Realization(A, B, C, D, dt=dt, provenance=provenance or {})


# --- (A, B, C, D) tuple behavior -----------

def test_tuple_unpacking_returns_the_matrices():
    r = _toy_realization()
    A, B, C, D = r
    assert A is r.A and B is r.B and C is r.C and D is r.D


def test_integer_indexing_matches_tuple_positions():
    r = _toy_realization()
    assert r[0] is r.A
    assert r[3] is r.D


def test_slicing_is_supported():
    r = _toy_realization()
    assert r[1:] == (r.B, r.C, r.D)


def test_len_and_iteration():
    r = _toy_realization()
    assert len(r) == 4
    assert list(r) == [r.A, r.B, r.C, r.D]


def test_simulate_accepts_a_realization():
    # mdof.simulate.simulate unpacks (A, B, C, D) = system
    r = _toy_realization(n=2, q=1, p=1)
    inputs = np.random.default_rng(0).standard_normal((1, 50))
    out = mdof.simulate.simulate(r, inputs)
    assert out.shape[1] == inputs.shape[1]


# --- stabilize(): immutability + provenance bookkeeping ------------

def test_stabilize_records_removed_modes_in_provenance():
    # A has one unstable discrete mode (|eig| > 1) and one
    # stable mode.
    A = np.diag([1.5, 0.5])
    r = Realization(A, np.ones((2, 1)), np.ones((1, 2)),
                    np.zeros((1, 1)), dt=0.02)
    rs = r.stabilize()

    assert rs.provenance["stabilized"] is True
    removed = rs.provenance["modes_removed"]
    assert 0 in removed["indices"]        # the unstable mode
    # dt is known, so periods are recorded alongside indices
    assert len(removed["periods"]) == len(removed["indices"])


def test_stabilize_leaves_original_unchanged():
    A = np.diag([1.5, 0.5])
    r = Realization(A, np.ones((2, 1)), np.ones((1, 2)),
                    np.zeros((1, 1)), dt=0.02)
    r.stabilize()
    # frozen + returns a new object: original A and provenance
    # untouched
    assert np.allclose(r.A, np.diag([1.5, 0.5]))
    assert "stabilized" not in r.provenance


def test_stabilize_without_dt_omits_periods():
    A = np.diag([1.5, 0.5])
    r = Realization(A, np.ones((2, 1)), np.ones((1, 2)),
                    np.zeros((1, 1)))
    rs = r.stabilize()
    assert "periods" not in rs.provenance["modes_removed"]


# --- dt safety -----------------------------------------------------

def test_modes_requires_a_timestep():
    r = _toy_realization(dt=None)
    with pytest.raises(ValueError):
        r.modes()


def test_modes_uses_the_objects_dt():
    # modes() must call system_modes with the object's own
    # (effective) dt.
    captured = {}
    r = _toy_realization(dt=0.05)

    import mdof.modal as modal
    orig = modal.system_modes

    def spy(realization, dt, *args, **kwargs):
        captured["dt"] = dt
        return {}

    modal.system_modes = spy
    try:
        r.modes()
    finally:
        modal.system_modes = orig

    assert captured["dt"] == pytest.approx(0.05)


def test_modes_rejects_decimation_argument():
    # dt already reflects decimation, so passing it into modes()
    # is an error.
    r = _toy_realization(dt=0.05)
    with pytest.raises(TypeError):
        r.modes(decimation=5)


# --- wiring through sysid --------------------------------------

def test_sysid_returns_a_realization_with_effective_dt():
    rng = np.random.default_rng(1)
    # Generate data from a known stable SISO system.
    A = np.array([[0.6, 0.2], [0.0, 0.5]])
    B = np.array([[1.0], [0.5]])
    C = np.array([[1.0, 0.0]])
    D = np.zeros((1, 1))
    inputs = rng.standard_normal((1, 400))
    outputs = mdof.simulate.simulate((A, B, C, D), inputs)

    dt_raw, decimation = 0.02, 2
    result = mdof.sysid(inputs, outputs, dt=dt_raw,
                        decimation=decimation, method="srim")

    assert isinstance(result, Realization)
    # still unpackable
    Ai, Bi, Ci, Di = result
    assert Ai.ndim == 2
    # effective dt reflects decimation
    assert result.dt == pytest.approx(dt_raw * decimation)
    assert result.provenance["method"] == "srim"
    assert result.provenance["decimation"] == decimation


def test_decimation_without_dt_warns():
    rng = np.random.default_rng(3)
    u = rng.standard_normal((1, 300))
    A = np.array([[0.6, 0.2], [0.0, 0.5]])
    B = np.array([[1.0], [0.5]])
    C = np.array([[1.0, 0.0]])
    D = np.zeros((1, 1))
    y = mdof.simulate.simulate((A, B, C, D), u)

    with pytest.warns(FutureWarning):
        mdof.sysid(u, y, decimation=2)


def test_decimation_with_dt_does_not_warn():
    rng = np.random.default_rng(4)
    u = rng.standard_normal((1, 300))
    A = np.array([[0.6, 0.2], [0.0, 0.5]])
    B = np.array([[1.0], [0.5]])
    C = np.array([[1.0, 0.0]])
    D = np.zeros((1, 1))
    y = mdof.simulate.simulate((A, B, C, D), u)

    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("error", FutureWarning)
        mdof.sysid(u, y, dt=0.01, decimation=2)  # must not raise


# --- modes() decimation correctness / deprecation ------------------

def _siso_data(seed, n=2000):
    rng = np.random.default_rng(seed)
    A = np.array([[0.6, 0.2], [0.0, 0.5]])
    B = np.array([[1.0], [0.5]])
    C = np.array([[1.0, 0.0]])
    D = np.zeros((1, 1))
    u = rng.standard_normal((1, n))
    y = mdof.simulate.simulate((A, B, C, D), u)
    return u, y


def test_top_level_modes_uses_effective_dt_under_decimation():
    # Regression: mdof.modes(..., decimation=k) must use the effective
    # (post-decimation) timestep, i.e. agree with Realization.modes().
    u, y = _siso_data(seed=0)
    dt, decimation = 0.01, 2

    P_modes, _ = mdof.modes(u, y, dt, decimation=decimation)
    R_modes = mdof.sysid(u, y, dt=dt, decimation=decimation).modes()
    P_obj = [1 / m["freq"]
             for m in R_modes.values()]

    assert sorted(P_modes) == pytest.approx(sorted(P_obj))


def test_top_level_modes_does_not_warn():
    # dt is threaded through, so neither system() nor
    # system_modes should emit the decimation FutureWarning.
    u, y = _siso_data(seed=1)
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("error", FutureWarning)
        mdof.modes(u, y, 0.01, decimation=2)


def test_system_modes_decimation_kwarg_is_deprecated():
    u, y = _siso_data(seed=2)
    r = mdof.sysid(u, y, dt=0.01)
    from mdof import modal
    with pytest.warns(FutureWarning):
        modal.system_modes(r.to_tuple(), 0.01, decimation=2)
