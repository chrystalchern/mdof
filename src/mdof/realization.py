"""
Core state-space realization object for ``mdof``.

A :class:`Realization` bundles the state-space coefficient matrices
``(A, B, C, D)`` together with the effective timestep ``dt`` of the
identified model and a ``provenance`` record describing how it was produced
(identification method, decimation factor, stabilization, etc.).

The object is **tuple-compatible**: it iterates, indexes, and
slices exactly like a plain ``(A, B, C, D)`` tuple. Usage such as::

    A,B,C,D = mdof.sysid(inputs, outputs)

is allowed.
"""

from dataclasses import dataclass, field, replace
from typing import Optional

import numpy as np


@dataclass(frozen=True)
class Realization:
    """
    A discrete-time state-space realization ``(A, B, C, D)``.

    :param A:           state transition matrix, shape
                        :math:`(n, n)`.
    :param B:           input matrix, shape :math:`(n, q)`.
    :param C:           output matrix, shape :math:`(p, n)`.
    :param D:           feedthrough matrix, shape :math:`(p, q)`.
    :param dt:          effective timestep of the realization, in
                        seconds. When decimated data is used for
                        identification, this is the *post-decimation*
                        timestep (``dt_raw * decimation``).
                        ``None`` if unknown.
    :param provenance:  free-form record of how the realization was
                        produced, e.g. ``{"method": "srim",
                        "decimation": 5, "stabilized": True}``.
    """

    A: np.ndarray
    B: np.ndarray
    C: np.ndarray
    D: np.ndarray
    dt: Optional[float] = None
    provenance: dict = field(default_factory=dict)

    # -- tuple compatibility -------------------------------------------------

    def _matrices(self):
        return (self.A, self.B, self.C, self.D)

    def __iter__(self):
        return iter(self._matrices())

    def __len__(self):
        return 4

    def __getitem__(self, index):
        # Supports both integer indexing (realization[0]) and slicing
        # (realization[1:]).
        return self._matrices()[index]

    def to_tuple(self):
        """Return the plain ``(A, B, C, D)`` tuple."""
        return self._matrices()

    # -- convenience methods -------------------------------------------------

    @property
    def order(self):
        """State-space order :math:`n` (number of states)."""
        return self.A.shape[0]

    def _require_dt(self, dt):
        dt = self.dt if dt is None else dt
        if dt is None:
            raise ValueError(
                "No timestep available. Either identify with a `dt` "
                "(e.g. mdof.sysid(inputs, outputs, dt=...)) or pass dt=... here."
            )
        return dt

    def modes(self, dt=None, **options):
        """
        Modal properties (frequencies, damping, modeshapes).

        :param dt: override the stored timestep, in seconds.
        :return: dictionary of modes, as returned by :func:`mdof.modal.system_modes`.
        """
        from mdof import modal
        dt = self._require_dt(dt)
        # `dt` already reflects any decimation, so `decimation` here would
        # double-count. Reject it loudly rather than silently ignore it.
        if "decimation" in options:
            raise TypeError(
                "`Realization.modes()` does not accept `decimation`: the "
                "realization's `dt` already reflects any decimation used "
                "during identification."
            )
        return modal.system_modes(self._matrices(), dt, **options)

    def predict(self, inputs, **options):
        """
        Simulate the output response for a given ``inputs`` array.

        :param inputs: input time history, shape :math:`(q, nt)`.
        :return: output response history, shape :math:`(p, nt)`.
        """
        from mdof.simulate import simulate
        return simulate(self._matrices(), inputs, **options)

    def stabilize(self, **options):
        """
        Return a new :class:`Realization` with unstable discrete modes filtered
        out (see :func:`mdof.validation.stabilize_discrete`).

        Provenance is updated with ``{"stabilized": True}``; the original object
        is left unchanged (this object is immutable).
        """
        from mdof.validation import stabilize_discrete
        A_stable = stabilize_discrete(A=self.A, **options)
        new_provenance = {**self.provenance, "stabilized": True}
        return replace(self, A=A_stable, provenance=new_provenance)

    def __repr__(self):
        p, n = self.C.shape
        q = self.B.shape[1] if self.B.ndim == 2 else 1
        dt = "unknown" if self.dt is None else f"{self.dt:g}s"
        prov = ", ".join(f"{k}={v}" for k, v in self.provenance.items())
        prov = f", {prov}" if prov else ""
        return (
            f"Realization(order={n}, inputs={q}, outputs={p}, "
            f"dt={dt}{prov})"
        )
