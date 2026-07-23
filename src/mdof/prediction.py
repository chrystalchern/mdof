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

def get_error(
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


# Backwards-compatible re-exports
from mdof.utilities.events import Event, Prediction, norm, _get_error
Realization = Event
