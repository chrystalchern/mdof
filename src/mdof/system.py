import warnings

from mdof import markov, realize
from .process import decimate
from .realization import Realization

def system(inputs, outputs, dt=None, **options):
    """
    State space system realization from input and output data.

    :param inputs:      input time history. dimensions:
                        :math:`(q,nt)`, where :math:`q` = number of
                        inputs, and :math:`nt` = number of timesteps
    :type inputs:       array
    
    :param outputs:     output response history.
                        dimensions: :math:`(p,nt)`, where
                        :math:`p` = number of outputs, and
                        :math:`nt` = number of timesteps
    :type outputs:      array
    
    :param method:      (optional) system identification method.
                        default is "srim", other options are
                        "okid-era", "okid-era-dc", and
                        "deterministic".
    :type method:       string
    
    :param dt:          (optional) timestep, in seconds.
                        Stored on the returned realization, adjusted
                        for any decimation.
    :type dt:           float
    
    :param decimation:  (optional) decimation factor. default: no
                        decimation
    :type decimation:   int

    :return:            state-space realization ``(A,B,C,D)`` of
                        :class:`mdof.realization.Realization`.
                        Unpacks and indexes as a tuple but carries
                        the effective ``dt`` and a provenance record.
    :rtype:             :class:`mdof.realization.Realization`
    """

    method = options.get("method", "srim")
    decimation = options.get("decimation", None)

    if method not in {
        "srim",
        "okid-era",
        "okid-era-dc",
        "deterministic"
    }: raise ValueError(f"Unknown method {method}")

    if decimation is not None and dt is None:
        warnings.warn(
            "`decimation` was set without passing `dt`. The "
            "identified realization changes the effective "
            "timestep to `dt * decimation`, but with no `dt` "
            "the resulting model cannot record its own "
            "timestep, so modal analysis may be silently off "
            "by the decimation factor. Pass `dt=...` so mdof "
            "can track it.",
            FutureWarning,
            stacklevel=2,
        )

    if decimation is not None:
        inputs = decimate(inputs, decimation=decimation)
        outputs = decimate(outputs, decimation=decimation)

    if method == "okid-era":
        Y = markov.okid(inputs, outputs, **options)
        realization = realize.era(Y, **options)

    if method == "okid-era-dc":
        Y = markov.okid(inputs, outputs, **options)
        realization = realize.era_dc(Y, **options)

    if method == "srim":
        realization = realize.srim(inputs, outputs, **options)

    if method == "deterministic":
        realization = realize.deterministic(inputs, outputs,
                                            **options)

    effective_dt = None
    if dt is not None:
        effective_dt = dt * (
            decimation if decimation is not None else 1)

    A,B,C,D = realization

    return Realization(
        A,B,C,D,
        dt=effective_dt,
        provenance={
            "method": method,
            "decimation": decimation,
            "dt_raw": dt,
        },
    )