from mdof import markov, realize, modal
from .numerics import decimate

def system(inputs, outputs, **options):
    """
    State space system realization from input and output data.

    :param inputs:      input time history. dimensions: :math:`(q,nt)`, where
                        :math:`q` = number of inputs, and :math:`nt` = number of timesteps
    :type inputs:       array
    :param outputs:     output response history.
                        dimensions: :math:`(p,nt)`, where :math:`p` = number of outputs, and
                        :math:`nt` = number of timesteps
    :type outputs:      array
    :param method:      system identification method. default is "srim", other options are "okid-era" and "okid-era-dc".
    :type method:       string, optional
    :param decimation:  decimation factor. default: 8
    :type decimation:   int, optional

    :return:            system realization in the form of state space coefficients ``(A,B,C,D)``
    :rtype:             tuple of arrays
    """
    method = options.get("method", "srim")
    decimation = options.get("decimation", None)

    if method not in {
        "srim",
        "okid-era",
        "okid-era-dc"
    }: raise ValueError(f"Unknown method {method}")

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

    return realization



def modes(inputs, outputs, dt, **options):
    """
    Modal identification using state space system realization from input and output data.

    :param inputs:      input time history. dimensions: :math:`(q,nt)`, where
                        :math:`q` = number of inputs, and :math:`nt` = number of timesteps
    :type inputs:       array
    :param outputs:     output response history.
                        dimensions: :math:`(p,nt)`, where :math:`p` = number of outputs, and
                        :math:`nt` = number of timesteps
    :type outputs:      array
    :param dt:          timestep.
    :type dt:           float
    :param method:      system identification method. default is "srim", other options are
                        "okid-era" and "okid-era-dc".
    :type method:       string, optional
    :param decimation:  decimation factor. default: 8
    :type decimation:   int, optional

    :return:            system modes, including natural frequencies, damping ratios, mode shapes,
                        condition numbers, and modal validation metrics EMAC and MPC.
    :rtype:             dictionary
    """
    realization = system(inputs, outputs, **options)

    return modal.system_modes(realization, dt, options["decimation"])

