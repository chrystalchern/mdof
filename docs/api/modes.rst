.. py:function:: ssid.modes(realization,**options)

    Determine the modal parameters of a given state space system realization.

    :param realization: realization in the form of state space coefficients ``A,B,C,D``
    :type realization: list
    :param dt: realization in the form of state space coefficients ``A,B,C,D``
    :type realization: list
    
    :return: ``freq`` natural frequencies, ``damp`` damping ratios, and ``mode`` mode shapes
    :rtype: dictionary
