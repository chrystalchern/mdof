
.. py:function:: ssid.modes(realization,**options)

    Generate a state space system realization from a given set of input and output data using a specified system identification method.

    :param realization: realization in the form of state space coefficients ``A,B,C,D``
    :type realization: list

    :return: ``freq`` natural frequencies, ``damp`` damping ratios, and ``mode`` mode shapes
    :rtype: dictionary
