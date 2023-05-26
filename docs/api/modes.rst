.. py:function:: ssid.modes.modes(realization)

    Determine the modal parameters of a given state space system realization.

    :param realization: realization in the form of state space coefficients ``(A,B,C,D)``
    :type realization: tuple
    :param dt: timestep
    :type realization: float
    
    :return: natural frequencies, damping ratios, mode shapes, and eigenvalue condition numbers for each mode in the ``freq``, ``damp``, ``modeshape``, and ``cnd`` fields, respectively.
    :rtype: dictionary
