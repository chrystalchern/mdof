
.. py:function:: ssid.system(output,input,**options)

    Generate a state space system realization from a given set of input and output data using a specified system identification method.

    :param output: output acceleration response history to an arbitrary input. dimensions: :math:`(p,nt)`, where :math:`p` = number of output channels, and :math:`nt` = number of timesteps
    :type output: array
    :param input: input acceleration time history. dimensions: :math:`(q,nt)`, where :math:`q` = number of input channels, and :math:`nt` = number of timesteps

    :return: realization in the form of state space coefficients ``A,B,C,D``
    :rtype: list[array]
