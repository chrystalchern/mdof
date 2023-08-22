
.. py:function:: ssid.system(method="srim",*args,**options)

    Generate a state space system realization from a given set of input and output data using a specified system identification method.

    :param method: system identification method. default is "srim", other options are "okid-era" and "okid-era-dc".
    :type method: string
    :param input: input time history data. dimensions: :math:`(q,nt)`, where :math:`q` = number of input channels, and :math:`nt` = number of timesteps
    :type input: array
    :param output: output response history data. dimensions: :math:`(p,nt)`, where :math:`p` = number of output channels, and :math:`nt` = number of timesteps
    :type output: array

    :return: system realization in the form of state space coefficients ``(A,B,C,D)``
    :rtype: tuple
