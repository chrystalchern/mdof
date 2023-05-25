``markov`` Module
=================

.. py:function:: ssid.markov.okid(input, output, m)

    Identify Markov parameters, or discrete impulse response data, for a given set of input and output data, using the Observer Kalman Identification Algorithm (OKID) (Juang, Phan, Horta, Longman, 1993).

    :param output: output acceleration response history to an arbitrary input. dimensions: :math:`(p,nt)`, where :math:`p` = number of output channels, and :math:`nt` = number of timesteps
    :type output: array
    :param input: input acceleration time history. dimensions: :math:`(q,nt)`, where :math:`q` = number of input channels, and :math:`nt` = number of timesteps
    :type input: array
    :param m: number of Markov paramters to compute
    :type m: int

    :return: the Markov parameters, with dimensions :math:`(p,q*nt)`
    :rtype: array

