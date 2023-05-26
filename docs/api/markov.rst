``markov`` Module
=================

.. py:function:: ssid.markov.okid(input,output,**options)

    Identify Markov parameters, or discrete impulse response data, for a given set of input and output data, using the Observer Kalman Identification Algorithm (OKID) (Juang, Phan, Horta, Longman, 1993).

    :param input: input time history data. dimensions: :math:`(q,nt)`, where :math:`q` = number of input channels, and :math:`nt` = number of timesteps
    :type input: array
    :param output: output response history data. dimensions: :math:`(p,nt)`, where :math:`p` = number of output channels, and :math:`nt` = number of timesteps
    :type output: array
    :param m: number of Markov parameters to compute. if None, uses the minimum of 300 and the number of timesteps in `input` and `output`.
    :type m: int

    :return: the Markov parameters, with dimensions :math:`(p,q,m+1)`
    :rtype: array

