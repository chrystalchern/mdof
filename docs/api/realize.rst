``realize`` Module
===================

.. py:function:: ssid.realize.era(Y,**options)

    Generate a state space system realization from set of Markov parameters, or discrete impulse response data, using the Eigensystem Realization Algorithm (ERA) (Juang and Pappa, 1985).

    :param Y: Markov parameters. dimensions: :math:`(p,q,nt)`, where :math:`p` = number of output channels, :math:`q` = number of input channels, and :math:`nt` = number of timesteps
    :type Y: array
    :param no: number of block rows in Hankel matrix = order of observability matrix
    :type no: int
    :param nc: number of block columns in Hankel matrix = order of controllability matrix
    :type nc: int
    :param r: reduced model order. default: minimum of 10 and ``no``/2
    :type r: int, optional

    :return: realization in the form of state space coefficients ``(A,B,C,D)``
    :rtype: tuple

.. py:function:: ssid.realize.era_dc(Y,**options)

    Generate a state space system realization from set of Markov parameters, or discrete impulse response data, using the Eigensystem Realization Algorithm with Data Correlations (ERA/DC) (Juang, Cooper, and Wright, 1988).

    :param Y: Markov parameters. dimensions: :math:`(p,q,nt)`, where :math:`p` = number of output channels, :math:`q` = number of input channels, and :math:`nt` = number of timesteps
    :type Y: array
    :param no: number of block rows in Hankel matrix = order of observability matrix
    :type no: int
    :param nc: number of block columns in Hankel matrix = order of controllability matrix
    :type nc: int
    :param a: (alpha) number of block rows in Hankel of correlation matrix. default: 0
    :type a: int, optional
    :param b: (beta) number of block columns in Hankel of correlation matrix. default: 0
    :type b: int, optional
    :param l: initial lag for data correlations. default: 0
    :type l: int, optional
    :param g: lags (gap) between correlation matrices. default: 1
    :type g: int, optional
    :param r: reduced model order. default: minimum of 10 and ``no``/2
    :type r: int, optional

    :return: realization in the form of state space coefficients ``(A,B,C,D)``
    :rtype: tuple

.. py:function:: ssid.realize.srim(input,output,**options)

    Generate a state space system realization from set of Markov parameters, or discrete impulse response data, for a given set of input and output data, using System Realization Using Information Matrix (SRIM) (Juang, 1996).
    
    :param input: input acceleration time history. dimensions: :math:`(q,nt)`, where :math:`q` = number of input channels, and :math:`nt` = number of timesteps
    :type input: array
    :param output: output acceleration response history to an arbitrary input. dimensions: :math:`(p,nt)`, where :math:`p` = number of output channels, and :math:`nt` = number of timesteps
    :type output: array
    :param no: number of steps used for identification (prediction horizon). default: minimum of 300 and number of timesteps in input and output data.
    :type no: int
    :param r: reduced model order. default: minimum of 10 and ``no``/2
    :type r: int
    :param full: if True, full SVD. default: True
    :type r: bool

    :return: realization in the form of state space coefficients ``(A,B,C,D)``
    :rtype: tuple