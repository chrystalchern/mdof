``realize`` Module
===================

.. py:function:: ssid.realize.era(Y,no=None,nc=None,r=None)

    Generate a state space system realization from set of Markov parameters, or discrete impulse response data, using the Eigensystem Realization Algorithm (ERA) (Juang and Pappa, 1985).

    :param Y: Markov parameters. dimensions: :math:`(p,q*nt)`, where :math:`p` = number of output channels, :math:`q` = number of input channels, and :math:`nt` = number of timesteps
    :type Y: array
    :param no: number of block rows in Hankel matrix = order of controllability matrix
    :type no: int
    :param nc: number of block columns in Hankel matrix = order of observability matrix
    :type nc: int
    :param r: reduced model order; default None
    :type r: int, optional

    :return: realization in the form of state space coefficients ``A,B,C,D``
    :rtype: list[array]

.. py:function:: ssid.realize.era_dc(Y,no=None,nc=None,a=0,b=0,l=0,g=1,r=None)

    Generate a state space system realization from set of Markov parameters, or discrete impulse response data, using the Eigensystem Realization Algorithm with Data Correlations (ERA/DC) (Juang, Cooper, and Wright, 1988).

    :param Y: Markov parameters. dimensions: :math:`(p,q*nt)`, where :math:`p` = number of output channels, :math:`q` = number of input channels, and :math:`nt` = number of timesteps
    :type Y: array
    :param no: number of block rows in Hankel matrix = order of controllability matrix
    :type no: int
    :param nc: number of block columns in Hankel matrix = order of observability matrix
    :type nc: int
    :param a: (alpha) number of block rows in Hankel of correlation matrix
    :type a: int, optional
    :param b: (beta) number of block columns in Hankel of correlation matrix
    :type b: int, optional
    :param l: initial lag for data correlations
    :type l: int, optional
    :param g: lags (gap) between correlation matrices
    :type g: int, optional
    :param r: reduced model order
    :type r: int, optional

    :return: realization in the form of state space coefficients ``A,B,C,D``
    :rtype: list[array]

.. py:function:: ssid.realize.srim(input, output, mro=100, orm=20)

    Generate a state space system realization from set of Markov parameters, or discrete impulse response data, for a given set of input and output data, using System Realization Using Information Matrix (SRIM) (Juang, 1996).
    
    :param input: input acceleration time history. dimensions: :math:`(q,nt)`, where :math:`q` = number of input channels, and :math:`nt` = number of timesteps
    :type input: array
    :param output: output acceleration response history to an arbitrary input. dimensions: :math:`(p,nt)`, where :math:`p` = number of output channels, and :math:`nt` = number of timesteps
    :type output: array
    :param mro: number of steps used for identification (prediction horizon)
    :type mro: int
    :param rmo: reduced model order
    :type rmo: int, optional

    :return: realization in the form of state space coefficients ``A,B,C,D``
    :rtype: list[array]