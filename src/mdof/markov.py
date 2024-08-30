import numpy as np
import warnings

def okid(inputs, outputs,**options):
    r"""
    Identify Markov parameters, or discrete impulse response data, for a given set of input and output data.
    Observer Kalman Identification Algorithm (OKID) [1]_.

    :param inputs:  input time history. dimensions: :math:`(q,nt)`, where
                    :math:`q` = number of inputs, and :math:`nt` = number of timesteps
    :type inputs:   array
    :param outputs: output response history.
                    dimensions: :math:`(p,nt)`, where :math:`p` = number of outputs, and
                    :math:`nt` = number of timesteps
    :type outputs:  array
    :param m:       (optional) number of Markov parameters to compute, excluding timestep zero. default: :math:`\min(300, nt)`
    :type m:        int
    :param rcond:   (optional) cut-off ratio for small singular values in pseudoinverse computation. default: `1e-15`
    :type m:        float

    :return: the Markov parameters, with dimensions :math:`(p,q,m+1)`
    :rtype: array

    References
    ----------
    .. [1]  Juang, J. N., Phan, M., Horta, L. G., & Longman, R. W. (1993). Identification of
            observer/Kalman filter Markov parameters-Theory and experiments. Journal of Guidance,
            Control, and Dynamics, 16(2), 320-329. (https://doi.org/10.2514/3.21006)
    """

    if len(inputs.shape) == 1:
        inputs = inputs[None,:]
    if len(outputs.shape) == 1:
        outputs = outputs[None,:]
    
    if inputs.shape[0] > inputs.shape[1]:
        warnings.warn("input data has more channels (dim 1) than timesteps (dim 2)")
    if outputs.shape[0] > outputs.shape[1]:
        warnings.warn("output data has more channels (dim 1) than timesteps (dim 2)")

    q,nt = inputs.shape
    p = outputs.shape[0]
    assert nt == outputs.shape[1]

    m = options.get("m", None)
    if m is None:
        m = min(1000,nt)

    # Form data matrix V
    V = np.zeros((m+1,nt,q+p))                    # m+1 block rows, each with nt columns and q+p rows each.  DIM: (m+1)x(nt)x(q+p)
    for i in range(m+1):                          # From block row 0 to block row m
        V[i,-(nt-i):,:q] = inputs[:,:nt-i].T      # Populate the first q rows of the ith block row, last nt-i columns, with the first nt-i inputs.  DIM: (nt-i)x(q)
        V[i,-(nt-i):,-p:] = outputs[:,:nt-i].T    # Populate the last p rows of the ith block row, last nt-i columns, with the first nt-i outputs.  DIM: (nt-i)x(p)
    V = V.transpose((1,0,2)).reshape((nt,(m+1)*(p+q))).T    # Transpose and reshape data matrix to stack the block rows in.  DIM: ((m+1)(q+p))x(nt)
    # V = np.concatenate((V[:q,:],V[q+p:,:]))     # Remove rows q+1 to q+p (TODO: see if it's necessary to remove rows, or if solving for Ybar can include these)
    V = np.delete(V, slice(q,q+p), axis=0)        # Remove rows q+1 to q+p (TODO: see if it's necessary to remove rows, or if solving for Ybar can include these)

    # Solve for observer Markov parameters Ybar
    Ybar = outputs @ np.linalg.pinv(V,rcond=options.get('rcond',1e-15))  # DIM: (p)x(q+m(q+p))
    assert Ybar.shape == (p,q+m*(q+p))
    
    # Isolate system Markov parameters
    H = np.zeros((p,q,m+1))
    D = Ybar[:,:q] # feed-through term (or D matrix) is the first block column
    H[:,:,0] = D

    Y = np.zeros((p,q,m))
    Ybar1 = np.zeros((p,q,m))
    Ybar2 = np.zeros((p,p,m))
    
    for i in range(m):
        Ybar1[:,:,i] = Ybar[:,q+(q+p)*i : q+(q+p)*i+q]
        Ybar2[:,:,i] = Ybar[:,q+(q+p)*i+q : q+(q+p)*(i+1)]
    
    Y[:,:,0] = Ybar1[:,:,0] + Ybar2[:,:,0] @ D
    for k in range(1,m):
        Y[:,:,k] = Ybar1[:,:,k] + Ybar2[:,:,k] @ D
        for i in range(k-1):
            Y[:,:,k] += Ybar2[:,:,i] @ Y[:,:,k-i-1]

    H[:,:,1:] = Y

    return H


from scipy.signal import lti, dlti, dlsim
def _dimpulse(system, x0=None, t=None, n=None):
    # scipy signal
    """
    Impulse response of discrete-time system.

    Parameters
    ----------
    system : tuple of array_like or instance of `dlti`
        A tuple describing the system.
        The following gives the number of elements in the tuple and
        the interpretation:

            * 1: (instance of `dlti`)
            * 3: (num, den, dt)
            * 4: (zeros, poles, gain, dt)
            * 5: (A, B, C, D, dt)

    x0 : array_like, optional
        Initial state-vector.  Defaults to zero.
    t : array_like, optional
        Time points.  Computed if not given.
    n : int, optional
        The number of time points to compute (if `t` is not given).

    Returns
    -------
    tout : ndarray
        Time values for the output, as a 1-D array.
    yout : tuple of ndarray
        Impulse response of system.  Each element of the tuple represents
        the output of the system based on an impulse in each input.

    See Also
    --------
    impulse, dstep, dlsim, cont2discrete

    Examples
    --------
    >>> import numpy as np
    >>> from scipy import signal
    >>> import matplotlib.pyplot as plt

    >>> butter = signal.dlti(*signal.butter(3, 0.5))
    >>> t, y = signal.dimpulse(butter, n=25)
    >>> plt.step(t, np.squeeze(y))
    >>> plt.grid()
    >>> plt.xlabel('n [samples]')
    >>> plt.ylabel('Amplitude')

    """
    # Convert system to dlti-StateSpace
    if isinstance(system, dlti):
        system = system._as_ss()
    elif isinstance(system, lti):
        raise AttributeError('dimpulse can only be used with discrete-time '
                             'dlti systems.')
    else:
        system = dlti(*system[:-1], dt=system[-1])._as_ss()

    # Default to 100 samples if unspecified
    if n is None:
        n = 100

    # If time is not specified, use the number of samples
    # and system dt
    if t is None:
        t = np.linspace(0, n * system.dt, n, endpoint=False)
    else:
        t = np.asarray(t)

    # For each input, implement a step change
    yout = None
    for i in range(0, system.inputs):
        u = np.zeros((t.shape[0], system.inputs))
        u[0, i] = 1.0

        one_output = dlsim(system, u, t=t, x0=x0)

        if yout is None:
            yout = (one_output[1],)
        else:
            yout = yout + (one_output[1],)

        tout = one_output[0]

    return tout, yout