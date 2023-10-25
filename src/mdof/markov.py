import numpy as np
import warnings

# inputs = input data. dimensions of inputs: q x nt, where nt = number of timesteps.
# outputs = response data due to input data. dimensions of outputs: p x nt.
# m = number of Markov parameters (impulse response timesteps), not including timestep zero, to solve for.
def okid(inputs,outputs,**options):
    """
    Identify Markov parameters, or discrete impulse response data, for a given set of input and output data.
    Observer Kalman Identification Algorithm (OKID) [1]_.

    :param inputs:  input time history. dimensions: :math:`(q,nt)`, where
                    :math:`q` = number of inputs, and :math:`nt` = number of timesteps
    :type inputs:   array
    :param outputs: output response history.
                    dimensions: :math:`(p,nt)`, where :math:`p` = number of outputs, and
                    :math:`nt` = number of timesteps
    :type outputs:  array
    :param m:       number of Markov parameters to compute. default: :math:`min(300, nt)`
    :type m:        int, optional

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
    Ybar = outputs @ np.linalg.pinv(V,rcond=10e-3)  # DIM: (p)x(q+m(q+p))
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