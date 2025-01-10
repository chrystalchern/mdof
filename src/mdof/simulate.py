import numpy as np
import scipy
import warnings


def simulate(system, inputs, w=None, backend=None):
    """
    Simulate the response for a discrete-time system from its system
    matrices ``(A,B,C,D)`` and an `inputs` array.
    """
    (A, B, C, D) = system
    if backend is None:
        return _simulate(system=(A, B, C, D), inputs=inputs, w=w)[1]
    elif backend=='scipy':
        return scipy.signal.dlsim(scipy.signal.dlti(A, B, C, D), inputs.T)[1].T
    elif backend=='control':
        return _simulate_control(realization=system, dt=1, inputs=inputs)
    else:
        warnings.warn("backend not recognized.")
        return
    

# in-house simulate
def _simulate(system, inputs, w=None):
    (A,B,C,D) = system
    p,n = C.shape
    if inputs.ndim == 1:
        inputs = inputs[None,:]
    elif inputs.ndim != 2:
        raise ValueError("inputs must be a 1D or 2D array.")
    
    if inputs.shape[0] > inputs.shape[1]:
        warnings.warn("There are more inputs than timesteps. Recommend double-checking the shape of the input array.")
        
    _,T = inputs.shape

    x = np.empty((n,T), dtype=complex)
    x[:,0] = np.zeros((n,))
    for t in range(T-1):
        x[:,t+1] = A@(x[:,t]) + B@(inputs[:,t])

    if w is None:
        w = inputs
    y = np.real(C@x + D@w)
    assert y.shape == (p,T)

    return x, y


from scipy.interpolate import make_interp_spline
from scipy.signal import lti, dlti, StateSpace
# scipy signal dlsim
def dlsim(system, u, t=None, x0=None):
    """
    Simulate output of a discrete-time linear system.

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

    u : array_like
        An input array describing the input at each time `t` (interpolation is
        assumed between given times).  If there are multiple inputs, then each
        column of the rank-2 array represents an input.
    t : array_like, optional
        The time steps at which the input is defined.  If `t` is given, it
        must be the same length as `u`, and the final value in `t` determines
        the number of steps returned in the output.
    x0 : array_like, optional
        The initial conditions on the state vector (zero by default).

    Returns
    -------
    tout : ndarray
        Time values for the output, as a 1-D array.
    yout : ndarray
        System response, as a 1-D array.
    xout : ndarray, optional
        Time-evolution of the state-vector.  Only generated if the input is a
        `StateSpace` system.

    See Also
    --------
    lsim, dstep, dimpulse, cont2discrete

    Examples
    --------
    A simple integrator transfer function with a discrete time step of 1.0
    could be implemented as:

    >>> import numpy as np
    >>> from scipy import signal
    >>> tf = ([1.0,], [1.0, -1.0], 1.0)
    >>> t_in = [0.0, 1.0, 2.0, 3.0]
    >>> u = np.asarray([0.0, 0.0, 1.0, 1.0])
    >>> t_out, y = signal.dlsim(tf, u, t=t_in)
    >>> y.T
    array([[ 0.,  0.,  0.,  1.]])

    """
    # We are looking for A, B, C, D, dt
    # ---------------------------------------------------------------------------

    # Convert system to dlti-StateSpace
    if isinstance(system, lti):
        raise AttributeError('dlsim can only be used with discrete-time dlti '
                             'systems.')
    elif not isinstance(system, dlti):
        system = dlti(*system[:-1], dt=system[-1])

    # Condition needed to ensure output remains compatible
    is_ss_input = isinstance(system, StateSpace)
    system = system._as_ss()

    A = system.A 
    B = system.B
    C = system.C 
    D = system.D
    dt = system.dt

    # We have A, B, C, D, dt
    # ---------------------------------------------------------------------------

    u = np.atleast_1d(u)

    if u.ndim == 1:
        u = np.atleast_2d(u).T

    if t is None:
        out_samples = len(u)
        stoptime = (out_samples - 1) * dt
    else:
        stoptime = t[-1]
        out_samples = int(np.floor(stoptime / dt)) + 1

    # Pre-build output arrays
    xout = np.zeros((out_samples, A.shape[0]))
    yout = np.zeros((out_samples, C.shape[0]))
    tout = np.linspace(0.0, stoptime, num=out_samples)

    # Check initial condition
    if x0 is None:
        xout[0, :] = np.zeros((A.shape[1],))
    else:
        xout[0, :] = np.asarray(x0)

    # Pre-interpolate inputs into the desired time steps
    if t is None:
        u_dt = u
    else:
        if len(u.shape) == 1:
            u = u[:, np.newaxis]

        u_dt = make_interp_spline(t, u, k=1)(tout)

    # Simulate the system
    for i in range(0, out_samples - 1):
        xout[i+1,:] = A@xout[i,:] + B@u_dt[i,:]
        yout[i,:]   = C@xout[i,:] + D@u_dt[i,:]

    # Last point
    yout[out_samples-1, :] = C@xout[out_samples-1,:] + D@u_dt[out_samples-1,:]

    if is_ss_input:
        return tout, yout, xout
    else:
        return tout, yout


from control import ss, forced_response
# control forced_response
def _simulate_control(realization, dt, inputs):
    out_pred = forced_response(ss(*realization,dt), U=inputs, squeeze=False, return_x=False).outputs
    return out_pred
