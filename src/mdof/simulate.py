def simulate((A,B,C,D), (u, w)):
    p,n = C.shape
    T, _ = u.shape
    # Initialize
    x = np.empty((T,n))

    x[0] = 0
    for t in xrange(T-1):
        x[t+1] = A.dot(x[t]) + B.dot(u[t])

    y = C.dot(x.T).T + D.dot(w.T).T

    return x, y
