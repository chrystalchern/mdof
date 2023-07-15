
def lsq_solver(options):
    import numpy.linalg
    return lambda *args: numpy.linalg.lstsq(*args, rcond=None)[0]

def lin_solver():
    pass

def svd_routine():
    import scipy
    def _svd(*args):
        U,S,V = scipy.linalg.svd(*args, lapack_driver="gesvd")
        return U,S,V.T.conj()
    return _svd

def form_observability():
    pass
