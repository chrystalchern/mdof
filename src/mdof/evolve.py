from . import numerics

def obsv2ac(Observability, no, p, **options):
    """
    Compute the system matrices ``A`` and ``C`` from the ``Observability`` matrix.
    """
    lsq_solve = numerics.lsq_solver(options.get("lsq", {}))
    A = lsq_solve(Observability[:(no-1)*p,:], Observability[p:no*p,:])
    C = Observability[:p,:]
    return (A,C)