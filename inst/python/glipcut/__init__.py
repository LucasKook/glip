"""glipcut: cut + reachability MIP formulation for glip (DAG, ADMG, LWF chain graphs).

Called from R (glip::cuts_optim(), glip::learn_graph(..., encoding = "cuts"))
as a subprocess: python -m glipcut <in.json> <out.json>. See solver.py for
the formulation; gurobipy is imported only when solving.
"""
__version__ = '0.1.0'


def solve(case, **overrides):
    from .solver import solve as _solve
    return _solve(case, **overrides)
