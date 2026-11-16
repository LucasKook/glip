"""Background knowledge with the semantics of the classic glip encodings.

Input (already validated and normalised by glip's .bk_validate() in R; 0/1
d x d matrices with zero diagonal, symmetric for *_bi and *_undirected):
forbidden, required (directed edges i -> j), forbidden_bi, required_bi
(ADMG), forbidden_undirected, required_undirected (chain graphs),
forbidden_anc, required_anc (directed paths; anterior paths for chain
graphs), edge_constraints (list of {coef, coef_bi, sense, rhs}) and
edge_penalty. In terms of glip's directed indicator x_ij (for chain graphs
x_ij = 1 iff i -> j or i - j, see classes.*.xexpr):
- DAG/ADMG: forbidden fixes x_ij = 0, required x_ij = 1; *_bi fix the
  bidirected indicator.
- Chain graphs (glip's .bk_chain()): required[i, j] -> x_ij = 1, x_ji = 0;
  required_undirected -> x_ij = x_ji = 1; forbidden[i, j] -> x_ij - x_ji <= 0;
  forbidden_undirected -> x_ij + x_ji <= 1.
- Ancestry: forbidden_anc[i, j]: no (directed / anterior) path from i to j,
  by reachability rows; required_anc[i, j]: a path, by a unit flow.
- edge_constraints: sum coef[i, j] x_ij + sum_{i<j} coef_bi[i, j] b_ij
  (sense) rhs.
- edge_penalty times every edge indicator (x_ij and b_ij), as in the classic
  encodings (so an undirected chain edge counts twice).
"""
from collections import deque

import gurobipy as gp
from gurobipy import GRB

MATRICES = ('forbidden', 'required', 'forbidden_bi', 'required_bi',
            'forbidden_undirected', 'required_undirected', 'forbidden_anc',
            'required_anc')


def normalise(K, G, d):
    """Drop empty entries; error for knowledge the class does not support."""
    K = dict(K or {})
    out = {}
    for k, v in K.items():
        if v is None or (k == 'edge_constraints' and not v) or (k == 'edge_penalty' and not v):
            continue
        if k in MATRICES:
            if k not in G.knowledge_keys:
                raise ValueError(f'`{k}` is not supported for class {G.name!r}')
            assert len(v) == d and all(len(r) == d for r in v), k
            M = [[int(v[a][b]) if a != b else 0 for b in range(d)] for a in range(d)]
            if any(M[a][b] for a in range(d) for b in range(d)):
                out[k] = M
        elif k == 'edge_constraints':
            cons = []
            for c in v:
                c = {f: c[f] for f in c if c[f] is not None}
                bad = set(c) - set(G.coef_keys) - {'sense', 'rhs'}
                if bad:
                    raise ValueError(f'edge constraint fields {sorted(bad)} not supported for class {G.name!r}')
                if c.get('sense') not in ('<=', '>=', '='):
                    raise ValueError('edge constraint sense must be <=, >= or =')
                cons.append(c)
            if cons:
                out[k] = cons
        elif k == 'edge_penalty':
            out[k] = float(v)
            if out[k] < 0:
                raise ValueError('edge_penalty must be non-negative')
        else:
            raise ValueError(f'unknown background knowledge `{k}`')
    return out


def _x(G, i, j, X):
    return gp.quicksum(X[e] for e in G.xexpr(i, j))


def _linear_terms(G, d, c):
    """(edge, coefficient) pairs of an edge constraint's left side."""
    terms = []
    cf = c.get('coef')
    if cf is not None:
        for a in range(d):
            for b in range(d):
                if a != b and cf[a][b] != 0:
                    terms += [(e, float(cf[a][b])) for e in G.xexpr(a, b)]
    cb = c.get('coef_bi')
    if cb is not None:
        for a in range(d):
            for b in range(a + 1, d):
                if cb[a][b] != 0:
                    terms.append((G.bi_var(a, b), float(cb[a][b])))
    return terms


def _fixings(G, d, K):
    """Linear rows (terms, sense, rhs) for edge fixings and forbidden edges."""
    rows = []
    F, R = K.get('forbidden'), K.get('required')
    pairs = [(a, b) for a in range(d) for b in range(d) if a != b]
    if G.name == 'chain':
        RU, FU = K.get('required_undirected'), K.get('forbidden_undirected')
        for a, b in pairs:
            if (R and R[a][b]) or (RU and RU[a][b]):
                rows.append(([(e, 1.0) for e in G.xexpr(a, b)], '>=', 1.0))
            if R and R[b][a]:
                rows.append(([(e, 1.0) for e in G.xexpr(a, b)], '<=', 0.0))
            if F and F[a][b]:
                rows.append(([(e, 1.0) for e in G.xexpr(a, b)] +
                             [(e, -1.0) for e in G.xexpr(b, a)], '<=', 0.0))
            if FU and FU[a][b] and a < b:
                rows.append(([(e, 1.0) for e in G.xexpr(a, b)] +
                             [(e, 1.0) for e in G.xexpr(b, a)], '<=', 1.0))
    else:
        for a, b in pairs:
            if F and F[a][b]:
                rows.append(([(e, 1.0) for e in G.xexpr(a, b)], '<=', 0.0))
            if R and R[a][b]:
                rows.append(([(e, 1.0) for e in G.xexpr(a, b)], '>=', 1.0))
        FB, RB = K.get('forbidden_bi'), K.get('required_bi')
        for a, b in pairs:
            if a < b and FB and FB[a][b]:
                rows.append(([(G.bi_var(a, b), 1.0)], '<=', 0.0))
            if a < b and RB and RB[a][b]:
                rows.append(([(G.bi_var(a, b), 1.0)], '>=', 1.0))
    for c in K.get('edge_constraints', []):
        rows.append((_linear_terms(G, d, c), c['sense'], float(c['rhs'])))
    return rows


def apply(m, X, G, d, K):
    """Add the knowledge to model m; return the penalty expression."""
    sense = {'<=': GRB.LESS_EQUAL, '>=': GRB.GREATER_EQUAL, '=': GRB.EQUAL}
    for terms, s, rhs in _fixings(G, d, K):
        if len(terms) == 1 and terms[0][1] == 1.0 and s != '=':
            x = X[terms[0][0]]  # single indicator: bound
            if s == '<=' and rhs == 0.0:
                x.UB = 0
                continue
            if s == '>=' and rhs == 1.0:
                x.LB = 1
                continue
        m.addLConstr(gp.quicksum(c * X[e] for e, c in terms), sense[s], rhs)
    arcs = G.anc_arcs(d)
    Fa = K.get('forbidden_anc')
    if Fa:
        for i in range(d):
            if not any(Fa[i]):
                continue
            a = [m.addVar(lb=1 if k == i else 0, ub=0 if Fa[i][k] else 1) for k in range(d)]
            for k, l, e in arcs:
                m.addConstr(a[l] >= a[k] + X[e] - 1)
    Ra = K.get('required_anc')
    if Ra:
        for i in range(d):
            for j in range(d):
                if i == j or not Ra[i][j]:
                    continue
                f = [m.addVar(lb=0, ub=1) for _ in arcs]
                for fv, (k, l, e) in zip(f, arcs):
                    m.addConstr(fv <= X[e])
                for v in range(d):
                    out = gp.quicksum(fv for fv, (k, l, e) in zip(f, arcs) if k == v)
                    inn = gp.quicksum(fv for fv, (k, l, e) in zip(f, arcs) if l == v)
                    m.addConstr(out - inn == (1 if v == i else -1 if v == j else 0))
    pen = K.get('edge_penalty', 0.0)
    if pen:
        return gp.quicksum(pen * c * X[e] for e, c in G.penalty_terms(d))
    return gp.LinExpr()


def _reach(G, d, val, i):
    out = {}
    for k, l, e in G.anc_arcs(d):
        if val[e] > .5:
            out.setdefault(k, []).append(l)
    seen, q = {i}, deque([i])
    while q:
        for l in out.get(q.popleft(), []):
            if l not in seen:
                seen.add(l)
                q.append(l)
    return seen


def satisfied(G, d, val, K, tol=1e-6):
    """Independent check of a 0/1 graph (edge values) against the knowledge."""
    ok = {'<=': lambda a, b: a <= b + tol, '>=': lambda a, b: a >= b - tol,
          '=': lambda a, b: abs(a - b) <= tol}
    for terms, s, rhs in _fixings(G, d, K):
        if not ok[s](sum(c * round(val[e]) for e, c in terms), rhs):
            return False
    Fa, Ra = K.get('forbidden_anc'), K.get('required_anc')
    for i in range(d):
        if (Fa and any(Fa[i])) or (Ra and any(Ra[i])):
            R = _reach(G, d, val, i)
            for j in range(d):
                if j != i and ((Fa and Fa[i][j] and j in R) or (Ra and Ra[i][j] and j not in R)):
                    return False
    return True


def penalty_value(G, d, val, K):
    pen = K.get('edge_penalty', 0.0)
    return pen * sum(c * round(val[e]) for e, c in G.penalty_terms(d)) if pen else 0.0
