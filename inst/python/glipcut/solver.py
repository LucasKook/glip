"""Cut + reachability MIP for exact GLIP learning.

Same problem as the classic glip encodings: choose a graph of the requested
class minimising sum_q w_q |connected_q(G) - p_q| over the CI decisions q =
(i, j, C) (p = 1: dependence). Only the MIP formulation differs.

Formulation (ported from the mu-GLIP solver `glipcut` in cometsTS and adapted
to m-separation, where a query connects if the target is reached with either
arrival mark):
- separated decisions (p = 0): closure rows per (source, C) over the states
  of the class automaton; z >= r_(j,h) for both marks h;
- dependent decisions (p = 1): y <= sum of the edges leaving a state set
  containing the source but no target state; static cut from the source,
  exact lazy cuts at integer solutions, root min cuts;
- joint cuts (option `joint`): y <= r_(j,0) + r_(j,1) when (source, C) has a
  closure block, and root min cuts with arc capacities min(x_e, r_head);
- subset consistency cuts (option `subset`): for every 3-node subset S, the
  queries inside S are the separations of the marginal (latent projection)
  of G on S, which is a 3-node model of the marginal class (see `local.py`):
  a static row sum w mismatch >= local minimum and root separation over the
  local polytope ("lp", lambda cuts);
- `orient`: m-connection is symmetric, so each query may use either endpoint
  as the source; "cover" chooses sources so that few closure blocks are
  needed and dependent queries share a closure block where possible;
- `prune`: closure states that reach no separated target are dropped;
- background knowledge and edge penalty with the semantics of the classic
  encodings (knowledge.py).
Acyclicity (directed edges; semidirected cycles for chain graphs) by static
2-cycle rows and lazy cycle cuts.

At an integer point z >= true connection (r is an upper bound on
reachability) and y <= true connection (lazy cuts), so every mismatch proxy
is >= the true mismatch and all subset/joint cuts are valid. At a time
limited incumbent the model objective may exceed the oracle loss; the oracle
value is reported (raw value kept as objval_model).
"""
import itertools
import math
import time

import gurobipy as gp
from gurobipy import GRB
import networkx as nx

from . import classes, knowledge
from .local import LOCAL_QUERIES, local_models

DEFAULTS = dict(time_limit=60.0, threads=1, seed=1, mip_gap=0.0, mip_gap_abs=0.0,
                nodecuts='root', subset='lp', joint=True, prune=False,
                orient='none', subset_rounds=50)
_STATUS = {getattr(GRB.Status, k): k for k in dir(GRB.Status) if k.isupper()}


def _orient(queries, how):
    """Return queries with chosen source endpoint (connection is symmetric)."""
    qs = [dict(q) for q in queries]
    if how == 'none':
        return qs
    if how != 'cover':
        raise ValueError(how)
    # greedy vertex cover of the separated pairs, per conditioning set
    sources = {}
    byC = {}
    for q in qs:
        if q['p'] == 0:
            byC.setdefault(tuple(sorted(q['C'])), []).append((q['i'], q['j']))
    for C, pairs in byC.items():
        left, chosen = set(map(frozenset, pairs)), set()
        while left:
            deg = {}
            for e in left:
                for v in e:
                    deg[v] = deg.get(v, 0) + 1
            v = max(sorted(deg), key=lambda u: deg[u])
            chosen.add(v)
            left = {e for e in left if v not in e}
        sources[C] = chosen
    for q in qs:
        S = sources.get(tuple(sorted(q['C'])), set())
        if q['i'] not in S and q['j'] in S:
            q['i'], q['j'] = q['j'], q['i']
    return qs


def solve(case, **overrides):
    wall0 = time.perf_counter()
    opts = dict(DEFAULTS)
    opts.update(case.get('options', {}) or {})
    opts.update(overrides)
    gclass = case.get('class', 'admg')
    G = classes.get(gclass)
    d = int(case['d'])
    K = knowledge.normalise(case.get('knowledge'), G, d)
    assert all(q['i'] != q['j'] and not {q['i'], q['j']} & set(q['C'])
               and q['p'] in (0, 1) and q['w'] >= 0 for q in case['queries'])
    queries = _orient([q for q in case['queries'] if q['w'] > 0], opts['orient'])

    env = gp.Env(empty=True)
    env.setParam('OutputFlag', 0)
    env.start()
    m = gp.Model('cuts', env=env)
    m.Params.Threads = int(opts['threads'])
    m.Params.Seed = int(opts['seed'])
    m.Params.TimeLimit = float(opts['time_limit'])
    m.Params.MIPGap = float(opts['mip_gap'])
    m.Params.MIPGapAbs = float(opts['mip_gap_abs'])
    for k, v in (opts.get('gurobi_params') or {}).items():
        m.setParam(k, v)
    m.Params.LazyConstraints = 1
    m.Params.PreCrush = 1

    X = {e: m.addVar(vtype=GRB.BINARY) for e in G.edges(d)}
    for x in X.values():
        x.BranchPriority = 1
    for row in G.static_rows(d):
        m.addConstr(gp.quicksum(X[e] for e in row) <= 1)
    penalty = knowledge.apply(m, X, G, d, K)

    start = case.get('start')
    sv = None
    warnings = []
    if start:
        sv = G.start_values(d, start)
        if G.violated_cycle(d, sv):
            warnings.append('The warm start is not a valid graph of the class and was ignored.')
            sv = None
        elif not knowledge.satisfied(G, d, sv, K):
            warnings.append('The warm start violates the background knowledge and was ignored.')
            sv = None
        else:
            for e, x in X.items():
                x.Start = sv[e]

    arcs_cache = {}

    def arcs_of(i, C):
        if (i, C) not in arcs_cache:
            arcs_cache[i, C] = G.state_arcs(d, set(C), i)
        return arcs_cache[i, C]

    stats = dict(lazy=0, user=0, cycle=0, closure_rows=0, joint_rows=0,
                 subset_static=0, subset_lp=0, subset_rounds=0)
    objective = gp.LinExpr() + penalty
    proxy = {}  # unordered query key -> list of (kind, var, w)

    # Separated decisions: closure blocks
    targets = {}
    for q in queries:
        if q['p'] == 0:
            key = q['i'], tuple(sorted(q['C']))
            targets.setdefault(key, set()).update(G.target_states(q['j']))
    R = {}
    for key, tg in targets.items():
        arcs = arcs_of(*key)
        states = {st for s, t, _ in arcs for st in (s, t) if st != 'S'}
        if opts['prune']:
            back = {}
            for s, t, _ in arcs:
                back.setdefault(t, []).append(s)
            keep, stack = set(tg) & states, list(set(tg) & states)
            while stack:
                for s in back.get(stack.pop(), []):
                    if s != 'S' and s not in keep:
                        keep.add(s)
                        stack.append(s)
            states = keep
        r = {st: m.addVar(lb=0, ub=1) for st in states}
        if sv is not None:
            seen = classes.reachable(arcs, lambda e: sv[e] > .5)
            for st, x in r.items():
                x.Start = int(st in seen)
        for s, t, e in arcs:
            if t not in r or (s != 'S' and s not in r):
                continue
            m.addConstr(r[t] >= X[e] if s == 'S' else r[t] >= r[s] + X[e] - 1)
            stats['closure_rows'] += 1
        R[key] = r
    dep = []
    for q in queries:
        key = q['i'], tuple(sorted(q['C']))
        conn0 = None if sv is None else G.connected_vals(d, sv, q['i'], q['j'], q['C'])
        if q['p'] == 0:
            z = m.addVar(lb=0, ub=1)
            for st in G.target_states(q['j']):
                if st in R[key]:
                    m.addConstr(z >= R[key][st])
            if conn0 is not None:
                z.Start = int(conn0)
            objective += q['w'] * z
            proxy.setdefault(classes.qkey(q), []).append(('z', z, q['w']))
        else:
            y = m.addVar(lb=0, ub=1)
            if conn0 is not None:
                y.Start = int(conn0)
            objective += q['w'] * (1 - y)
            arcs = arcs_of(*key)
            m.addConstr(y <= gp.quicksum(X[e] for e in {e for s, t, e in arcs if s == 'S'}))
            tset = set(G.target_states(q['j']))
            if opts['joint'] and key in R:
                have = [R[key][st] for st in tset if st in R[key]]
                if len(have) == len(tset):
                    m.addConstr(y <= gp.quicksum(have))
                    stats['joint_rows'] += 1
            dep.append((y, arcs, tset, key))
            proxy.setdefault(classes.qkey(q), []).append(('y', y, q['w']))
    m.setObjective(objective)

    # Subset consistency cuts over 3-node subsets
    subsets = []
    if opts['subset'] != 'none' and d >= 3:
        models = local_models(gclass)
        for S in itertools.combinations(range(d), 3):
            items = []
            for lq, (a, b, C) in enumerate(LOCAL_QUERIES):
                k = classes.qkey(dict(i=S[a], j=S[b], C=[S[c] for c in C]))
                for kind, var, w in proxy.get(k, []):
                    items.append((lq, kind, var, w, 0 if kind == 'z' else 1))
            if not items:
                continue
            rows = sorted({tuple(int(M[lq] != p) for lq, _, _, _, p in items) for M in models})
            exprs = [var if kind == 'z' else 1 - var for _, kind, var, _, _ in items]
            w = [it[3] for it in items]
            mS = min(sum(wi * ri for wi, ri in zip(w, row)) for row in rows)
            if mS > 1e-9:
                m.addConstr(gp.quicksum(wi * e for wi, e in zip(w, exprs)) >= mS)
                stats['subset_static'] += 1
            lp = None
            if opts['subset'] == 'lp':
                lpm = gp.Model(env=env)
                lpm.Params.Threads = 1
                lam = lpm.addVars(len(items), lb=0, ub=1)
                t = lpm.addVar(lb=-GRB.INFINITY)
                for row in rows:
                    lpm.addConstr(t <= gp.quicksum(lam[k] * row[k] for k in range(len(items))))
                lpm.ModelSense = GRB.MAXIMIZE
                lpm.update()
                lp = (lpm, lam, t)
            subsets.append((items, exprs, lp))

    rlist = [((key, st), v) for key, rb in R.items() for st, v in rb.items()] if opts['joint'] else []
    errors = []

    def subset_separation(model):
        vals = model.cbGetNodeRel([it[2] for items, _, _ in subsets for it in items])
        pos = 0
        for items, exprs, (lpm, lam, t) in subsets:
            v = []
            for it in items:
                v.append(vals[pos] if it[1] == 'z' else 1 - vals[pos])
                pos += 1
            lpm.setObjective(t - gp.quicksum(lam[k] * v[k] for k in range(len(items))))
            lpm.optimize()
            if lpm.Status == GRB.OPTIMAL and lpm.ObjVal > 1e-4:
                lv = [lam[k].X for k in range(len(items))]
                model.cbCut(gp.quicksum(lv[k] * exprs[k] for k in range(len(items))
                                        if lv[k] > 1e-9) >= t.X)
                stats['subset_lp'] += 1

    def cb(model, where):
        try:
            if where == GRB.Callback.MIPSOL:
                val = dict(zip(X, model.cbGetSolution(list(X.values()))))
                cyc = G.violated_cycle(d, val)
                if cyc:
                    model.cbLazy(gp.quicksum(X[e] for e in cyc) <= len(cyc) - 1)
                    stats['cycle'] += 1
                ys = model.cbGetSolution([dd[0] for dd in dep]) if dep else []
                for (y, arcs, tset, _), yy in zip(dep, ys):
                    if yy <= 1e-6:
                        continue
                    S = classes.reachable(arcs, lambda e: val[e] > .5)
                    if not (S & tset):
                        edges = {e for s, t, e in arcs if s in S and t not in S}
                        model.cbLazy(y <= gp.quicksum(X[e] for e in edges))
                        stats['lazy'] += 1
            elif where == GRB.Callback.MIPNODE:
                if model.cbGet(GRB.Callback.MIPNODE_STATUS) != GRB.OPTIMAL:
                    return
                root = model.cbGet(GRB.Callback.MIPNODE_NODCNT) == 0
                if subsets and opts['subset'] == 'lp' and root and \
                        stats['subset_rounds'] < opts['subset_rounds']:
                    stats['subset_rounds'] += 1
                    subset_separation(model)
                if opts['nodecuts'] == 'none' or (opts['nodecuts'] == 'root' and not root) or not dep:
                    return
                val = dict(zip(X, model.cbGetNodeRel(list(X.values()))))
                ys = model.cbGetNodeRel([dd[0] for dd in dep])
                rval = {}
                if rlist:
                    rval = dict(zip([k for k, _ in rlist], model.cbGetNodeRel([v for _, v in rlist])))
                for (y, arcs, tset, key), yy in zip(dep, ys):
                    if yy <= 1e-5:
                        continue
                    rb = R.get(key, {}) if opts['joint'] else {}
                    H = nx.DiGraph()
                    H.add_nodes_from(['S', 'T'])
                    terms = []
                    for s, t, e in arcs:
                        cap = max(val[e], 0)
                        use_r = t in rb and rval[key, t] < cap
                        if use_r:
                            cap = max(rval[key, t], 0)
                        terms.append((s, t, e, use_r))
                        if H.has_edge(s, t):
                            H[s][t]['capacity'] += cap
                        else:
                            H.add_edge(s, t, capacity=cap)
                    for st in tset:
                        H.add_edge(st, 'T', capacity=1e6)
                    cut, (S, _) = nx.minimum_cut(H, 'S', 'T')
                    if cut < yy - 1e-5:
                        edges = {e for s, t, e, u in terms if s in S and t not in S and not u}
                        heads = {t for s, t, e, u in terms if s in S and t not in S and u}
                        rhs = sum(val[e] for e in edges) + sum(rval[key, t] for t in heads)
                        if rhs < yy - 1e-5:
                            model.cbCut(y <= gp.quicksum(X[e] for e in edges) +
                                        gp.quicksum(rb[t] for t in heads))
                            stats['user'] += 1
        except Exception as ex:
            errors.append(repr(ex))
            model.terminate()

    m.update()
    build = time.perf_counter() - wall0
    m.optimize(cb)
    if errors:
        raise RuntimeError(errors)
    has = m.SolCount > 0
    bound = m.ObjBound if (has or m.Status == GRB.OPTIMAL) and math.isfinite(m.ObjBound) else None
    res = dict(status=_STATUS.get(m.Status, str(m.Status)), objval_model=m.ObjVal if has else None,
               objbound=bound, runtime=m.Runtime, build=build, nodes=m.NodeCount,
               nvars=m.NumVars, nconstrs=m.NumConstrs, cuts=stats, options=opts,
               version=gp.gurobi.version(), seed=opts['seed'], threads=opts['threads'],
               **{'class': gclass})
    if has:
        val = {e: x.X for e, x in X.items()}
        if G.violated_cycle(d, val):
            raise RuntimeError('Returned graph is not a valid graph of the class')
        if not knowledge.satisfied(G, d, val, K):
            raise RuntimeError('Returned graph violates the background knowledge')
        loss = sum(q['w'] * (G.connected_vals(d, val, q['i'], q['j'], q['C']) != bool(q['p']))
                   for q in queries)
        obj = loss + knowledge.penalty_value(G, d, val, K)
        if obj > res['objval_model'] + 1e-5 * (1 + abs(obj)):
            raise RuntimeError(('Model objective below the oracle objective', obj, res['objval_model']))
        if m.Status == GRB.OPTIMAL and abs(obj - res['objval_model']) > 1e-5 * (1 + abs(obj)):
            raise RuntimeError(('Optimal objective differs from the oracle objective', obj, res['objval_model']))
        res['objval'] = obj
        res['oracle_loss'] = loss
        res['gap'] = None if bound is None else (abs(obj - bound) / abs(obj) if abs(obj) > 1e-10 else 0.0)
        res.update(G.matrices(d, val))
        res['connected'] = [int(G.connected_vals(d, val, q['i'], q['j'], q['C'])) for q in case['queries']]
    else:
        res['objval'] = res['gap'] = None
    res['warnings'] = warnings
    res['wall'] = time.perf_counter() - wall0
    m.dispose()
    env.dispose()
    return res

