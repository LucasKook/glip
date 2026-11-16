"""Graph classes for cutsolver.py: edge variables, state automata, validity.

Edge keys: ('d', a, b) directed a -> b, ('b', a, b) bidirected (a < b).
A state is (node, mark) with mark 1 = arrived with an arrowhead at node and
mark 0 = arrived with a tail; 'S' is the source.
"""
import itertools
from collections import deque

import networkx as nx



def admg_state_arcs(d, C, i):
    """m-connection states (node, arrival mark); validated against
    is_m_separated() (2026-11-13: 10,500 queries, 120 random ADMGs)."""
    arcs = []
    for u in range(d):
        if u == i:
            continue
        arcs.extend([('S', (u, 1), ('d', i, u)),
                     ('S', (u, 0), ('d', u, i)),
                     ('S', (u, 1), ('b', min(i, u), max(i, u)))])
    for v in range(d):
        for h in (0, 1):
            head_ok = (h == 1 and v in C) or (h == 0 and v not in C)
            for u in range(d):
                if u == v:
                    continue
                if v not in C:
                    arcs.append(((v, h), (u, 1), ('d', v, u)))
                if head_ok:
                    arcs.extend([((v, h), (u, 0), ('d', u, v)),
                                 ((v, h), (u, 1), ('b', min(u, v), max(u, v)))])
    return arcs


def reachable(arcs, present):
    out = {}
    for s, t, e in arcs:
        if present(e):
            out.setdefault(s, []).append(t)
    seen = {'S'}
    queue = deque(['S'])
    while queue:
        for t in out.get(queue.popleft(), []):
            if t not in seen:
                seen.add(t)
                queue.append(t)
    return seen


def qkey(q):
    i, j = q['i'], q['j']
    return (min(i, j), max(i, j), tuple(sorted(q['C'])))


def _directed_cycle(d, val):
    H = nx.DiGraph()
    H.add_nodes_from(range(d))
    H.add_edges_from((a, b) for (k, a, b), v in val.items() if k == 'd' and v > .5)
    try:
        return [('d', a, b) for a, b in nx.find_cycle(H)]
    except nx.NetworkXNoCycle:
        return []


class ADMG:
    """Acyclic directed mixed graphs under m-separation."""
    name = 'admg'

    @staticmethod
    def edges(d):
        return [('d', a, b) for a in range(d) for b in range(d) if a != b] + \
               [('b', a, b) for a in range(d) for b in range(a + 1, d)]

    @staticmethod
    def static_rows(d):
        return [[('d', a, b), ('d', b, a)] for a in range(d) for b in range(a + 1, d)]

    @staticmethod
    def start_values(d, start):
        A, B = start['M1'], start['M2']
        assert len(A) == d and all(len(r) == d for r in A)
        assert all(A[a][a] == 0 and B[a][a] == 0 for a in range(d))
        assert all(B[a][b] == B[b][a] for a in range(d) for b in range(d))
        return {e: int(A[e[1]][e[2]] if e[0] == 'd' else B[e[1]][e[2]]) for e in ADMG.edges(d)}

    violated_cycle = staticmethod(_directed_cycle)

    @staticmethod
    def state_arcs(d, C, i):
        return admg_state_arcs(d, C, i)

    @staticmethod
    def target_states(j):
        return [(j, 0), (j, 1)]

    @classmethod
    def connected_vals(cls, d, val, i, j, C):
        S = reachable(cls.state_arcs(d, set(C), i), lambda e: val[e] > .5)
        return any(st in S for st in cls.target_states(j))

    @staticmethod
    def matrices(d, val):
        M1, M2 = [[0] * d for _ in range(d)], [[0] * d for _ in range(d)]
        for (k, a, b), v in val.items():
            x = int(round(v))
            if k == 'd':
                M1[a][b] = x
            else:
                M2[a][b] = M2[b][a] = x
        return dict(M1=M1, M2=M2)

    # Background knowledge (semantics of the classic glip encodings)
    knowledge_keys = ('forbidden', 'required', 'forbidden_bi', 'required_bi',
                      'forbidden_anc', 'required_anc')
    coef_keys = ('coef', 'coef_bi')

    @staticmethod
    def xexpr(i, j):
        """Edge variables summing to glip's directed indicator x_ij."""
        return [('d', i, j)]

    @staticmethod
    def bi_var(i, j):
        return ('b', min(i, j), max(i, j))

    @classmethod
    def penalty_terms(cls, d):
        """(edge, coefficient): edge_penalty applies to every edge indicator."""
        return [(e, 1.0) for e in cls.edges(d)]

    @staticmethod
    def anc_arcs(d):
        """Arcs (k, l, edge) along which k is an ancestor of l (directed paths)."""
        return [(k, l, ('d', k, l)) for k in range(d) for l in range(d) if k != l]

    @classmethod
    def enumerate(cls, d):
        """All graphs of the class on d nodes, as edge-value dicts."""
        E = cls.edges(d)
        for bits in itertools.product((0, 1), repeat=len(E)):
            val = dict(zip(E, bits))
            if any(sum(val[e] for e in row) > 1 for row in cls.static_rows(d)):
                continue
            if cls.violated_cycle(d, val):
                continue
            yield val


class DAG(ADMG):
    """DAGs under d-separation: the ADMG automaton without bidirected arcs."""
    name = 'dag'
    knowledge_keys = ('forbidden', 'required', 'forbidden_anc', 'required_anc')
    coef_keys = ('coef',)

    @staticmethod
    def edges(d):
        return [('d', a, b) for a in range(d) for b in range(d) if a != b]

    @staticmethod
    def start_values(d, start):
        A = start['M1']
        assert len(A) == d and all(len(r) == d for r in A) and all(A[a][a] == 0 for a in range(d))
        assert not any(v for row in start.get('M2', [[0]]) for v in row), 'DAG start has bidirected edges'
        return {e: int(A[e[1]][e[2]]) for e in DAG.edges(d)}

    @staticmethod
    def state_arcs(d, C, i):
        return [a for a in admg_state_arcs(d, C, i) if a[2][0] == 'd']

    @staticmethod
    def matrices(d, val):
        M1 = [[0] * d for _ in range(d)]
        for (k, a, b), v in val.items():
            M1[a][b] = int(round(v))
        return dict(M1=M1, M2=[[0] * d for _ in range(d)])


class Chain(ADMG):
    """LWF chain graphs under c-separation (glip's csep()).

    Edge keys ('d', a, b) for a -> b and ('u', a, b) for a - b (a < b).
    Route criterion (Studeny 1998): a route is active given C iff every
    head-to-head section (maximal undirected subroute with arrowheads into
    both ends) contains a node in C and every other section contains none.
    States: (v, 0) = in a section entered with a tail (or the source
    section), no node of C so far; (v, 1, c) = in a section entered with an
    arrowhead, c = 1 iff the section so far contains a node of C. The source
    state (i, 0) is 'S'. Targets: (j, 0) and (j, 1, 0) (end sections are not
    head-to-head, so must avoid C). Validity: no semidirected cycle.
    """
    name = 'chain'
    knowledge_keys = ('forbidden', 'required', 'forbidden_undirected',
                      'required_undirected', 'forbidden_anc', 'required_anc')
    coef_keys = ('coef',)

    @staticmethod
    def xexpr(i, j):
        """glip's chain indicator x_ij = 1 iff i -> j or i - j."""
        return [('d', i, j), ('u', min(i, j), max(i, j))]

    @classmethod
    def penalty_terms(cls, d):
        """The penalty applies to every x_ij, so an undirected edge counts twice."""
        return [(e, 2.0 if e[0] == 'u' else 1.0) for e in cls.edges(d)]

    @staticmethod
    def anc_arcs(d):
        """Anterior paths: directed edges forward and undirected edges both ways."""
        out = []
        for k in range(d):
            for l in range(d):
                if k != l:
                    out += [(k, l, ('d', k, l)), (k, l, ('u', min(k, l), max(k, l)))]
        return out

    @staticmethod
    def edges(d):
        return [('d', a, b) for a in range(d) for b in range(d) if a != b] + \
               [('u', a, b) for a in range(d) for b in range(a + 1, d)]

    @staticmethod
    def static_rows(d):
        return [[('d', a, b), ('d', b, a), ('u', a, b)] for a in range(d) for b in range(a + 1, d)]

    @staticmethod
    def start_values(d, start):
        G = start['G'] if 'G' in start else start['M1']
        assert len(G) == d and all(len(r) == d for r in G) and all(G[a][a] == 0 for a in range(d))
        val = {}
        for a in range(d):
            for b in range(d):
                if a != b:
                    val['d', a, b] = int(G[a][b] == 1 and G[b][a] == 0)
                if a < b:
                    val['u', a, b] = int(G[a][b] == 1 and G[b][a] == 1)
        return val

    @staticmethod
    def violated_cycle(d, val):
        """Edge variables of a semidirected cycle, or []."""
        H = nx.DiGraph()
        H.add_nodes_from(range(d))
        dir_edges = []
        for (k, a, b), v in val.items():
            if v > .5:
                if k == 'd':
                    H.add_edge(a, b, var=('d', a, b))
                    dir_edges.append((a, b))
                else:
                    H.add_edge(a, b, var=('u', a, b))
                    H.add_edge(b, a, var=('u', a, b))
        for a, b in dir_edges:
            try:
                path = nx.shortest_path(H, b, a)
            except nx.NetworkXNoPath:
                continue
            cyc = {('d', a, b)} | {H[u][w]['var'] for u, w in zip(path, path[1:])}
            return sorted(cyc)
        return []

    @staticmethod
    def state_arcs(d, C, i):
        C = set(C)
        states = [(v, 0) for v in range(d) if v not in C] + \
                 [(v, 1, c) for v in range(d) for c in (0, 1) if c or v not in C]
        arcs = []
        for s in states:
            v = s[0]
            src = 'S' if s == (i, 0) else s
            tail_entered = len(s) == 2
            c = 0 if tail_entered else s[2]
            for w in range(d):
                if w == v:
                    continue
                wc = int(w in C)
                if c == 0:  # leave with a tail at v: non-collider, section free of C
                    arcs.append((src, (w, 1, wc), ('d', v, w)))
                # leave with an arrowhead at v: collider iff entered with an arrowhead
                if (tail_entered or c == 1) and not wc:
                    arcs.append((src, (w, 0), ('d', w, v)))
                u = ('u', min(v, w), max(v, w))
                if tail_entered:
                    if not wc:
                        arcs.append((src, (w, 0), u))
                else:
                    arcs.append((src, (w, 1, int(c or wc)), u))
        return [(s, t, e) for s, t, e in arcs if t != (i, 0)]

    @staticmethod
    def target_states(j):
        return [(j, 0), (j, 1, 0)]

    @staticmethod
    def matrices(d, val):
        G = [[0] * d for _ in range(d)]
        for (k, a, b), v in val.items():
            if round(v):
                G[a][b] = 1
                if k == 'u':
                    G[b][a] = 1
        return dict(G=G, M1=G, M2=[[0] * d for _ in range(d)])


class DG(DAG):
    """Directed graphs (cycles, including 2-cycles, allowed) under
    d-separation, as in glip's DG modes: the DAG automaton without any
    acyclicity constraint (the walk criterion is the same in cyclic graphs)."""
    name = 'dg'

    @staticmethod
    def static_rows(d):
        return []

    @staticmethod
    def violated_cycle(d, val):
        return []


class DMG(ADMG):
    """Directed mixed graphs (cycles allowed) under m-separation, as in
    glip's DMG modes: the ADMG automaton without acyclicity."""
    name = 'dmg'

    @staticmethod
    def static_rows(d):
        return []

    @staticmethod
    def violated_cycle(d, val):
        return []


CLASSES = {'admg': ADMG, 'dag': DAG, 'chain': Chain, 'dg': DG, 'dmg': DMG}


def get(name):
    if name not in CLASSES:
        raise ValueError(f'Unsupported class {name!r}; available: {sorted(CLASSES)}')
    return CLASSES[name]
