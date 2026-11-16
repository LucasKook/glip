"""Local 3-node separation models for the subset consistency cuts.

LOCAL_QUERIES: the 6 queries of a 3-node graph (unordered pair, C empty or
the third node). local_models(cls): the distinct connection vectors of all
3-node graphs of the *marginal class*: the separations among a 3-node subset
S of a graph G equal those of the latent projection of G on S. ADMGs are
closed under latent projection (Richardson 2003), so for ADMGs the marginal
class is the ADMGs themselves; for DAGs and chain graphs the marginals are
again ADMGs (DAGs) or need a separate enumeration (see MARGINAL below).
Tested against the R oracle in tests/testthat/test-cuts-oracle.R.
"""
from functools import lru_cache

from . import classes

LOCAL_QUERIES = [(0, 1, ()), (0, 1, (2,)), (0, 2, ()), (0, 2, (1,)), (1, 2, ()), (1, 2, (0,))]

# class -> class whose 3-node members give a superset of the local models.
# Every DAG, ADMG and LWF chain graph has a faithful Gaussian distribution
# (DAG: Spirtes et al.; ADMG: via the canonical latent DAG; LWF CG: Pena 2011),
# so its restriction to S is the CI model of a regular Gaussian on 3
# variables. Those are exactly 11 models (all marginal correlations nonzero:
# at most one vanishing partial correlation, 4; exactly one zero marginal: 3;
# two: 3; none: 1), and the 3-node ADMG models are the same 11.
# Cyclic DGs and DMGs (linear cyclic models with d-separation, Spirtes 1995;
# bidirected edges as latent common causes) give the same 11 models; checked
# on all DGs and DMGs on 3 and 4 nodes and random larger ones.
MARGINAL = {'admg': 'admg', 'dag': 'admg', 'chain': 'admg', 'dg': 'admg', 'dmg': 'admg'}


@lru_cache(maxsize=None)
def local_models(cls):
    G = classes.get(MARGINAL[cls])
    models = set()
    for val in G.enumerate(3):
        models.add(tuple(int(G.connected_vals(3, val, i, j, list(C))) for i, j, C in LOCAL_QUERIES))
    return tuple(sorted(models))
