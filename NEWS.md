# glip 0.0-1.9000

## New features

* Background knowledge as hard constraints for all encodings and through
  `learn_graph()`:
  * `forbidden`, `required`: forbidden and required directed edges (all
    modes).
  * `forbidden_bi`, `required_bi`: forbidden and required bidirected edges
    (ADMG and DMG modes).
  * `forbidden_undirected`, `required_undirected`: forbidden and required
    undirected edges (chain graph modes).
  * `forbidden_anc`, `required_anc`: forbidden and required ancestral
    relations (directed paths; anterior paths for chain graphs).
  * `edge_constraints`: general linear constraints on the edge indicators
    (`coef` for directed edges, and `coef_bi` for bidirected edges in ADMG
    and DMG modes).
  * `edge_penalty`: objective coefficient per edge, e.g., to break ties in
    favour of sparse graphs.

  The constraints are imposed as bounds or additional rows after the cached
  constraint matrices are loaded, so caches are unaffected. Arguments for a
  type of edge that the class of graphs does not contain (e.g.,
  `forbidden_bi` for DAGs) raise an error.

* `learn_graph()` passes the background knowledge to the PC/FCI warm start
  as far as possible: pairs that cannot or must be adjacent as `fixedGaps`
  and `fixedEdges` (FCI then uses `pcalg::fci()` instead of
  `pcalg::fciPlus()`), and implied orientations of the PC result via
  `pcalg::addBgKnowledge()`. See `?learn_graph` for what cannot be passed.

* If the solver stops without a feasible solution (infeasible constraints,
  or early termination via `gurobi_args`), the optimizers return
  `graph = NULL` together with the solver output instead of failing.

## Bug fixes

* `compute_largest_cg()` failed for chain graphs with directed edges.
* `dcon_optim()` now matches its `mode` argument.
* `dcon_lean_optim()` and `dcon_optim()` include `mode` in the names of the
  cache files, so DAG and DG modes no longer share cached constraint
  matrices (which include the acyclicity constraints for DAGs only).
  Existing caches for these functions are rebuilt on first use.
* `tests/testthat.R` loads and checks glip.
* `learn_graph()`, `falsify_graph()` and the oracle tests failed with
  `parallel = TRUE` unless future.apply was attached, because the package
  imported nothing in its NAMESPACE. Functions from methods, stats, utils,
  future, future.apply, comets and pcalg are now imported explicitly.
* With `comets = FALSE` and a test other than `"gaussCItest"`, a missing
  bnlearn package now raises an error. Previously, every test failed
  silently (with a warning) and returned a p-value close to one.
* Dependencies: bnlearn is now suggested, graph, methods and utils are
  imported, and the unused Formula, mvtnorm and stringr were removed.
* The methods of `prf1()` now accept `...` like the generic.
* `dcon_optim()`, `dag_optim()`, `admg_lean_optim()`, `admg_optim()` and
  `chain_optim()` (modes `"dag-dcon"`, `"dag-dc"`, `"dg-dc"`,
  `"admg-lean"`, `"dmg-lean"`, `"admg-dc"`, `"dmg-dc"`, `"chain-dcon"`)
  failed for `max_size = 0`.
* `mode = "dg-dcon"` is now supported by `learn_graph()` and the internal
  helpers; `dcon_optim()` already supported it.
* `learn_graph()` checks `mode` and lists the supported modes in the error
  and the documentation. The warning that learning relies on d-separation
  now applies to all modes for DGs and DMGs.
* Without the gurobi package, the optimizers now stop immediately with an
  informative error instead of warning, building the model, and failing
  with "object 'sol' not found".
* `R CMD check` passes without warnings or notes.
* `learn_graph()` with `mode = "dag-dc"` or `"dg-dc"` failed because it
  passed an FCI (ADMG) warm start to a DAG optimizer; all DAG and DG modes
  now use a PC warm start. No warm start is computed for `"chain-dcon"`,
  whose optimizer does not use one.
