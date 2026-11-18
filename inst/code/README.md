# Experiments of the paper

This folder contains the code for the experiments in Kook and Mogensen
(2026), *Exact Graph Learning via Integer Programming* (arXiv:2601.20589).
Every experiment is split into small, independent tasks (one data set
each), so the computations can run in parallel and be resumed after an
interruption.

| Experiment | Paper | Code |
| --- | --- | --- |
| `asp`: timing against the ASP method of Hyttinen et al. (2014) | Sec. 5.1, Fig. 7; App. D.1 | `exp-asp.R` |
| `simulation`: simulated DAGs and ADMGs | Sec. 5.2, Figs. 8, 9, Table 3; App. D.2 | `exp-simulation.R` |
| `benchmark`: benchmark data sets | Sec. 5.2.2, Table 4; App. D.3, Table 6 | `exp-benchmark.R` |
| `chain`: chain graphs | App. D.4, Figs. 15, 16, Table 7 | `exp-chain.R` |
| Example 1 | Sec. 1.1, Fig. 1 | `intro-example.R` |

## Requirements

- A clone of this repository. The scripts load **glip** from the clone
  (via **pkgload**), so its dependencies and Gurobi with the **gurobi** R
  package must be installed (see the main README).
- R packages: pkgload, pcalg, bnlearn, dagitty, igraph, mvtnorm, reticulate,
  deal, combinat, hash, graph, gRbase, ggm, tidyverse, scales, knitr,
  kableExtra, and lcd 0.7-3, which is archived on CRAN:
  `bash inst/code/install-lcd.sh <R library>` installs it (patched for
  current versions of R).
- Python with **dagma** and **CausalDisco** for the competitors NOTEARS and
  R2SORT (`pip install dagma CausalDisco`); point `GLIP_COMP_VENV` to the
  virtualenv if it is not the default Python of **reticulate**.
- The ASP comparison uses the code of Hyttinen et al. (2014) with the clingo
  binary contained in `inst/asp/hyttinen2014uai_ver6.zip` (Linux, x86-64);
  see `inst/asp-glip/README.md`. On other systems, set `GLIP_CLINGO` to a
  clingo executable (e.g. from `pip install clingo`); the ASP timings then
  depend on that version.

## Running the experiments

From the repository root:

```sh
SMOKE=10 bash inst/scripts/reproduce.sh results/smoke  # check the installation (a few minutes)
bash inst/scripts/reproduce.sh results/paper           # all experiments
bash inst/scripts/simulation.sh results/paper          # a single experiment (also asp, benchmark, chain)
bash inst/scripts/intro-example.sh                     # Example 1 (small)
```

`reproduce.sh` uses all cores by default (`WORKERS=<n>` to change this),
with one Gurobi thread per task. The time limits are those of the paper
(600 s per fit for `asp`, 300 s otherwise). If every fit hits its time
limit, all experiments together take about 550 CPU hours. Finished tasks
are skipped, so the script can be started again after an interruption.

The results are written to `<out>/collected/*.csv` (one file per
experiment), the figures to `<out>/figures` and the tables to
`<out>/tables`. To recreate the figures and tables from collected results:

```sh
Rscript inst/code/vis.R <out>/collected <figures> <tables>
```

The individual steps of `reproduce.sh`, e.g. for a cluster:

```sh
Rscript inst/code/tasks.R plan <out>                    # list of tasks, <out>/tasks.csv
Rscript inst/code/tasks.R prebuild <out>                # encoding caches (run once, before the tasks)
Rscript inst/code/tasks.R run <out> <chunk> <n_chunks>  # every n_chunks-th task; run chunks 1..n_chunks in parallel
Rscript inst/code/tasks.R collect <out>                 # <out>/collected/*.csv
Rscript inst/code/tasks.R retry <out> <task> ...        # run failed tasks again
```

Each task writes its inputs, fits and results to `<out>/<task>/`; a failed
task writes the error message to `<out>/<task>/error.txt`. The versions of
R, the packages, Gurobi, Python libraries and clingo are stored with the
results.

## Details

- Data generation and seeds are those of the paper, so the tasks generate
  the data sets of the paper. The benchmark experiment uses
  min{C(d, d_max), 50} node subsets per data set.
- GLIP uses the min-length encodings of the paper (`learn_graph()` with
  `encoding = "minlength"`). The package also contains a cut encoding
  (`encoding = "cuts"`), which is not part of the paper and is not run by
  default. To fit it in addition, set `GLIP_ENCODINGS=minlength,cuts` and
  `GLIP_PYTHON` to a Python with the packages of
  `inst/python/requirements.txt`; the figures then show both encodings and
  the tables compare them.
- Timings and the fractions of certified optimal solutions depend on the
  hardware, the Gurobi version and the number of threads.
- Environment variables: `GLIP_ENCODINGS`, `GLIP_COMP_VENV`, `GLIP_PYTHON`,
  `GLIP_CLINGO`, `GLIP_CACHE` (cache directory of the encodings, default
  `.cache-glip` in the repository root) and `GLIP_ASP_HOME` (where the ASP
  code is unzipped, default a temporary directory).
