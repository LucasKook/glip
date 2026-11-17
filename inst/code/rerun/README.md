# Rerun of the paper experiments

Rerun of all experiments of the GLIP paper (arXiv:2601.20589, v2) with the
new min-length encodings and the cut encoding, as resumable tasks for a
machine with many cores (one Gurobi thread per solve).

| Experiment | Paper | Task function | Grid |
| --- | --- | --- | --- |
| asp | Sec. 5.1, Fig. 7; App. D.1 | `exp-asp.R` | DAG, ADMG; full d 3 to 9, weak k = 1 d 3 to 14; n 10,000; 20 reps; 600 s |
| simulation | Sec. 5.2, Figs. 8, 9, Table 3; App. D.2 | `exp-simulation.R` | DAG, ADMG; full d 3 to 9, weak k = 1, 2 d 3 to 12; n 400, 10,000; 30 reps; 300 s |
| benchmark | Sec. 5.2.2, Table 4; App. D.3, Table 6 | `exp-benchmark.R` | ADMG; alarm, asia, child, hepar2, sachs; d_max 6, 8; k = 1, 2, d - 2; min(choose(d, d_max), 50) subsets; 300 s |
| chain | App. D.4, Figs. 15, 16 | `exp-chain.R` | LWF chain graphs, d 3 to 10, full; n 1,000, 100,000; 30 reps; 300 s |

Every GLIP fit is run with both encodings (`minlength`: `dcon_lean_optim()`,
`admg_lean_optim(mode = "admg-lean")`, `chain_lean_optim()`; `cuts`:
`cuts_optim()`) on the same tests, decisions, weights and warm start, with
`Threads = 1`. Data generation and seeds follow the original drivers in
`inst/code/run-*.R`. Differences from the paper drivers: PC/FCI warm starts
on weak banks use `gaussCItest` on the data instead of the full d - 2 test
bank; the benchmark uses all or 50 subsets; the competitors are run once per
data set.

## Usage (from the repository root)

```sh
Rscript inst/code/rerun/tasks.R plan <out>             # <out>/tasks.csv
seq 1 32 | xargs -P 32 -I{} Rscript inst/code/rerun/tasks.R run <out> {} 32
Rscript inst/code/rerun/tasks.R collect <out>          # <out>/collected/*.csv
Rscript inst/code/rerun/tasks.R smoke <out> 10         # one task per experiment
Rscript inst/code/rerun/vis.R <out>/collected <figures> <tables>
Rscript inst/code/rerun/intro-example.R <outdir>       # Example 1 (local)
```

`vis.R` writes the figures and tables of the paper from the collected CSVs
(both encodings as methods; certified = status OPTIMAL) and a comparison of
the two encodings (`*-encodings.csv`, paired `*-encodings-paired.csv`).

Server (UCloud, drive `/work/glip`): `server/deploy-setup.sh <port>` packs
the working tree and runs `server/setup.sh` on a job (R, Gurobi 13.0.3,
licence check, R packages, lcd, Python environments, clingo, package tests,
plan, cache prebuild, smoke); `server/deploy-run.sh <port> <code_dir>`
starts `server/run.sh` (32 single-thread workers, resumable, ends the job
after collecting); `server/copy-back.sh <port>` copies the results with an
MD5 check to `inst/results/rerun-2026-11/raw`.

Each task writes `<out>/<id>/` (inputs, one rds per fit, `result.rds`);
finished tasks and fits are skipped, so `run` can be restarted at any time.
Errors are written to `<out>/<id>/error.txt`; delete the task directory to
rerun it. Versions of R, packages, Gurobi, gurobipy, the Python competitors
and clingo are attached to every result.

Environment: `GLIP_PYTHON` (Python with gurobipy for the cut encoding),
`GLIP_COMP_VENV` (virtualenv with `dagma` and `CausalDisco`), `GLIP_CACHE`
(cache of the min-length encodings), `GLIP_ASP_HOME`, `GLIP_CLINGO` (see
`inst/asp-glip/README.md`). R packages: pcalg, bnlearn, dagitty, igraph,
mvtnorm, reticulate, gurobi, deal, combinat, hash, graph, gRbase, ggm and lcd
(archived on CRAN: `install-lcd.sh`).
