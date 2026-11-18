#!/usr/bin/env bash
# Reproduce the experiments of the paper on one machine (see
# inst/code/README.md): plan the tasks, build the encoding caches, run the
# tasks in parallel (one Gurobi thread each), collect the results and write
# the figures and tables. Finished tasks are skipped, so the script can be
# started again after an interruption.
#
# Usage (repository root):
#   bash inst/scripts/reproduce.sh <out> [exp]
#     exp: asp, simulation, benchmark, chain or a comma-separated list
#          (default: all four)
# Environment:
#   WORKERS         parallel tasks (default: number of cores)
#   SMOKE           if set (seconds), only one small task per experiment with
#                   this time limit, to check the installation
#   GLIP_ENCODINGS  "minlength" (default, as in the paper) or
#                   "minlength,cuts" to add the cut encoding of the package
#   GLIP_COMP_VENV, GLIP_PYTHON, GLIP_CACHE, GLIP_CLINGO: see inst/code/README.md
# Results: <out>/collected/*.csv, <out>/figures, <out>/tables
set -uo pipefail
out="${1:?usage: reproduce.sh <out> [asp,simulation,benchmark,chain]}"
exp="${2:-}"
WORKERS="${WORKERS:-$(getconf _NPROCESSORS_ONLN)}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
[ -f inst/code/tasks.R ] || { echo "run from the repository root"; exit 1; }
if [ -n "${SMOKE:-}" ]; then
  RR_SMOKE_EXP="$exp" Rscript inst/code/tasks.R smoke "$out" "$SMOKE" || exit 1
else
  [ -f "$out/tasks.csv" ] || Rscript inst/code/tasks.R plan "$out" || exit 1
  Rscript inst/code/tasks.R prebuild "$out" > "$out/prebuild.log" 2>&1 ||
    { tail "$out/prebuild.log"; exit 1; }
  mkdir -p "$out/logs"
  echo "== start $(date) workers=$WORKERS exp=${exp:-all}"
  seq 1 "$WORKERS" | xargs -P "$WORKERS" -I{} sh -c \
    "Rscript inst/code/tasks.R run '$out' {} $WORKERS $exp >> '$out/logs/chunk-{}.log' 2>&1 || echo 'chunk {} exit '\$?"
  echo "== tasks done $(date)"
  Rscript inst/code/tasks.R collect "$out"
fi
n_err=$(find "$out" -name error.txt | wc -l | tr -d ' ')
[ "$n_err" = 0 ] || echo "failed tasks ($n_err): see <out>/<task>/error.txt; retry with: Rscript inst/code/tasks.R retry $out <task> ..."
Rscript inst/code/vis.R "$out/collected" "$out/figures" "$out/tables"
