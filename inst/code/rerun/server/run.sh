#!/usr/bin/env bash
# Server (32 vCPU job): run all tasks of the rerun with WORKERS single-thread
# workers, collect, and end the job END_DELAY seconds later (END_JOB=0 keeps
# it). Resumable: finished tasks and fits are skipped.
# Usage: setsid nohup bash /work/glip/run.sh <code_dir_name> > /work/glip/run.log 2>&1 &
set -uo pipefail
base=/work/glip
code="$base/$1"
WORKERS="${WORKERS:-32}"; CHUNKS="${CHUNKS:-$WORKERS}"; EXP="${EXP:-}"
"$base/install-r.sh" > /tmp/install-r.log 2>&1 || { echo "install-r failed"; exit 1; }
. "$base/env.sh"
cd "$code"
out="$base/runs/rerun-2026-11"
[ -f "$out/tasks.csv" ] || { echo "no tasks.csv in $out"; exit 1; }
mkdir -p "$out/logs"
echo "== start $(date) workers=$WORKERS chunks=$CHUNKS exp=${EXP:-all} tasks=$(($(wc -l < "$out/tasks.csv") - 1))"
seq 1 "$CHUNKS" | xargs -P "$WORKERS" -I{} sh -c \
  "Rscript inst/code/rerun/tasks.R run '$out' {} $CHUNKS $EXP >> '$out/logs/chunk-{}.log' 2>&1 || echo 'chunk {} exit '\$?"
echo "== tasks done $(date)"
Rscript inst/code/rerun/tasks.R collect "$out"
echo "== all done $(date)"
if [ "${END_JOB:-1}" = 1 ]; then sleep "${END_DELAY:-1800}"; sync; pkill -x ttyd; fi
