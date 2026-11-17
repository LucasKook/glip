#!/usr/bin/env bash
# Local: start run.sh on a running job (setup.sh must have finished).
# Usage (repository root): WORKERS=32 bash inst/code/rerun/server/deploy-run.sh <ssh_port> <code_dir_name>
set -euo pipefail
port="${1:?ssh port}"; name="${2:?code dir name (see setup.log)}"
base=/work/glip
scp -P "$port" -o ConnectTimeout=20 -o StrictHostKeyChecking=accept-new \
  inst/code/rerun/server/run.sh "ucloud@ssh.cloud.sdu.dk:$base/run.sh"
ssh -p "$port" -o ConnectTimeout=20 ucloud@ssh.cloud.sdu.dk bash -s -- "$name" "${WORKERS:-32}" "${END_JOB:-1}" <<'REMOTE'
set -euo pipefail
name=$1; workers=$2; endjob=$3; base=/work/glip
[ -d "$base/$name" ] || { echo "code $name missing"; exit 1; }
[ -f "$base/env.sh" ] || { echo "env.sh missing (run setup first)"; exit 1; }
[ -f "$base/runs/rerun-2026-11/tasks.csv" ] || { echo "tasks.csv missing"; exit 1; }
if pgrep -f "rerun/tasks.R run" > /dev/null; then echo "rerun already running"; exit 1; fi
[ ! -e "$base/run.log" ] || mv "$base/run.log" "$base/run.log.$(date +%s)"
WORKERS=$workers END_JOB=$endjob setsid nohup bash "$base/run.sh" "$name" > "$base/run.log" 2>&1 < /dev/null &
sleep 60
tail -3 "$base/run.log"
echo "finished tasks: $(find "$base/runs/rerun-2026-11" -name result.rds | wc -l)"
REMOTE
echo "status: ssh -p $port ucloud@ssh.cloud.sdu.dk tail /work/glip/run.log"
