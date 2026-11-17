#!/usr/bin/env bash
# Local: copy the rerun results (task folders, collected CSVs, logs, tasks.csv)
# from /work/glip/runs/rerun-2026-11 to inst/results/rerun-2026-11/raw and
# check the MD5 sum of the archive.
# Usage (repository root): bash inst/code/rerun/server/copy-back.sh <ssh_port>
set -euo pipefail
port="${1:?ssh port}"
dest=inst/results/rerun-2026-11/raw
mkdir -p "$dest"
stamp=$(date +%Y%m%d-%H%M)
ssh -p "$port" -o ConnectTimeout=20 ucloud@ssh.cloud.sdu.dk bash -s -- "$stamp" <<'REMOTE'
set -euo pipefail
cd /work/glip/runs
tar -czf "/work/glip/rerun-2026-11-$1.tar.gz" rerun-2026-11
md5sum "/work/glip/rerun-2026-11-$1.tar.gz"
REMOTE
remote_md5=$(ssh -p "$port" ucloud@ssh.cloud.sdu.dk md5sum "/work/glip/rerun-2026-11-$stamp.tar.gz" | cut -d' ' -f1)
scp -P "$port" "ucloud@ssh.cloud.sdu.dk:/work/glip/rerun-2026-11-$stamp.tar.gz" "$dest/"
local_md5=$(md5 -q "$dest/rerun-2026-11-$stamp.tar.gz" 2>/dev/null || md5sum "$dest/rerun-2026-11-$stamp.tar.gz" | cut -d' ' -f1)
[ "$remote_md5" = "$local_md5" ] || { echo "MD5 mismatch"; exit 1; }
echo "$local_md5  rerun-2026-11-$stamp.tar.gz" >> "$dest/MD5SUMS"
(cd "$dest" && tar -xzf "rerun-2026-11-$stamp.tar.gz")
echo "copied and verified: $dest/rerun-2026-11-$stamp.tar.gz"
