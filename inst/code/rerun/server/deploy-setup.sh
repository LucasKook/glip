#!/usr/bin/env bash
# Local: pack the working tree (tracked + untracked files that are not
# ignored, plus the benchmark data sets in inst/data/datasets), copy it to
# /work/glip/releases on a running UCloud job and start setup.sh there.
# Usage (repository root): bash inst/code/rerun/server/deploy-setup.sh <ssh_port>
set -euo pipefail
port="${1:-}"; [ -n "$port" ] || [ "${DRY:-0}" = 1 ] || { echo "usage: deploy-setup.sh <ssh_port>"; exit 1; }
base=/work/glip
ssh_cmd=(ssh -p "$port" -o ConnectTimeout=20 -o StrictHostKeyChecking=accept-new ucloud@ssh.cloud.sdu.dk)
scp_cmd=(scp -P "$port" -o ConnectTimeout=20 -o StrictHostKeyChecking=accept-new)
head=$(git rev-parse --short HEAD)
tmpdir=$(mktemp -d)
{ git ls-files; git ls-files --others --exclude-standard; ls inst/data/datasets/*.txt; echo inst/code/test-one-thread.R; } |
  sort -u | grep -v '__pycache__' | while read -r f; do [ -f "$f" ] && echo "$f"; done > "$tmpdir/files.txt"
hash=$(tar -cf - -T "$tmpdir/files.txt" | shasum -a 256 | cut -c1-8)
dirty=$(git status --porcelain --untracked-files=no | wc -l | tr -d ' ')
name="glip-$head-$hash"
mkdir -p "$tmpdir/stage/glip"
rsync -a --files-from="$tmpdir/files.txt" ./ "$tmpdir/stage/glip/"
echo "$head (tracked files modified: $dirty, content hash $hash)" > "$tmpdir/stage/glip/REVISION"
COPYFILE_DISABLE=1 tar -czf "$tmpdir/$name.tar.gz" -C "$tmpdir/stage" glip
md5=$(md5 -q "$tmpdir/$name.tar.gz" 2>/dev/null || md5sum "$tmpdir/$name.tar.gz" | cut -d' ' -f1)
echo "$name md5 $md5 ($(wc -l < "$tmpdir/files.txt" | tr -d ' ') files)"
if [ "${DRY:-0}" = 1 ]; then echo "dry run: $tmpdir/$name.tar.gz"; exit 0; fi
"${ssh_cmd[@]}" "mkdir -p $base/releases"
"${scp_cmd[@]}" "$tmpdir/$name.tar.gz" "ucloud@ssh.cloud.sdu.dk:$base/releases/"
"${scp_cmd[@]}" inst/code/rerun/server/setup.sh inst/code/rerun/server/run.sh "ucloud@ssh.cloud.sdu.dk:$base/"
"${ssh_cmd[@]}" bash -s -- "$name" "$md5" "$base" <<'REMOTE'
set -euo pipefail
name=$1; md5=$2; base=$3
cd "$base"
[ "$(md5sum "releases/$name.tar.gz" | cut -d' ' -f1)" = "$md5" ] || { echo "md5 mismatch"; exit 1; }
if [ ! -e "$name" ]; then x=$(mktemp -d); tar -xzf "releases/$name.tar.gz" -C "$x"; mv "$x/glip" "$name"; fi
chmod +x setup.sh run.sh
[ ! -e setup.log ] || mv setup.log "setup.log.$(date +%s)"
setsid nohup bash "$base/setup.sh" "$name" > "$base/setup.log" 2>&1 < /dev/null &
echo "setup started for $name (log $base/setup.log)"
REMOTE
echo "follow: ssh -p $port ucloud@ssh.cloud.sdu.dk tail -f $base/setup.log"
