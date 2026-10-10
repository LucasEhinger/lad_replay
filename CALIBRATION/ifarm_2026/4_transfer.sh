#!/bin/bash
# Step 4 (subMIT): move finished replays from the ifarm /volatile to subMIT with Globus, check them, and delete them
# on the ifarm. Uses only the Globus CLI (no ssh). A replay counts as finished when its report file exists.
#   [CHUNK_GB=300] ./4_transfer.sh [dest_dir]        one pass; run it again (or from a loop) while the workflow runs
set -euo pipefail
export PATH=$HOME/.local/bin:$PATH
JLAB=b0fca1ad-f485-4a00-8fcd-bca0b93a2a1c          # jlab#gw1
SUBMIT=dc702c96-b104-452c-bb45-db43ef85f390        # SubMIT Globus
SRC=/expphy/volatile/hallc/c-lad/ehingerl/lad_replay_2026
DEST=${1:-/ceph/submit/data/user/e/ehingerl/hallc/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/replay_2026}
FALLBACK=/ceph/submit/data/group/lad/lad_replay_2026
CHUNK_GB=${CHUNK_GB:-300}                         # per Globus task; files are deleted on the ifarm after each chunk
# switch to the group area when the user area has less than 1 TB left
quota_free=$(( $(getfattr --only-values -n ceph.quota.max_bytes /ceph/submit/data/user/e/ehingerl) - $(getfattr --only-values -n ceph.dir.rbytes /ceph/submit/data/user/e/ehingerl) ))
if [ "$quota_free" -lt 1000000000000 ]; then DEST=$FALLBACK; echo "user area has < 1 TB left; using $DEST"; fi
mkdir -p "$DEST" "$DEST/reports"
tmp=$(mktemp -d)
globus ls -F json "$JLAB:$SRC/REPORT_OUTPUT/LAD_COIN/PRODUCTION/" | python3 -c '
import json,sys,re
for e in json.load(sys.stdin)["DATA"]:
    m=re.match(r"replayReport_LAD_coin_production_(\d+)_(\d+)_(\d+)_(-?\d+)\.report$",e["name"])
    if m: print(e["name"], "LAD_COIN_%s_%s_%s_%s.root"%m.groups())' > $tmp/done.txt
globus ls -F json "$JLAB:$SRC/ROOTfiles/LAD_COIN/PRODUCTION/" | python3 -c '
import json,sys
for e in json.load(sys.stdin)["DATA"]: print(e["name"], e["size"])' > $tmp/roots.txt
: > $tmp/all.txt
while read -r rep root; do
  size=$(awk -v f="$root" '$1==f{print $2}' $tmp/roots.txt)
  [ -n "$size" ] || continue
  echo "$root $rep $size" >> $tmp/all.txt
done < $tmp/done.txt
n=$(wc -l < $tmp/all.txt)
if [ "$n" -eq 0 ]; then echo "nothing finished to move"; rm -rf $tmp; exit 0; fi
echo "moving $n finished replays ($(awk '{s+=$2} END{printf "%.1f GB", s/1e9}' $tmp/roots.txt) on the ifarm in total), in chunks of <= $CHUNK_GB GB"
# chunks of at most CHUNK_GB, each transferred, checked and deleted on the ifarm before the next, so space frees as it goes
awk -v max=$((CHUNK_GB * 1000000000)) '{if (s > 0 && s + $3 > max) {c++; s = 0} s += $3; print > (dir "/chunk_" c ".txt")}' dir=$tmp c=0 $tmp/all.txt
for ch in $(ls $tmp/chunk_*.txt | sort -V); do
  : > $tmp/batch.txt
  while read -r root rep size; do
    echo "$SRC/ROOTfiles/LAD_COIN/PRODUCTION/$root $DEST/$root" >> $tmp/batch.txt
    echo "$SRC/REPORT_OUTPUT/LAD_COIN/PRODUCTION/$rep $DEST/reports/$rep" >> $tmp/batch.txt
  done < $ch
  task=$(globus transfer "$JLAB" "$SUBMIT" --batch $tmp/batch.txt --sync-level checksum --verify-checksum \
         --label "lad_replay_2026 $(date +%m%d_%H%M)" --jmespath task_id -F unix)
  echo "task $task: $(wc -l < $ch) files, $(awk '{s+=$3} END{printf "%.0f GB", s/1e9}' $ch)"
  globus task wait "$task" --timeout 86400 --polling-interval 60
  [ "$(globus task show "$task" --jmespath status -F unix)" = SUCCEEDED ] || { echo "task $task did not succeed" >&2; exit 1; }
  # check sizes, then delete on the ifarm
  : > $tmp/del.txt
  while read -r root rep size; do
    have=$(stat -c %s "$DEST/$root" 2>/dev/null || echo 0)
    if [ "$size" = "$have" ]; then echo "$SRC/ROOTfiles/LAD_COIN/PRODUCTION/$root" >> $tmp/del.txt; else echo "size mismatch $root: $size vs $have" >&2; fi
  done < $ch
  if [ -s $tmp/del.txt ]; then
    dtask=$(globus delete "$JLAB" --batch $tmp/del.txt --label "lad_replay_2026 cleanup" --jmespath task_id -F unix)
    globus task wait "$dtask" --timeout 7200 --polling-interval 30
    echo "deleted $(wc -l < $tmp/del.txt) files on the ifarm (task $dtask)"
  fi
done
rm -rf $tmp
