#!/bin/bash
# Step 2 (ifarm): job lists from the tape stubs and one jcache request.
#   ./2_stage.sh runs_allseg.txt allseg 1     every segment, one segment per job
#   ./2_stage.sh runs_seg0.txt   seg0   0     first segment only
set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/env.sh"
cd "$IF"
runs=$1; out=$2; nseg=${3:-1}
./make_lists.sh "$runs" "$out" "$nseg"
[ -s "${out}_missing.txt" ] && echo "Runs without raw data: $(tr '\n' ' ' < "${out}_missing.txt")"
jcache get -D 21 -e "$EMAIL" $(cat "${out}_files.txt")
echo "Requested $(wc -l < "${out}_files.txt") files (pinned 21 days). When jcache's e-mail arrives, run 3_submit.sh $out <workflow>."
