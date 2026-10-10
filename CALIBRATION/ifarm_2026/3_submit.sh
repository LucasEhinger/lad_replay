#!/bin/bash
# Step 3 (ifarm): check that every raw file of a list is on /cache, then build and start the swif2 workflow.
#   ./3_submit.sh allseg lad_2026_allseg
set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/env.sh"
cd "$IF"
out=$1; wf=$2
TARBALL=$(cat "$TARBALL_FILE")
[ -f "$TARBALL" ] || { echo "tarball $TARBALL missing; run 1_setup.sh" >&2; exit 1; }
missing=0
while read -r f; do [ -f "${f/#\/mss/\/cache}" ] || missing=$((missing + 1)); done < "${out}_files.txt"
if [ "$missing" -gt 0 ]; then
  echo "$missing of $(wc -l < "${out}_files.txt") raw files are not on /cache yet; try again later." >&2
  exit 1
fi
cd "$HCSWIF"
./hcswif.py --mode replay --spectrometer LAD_COIN --account hallc --all_segs true --time 86400 \
  --run file "$IF/${out}_runlist.dat" --specify_replay "$TARBALL" --name "$wf"
swif2 import -file "jsons/$wf.json"
# swif2 leaves a new workflow undispatched until `swif2 run` is called twice, and a pair issued right after the
# import is sometimes not enough (2026-10-04, 2026-10-09): repeat the pair until jobs are dispatched.
for i in 1 2 3 4 5; do
  swif2 run "$wf"; swif2 run "$wf"
  [ "$(swif2 status "$wf" | awk '$1=="undispatched"{print $3}')" = 0 ] && break
  sleep 60
done
swif2 status "$wf" | grep -E '^(undispatched|dispatched) '
swif2 notify -workflow "$wf" -email "$EMAIL" -when done
echo "Submitted $wf ($(grep -vc '^#' "$IF/${out}_runlist.dat") jobs). Progress: swif2 status $wf"
