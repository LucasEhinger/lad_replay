#!/bin/bash
# Step 3: check every raw file is on /cache, then build and start the swif2 workflow.
set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/rf_env.sh"
cd "$IFARM_RF"

missing=0
while read -r f; do
  [ -f "${f/#\/mss/\/cache}" ] || missing=$((missing + 1))
done < rf_seg0_files.txt
if [ "$missing" -gt 0 ]; then
  echo "$missing of $(wc -l < rf_seg0_files.txt) raw files are not on /cache yet; try again later." >&2
  exit 1
fi
echo "All $(wc -l < rf_seg0_files.txt) raw files are on /cache."

cd "$HCSWIF"
./hcswif.py --mode replay --spectrometer LAD_COIN --account hallc --all_segs true --time 86400 \
  --run file "$IFARM_RF/rf_seg0_runlist.dat" --specify_replay "$TARBALL" --name "$WF"
swif2 import -file "jsons/$WF.json"
swif2 run "$WF"
swif2 notify -workflow "$WF" -email "$EMAIL" -when done
echo "Submitted $WF. Progress: swif2 status $WF. When it is done, run rf_4_check.sh."
