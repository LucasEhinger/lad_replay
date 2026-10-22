#!/bin/bash
# Step 2: list the first segment of every beam run and request it from tape (one jcache request).
set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/rf_env.sh"
cd "$IFARM_RF"

./make_rf_lists.sh rf_beam_runs.txt rf_seg0
echo "Runs without raw data: $(tr '\n' ' ' < rf_seg0_missing.txt)"
jcache get -e "$EMAIL" $(cat rf_seg0_files.txt)
echo "Requested $(wc -l < rf_seg0_files.txt) files. When jcache's email arrives (or rf_3_submit.sh"
echo "reports nothing missing), run rf_3_submit.sh."
