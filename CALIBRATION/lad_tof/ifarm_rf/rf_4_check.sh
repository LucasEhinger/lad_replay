#!/bin/bash
# Step 4: per-run RF/photon histograms from the replay output (one Slurm job).
set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/rf_env.sh"
cd "$RF/CALIBRATION/lad_tof"

ls "$OUT"/ROOTfiles/LAD_COIN/PRODUCTION/LAD_COIN_*_0_0_-1.root > rf_seg0_rootfiles.dat
echo "$(wc -l < rf_seg0_rootfiles.dat) replay files (expected $(wc -l < "$IFARM_RF/rf_seg0_runlist.dat"))"
sbatch -A hallc -p production -c 8 --mem=16G -t 24:00:00 -o rf_check_seg0.log --job-name rf_check \
  --wrap "source /etc/profile; module use /cvmfs/oasis.opensciencegrid.org/jlab/scicomp/sw/el9/modulefiles; \
module load root/6.30.04-gcc11.4.0; \
root -l -b -q 'lad_rf_offset_check.C+(\"rf_seg0_rootfiles.dat\",\"rf_check_seg0.root\",8)'"
echo "When the job finishes, copy $RF/CALIBRATION/lad_tof/rf_check_seg0.root to subMIT:"
echo "  /ceph/submit/data/user/e/ehingerl/hallc/lad_tof_scratch/rf_check/"
