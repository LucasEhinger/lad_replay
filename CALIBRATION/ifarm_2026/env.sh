# Paths shared by the ifarm_2026 scripts (sourced, not run). Edit here if your layout differs.
SW=/work/hallc/c-lad/$USER/software
REPLAY=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)        # this lad_replay checkout (branch lad-replay-2026)
IF=$REPLAY/CALIBRATION/ifarm_2026
OUT=/volatile/hallc/c-lad/$USER/lad_replay_2026                   # replay output (ROOTfiles, REPORT_OUTPUT)
HCSWIF=$SW/hcswif_LAD
LADLIB=$HOME/hallc/software/LADlib                                # the build hcswif_LAD/setup.sh loads
LADLIB_BRANCH=lad-replay-2026
EMAIL=ehingerl@mit.edu
TARBALL_FILE=$IF/tarball_path.txt                                 # written by 1_setup.sh
