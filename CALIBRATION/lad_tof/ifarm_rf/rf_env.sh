# Paths shared by the rf_*.sh scripts (sourced, not run). Edit here if your layout differs.
SW=/work/hallc/c-lad/$USER/software
RF=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)      # the unpacked lad_replay_rf
IFARM_RF=$RF/CALIBRATION/lad_tof/ifarm_rf
OUT=/volatile/hallc/c-lad/$USER/lad_replay_rf                   # replay output (ROOTfiles, REPORT_OUTPUT)
TARBALL=$SW/lad_replay_versions/lad_replay_rf_2026-10-03.tar.gz  # what the swif jobs unpack
HCSWIF=$SW/hcswif_LAD
LADLIB=$HOME/hallc/software/LADlib                               # the build hcswif_LAD/setup.sh loads
WF=lad_rf_seg0                                                   # swif2 workflow name
EMAIL=ehingerl@mit.edu
LADLIB_BASE=49d4e2f                                              # origin commit the patches apply to
