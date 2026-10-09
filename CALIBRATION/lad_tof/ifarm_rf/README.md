# RF offset study: replays on the JLab ifarm

Goal: per-run `l_rf_offset` for every beam run, from the RF - vertex time and from the
LAD photon flash. Step A (first segment of every run) gives the RF - vertex phase of every
run and a photon-flash check of every jump. Step B (all segments) gives per-run
photon-flash values at the ~0.05 ns level; it is the full production replay.

Files here:
- `rf_beam_runs.txt`: the 1447 beam runs in `CALIBRATION/files/LAD runs.xlsx`, without cosmics,
  junk runs and rows without a target.
- `rf_env.sh`, `rf_1_setup.sh` ... `rf_4_check.sh`: the steps below; `make_rf_lists.sh` writes the jcache and hcswif lists from the tape stubs.
- `LADlib_patches/`: LADlib commits not yet on GitHub (per-spectrometer ToF offset, vertex z
  in the ToF). The replay must be run with them.
- `../lad_rf_offset_check.C` + `../lad_tof_offset.h`: per-run histograms from the replay output.

## How to run

On an ifarm node:

```bash
cd /work/hallc/c-lad/$USER/software
mkdir lad_replay_rf && tar -xzf lad_replay_rf_2026-10-03.tar.gz -C lad_replay_rf
cd lad_replay_rf/CALIBRATION/lad_tof/ifarm_rf

./rf_1_setup.sh     # output links to /volatile, swif tarball, LADlib patches + build (once)
./rf_2_stage.sh     # run lists + one jcache request for all first segments
#   ... wait for jcache's email ...
./rf_3_submit.sh    # refuses until every file is on /cache; then hcswif.py + swif2 import/run/notify
#   ... wait for swif2's email (swif2 status lad_rf_seg0) ...
./rf_4_check.sh     # Slurm job: lad_rf_offset_check.C -> rf_check_seg0.root, copy that to subMIT
```

Paths and names (output directory, tarball, workflow name, e-mail) are in `rf_env.sh`.
Each script stops at the first error. Replay output:
`/volatile/hallc/c-lad/$USER/lad_replay_rf/ROOTfiles/LAD_COIN/PRODUCTION/LAD_COIN_<run>_0_0_-1.root`.

Step B (all segments) would use `./make_rf_lists.sh rf_beam_runs.txt rf_all 5` (jobs of 5
segments) with its own workflow name; it is the full production replay (several hundred TB).
