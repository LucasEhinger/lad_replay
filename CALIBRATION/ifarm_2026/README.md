# 2026 replay on the JLab ifarm

Full replay with the 2026 LAD calibration (lad_replay and LADlib branch `lad-replay-2026`):
2026 walk, per-HV-period bar timing and proton-scale energy, 10-row RF table, 1D-cluster
tracking, GEM zero suppression at 5 sigma, 100-1 disabled, GEM space points in the output.

Runs (`runs_allseg.txt`, every segment, one segment per job):
- the 5-pass carbon program (`CALIBRATION/files/run-lists/all_C3_runlist.dat`: carbon, optics, empty target),
- 3-pass carbon, optics and dummy,
- LD2 runs 22811, 23110 and 23111, next to carbon runs, for an LD2/carbon comparison.

`runs_seg0.txt`: first segments of LD2 runs in the HV periods without carbon (timing/energy checks).

## Steps

```bash
# ifarm (bash): LADlib in ~/hallc/software/LADlib on branch lad-replay-2026, built with ./build.sh -c
cd /work/hallc/c-lad/$USER/software/lad_replay_2026/CALIBRATION/ifarm_2026
./1_setup.sh                          # output links -> /volatile/.../lad_replay_2026, pinned swif tarball
./2_stage.sh runs_allseg.txt allseg 1 # job lists + jcache request (pinned 21 days)
./2_stage.sh runs_seg0.txt seg0 0
#   ... jcache e-mail ...
./3_submit.sh allseg lad_2026_allseg  # hcswif.py + swif2 import/run/notify
./3_submit.sh seg0 lad_2026_seg0
# subMIT, while the workflows run (keeps the ifarm /volatile use low):
./4_transfer.sh                       # Globus: finished replays -> subMIT, size check, delete on the ifarm
```

Output: `LAD_COIN_<run>_<seg>_<seg>_-1.root` (one file per segment) in
`/ceph/submit/data/user/e/ehingerl/hallc/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/replay_2026/`
(the group area `/ceph/submit/data/group/lad/lad_replay_2026` when the user area has less than 1 TB left),
reports in its `reports/` subdirectory.
