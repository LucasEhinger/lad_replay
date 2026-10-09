#!/bin/bash
# Step 1 (ifarm, once): check LADlib, make the output links and the pinned swif tarball.
set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/env.sh"

echo "== LADlib $LADLIB"
b=$(git -C "$LADLIB" rev-parse --abbrev-ref HEAD)
[ "$b" = "$LADLIB_BRANCH" ] || { echo "LADlib is on $b, not $LADLIB_BRANCH" >&2; exit 1; }
[ -z "$(git -C "$LADLIB" status --porcelain --untracked-files=no)" ] || { echo "LADlib has local changes" >&2; exit 1; }
[ "$LADLIB/install/lib64/libLAD.so" -nt "$LADLIB/.git/HEAD" ] || echo "warning: libLAD.so older than the checkout; rebuild with ./build.sh -c"
echo "   $(git -C "$LADLIB" log --oneline -1)"

echo "== lad_replay $REPLAY ($(git -C "$REPLAY" log --oneline -1))"
[ -z "$(git -C "$REPLAY" status --porcelain --untracked-files=no)" ] || { echo "lad_replay has local changes" >&2; exit 1; }

echo "== output links -> $OUT"
mkdir -p "$OUT/ROOTfiles/LAD_COIN/PRODUCTION" "$OUT/REPORT_OUTPUT/LAD_COIN/PRODUCTION"
for d in ROOTfiles REPORT_OUTPUT; do
  [ -e "$REPLAY/$d" ] || [ -L "$REPLAY/$d" ] || ln -s "$OUT/$d" "$REPLAY/$d"
  echo "   $REPLAY/$d -> $(readlink "$REPLAY/$d")"
done

sha=$(git -C "$REPLAY" rev-parse --short HEAD)
lsha=$(git -C "$LADLIB" rev-parse --short HEAD)
TARBALL=$SW/lad_replay_versions/lad_replay_2026_${sha}_LADlib_${lsha}.tar.gz
echo "== swif tarball $TARBALL"
mkdir -p "$(dirname "$TARBALL")"
tar -czf "$TARBALL" -C "$REPLAY" --exclude=./.git --exclude='./CALIBRATION/ifarm_2026/*_files.txt' \
    --exclude='./CALIBRATION/ifarm_2026/*_runlist.dat' .
echo "$TARBALL" > "$TARBALL_FILE"
ls -lh "$TARBALL"
echo "Step 1 done. Next: 2_stage.sh"
