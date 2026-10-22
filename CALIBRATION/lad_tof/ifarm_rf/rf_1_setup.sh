#!/bin/bash
# Step 1: output links, swif tarball, LADlib with the two patches. Run once, interactively.
set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/rf_env.sh"

echo "== output links -> $OUT"
mkdir -p "$OUT/ROOTfiles/LAD_COIN/PRODUCTION" "$OUT/REPORT_OUTPUT/LAD_COIN/PRODUCTION"
for d in ROOTfiles REPORT_OUTPUT; do
  [ -e "$RF/$d" ] || [ -L "$RF/$d" ] || ln -s "$OUT/$d" "$RF/$d"
  echo "   $RF/$d -> $(readlink "$RF/$d")"
done

echo "== swif tarball $TARBALL"
mkdir -p "$(dirname "$TARBALL")"
tar -czf "$TARBALL" -C "$RF" --exclude='./CALIBRATION/lad_tof/ifarm_rf/rf_seg0_*' .
ls -lh "$TARBALL"

echo "== LADlib in $LADLIB"
cd "$LADLIB"
if [ -d .git/rebase-apply ]; then
  echo "LADlib has an unfinished git am (.git/rebase-apply); run 'git am --quit' there after checking git status" >&2
  exit 1
fi
# no pipe into grep -q here: with pipefail, git log can die of SIGPIPE and the check fails
if [[ "$(git log --format=%s -20)" == *"LADKINE: vertex z in the ToF path length"* ]]; then
  echo "   patches already applied"
else
  if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
    echo "LADlib has uncommitted changes; commit or stash them first" >&2
    exit 1
  fi
  git fetch
  git checkout lad-1Dcluster-tracking-upstream
  git pull --ff-only
  git merge-base --is-ancestor "$LADLIB_BASE" HEAD || { echo "LADlib is not at/after $LADLIB_BASE" >&2; exit 1; }
  git am "$IFARM_RF"/LADlib_patches/*.patch
fi
git log --oneline -3

echo "== build (environment from $HCSWIF/setup.sh)"
set +eu # setup.sh has commands that return non-zero; it is not written for set -e
source "$HCSWIF/setup.sh" >/dev/null
set -eu
command -v hcana >/dev/null || { echo "hcana not found after sourcing $HCSWIF/setup.sh" >&2; exit 1; }
cd "$LADLIB"
./build.sh -j 8
ls -l install/lib64/libLAD.so
echo "Step 1 done. Next: rf_2_stage.sh"
