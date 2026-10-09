#!/bin/bash
# Run on the JLab ifarm. From a list of run numbers, write
#   <out>_files.txt    raw files to stage with jcache (/mss paths)
#   <out>_runlist.dat  hcswif.py run list: "run seg_end seg_start run_type 0"
#   <out>_missing.txt  runs with no raw file on tape
# The raw file prefix (run type) and the last segment are read from the tape stubs in
# /mss/hallc/c-lad/raw, so the files do not need to be cached first.
#
# Usage: ./make_lists.sh <runs.txt> <out_prefix> [segs_per_job]
#   segs_per_job = 0 (default): first segment (.dat.0) of each run only.
#   segs_per_job = N > 0:       every segment, in jobs of N segments each.
#
# run_type indices are those of hcswif.py / make_lad_runlist.sh.

RAW=${RAW_DIR:-/mss/hallc/c-lad/raw}
run_types=("lad_Production" "lad_Production_noGEM" "lad_LADwGEMwROC2" "lad_GEMonly" "lad_LADonly"
           "lad_SHMS_HMS" "lad_SHMS" "lad_HMS")

runs_file=$1
out=$2
nseg=${3:-0}
if [ -z "$runs_file" ] || [ -z "$out" ]; then
  echo "Usage: $0 <runs.txt> <out_prefix> [segs_per_job]" >&2
  exit 1
fi

: >"${out}_files.txt"
: >"${out}_runlist.dat"
: >"${out}_missing.txt"
for run in $(grep -v "^#" "$runs_file" | sort -un); do
  type=-1
  for i in "${!run_types[@]}"; do
    if ls "$RAW/${run_types[$i]}_${run}.dat."* >/dev/null 2>&1; then
      type=$i
      break
    fi
  done
  if [ "$type" -lt 0 ]; then
    echo "$run" >>"${out}_missing.txt"
    continue
  fi
  stem="$RAW/${run_types[$type]}_${run}"
  last=$(ls "$stem".dat.* | sed -n 's/.*\.dat\.\([0-9][0-9]*\)$/\1/p' | sort -n | tail -1)
  if [ "$nseg" -le 0 ]; then
    echo "$stem.dat.0" >>"${out}_files.txt"
    echo "$run 0 0 $type 0" >>"${out}_runlist.dat"
  else
    for s in $(seq 0 "$last"); do echo "$stem.dat.$s"; done >>"${out}_files.txt"
    for ((start = 0; start <= last; start += nseg)); do
      end=$((start + nseg - 1))
      [ "$end" -gt "$last" ] && end=$last
      echo "$run $end $start $type 0" >>"${out}_runlist.dat"
    done
  fi
done
echo "$(wc -l <"${out}_runlist.dat") jobs, $(wc -l <"${out}_files.txt") raw files," \
     "$(wc -l <"${out}_missing.txt") runs without raw data (${out}_missing.txt)"
