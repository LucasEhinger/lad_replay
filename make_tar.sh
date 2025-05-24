#!/bin/bash

# This will have to be changed for your specific directory. Is currently set up for Lucas's
TAR_PATH="/work/hallc/c-lad/ehingerl/software/lad_replay.tar.gz"
if [ -e "$TAR_PATH" ]; then
  i=1
  while [ -e "/work/hallc/c-lad/ehingerl/software/lad_replay.tar.gz.$i" ]; do
    i=$((i+1))
  done
  TAR_PATH="/work/hallc/c-lad/ehingerl/software/lad_replay.tar.gz.$i"
fi

cd /work/hallc/c-lad/ehingerl/software/lad_replay
tar --exclude='*.root' -czf ../lad_replay.tar.gz .
cd -