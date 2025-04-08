#!/bin/bash
#This script can be used to automatically run the plotW script on any number of events or runs

#Define the list of run numbers and the number of events
run_numbers=(22049 22050 22052 22098 22099 22100 22101 22102 22103 22106 22107 22108 22109 22110 22110)

nEvts=15000 

cd /home/cdaq/lad-2024/lad_replay
# Loop through each run number
# Run hcana for all run numbers in parallel
for runnum in "${run_numbers[@]}"; do
  hcana -l -q "SCRIPTS/LAD_COIN/PRODUCTION/replay_production_lad_spec.C(${runnum},${nEvts},1)" &
done
# Wait for all background processes to finish
wait

# Create the directory if it doesn't exist
mkdir -p CALIBRATION/delta_scan/15k_autoruns
# After processing all run numbers, call the root script for each run number and both P and H
for runnum in "${run_numbers[@]}"; do
  root -q "CALIBRATION/delta_scan/plot_W.C(\"ROOTfiles/CALIB/LAD_COIN_${runnum}_${nEvts}.root\",${runnum},\"P\")"
  root -q "CALIBRATION/delta_scan/plot_W.C(\"ROOTfiles/CALIB/LAD_COIN_${runnum}_${nEvts}.root\",${runnum},\"H\")"
  # Find the two most recent files created in the last 5 minutes
  recent_files=$(find CALIBRATION/delta_scan/ -maxdepth 1 -name "*${runnum}*.pdf" -mmin -5 -type f | sort -t '\0' -n | tail -n 2)
  if [ -n "$recent_files" ]; then
    mv $recent_files CALIBRATION/delta_scan/15k_autoruns/
  fi
done

# Merge all PDF files in 15k_autoruns into one for P and one for H
cd CALIBRATION/delta_scan/15k_autoruns

# Merge PDFs for P
if ls P*.pdf 1>/dev/null 2>&1; then
  pdfunite $(ls P*.pdf | sort -t '_' -k2,2n) all_delta_scan_plots_P.pdf
else
  echo "No P*.pdf files found to merge."
fi

# Merge PDFs for H
if ls H*.pdf 1>/dev/null 2>&1; then
  pdfunite $(ls H*.pdf | sort -t '_' -k2,2n) all_delta_scan_plots_H.pdf
else
  echo "No H*.pdf files found to merge."
fi