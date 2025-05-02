#! /bin/bash

# Which spectrometer are we analyzing.
spec="lad_coin"
SPEC="LAD_COIN"

# What is the last run number for the spectrometer.
# The pre-fix zero must be stripped because ROOT is ... well ROOT
# lastRun=$( \
#     ls raw/"${spec}"_all_*.dat raw/../raw.copiedtotape/"${spec}"_all_*.dat -R 2>/dev/null | perl -ne 'if(/0*(\d+)/) {print "$1\n"}' | sort -n | tail -1 \
# )


runNum=$1
if ! [[ "$runNum" =~ ^[0-9]+$ ]]; then
  runNum=$(find ./raw/ -type f -name "*.dat.*" \
  | grep -Eo '[0-9]{5}\.dat\.' \
  | grep -Eo '[0-9]{5}' \
  | sort -n | tail -1)
fi

#Currently not pointing to cache. Can change later if necessary 
# /cache/hallc/c-lad/raw/lad_Production_*.dat.* (or just ./cache soft link)
# Find the most recent run number

# Define an array of directories to search
# Define an array of directories to search
directories=(./raw/ ./cache/) #  ./cache/ currently causing problems

for dir in "${directories[@]}"; do
  lastRunFile=$(ls "${dir}"/*"${runNum}.dat.0" 2>/dev/null | head -n 1)
  if [[ -n $lastRunFile ]]; then
    break
  fi
done


#/volatile/hallc/c-lad/ehingerl/raw_data/LAD_cosmic/
# Determine the run_type based on the lastRun value
case "$lastRunFile" in
  *Production_noGEM*) run_type=1 ;;
  *Production*) run_type=0 ;;
  *LADwGEMwROC2*) run_type=2 ;;
  *GEMonly*) run_type=3 ;;
  *LADonly*) run_type=4 ;;
  *SHMS_HMS*) run_type=5 ;;
  *SHMS*) run_type=6 ;;
  *HMS*) run_type=7 ;;
  *) run_type=-1 ;; # Default case if no match is found
esac

# Which run to analyze.

numEvents=$2

if [ -z "$numEvents" ]; then
  numEvents=50000
fi

numEventsk=$((numEvents / 1000))
# Which scripts to run.
script="SCRIPTS/${SPEC}/PRODUCTION/replay_production_lad_spec.C"
#config="CONFIG/${SPEC}/PRODUCTION/${spec}_production.cfg"
#expertConfig="CONFIG/${SPEC}/PRODUCTION/${spec}_production_expert.cfg"

# Define some useful directories
rootFileDir="./ROOTfiles/${SPEC}/${spec}50k"
monRootDir="./HISTOGRAMS/${SPEC}/ROOT"
monPdfDir="./HISTOGRAMS/${SPEC}/PDF"
reportFileDir="./REPORT_OUTPUT/${SPEC}"
reportMonDir="./UTIL_OL/REP_MON"
reportMonOutDir="./MON_OUTPUT/${SPEC}/REPORT"

# Name of the report monitoring file
reportMonFile="summary_output_${runNum}.txt"

# Which commands to run.
runHcana="hcana -q \"${script}(${runNum}, ${numEvents},${run_type})\""
#runHcana="/home/cdaq/cafe-2022/hcana/hcana -q \"${script}(${runNum}, ${numEvents})\""
# The original onlineGUI commands (single configuration):
# runOnlineGUI="panguin -f ${config} -r ${runNum}"
# saveOnlineGUI="panguin -f ${config} -r ${runNum} -P"
# saveExpertOnlineGUI="./online -f ${expertConfig} -r ${runNum} -P"
# Replace these with a multi-config loop below.

runReportMon="./${reportMonDir}/reportSummary.py ${runNum} ${numEvents} ${spec} singles"
openReportMon="emacs ${reportMonOutDir}/${reportMonFile}"

# Name of the replay ROOT file
replayFile="${spec}_replay_production_${runNum}"
rootFile="${replayFile}_${numEvents}.root"
latestRootFile="${rootFileDir}/${replayFile}_latest.root"

# Names of the monitoring file
monRootFile="${spec}_production_${runNum}.root"
monPdfFile="${spec}_production_${runNum}.pdf"
#monExpertPdfFile="${spec}_production_expert_${runNum}.pdf"
latestMonRootFile="${monRootDir}/${spec}_production_latest.root"
latestMonPdfFile="${monPdfDir}/${spec}_production_latest.pdf"

# Where to put log
reportFile="${reportFileDir}/replay_${spec}_production_${runNum}_${numEvents}.report"
summaryFile="${reportFileDir}/summary_production_${runNum}_${numEvents}.txt"

# What is base name of panguin output.
outFile="${spec}_production_${runNum}"
#outExpertFile="${spec}_production_expert_${runNum}"

outFileMonitor="output.txt"

# Replay out files
replayReport="${reportFileDir}/REPLAY_REPORT/replayReport_${spec}_production_${runNum}_${numEvents}.txt"

#############################################
# Define arrays for multiple GUI configurations.
# Each index corresponds to:
#   1. A tag for output naming.
#   2. The normal GUI configuration file.
#   3. The expert GUI configuration file.
# gui_tags=("lad_coin" "lad_kin" "shms" "hms")
gui_tags=("lad_gem" "lad_coin" "shms" "hms")
# gui_tags=("shms" "hms")

gui_configs=(
  "CONFIG/LAD/PRODUCTION/lad_gem.cfg"
  "CONFIG/LAD/PRODUCTION/lad_coin_production.cfg"
  # "CONFIG/LAD/PRODUCTION/lad_coin_kin.cfg"
  "CONFIG/SHMS/PRODUCTION/shms_production.cfg"
  "CONFIG/HMS/PRODUCTION/hms_production.cfg"
)

expert_configs=(
  "CONFIG/LAD/PRODUCTION/lad_gem.cfg"
  "CONFIG/LAD/PRODUCTION/lad_coin_production.cfg"
  # "CONFIG/LAD/PRODUCTION/lad_coin_kin.cfg"
  "CONFIG/SHMS/PRODUCTION/shms_production.cfg"
  "CONFIG/HMS/PRODUCTION/hms_production.cfg"
)
#############################################

# Start analysis and monitoring plots.
{
  echo ""
  echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="
  echo ""
  date
  echo ""
  echo "Running ${SPEC} analysis on the run ${runNum}:"
  echo " -> SCRIPT:  ${script}"
  echo " -> RUN:     ${runNum}"
  echo " -> NEVENTS: ${numEvents}"
  echo " -> COMMAND: ${runHcana}"
  echo ""
  echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="

  sleep 2
  eval ${runHcana}

  # Link the ROOT file to latest for online monitoring
  #Need to match ${rootFile} to the output name format of the coin_replay
  #ln -fs ${rootFile} ${latestRootFile} #
  ln -fs "../../CALIB/${SPEC}_${runNum}_${numEvents}.root" ${latestRootFile}
  #ln -fs "../../LADC_COIN/PRODUCTION/${SPEC}_production_hall_${runNum}_${numEvents}.root" ${latestRootFile}
  #ln -fs "../../LADC_COIN/CALIBRATION/${SPEC}_calibration_hall_${runNum}_${numEvents}.root" ${latestRootFile}
  sleep 2
  function yes_or_no() {
    while true; do
      read -p "$* [y/n]: " yn
      case $yn in
      [Yy]*) return 0 ;;
      [Nn]*)
        echo "No entered"
        return 1
        ;;
      esac
    done
  }

  echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="
  echo ""
  # echo "Do you want to generate expert LAD detector (gem + hodoscope) plots?"
  # yes_or_no "Do you want to generate expert LAD detector (gem + hodoscope) plots?" &&
  # {
    echo "Generating expert LAD detector plots."
    # This script runs a ROOT macro to process the latest ROOT file and generate histograms.
    # The macro "lad_histos.C" is executed with two arguments:
    # - The first argument (${latestRootFile}) specifies the latest ROOT file to process.
    # - The second argument (0) indicates that the histograms are generated for both HMS and SHMS LAD.
    root -l -b -q "macros/LAD/lad_histos_MT.C(\"${latestRootFile}\",0,${numEvents})" 
    # Currently on generating for 1k events. Will have to come up with a faster way to make these histograms.
  # }



  echo ""
  echo ""
  echo ""
  echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="
  echo ""
  echo "Running onlineGUI for analyzed ${SPEC} run ${runNum}:"
  echo " -> Using multiple configuration files"
  echo " -> RUN:     ${runNum}"
  echo "-------------------------------------------------------------"
  
  sleep 1
  # Loop over each GUI configuration.
  # Initialize an empty array to hold the PDF filenames.
  mergedPDFs=()
  for i in "${!gui_tags[@]}"; do
    tag="${gui_tags[$i]}"
    config="${gui_configs[$i]}"
    expertConfig="${expert_configs[$i]}"
    # Define output name including the tag.
    outFile="${spec}_production_${runNum}_${tag}"
    outExpertFile="summaryPlots_${runNum}_${expertConfig##*/}"
    outExpertFile="${outExpertFile%.cfg}"

    echo ""
    echo "Processing GUI for configuration: ${tag}"
    echo " -> CONFIG:  ${config}"
    echo " -> EXPERT CONFIG: ${expertConfig}"
    echo " -> RUN:     ${runNum}"
    echo "-------------------------------------------------------------"

    sleep 2
    cd onlineGUI || exit 1
    yes_or_no "Do you want to view plots for ${tag}?" &&
    {
      # Run the normal GUI command.
      panguin -f "${config}" -r "${runNum}"
    }
    # Run the expert GUI command (-P flag).
    panguin -f "${expertConfig}" -r "${runNum}" -P

    # Display current directory and output file info.
    pwd
    echo " -> outExpertFile: ${outExpertFile}"
    echo "../HISTOGRAMS/${SPEC}/PDF/${outFile}.pdf"

    # Move the resulting expert PDF to the appropriate directory with the tag in its name.
    #monExpertPdfFile="../HISTOGRAMS/${SPEC}/PDF/${outFile}_expert.pdf"
    echo "Moving Expert PDF to ${monExpertPdfFile}"
    monExpertPdfFile="$(readlink -f "../HISTOGRAMS/${SPEC}/PDF/${outFile}_expert.pdf")"
    mv "${outExpertFile}.pdf" ${monExpertPdfFile}

    mergedPDFs+=("${monExpertPdfFile}")
    cd .. || exit 1
  done
  #merge the pdfs from difference specs
  echo "Merging PDF file from all specs ${mergedPDFs[@]}"
  # Initialize an array for valid PDFs (with separators)
  validPDFs=()

  # Loop over each generated expert PDF in your mergedPDFs array
  for pdf in "${mergedPDFs[@]}"; do
    # Check if the PDF is valid
    if pdfinfo "$pdf" >/dev/null 2>&1; then
      # Create a separator page with a title based on the PDF name.
      label=$(basename "$pdf" _expert.pdf)
      # Create a temporary separator page PDF.
      separator="./HISTOGRAMS/LAD_COIN/PDF/tmp/separator_${label}.pdf"
      convert -size 1600x800 xc:white -fill black -gravity center -pointsize 36 "label:${label}" "$separator"
      # Append the separator and the valid PDF to the list.
      validPDFs+=("$separator")
      validPDFs+=("$pdf")
    else
      echo "Skipping damaged PDF: $pdf"
    fi
  done

  # Now merge all valid PDFs (which include separator pages) into one file.
  pdfunite "${validPDFs[@]}" "${latestMonPdfFile}"
  rm ./HISTOGRAMS/LAD_COIN/PDF/tmp/separator_*.pdf
  #pdfunite "${mergedPDFs[@]}" "${latestMonPdfFile}"
  # Update the latest PDF link (using the expert PDF of the last configuration processed).
  #ln -fs ${latestMonPdFile} ${latestMonPdfFile}

  ###########################################################
  # function used to prompt user for questions
  # post pdfs in hclog
   yes_or_no "Upload these plots to logbook HCLOG? " && {
    read -p "Enter a text body for the log entry (or leave blank): " logCaption
    echo "$logCaption" > caption.txt
    /site/ace/certified/apps/bin/logentry \
      -cert /home/cdaq/.elogcert \
      -t "${numEventsk}k replay plots for run ${runNum}" \
      -e cdaq \
      -l HCLOG \
      -a ${latestMonPdfFile} \
      -b "caption.txt"

      rm -rf "caption.txt"
  }

  #    /home/cdaq/bin/hclog \
  #    --logbook "HCLOG" \
  #    --tag Autolog \
  #    --title ${events}" replay plots for run ${runnum} TEST TEST TEST" \
  #    --attach ${latestMonPdfFile}
  ###########################################################
  #    --title ${events}" replay plots for run ${runnum}" \

  echo ""
  echo ""
  echo ""
  echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="
  echo ""
  echo "Done analyzing ${SPEC} run ${runNum}."
  echo ""
  echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="

  sleep 2

  echo ""
  echo ""
  echo ""
  echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="
  echo ""
  echo "(Not) Generating report file monitoring data file ${SPEC} run ${runNum}."
  echo ""
  echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="

  #  eval ${runReportMon}
  #  mv "${outFileMonitor}" "${reportMonOutDir}/${reportMonFile}"
  #  eval ${openReportMon}

  # sleep 2

  # echo ""
  # echo ""
  # echo ""
  # echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="
  # echo ""
  # echo "Taking hextuple ratios in  ${latestRootFile}"
  # echo ""
  # echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="

  #if [[ "${spec}" == "shms" ]]; then
  #  ./hcana -l << EOF
  #  .L ${shmsCounter}
  #  run_el_counter_${spec}("${latestRootFile}");
  #EOF
  #fi
  #if [[ "${spec}" == "hms" ]]; then
  #  ./hcana -l << EOF
  #  .L ${hmsCounter}
  #  run_el_counter_${spec}("${latestRootFile}");
  #EOF
  #fi

  sleep 2

  echo ""
  echo ""
  echo ""
  echo "-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|"
  echo ""
  echo "Done!"
  echo ""
  echo "-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|"
  echo ""
  echo "Useful links:"
  echo " -> ROOT file: ${rootFileDir}/${rootFile}"
  echo " -> Histogram (ROOT) file: ${monRootDir}/${monRootFile}"
  echo " -> Histogram (PDF) file: ${monPdfDir}/${monPdfFile}"
  echo " -> Raw EVIO file: ${lastRunFile}"
  echo ""

} 2>&1 | tee "${replayReport}"
