#!/bin/bash

currentRun=$1
if [ -z "$currentRun" ]; then
  echo "Usage: $0 <run_number>"
  exit 1
fi

base_out=/home/hydra/hydra_in/

declare -A dir_map=(
  [HMS]="HMS"
  [SHMS]="SHMS"
  [LAD]="LAD"
  [LAD_GEM]="LAD_GEM"
)

declare -A plots_by_dir

plots_by_dir[HMS]="
HMS_Hodoscope_ADC_Occ_Mult
HMS_Hodoscope_TDC_Occ_Mult
HMS_Hodoscope_Pedestals
HMS_Hodoscope_Pedestal_Monitoring
HMS_Drift_Chamber_Wire_Maps_All_Hits
HMS_Drift_Chamber_Drift_Distance
HMS_Drift_Chamber_Drift_Time
HMS_Focal_Plane
HMS_Target_Quantities
HMS_Cherenkov_Occ_Mult_Ped
HMS_Cherenkov_Pedestal_Monitoring
HMS_Cherenkov_NPE
HMS_Calorimeter_Occupancy
HMS_Calorimeter_Multiplicity
HMS_Calorimeter_Pedestals
HMS_Calorimeter_Pedestal_Monitoring
HMS_Drift_Chamber_Reference_Times
HMS_Trigger_Reference_Times
HMS_Raw_Fast_Raster_No_Track
HMS_Raw_Fast_Raster_Track_Cut
HMS_Fast_Raster_Position_No_Track
HMS_Fast_Raster_Position_Track_Cut
HMS_EPICS_BPM
HMS_PID
HMS_Trigger_Pedestal_Tracking
"

plots_by_dir[SHMS]="
SHMS_Hodoscope_ADC_Occ_Mult
SHMS_Hodoscope_TDC_Occ_Mult
SHMS_Hodoscope_Pedestals
SHMS_Hodoscope_Pedestal_Monitoring
SHMS_Drift_Chamber_Wire_Maps_All_Hits
SHMS_Drift_Chamber_Drift_Distance
SHMS_Drift_Chamber_Drift_Time
SHMS_Focal_Plane
SHMS_Target_Quantites
SHMS_Cherenkov_Occ_Mult
SHMS_Cherenkov_Pedestals
SHMS_Cherenkov_Pedestal_Monitoring
SHMS_Cherenkov_NPE
SHMS_Calorimeter_Occ_Mult
SHMS_Calorimeter_Pedestals
SHMS_Pre_Shower_Pedestal_Monitoring
SHMS_Calorimeter_Pedestal_Monitoring
SHMS_Drift_Chamber_Reference_Times
SHMS_Trigger_Reference_Times
SHMS_Raw_Fast_Raster_No_Track
SHMS_Raw_Fast_Raster_Track_Cut
SHMS_Fast_Raster_Position_No_Track
SHMS_Fast_Raster_Positions_Track_Cuts
SHMS_EPICS_BPM
SHMS_PID
SHMS_Trigger_Pedestal_Tracking
"

plots_by_dir[LAD]="
LAD_HMS_Hodoscope_ADC_Occ_Mult
LAD_HMS_Hodoscope_TDC_Occ_Mult
LAD_HMS_Hodoscope_Pedestals
LAD_HMS_Hodoscope_Pedestal_Monitoring
LAD_HMS_Hodoscope_PMT_Time
LAD_SHMS_Hodoscope_ADC_Occ_Mult
LAD_SHMS_Hodoscope_TDC_Occ_Mult
LAD_SHMS_Hodoscope_Pedestals
LAD_SHMS_Hodoscope_Pedestal_Monitoring
LAD_SHMS_Hodoscope_PMT_Time
LAD_1D_Hodoscope_Time_PageOne
LAD_1D_Hodoscope_Time_PageTwo
LAD_Hodoscope_ADC_PageOne
LAD_1D_Hodoscope_ADC_PageTwo
LAD_Full_Hit_Multiplicity
LAD_Full_Hit_Occupancy
"

plots_by_dir[LAD_GEM]="
GEM_HMS_Layers_with_hits
GEM_HMS_Strip_and_cluster_multiplicities
GEM_HMS_Strip_Max_Sum
GEM_HMS_Raw_occupancies_by_layer
GEM_HMS_Cluster_counts
GEM_HMS_Clustering_Results_PageOne
GEM_HMS_Clustering_Results_PageTwo
GEM_HMS_Max_Cluster
GEM_HMS_Cluster_Sum
GEM_HMS_Hit_Multiplicity
GEM_HMS_Hit_Map
GEM_HMS_Tracks
GEM_HMS_Strip_hitmap
GEM_SHMS_Layers_with_hits
GEM_SHMS_Strip_and_cluster_multiplicities
GEM_SHMS_Strip_Max_Sum
GEM_SHMS_Raw_occupancies_by_layer
GEM_SHMS_Cluster_counts
GEM_SHMS_Clustering_Results_PageOne
GEM_SHMS_Clustering_Results_PageTwo
GEM_SHMS_Max_Cluster
GEM_SHMS_Cluster_Sum_Clust
GEM_SHMS_Hit_Multiplicity
GEM_SHMS_Hit_Map
GEM_SHMS_Tracks
GEM_SHMS_Strip_hitmap
"

# Process each directory
for dir in "${!plots_by_dir[@]}"; do
  input_dir="../HISTOGRAMS/LAD_COIN/${currentRun}/${dir}"
  output_dir="${base_out}/${currentRun}"

  if [ ! -d "$output_dir" ]; then
    echo "Creating $output_dir"
    mkdir -p "$output_dir"
    # The output_dir has g+ws permission so hydra can clean dir.
    chmod 2775 "$output_dir"
  fi

  echo "Copying plots from $input_dir to $output_dir"

  i=1
  while IFS= read -r plot_name; do
    [[ -z "$plot_name" ]] && continue  # skip empty lines

    index=$(printf "%02d" $i)
    src="${input_dir}/tmp_${index}.png"
    dst="${output_dir}/${plot_name}.png"

    if [ -f "$src" ]; then
      cp "$src" "$dst"
      echo "Copied $src -> $dst"
    else
      echo "Warning: Missing file $src"
    fi

    i=$((i + 1))
  done <<< "${plots_by_dir[$dir]}"
done
