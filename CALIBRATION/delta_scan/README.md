HMS/SHMS Delta Scan Calibration
============================================
Storing scripts here for running the delta scan. Scripts don't do much other than make plots.




Running code
---------------
NOTE: P: SHMS,  H: HMS, 

From lad_replay/, do
`root -q 'CALIBRATION/delta_scan/plot_W.C("ROOTfiles/CALIB/LAD_COIN_22050_50000.root",22050,"P")'`

You'll also have to update the 'DBASE/LAD_COIN/standard.kinematics' file with the correct shms (p) and hms (h) momenta and angles. The LAD runs start from the top (one of the top 10 configurations). Copy the most recent config block (~8 lines), paste it, and change the run numbers and values.

Where you specify the input rootfile, run number (used to grab run info from standard_kinematics), and spectrometer.
Standard database and file save location (CALIBRATION/delta_scan) are hard coded

`run_all.sh script`
This script automatically runs the replay on a specified number of events, and then runs the W-plotting script on that. It then copies all pdfs into a new directory and merges them.