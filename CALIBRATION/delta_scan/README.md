HMS/SHMS Delta Scan Calibration
============================================
Storing scripts here for running the delta scan. Scripts don't do much other than make plots.




Running code
---------------
NOTE: P: SHMS,  H: HMS, 

From lad_replay/, do
`root -q 'CALIBRATION/delta_scan/plot_W.C("ROOTfiles/CALIB/LAD_COIN_22050_50000.root",22050,"P")'`

Where you specify the input rootfile, run number (used to grab run info from standard_kinematics), and spectrometer.
Standard database and file save location (CALIBRATION/delta_scan) are hard coded