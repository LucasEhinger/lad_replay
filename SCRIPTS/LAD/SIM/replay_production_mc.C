// #include "../LAD_link_defs.h"
void replay_production_mc(Int_t RunNumber = 0, Int_t MaxEvent = 0) {

  // Get RunNumber and MaxEvent if not provided.
  if (RunNumber == 0) {
    cout << "Enter a Run Number (-1 to exit): ";
    cin >> RunNumber;
    if (RunNumber <= 0)
      return;
  }
  if (MaxEvent == 0) {
    cout << "\nNumber of Events to analyze: ";
    cin >> MaxEvent;
    if (MaxEvent == 0) {
      cerr << "...Invalid entry\n";
      exit;
    }
  }

  // Load Global parameters
  // Add variables to global list.
  gHcParms->Define("gen_run_number", "Run Number", RunNumber);
  gHcParms->AddString("g_ctp_database_filename", "DBASE/LAD/standard.database");
  gHcParms->Load(gHcParms->GetString("g_ctp_database_filename"), RunNumber);
  gHcParms->Load(gHcParms->GetString("g_ctp_parm_filename"));
  gHcParms->Load(gHcParms->GetString("g_ctp_kinematics_filename"), RunNumber);

  // Load the Hall C detector map
  gHcDetectorMap = new THcDetectorMap();
  gHcDetectorMap->Load("MAPS/LAD/DETEC/HODO/lhodo_mc.map");

  // Add LAD detector
  THcLADSpectrometer *LAD = new THcLADSpectrometer("L", "LAD");
  gHaApps->Add(LAD);

  THcLADHodoscope *lhod = new THcLADHodoscope("hod", "LAD Hodoscope");
  LAD->AddDetector(lhod);

  THcLADGEM *gem = new THcLADGEM("gem", "gem");
  LAD->AddDetector(gem);

  THcAnalyzer *analyzer = new THcAnalyzer;

  // Set the decoder to use simulation input
  THcInterface::SetDecoder(LADSimDecoder::Class());

  THaEvent *event = new THaEvent;

  // TString run_file = "../libLADdig/test_scripts/lad_hodo_gem_sim.root";
  TString run_file = Form("/volatile/hallc/c-lad/ehingerl/G4_LAD/carlos_proton/dig/ScanLAD_proton_%dMeV_10k_20240205_dig.root", RunNumber);
  THaRunBase *run = new LADSimFile(run_file.Data(), "lad", "");

  // Set to read in Hall C run database parameters
  // run->SetRunParamClass("THcRunParameters");

  run->SetEventRange(1, MaxEvent); // Physics Event number, does not include scaler or control events.
  // run->SetNscan(1);
  run->SetDataRequired(0);
  // run->SetDate("2022-10-06 00:00:00");
  run->Print();

  // Define the analysis parameters
  analyzer->SetCountMode(2); // 0 = counter is # of physics triggers
                             // 1 = counter is # of all decode reads
                             // 2 = counter is event number

  analyzer->SetEvent(event);


  // Define crate map
  analyzer->SetCrateMapFileName("MAPS/db_cratemap_mc.dat");
  // Define output ROOT file

  analyzer->SetOutFile(Form("ROOTfiles/full_mc_%dMeV_%d.root", RunNumber, MaxEvent));
  // // Define output DEF-file
  analyzer->SetOdefFile("DEF-files/LAD/PRODUCTION/lstackana_production_all.def");
  // // Define cuts file
  analyzer->SetCutFile("DEF-files/LAD/PRODUCTION/CUTS/lstackana_production_cuts.def"); // optional
  // File to record cuts accounting information for cuts
  // analyzer->SetSummaryFile(Form("REPORT_OUTPUT/HMS/PRODUCTION/summary_production_%d_%d.report", RunNumber,
  // MaxEvent));    // optional Start the actual analysis.
  analyzer->Process(run);
  // Create report file from template.
  // analyzer->PrintReport("TEMPLATES/LAD/PRODUCTION/lstackana_production.template",
  //                       Form("REPORT_OUTPUT/HMS/PRODUCTION/replay_hms_production_%d_%d.report", RunNumber,
  //                       MaxEvent));
}