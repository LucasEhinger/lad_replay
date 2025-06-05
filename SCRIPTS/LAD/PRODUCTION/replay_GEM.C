// Only replay's LAD GEM detector
// Currently set up to replay GEM data from the test lab
// #include "../../LAD/LAD_link_defs.h". Leave this line commented. Used for debugging purposes.
void load_GEM_CM_PED(int runNumber) {
  std::vector<int> ped_cm_runs;
  int ped_cm_runs_count = 0;

  // Open the file
  std::ifstream infile("PARAM/LAD/GEM/lgem_cm_ped_runs.param");
  if (infile) {
    std::string content((std::istreambuf_iterator<char>(infile)), std::istreambuf_iterator<char>());
    std::stringstream ss(content);
    std::string token;
    while (std::getline(ss, token, ',')) {
      std::stringstream token_ss(token);
      int value;
      if (token_ss >> value) {
        ped_cm_runs.push_back(value);
      }
    }
    ped_cm_runs_count = ped_cm_runs.size();
  } else {
    std::cerr << "Error: Could not open PARAM/LAD/GEM/lgem_cm_ped_runs.param" << std::endl;
    ped_cm_runs_count = 0;
  }
  std::sort(ped_cm_runs.begin(), ped_cm_runs.end());

  int ped_cm_file_num = ped_cm_runs[ped_cm_runs_count - 1]; // Default to the last run number in the list
  for (int i = ped_cm_runs_count - 1; i > 0; --i) {
    if (ped_cm_runs[i] < runNumber) {
      ped_cm_file_num = ped_cm_runs[i];
      break;
    }
  }

  gHcParms->AddString("lgem_pedfile", Form("PARAM/LAD/GEM/PED/gem_ped_%d.dat", ped_cm_file_num));
  gHcParms->AddString("lgem_cmfile", Form("PARAM/LAD/GEM/CM/CommonModeRange_%d.txt", ped_cm_file_num));
  return;
}

void replay_GEM(int RunNumber = 1559, int MaxEvent = 100) {

  // datafile name pattern
  const char *RunFileNamePattern = "ssp_gem_%d.evio.0";
  vector<TString> pathList;
  pathList.push_back("./");
  pathList.push_back("./ROOTfiles/GEM_EVIO_COSMIC/");

  const char *ROOTFileNamePattern = "ROOTfiles/lad_gem_%d_%d.root";
  TString ROOTFileName = Form(ROOTFileNamePattern, RunNumber, MaxEvent);

  // Params
  gHcParms->Define("gen_run_number", "Run Number", RunNumber);
  gHcParms->AddString("g_ctp_database_filename", "DBASE/LAD/standard.database");
  gHcParms->Load(gHcParms->GetString("g_ctp_database_filename"), RunNumber);
  gHcParms->Load(gHcParms->GetString("g_ctp_parm_filename"));
  gHcParms->Load(gHcParms->GetString("g_ctp_kinematics_filename"), RunNumber);
  // gHcParms->Load("PARAM/LAD/GEM/lgem_geom.param");
  // gHcParms->Load("DB_LAD/lgem_chan.map");
  // gHcParms->Load("DB_LAD/lgem_cuts.param");

  // Load correct GEM common mode and pedestal files
  load_GEM_CM_PED(RunNumber);
  // Detector map
  //  gHcDetectorMap = new THcDetectorMap();
  //  gHcDetectorMap->Load("detector.map");

  THcLADSpectrometer *lad = new THcLADSpectrometer("L", "LAD");
  gHaApps->Add(lad);

  THcLADGEM *gem = new THcLADGEM("gem", "gem");
  lad->AddDetector(gem);

  THcAnalyzer *analyzer = new THcAnalyzer;
  THaEvent *event       = new THaEvent;

  //  THcRun* run = new THcRun(pathList, Form(RunFileNamePattern, RunNumber) );
  THaRunBase *run = new THaRun(pathList, Form(RunFileNamePattern, RunNumber));

  run->SetEventRange(1, MaxEvent);
  // run->SetDataRequired(THaRunBase::kDate|THaRunBase::kRunNumber);
  run->SetDataRequired(0);
  run->Print();

  analyzer->SetEvent(event);
  analyzer->SetCountMode(2); // 2 = counter is event number

  analyzer->SetCrateMapFileName("MAPS/db_cratemap_gem_testlab.dat");
  analyzer->SetOutFile(ROOTFileName.Data());
  analyzer->SetOdefFile("DEF-files/LAD/PRODUCTION/lstackana_production_gem.def");
   analyzer->SetCutFile( "DEF-files/LAD/PRODUCTION/CUTS/lstackana_production_cuts_gem.def" );
  //  analyzer->SetCrateMapFileName("cratemap.dat");

  analyzer->Process(run);
}
