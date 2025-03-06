// #include "../../LAD/LAD_link_defs.h"
void replay_SSP_SIM(int RunNumber = 111, int MaxEvent = 100) {


  // Params
  gHcParms->Define("gen_run_number", "Run Number", RunNumber);
  gHcParms->AddString("g_ctp_database_filename", "DBASE/LAD/standard.database");
  gHcParms->Load(gHcParms->GetString("g_ctp_database_filename"), RunNumber);
  gHcParms->Load(gHcParms->GetString("g_ctp_parm_filename"));
  gHcParms->Load(gHcParms->GetString("g_ctp_kinematics_filename"), RunNumber);
  // gHcParms->Load("PARAM/LAD/GEM/lgem_geom.param");
  // gHcParms->Load("MAPS/LAD/DETEC/GEM/lgem_chan_mc.map");
  // gHcParms->Load("DB_LAD/lgem_cuts.param");

  // Detector map
  // gHcDetectorMap = new THcDetectorMap();
  // gHcDetectorMap->Load("MAPS/LAD/DETEC/GEM/lgem_chan_mc.map");//Think we can use the same one as for real data, and just route sim to same maping

  THcLADSpectrometer *lad = new THcLADSpectrometer("L", "LAD");
  gHaApps->Add(lad);

  THcLADGEM *gem = new THcLADGEM("gem", "gem");
  lad->AddDetector(gem);

  // THcConfigEvtHandler *ev125 = new THcConfigEvtHandler("HC", "Config Event type 125");
  // gHaEvtHandlers->Add(ev125);

  THcAnalyzer *analyzer = new THcAnalyzer;

  // Set the decoder to use simulation input
  THcInterface::SetDecoder(LADSimDecoder::Class());

  THaEvent *event = new THaEvent;


  TString run_file = "../libLADdig/test_scripts/lad_hodo_gem_sim.root";
  // TString run_file = "../libLADdig/test_scripts/lad_hodo_gem_sim_strip_test.root";
  THaRunBase *run  = new LADSimFile(run_file.Data(), "lad", "");
  // run->SetRunParamClass("THcRunParameters");
  
  run->SetEventRange(1, MaxEvent);
  // run->SetDataRequired(THaRunBase::kDate | THaRunBase::kRunNumber);
  run->SetDataRequired(0);
  run->Print();

  analyzer->SetEvent(event);
  analyzer->SetCountMode(2); // 2 = counter is event number

  analyzer->SetCrateMapFileName("MAPS/db_cratemap_mc.dat");
  TString outfile = Form("ROOTfiles/gem_mc_%d.root", RunNumber);
  analyzer->SetOutFile(outfile);
  analyzer->SetOdefFile("DEF-files/LAD/PRODUCTION/GEM/test.def");
  //  analyzer->SetCutFile( "test_cuts.def" );
  //  analyzer->SetCrateMapFileName("cratemap.dat");

  analyzer->Process(run);
}
