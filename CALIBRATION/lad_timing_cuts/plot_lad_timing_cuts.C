#include "../../SCRIPTS/LAD/LAD_link_defs.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include <iostream>
#include <string>

const static int MAXHITS                   = 1000;
const static int MAXPADDLES                = 11;
const static int MAXPLANES                 = 6;
const static string plane_names[MAXPLANES] = {"000", "001", "100", "101", "200", "REFBAR"};
const static double TDC_TO_TIME            = 0.09766; // ns
const static int NBINS_TDC                 = 500;
const static double TDC_MIN                = 0.0;
const static double TDC_MAX                = 45000.0 * TDC_TO_TIME;
const static double NBINS_ADC_TDC          = 50;
const static double ADC_TDC_MIN            = -50.0;
const static double ADC_TDC_MAX            = 50.0;

void plot_lad_timing_cuts(TString infile, int RunNumber, TString outfile = "move_me.root", const char* spec_prefix = "H") {
  gROOT->SetBatch(kTRUE); // do not display plots
  gHcParms->Define("gen_run_number", "Run Number", RunNumber);
  gHcParms->AddString("g_ctp_database_filename", "DBASE/LAD_COIN/standard.database");
  gHcParms->Load(gHcParms->GetString("g_ctp_database_filename"), RunNumber);
  gHcParms->Load(gHcParms->GetString("g_ctp_parm_filename"));
  gHcParms->Load(gHcParms->GetString("g_ctp_kinematics_filename"), RunNumber);
  // Load parameters for SHMS trigger configuration
  gHcParms->Load("PARAM/TRIG/tshms.param");
  gHcParms->Load("PARAM/TRIG/thms.param");


  const int nrefcut     = 1;
  Int_t *flad_tdcrefcut = new Int_t[nrefcut];
  Int_t *flad_adcrefcut = new Int_t[nrefcut];

  Int_t fladhodo_tdcmin              = 0;
  Int_t fladhodo_tdcmax              = 0;
  UInt_t fHodoScin                   = MAXPLANES * MAXPADDLES; // This could be determined from hcana DB
  Double_t *fHodoTopAdcTimeWindowMin = new Double_t[fHodoScin];
  Double_t *fHodoTopAdcTimeWindowMax = new Double_t[fHodoScin];
  Double_t *fHodoBtmAdcTimeWindowMin = new Double_t[fHodoScin];
  Double_t *fHodoBtmAdcTimeWindowMax = new Double_t[fHodoScin];
  Double_t *fADC_TDC_Offset          = new Double_t[MAXPLANES];
  Double_t *fTDC_Offset              = new Double_t[MAXPLANES];

  const char *prefix = "l";
  DBRequest windowList[] = {
      {"hodo_tdcrefcut", flad_tdcrefcut, kInt, 1, 1},
      {"hodo_adcrefcut", flad_adcrefcut, kInt, 1, 1},
      {"ladhodo_tdc_min", &fladhodo_tdcmin, kInt, 1, 1},
      {"ladhodo_tdc_max", &fladhodo_tdcmax, kInt, 1, 1},
      {"ladhodo_TopAdcTimeWindowMin", fHodoTopAdcTimeWindowMin, kDouble, fHodoScin, 1},
      {"ladhodo_TopAdcTimeWindowMax", fHodoTopAdcTimeWindowMax, kDouble, fHodoScin, 1},
      {"ladhodo_BtmAdcTimeWindowMin", fHodoBtmAdcTimeWindowMin, kDouble, fHodoScin, 1},
      {"ladhodo_BtmAdcTimeWindowMax", fHodoBtmAdcTimeWindowMax, kDouble, fHodoScin, 1},
      {0},
  };

  // gHcParms->LoadParmValues((DBRequest*)&trigRefList, trig_prefix);
  gHcParms->LoadParmValues((DBRequest *)&windowList, prefix);

  std::string spec_prefix_lower = spec_prefix;
  std::transform(spec_prefix_lower.begin(), spec_prefix_lower.end(), spec_prefix_lower.begin(), ::tolower);
  DBRequest list[] = {
      {"ladhodo_tdc_offset", fTDC_Offset, kDouble, (UInt_t)MAXPLANES, 1},
      {"ladhodo_adc_tdc_offset", fADC_TDC_Offset, kDouble, (UInt_t)MAXPLANES, 1},
      {0},
  };
  gHcParms->LoadParmValues((DBRequest *)&list, spec_prefix_lower.c_str());


  // cout << fdc_tdcrefcut[0] << endl;
  // cout << fhodo_tdcrefcut[0] << endl;
  // cout << fhodo_adcrefcut[0] << endl;
  // cout << fngcer_adcrefcut[0] << endl;
  // cout << fhgcer_adcrefcut[0] << endl;
  // cout << faero_adcrefcut[0] << endl;
  // cout << fcal_adcrefcut[0] << endl;

  gStyle->SetOptStat(0);
  gROOT->SetBatch(kTRUE); // do not display plots
  const int n_pdc_refs = 10;
  const int n_pT_refs  = 2;
  TFile *f1            = new TFile(infile, "UPDATE");
  if (f1->IsZombie()) {
    cout << "Cannot find : " << infile << endl;
    return;
  }
  TTree *tree = nullptr;
  f1->GetObject("T", tree);
  if (!tree) {
    std::cerr << "Failed to retrieve tree 'T' from file: " << infile << std::endl;
    return;
  }
  TFile *f2 = new TFile(outfile, "RECREATE");

  TCanvas *canvases[MAXPLANES];
  for (int i = 0; i < MAXPLANES; i++) {
    canvases[i] = new TCanvas(Form("c_plane_%d", i + 1), Form("Canvas for Plane %d", i + 1));
    canvases[i]->Divide(1, 2);
  }

  vector<TH1F *> pladhod_time_raw_top;
  vector<TH1F *> pladhod_time_raw_btm;

  for (int i = 0; i < MAXPLANES; i++) {
    TString branchNameTop = Form("%s.ladhod.%s.TopTdcTime", spec_prefix, plane_names[i].c_str());
    TString branchNameBtm = Form("%s.ladhod.%s.BtmTdcTime", spec_prefix, plane_names[i].c_str());
    TString histNameTop   = Form("h_plane_%s_top", plane_names[i].c_str());
    TString histNameBtm   = Form("h_plane_%s_btm", plane_names[i].c_str());
    TString histTitleTop  = Form("Plane %s Top TDC Time Raw", plane_names[i].c_str());
    TString histTitleBtm  = Form("Plane %s Bottom TDC Time Raw", plane_names[i].c_str());

    TH1F *histTop = new TH1F(histNameTop, histTitleTop, NBINS_TDC, TDC_MIN, TDC_MAX);
    TH1F *histBtm = new TH1F(histNameBtm, histTitleBtm, NBINS_TDC, TDC_MIN, TDC_MAX);

    tree->Draw(Form("((%s+ %f)*%f )>>%s", branchNameTop.Data(), fTDC_Offset[i], TDC_TO_TIME, histNameTop.Data()), "",
               "goff");
    tree->Draw(Form("((%s+ %f)*%f )>>%s", branchNameBtm.Data(), fTDC_Offset[i], TDC_TO_TIME, histNameBtm.Data()), "",
               "goff");

    pladhod_time_raw_top.push_back(histTop);
    pladhod_time_raw_btm.push_back(histBtm);

    if (histTop->GetEntries() == 0) {
      std::cerr << "No entries for histogram: " << histNameTop << " (branch: " << branchNameTop << ")" << std::endl;
    }
    if (histBtm->GetEntries() == 0) {
      std::cerr << "No entries for histogram: " << histNameBtm << " (branch: " << branchNameBtm << ")" << std::endl;
    }
  }

  for (int i = 0; i < MAXPLANES; i++) {
    canvases[i]->cd(1);
    pladhod_time_raw_top[i]->Draw();

    TLine *lineTopMin = new TLine(fladhodo_tdcmin * TDC_TO_TIME, 0, fladhodo_tdcmin * TDC_TO_TIME,
                                  pladhod_time_raw_top[i]->GetMaximum());
    TLine *lineTopMax = new TLine(fladhodo_tdcmax * TDC_TO_TIME, 0, fladhodo_tdcmax * TDC_TO_TIME,
                                  pladhod_time_raw_top[i]->GetMaximum());
    lineTopMin->SetLineColor(kRed);
    lineTopMax->SetLineColor(kRed);
    lineTopMin->SetLineStyle(2);
    lineTopMax->SetLineStyle(2);
    lineTopMin->Draw("same");
    lineTopMax->Draw("same");

    canvases[i]->cd(2);
    pladhod_time_raw_btm[i]->Draw();
    TLine *lineBtmMin = new TLine(fladhodo_tdcmin * TDC_TO_TIME, 0, fladhodo_tdcmin * TDC_TO_TIME,
                                  pladhod_time_raw_btm[i]->GetMaximum());
    TLine *lineBtmMax = new TLine(fladhodo_tdcmax * TDC_TO_TIME, 0, fladhodo_tdcmax * TDC_TO_TIME,
                                  pladhod_time_raw_btm[i]->GetMaximum());
    lineBtmMin->SetLineColor(kRed);
    lineBtmMax->SetLineColor(kRed);
    lineBtmMin->SetLineStyle(2);
    lineBtmMax->SetLineStyle(2);
    lineBtmMin->Draw("same");
    lineBtmMax->Draw("same");

    canvases[i]->cd(1);
    TLatex *latexTop = new TLatex();
    latexTop->SetNDC();
    latexTop->SetTextSize(0.03);
    latexTop->SetTextAlign(13);
    latexTop->DrawLatex(0.7, 0.85, Form("Low Cut: %.2f ns", fladhodo_tdcmin * TDC_TO_TIME));
    latexTop->DrawLatex(0.7, 0.8, Form("High Cut: %.2f ns", fladhodo_tdcmax * TDC_TO_TIME));

    canvases[i]->cd(2);
    TLatex *latexBtm = new TLatex();
    latexBtm->SetNDC();
    latexBtm->SetTextSize(0.03);
    latexBtm->SetTextAlign(13);
    latexBtm->DrawLatex(0.7, 0.85, Form("Low Cut: %.2f ns", fladhodo_tdcmin * TDC_TO_TIME));
    latexBtm->DrawLatex(0.7, 0.8, Form("High Cut: %.2f ns", fladhodo_tdcmax * TDC_TO_TIME));
  }

  for (int i = 0; i < MAXPLANES; i++) {
    canvases[i]->Write();
  }

  // Timing Window Cuts

  TCanvas *c_timing_window[MAXPLANES];
  for (int i = 0; i < MAXPLANES; i++) {
    c_timing_window[i] = new TCanvas(Form("c_timing_window_plane_%d", i + 1),
                                     Form("Timing Window Cuts for Plane %s", plane_names[i].c_str()), 800, 600);
    c_timing_window[i]->Divide(3, 4);
  }
  for (int plane = 0; plane < MAXPLANES; plane++) {
    for (int paddle = 0; paddle < MAXPADDLES; paddle++) {
      c_timing_window[plane]->cd(paddle + 1);

      TString branchNameTopAdcTdcDiff =
          Form("%s.ladhod.%s.GoodTopAdcTdcDiffTime[%d]", spec_prefix, plane_names[plane].c_str(), paddle);
      TString branchNameBtmAdcTdcDiff =
          Form("%s.ladhod.%s.GoodBtmAdcTdcDiffTime[%d]", spec_prefix, plane_names[plane].c_str(), paddle);

      TString histNameTop  = Form("h_plane_%d_paddle_%d_top_diff", plane, paddle);
      TString histTitleTop = Form("Plane %d Paddle %d Top TDC - ADC Time", plane, paddle);
      TString histNameBtm  = Form("h_plane_%d_paddle_%d_btm_diff", plane, paddle);
      TString histTitleBtm = Form("Plane %d Paddle %d Bottom TDC - ADC Time", plane, paddle);

      TH1F *histTop = new TH1F(histNameTop, histTitleTop, NBINS_ADC_TDC, ADC_TDC_MIN, ADC_TDC_MAX);
      TH1F *histBtm = new TH1F(histNameBtm, histTitleBtm, NBINS_ADC_TDC, ADC_TDC_MIN, ADC_TDC_MAX);

      tree->Draw(Form("%s>>%s", branchNameTopAdcTdcDiff.Data(), histNameTop.Data()), "", "goff");
      tree->Draw(Form("%s>>%s", branchNameBtmAdcTdcDiff.Data(), histNameBtm.Data()), "", "goff");

      histTop->SetLineColor(kBlue);
      histBtm->SetLineColor(kGreen);

      histTop->Draw();
      histBtm->Draw("same");

      int idx = MAXPLANES * paddle + plane;
      TLine *lineTopMin = new TLine(fHodoTopAdcTimeWindowMin[idx], 0,
                                    fHodoTopAdcTimeWindowMin[idx], histTop->GetMaximum());
      TLine *lineTopMax = new TLine(fHodoTopAdcTimeWindowMax[idx], 0,
                                    fHodoTopAdcTimeWindowMax[idx], histTop->GetMaximum());
      TLine *lineBtmMin = new TLine(fHodoBtmAdcTimeWindowMin[idx], 0,
                                    fHodoBtmAdcTimeWindowMin[idx], histBtm->GetMaximum());
      TLine *lineBtmMax = new TLine(fHodoBtmAdcTimeWindowMax[idx], 0,
                                    fHodoBtmAdcTimeWindowMax[idx], histBtm->GetMaximum());

      lineTopMin->SetLineColor(kRed);
      lineTopMax->SetLineColor(kRed);
      lineBtmMin->SetLineColor(kMagenta);
      lineBtmMax->SetLineColor(kMagenta);

      lineTopMin->SetLineStyle(2);
      lineTopMax->SetLineStyle(2);
      lineBtmMin->SetLineStyle(2);
      lineBtmMax->SetLineStyle(2);

      lineTopMin->Draw("same");
      lineTopMax->Draw("same");
      lineBtmMin->Draw("same");
      lineBtmMax->Draw("same");

      TLatex *latexTop = new TLatex();
      latexTop->SetNDC();
      latexTop->SetTextSize(0.03);
      latexTop->SetTextAlign(13);
      latexTop->DrawLatex(0.7, 0.85, Form("Low Cut: %.2f ns", fHodoTopAdcTimeWindowMin[plane * MAXPADDLES + paddle]));
      latexTop->DrawLatex(0.7, 0.8, Form("High Cut: %.2f ns", fHodoTopAdcTimeWindowMax[plane * MAXPADDLES + paddle]));
      TLatex *latexBtm = new TLatex();
      latexBtm->SetNDC();
      latexBtm->SetTextSize(0.03);
      latexBtm->SetTextAlign(13);
      latexBtm->DrawLatex(0.7, 0.85, Form("Low Cut: %.2f ns", fHodoBtmAdcTimeWindowMin[plane * MAXPADDLES + paddle]));
      latexBtm->DrawLatex(0.7, 0.8, Form("High Cut: %.2f ns", fHodoBtmAdcTimeWindowMax[plane * MAXPADDLES + paddle]));

      histTop->GetXaxis()->SetTitle("TDC - ADC Time (ns)");
      histTop->GetYaxis()->SetTitle("Counts");
      histBtm->GetXaxis()->SetTitle("TDC - ADC Time (ns)");
      histBtm->GetYaxis()->SetTitle("Counts");
    }
    c_timing_window[plane]->Write();
  }

  // Clean up
  for (int i = 0; i < MAXPLANES; i++) {
    delete canvases[i];
    delete pladhod_time_raw_top[i];
    delete pladhod_time_raw_btm[i];
  }
  delete[] flad_tdcrefcut;
  delete[] flad_adcrefcut;
  delete[] fHodoTopAdcTimeWindowMin;
  delete[] fHodoTopAdcTimeWindowMax;
  delete[] fHodoBtmAdcTimeWindowMin;
  delete[] fHodoBtmAdcTimeWindowMax;
  f1->Close();
  delete f1;
  // f2->Write();
  f2->Close();
  delete f2;

  return;
}

// void plot_lad_timing_cuts(TString infile, int RunNumber, TString outfile = "move_me.root") {

//   // run_shms_reference_time_setup(infile, RunNumber, outfile, "coin");
//   plot_lad_timing_cuts(infile, RunNumber, outfile, "shms");
//   return;
// }
