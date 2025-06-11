/*HMS Hodo Calibration script to determine scin. prop velocity, cable time offsets per paddle,
  and time offsets between paddles in different planes
  Author: Carlos Yero
  Dated: June 5, 2018
*/

#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"
#include "TTree.h"
#include <TDecompSVD.h>
#include <TF1.h>
#include <TGraph.h>
#include <TH2.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TNtuple.h>
#include <TObjArray.h>
#include <TPaveLabel.h>
#include <TPolyLine.h>
#include <TProfile.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TVectorD.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace std;

const static double TdcTimeTWCorr_MAX = 10000.0; // Was originally 100.
const static int minADC_PulseAmp = 50; // Minimum ADC Pulse Amplitude to consider for Time Walk correction

// void fitHodoCalib(TString filename, Int_t runNUM, Bool_t cosmic_flag = kFALSE) {
void fitHodoCalib(Int_t runNUM) {

  gROOT->SetBatch(kTRUE);
  // TString filename = Form("../../ROOTfiles/COSMICS/LAD_wREF_cosmic_hall_%d_-1.root", runNUM);
  TString filename = Form("../../ROOTfiles/COSMICS/LAD_wGEM_cosmic_hall_%d_-1.root", runNUM);

  gStyle->SetOptFit();
  // gROOT->SetBatch(kTRUE); // do not display plots

  Int_t evtNUM = 30000;

  TFile *data_file = new TFile(filename, "READ");
  TTree *T         = (TTree *)data_file->Get("T");

  /******Define Fixed Quantities********/
  static const Int_t NPLANES                       = 6;
  static const Int_t MAX_PADDLES                   = 11;
  static const Int_t SIDES                         = 2;
  TString spec                                     = "L";
  TString det                                      = "ladhod";
  TString pl_names[NPLANES]                        = {"000", "001", "100", "101", "200", "REFBAR"};
  TString side_names[2]                            = {"GoodTop", "GoodBtm"};
  TString nsign[2]                                 = {"+", "-"};
  Int_t maxPMT[NPLANES]                            = {11, 11, 11, 11, 11, 1};
  Int_t refBar                                     = 0;
  Int_t refPlane                                   = 0;
  Double_t lladhodo_velArr[NPLANES][MAX_PADDLES]   = {0.0}; // store hhodo velocity parameters (1/slope of the line fit)
  Double_t lladhodo_cableArr[NPLANES][MAX_PADDLES] = {0.0}; // store hhodo cableLength differences (y-int of line fit)
  Double_t lladhodo_LCoeff[NPLANES][MAX_PADDLES]   = {0.0}; // Variables to write out LCoeff. parameter file
  Double_t lladhodo_sigArr[NPLANES][MAX_PADDLES]   = {0.0}; // store hhodo sigma parameters
  Double_t vp                                      = 30.0;  // speed of light [cm/ns]
  Double_t paddle_vel_fit_max_threshold            = 0.1;   // 0-1 range, 0.1 = 10% threshold
  Double_t barlength[NPLANES][MAX_PADDLES]         = {
      {387.5, 393.3, 399.0, 404.8, 410.5, 416.3, 422.0, 427.8, 433.6, 439.3, 445.1},
      {387.5, 393.3, 399.0, 404.8, 410.5, 416.3, 422.0, 427.8, 433.6, 439.3, 445.1},
      {387.5, 393.3, 399.0, 404.8, 410.5, 416.3, 422.0, 427.8, 433.6, 439.3, 445.1},
      {387.5, 393.3, 399.0, 404.8, 410.5, 416.3, 422.0, 427.8, 433.6, 439.3, 445.1},
      {387.5, 393.3, 399.0, 404.8, 410.5, 416.3, 422.0, 427.8, 433.6, 439.3, 445.1},
      {387.5, 393.3, 399.0, 404.8, 410.5, 416.3, 422.0, 427.8, 433.6, 439.3, 445.1}};

  static const Int_t TDC_T_NBINS         = 1000;
  static const Double_t TDC_T_MIN        = -60.0;
  static const Double_t TDC_T_MAX        = 60.0;
  static const Int_t ADC_TDC_diff_NBINS  = 100;
  static const Double_t ADC_TDC_diff_MIN = -360.0;
  static const Double_t ADC_TDC_diff_MAX = -330.0;
  static const Int_t DiffTime_NBINS      = 1000;
  static const Double_t DiffTime_MIN     = -40.0;
  static const Double_t DiffTime_MAX     = 40.0;
  /******Define Leafs to be read from TTree******/

  //---Names---
  TString base;
  TString nTdcTimeUnCorr;
  TString nTdcTimeTWCorr;
  TString nAdcPulseTime;
  TString nAdcPulseAmp;

  TString nTrackXPos;
  TString nTrackYPos;
  TString nDiffTWCorr;

  //---Variables---
  Double_t TdcTimeUnCorr[NPLANES][SIDES][MAX_PADDLES];
  Double_t TdcTimeTWCorr[NPLANES][SIDES][MAX_PADDLES];
  Double_t AdcPulseTime[NPLANES][SIDES][MAX_PADDLES];
  Double_t AdcPulseAmp[NPLANES][SIDES][MAX_PADDLES];

  Double_t hcal_etrkNorm;
  Double_t hcer_npeSum;
  Double_t hdc_ntrack;
  Double_t hod_nhits[NPLANES];
  Double_t beta;

  /******Define Matrices/Vectors and related *******/
  // Int_t good_pad[PLANES]; // keep track of good paddle hits
  // Int_t ngood = 0;

  // static const Int_t npar = 55 - 1; // reference paddle ??? is fixed (so we actually have 54)
  // static const Int_t ROW  = 6 * evtNUM;
  // static const Int_t COL  = npar;

  // Int_t row1, row2, row3, row4, row5, row6; // keep track of ith row element
  // Int_t cnt;                                // keep track of good plane hits

  // /*******Define Canvas and Histograms*******/

  // //----Canvas----
  TCanvas *TWAvg_canv[NPLANES];
  TCanvas *TWAvg_canv_2D[NPLANES];
  TCanvas *TWDiff_canv[NPLANES];

  TCanvas *TWUnCorr_canv[NPLANES][SIDES];
  TCanvas *TWCorr_canv[NPLANES][SIDES];

  TCanvas *TWCorr_BarLCoef[NPLANES];

  // //----Histograms----
  TH1F *h1Hist_TWAvg[NPLANES][MAX_PADDLES];     // (TWCorr_Top + TWCorr_Btm) / 2       <-------
  TH1F *h1Hist_TWAvg_CUT[NPLANES][MAX_PADDLES]; //<------
  TH1F *h1Hist_TWDiff[NPLANES][MAX_PADDLES];    // (TWCorr_Top - TWCorr_Btm) / 2       <-------

  TH2F *h2Hist_TW_UnCorr[NPLANES][SIDES][MAX_PADDLES]; // Time-Walk Uncorrected vs. ADC Pulse Amp Hist
  TH2F *h2Hist_TW_Corr[NPLANES][SIDES][MAX_PADDLES];   // Time-Walk Corrected vs. ADC Pulse Amp Hist

  TH1F *h1Hist_TWCorr_BarLCoef[NPLANES][MAX_PADDLES]; // Time-Walk Corrected vs. ADC Pulse Amp Hist
  // /*******Define Fit Functions and Related Variables*******/

  Double_t Mean;   // variable to get Mean to make a 3sig cut
  Double_t StdDev; // variable to get satndar deviation to make a 3sig cut
  Double_t nSig;   // multiple of Sigma used for sigmaCut

  /********Initialize HISTOS and GET TTREE VARIABLES*********/

  // Loop over hodo planes
  for (Int_t npl = 0; npl < NPLANES; npl++) {

    // Loop over hodo side
    for (Int_t side = 0; side < SIDES; side++) {

      // Loop over hodo PMTs
      for (Int_t ipmt = 0; ipmt < maxPMT[npl]; ipmt++) {

        // Initialize Histograms
        h2Hist_TW_UnCorr[npl][side][ipmt] =
            new TH2F(Form("TW_UnCorr PMT %s%d%s", pl_names[npl].Data(), ipmt + 1, nsign[side].Data()),
                     Form("PMT %s%d%s: UnCorr. (TDC - ADC) Pulse Time vs. ADC Pulse Amplitude ", pl_names[npl].Data(),
                          ipmt + 1, nsign[side].Data()),
                     600, 0, 420, ADC_TDC_diff_NBINS, ADC_TDC_diff_MIN, ADC_TDC_diff_MAX);
        h2Hist_TW_Corr[npl][side][ipmt] =
            new TH2F(Form("TW_Corr PMT %s%d%s", pl_names[npl].Data(), ipmt + 1, nsign[side].Data()),
                     Form("PMT %s%d%s: Corr. (TDC - ADC) Pulse Time vs. ADC Pulse Amplitude ", pl_names[npl].Data(),
                          ipmt + 1, nsign[side].Data()),
                     600, 0, 420, ADC_TDC_diff_NBINS, ADC_TDC_diff_MIN, ADC_TDC_diff_MAX);

        h2Hist_TW_UnCorr[npl][side][ipmt]->GetYaxis()->SetTitle("Time Walk UnCorr.(TDC - ADC) Pulse Time (ns)");
        h2Hist_TW_UnCorr[npl][side][ipmt]->GetXaxis()->SetTitle("ADC Pulse Amplitude (mV)");
        h2Hist_TW_UnCorr[npl][side][ipmt]->GetXaxis()->CenterTitle();

        h2Hist_TW_Corr[npl][side][ipmt]->GetYaxis()->SetTitle("Time Walk Corr.(TDC - ADC) Pulse Time (ns)");

        h2Hist_TW_Corr[npl][side][ipmt]->GetXaxis()->SetTitle("ADC Pulse Amplitude (mV)");
        h2Hist_TW_Corr[npl][side][ipmt]->GetXaxis()->CenterTitle();

        if (side == 0) // require ONLY one side, since a time diff between two pmts at each end is taken
        {
          h1Hist_TWAvg[npl][ipmt] =
              new TH1F(Form("Avg. Time: Paddle %s%d", pl_names[npl].Data(), ipmt + 1),
                       Form("Paddle %s%d: Time-Walk Corrected Average Time", pl_names[npl].Data(), ipmt + 1),
                       TDC_T_NBINS, TDC_T_MIN, TDC_T_MAX);

          h1Hist_TWAvg_CUT[npl][ipmt] =
              new TH1F(Form("Avg. Time CUT: Paddle %s%d", pl_names[npl].Data(), ipmt + 1),
                       Form("Paddle %s%d: Time-Walk Corrected Average (CUT)", pl_names[npl].Data(), ipmt + 1),
                       TDC_T_NBINS, TDC_T_MIN, TDC_T_MAX);

          h1Hist_TWDiff[npl][ipmt] = new TH1F(
              Form("Diff. Time: Paddle %s%d", pl_names[npl].Data(), ipmt + 1),
              Form("Paddle %s%d: Time-Walk Corrected Difference Time (Top - Btm)", pl_names[npl].Data(), ipmt + 1),
              DiffTime_NBINS, DiffTime_MIN, DiffTime_MAX);

          // Set Axis Titles
          h1Hist_TWAvg[npl][ipmt]->GetXaxis()->SetTitle("Time-Walk Corr. TDC Average Paddle Time (ns)");
          h1Hist_TWAvg_CUT[npl][ipmt]->GetXaxis()->SetTitle("Time-Walk Corr. TDC Average Paddle Time (ns)");

          h1Hist_TWAvg[npl][ipmt]->GetXaxis()->CenterTitle();
          h1Hist_TWAvg_CUT[npl][ipmt]->GetXaxis()->CenterTitle();

          h2Hist_TW_UnCorr[npl][side][ipmt]->GetYaxis()->SetTitle("Time Walk UnCorr.(TDC - ADC) Pulse Time (ns)");
          h2Hist_TW_UnCorr[npl][side][ipmt]->GetXaxis()->SetTitle("ADC Pulse Amplitude (mV)");
        } // end require SINGLE side requirement

        //----Define TTree Leaf Names-----
        base = spec + "." + det + "." + pl_names[npl];

        nTdcTimeUnCorr = base + "." + side_names[side] + "TdcTimeUnCorr";
        nTdcTimeTWCorr = base + "." + side_names[side] + "TdcTimeWalkCorr";
        nAdcPulseTime  = base + "." + side_names[side] + "AdcPulseTime";
        nAdcPulseAmp   = base + "." + side_names[side] + "AdcPulseAmp";

        //------Set Branch Address-------
        T->SetBranchAddress(nTdcTimeUnCorr, &TdcTimeUnCorr[npl][side]);
        T->SetBranchAddress(nTdcTimeTWCorr, &TdcTimeTWCorr[npl][side]);
        T->SetBranchAddress(nAdcPulseTime, &AdcPulseTime[npl][side]);
        T->SetBranchAddress(nAdcPulseAmp, &AdcPulseAmp[npl][side]);

      } // end loop over hodo PMTs

    } // end loop over hodo side

  } // end loop over hodo planes

  //**************************************************************//
  // FIRST PASS OF EVENT LOOP (Get the StdDev) of (TDC+ + TDC-)/2 //
  //**************************************************************//
  cout << "Initializing 1st Pass of Event Loop: " << endl;

  Long64_t nentries = T->GetEntries();

  // Loop over all entries
  for (Long64_t i = 0; i < nentries; i++) {
    T->GetEntry(i);

    // Loop over hodo planes
    for (Int_t npl = 0; npl < NPLANES; npl++) {

      // Loop over pmt
      for (Int_t ipmt = 0; ipmt < maxPMT[npl]; ipmt++) {

        if (AdcPulseAmp[npl][0][ipmt] < minADC_PulseAmp || AdcPulseAmp[npl][1][ipmt] < minADC_PulseAmp) {
          continue; // skip if ADC Pulse Amplitude is below threshold
        }
        if (TdcTimeTWCorr[npl][0][ipmt] < TdcTimeTWCorr_MAX && TdcTimeTWCorr[npl][1][ipmt] < TdcTimeTWCorr_MAX) {
          // Fill Average TW Corr TDC Time
          h1Hist_TWAvg[npl][ipmt]->Fill((TdcTimeTWCorr[npl][0][ipmt] + TdcTimeTWCorr[npl][1][ipmt]) / 2.);

        } // end time cut

      } // end pmt loop

    } // end plane loop

    cout << std::setprecision(2) << double(i) / nentries * 100. << "  % " << std::flush << "\r";

  } // end loop over entries

  // Set cut on Sigma,
  nSig = 1;

  //************************************//
  //    SECOND PASS OF EVENT LOOP       //
  //************************************//

  cout << "Initializing 2nd Pass of Event Loop: " << endl;

  // Loop over all entries
  for (Long64_t i = 0; i < nentries; i++) {

    T->GetEntry(i);

    // Loop over hodo planes
    for (Int_t npl = 0; npl < NPLANES; npl++) {

      // Loop over plane side
      for (Int_t side = 0; side < SIDES; side++) {

        // Loop over pmt
        for (Int_t ipmt = 0; ipmt < maxPMT[npl]; ipmt++) {

          // Get Standard deviation from initial entry fill
          StdDev = h1Hist_TWAvg[npl][ipmt]->GetStdDev();
          Mean   = h1Hist_TWAvg[npl][ipmt]->GetMean();

          // FIll Uncorrected/Corrected Time Walk Histos
          h2Hist_TW_UnCorr[npl][side][ipmt]->Fill(AdcPulseAmp[npl][side][ipmt],
                                                  TdcTimeUnCorr[npl][side][ipmt] - AdcPulseTime[npl][side][ipmt]);
          h2Hist_TW_Corr[npl][side][ipmt]->Fill(AdcPulseAmp[npl][side][ipmt],
                                                TdcTimeTWCorr[npl][side][ipmt] - AdcPulseTime[npl][side][ipmt]);

          // Add Time Cuts to get rid of kBig - kBig values, which yielded high evt density at zero
          if (TdcTimeTWCorr[npl][0][ipmt] < TdcTimeTWCorr_MAX && TdcTimeTWCorr[npl][1][ipmt] < TdcTimeTWCorr_MAX) {

            if (side == 0) // require only one side, as a time diff between the two ends of a paddle is take
            {
              if ((((TdcTimeTWCorr[npl][0][ipmt] + TdcTimeTWCorr[npl][1][ipmt]) / 2.) > (Mean - nSig * StdDev)) &&
                  (((TdcTimeTWCorr[npl][0][ipmt] + TdcTimeTWCorr[npl][1][ipmt]) / 2.) < (Mean + nSig * StdDev))) {

                h1Hist_TWAvg_CUT[npl][ipmt]->Fill((TdcTimeTWCorr[npl][0][ipmt] + TdcTimeTWCorr[npl][1][ipmt]) / 2.);
                h1Hist_TWDiff[npl][ipmt]->Fill((TdcTimeTWCorr[npl][1][ipmt] - TdcTimeTWCorr[npl][0][ipmt]) / 2.);

              } // end 3SIGMA CUT of TW Corr Time

            } // end require ONLY single side

          } // end time cuts

        } // end pmt loop

      } // end side loop

    } // end plane loop

    cout << std::setprecision(2) << double(i) / nentries * 100. << "  % " << std::flush << "\r";

  } // end entry loop

  /***********DRAW HISTOGRAMS TO CANVAS***************/
  for (Int_t npl = 0; npl < NPLANES; npl++) {

    // Create Canvas to store TW-Corr Time/Dist vs. trk position
    TWAvg_canv[npl] = new TCanvas(Form("TWAvg_%d", npl), Form("TWAvg, plane %s", pl_names[npl].Data()), 1000, 700);
    TWAvg_canv_2D[npl] =
        new TCanvas(Form("TWAvg2D_%d", npl), Form("TWAvg2D, plane %s", pl_names[npl].Data()), 1000, 700);
    TWDiff_canv[npl] = new TCanvas(Form("TWDiff_%d", npl), Form("TWDiff, plane %s", pl_names[npl].Data()), 1000, 700);

    if (npl < 5) {
      TWAvg_canv[npl]->Divide(4, 3);
      TWDiff_canv[npl]->Divide(4, 3);
      TWAvg_canv_2D[npl]->Divide(4, 3);
    }

    // Loop over plane side
    for (Int_t side = 0; side < SIDES; side++) {
      // Create Canvas
      TWUnCorr_canv[npl][side] =
          new TCanvas(Form("TWUnCorrCanv%d%d", npl, side),
                      Form("plane %s_%s", pl_names[npl].Data(), side_names[side].Data()), 1000, 700);
      TWCorr_canv[npl][side] =
          new TCanvas(Form("TWCorrCanv%d%d", npl, side),
                      Form("plane %s_%s", pl_names[npl].Data(), side_names[side].Data()), 1000, 700);

      // Divide Canvas
      if (npl < 5) {
        TWUnCorr_canv[npl][side]->Divide(4, 3);
        TWCorr_canv[npl][side]->Divide(4, 3);
      }

      // Loop over pmt
      Int_t fit_status;

      for (Int_t ipmt = 0; ipmt < maxPMT[npl]; ipmt++) {

        TWUnCorr_canv[npl][side]->cd(ipmt + 1);
        h2Hist_TW_UnCorr[npl][side][ipmt]->Draw("colz");

        TWCorr_canv[npl][side]->cd(ipmt + 1);
        h2Hist_TW_Corr[npl][side][ipmt]->Draw("colz");
        fit_status = -1;

        // Require ONLY one side
        if (side == 0) {

          lladhodo_cableArr[npl][ipmt] = h1Hist_TWDiff[npl][ipmt]->GetMean();
          lladhodo_sigArr[npl][ipmt]   = h1Hist_TWDiff[npl][ipmt]->GetStdDev();

          // Calculate the value at which h1Hist_TWDiff drops below 0.1 times its peak for both ends
          double peak      = h1Hist_TWDiff[npl][ipmt]->GetMaximum();
          double peak_bin  = h1Hist_TWDiff[npl][ipmt]->GetMaximumBin();
          double threshold = paddle_vel_fit_max_threshold * peak;
          int bin_left     = peak_bin;
          int bin_right    = peak_bin;

          // Find the bin to the left where the value drops below the threshold
          while (bin_left > 1 && h1Hist_TWDiff[npl][ipmt]->GetBinContent(bin_left) > threshold) {
            bin_left--;
          }
          double value_below_threshold_btm = h1Hist_TWDiff[npl][ipmt]->GetBinCenter(bin_left);

          // Find the bin to the right where the value drops below the threshold
          while (bin_right < h1Hist_TWDiff[npl][ipmt]->GetNbinsX() &&
                 h1Hist_TWDiff[npl][ipmt]->GetBinContent(bin_right) > threshold) {
            bin_right++;
          }
          double value_below_threshold_top = h1Hist_TWDiff[npl][ipmt]->GetBinCenter(bin_right);

          lladhodo_velArr[npl][ipmt] = barlength[npl][ipmt] / (value_below_threshold_top - value_below_threshold_btm);

          // Draw 1D TWAvg Histos
          TWAvg_canv[npl]->cd(ipmt + 1);
          h1Hist_TWAvg[npl][ipmt]->SetLineColor(kBlack);
          h1Hist_TWAvg_CUT[npl][ipmt]->SetLineColor(kRed);
          h1Hist_TWAvg[npl][ipmt]->Draw();
          h1Hist_TWAvg_CUT[npl][ipmt]->Draw("same");

          // Draw 1D TWDiff Histos
          TWDiff_canv[npl]->cd(ipmt + 1);
          h1Hist_TWDiff[npl][ipmt]->Draw();

          // Draw a line at lladhodo_cableArr[npl][ipmt]
          TLine *line_cable = new TLine(lladhodo_cableArr[npl][ipmt], 0, lladhodo_cableArr[npl][ipmt],
                                        h1Hist_TWDiff[npl][ipmt]->GetMaximum());
          line_cable->SetLineColor(kBlack);
          line_cable->SetLineStyle(2);
          line_cable->Draw("same");

          // Draw lines at bin_left and bin_right
          TLine *line_left = new TLine(value_below_threshold_btm, 0, value_below_threshold_btm,
                                       h1Hist_TWDiff[npl][ipmt]->GetMaximum());
          line_left->SetLineColor(kRed);
          line_left->SetLineStyle(2);
          line_left->Draw("same");

          TLine *line_right = new TLine(value_below_threshold_top, 0, value_below_threshold_top,
                                        h1Hist_TWDiff[npl][ipmt]->GetMaximum());
          line_right->SetLineColor(kRed);
          line_right->SetLineStyle(2);
          line_right->Draw("same");
        } // end single SIDE requirement

      } // end pmt loop

    } // end side loop

  } // end plane loop

  /***********DRAW CANVAS FOR lladhodo_cableArr***************/
  TCanvas *cableArr_canv = new TCanvas("cableArr_canv", "Cable Time Offsets per Paddle", 1000, 700);
  cableArr_canv->Divide(3, 2);

  for (Int_t npl = 0; npl < NPLANES; npl++) {
    TGraph *g_CableArr = new TGraph(maxPMT[npl]);
    g_CableArr->SetTitle(Form("Cable Time Offsets per Paddle for Plane %s", pl_names[npl].Data()));
    g_CableArr->GetXaxis()->SetTitle("Paddle Number");
    g_CableArr->GetYaxis()->SetTitle("Cable Time Offset (ns)");
    g_CableArr->SetMarkerStyle(20);    // Set marker style to points
    g_CableArr->SetMarkerColor(kBlue); // Set marker color

    double minY = *min_element(lladhodo_cableArr[npl], lladhodo_cableArr[npl] + maxPMT[npl]);
    double maxY = *max_element(lladhodo_cableArr[npl], lladhodo_cableArr[npl] + maxPMT[npl]);
    g_CableArr->GetYaxis()->SetRangeUser(minY - 1, maxY + 1); // Adjust y range

    for (Int_t ipmt = 0; ipmt < maxPMT[npl]; ipmt++) {
      g_CableArr->SetPoint(ipmt, ipmt + 1, lladhodo_cableArr[npl][ipmt]);
    }
    cableArr_canv->cd(npl + 1);
    g_CableArr->Draw("AP"); // Draw only points

    // Draw a dashed line at y = 0
    TLine *line = new TLine(0, 0, maxPMT[npl] + 1, 0);
    line->SetLineColor(kBlack);
    line->SetLineStyle(2); // Dashed line
    line->Draw("same");
  }

  /***********DRAW CANVAS FOR lladhodo_velArr***************/
  TCanvas *velArr_canv = new TCanvas("velArr_canv", "Propagation Velocities per Paddle", 1000, 700);
  velArr_canv->Divide(3, 2);

  for (Int_t npl = 0; npl < NPLANES; npl++) {
    TGraph *g_VelArr = new TGraph(maxPMT[npl]);
    g_VelArr->SetTitle(Form("Propagation Velocities per Paddle for Plane %s", pl_names[npl].Data()));
    g_VelArr->GetXaxis()->SetTitle("Paddle Number");
    g_VelArr->GetYaxis()->SetTitle("Propagation Velocity (cm/ns)");
    g_VelArr->SetMarkerStyle(20);   // Set marker style to points
    g_VelArr->SetMarkerColor(kRed); // Set marker color

    double minY = *min_element(lladhodo_velArr[npl], lladhodo_velArr[npl] + maxPMT[npl]);
    double maxY = *max_element(lladhodo_velArr[npl], lladhodo_velArr[npl] + maxPMT[npl]);
    g_VelArr->GetYaxis()->SetRangeUser(minY - 1, maxY + 1); // Adjust y range

    for (Int_t ipmt = 0; ipmt < maxPMT[npl]; ipmt++) {
      g_VelArr->SetPoint(ipmt, ipmt + 1, lladhodo_velArr[npl][ipmt]);
    }
    velArr_canv->cd(npl + 1);
    g_VelArr->Draw("AP"); // Draw only points
  }

  /************WRITE FIT RESULTS TO PARAMETER FILE***************/

  ofstream outPARAM;
  outPARAM.open(Form("../../PARAM/LAD/HODO/ladhodo_Vpcalib_%d.param", runNUM));

  outPARAM << "; LAD Hodoscope Parameter File Containing propagation velocities per paddle " << endl;
  outPARAM << "; and signal cable time diff. offsets per paddle " << endl;
  outPARAM << Form("; Run %d ", runNUM) << endl;
  outPARAM << " " << endl;
  outPARAM << " " << endl;

  outPARAM << ";Propagation Velocity Per Paddle" << endl;
  outPARAM << "; ";
  for (Int_t npl = 0; npl < NPLANES; npl++) {
    outPARAM << setw(16) << pl_names[npl] << " ";
  }
  outPARAM << endl;
  outPARAM << "lladhodo_velFit = ";

  //--------Write Velocity Parameters to Param File-----
  for (Int_t ipmt = 0; ipmt < MAX_PADDLES; ipmt++) {
    for (Int_t npl = 0; npl < NPLANES; npl++) {
      if (npl == 0) {
        if (ipmt == 0) {
          outPARAM << lladhodo_velArr[npl][ipmt];
        } else {
          outPARAM << setw(26) << lladhodo_velArr[npl][ipmt];
        }
      } else {
        outPARAM << ", " << setw(15) << lladhodo_velArr[npl][ipmt];
      }
    }
    outPARAM << fixed << endl;
  }

  outPARAM << " " << endl;
  outPARAM << " " << endl;
  outPARAM << " " << endl;

  outPARAM << ";PMTs Signal Cable Time Diff. Per Paddle" << endl;
  outPARAM << "; ";
  for (Int_t npl = 0; npl < NPLANES; npl++) {
    outPARAM << setw(16) << pl_names[npl] << " ";
  }
  outPARAM << endl;
  outPARAM << "lladhodo_cableFit = ";

  // Write Cable Length Time Diff. Parameters to Param File
  for (Int_t ipmt = 0; ipmt < MAX_PADDLES; ipmt++) {
    for (Int_t npl = 0; npl < NPLANES; npl++) {
      if (npl == 0) {
        if (ipmt == 0) {
          outPARAM << lladhodo_cableArr[npl][ipmt];
        } else {
          outPARAM << setw(26) << lladhodo_cableArr[npl][ipmt];
        }
      } else {
        outPARAM << ", " << setw(15) << lladhodo_cableArr[npl][ipmt];
      }
    }
    outPARAM << fixed << endl;
  }

  outPARAM << " " << endl;
  outPARAM << " " << endl;
  outPARAM << " " << endl;

  outPARAM << ";PMTs Time Diff. Sigma Parameters" << endl;
  outPARAM << "; ";
  for (Int_t npl = 0; npl < NPLANES; npl++) {
    outPARAM << setw(16) << pl_names[npl] << " ";
  }
  outPARAM << endl;
  outPARAM << "lladhodo_TopSigma = ";

  // Write Sigma Parameters to file
  for (Int_t ipmt = 0; ipmt < MAX_PADDLES; ipmt++) {
    for (Int_t npl = 0; npl < NPLANES; npl++) {
      if (npl == 0) {
        if (ipmt == 0) {
          outPARAM << lladhodo_sigArr[npl][ipmt];
        } else {
          outPARAM << setw(26) << lladhodo_sigArr[npl][ipmt];
        }
      } else {
        outPARAM << ", " << setw(15) << lladhodo_sigArr[npl][ipmt];
      }
    }
    outPARAM << fixed << endl;
  }

  outPARAM << " " << endl;
  outPARAM << " " << endl;
  outPARAM << "lladhodo_BtmSigma = ";

  // Write Sigma Parameters to file

  cout << "FINISHED Getting Vp and Cable Fits . . . " << endl;
  cout << "Starting the code to fit Hodo Matrix . . . " << endl;

  /**********BEGIN CODE TO FIT HODO MATRIX**************/

  // Initialize h1Hist_TWCorr_BarLCoef histograms
  for (Int_t npl = 0; npl < NPLANES; npl++) {
    for (Int_t ipmt = 0; ipmt < maxPMT[npl]; ipmt++) {
      h1Hist_TWCorr_BarLCoef[npl][ipmt] =
          new TH1F(Form("TWCorr_BarLCoef_%s%d", pl_names[npl].Data(), ipmt + 1),
                   Form("TWCorr_BarLCoef, plane %s, paddle %d", pl_names[npl].Data(), ipmt + 1), DiffTime_NBINS,
                   DiffTime_MIN, DiffTime_MAX);
    }
  }
  Double_t good_TW_top, good_TW_btm, good_TW_top_ref, good_TW_btm_ref;
  for (Long64_t i = 0; i < nentries; i++) {
    T->GetEntry(i);
    // Loop over hodo planes
    for (Int_t npl = 0; npl < NPLANES; npl++) {

      // Loop over pmt
      for (Int_t ipmt = 0; ipmt < maxPMT[npl]; ipmt++) {
        if (TdcTimeTWCorr[npl][0][ipmt] < TdcTimeTWCorr_MAX && TdcTimeTWCorr[npl][1][ipmt] < TdcTimeTWCorr_MAX) {

          good_TW_top = TdcTimeTWCorr[npl][0][ipmt];
          good_TW_btm = TdcTimeTWCorr[npl][1][ipmt] - 2 * lladhodo_cableArr[npl][ipmt]; // IMPORTANT: Apply cable time
          // correction obtained from fits
          good_TW_top_ref = TdcTimeTWCorr[refPlane][0][refBar];
          good_TW_btm_ref = TdcTimeTWCorr[refPlane][1][refBar] - 2 * lladhodo_cableArr[refPlane][refBar];
          if (good_TW_top_ref == 0) {
            continue;
          }
          h1Hist_TWCorr_BarLCoef[npl][ipmt]->Fill(((good_TW_top_ref + good_TW_btm_ref) - (good_TW_top + good_TW_btm)) /
                                                  2);

        } // end time cut

      } // end pmt loop

    } // end plane loop
    cout << std::setprecision(2) << double(i) / nentries * 100. << "  % " << std::flush << "\r";
  } // end event loop

  // Loop through h1Hist_TWCorr_BarLCoef and fill the mean into lladhodo_LCoeff
  double refBarTime = h1Hist_TWCorr_BarLCoef[refPlane][refBar]->GetMean();
  for (Int_t npl = 0; npl < NPLANES; npl++) {
    for (Int_t ipmt = 0; ipmt < maxPMT[npl]; ipmt++) {
      // if (npl == refPlane && ipmt == refBar) {
      //   continue;
      // }
      lladhodo_LCoeff[npl][ipmt] = (h1Hist_TWCorr_BarLCoef[npl][ipmt]->GetMean());
    }
  }
  // Create Canvas to store TW-Corr Time/Dist vs. trk position
  for (Int_t npl = 0; npl < NPLANES; npl++) {

    // Create Canvas
    TWCorr_BarLCoef[npl] = new TCanvas(Form("TWCorr_BarLCoef_%d", npl),
                                       Form("TWCorr_BarLCoef, plane %s", pl_names[npl].Data()), 1000, 700);

    // Divide Canvas
    if (npl < 5) {
      TWCorr_BarLCoef[npl]->Divide(4, 3);
    }

    // Loop over pmt
    for (Int_t ipmt = 0; ipmt < maxPMT[npl]; ipmt++) {

      TWCorr_BarLCoef[npl]->cd(ipmt + 1);
      h1Hist_TWCorr_BarLCoef[npl][ipmt]->Draw();

      lladhodo_LCoeff[npl][ipmt] =
          h1Hist_TWCorr_BarLCoef[npl][ipmt]->GetMean() - h1Hist_TWCorr_BarLCoef[refPlane][refBar]->GetMean();
    } // end pmt loop

  } // end plane loop

  /***********DRAW CANVAS FOR lladhodo_LCoeff***************/
  TCanvas *LCoeff_canv = new TCanvas("LCoeff_canv", "Timing Corrections per Paddle", 1000, 700);
  LCoeff_canv->Divide(3, 2);

  for (Int_t npl = 0; npl < NPLANES; npl++) {
    TGraph *g_LCoeff = new TGraph(maxPMT[npl]);
    g_LCoeff->SetTitle(Form("Timing Corrections per Paddle for Plane %s", pl_names[npl].Data()));
    g_LCoeff->GetXaxis()->SetTitle("Paddle Number");
    g_LCoeff->GetYaxis()->SetTitle("Timing Correction (ns)");
    g_LCoeff->SetMarkerStyle(20);   // Set marker style to points
    g_LCoeff->SetMarkerColor(kRed); // Set marker color

    double minY = *min_element(lladhodo_LCoeff[npl], lladhodo_LCoeff[npl] + maxPMT[npl]);
    double maxY = *max_element(lladhodo_LCoeff[npl], lladhodo_LCoeff[npl] + maxPMT[npl]);
    g_LCoeff->GetYaxis()->SetRangeUser(minY - 1, maxY + 1); // Adjust y range

    for (Int_t ipmt = 0; ipmt < maxPMT[npl]; ipmt++) {
      g_LCoeff->SetPoint(ipmt, ipmt + 1, lladhodo_LCoeff[npl][ipmt]);
    }
    LCoeff_canv->cd(npl + 1);
    g_LCoeff->Draw("AP"); // Draw only points
  }

  // Write Fit Results to Parameter File

  outPARAM << "" << endl;
  outPARAM << "" << endl;
  outPARAM << "" << endl;
  outPARAM << ";Timing Corrections Per Paddle, where 000 Paddle 00 has been set as the reference paddle" << endl;
  outPARAM << "; ";
  for (Int_t npl = 0; npl < NPLANES; npl++) {
    outPARAM << setw(16) << pl_names[npl] << " ";
  }
  outPARAM << endl;
  outPARAM << "lladhodo_LCoeff = ";

  // Write Lambda Time Coeff. Parameters to Param File
  for (Int_t ipmt = 0; ipmt < MAX_PADDLES; ipmt++) {
    for (Int_t npl = 0; npl < NPLANES; npl++) {
      if (npl == 0) {
        if (ipmt == 0) {
          outPARAM << lladhodo_LCoeff[npl][ipmt];
        } else {
          outPARAM << setw(26) << lladhodo_LCoeff[npl][ipmt];
        }
      } else {
        outPARAM << ", " << setw(15) << lladhodo_LCoeff[npl][ipmt];
      }
    }
    outPARAM << fixed << endl;
  }
  cout << "FINISHED Fitting Hodo Matrix . . . " << endl;
  cout << "Parameter File Created: " << Form("../../PARAM/LAD/HODO/ladhodo_Vpcalib_%d.param", runNUM) << endl;

  // Write Histograms to ROOT file
  // Create output root file where histograms will be stored
  TFile *outROOT = new TFile(Form("HodoCalibPlots_%d.root", runNUM), "recreate");
  // Create directories for each type of canvas
  outROOT->mkdir("TWAvg");
  outROOT->mkdir("TWAvg2D");
  outROOT->mkdir("TWDiff");
  outROOT->mkdir("TWCorr_BarLCoef");
  outROOT->mkdir("TWUnCorr");
  outROOT->mkdir("TWCorr");
  outROOT->mkdir("Param");

  // Loop through all canvases and save each canvas in the appropriate directory
  for (Int_t npl = 0; npl < NPLANES; npl++) {
    outROOT->cd("TWAvg");
    TWAvg_canv[npl]->Write();
    outROOT->cd("TWAvg2D");
    TWAvg_canv_2D[npl]->Write();
    outROOT->cd("TWDiff");
    TWDiff_canv[npl]->Write();
    outROOT->cd("TWCorr_BarLCoef");
    TWCorr_BarLCoef[npl]->Write();
    for (Int_t side = 0; side < SIDES; side++) {
      outROOT->cd("TWUnCorr");
      TWUnCorr_canv[npl][side]->Write();
      outROOT->cd("TWCorr");
      TWCorr_canv[npl][side]->Write();
    }
  }
  outROOT->cd("Param");
  cableArr_canv->Write();
  velArr_canv->Write();
  LCoeff_canv->Write();
  outROOT->Write();
  outROOT->Close();
}
