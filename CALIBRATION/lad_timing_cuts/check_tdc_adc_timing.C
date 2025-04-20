#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TTree.h>
#include <iostream>

const static int NPLANES  = 6;
const static int NBARS    = 11;
const static int MAX_DATA = 105;

const static double adc2ns       = 0.0625;  // ns
const static double tdc2ns       = 0.09766; // ns
const static int NBINS_DIFF      = 500;
const static int DIFF_MIN        = -3087;
const static int DIFF_MAX        = -3070; // ns
const static int NBINS_DIFF_RAW  = 50;
const static int DIFF_RAW_MIN    = 1560;
const static int DIFF_RAW_MAX    = 1620; // ns
const static int NBINS_DIFF_GOOD = 100;
const static int DIFF_GOOD_MIN   = -20; // ns
const static int DIFF_GOOD_MAX   = 50;  // ns
const static int NBINS_AMP       = 100;
const static int AMP_MIN         = 0;
const static int AMP_MAX         = 1000; // mV

const char *spec_prefix  = "H";
const static int nSides  = 2;
const char *side[nSides] = {"Top", "Btm"};

const static int kBig = 1000000;

const static char *plane_names[NPLANES] = {"000", "001", "100", "101", "200", "REFBAR"};
void check_tdc_adc_timing(TString inputFileName = "ROOTfiles/CALIB/LAD_COIN_22300_300000.root",
                          TString outfileName = "move_me.root") {
  // Set batch mode to suppress graphical output
  gROOT->SetBatch(kTRUE);

  TFile *file           = TFile::Open(inputFileName, "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "Error: Cannot open the ROOT file!" << std::endl;
    return;
  }

  // Get the tree
  TTree *tree = dynamic_cast<TTree *>(file->Get("T"));
  if (!tree) {
    std::cerr << "Error: Cannot find the tree 'T' in the file!" << std::endl;
    file->Close();
    return;
  }

  // Create histograms
  TH1D hist_adc_tdc_diff[nSides][NPLANES][NBARS];
  TH1D hist_adc_tdc_diff_raw[nSides][NPLANES][NBARS];
  TH1D hist_adc_tdc_diff_TWUnCorr[nSides][NPLANES][NBARS];
  TH1D hist_adc_tdc_diff_Corr[nSides][NPLANES][NBARS];
  TH1D hist_adc_tdc_diff_TWCorr[nSides][NPLANES][NBARS];

  TH2D hist_adc_tdc_diff_amp[nSides][NPLANES][NBARS];
  TH2D hist_adc_tdc_diff_amp_raw[nSides][NPLANES][NBARS];
  TH2D hist_adc_tdc_diff_TWCorr_amp[nSides][NPLANES][NBARS];
  TH2D hist_adc_tdc_diff_TWUnCorr_amp[nSides][NPLANES][NBARS];
  TH2D hist_adc_tdc_diff_Corr_amp[nSides][NPLANES][NBARS];

  for (int plane = 0; plane < NPLANES; ++plane) {
    for (int i = 0; i < NBARS; ++i) {
      for (int side_idx = 0; side_idx < nSides; ++side_idx) {
        hist_adc_tdc_diff[side_idx][plane][i] = TH1D(
            Form("hist_adc_tdc_diff_%s_plane%d_bar%d", side[side_idx], plane, i),
            Form("Plane %d %s: TdcTime - AdcPulseTime (Bar %d);Time Difference (ns);Entries", plane, side[side_idx], i),
            NBINS_DIFF, DIFF_MIN, DIFF_MAX);
        hist_adc_tdc_diff_raw[side_idx][plane][i] =
            TH1D(Form("hist_diff_raw_%s_plane%d_bar%d", side[side_idx], plane, i),
                 Form("Plane %d %s: TdcTimeRaw - AdcPulseTimeRaw (Bar %d);Time Difference (ns);Entries", plane,
                      side[side_idx], i),
                 NBINS_DIFF_RAW, DIFF_RAW_MIN, DIFF_RAW_MAX);
        hist_adc_tdc_diff_amp[side_idx][plane][i] =
            TH2D(Form("hist_adc_tdc_diff_amp_%s_plane%d_bar%d", side[side_idx], plane, i),
                 Form("Plane %d %s: TdcTime - AdcPulseTime (Bar %d);ADC Amplitude (mV);Time Difference (ns)", plane,
                      side[side_idx], i),
                 NBINS_AMP, AMP_MIN, AMP_MAX, NBINS_DIFF, DIFF_MIN, DIFF_MAX);
        hist_adc_tdc_diff_amp_raw[side_idx][plane][i] =
            TH2D(Form("hist_adc_tdc_diff_amp_raw_%s_plane%d_bar%d", side[side_idx], plane, i),
                 Form("Plane %d %s: TdcTimeRaw - AdcPulseTimeRaw (Bar %d);ADC Amplitude (mV);Time Difference (ns)",
                      plane, side[side_idx], i),
                 NBINS_AMP, AMP_MIN, AMP_MAX, NBINS_DIFF_RAW, DIFF_RAW_MIN, DIFF_RAW_MAX);
        hist_adc_tdc_diff_TWCorr_amp[side_idx][plane][i] = TH2D(
            Form("hist_tdc_TWCorr_%s_plane%d_bar%d", side[side_idx], plane, i),
            Form("Plane %d %s: TdcTimeWalkCorr (Bar %d);Amplitude (mV);Time Difference (ns)", plane, side[side_idx], i),
            NBINS_AMP, AMP_MIN, AMP_MAX, NBINS_DIFF_GOOD, DIFF_GOOD_MIN, DIFF_GOOD_MAX);
        hist_adc_tdc_diff_TWUnCorr_amp[side_idx][plane][i] = TH2D(
            Form("hist_tdc_diff_TWUnCorr_%s_plane%d_bar%d", side[side_idx], plane, i),
            Form("Plane %d %s: TdcTimeUnCorr (Bar %d);Amplitude (mV);Time Difference (ns)", plane, side[side_idx], i),
            NBINS_AMP, AMP_MIN, AMP_MAX, NBINS_DIFF_GOOD, DIFF_GOOD_MIN, DIFF_GOOD_MAX);
        hist_adc_tdc_diff_Corr_amp[side_idx][plane][i] =
            TH2D(Form("hist_tdc_diff_Corr_%s_plane%d_bar%d", side[side_idx], plane, i),
                 Form("Plane %d %s: TdcTimeCorr - AdcPulseTime (Bar %d);Amplitude (mV);Time Difference (ns)", plane,
                      side[side_idx], i),
                 NBINS_AMP, AMP_MIN, AMP_MAX, NBINS_DIFF_GOOD, DIFF_GOOD_MIN, DIFF_GOOD_MAX);

        hist_adc_tdc_diff_TWUnCorr[side_idx][plane][i] =
            TH1D(Form("hist_adc_tdc_diff_TWUnCorr_%s_plane%d_bar%d", side[side_idx], plane, i),
                 Form("Plane %d %s: TdcTimeUnCorr - AdcPulseTime (Bar %d);Time Difference (ns);Entries", plane,
                      side[side_idx], i),
                 NBINS_DIFF_GOOD, DIFF_GOOD_MIN, DIFF_GOOD_MAX);
        hist_adc_tdc_diff_Corr[side_idx][plane][i] =
            TH1D(Form("hist_adc_tdc_diff_Corr_%s_plane%d_bar%d", side[side_idx], plane, i),
                 Form("Plane %d %s: TdcTimeCorr - AdcPulseTime (Bar %d);Time Difference (ns);Entries", plane,
                      side[side_idx], i),
                 NBINS_DIFF_GOOD, DIFF_GOOD_MIN, DIFF_GOOD_MAX);
        hist_adc_tdc_diff_TWCorr[side_idx][plane][i] =
            TH1D(Form("hist_adc_tdc_diff_TWCorr_%s_plane%d_bar%d", side[side_idx], plane, i),
                 Form("Plane %d %s: TdcTimeWalkCorr - AdcPulseTime (Bar %d);Time Difference (ns);Entries", plane,
                      side[side_idx], i),
                 NBINS_DIFF_GOOD, DIFF_GOOD_MIN, DIFF_GOOD_MAX);
      }
    }
  }

  // Variables for branches
  Double_t GoodAdcTdcDiffTime[nSides][NPLANES][MAX_DATA]  = {0};
  Double_t GoodTdcTimeUnCorr[nSides][NPLANES][MAX_DATA]   = {0};
  Double_t GoodTdcTimeWalkCorr[nSides][NPLANES][MAX_DATA] = {0};
  Double_t GoodTdcTimeCorr[nSides][NPLANES][MAX_DATA]     = {0};
  Double_t GoodAdcPulseAmp[nSides][NPLANES][MAX_DATA]     = {0};
  Double_t GoodAdcPulseTime[nSides][NPLANES][MAX_DATA]    = {0};
  Double_t AdcPulseTimeRaw[nSides][NPLANES][MAX_DATA]     = {0};
  Double_t AdcPulseTime[nSides][NPLANES][MAX_DATA]        = {0};
  Double_t AdcPulseAmp[nSides][NPLANES][MAX_DATA]         = {0};
  Double_t TdcTimeRaw[nSides][NPLANES][MAX_DATA]          = {0};
  Double_t TdcTime[nSides][NPLANES][MAX_DATA]             = {0};
  Double_t AdcCounter[nSides][NPLANES][MAX_DATA]          = {0};
  Double_t TdcCounter[nSides][NPLANES][MAX_DATA]          = {0};
  Int_t NDataAdc[nSides][NPLANES]                         = {0};
  Int_t NDataTdc[nSides][NPLANES]                         = {0};

  // Set branch addresses
  // tree->SetBranchAddress("H.ladhod.000.GoodTopAdcTdcDiffTime", GoodTopAdcTdcDiffTime);
  for (int i_side = 0; i_side < nSides; i_side++) {
    for (int i_plane = 0; i_plane < NPLANES; i_plane++) {
      tree->SetBranchAddress(Form("%s.ladhod.%s.%sAdcPulseTimeRaw", spec_prefix, plane_names[i_plane], side[i_side]),
                             AdcPulseTimeRaw[i_side][i_plane]);
      tree->SetBranchAddress(Form("%s.ladhod.%s.%sAdcPulseTime", spec_prefix, plane_names[i_plane], side[i_side]),
                             AdcPulseTime[i_side][i_plane]);
      tree->SetBranchAddress(Form("%s.ladhod.%s.%sAdcPulseAmp", spec_prefix, plane_names[i_plane], side[i_side]),
                             AdcPulseAmp[i_side][i_plane]);
      tree->SetBranchAddress(Form("%s.ladhod.%s.%sTdcTimeRaw", spec_prefix, plane_names[i_plane], side[i_side]),
                             TdcTimeRaw[i_side][i_plane]);
      tree->SetBranchAddress(Form("%s.ladhod.%s.%sTdcTime", spec_prefix, plane_names[i_plane], side[i_side]),
                             TdcTime[i_side][i_plane]);
      tree->SetBranchAddress(Form("%s.ladhod.%s.%sAdcCounter", spec_prefix, plane_names[i_plane], side[i_side]),
                             AdcCounter[i_side][i_plane]);
      tree->SetBranchAddress(Form("%s.ladhod.%s.%sTdcCounter", spec_prefix, plane_names[i_plane], side[i_side]),
                             TdcCounter[i_side][i_plane]);
      tree->SetBranchAddress(Form("Ndata.%s.ladhod.%s.%sAdcPulseTime", spec_prefix, plane_names[i_plane], side[i_side]),
                             &NDataAdc[i_side][i_plane]);
      tree->SetBranchAddress(Form("Ndata.%s.ladhod.%s.%sTdcTime", spec_prefix, plane_names[i_plane], side[i_side]),
                             &NDataTdc[i_side][i_plane]);
      tree->SetBranchAddress(Form("%s.ladhod.%s.Good%sTdcTimeUnCorr", spec_prefix, plane_names[i_plane], side[i_side]),
                             GoodTdcTimeUnCorr[i_side][i_plane]);
      tree->SetBranchAddress(
          Form("%s.ladhod.%s.Good%sTdcTimeWalkCorr", spec_prefix, plane_names[i_plane], side[i_side]),
          GoodTdcTimeWalkCorr[i_side][i_plane]);
      tree->SetBranchAddress(Form("%s.ladhod.%s.Good%sTdcTimeCorr", spec_prefix, plane_names[i_plane], side[i_side]),
                             GoodTdcTimeCorr[i_side][i_plane]);
      tree->SetBranchAddress(Form("%s.ladhod.%s.Good%sAdcTdcDiffTime", spec_prefix, plane_names[i_plane], side[i_side]),
                             GoodAdcTdcDiffTime[i_side][i_plane]);
      tree->SetBranchAddress(Form("%s.ladhod.%s.Good%sAdcPulseAmp", spec_prefix, plane_names[i_plane], side[i_side]),
                             GoodAdcPulseAmp[i_side][i_plane]);
      tree->SetBranchAddress(Form("%s.ladhod.%s.Good%sAdcPulseTime", spec_prefix, plane_names[i_plane], side[i_side]),
                             GoodAdcPulseTime[i_side][i_plane]);
    }
  }

  // Loop over entries
  Long64_t nEntries = tree->GetEntries();
  for (Long64_t i = 0; i < nEntries; ++i) {
    // for (Long64_t i = 0; i < 5000; ++i) {

    tree->GetEntry(i);
    for (int i_side = 0; i_side < nSides; i_side++) {
      for (int i_plane = 0; i_plane < NPLANES; i_plane++) {
        // Check if the number of ADC and TDC data exceeds MAX_DATA
        if (NDataAdc[i_side][i_plane] > MAX_DATA || NDataTdc[i_side][i_plane] > MAX_DATA) {
          std::cerr << "Error: NDataAdc or NDataTdc exceeds MAX_DATA. Skipping entry " << i << "." << std::endl;
          continue;
        }

        // Fill GoodTDCTimes - Good ADC Time
        for (int bar = 0; bar < NBARS; ++bar) {
          if (GoodAdcPulseTime[i_side][i_plane][bar] < kBig) {
            hist_adc_tdc_diff_TWUnCorr[i_side][i_plane][bar].Fill(GoodTdcTimeUnCorr[i_side][i_plane][bar] -
                                                                  GoodAdcPulseTime[i_side][i_plane][bar]);
            hist_adc_tdc_diff_Corr[i_side][i_plane][bar].Fill(GoodTdcTimeCorr[i_side][i_plane][bar] -
                                                              GoodAdcPulseTime[i_side][i_plane][bar]);

            hist_adc_tdc_diff_TWCorr[i_side][i_plane][bar].Fill(GoodTdcTimeWalkCorr[i_side][i_plane][bar] -
                                                                GoodAdcPulseTime[i_side][i_plane][bar]);
            hist_adc_tdc_diff_TWCorr_amp[i_side][i_plane][bar].Fill(GoodAdcPulseAmp[i_side][i_plane][bar],
                                                                    GoodTdcTimeWalkCorr[i_side][i_plane][bar] -
                                                                        GoodAdcPulseTime[i_side][i_plane][bar]);
            hist_adc_tdc_diff_TWUnCorr_amp[i_side][i_plane][bar].Fill(GoodAdcPulseAmp[i_side][i_plane][bar],
                                                                      GoodTdcTimeUnCorr[i_side][i_plane][bar] -
                                                                          GoodAdcPulseTime[i_side][i_plane][bar]);
            hist_adc_tdc_diff_Corr_amp[i_side][i_plane][bar].Fill(GoodAdcPulseAmp[i_side][i_plane][bar],
                                                                  GoodTdcTimeCorr[i_side][i_plane][bar] -
                                                                      GoodAdcPulseTime[i_side][i_plane][bar]);
          }
        }
        // Loop through Top ADC and TDC data
        for (int j = 0; j < NDataAdc[i_side][i_plane]; ++j) {
          for (int k = 0; k < NDataTdc[i_side][i_plane]; ++k) {
            // Check if the bar indices match
            int bar_tdc = TdcCounter[i_side][i_plane][k];
            int bar_adc = AdcCounter[i_side][i_plane][j];
            if (bar_tdc != bar_adc) {
              continue;
            }
            // Calculate the difference between AdcPulseTimeRaw and TdcTimeRaw
            double diff_raw = TdcTimeRaw[i_side][i_plane][k] * tdc2ns - AdcPulseTimeRaw[i_side][i_plane][j] * adc2ns;
            double diff     = TdcTime[i_side][i_plane][k] * tdc2ns - AdcPulseTime[i_side][i_plane][j];
            // Fill the corresponding histogram based on the bar index
            int bar_index = AdcCounter[i_side][i_plane][j]; // Assuming bar index is stored in AdcCounter

            hist_adc_tdc_diff[i_side][i_plane][bar_index - 1].Fill(diff);
            hist_adc_tdc_diff_raw[i_side][i_plane][bar_index - 1].Fill(diff_raw);
            hist_adc_tdc_diff_amp[i_side][i_plane][bar_index - 1].Fill(AdcPulseAmp[i_side][i_plane][j], diff);
            hist_adc_tdc_diff_amp_raw[i_side][i_plane][bar_index - 1].Fill(AdcPulseAmp[i_side][i_plane][j], diff_raw);
          }
        }
      }
    }
    //   // Display progress
    std::cout << "\rProcessing: " << std::fixed << std::setprecision(2) << (100.0 * (i + 1) / nEntries) << "% completed"
              << std::flush;
  }

  std::cout << std::endl;
  std::cout << "Processing completed." << std::endl;
  // Save histograms to a ROOT file

  TFile *outputFile = new TFile(outfileName, "RECREATE");
  if (!outputFile || outputFile->IsZombie()) {
    std::cerr << "Error: Cannot create the output ROOT file!" << std::endl;
    return;
  }
  // Create directories for each histogram type
  TDirectory *dir_diff              = outputFile->mkdir("adc_tdc_diff");
  TDirectory *dir_diff_raw          = outputFile->mkdir("adc_tdc_diff_raw");
  TDirectory *dir_diff_amp          = outputFile->mkdir("adc_tdc_diff_amp");
  TDirectory *dir_diff_amp_raw      = outputFile->mkdir("adc_tdc_diff_amp_raw");
  TDirectory *dir_diff_TWComp       = outputFile->mkdir("adc_tdc_diff_TWComp");
  TDirectory *dir_TWCorr_amp        = outputFile->mkdir("adc_tdc_TWCorr_amp");
  TDirectory *dir_diff_TWUnCorr_amp = outputFile->mkdir("adc_tdc_diff_TWUnCorr_amp");
  TDirectory *dir_diff_Corr_amp     = outputFile->mkdir("adc_tdc_diff_Corr_amp");

  for (int plane = 0; plane < NPLANES; ++plane) {
    // Create subdirectories for each plane
    for (int side_idx = 0; side_idx < nSides; ++side_idx) {
      TCanvas *canvas_diff         = new TCanvas(Form("canvas_diff_plane_%d_%s", plane, side[side_idx]),
                                                 Form("ADC-TDC Diff %s", side[side_idx]), 1200, 800);
      TCanvas *canvas_diff_raw     = new TCanvas(Form("canvas_diff_raw_plane_%d_%s", plane, side[side_idx]),
                                                 Form("ADC-TDC Diff Raw %s", side[side_idx]), 1200, 800);
      TCanvas *canvas_diff_amp     = new TCanvas(Form("canvas_diff_amp_plane_%d_%s", plane, side[side_idx]),
                                                 Form("ADC-TDC Diff Amp %s", side[side_idx]), 1200, 800);
      TCanvas *canvas_diff_amp_raw = new TCanvas(Form("canvas_diff_amp_raw_plane_%d_%s", plane, side[side_idx]),
                                                 Form("ADC-TDC Diff Amp Raw %s", side[side_idx]), 1200, 800);
      TCanvas *canvas_diff_TW_Comp = new TCanvas(Form("canvas_diff_TWComp_plane_%d_%s", plane, side[side_idx]),
                                                 Form("ADC-TDC Diff %s", side[side_idx]), 1200, 800);
      TCanvas *canvas_TWCorr_amp   = new TCanvas(Form("canvas_TWCorr_amp_plane_%d_%s", plane, side[side_idx]),
                                                 Form("ADC-TDC TWCorr Amp %s", side[side_idx]), 1200, 800);
      TCanvas *canvas_diff_TWUnCorr_amp =
          new TCanvas(Form("canvas_diff_TWUnCorr_amp_plane_%d_%s", plane, side[side_idx]),
                      Form("ADC-TDC Diff TWUnCorr Amp %s", side[side_idx]), 1200, 800);
      TCanvas *canvas_diff_Corr_amp = new TCanvas(Form("canvas_diff_Corr_amp_plane_%d_%s", plane, side[side_idx]),
                                                  Form("ADC-TDC Diff Corr Amp %s", side[side_idx]), 1200, 800);

      int nCols = 4;                           // Number of columns for subplots
      int nRows = (NBARS + nCols - 1) / nCols; // Calculate rows based on bars and columns

      canvas_diff->Divide(nCols, nRows);
      canvas_diff_raw->Divide(nCols, nRows);
      canvas_diff_amp->Divide(nCols, nRows);
      canvas_diff_amp_raw->Divide(nCols, nRows);
      canvas_diff_TW_Comp->Divide(nCols, nRows);
      canvas_TWCorr_amp->Divide(nCols, nRows);
      canvas_diff_TWUnCorr_amp->Divide(nCols, nRows);
      canvas_diff_Corr_amp->Divide(nCols, nRows);

      for (int bar = 0; bar < NBARS; ++bar) {
        canvas_diff->cd(bar + 1);
        hist_adc_tdc_diff[side_idx][plane][bar].Draw();

        canvas_diff_raw->cd(bar + 1);
        hist_adc_tdc_diff_raw[side_idx][plane][bar].Draw();

        canvas_diff_amp->cd(bar + 1);
        hist_adc_tdc_diff_amp[side_idx][plane][bar].Draw("COLZ");

        canvas_diff_amp_raw->cd(bar + 1);
        hist_adc_tdc_diff_amp_raw[side_idx][plane][bar].Draw("COLZ");

        canvas_diff_TW_Comp->cd(bar + 1);
        hist_adc_tdc_diff_TWUnCorr[side_idx][plane][bar].SetLineColor(kRed);
        hist_adc_tdc_diff_TWUnCorr[side_idx][plane][bar].Draw();

        hist_adc_tdc_diff_Corr[side_idx][plane][bar].SetLineColor(kBlue);
        hist_adc_tdc_diff_Corr[side_idx][plane][bar].Draw("SAME");

        hist_adc_tdc_diff_TWCorr[side_idx][plane][bar].SetLineColor(kGreen);
        hist_adc_tdc_diff_TWCorr[side_idx][plane][bar].Draw("SAME");

        // Adjust y-axis range to fit all graphs
        double max_y = std::max({hist_adc_tdc_diff_TWUnCorr[side_idx][plane][bar].GetMaximum(),
                                 hist_adc_tdc_diff_Corr[side_idx][plane][bar].GetMaximum(),
                                 hist_adc_tdc_diff_TWCorr[side_idx][plane][bar].GetMaximum()});
        hist_adc_tdc_diff_TWUnCorr[side_idx][plane][bar].SetMaximum(max_y * 1.1);

        // Add a legend
        TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
        legend->AddEntry(&hist_adc_tdc_diff_TWUnCorr[side_idx][plane][bar], "TWUnCorr", "l");
        legend->AddEntry(&hist_adc_tdc_diff_Corr[side_idx][plane][bar], "Corr", "l");
        legend->AddEntry(&hist_adc_tdc_diff_TWCorr[side_idx][plane][bar], "TWCorr", "l");
        legend->Draw();

        canvas_TWCorr_amp->cd(bar + 1);
        hist_adc_tdc_diff_TWCorr_amp[side_idx][plane][bar].Draw("COLZ");

        canvas_diff_TWUnCorr_amp->cd(bar + 1);
        hist_adc_tdc_diff_TWUnCorr_amp[side_idx][plane][bar].Draw("COLZ");

        canvas_diff_Corr_amp->cd(bar + 1);
        hist_adc_tdc_diff_Corr_amp[side_idx][plane][bar].Draw("COLZ");
      }

      // Write canvases to the output file
      dir_diff->cd();
      canvas_diff->Write();
      dir_diff_raw->cd();
      canvas_diff_raw->Write();
      dir_diff_amp->cd();
      canvas_diff_amp->Write();
      dir_diff_amp_raw->cd();
      canvas_diff_amp_raw->Write();
      dir_diff_TWComp->cd();
      canvas_diff_TW_Comp->Write();
      dir_TWCorr_amp->cd();
      canvas_TWCorr_amp->Write();
      dir_diff_TWUnCorr_amp->cd();
      canvas_diff_TWUnCorr_amp->Write();
      dir_diff_Corr_amp->cd();
      canvas_diff_Corr_amp->Write();

      // Clean up
      delete canvas_diff;
      delete canvas_diff_raw;
      delete canvas_diff_amp;
      delete canvas_diff_amp_raw;
      delete canvas_diff_TW_Comp;
      delete canvas_TWCorr_amp;
      delete canvas_diff_TWUnCorr_amp;
      delete canvas_diff_Corr_amp;
    }
  }

  outputFile->Close();

  // // Clean up
  // file->Close();
  // delete file;
  // for (int i = 0; i < 11; ++i) {
  //   delete hist_adc_tdc_diff[i];
  //   delete hist_adc_tdc_diff_raw[i];
  //   delete hist_adc_tdc_diff_amp[i];
  //   delete hist_adc_tdc_diff_amp_raw[i];
  // }
}