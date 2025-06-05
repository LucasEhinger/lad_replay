

#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include <string>

using namespace std;

const int nPlanes                     = 6;
const std::string planeNames[nPlanes] = {"000", "001", "100", "101", "200", "REFBAR"};
const int nSides                      = 2;
const int maxPaddles[nPlanes]         = {11, 11, 11, 11, 11, 1};

const double Tdc2Time = 0.09766;

struct hist_param {
  int nBins;
  double xMin;
  double xMax;
};

hist_param tdc_time = {100, 460, 490}; // Histogram parameters for TDC time

int check_LAD_timing() {
  gROOT->SetBatch(kTRUE);
  string inputFileName =
      "/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/COSMICS/LAD_wGEM_cosmic_hall_22831_-1.root";
  TFile *file = TFile::Open(inputFileName.c_str(), "READ");
  if (!file || file->IsZombie()) {
    printf("Error opening file\n");
    return 1;
  }

  TTree *T = (TTree *)file->Get("T");
  if (!T) {
    printf("Error: Tree 'T' not found in file\n");
    file->Close();
    return 1;
  }

  // Define branches
  Double_t TdcTimeCorr[nPlanes][nSides][maxPaddles[0]];
  Double_t RefTime;
  Double_t PhotodiodeTimeRaw;

  // Set branch addresses
  T->SetBranchAddress("L.ladhod.000.TopTdcRefTime", &RefTime);
  T->SetBranchAddress("T.shms.photodiodeLAD_tdcTimeRaw", &PhotodiodeTimeRaw);
  for (int npl = 0; npl < nPlanes; npl++) {
    for (int side = 0; side < nSides; side++) {
      TString branchName = "L.ladhod." + planeNames[npl] + "." + (side == 0 ? "GoodTop" : "GoodBtm") + "TdcTimeCorr";
      T->SetBranchAddress(branchName, TdcTimeCorr[npl][side]);
    }
  }

  // Define Histograms
  TH1F *hTdcTimeCorr[nPlanes][nSides][maxPaddles[0]];
  TH2F *hTdcTimeCorr2D[nPlanes][nSides];

  for (int npl = 0; npl < nPlanes; npl++) {
    for (int side = 0; side < nSides; side++) {
      for (int paddle = 0; paddle < maxPaddles[npl]; paddle++) {
        TString histName = Form("hTdcTimeCorr_%s_%s_%d", planeNames[npl].c_str(), (side == 0 ? "Top" : "Btm"), paddle);
        hTdcTimeCorr[npl][side][paddle] = new TH1F(histName, histName, tdc_time.nBins, tdc_time.xMin, tdc_time.xMax);
      }
      TString hist2DName        = Form("hTdcTimeCorr2D_%s_%s", planeNames[npl].c_str(), (side == 0 ? "Top" : "Btm"));
      hTdcTimeCorr2D[npl][side] = new TH2F(hist2DName, hist2DName, maxPaddles[npl], -0.5, maxPaddles[npl] - 0.5,
                                           tdc_time.nBins, tdc_time.xMin, tdc_time.xMax);
    }
  }

  // Loop over events
  Long64_t nEntries = T->GetEntries();
  for (Long64_t i = 0; i < nEntries; i++) {
    T->GetEntry(i);

    // Fill histograms
    for (int npl = 0; npl < nPlanes; npl++) {
      for (int side = 0; side < nSides; side++) {
        for (int paddle = 0; paddle < maxPaddles[npl]; paddle++) {
          if (TdcTimeCorr[npl][side][paddle] > 0 &&
              TdcTimeCorr[npl][side][paddle] < 10000000) { // Check for valid TDC time
            double correctedTime = TdcTimeCorr[npl][side][paddle] - (PhotodiodeTimeRaw - RefTime) * Tdc2Time;
            hTdcTimeCorr[npl][side][paddle]->Fill(correctedTime);
            hTdcTimeCorr2D[npl][side]->Fill(paddle, correctedTime);
          }
        }
      }
    }

    if (i % (nEntries / 100 == 0 ? 1 : nEntries / 100) == 0) {
      printf("\rProcessing event %lld / %lld (%.1f%%)", i, nEntries, 100.0 * i / nEntries);
      fflush(stdout);
    }
  }

  TString outFileName = "lad_timing_histograms.root";
  TFile *outFile      = new TFile(outFileName, "RECREATE");
  TDirectory *dir1D = outFile->mkdir("1D_Histograms");
  TDirectory *dir2D = outFile->mkdir("2D_Histograms");

  for (int npl = 0; npl < nPlanes; npl++) {
    for (int side = 0; side < nSides; side++) {
      TString canvasName = Form("cPlane_%s_%s", planeNames[npl].c_str(), (side == 0 ? "Top" : "Btm"));
      int npaddles       = maxPaddles[npl];
      TCanvas *c         = new TCanvas(canvasName, canvasName, 200 * npaddles, 400);
      if (npaddles > 1) {
        c->Divide(3, 4);
      }

      for (int paddle = 0; paddle < npaddles; paddle++) {
        c->cd(paddle + 1);
        gPad->SetLogy();
        hTdcTimeCorr[npl][side][paddle]->Draw();
      }
      dir1D->cd();
      c->Write();
      dir2D->cd();
      hTdcTimeCorr2D[npl][side]->Write();
    }
  }

  outFile->Write();
  outFile->Close();

  file->Close();
  return 0;
}