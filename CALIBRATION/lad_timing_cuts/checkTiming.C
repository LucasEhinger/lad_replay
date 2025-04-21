// Macro to perform time-walk fits and extract the calibration parameters
// Author: Eric Pooser, pooser@jlab.org
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"
#include "TTree.h"
#include <TF1.h>
#include <TGraph.h>
#include <TH2.h>
#include <TLegend.h>
#include <TMath.h>
#include <TMultiGraph.h>
#include <TNtuple.h>
#include <TObjArray.h>
#include <TPaveLabel.h>
#include <TPaveText.h>
#include <TPolyLine.h>
#include <TProfile.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <fstream>
#include <iostream>

// Declare ROOT files
TFile *histoFile, *outFile;

// Declare constants
static const UInt_t nPlanes  = 6;
static const UInt_t nSides   = 2;
static const UInt_t nBarsMax = 11;

static const Double_t tdcThresh      = 120.0; // 30 mV in units of FADC channels
static const Double_t twFitRangeLow  = 20.0;
static const Double_t twFitRangeHigh = 600.0;
static const Double_t c0twParInit    = 1.0;
static const Double_t c1twParInit    = 1.0;

// Parameter values to be written to param file
Double_t c1[nPlanes][nSides][nBarsMax] = {0.};
Double_t c2[nPlanes][nSides][nBarsMax] = {0.};

static const Double_t fontSize     = 0.05;
static const Double_t yTitleOffset = 0.75;
static const Double_t markerSize   = 2.0;
static const Double_t minScale     = 0.75;
static const Double_t maxScale     = 0.75;

static const UInt_t lineWidth = 4;
static const UInt_t lineStyle = 7;

static const UInt_t nbars[nPlanes]       = {11, 11, 11, 11, 11, 1};
static const TString planeNames[nPlanes] = {"000", "001", "100", "101", "200", "REFBAR"};
static const TString sideNames[nSides]   = {"Top", "Btm"};

static const int NBINS_TDC = 100;
static const int TDC_MIN   = -180;
static const int TDC_MAX   = -120;

void checkTiming(const TString &inputFileName, const TString &outputFileName) {
  // Set ROOT to batch mode
  gROOT->SetBatch(kTRUE);
  // Open the ROOT file
  TFile *file = TFile::Open(inputFileName);
  if (!file || file->IsZombie()) {
    std::cerr << "Error: Cannot open file " << inputFileName << std::endl;
    return;
  }

  // Access the TTree
  TTree *T = (TTree *)file->Get("T");
  if (!T) {
    std::cerr << "Error: Cannot find TTree 'T' in file " << inputFileName << std::endl;
    file->Close();
    return;
  }

  // Create a new ROOT file to save the histograms
  TFile *outputFile = new TFile(outputFileName, "RECREATE");
  if (!outputFile || outputFile->IsZombie()) {
    std::cerr << "Error: Cannot create output file '" << outputFileName << "'" << std::endl;
    file->Close();
    return;
  }

  // Create directories for Corrected vs Uncorrected and subdirectories for each plane
  outputFile->mkdir("Corrected_vs_Uncorrected");
  // Create directories for 2D Corrected and 2D Uncorrected
  outputFile->mkdir("2D_Corrected");
  outputFile->mkdir("2D_Uncorrected");

  // Loop over planes and bars
  for (UInt_t plane = 0; plane < nPlanes; ++plane) {
    std::vector<TH2F *> hist2DCorrBarsTop;
    std::vector<TH2F *> hist2DUnCorrBarsTop;
    std::vector<TH2F *> hist2DCorrBarsBtm;
    std::vector<TH2F *> hist2DUnCorrBarsBtm;

    for (UInt_t bar = 0; bar < nbars[plane]; ++bar) {

      ////////////////////////////////////////////////////////////////////////////////////////////
      // Create a canvas for the current plane if it doesn't exist
      TCanvas *cTopPlane = (TCanvas *)gROOT->FindObject(Form("cTopPlane_%d", plane));
      if (!cTopPlane) {
        cTopPlane = new TCanvas(Form("cTopPlane_%d", plane), Form("Top Plane %d", plane), 1200, 800);
        cTopPlane->Divide(4, 3); // Divide canvas into subpads for each bar
      }

      // Draw corrected and uncorrected times on the same histogram for Top
      TString histNameTopCorr   = Form("GoodTopTdcTimeCorr_plane%d_bar%d", plane, bar);
      TString histNameTopUnCorr = Form("GoodTopTdcTimeUnCorr_plane%d_bar%d", plane, bar);
      TString drawCmdTopCorr = Form("L.ladhod.%s.GoodTopTdcTimeCorr[%d]>>%s(%d, %d, %d)", planeNames[plane].Data(), bar,
                histNameTopCorr.Data(), NBINS_TDC, TDC_MIN, TDC_MAX);
      TString drawCmdTopUnCorr = Form("L.ladhod.%s.GoodTopTdcTimeUnCorr[%d]>>%s(%d, %d, %d)", planeNames[plane].Data(),
                  bar, histNameTopUnCorr.Data(), NBINS_TDC, TDC_MIN, TDC_MAX);

      T->Draw(drawCmdTopCorr, "", "goff");
      T->Draw(drawCmdTopUnCorr, "", "goff");

      TH1F *histTopCorr   = (TH1F *)gROOT->FindObject(histNameTopCorr);
      TH1F *histTopUnCorr = (TH1F *)gROOT->FindObject(histNameTopUnCorr);

      if (histTopCorr && histTopUnCorr) {
        cTopPlane->cd(bar + 1); // Select the subpad corresponding to the current bar
        histTopCorr->SetLineColor(kRed);
        histTopUnCorr->SetLineColor(kBlue);
        histTopCorr->SetTitle(Form("Plane %d, Bar %d", plane, bar));
        histTopCorr->GetXaxis()->SetTitle("Time (ns)");
        histTopCorr->GetYaxis()->SetTitle("Counts");
        histTopCorr->Draw();
        histTopUnCorr->Draw("SAME");

        TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
        legend->AddEntry(histTopCorr, "Corrected", "l");
        legend->AddEntry(histTopUnCorr, "Uncorrected", "l");
        legend->Draw();
      }

      // Save the canvas for the current plane after processing all bars
      if (bar == nbars[plane] - 1) {
        outputFile->cd("Corrected_vs_Uncorrected");
        cTopPlane->Write();
      }

      // Draw corrected and uncorrected times for all bars on the same canvas for Btm
      TCanvas *cBtmPlane = (TCanvas *)gROOT->FindObject(Form("cBtmPlane_%d", plane));
      if (!cBtmPlane) {
        cBtmPlane = new TCanvas(Form("cBtmPlane_%d", plane), Form("Btm Plane %d", plane), 1200, 800);
        cBtmPlane->Divide(4, 3); // Divide canvas into subpads for each bar
      }

      TString histNameBtmCorr   = Form("GoodBtmTdcTimeCorr_plane%d_bar%d", plane, bar);
      TString histNameBtmUnCorr = Form("GoodBtmTdcTimeUnCorr_plane%d_bar%d", plane, bar);
      TString drawCmdBtmCorr = Form("L.ladhod.%s.GoodBtmTdcTimeCorr[%d]>>%s(%d, %d, %d)", planeNames[plane].Data(), bar,
                histNameBtmCorr.Data(), NBINS_TDC, TDC_MIN, TDC_MAX);
      TString drawCmdBtmUnCorr = Form("L.ladhod.%s.GoodBtmTdcTimeUnCorr[%d]>>%s(%d, %d, %d)", planeNames[plane].Data(),
                  bar, histNameBtmUnCorr.Data(), NBINS_TDC, TDC_MIN, TDC_MAX);

      T->Draw(drawCmdBtmCorr, "", "goff");
      T->Draw(drawCmdBtmUnCorr, "", "goff");

      TH1F *histBtmCorr   = (TH1F *)gROOT->FindObject(histNameBtmCorr);
      TH1F *histBtmUnCorr = (TH1F *)gROOT->FindObject(histNameBtmUnCorr);

      if (histBtmCorr && histBtmUnCorr) {
        cBtmPlane->cd(bar + 1); // Select the subpad corresponding to the current bar
        histBtmCorr->SetLineColor(kRed);
        histBtmUnCorr->SetLineColor(kBlue);
        histBtmCorr->SetTitle(Form("Plane %d, Bar %d", plane, bar));
        histBtmCorr->GetXaxis()->SetTitle("Time (ns)");
        histBtmCorr->GetYaxis()->SetTitle("Counts");
        histBtmCorr->Draw();
        histBtmUnCorr->Draw("SAME");

        TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
        legend->AddEntry(histBtmCorr, "Corrected", "l");
        legend->AddEntry(histBtmUnCorr, "Uncorrected", "l");
        legend->Draw();
      }

      // Save the canvas for the current plane after processing all bars
      if (bar == nbars[plane] - 1) {
        outputFile->cd("Corrected_vs_Uncorrected");
        cBtmPlane->Write();
      }
      ////////////////////////////////////////////////////////////////////////////////////////////

      // Create 2D histograms for corrected and uncorrected times for each bar (Top)
      TString histName2DCorrBarTop   = Form("GoodTopTdcTimeCorr_2D_plane%d_bar%d", plane, bar);
      TString histName2DUnCorrBarTop = Form("GoodTopTdcTimeUnCorr_2D_plane%d_bar%d", plane, bar);
      TString drawCmd2DCorrBarTop =
          Form("L.ladhod.%s.GoodTopTdcTimeCorr[%d]:%d>>%s(%d, 0, %d, %d, %d, %d)", planeNames[plane].Data(), bar, bar,
               histName2DCorrBarTop.Data(), nbars[plane], nbars[plane], NBINS_TDC, TDC_MIN, TDC_MAX);
      TString drawCmd2DUnCorrBarTop =
          Form("L.ladhod.%s.GoodTopTdcTimeUnCorr[%d]:%d>>%s(%d, 0, %d, %d, %d, %d)", planeNames[plane].Data(), bar, bar,
               histName2DUnCorrBarTop.Data(), nbars[plane], nbars[plane], NBINS_TDC, TDC_MIN, TDC_MAX);

      T->Draw(drawCmd2DCorrBarTop, "", "goff");
      T->Draw(drawCmd2DUnCorrBarTop, "", "goff");

      TH2F *hist2DCorrBarTop   = (TH2F *)gROOT->FindObject(histName2DCorrBarTop);
      TH2F *hist2DUnCorrBarTop = (TH2F *)gROOT->FindObject(histName2DUnCorrBarTop);

      if (hist2DCorrBarTop) {
        hist2DCorrBarsTop.push_back(hist2DCorrBarTop);
      }
      if (hist2DUnCorrBarTop) {
        hist2DUnCorrBarsTop.push_back(hist2DUnCorrBarTop);
      }

      // Create 2D histograms for corrected and uncorrected times for each bar (Btm)
      TString histName2DCorrBarBtm   = Form("GoodBtmTdcTimeCorr_2D_plane%d_bar%d", plane, bar);
      TString histName2DUnCorrBarBtm = Form("GoodBtmTdcTimeUnCorr_2D_plane%d_bar%d", plane, bar);
      TString drawCmd2DCorrBarBtm =
          Form("L.ladhod.%s.GoodBtmTdcTimeCorr[%d]:%d>>%s(%d, 0, %d, %d, %d, %d)", planeNames[plane].Data(), bar, bar,
               histName2DCorrBarBtm.Data(), nbars[plane], nbars[plane], NBINS_TDC, TDC_MIN, TDC_MAX);
      TString drawCmd2DUnCorrBarBtm =
          Form("L.ladhod.%s.GoodBtmTdcTimeUnCorr[%d]:%d>>%s(%d, 0, %d, %d, %d, %d)", planeNames[plane].Data(), bar, bar,
               histName2DUnCorrBarBtm.Data(), nbars[plane], nbars[plane], NBINS_TDC, TDC_MIN, TDC_MAX);

      T->Draw(drawCmd2DCorrBarBtm, "", "goff");
      T->Draw(drawCmd2DUnCorrBarBtm, "", "goff");

      TH2F *hist2DCorrBarBtm   = (TH2F *)gROOT->FindObject(histName2DCorrBarBtm);
      TH2F *hist2DUnCorrBarBtm = (TH2F *)gROOT->FindObject(histName2DUnCorrBarBtm);

      if (hist2DCorrBarBtm) {
        hist2DCorrBarsBtm.push_back(hist2DCorrBarBtm);
      }
      if (hist2DUnCorrBarBtm) {
        hist2DUnCorrBarsBtm.push_back(hist2DUnCorrBarBtm);
      }
    }

    // Merge all individual bar histograms into one for corrected times (Top)
    if (!hist2DCorrBarsTop.empty()) {
      TString histName2DCorrTop = Form("GoodTopTdcTimeCorr_2D_plane%d", plane);
      TH2F *hist2DCorrTop       = (TH2F *)hist2DCorrBarsTop[0]->Clone(histName2DCorrTop);
      for (size_t i = 1; i < hist2DCorrBarsTop.size(); ++i) {
        hist2DCorrTop->Add(hist2DCorrBarsTop[i]);
      }

      TCanvas *c2DCorrTop =
          new TCanvas(Form("c2DCorrTop_plane%d", plane), Form("Corrected Times 2D Top Plane %d", plane), 800, 600);
      hist2DCorrTop->SetTitle(Form("Corrected Times for Top Plane %d", plane));
      hist2DCorrTop->GetXaxis()->SetTitle("Bar Number");
      hist2DCorrTop->GetYaxis()->SetTitle("Time (ns)");
      hist2DCorrTop->Draw("COLZ");

      // Save merged histogram and canvas to the output file
      outputFile->cd("2D_Corrected");
      c2DCorrTop->Write();
    }

    // Merge all individual bar histograms into one for uncorrected times (Top)
    if (!hist2DUnCorrBarsTop.empty()) {
      TString histName2DUnCorrTop = Form("GoodTopTdcTimeUnCorr_2D_plane%d", plane);
      TH2F *hist2DUnCorrTop       = (TH2F *)hist2DUnCorrBarsTop[0]->Clone(histName2DUnCorrTop);
      for (size_t i = 1; i < hist2DUnCorrBarsTop.size(); ++i) {
        hist2DUnCorrTop->Add(hist2DUnCorrBarsTop[i]);
      }

      TCanvas *c2DUnCorrTop =
          new TCanvas(Form("c2DUnCorrTop_plane%d", plane), Form("Uncorrected Times 2D Top Plane %d", plane), 800, 600);
      hist2DUnCorrTop->SetTitle(Form("Uncorrected Times for Top Plane %d", plane));
      hist2DUnCorrTop->GetXaxis()->SetTitle("Bar Number");
      hist2DUnCorrTop->GetYaxis()->SetTitle("Time (ns)");
      hist2DUnCorrTop->Draw("COLZ");

      // Save merged histogram and canvas to the output file
      outputFile->cd("2D_Uncorrected");
      c2DUnCorrTop->Write();
    }

    // Merge all individual bar histograms into one for corrected times (Btm)
    if (!hist2DCorrBarsBtm.empty()) {
      TString histName2DCorrBtm = Form("GoodBtmTdcTimeCorr_2D_plane%d", plane);
      TH2F *hist2DCorrBtm       = (TH2F *)hist2DCorrBarsBtm[0]->Clone(histName2DCorrBtm);
      for (size_t i = 1; i < hist2DCorrBarsBtm.size(); ++i) {
        hist2DCorrBtm->Add(hist2DCorrBarsBtm[i]);
      }

      TCanvas *c2DCorrBtm =
          new TCanvas(Form("c2DCorrBtm_plane%d", plane), Form("Corrected Times 2D Btm Plane %d", plane), 800, 600);
      hist2DCorrBtm->SetTitle(Form("Corrected Times for Btm Plane %d", plane));
      hist2DCorrBtm->GetXaxis()->SetTitle("Bar Number");
      hist2DCorrBtm->GetYaxis()->SetTitle("Time (ns)");
      hist2DCorrBtm->Draw("COLZ");

      // Save merged histogram and canvas to the output file
      outputFile->cd("2D_Corrected");
      c2DCorrBtm->Write();
    }

    // Merge all individual bar histograms into one for uncorrected times (Btm)
    if (!hist2DUnCorrBarsBtm.empty()) {
      TString histName2DUnCorrBtm = Form("GoodBtmTdcTimeUnCorr_2D_plane%d", plane);
      TH2F *hist2DUnCorrBtm       = (TH2F *)hist2DUnCorrBarsBtm[0]->Clone(histName2DUnCorrBtm);
      for (size_t i = 1; i < hist2DUnCorrBarsBtm.size(); ++i) {
        hist2DUnCorrBtm->Add(hist2DUnCorrBarsBtm[i]);
      }

      TCanvas *c2DUnCorrBtm =
          new TCanvas(Form("c2DUnCorrBtm_plane%d", plane), Form("Uncorrected Times 2D Btm Plane %d", plane), 800, 600);
      hist2DUnCorrBtm->SetTitle(Form("Uncorrected Times for Btm Plane %d", plane));
      hist2DUnCorrBtm->GetXaxis()->SetTitle("Bar Number");
      hist2DUnCorrBtm->GetYaxis()->SetTitle("Time (ns)");
      hist2DUnCorrBtm->Draw("COLZ");

      // Save merged histogram and canvas to the output file
      outputFile->cd("2D_Uncorrected");
      c2DUnCorrBtm->Write();
    }
  }
  ////////////////////////////////////////////////////////////////////////////////////////////

  // Close the output file
  outputFile->Close();

  // Close the file
  file->Close();
}
