#include "TFile.h"
#include "TH1D.h"
#include <iostream>

void ped_tracking_ladonly(TString golden_file = "", TString detector = "",
                           TString spect = "", Double_t polarity = 0) {

  if (golden_file == "") {
    cout << "Enter golden run root file name: " << endl;
    cin >> golden_file;
  }
  if (detector == "") {
    cout << "Enter detector prefix (hodo_1x (etc.), hgcer, aero, cal_prshwr, "
            "cal_shwr, ngcer, cal_hA (etc.)): "
         << endl;
    cin >> detector;
  }
  if (spect == "") {
    cout << "Enter a spectrometer (p or h): " << endl;
    cin >> spect;
  }
  if (polarity == 0) {
    cout << "Enter a detector polarity (top = 1, btm = 2): " << endl;
    cin >> polarity;
  }

  TString histname = Form("%s%s", spect.Data(), detector.Data());
  if (histname.Contains("ladhod") && histname.Contains("000") && polarity == 1)
    histname = Form("%s%s", histname.Data(), "_good_pped_vs_pmt_top");
  if (histname.Contains("ladhod") && histname.Contains("000") && polarity == 2)
    histname = Form("%s%s", histname.Data(), "_good_pped_vs_pmt_btm");
  if (histname.Contains("ladhod") && histname.Contains("001") && polarity == 1)
    histname = Form("%s%s", histname.Data(), "_good_pped_vs_pmt_top");
  if (histname.Contains("ladhod") && histname.Contains("001") && polarity == 2)
    histname = Form("%s%s", histname.Data(), "_good_pped_vs_pmt_btm");
  if (histname.Contains("ladhod") && histname.Contains("100") && polarity == 1)
    histname = Form("%s%s", histname.Data(), "_good_pped_vs_pmt_top");
  if (histname.Contains("ladhod") && histname.Contains("100") && polarity == 2)
    histname = Form("%s%s", histname.Data(), "_good_pped_vs_pmt_btm");
  if (histname.Contains("ladhod") && histname.Contains("101") && polarity == 1)
    histname = Form("%s%s", histname.Data(), "_good_pped_vs_pmt_top");
  if (histname.Contains("ladhod") && histname.Contains("101") && polarity == 2)
    histname = Form("%s%s", histname.Data(), "_good_pped_vs_pmt_btm");
  if (histname.Contains("ladhod") && histname.Contains("200") && polarity == 1)
    histname = Form("%s%s", histname.Data(), "_good_pped_vs_pmt_top");
  if (histname.Contains("ladhod") && histname.Contains("200") && polarity == 2)
    histname = Form("%s%s", histname.Data(), "_good_pped_vs_pmt_btm");
  if (histname.Contains("aero") && polarity == 1)
    histname = Form("%s%s", histname.Data(), "_good_pped_vs_pmt_top");
  if (histname.Contains("aero") && polarity == 2)
    histname = Form("%s%s", histname.Data(), "_good_pped_vs_pmt_btm");
  if (histname.Contains("hcer"))
    histname = Form("%s%s", histname.Data(), "_good_pped_vs_pmt");
  if (histname.Contains("hgcer"))
    histname = Form("%s%s", histname.Data(), "_good_pped_vs_pmt");
  if (histname.Contains("ngcer"))
    histname = Form("%s%s", histname.Data(), "_good_pped_vs_pmt");
  if (histname.Contains("_prshwr") && polarity == 1)
    histname = Form("%s%s", histname.Data(), "_good_pped_vs_pmt_top");
  if (histname.Contains("_prshwr") && polarity == 2)
    histname = Form("%s%s", histname.Data(), "_good_pped_vs_pmt_btm");
  if (histname.Contains("_shwr"))
    histname = Form("%s%s", histname.Data(), "_good_pped_vs_pmt");
  if (histname.Contains("hcal_h") && polarity == 1)
    histname = Form("%s%s", histname.Data(), "_good_pped_vs_pmt_top");
  if (histname.Contains("hcal_h") && polarity == 2)
    histname = Form("%s%s", histname.Data(), "_good_pped_vs_pmt_btm");

  TH2F* H1_ped_vs_pmt;
  TH2F* H2_ped_vs_pmt;

  TString protorootpath = Form("%s", gDirectory->GetPath());

  H2_ped_vs_pmt = (TH2F*)gDirectory->Get(histname.Data());

  TFile* f1 = new TFile(golden_file, "READ");
  if (f1->IsZombie()) {
    cout << "Cannot find : " << golden_file << endl;
    return;
  }

  f1->GetObject(histname.Data(), H1_ped_vs_pmt);

  gDirectory->cd(Form("%s", protorootpath.Data()));

  TH1D* H1_pmt;
  TH1D* H2_pmt;

  H1_pmt = (TH1D*)H1_ped_vs_pmt->ProjectionX("H1_pmt", 1,
                                             H1_ped_vs_pmt->GetSize() - 2);
  H2_pmt = (TH1D*)H2_ped_vs_pmt->ProjectionX("H2_pmt", 1,
                                             H2_ped_vs_pmt->GetSize() - 2);

  TH1D* H1_ped[H1_pmt->GetSize() - 2];
  TH1D* H2_ped[H2_pmt->GetSize() - 2];

  for (Int_t ipmt = 0; ipmt < (H1_pmt->GetSize() - 2); ipmt++) {
    H1_ped[ipmt] = (TH1D*)H1_ped_vs_pmt->ProjectionY(
        Form("H1_ped_pmt%d", ipmt + 1), ipmt + 1, ipmt + 1);
    H2_ped[ipmt] = (TH1D*)H2_ped_vs_pmt->ProjectionY(
        Form("H2_ped_pmt%d", ipmt + 1), ipmt + 1, ipmt + 1);
  }

  Double_t H1_ped_peak[H1_pmt->GetSize() - 2];
  Double_t H2_ped_peak[H2_pmt->GetSize() - 2];
  for (Int_t ipmt = 0; ipmt < (H1_pmt->GetSize() - 2); ipmt++) {

    if (H1_ped[ipmt]->GetEntries() > 25) {
      TSpectrum* s = new TSpectrum(1);
      gSystem->RedirectOutput("/dev/null", "a");
      s->Search(H1_ped[ipmt], 1.0, "nobackground&&nodraw", 0.001);
      gSystem->RedirectOutput(0);
      TList* functions = H1_ped[ipmt]->GetListOfFunctions();
      TPolyMarker* pm = (TPolyMarker*)functions->FindObject("TPolyMarker");
      if (pm) {
        H1_ped_peak[ipmt] = *pm->GetX();
      } else {
        H1_ped_peak[ipmt] = 1e+38;
      }
    } else {
      H1_ped_peak[ipmt] = 1e+38;
    }

    if (H2_ped[ipmt]->GetEntries() > 25) {
      TSpectrum* s = new TSpectrum(1);
      gSystem->RedirectOutput("/dev/null", "a");
      s->Search(H2_ped[ipmt], 1.0, "nobackground&&nodraw", 0.001);
      gSystem->RedirectOutput(0);
      TList* functions = H2_ped[ipmt]->GetListOfFunctions();
      TPolyMarker* pm = (TPolyMarker*)functions->FindObject("TPolyMarker");
      if (pm) {
        H2_ped_peak[ipmt] = *pm->GetX();
      } else {
        H2_ped_peak[ipmt] = 1e+38;
      }
    } else {
      H2_ped_peak[ipmt] = 1e+38;
    }
  }

  TF1* Gaussian =
      new TF1("Gaussian", "[0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))");
  Gaussian->SetParLimits(0, 0, 1000);
  Gaussian->SetParLimits(2, 0, 2);

  for (Int_t ipmt = 0; ipmt < (H1_pmt->GetSize() - 2); ipmt++) {
    Gaussian->SetParameter(1, H1_ped_peak[ipmt]);
    gSystem->RedirectOutput("/dev/null", "a");
    H1_ped[ipmt]->Fit(Gaussian, "QMN");
    gSystem->RedirectOutput(0);
    H1_pmt->SetBinContent(ipmt + 1, Gaussian->GetParameter(1));
    H1_pmt->SetBinError(ipmt + 1, Gaussian->GetParameter(2));
    Gaussian->SetParameter(1, H2_ped_peak[ipmt]);
    gSystem->RedirectOutput("/dev/null", "a");
    H2_ped[ipmt]->Fit(Gaussian, "QMN");
    gSystem->RedirectOutput(0);
    H2_pmt->SetBinContent(ipmt + 1, Gaussian->GetParameter(1));
    H2_pmt->SetBinError(ipmt + 1, Gaussian->GetParameter(2));
  }

  gSystem->RedirectOutput("/dev/null", "a");
  TH1D* Ped_Difference =
      new TH1D("Ped_Difference",
               Form("%s %s;PMT Number;  (Golden - Present) Pedestal Mean (mV)",
                    detector.Data(), (polarity == 1) ? "+" : "-"),
               (H1_pmt->GetSize() - 2), 0.5, (H1_pmt->GetSize() - 2) + 0.5);
  TH1D* Ped_Difference_Zero =
      new TH1D("Ped_Difference_Zero",
               Form("%s %s;PMT Number;  (Golden - Present) Pedestal Mean (mV)",
                    detector.Data(), (polarity == 1) ? "+" : "-"),
               (H1_pmt->GetSize() - 2), 0.5, (H1_pmt->GetSize() - 2) + 0.5);
  gSystem->RedirectOutput(0);

  Double_t histmaxdiff = 0.;
  for (Int_t ipmt = 0; ipmt < (H1_pmt->GetSize() - 2); ipmt++) {
    Double_t peddiff =
        (H1_pmt->GetBinContent(ipmt + 1) - H2_pmt->GetBinContent(ipmt + 1));
    if (TMath::Abs(peddiff) == 1e+38) {
      if (H1_pmt->GetBinContent(ipmt + 1) == 1e+38) {
        peddiff = -H2_pmt->GetBinContent(ipmt + 1);
      }
      if (H2_pmt->GetBinContent(ipmt + 1) == 1e+38)
        peddiff = H1_pmt->GetBinContent(ipmt + 1);
      Ped_Difference_Zero->SetBinContent(ipmt + 1, peddiff);
      Ped_Difference_Zero->SetBinError(
          ipmt + 1, sqrt(pow(H1_pmt->GetBinError(ipmt + 1), 2) +
                         pow(H2_pmt->GetBinError(ipmt + 1), 2)));
    } else {
      Ped_Difference_Zero->SetBinContent(ipmt + 1, 666e6);
    }
    Ped_Difference->SetBinContent(ipmt + 1, peddiff);
    if (TMath::Abs(peddiff) > histmaxdiff)
      histmaxdiff = TMath::Abs(peddiff);
    Ped_Difference->SetBinError(ipmt + 1,
                                sqrt(pow(H1_pmt->GetBinError(ipmt + 1), 2) +
                                     pow(H2_pmt->GetBinError(ipmt + 1), 2)));
  }
  histmaxdiff = histmaxdiff * 1.2;
  if (histmaxdiff < 15)
    histmaxdiff = 15;
  // Start Drawing
  gStyle->SetOptStat(0);
  Ped_Difference->SetAxisRange(-histmaxdiff, histmaxdiff, "Y");
  Ped_Difference->SetMarkerStyle(8);
  Ped_Difference->SetMarkerSize(1);
  Ped_Difference->DrawClone("PE1");
  Ped_Difference_Zero->SetAxisRange(-histmaxdiff, histmaxdiff, "Y");
  Ped_Difference_Zero->SetMarkerStyle(8);
  Ped_Difference_Zero->SetMarkerColor(kGray);
  Ped_Difference_Zero->SetLineColor(kGray);
  Ped_Difference_Zero->SetMarkerSize(1);
  Ped_Difference_Zero->DrawClone("PE1same");
  gPad->Update();
  TLine* Lower_Limit =
      new TLine(gPad->GetUxmin(), -3.5, gPad->GetUxmax(), -3.5);
  Lower_Limit->SetLineColor(kRed);
  Lower_Limit->SetLineWidth(3);
  Lower_Limit->SetLineStyle(2);
  Lower_Limit->Draw();
  TLine* Upper_Limit =
      new TLine(gPad->GetUxmin(), +3.5, gPad->GetUxmax(), +3.5);
  Upper_Limit->SetLineColor(kRed);
  Upper_Limit->SetLineWidth(3);
  Upper_Limit->SetLineStyle(2);
  Upper_Limit->Draw();

  delete Ped_Difference;
  delete Ped_Difference_Zero;
}
