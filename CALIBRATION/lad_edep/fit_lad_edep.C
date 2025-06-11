#include <TCanvas.h>
#include <TF2.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2D.h>
#include <TProfile.h>
#include <string>

using namespace std;

const char spec_prefix = 'H'; // Spectrometer to replay

const int N_PLANES                 = 5;
const int N_PADDLES                = 11;
const string plane_names[N_PLANES] = {"000", "001", "100", "101", "200"};
const int N_SIDES                  = 2;
const string side_names[N_SIDES]   = {"Top", "Btm"};

const double fit_func_xmin = 2.2; // Minimum x value for fit function
const double fit_func_xmax = 6.0; // Maximum x value for fit function

double edep_fit_func_profile(double *x, double *par) {
  // x[0]: x, x[1]: y
  // par[0]: amplitude, par[1]: x-mean, par[2]: x-sigma, par[3]: y-mean, par[4]: y-sigma, par[5]: x0, par[6]: constant
  // 2D Gaussian for x < x0, else constant
  // Only one free parameter: par[0] is the y-axis scaling factor
  double y;
  if (x[0] < 3.5) {
    // Line between (1.5,40) and (3.5,180)
    double x1 = 1.5, y1 = 40.0;
    double x2 = 3.5, y2 = 180.0;
    double m = (y2 - y1) / (x2 - x1);
    double b = y1 - m * x1;
    y        = m * x[0] + b;
  } else {
    // Line between (3.5,180) and (6,20)
    double x1 = 3.5, y1 = 180.0;
    double x2 = 6.0, y2 = 20.0;
    double m = (y2 - y1) / (x2 - x1);
    double b = y1 - m * x1;
    y        = m * x[0] + b;
  }
  return par[0] * y;
}

double edep_fit_func_profile_quad(double *x, double *par) {
  double y;
  // Fit two quadratics: one for x < x_mid, one for x >= x_mid
  // par[0]: scaling
  // par[1]: a_left, par[2]: b_left, par[3]: c_left
  // par[4]: a_right, par[5]: b_right, par[6]: c_right
  double x_mid = 3.3;
  double x1, y1, x2, y2, a;
  if (x[0] < x_mid) {
    // Quadratic passing through (1,20) and (3.3,190)
    // One free parameter
    x1 = 1.0, y1 = 20.0;
    x2 = 3.3, y2 = 190.0;
    a = par[0]; // Scaling factor
  } else {
    // Quadratic passing through (3.3,190) and (7,20)
    x1 = 3.3, y1 = 190.0;
    x2 = 6.0, y2 = 40.0;
    a = par[1]; // Scaling factor
  }
  // double a = par[0];
  double b_tmp = (y2 - y1 - a * (x2 * x2 - x1 * x1)) / (x2 - x1);
  y            = a * (x[0] * x[0] - x1 * x1) + b_tmp * (x[0] - x1) + y1;

  return y;
}

void fit_lad_edep() {

  TFile *file = TFile::Open("lad_edep_plots_22609_H.root", "READ");
  if (!file || file->IsZombie()) {
    printf("Error: Cannot open ROOT file.\n");
    return;
  }

  TFile *outfile = TFile::Open("fit_lad_edep_output.root", "RECREATE");
  if (!outfile || outfile->IsZombie()) {
    printf("Error: Cannot create output ROOT file.\n");
    file->Close();
    return;
  }

  TCanvas *canvas = (TCanvas *)file->Get("KIN/TDC_DIFF_VS_ADC_AMP/c_TDC_DIFF_VS_ADC_AMP_plane_001");
  if (!canvas) {
    printf("Error: Cannot find canvas 'KIN/TDC_DIFF_VS_ADC_AMP'.\n");
    file->Close();
    return;
  }

  TList *primitives = canvas->GetListOfPrimitives();
  if (!primitives || primitives->GetSize() < 11) {
    printf("Error: Not enough histograms on the canvas.\n");
    file->Close();
    return;
  }

  TH2D *hist_arr[N_PLANES][N_PADDLES] = {nullptr};

  TH2D *hist = nullptr;
  int count  = 0;
  for (TObject *obj : *primitives) {
    // cout << "Object name: " << obj->GetName() << endl;
    if (obj->InheritsFrom("TPad")) {
      TPad *pad            = (TPad *)obj;
      TList *padPrimitives = pad->GetListOfPrimitives();
      for (TObject *padObj : *padPrimitives) {
        if (padObj->InheritsFrom("TH2")) {
          // cout << "Pad object name: " << padObj->GetName() << endl;
          // cout << "Pad object class: " << padObj->ClassName() << endl;
          hist_arr[0][count++] = (TH2D *)padObj;
          // Now you can use 'hist' as before
        }
      }
    }
  }

  TCanvas *c1 = new TCanvas("c1", "Multiple Pads", 800, 600);
  c1->Divide(4, 3);
  c1->Update();

  int pad_idx = 0;
  for (int i = 0; i < 1; ++i) {
    for (int j = 0; j < N_PADDLES; ++j) {
      TH2D *hist = hist_arr[i][j];
      if (!hist)
        continue;
      c1->cd(++pad_idx);
      TProfile *profile = new TProfile("profile", "Profile above line", hist->GetNbinsX(), hist->GetXaxis()->GetXmin(),
                                       hist->GetXaxis()->GetXmax());

      for (int binx = 1; binx <= hist->GetNbinsX(); ++binx) {
        double x      = hist->GetXaxis()->GetBinCenter(binx);
        double y_line = -64.5 * (x - 2) + 129;
        for (int biny = 1; biny <= hist->GetNbinsY(); ++biny) {
          double y = hist->GetYaxis()->GetBinCenter(biny);
          if (y > y_line) {
            double content = hist->GetBinContent(binx, biny);
            // Fill the profile with y, weighted by the bin content
            for (int k = 0; k < int(content); ++k) {
              profile->Fill(x, y);
            }
          }
        }
      }

      hist->Draw("colz");
      cout << "Drawing histogram: " << hist->GetName() << " Number " << pad_idx << endl;
      if (profile) {
        // profile->Draw("same");
        TF1 *edep_func = new TF1("edep_fit_func_profile_quad", edep_fit_func_profile_quad, fit_func_xmin, fit_func_xmax, 2);
        // TF1 *edep_func = new TF1("edep_fit_func_profile", edep_fit_func_profile, fit_func_xmin, fit_func_xmax, 1);
        // Initial guess for scaling
        profile->Fit(edep_func, "QNR");
        edep_func->SetLineColor(kRed);
        edep_func->SetLineWidth(2);
        edep_func->Draw("same");
        // profile->SetLineColor(kRed);
        // profile->Draw("same");
      }
      delete profile; // Clean up the profile after drawing
      if (pad_idx >= 12)
        break;
    }
    if (pad_idx >= 12)
      break;
  }

  // // Create the TF2 using the C++ method
  // TF2 *fitFunc = new TF2("fitFunc", edep_fit_func, fit_func_xmin, fit_func_xmax, hist->GetYaxis()->GetXmin(),
  //                        hist->GetYaxis()->GetXmax(), 7);

  // // Initial parameters: amplitude, x-mean, x-sigma, y-mean, y-sigma, x0, constant
  // fitFunc->SetParameters(hist->GetMaximum(),                    // [0]: amplitude
  //                        hist->GetMean(1),                      // [1]: x-mean
  //                        hist->GetRMS(1),                       // [2]: x-sigma
  //                        hist->GetMean(2),                      // [3]: y-mean
  //                        hist->GetRMS(2),                       // [4]: y-sigma
  //                        (fit_func_xmin + fit_func_xmax) / 2.0, // [5]: x0 (piecewise boundary)
  //                        0.0                                    // [6]: constant value for x >= x0
  // );
  // hist->Fit(fitFunc, "R");
  // fitFunc->Draw("same surf1");
  // // Now 'hist' points to the 11th histogram on the canvas

  outfile->cd();
  // c1->Update();
  c1->Write();
  outfile->Close();
  file->Close();
  delete c1;   // Clean up the canvas
  delete file; // Clean up the input file
  return;
}