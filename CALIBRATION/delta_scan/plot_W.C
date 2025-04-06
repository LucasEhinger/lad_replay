#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLine.h>
#include <TTree.h>
#include <TROOT.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

// Function to parse the standard kinematics file and retrieve values
std::map<std::string, double> getKinematics(const char *filename, int runnumber, const char *prefix) {
  // Convert prefix to lowercase
  std::string lower_prefix = prefix;
  std::transform(lower_prefix.begin(), lower_prefix.end(), lower_prefix.begin(), ::tolower);
  prefix = lower_prefix.c_str();

  std::ifstream infile(filename);
  if (!infile.is_open()) {
    throw std::runtime_error("Error: Could not open kinematics file " + std::string(filename));
  }

  std::string line;
  std::map<std::string, double> result;
  bool found = false;

  while (std::getline(infile, line)) {
    // Ignore lines starting with '#'
    if (line.empty() || line[0] == '#') {
      continue;
    }

    // Check if the line specifies a run number range
    if (line.find('-') != std::string::npos && line.find('=') == std::string::npos) {
      std::istringstream iss(line);
      int startRun, endRun;
      char dash;
      iss >> startRun >> dash >> endRun;
      // Check if the run number is in the range
      if (runnumber >= startRun && runnumber <= endRun) {
        found = true;
        continue;
      }
    }

    // If the correct run range is found, parse the variables
    if (found && !line.empty()) {
      std::istringstream iss(line);
      std::string key;
      double value;
      char equals;

      if (iss >> key >> equals >> value && equals == '=') {
        result[key] = value;
      }
    }

    // Stop parsing if a new run range is encountered
    if (found && line.find('-') != std::string::npos && line.find('=') == std::string::npos) {
      break;
    }
  }

  infile.close();

  if (!found) {
    throw std::runtime_error("Error: Run number " + std::to_string(runnumber) + " not found in kinematics file.");
  }

  // Retrieve the required variables based on the prefix
  std::string theta_key    = std::string(prefix) + "theta_lab";
  std::string pcentral_key = std::string(prefix) + "pcentral";

  if (result.find(theta_key) == result.end() || result.find(pcentral_key) == result.end()) {
    throw std::runtime_error("Error: Required keys not found for prefix " + std::string(prefix));
  }

  return {{"theta_lab", result[theta_key]}, {"pcentral", result[pcentral_key]}};
}

void plot_W(const char *filename, const int runnum, const char *prefix = "P") {
  // Set ROOT to batch mode
  gROOT->SetBatch(kTRUE);
  // Adjustable constants
  const double W_min          = 0.0;
  const double W_max          = 2.0;
  const double W_nBins        = 100;
  const double dp_min         = -30.0;
  const double dp_max         = 30.0;
  const double dp_nBins       = 100;
  const double W_integral_min = 0.8; // Integration range for W
  const double W_integral_max = 1.05;

  // Convert prefix to uppercase
  std::string upper_prefix = prefix;
  std::transform(upper_prefix.begin(), upper_prefix.end(), upper_prefix.begin(), ::toupper);
  prefix = upper_prefix.c_str();
  // Open the ROOT file
  TFile *file = TFile::Open(filename, "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "Error: Could not open file " << filename << std::endl;
    return;
  }

  // Get the TTree
  TTree *tree = (TTree *)file->Get("T");
  if (!tree) {
    std::cerr << "Error: Could not find tree 'T' in file " << filename << std::endl;
    file->Close();
    return;
  }

  // Create histograms
  TH1D *hW    = new TH1D("hW", "W Distribution;W;Counts", W_nBins, W_min, W_max);
  TH1D *hdp   = new TH1D("hdp", "dp Distribution;dp (%);Counts", dp_nBins, dp_min, dp_max);
  TH2D *hW_dp = new TH2D("hW_dp", "W vs dp;W;dp (%)", W_nBins, W_min, W_max, dp_nBins, dp_min, dp_max);

  // Construct variable names with the prefix
  std::string W_var  = std::string(prefix) + ".kin.W";
  std::string dp_var = std::string(prefix) + ".gtr.dp";

  // Fill histograms
  tree->Draw((W_var + ">>hW").c_str(), "", "goff");
  tree->Draw((dp_var + ">>hdp").c_str(), "", "goff");
  tree->Draw((dp_var + ":" + W_var + ">>hW_dp").c_str(), "", "goff");

  // Integrate W over the specified range
  int bin_min     = hW->FindBin(W_integral_min);
  int bin_max     = hW->FindBin(W_integral_max);
  double integral = hW->Integral(bin_min, bin_max);
  std::cout << "Integral of W from " << W_integral_min << " to " << W_integral_max << " is " << integral << std::endl;

  // Draw histograms
  TCanvas *c1 = new TCanvas("c1", "Plots", 1200, 800);
  c1->Divide(2, 2);

  c1->cd(1);
  hW->Draw();
  // Highlight the integration region
  TLine *line_min = new TLine(W_integral_min, 0, W_integral_min, hW->GetMaximum());
  TLine *line_max = new TLine(W_integral_max, 0, W_integral_max, hW->GetMaximum());
  line_min->SetLineColor(kRed);
  line_max->SetLineColor(kRed);
  line_min->SetLineStyle(2);
  line_max->SetLineStyle(2);
  line_min->Draw("same");
  line_max->Draw("same");

  c1->cd(2);
  hdp->Draw();

  c1->cd(3);
  hW_dp->Draw("COLZ");

  c1->cd(4);
  // Retrieve kinematics information
  TLatex latex;
  latex.SetTextSize(0.05); // Increase font size
  std::string spectrometer = (std::string(prefix) == "P") ? "SHMS" : "HMS";
  latex.DrawLatex(0.1, 0.8, Form("Spectrometer: %s", spectrometer.c_str()));
  latex.DrawLatex(0.1, 0.7, Form("Run Number: %d", runnum));
  try {
    auto kinematics  = getKinematics("DBASE/LAD_COIN/standard.kinematics", runnum, prefix);
    double theta_lab = kinematics["theta_lab"];
    double pcentral  = kinematics["pcentral"];

    // Display the angle, momenta, and integral on the canvas
    latex.DrawLatex(0.1, 0.6, Form("%s Theta_lab: %.2f deg", prefix, theta_lab));
    latex.DrawLatex(0.1, 0.5, Form("%s Pcentral: %.2f GeV/c", prefix, pcentral));
    latex.DrawLatex(0.1, 0.4, Form("W Integral: %d events", static_cast<int>(integral)));
  } catch (const std::exception &e) {
    std::cerr << "Error retrieving kinematics: " << e.what() << std::endl;
  }

  // Construct the output filename
  auto kinematics = getKinematics("DBASE/LAD_COIN/standard.kinematics", runnum, prefix);
  double theta_lab = kinematics["theta_lab"];
  double pcentral = kinematics["pcentral"];
  std::string outfilename = "CALIBRATION/delta_scan/" + std::string(prefix) + "_delta_scan_p" + std::to_string(static_cast<int>(pcentral * 1000)) + "_theta" + std::to_string(static_cast<int>(theta_lab * 10)) + "_run" + std::to_string(runnum) + ".pdf";
  
  // Save the canvas
  c1->SaveAs(outfilename.c_str());

  // Clean up
  delete hW;
  delete hdp;
  delete hW_dp;
  delete line_min;
  delete line_max;
  delete c1;
  file->Close();
  delete file;
}