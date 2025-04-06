// Created by ehingerl on 04/01/25.
// Produces low level histograms for LAD Hodo and GEMs. Had to use this, as DEF-Files failed to make certain histograms.
// 1D & 2D histograms can be defined using the scructs. The script will loop over all leaves, and evaluate cuts on a
// per-leaf basis
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2.h>
#include <TLeaf.h>
#include <TObjString.h>
#include <TString.h>
#include <TTree.h>
#include <functional> // Ensure std::function is included
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

const int NDDATA_MAX = 10000; // Maximum number of data points to read
struct Histo1DCommand {
  TString name;       // Histogram name
  TString title;      // Histogram title
  TString expression; // Expression to fill the histogram
  int bins;           // Number of bins
  double xMin;        // Minimum x-axis value
  double xMax;        // Maximum x-axis value
  TString cut;        // Cut expression
};
struct Histo2DCommand {
  TString name;         // Histogram name
  TString title;        // Histogram title
  TString expression_x; // Expression for x-axis
  TString expression_y; // Expression for y-axis
  int binsX;            // Number of bins in x-axis
  double xMin;          // Minimum x-axis value
  double xMax;          // Maximum x-axis value
  int binsY;            // Number of bins in y-axis
  double yMin;          // Minimum y-axis value
  double yMax;          // Maximum y-axis value
  TString cut;          // Cut expression
};

// Define a structure to store cut information
struct CutInfo {
  std::string var1; // First variable or branch name
  std::string op;   // Operator (e.g., ==, !=, >, <, >=, <=)
  std::string var2; // Second variable or branch name or constant value
};

// Function to parse a cut expression into CutInfo objects
std::vector<CutInfo> parseCutExpression(const TString &cutExpression) {
  std::vector<CutInfo> cuts;
  std::string cutExpr = cutExpression.Data();
  std::vector<std::string> tokens;
  tokens.reserve(10); // Reserve space for tokens to avoid multiple allocations
  size_t start = 0, end = 0;

  // Tokenize the cut expression based on "&&" or "||"
  while ((end = cutExpr.find_first_of("&|", start)) != std::string::npos) {
    if (end > start) {
      tokens.push_back(cutExpr.substr(start, end - start));
    }
    start = cutExpr.find_first_not_of("&|", end);
  }
  if (start < cutExpr.size()) {
    tokens.push_back(cutExpr.substr(start));
  }

  for (const auto &token : tokens) {
    std::string trimmedToken = token;
    trimmedToken.erase(remove(trimmedToken.begin(), trimmedToken.end(), ' '), trimmedToken.end()); // Remove spaces

    // Parse the token into var1, operator, and var2
    std::string var1, op, var2;
    size_t opPos = trimmedToken.find("==");
    if (opPos == std::string::npos)
      opPos = trimmedToken.find("!=");
    if (opPos == std::string::npos)
      opPos = trimmedToken.find(">=");
    if (opPos == std::string::npos)
      opPos = trimmedToken.find("<=");
    if (opPos == std::string::npos)
      opPos = trimmedToken.find(">");
    if (opPos == std::string::npos)
      opPos = trimmedToken.find("<");

    if (opPos != std::string::npos) {
      var1 = trimmedToken.substr(0, opPos);
      if (trimmedToken.substr(opPos, 2) == "==" || trimmedToken.substr(opPos, 2) == "!=" ||
          trimmedToken.substr(opPos, 2) == ">=" || trimmedToken.substr(opPos, 2) == "<=") {
        op = trimmedToken.substr(opPos, 2);
      } else {
        op = trimmedToken.substr(opPos, 1); // Handle single-character operators like > or <
      }
      var2 = trimmedToken.substr(opPos + op.length());
      CutInfo cut;
      cut.var1 = var1;
      cut.op   = op;
      cut.var2 = var2;
      cuts.push_back(cut);
    } else {
      std::cerr << "Error: Invalid cut token format: " << trimmedToken << std::endl;
    }
  }
  return cuts;
}

bool getValues(const TString &var, std::vector<double> &values,
               const std::map<std::string, std::vector<Double_t>> &branchDataDouble,
               const std::map<std::string, std::vector<Int_t>> &branchDataInt, int nEntries_Evt) {
  if (nEntries_Evt < 0)
    return false;
  if (branchDataDouble.find(var.Data()) != branchDataDouble.end()) {
    values = branchDataDouble.at(var.Data());
  } else if (branchDataInt.find(var.Data()) != branchDataInt.end()) {
    for (const auto &val : branchDataInt.at(var.Data())) {
      values.push_back(static_cast<double>(val));
    }
  } else if (var.IsFloat() || var.IsDigit()) {
    values.insert(values.end(), nEntries_Evt, var.Atof());
  } else {
    std::cerr << "Error: Could not find or convert " << var << std::endl;
    return false;
  }
  return true;
}

void processCutExpression(const TString &cutExpression, const std::function<void(const TString &)> &processBranch) {
  if (!cutExpression.IsNull()) {
    auto cutInfos = parseCutExpression(cutExpression);
    for (const auto &cut : cutInfos) {
      processBranch(cut.var1.c_str());
      processBranch(cut.var2.c_str());
    }
  }
}

// void gem_histos(const char *inputFileName, const char *treeName) {
void lad_histos_expert(
    const char *inputFileName = "/volatile/hallc/c-lad/ehingerl/ROOTfiles/COSMICS/LAD_wGEM_cosmic_hall_296_-1.root",
    bool lad_only = true, int num_evts = -1) {
  const char *treeName = "T";
  // Open the ROOT file and get the tree
  TFile *inputFile = TFile::Open(inputFileName, "UPDATE");
  if (!inputFile || inputFile->IsZombie()) {
    std::cerr << "Error: Could not open file " << inputFileName << std::endl;
    return;
  }

  TTree *tree = (TTree *)inputFile->Get(treeName);
  if (!tree) {
    std::cerr << "Error: Could not find tree " << treeName << " in file " << inputFileName << std::endl;
    inputFile->Close();
    return;
  }

  // Define the list of histogram commands
  std::vector<TString> prefixes;
  if (lad_only) {
    prefixes.push_back("L");
  } else {
    prefixes.push_back("H");
    prefixes.push_back("P");
  }
  std::vector<Histo1DCommand> hist1DCommands;
  std::vector<Histo2DCommand> hist2DCommands;

  for (const auto &prefix : prefixes) {
    TString prefix_lower = prefix;
    prefix_lower.ToLower();
    // Add 1D histogram commands for each prefix
    for (const TString &plane : {"000", "001", "100", "101", "200"}) {
      for (int bar = 0; bar < 11; ++bar) {
      TString eventType = (prefix_lower == "p") ? "1" : "2";
      hist1DCommands.push_back(
        {prefix_lower + "ladhod_" + plane + "_fullhit_time_" + TString::Format("%d", bar),
         "SHMS-LAD " + plane + " Full Hit Time Paddle " + TString::Format("%d", bar + 1) + "; Time [ns]; Counts",
         prefix + ".ladhod." + plane + ".HodoHitTime", 200, -3100, -800,
         "g.evtyp==" + eventType + "&&" + prefix + ".ladhod." + plane + ".HodoHitPaddleNum==" + TString::Format("%d", bar + 1)});
      }
    }

    // Add 2D histogram commands for each prefix
    
    prefix_lower.ToLower();
    for (const TString &plane : {"000", "001", "100", "101", "200"}) {
      for (int bar = 0; bar < 11; ++bar) {
      TString eventType = (prefix_lower == "p") ? "1" : "2";
      hist2DCommands.insert(
        hist2DCommands.end(),
        {{prefix_lower + "ladhod2D_" + plane + "_fullhit_time_v_ADC_" + TString::Format("%d", bar),
          "SHMS-LAD " + plane + " Full Hit Time Paddle " + TString::Format("%d", bar + 1) + "; Time [ns]; ADC Counts",
          prefix + ".ladhod." + plane + ".HodoHitTime", prefix + ".ladhod." + plane + ".HodoHitEdep", 
          200, -3100, -800, 200, 0, 300, 
          "g.evtyp==" + eventType + "&&" + prefix + ".ladhod." + plane + ".HodoHitPaddleNum==" + TString::Format("%d", bar + 1)}});
      }
    }
  }
  int nHistograms = hist1DCommands.size() + hist2DCommands.size();
  int histCount   = 0; // Counter for histograms

  // Loop over the histogram commands and create histograms
  Long64_t nEntries = tree->GetEntries();
  // if (num_evts > 0 && num_evts < nEntries) {
  //   nEntries = num_evts;
  // }

  // Declare the map to store 1D histograms
  static std::map<TString, TH1D *> histograms1D;
  static std::map<TString, TH2D *> histograms2D;

  // Declare maps to store branch data
  std::map<std::string, std::vector<Double_t>> branchDataDouble;
  std::map<std::string, std::vector<Int_t>> branchDataInt;

  // Loop through histogram expressions and add each tree name to branchDataDouble or branchDataInt
  auto processBranch = [&](const TString &branchName) {
    if (branchDataDouble.find(branchName.Data()) == branchDataDouble.end() &&
        branchDataInt.find(branchName.Data()) == branchDataInt.end()) {
      TBranch *branch = tree->GetBranch(branchName.Data());
      if (branch) {
        TLeaf *leaf = branch->GetLeaf(branchName.Data());
        if (leaf) {
          TString leafType = leaf->GetTypeName();
          if (leafType == "Double_t" || leafType == "Float_t") {
            branchDataDouble[branchName.Data()] = std::vector<Double_t>(NDDATA_MAX);
          } else if (leafType == "Int_t") {
            branchDataInt[branchName.Data()] = std::vector<Int_t>(NDDATA_MAX);
          } else {
            std::cerr << "Error: Unsupported leaf type " << leafType << " for branch " << branchName << std::endl;
          }
        }
        // delete leaf;
      }
      // delete branch;
    }
  };

  // Process 1D histogram expressions and cuts
  for (const auto &cmd : hist1DCommands) {
    if (!cmd.expression.IsNull()) {
      processBranch(cmd.expression);
    }
    processCutExpression(cmd.cut, processBranch);
  }

  // Process 2D histogram expressions and cuts
  for (const auto &cmd : hist2DCommands) {
    if (!cmd.expression_x.IsNull()) {
      processBranch(cmd.expression_x);
    }
    if (!cmd.expression_y.IsNull()) {
      processBranch(cmd.expression_y);
    }
    processCutExpression(cmd.cut, processBranch);
  }

  // Set branch addresses for all branches
  for (auto &branchPair : branchDataDouble) {
    tree->SetBranchAddress(branchPair.first.c_str(), branchPair.second.data());
  }
  for (auto &branchPair : branchDataInt) {
    tree->SetBranchAddress(branchPair.first.c_str(), branchPair.second.data());
  }
  std::map<std::string, int> nEntriesMap;

  for (auto &branchPair : branchDataDouble) {
    TString nEntriesBranchName = "Ndata." + TString(branchPair.first.c_str());
    TBranch *nEntriesBranch    = tree->GetBranch(nEntriesBranchName.Data());
    if (nEntriesBranch) {
      nEntriesMap[branchPair.first] = 0; // Initialize nEntries for this branch
      tree->SetBranchAddress(nEntriesBranchName.Data(), &nEntriesMap[branchPair.first]);
    } else {
      nEntriesMap[branchPair.first] = 1;
    }
  }

  for (auto &branchPair : branchDataInt) {
    TString nEntriesBranchName = "Ndata." + TString(branchPair.first.c_str());
    TBranch *nEntriesBranch    = tree->GetBranch(nEntriesBranchName.Data());
    if (nEntriesBranch) {
      nEntriesMap[branchPair.first] = 0; // Initialize nEntries for this branch
      tree->SetBranchAddress(nEntriesBranchName.Data(), &nEntriesMap[branchPair.first]);
    }
  }

  for (Long64_t i = 0; i < nEntries; ++i) {
    // Clear all vectors
    for (auto &branchPair : branchDataDouble) {
      branchPair.second.clear();
      branchPair.second.resize(NDDATA_MAX, 0);
    }
    for (auto &branchPair : branchDataInt) {
      branchPair.second.clear();
      branchPair.second.resize(NDDATA_MAX, 0);
    }
    // Get the entry from the tree
    tree->GetEntry(i);
    // if (nEntries_Evt < 1)
    //   continue;

    // Loop over the 1D histogram commands
    for (const auto &cmd : hist1DCommands) {
      // Create the histogram if it doesn't exist
      if (histograms1D.find(cmd.name) == histograms1D.end()) {
        histograms1D[cmd.name] = new TH1D(cmd.name, cmd.title, cmd.bins, cmd.xMin, cmd.xMax);
      }
      TH1D *hist = histograms1D[cmd.name];

      // Parse the cut expression into CutInfo objects
      std::vector<CutInfo> cutInfos;
      if (!cmd.cut.IsNull()) {
        cutInfos = parseCutExpression(cmd.cut);
      }

      // Check if the cut expression is valid
      std::vector<bool> passCut(nEntriesMap[cmd.expression.Data()], true);
      for (const auto &cut : cutInfos) {
        std::vector<double> values1, values2;

        // Get the values of var1
        if (!getValues(cut.var1, values1, branchDataDouble, branchDataInt, nEntriesMap[cmd.expression.Data()])) {
          passCut = {false};
          break;
        }

        // Get the values of var2
        if (!getValues(cut.var2, values2, branchDataDouble, branchDataInt, nEntriesMap[cmd.expression.Data()])) {
          passCut = {false};
          break;
        }

        // Compare values based on the operator
        for (size_t j = 0; j < nEntriesMap[cmd.expression.Data()]; ++j) {
          if (cut.op == "==") {
            passCut[j] = (values1[j] == values2[j]);
          } else if (cut.op == "!=") {
            passCut[j] = (values1[j] != values2[j]);
          } else if (cut.op == ">") {
            passCut[j] = (values1[j] > values2[j]);
          } else if (cut.op == "<") {
            passCut[j] = (values1[j] < values2[j]);
          } else if (cut.op == ">=") {
            passCut[j] = (values1[j] >= values2[j]);
          } else if (cut.op == "<=") {
            passCut[j] = (values1[j] <= values2[j]);
          } else {
            std::cerr << "Error: Unknown operator " << cut.op << std::endl;
            passCut[j] = false;
            break;
          }
        }
      }

      std::vector<double> valueVec;
      if (!getValues(cmd.expression, valueVec, branchDataDouble, branchDataInt, nEntriesMap[cmd.expression.Data()])) {
        std::cerr << "Error: Could not retrieve value for expression " << cmd.expression << std::endl;
        continue;
      }
      for (size_t j = 0; j < nEntriesMap[cmd.expression.Data()]; ++j) {
        if (passCut[j]) {
          hist->Fill(valueVec[j]);
        }
      }
    }

    // Loop over the 2D histogram commands
    for (const auto &cmd : hist2DCommands) {
      // Create the 2D histogram if it doesn't exist
      if (histograms2D.find(cmd.name) == histograms2D.end()) {
        histograms2D[cmd.name] =
            new TH2D(cmd.name, cmd.title, cmd.binsX, cmd.xMin, cmd.xMax, cmd.binsY, cmd.yMin, cmd.yMax);
      }
      TH2D *hist2D = histograms2D[cmd.name];

      // Parse the cut expression into CutInfo objects
      std::vector<CutInfo> cutInfos;
      if (!cmd.cut.IsNull()) {
        cutInfos = parseCutExpression(cmd.cut);
      }

      // Check if the cut expression is valid
      std::vector<bool> passCut(nEntriesMap[cmd.expression_x.Data()], true);
      for (const auto &cut : cutInfos) {
        std::vector<double> values1, values2;

        // Get the values of var1
        if (!getValues(cut.var1, values1, branchDataDouble, branchDataInt, nEntriesMap[cmd.expression_x.Data()])) {
          passCut = {false};
          break;
        }

        // Get the values of var2
        if (!getValues(cut.var2, values2, branchDataDouble, branchDataInt, nEntriesMap[cmd.expression_x.Data()])) {
          passCut = {false};
          break;
        }

        // Compare values based on the operator
        for (size_t j = 0; j < nEntriesMap[cmd.expression_x.Data()]; ++j) {
          if (cut.op == "==") {
            passCut[j] = (values1[j] == values2[j]);
          } else if (cut.op == "!=") {
            passCut[j] = (values1[j] != values2[j]);
          } else if (cut.op == ">") {
            passCut[j] = (values1[j] > values2[j]);
          } else if (cut.op == "<") {
            passCut[j] = (values1[j] < values2[j]);
          } else if (cut.op == ">=") {
            passCut[j] = (values1[j] >= values2[j]);
          } else if (cut.op == "<=") {
            passCut[j] = (values1[j] <= values2[j]);
          } else {
            std::cerr << "Error: Unknown operator " << cut.op << std::endl;
            passCut[j] = false;
            break;
          }
        }
      }

      std::vector<double> valueVecX, valueVecY;
      if (!getValues(cmd.expression_x, valueVecX, branchDataDouble, branchDataInt,
                     nEntriesMap[cmd.expression_x.Data()])) {
        std::cerr << "Error: Could not retrieve value for expression " << cmd.expression_x << std::endl;
        continue;
      }
      if (!getValues(cmd.expression_y, valueVecY, branchDataDouble, branchDataInt,
                     nEntriesMap[cmd.expression_y.Data()])) {
        std::cerr << "Error: Could not retrieve value for expression " << cmd.expression_y << std::endl;
        continue;
      }
      for (size_t j = 0; j < nEntriesMap[cmd.expression_x.Data()]; ++j) {
        if (passCut[j]) {
          hist2D->Fill(valueVecX[j], valueVecY[j]);
        }
      }
    }

    // Print progress
    if (i % 10 == 0 || i == nEntries - 1) {
      double progress = ((double)i) / nEntries * 100.0;
      std::cout << "\rProgress: " << std::fixed << std::setprecision(1) << progress << "%" << std::flush;
    }
  }

  // Write and delete histograms
  for (auto &histPair : histograms1D) {
    histPair.second->Write();
    delete histPair.second;
  }
  for (auto &histPair : histograms2D) {
    histPair.second->Write();
    delete histPair.second;
  }
  // Close the file
  std::cout << "\nWriting histograms to file..." << std::endl;
  inputFile->Write();
  inputFile->Close();
  std::cout << "All histograms written to file." << std::endl;
}