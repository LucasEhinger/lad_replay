#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2.h>
#include <TLeaf.h>
#include <TObjString.h>
#include <TString.h>
#include <TTree.h>
#include <iostream>
#include <vector>

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

// void gem_histos(const char *inputFileName, const char *treeName) {
void gem_histos() {
  const char *inputFileName = "/volatile/hallc/c-lad/ehingerl/ROOTfiles/COSMICS/LAD_wGEM_cosmic_hall_296_-1.root";
  const char *treeName      = "T";
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
  std::vector<Histo1DCommand> histCommands = {
      {"h1_gem_clustWidthU_0", "X cluster size Layer 0; cluster size", "L.gem.clust.nstrip", 10, 0.5, 10.5,
       "L.gem.clust.axis<1&&L.gem.clust.layer<1"},
      {"h1_gem_clustWidthV_0", "Y cluster size Layer 0; cluster size", "L.gem.clust.nstrip", 10, 0.5, 10.5,
       "L.gem.clust.axis>0&&L.gem.clust.layer<1"},
      {"h1_gem_clustWidthU_1", "X cluster size Layer 1; cluster size", "L.gem.clust.nstrip", 10, 0.5, 10.5,
       "L.gem.clust.axis<1&&L.gem.clust.layer>0"},
      {"h1_gem_clustWidthV_1", "Y cluster size Layer 1; cluster size", "L.gem.clust.nstrip", 10, 0.5, 10.5,
       "L.gem.clust.axis>0&&L.gem.clust.layer>0"},
      {"h1_gem_clustSampMaxU_0", "Peak time sample Layer 0; Peak time sample (X strip)", "L.gem.clust.maxsamp", 6, -0.5,
       5.5, "L.gem.clust.axis<1&&L.gem.clust.layer<1"},
      {"h1_gem_clustSampMaxV_0", "Peak time sample Layer 0; Peak time sample (Y strip)", "L.gem.clust.maxsamp", 6, -0.5,
       5.5, "L.gem.clust.axis>0&&L.gem.clust.layer<1"},
      {"h1_gem_clustSampMaxU_1", "Peak time sample Layer 1; Peak time sample (X strip)", "L.gem.clust.maxsamp", 6, -0.5,
       5.5, "L.gem.clust.axis<1&&L.gem.clust.layer>0"},
      {"h1_gem_clustSampMaxV_1", "Peak time sample Layer 1; Peak time sample (Y strip)", "L.gem.clust.maxsamp", 6, -0.5,
       5.5, "L.gem.clust.axis>0&&L.gem.clust.layer>0"},
      {"h1_gem_clustADCMaxU_0", "X cluster Max strip ADC Layer 0; MAX strip ADC", "L.gem.clust.maxadc", 1500, 0, 1500,
       "L.gem.clust.axis<1&&L.gem.clust.layer<1"},
      {"h1_gem_clustADCMaxV_0", "Y cluster Max strip ADC Layer 0; MAX strip ADC", "L.gem.clust.maxadc", 1500, 0, 1500,
       "L.gem.clust.axis>0&&L.gem.clust.layer<1"},
      {"h1_gem_clustADCMaxU_1", "X cluster Max strip ADC Layer 1; MAX strip ADC", "L.gem.clust.maxadc", 1500, 0, 1500,
       "L.gem.clust.axis<1&&L.gem.clust.layer>0"},
      {"h1_gem_clustADCMaxV_1", "Y cluster Max strip ADC Layer 1; MAX strip ADC", "L.gem.clust.maxadc", 1500, 0, 1500,
       "L.gem.clust.axis>0&&L.gem.clust.layer>0"},
      {"h1_gem_clustADCSumU_0", "X cluster ADC sum Layer 0; ADC sum", "L.gem.clust.adc", 1500, 0, 3000,
       "L.gem.clust.axis<1&&L.gem.clust.layer<1"},
      {"h1_gem_clustADCSumV_0", "Y cluster ADC sum Layer 0; ADC sum", "L.gem.clust.adc", 1500, 0, 3000,
       "L.gem.clust.axis>0&&L.gem.clust.layer<1"},
      {"h1_gem_clustADCSumU_1", "X cluster ADC sum Layer 1; ADC sum", "L.gem.clust.adc", 1500, 0, 3000,
       "L.gem.clust.axis<1&&L.gem.clust.layer>0"},
      {"h1_gem_clustADCSumV_1", "Y cluster ADC sum Layer 1; ADC sum", "L.gem.clust.adc", 1500, 0, 3000,
       "L.gem.clust.axis>0&&L.gem.clust.layer>0"},
      {"h1_gem_clustTimeMeanU_0", "Cluster time 0", "L.gem.clust.adc", 6, -0.5, 5.5,
       "L.gem.clust.axis==0&&L.gem.clust.layer==0"},
      {"h1_gem_clustTimeMeanV_0", "Cluster time 0", "L.gem.clust.adc", 6, -0.5, 5.5,
       "L.gem.clust.axis==1&&L.gem.clust.layer==0"},
      {"h1_gem_clustTimeMeanU_1", "Cluster time 1", "L.gem.clust.adc", 6, -0.5, 5.5,
       "L.gem.clust.axis==0&&L.gem.clust.layer==1"},
      {"h1_gem_clustTimeMeanV_1", "Cluster time 1", "L.gem.clust.adc", 6, -0.5, 5.5,
       "L.gem.clust.axis==1&&L.gem.clust.layer==1"}
      // Add more histogram commands here as needed
  };

  // All TH2D already included in DEF file
  std::vector<Histo2DCommand> hist2DCommands = {};
  //     {"h2_gem_stripU_adc_0", "Layer 0 Y Strip vs ADC Strip Sum; Strip; ADC Sum", "L.gem.m0.strip.istrip",
  //      "L.gem.m0.strip.ADCsum", 1536, -0.5, 1535.5, 100, 0, 8000, "L.gem.m0.strip.IsU"},
  //     {"h2_gem_stripV_adc_0", "Layer 0 X Strip vs ADC Strip Sum; Strip; ADC Sum", "L.gem.m0.strip.istrip",
  //      "L.gem.m0.strip.ADCsum", 3072, -0.5, 3071.5, 100, 0, 8000, "L.gem.m0.strip.IsV"},
  //     {"h2_gem_stripU_adc_1", "Layer 1 Y Strip vs ADC Strip Sum; Strip; ADC Sum", "L.gem.m1.strip.istrip",
  //      "L.gem.m1.strip.ADCsum", 1536, -0.5, 1535.5, 100, 0, 8000, "L.gem.m0.strip.IsU"},
  //     {"h2_gem_stripV_adc_1", "Layer 1 X Strip vs ADC Strip Sum; Strip; ADC Sum", "L.gem.m1.strip.istrip",
  //      "L.gem.m1.strip.ADCsum", 3072, -0.5, 3071.5, 100, 0, 8000, "L.gem.m0.strip.IsV"}
  // };

  int nHistograms = histCommands.size() + hist2DCommands.size();
  int histCount = 0; // Counter for histograms
  // Loop over the histogram commands and create histograms
  for (const auto &cmd : histCommands) {
    // Create the histogram
    TH1D *hist = new TH1D(cmd.name, cmd.title, cmd.bins, cmd.xMin, cmd.xMax);

    // Determine the type of "expression" in the ROOT tree
    // Get the branch and leaf for the histogram expression
    auto branch = tree->GetBranch(cmd.expression.Data());
    if (!branch) {
      std::cerr << "Error: Could not find branch " << cmd.expression << " in tree " << treeName << std::endl;
      continue;
    }
    TLeaf *leaf = branch->GetLeaf(cmd.expression.Data());
    if (!leaf) {
      std::cerr << "Error: Could not find leaf " << cmd.expression << " in branch " << branch->GetName() << std::endl;
      continue;
    }
    // Loop over the entries in the tree and fill the histogram
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
      tree->GetEntry(i);

      // Check if the cut expression is valid
      for (int j = 0; j < leaf->GetLen(); ++j) {
        bool passCut = true; // Initialize passCut to true
        // Parse the cut expression and get branches/leaves for any variables used in the cut

        if (!cmd.cut.IsNull()) {
          TObjArray *cutTokens = cmd.cut.Tokenize("&&||"); // Tokenize the cut expression
          for (int i = 0; i < cutTokens->GetEntries(); ++i) {
            TString token = ((TObjString *)cutTokens->At(i))->GetString();
            token.ReplaceAll(" ", ""); // Remove any spaces

            // Parse the token into var1, operator, and var2
            TString var1, op, var2;
            Ssiz_t opPos = token.Index("==");
            if (opPos == kNPOS)
              opPos = token.Index("!=");
            if (opPos == kNPOS)
              opPos = token.Index(">=");
            if (opPos == kNPOS)
              opPos = token.Index("<=");
            if (opPos == kNPOS)
              opPos = token.Index(">");
            if (opPos == kNPOS)
              opPos = token.Index("<");

            if (opPos != kNPOS) {
              var1 = token(0, opPos);
              if (token(opPos, 2) == "==") { // Explicitly check for "=="
                op = "==";
              } else if (token(opPos, 2) == "!=") {
                op = "!=";
              } else if (token(opPos, 2) == ">=") {
                op = ">=";
              } else if (token(opPos, 2) == "<=") {
                op = "<=";
              } else {
                op = token(opPos, 1); // Handle single-character operators like > or <
              }
              var2 = token(opPos + op.Length(), token.Length() - opPos - op.Length());
            } else {
              std::cerr << "Error: Invalid cut token format: " << token << std::endl;
              passCut = false;
              break;
            }

            // Get the value of var1
            double value1   = 0.0;
            TLeaf *cutLeaf1 = tree->GetLeaf(var1.Data());
            if (cutLeaf1) {
              value1 = cutLeaf1->GetValue(j);
            } else {
              // If var1 is not a leaf, try converting it to a double
              if (var1.IsFloat() || var1.IsDigit()) {
                value1 = var1.Atof();
              } else {
                std::cerr << "Error: Could not find leaf or convert " << var1 << " in tree " << treeName << std::endl;
                passCut = false;
                break;
              }
            }

            // Get the value of var2
            double value2   = 0.0;
            TLeaf *cutLeaf2 = tree->GetLeaf(var2.Data());
            if (cutLeaf2) {
              value2 = cutLeaf2->GetValue(j);
            } else {
              // If var2 is not a leaf, try converting it to a double
              if (var2.IsFloat() || var2.IsDigit()) {
                value2 = var2.Atof();
              } else {
                std::cerr << "Error: Could not find leaf or convert " << var2 << " in tree " << treeName << std::endl;
                passCut = false;
                break;
              }
            }

            // Evaluate the condition
            if (op == "==") {
              passCut = (value1 == value2);
            } else if (op == "!=") {
              passCut = (value1 != value2);
            } else if (op == ">") {
              passCut = (value1 > value2);
            } else if (op == "<") {
              passCut = (value1 < value2);
            } else if (op == ">=") {
              passCut = (value1 >= value2);
            } else if (op == "<=") {
              passCut = (value1 <= value2);
            } else {
              std::cerr << "Error: Unknown operator " << op << std::endl;
              passCut = false;
              break;
            }

            // If any condition fails, break out of the loop
            if (!passCut) {
              break;
            }
          }
          delete cutTokens;
        }

        if (passCut) {
          // Array (double[] or int[])

          double value = leaf->GetValue(j);
          hist->Fill(value);
        }
      }
      // Print progress
      if (i % 10 == 0 || i == nEntries - 1) {
        double progress = (static_cast<double>(i + 1 + histCount * nEntries) / (nHistograms * nEntries)) * 100.0;
        std::cout << "\rProgress: " << progress << "%" << std::flush;
      }
    }
    histCount++;
    hist->Write();
    delete hist;
  }

  // Loop over the 2D histogram commands and create histograms
  for (const auto &cmd : hist2DCommands) {
    // Create the 2D histogram
    TH2D *hist2D = new TH2D(cmd.name, cmd.title, cmd.binsX, cmd.xMin, cmd.xMax, cmd.binsY, cmd.yMin, cmd.yMax);

    // Get the branches and leaves for the x and y expressions
    auto branchX = tree->GetBranch(cmd.expression_x.Data());
    auto branchY = tree->GetBranch(cmd.expression_y.Data());
    if (!branchX || !branchY) {
      std::cerr << "Error: Could not find branches for expressions " << cmd.expression_x << " or " << cmd.expression_y
                << " in tree " << treeName << std::endl;
      continue;
    }

    TLeaf *leafX = branchX->GetLeaf(cmd.expression_x.Data());
    TLeaf *leafY = branchY->GetLeaf(cmd.expression_y.Data());
    if (!leafX || !leafY) {
      std::cerr << "Error: Could not find leaves for expressions " << cmd.expression_x << " or " << cmd.expression_y
                << " in branches " << branchX->GetName() << " or " << branchY->GetName() << std::endl;
      continue;
    }

    // Loop over the entries in the tree and fill the 2D histogram
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
      tree->GetEntry(i);

      // Check if the cut expression is valid
      for (int j = 0; j < leafX->GetLen(); ++j) {
        bool passCut = true;

        if (!cmd.cut.IsNull()) {
          TObjArray *cutTokens = cmd.cut.Tokenize("&&||");
          for (int k = 0; k < cutTokens->GetEntries(); ++k) {
            TString token = ((TObjString *)cutTokens->At(k))->GetString();
            token.ReplaceAll(" ", "");

            TString var1, op, var2;
            Ssiz_t opPos = token.Index("==");
            if (opPos == kNPOS)
              opPos = token.Index("!=");
            if (opPos == kNPOS)
              opPos = token.Index(">=");
            if (opPos == kNPOS)
              opPos = token.Index("<=");
            if (opPos == kNPOS)
              opPos = token.Index(">");
            if (opPos == kNPOS)
              opPos = token.Index("<");

            if (opPos != kNPOS) {
              var1 = token(0, opPos);
              if (token(opPos, 2) == "==") { // Explicitly check for "=="
                op = "==";
              } else if (token(opPos, 2) == "!=") {
                op = "!=";
              } else if (token(opPos, 2) == ">=") {
                op = ">=";
              } else if (token(opPos, 2) == "<=") {
                op = "<=";
              } else {
                op = token(opPos, 1); // Handle single-character operators like > or <
              }
              var2 = token(opPos + op.Length(), token.Length() - opPos - op.Length());

            } else {
              std::cerr << "Error: Invalid cut token format: " << token << std::endl;
              passCut = false;
              break;
            }

            double value1   = 0.0;
            TLeaf *cutLeaf1 = tree->GetLeaf(var1.Data());
            if (cutLeaf1) {
              value1 = cutLeaf1->GetValue(j);
            } else {
              if (var1.IsFloat() || var1.IsDigit()) {
                value1 = var1.Atof();
              } else {
                std::cerr << "Error: Could not find leaf or convert " << var1 << " in tree " << treeName << std::endl;
                passCut = false;
                break;
              }
            }

            double value2   = 0.0;
            TLeaf *cutLeaf2 = tree->GetLeaf(var2.Data());
            if (cutLeaf2) {
              value2 = cutLeaf2->GetValue(j);
            } else {
              if (var2.IsFloat() || var2.IsDigit()) {
                value2 = var2.Atof();
              } else {
                std::cerr << "Error: Could not find leaf or convert " << var2 << " in tree " << treeName << std::endl;
                passCut = false;
                break;
              }
            }

            if (op == "==") {
              passCut = (value1 == value2);
            } else if (op == "!=") {
              passCut = (value1 != value2);
            } else if (op == ">") {
              passCut = (value1 > value2);
            } else if (op == "<") {
              passCut = (value1 < value2);
            } else if (op == ">=") {
              passCut = (value1 >= value2);
            } else if (op == "<=") {
              passCut = (value1 <= value2);
            } else {
              std::cerr << "Error: Unknown operator " << op << std::endl;
              passCut = false;
              break;
            }

            if (!passCut) {
              break;
            }
          }
          delete cutTokens;
        }

        if (passCut) {
          double valueX = leafX->GetValue(j);
          double valueY = leafY->GetValue(j);
          hist2D->Fill(valueX, valueY);
        }
      }

      // Print progress
      if (i % 10 == 0 || i == nEntries - 1) {
        double progress = (static_cast<double>(i + 1 + histCount * nEntries) / (nHistograms * nEntries)) * 100.0;
        std::cout << "\rProgress: " << progress << "%" << std::flush;
      }
    }
    histCount++;
    hist2D->Write();
    delete hist2D;
  }

  // Close the file
  inputFile->Close();
}