// lad_tof_fast.C
//
// Multi-threaded LAD ToF analysis using RDataFrame.
//
// Reads a .dat file (one ROOT file path per line) into a TChain (tree "T"),
// books histograms per plane/paddle for tof, ypos, hittime, edep and
// edep_amp, and writes canvases (one subdirectory per variable) into the
// output ROOT file.
//
// Branch convention (hcana):
//   P.ladhod.goodhit_<var>_0    -> arrays for planes 0, 2, 4
//   P.ladhod.goodhit_<var>_1    -> arrays for planes 1, 3
//   Ndata.P.ladhod.goodhit_<var>_{0,1}  gives the array length per event.
//
// Per variable, the output ROOT file contains a subdirectory holding:
//   * 5 canvases (one per plane) each divided into 11 paddle pads.
//   * 1 summary canvas with 6 subpads: 5 plane sums + 1 grand total.
//
// Usage:
//   root -l -b -q 'lad_tof_fast.C("input.dat","lad_tof_fast.root")'
//

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <RVersion.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TROOT.h>

// Built-in IMT-safe RDataFrame progress bar exists from ROOT 6.30 onwards.
// The helper header has lived at two paths across versions; probe both.
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 30, 0)
#if __has_include(<ROOT/RDFHelpers.hxx>)
#include <ROOT/RDFHelpers.hxx>
#define LAD_HAS_RDF_PROGRESSBAR 1
#elif __has_include(<ROOT/RDF/RDFHelpers.hxx>)
#include <ROOT/RDF/RDFHelpers.hxx>
#define LAD_HAS_RDF_PROGRESSBAR 1
#endif
#endif

#include <atomic>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

// =====================================================================
// Histogram binning configuration -- edit these as needed
// =====================================================================
const int NBINS_TOF   = 200;
const double XMIN_TOF = 0.0;
const double XMAX_TOF = 100.0;

const int NBINS_YPOS   = 200;
const double XMIN_YPOS = -100.0;
const double XMAX_YPOS = 100.0;

const int NBINS_HITTIME   = 200;
const double XMIN_HITTIME = 1550.0;
const double XMAX_HITTIME = 2050.0;

const int NBINS_EDEP   = 200;
const double XMIN_EDEP = 0.0;
const double XMAX_EDEP = 20.0;

const int NBINS_EDEP_AMP   = 200;
const double XMIN_EDEP_AMP = 0.0;
const double XMAX_EDEP_AMP = 1.0;

const int NBINS_PUNCHTHROUGH   = 150;
const double XMIN_PUNCHTHROUGH = -5.0;
const double XMAX_PUNCHTHROUGH = 10.0;

const int N_PLANES  = 5;  // planes 0..4
const int N_PADDLES = 11; // paddles 0..10

// Default I/O paths (overridable via function arguments)
const char *DEFAULT_DAT_FILE = "files/all_C3_test.dat";
const char spec_prefix       = 'P'; // Spectrometer to replay
const char *DEFAULT_OUT_FILE = Form("files/root_fast/timing_C3_22745-23590_%c.root", spec_prefix);

// =====================================================================

struct VarConfig {
  std::string name;
  int nbins;
  double xmin;
  double xmax;
};

void lad_tof_fast(const char *dat_file = DEFAULT_DAT_FILE, const char *out_file = DEFAULT_OUT_FILE) {

  // Force batch mode: no graphical windows pop up while drawing canvases.
  gROOT->SetBatch(kTRUE);

  // Detach all TH1s from gDirectory so the output file only stores the
  // canvases we explicitly Write(); DrawCopy() embeds clones in the pads.
  TH1::AddDirectory(kFALSE);
  ROOT::EnableImplicitMT();

  // ----------------------------------------------------------------
  // 1. Build TChain from the .dat list
  // ----------------------------------------------------------------
  TChain chain("T");
  std::ifstream fin(dat_file);
  if (!fin.is_open()) {
    std::cerr << "[lad_tof_fast] ERROR: cannot open dat file '" << dat_file << "'\n";
    return;
  }
  int nfiles = 0;
  std::string line;
  while (std::getline(fin, line)) {
    // trim leading/trailing whitespace
    const size_t a = line.find_first_not_of(" \t\r\n");
    if (a == std::string::npos)
      continue;
    const size_t b   = line.find_last_not_of(" \t\r\n");
    std::string path = line.substr(a, b - a + 1);
    if (path.empty() || path[0] == '#')
      continue;
    chain.Add(path.c_str());
    ++nfiles;
  }
  std::cout << "[lad_tof_fast] Chained " << nfiles << " ROOT files. Total entries: " << chain.GetEntries() << "\n";
  if (chain.GetEntries() == 0) {
    std::cerr << "[lad_tof_fast] WARNING: chain is empty; aborting.\n";
    return;
  }

  // ----------------------------------------------------------------
  // 2. Build RDataFrame and alias the dotted branch names so we can
  //    refer to them in JIT-compiled Define() expressions.
  // ----------------------------------------------------------------
  ROOT::RDataFrame rdf(chain);
  ROOT::RDF::RNode df = rdf;

  // ---- Progress reporting --------------------------------------------
  // Modern ROOT: built-in, IMT-safe, prints a live percentage + ETA bar.
  // Older ROOT: hook OnPartialResult onto a Count() result. Count is
  // lazy, runs in the same event loop as the histograms, and its
  // partial result is the thread-safe running total of processed
  // entries -- exactly what we want for a progress fraction.
  ROOT::RDF::RResultPtr<ULong64_t> progress_count; // kept alive for the loop
#ifdef LAD_HAS_RDF_PROGRESSBAR
  ROOT::RDF::Experimental::AddProgressBar(df);
#else
  {
    const ULong64_t total        = chain.GetEntries();
    const ULong64_t report_every = std::max<ULong64_t>(1ULL, total / 200ULL); // ~0.5% steps
    auto t0        = std::make_shared<std::chrono::steady_clock::time_point>(std::chrono::steady_clock::now());
    progress_count = df.Count();
    // OnPartialResult is invoked thread-safely (never concurrently
    // with itself) every `report_every` processed entries.
    progress_count.OnPartialResult(report_every, [total, t0](ULong64_t partial) {
      const double frac = total ? double(partial) / double(total) : 0.0;
      const auto now    = std::chrono::steady_clock::now();
      const double secs = std::chrono::duration<double>(now - *t0).count();
      const double eta  = (frac > 0.0) ? secs * (1.0 / frac - 1.0) : 0.0;
      std::fprintf(stderr,
                   "\r[lad_tof_fast] %6.2f%%  (%llu / %llu)  "
                   "elapsed %5.1fs  eta %5.1fs   ",
                   100.0 * frac, (unsigned long long)partial, (unsigned long long)total, secs, eta);
      std::fflush(stderr);
    });
  }
#endif
  // --------------------------------------------------------------------

  for (const std::string &side : {std::string("0"), std::string("1")}) {
    df = df.Alias("plane_" + side, "P.ladhod.goodhit_plane_" + side);
    df = df.Alias("paddle_" + side, "P.ladhod.goodhit_paddle_" + side);
    df = df.Alias("hittime_" + side, "P.ladhod.goodhit_hittime_" + side);
    df = df.Alias("edep_" + side, "P.ladhod.goodhit_hitedep_" + side);
    df = df.Alias("edep_amp_" + side, "P.ladhod.goodhit_hitedep_amp_" + side);
    df = df.Alias("tof_" + side, "P.ladhod.goodhit_hit_tof_" + side);
    df = df.Alias("ypos_" + side, "P.ladhod.goodhit_hit_ypos_" + side);
  }

  // Histogram sets to produce
  const std::vector<VarConfig> vars = {
      {"tof", NBINS_TOF, XMIN_TOF, XMAX_TOF},
      {"ypos", NBINS_YPOS, XMIN_YPOS, XMAX_YPOS},
      {"hittime", NBINS_HITTIME, XMIN_HITTIME, XMAX_HITTIME},
      {"edep", NBINS_EDEP, XMIN_EDEP, XMAX_EDEP},
      {"edep_amp", NBINS_EDEP_AMP, XMIN_EDEP_AMP, XMAX_EDEP_AMP},
  };

  // ----------------------------------------------------------------
  // 3. Define masked array columns:
  //      <var>_p<plane>_b<paddle>  -- one paddle inside one plane
  //      <var>_p<plane>_sum        -- all paddles of one plane
  //      <var>_total               -- concat of _0 and _1 (all planes)
  //
  // Planes 0,2,4 live in the "_0" arrays, planes 1,3 in the "_1" arrays.
  // ----------------------------------------------------------------
  for (int plane = 0; plane < N_PLANES; ++plane) {
    const std::string side       = (plane % 2 == 0) ? "0" : "1";
    const std::string plane_col  = "plane_" + side;
    const std::string paddle_col = "paddle_" + side;

    // per-paddle masks
    for (int paddle = 0; paddle < N_PADDLES; ++paddle) {
      const std::string mask =
          "(" + plane_col + "==" + std::to_string(plane) + " && " + paddle_col + "==" + std::to_string(paddle) + ")";
      for (const auto &v : vars) {
        const std::string col = v.name + "_p" + std::to_string(plane) + "_b" + std::to_string(paddle);
        const std::string src = v.name + "_" + side;
        df                    = df.Define(col, src + "[" + mask + "]");
      }
    }

    // per-plane (all paddles) mask
    const std::string mask_plane = plane_col + "==" + std::to_string(plane);
    for (const auto &v : vars) {
      const std::string col = v.name + "_p" + std::to_string(plane) + "_sum";
      const std::string src = v.name + "_" + side;
      df                    = df.Define(col, src + "[" + mask_plane + "]");
    }
  }

  // grand total -- concat of _0 and _1 arrays for each variable
  for (const auto &v : vars) {
    const std::string col = v.name + "_total";
    df                    = df.Define(col, "ROOT::VecOps::Concatenate(" + v.name + "_0, " + v.name + "_1)");
  }

  // ----------------------------------------------------------------
  // 3b. Define punch-through time columns (hittime_1 - hittime_0)
  // ----------------------------------------------------------------
  // For each paddle, compute time difference between layer 1 and layer 0 hits
  for (int paddle = 0; paddle < N_PADDLES; ++paddle) {
    const std::string pt_col = "punchthrough_b" + std::to_string(paddle);
    const double pad_val     = paddle; // Capture specific value as double

    df = df.Define(pt_col,
                   [pad_val](ROOT::VecOps::RVec<double> pd0, ROOT::VecOps::RVec<double> ht0,
                             ROOT::VecOps::RVec<double> pd1, ROOT::VecOps::RVec<double> ht1) {
                     ROOT::VecOps::RVec<double> result;
                     // Get all hittimes for this paddle in layer 0
                     auto ht0_at_paddle = ht0[pd0 == pad_val];
                     // Get all hittimes for this paddle in layer 1
                     auto ht1_at_paddle = ht1[pd1 == pad_val];
                     // Calculate difference for each pair
                     for (size_t i = 0; i < ht0_at_paddle.size(); ++i) {
                       for (size_t j = 0; j < ht1_at_paddle.size(); ++j) {
                         result.push_back(ht1_at_paddle[j] - ht0_at_paddle[i]);
                       }
                     }
                     return result;
                   },
                   {"paddle_0", "hittime_0", "paddle_1", "hittime_1"});
  }

  // Define edep vs punchthrough columns (edep values paired with punchthrough time)
  for (int paddle = 0; paddle < N_PADDLES; ++paddle) {
    const double pad_val = paddle;

    // For side 0: edep values paired with punchthrough
    const std::string edep0_pt_col = "edep_0_vs_pt_b" + std::to_string(paddle);
    df                             = df.Define(edep0_pt_col,
                                               [pad_val](ROOT::VecOps::RVec<double> pd0, ROOT::VecOps::RVec<double> edep0,
                             ROOT::VecOps::RVec<double> pd1, ROOT::VecOps::RVec<double> ht0,
                             ROOT::VecOps::RVec<double> ht1) {
                     ROOT::VecOps::RVec<std::pair<double, double>> result;
                     auto edep0_at_paddle = edep0[pd0 == pad_val];
                     auto ht0_at_paddle   = ht0[pd0 == pad_val];
                     auto ht1_at_paddle   = ht1[pd1 == pad_val];
                     for (size_t i = 0; i < edep0_at_paddle.size(); ++i) {
                       for (size_t j = 0; j < ht1_at_paddle.size(); ++j) {
                         result.push_back({edep0_at_paddle[i], ht1_at_paddle[j] - ht0_at_paddle[i]});
                       }
                     }
                     return result;
                   },
                                               {"paddle_0", "edep_0", "paddle_1", "hittime_0", "hittime_1"});

    // For side 1: edep values paired with punchthrough
    const std::string edep1_pt_col = "edep_1_vs_pt_b" + std::to_string(paddle);
    df                             = df.Define(edep1_pt_col,
                                               [pad_val](ROOT::VecOps::RVec<double> pd1, ROOT::VecOps::RVec<double> edep1,
                             ROOT::VecOps::RVec<double> pd0, ROOT::VecOps::RVec<double> ht0,
                             ROOT::VecOps::RVec<double> ht1) {
                     ROOT::VecOps::RVec<std::pair<double, double>> result;
                     auto edep1_at_paddle = edep1[pd1 == pad_val];
                     auto ht0_at_paddle   = ht0[pd0 == pad_val];
                     auto ht1_at_paddle   = ht1[pd1 == pad_val];
                     for (size_t j = 0; j < edep1_at_paddle.size(); ++j) {
                       for (size_t i = 0; i < ht0_at_paddle.size(); ++i) {
                         result.push_back({edep1_at_paddle[j], ht1_at_paddle[j] - ht0_at_paddle[i]});
                       }
                     }
                     return result;
                   },
                                               {"paddle_1", "edep_1", "paddle_0", "hittime_0", "hittime_1"});

    // For edep_amp: side 0
    const std::string edepamp0_pt_col = "edepamp_0_vs_pt_b" + std::to_string(paddle);
    df                                = df.Define(edepamp0_pt_col,
                                                  [pad_val](ROOT::VecOps::RVec<double> pd0, ROOT::VecOps::RVec<double> edepamp0,
                             ROOT::VecOps::RVec<double> pd1, ROOT::VecOps::RVec<double> ht0,
                             ROOT::VecOps::RVec<double> ht1) {
                     ROOT::VecOps::RVec<std::pair<double, double>> result;
                     auto edepamp0_at_paddle = edepamp0[pd0 == pad_val];
                     auto ht0_at_paddle      = ht0[pd0 == pad_val];
                     auto ht1_at_paddle      = ht1[pd1 == pad_val];
                     for (size_t i = 0; i < edepamp0_at_paddle.size(); ++i) {
                       for (size_t j = 0; j < ht1_at_paddle.size(); ++j) {
                         result.push_back({edepamp0_at_paddle[i], ht1_at_paddle[j] - ht0_at_paddle[i]});
                       }
                     }
                     return result;
                   },
                                                  {"paddle_0", "edep_amp_0", "paddle_1", "hittime_0", "hittime_1"});

    // For edep_amp: side 1
    const std::string edepamp1_pt_col = "edepamp_1_vs_pt_b" + std::to_string(paddle);
    df                                = df.Define(edepamp1_pt_col,
                                                  [pad_val](ROOT::VecOps::RVec<double> pd1, ROOT::VecOps::RVec<double> edepamp1,
                             ROOT::VecOps::RVec<double> pd0, ROOT::VecOps::RVec<double> ht0,
                             ROOT::VecOps::RVec<double> ht1) {
                     ROOT::VecOps::RVec<std::pair<double, double>> result;
                     auto edepamp1_at_paddle = edepamp1[pd1 == pad_val];
                     auto ht0_at_paddle      = ht0[pd0 == pad_val];
                     auto ht1_at_paddle      = ht1[pd1 == pad_val];
                     for (size_t j = 0; j < edepamp1_at_paddle.size(); ++j) {
                       for (size_t i = 0; i < ht0_at_paddle.size(); ++i) {
                         result.push_back({edepamp1_at_paddle[j], ht1_at_paddle[j] - ht0_at_paddle[i]});
                       }
                     }
                     return result;
                   },
                                                  {"paddle_1", "edep_amp_1", "paddle_0", "hittime_0", "hittime_1"});
  }

  // Define tof vs edep columns (edep values paired with tof)
  for (int paddle = 0; paddle < N_PADDLES; ++paddle) {
    const double pad_val = paddle;

    // tof vs edep for side 0
    const std::string tof_edep0_col = "tof_edep_0_b" + std::to_string(paddle);
    df                              = df.Define(
        tof_edep0_col,
        [pad_val](ROOT::VecOps::RVec<double> pd0, ROOT::VecOps::RVec<double> tof0, ROOT::VecOps::RVec<double> edep0) {
          ROOT::VecOps::RVec<std::pair<double, double>> result;
          auto tof0_at_paddle  = tof0[pd0 == pad_val];
          auto edep0_at_paddle = edep0[pd0 == pad_val];
          for (size_t i = 0; i < tof0_at_paddle.size(); ++i) {
            if (i < edep0_at_paddle.size()) {
              result.push_back({tof0_at_paddle[i], edep0_at_paddle[i]});
            }
          }
          return result;
        },
        {"paddle_0", "tof_0", "edep_0"});

    // tof vs edep for side 1
    const std::string tof_edep1_col = "tof_edep_1_b" + std::to_string(paddle);
    df                              = df.Define(
        tof_edep1_col,
        [pad_val](ROOT::VecOps::RVec<double> pd1, ROOT::VecOps::RVec<double> tof1, ROOT::VecOps::RVec<double> edep1) {
          ROOT::VecOps::RVec<std::pair<double, double>> result;
          auto tof1_at_paddle  = tof1[pd1 == pad_val];
          auto edep1_at_paddle = edep1[pd1 == pad_val];
          for (size_t i = 0; i < tof1_at_paddle.size(); ++i) {
            if (i < edep1_at_paddle.size()) {
              result.push_back({tof1_at_paddle[i], edep1_at_paddle[i]});
            }
          }
          return result;
        },
        {"paddle_1", "tof_1", "edep_1"});

    // tof vs edep_amp for side 0
    const std::string tof_edepamp0_col = "tof_edepamp_0_b" + std::to_string(paddle);
    df                                 = df.Define(tof_edepamp0_col,
                                                   [pad_val](ROOT::VecOps::RVec<double> pd0, ROOT::VecOps::RVec<double> tof0,
                             ROOT::VecOps::RVec<double> edepamp0) {
                     ROOT::VecOps::RVec<std::pair<double, double>> result;
                     auto tof0_at_paddle     = tof0[pd0 == pad_val];
                     auto edepamp0_at_paddle = edepamp0[pd0 == pad_val];
                     for (size_t i = 0; i < tof0_at_paddle.size(); ++i) {
                       if (i < edepamp0_at_paddle.size()) {
                         result.push_back({tof0_at_paddle[i], edepamp0_at_paddle[i]});
                       }
                     }
                     return result;
                   },
                                                   {"paddle_0", "tof_0", "edep_amp_0"});

    // tof vs edep_amp for side 1
    const std::string tof_edepamp1_col = "tof_edepamp_1_b" + std::to_string(paddle);
    df                                 = df.Define(tof_edepamp1_col,
                                                   [pad_val](ROOT::VecOps::RVec<double> pd1, ROOT::VecOps::RVec<double> tof1,
                             ROOT::VecOps::RVec<double> edepamp1) {
                     ROOT::VecOps::RVec<std::pair<double, double>> result;
                     auto tof1_at_paddle     = tof1[pd1 == pad_val];
                     auto edepamp1_at_paddle = edepamp1[pd1 == pad_val];
                     for (size_t i = 0; i < tof1_at_paddle.size(); ++i) {
                       if (i < edepamp1_at_paddle.size()) {
                         result.push_back({tof1_at_paddle[i], edepamp1_at_paddle[i]});
                       }
                     }
                     return result;
                   },
                                                   {"paddle_1", "tof_1", "edep_amp_1"});
  }

  // Define x,y coordinate columns for 2D histograms from pair columns
  // edep vs punchthrough
  for (int paddle = 0; paddle < N_PADDLES; ++paddle) {
    const std::string pair_col0 = "edep_0_vs_pt_b" + std::to_string(paddle);
    const std::string x_col0    = "pt_0_b" + std::to_string(paddle);
    const std::string y_col0    = "edep_0_b_pt" + std::to_string(paddle);
    df                          = df.Define(x_col0,
                                            [](const ROOT::VecOps::RVec<std::pair<double, double>> &v) {
                     ROOT::VecOps::RVec<double> x;
                     for (const auto &p : v)
                       x.push_back(p.second);
                     return x;
                   },
                                            {pair_col0});
    df                          = df.Define(y_col0,
                                            [](const ROOT::VecOps::RVec<std::pair<double, double>> &v) {
                     ROOT::VecOps::RVec<double> y;
                     for (const auto &p : v)
                       y.push_back(p.first);
                     return y;
                   },
                                            {pair_col0});

    const std::string pair_col1 = "edep_1_vs_pt_b" + std::to_string(paddle);
    const std::string x_col1    = "pt_1_b" + std::to_string(paddle);
    const std::string y_col1    = "edep_1_b_pt" + std::to_string(paddle);
    df                          = df.Define(x_col1,
                                            [](const ROOT::VecOps::RVec<std::pair<double, double>> &v) {
                     ROOT::VecOps::RVec<double> x;
                     for (const auto &p : v)
                       x.push_back(p.second);
                     return x;
                   },
                                            {pair_col1});
    df                          = df.Define(y_col1,
                                            [](const ROOT::VecOps::RVec<std::pair<double, double>> &v) {
                     ROOT::VecOps::RVec<double> y;
                     for (const auto &p : v)
                       y.push_back(p.first);
                     return y;
                   },
                                            {pair_col1});

    const std::string pair_colamp0 = "edepamp_0_vs_pt_b" + std::to_string(paddle);
    const std::string x_colamp0    = "ptamp_0_b" + std::to_string(paddle);
    const std::string y_colamp0    = "edepamp_0_b_pt" + std::to_string(paddle);
    df                             = df.Define(x_colamp0,
                                               [](const ROOT::VecOps::RVec<std::pair<double, double>> &v) {
                     ROOT::VecOps::RVec<double> x;
                     for (const auto &p : v)
                       x.push_back(p.second);
                     return x;
                   },
                                               {pair_colamp0});
    df                             = df.Define(y_colamp0,
                                               [](const ROOT::VecOps::RVec<std::pair<double, double>> &v) {
                     ROOT::VecOps::RVec<double> y;
                     for (const auto &p : v)
                       y.push_back(p.first);
                     return y;
                   },
                                               {pair_colamp0});

    const std::string pair_colamp1 = "edepamp_1_vs_pt_b" + std::to_string(paddle);
    const std::string x_colamp1    = "ptamp_1_b" + std::to_string(paddle);
    const std::string y_colamp1    = "edepamp_1_b_pt" + std::to_string(paddle);
    df                             = df.Define(x_colamp1,
                                               [](const ROOT::VecOps::RVec<std::pair<double, double>> &v) {
                     ROOT::VecOps::RVec<double> x;
                     for (const auto &p : v)
                       x.push_back(p.second);
                     return x;
                   },
                                               {pair_colamp1});
    df                             = df.Define(y_colamp1,
                                               [](const ROOT::VecOps::RVec<std::pair<double, double>> &v) {
                     ROOT::VecOps::RVec<double> y;
                     for (const auto &p : v)
                       y.push_back(p.first);
                     return y;
                   },
                                               {pair_colamp1});
  }

  // tof vs edep
  for (int paddle = 0; paddle < N_PADDLES; ++paddle) {
    const std::string pair_col0 = "tof_edep_0_b" + std::to_string(paddle);
    const std::string x_col0    = "tof_0_b" + std::to_string(paddle);
    const std::string y_col0    = "edep_0_b_tof" + std::to_string(paddle);
    df                          = df.Define(x_col0,
                                            [](const ROOT::VecOps::RVec<std::pair<double, double>> &v) {
                     ROOT::VecOps::RVec<double> x;
                     for (const auto &p : v)
                       x.push_back(p.first);
                     return x;
                   },
                                            {pair_col0});
    df                          = df.Define(y_col0,
                                            [](const ROOT::VecOps::RVec<std::pair<double, double>> &v) {
                     ROOT::VecOps::RVec<double> y;
                     for (const auto &p : v)
                       y.push_back(p.second);
                     return y;
                   },
                                            {pair_col0});

    const std::string pair_col1 = "tof_edep_1_b" + std::to_string(paddle);
    const std::string x_col1    = "tof_1_b" + std::to_string(paddle);
    const std::string y_col1    = "edep_1_b_tof" + std::to_string(paddle);
    df                          = df.Define(x_col1,
                                            [](const ROOT::VecOps::RVec<std::pair<double, double>> &v) {
                     ROOT::VecOps::RVec<double> x;
                     for (const auto &p : v)
                       x.push_back(p.first);
                     return x;
                   },
                                            {pair_col1});
    df                          = df.Define(y_col1,
                                            [](const ROOT::VecOps::RVec<std::pair<double, double>> &v) {
                     ROOT::VecOps::RVec<double> y;
                     for (const auto &p : v)
                       y.push_back(p.second);
                     return y;
                   },
                                            {pair_col1});

    const std::string pair_colamp0 = "tof_edepamp_0_b" + std::to_string(paddle);
    const std::string x_colamp0    = "tofamp_0_b" + std::to_string(paddle);
    const std::string y_colamp0    = "edepamp_0_b_tof" + std::to_string(paddle);
    df                             = df.Define(x_colamp0,
                                               [](const ROOT::VecOps::RVec<std::pair<double, double>> &v) {
                     ROOT::VecOps::RVec<double> x;
                     for (const auto &p : v)
                       x.push_back(p.first);
                     return x;
                   },
                                               {pair_colamp0});
    df                             = df.Define(y_colamp0,
                                               [](const ROOT::VecOps::RVec<std::pair<double, double>> &v) {
                     ROOT::VecOps::RVec<double> y;
                     for (const auto &p : v)
                       y.push_back(p.second);
                     return y;
                   },
                                               {pair_colamp0});

    const std::string pair_colamp1 = "tof_edepamp_1_b" + std::to_string(paddle);
    const std::string x_colamp1    = "tofamp_1_b" + std::to_string(paddle);
    const std::string y_colamp1    = "edepamp_1_b_tof" + std::to_string(paddle);
    df                             = df.Define(x_colamp1,
                                               [](const ROOT::VecOps::RVec<std::pair<double, double>> &v) {
                     ROOT::VecOps::RVec<double> x;
                     for (const auto &p : v)
                       x.push_back(p.first);
                     return x;
                   },
                                               {pair_colamp1});
    df                             = df.Define(y_colamp1,
                                               [](const ROOT::VecOps::RVec<std::pair<double, double>> &v) {
                     ROOT::VecOps::RVec<double> y;
                     for (const auto &p : v)
                       y.push_back(p.second);
                     return y;
                   },
                                               {pair_colamp1});
  }

  // ----------------------------------------------------------------
  // 3c. Define additional columns for new plots
  // ----------------------------------------------------------------

  // Punch-through for specific plane pairs: 0-1 and 2-3
  for (int paddle = 0; paddle < N_PADDLES; ++paddle) {
    const double pad_val = paddle;
    
    // Punchthrough between planes 0 and 1
    const std::string pt_01_col = "punchthrough_01_b" + std::to_string(paddle);
    df = df.Define(pt_01_col, [pad_val](ROOT::VecOps::RVec<double> plane0, ROOT::VecOps::RVec<double> paddle0,
                                        ROOT::VecOps::RVec<double> ht0,
                                        ROOT::VecOps::RVec<double> plane1, ROOT::VecOps::RVec<double> paddle1,
                                        ROOT::VecOps::RVec<double> ht1) {
      ROOT::VecOps::RVec<double> result;
      // Get hits in plane 0 at this paddle
      for (size_t i = 0; i < plane0.size(); ++i) {
        if (plane0[i] == 0.0 && paddle0[i] == pad_val) {
          // Get hits in plane 1 at this paddle
          for (size_t j = 0; j < plane1.size(); ++j) {
            if (plane1[j] == 1.0 && paddle1[j] == pad_val) {
              result.push_back(ht1[j] - ht0[i]);
            }
          }
        }
      }
      return result;
    }, {"plane_0", "paddle_0", "hittime_0", "plane_1", "paddle_1", "hittime_1"});
    
    // Punchthrough between planes 2 and 3
    const std::string pt_23_col = "punchthrough_23_b" + std::to_string(paddle);
    df = df.Define(pt_23_col, [pad_val](ROOT::VecOps::RVec<double> plane0, ROOT::VecOps::RVec<double> paddle0,
                                        ROOT::VecOps::RVec<double> ht0,
                                        ROOT::VecOps::RVec<double> plane1, ROOT::VecOps::RVec<double> paddle1,
                                        ROOT::VecOps::RVec<double> ht1) {
      ROOT::VecOps::RVec<double> result;
      // Get hits in plane 2 at this paddle (plane 2 is in the _0 array)
      for (size_t i = 0; i < plane0.size(); ++i) {
        if (plane0[i] == 2.0 && paddle0[i] == pad_val) {
          // Get hits in plane 3 at this paddle (plane 3 is in the _1 array)
          for (size_t j = 0; j < plane1.size(); ++j) {
            if (plane1[j] == 3.0 && paddle1[j] == pad_val) {
              result.push_back(ht1[j] - ht0[i]);
            }
          }
        }
      }
      return result;
    }, {"plane_0", "paddle_0", "hittime_0", "plane_1", "paddle_1", "hittime_1"});
  }

  // Front veto: hittime and tof for back (_1) with no hits in front (_0)
  df = df.Define("hittime_1_front_veto", 
                 [](ROOT::VecOps::RVec<double> pd0, ROOT::VecOps::RVec<double> pd1, 
                    ROOT::VecOps::RVec<double> ht1) {
    ROOT::VecOps::RVec<double> result;
    // For each hit in back, check if there's a hit in front at same paddle
    for (size_t j = 0; j < pd1.size(); ++j) {
      bool has_front_hit = false;
      for (size_t i = 0; i < pd0.size(); ++i) {
        if (pd0[i] == pd1[j]) {
          has_front_hit = true;
          break;
        }
      }
      // Only include if no front hit at this paddle
      if (!has_front_hit) {
        result.push_back(ht1[j]);
      }
    }
    return result;
  }, {"paddle_0", "paddle_1", "hittime_1"});

  df = df.Define("tof_1_front_veto", 
                 [](ROOT::VecOps::RVec<double> pd0, ROOT::VecOps::RVec<double> pd1, 
                    ROOT::VecOps::RVec<double> tof1) {
    ROOT::VecOps::RVec<double> result;
    for (size_t j = 0; j < pd1.size(); ++j) {
      bool has_front_hit = false;
      for (size_t i = 0; i < pd0.size(); ++i) {
        if (pd0[i] == pd1[j]) {
          has_front_hit = true;
          break;
        }
      }
      if (!has_front_hit) {
        result.push_back(tof1[j]);
      }
    }
    return result;
  }, {"paddle_0", "paddle_1", "tof_1"});

  // Photon path length corrected columns (subtract 55.52563)
  const double path_corr = 55.52563;
  
  df = df.Define("tof_0_corrected", 
                 [path_corr](ROOT::VecOps::RVec<double> tof0) {
    ROOT::VecOps::RVec<double> result;
    for (const auto& v : tof0) result.push_back(v - path_corr);
    return result;
  }, {"tof_0"});

  df = df.Define("tof_1_corrected", 
                 [path_corr](ROOT::VecOps::RVec<double> tof1) {
    ROOT::VecOps::RVec<double> result;
    for (const auto& v : tof1) result.push_back(v - path_corr);
    return result;
  }, {"tof_1"});

  df = df.Define("hittime_0_corrected", 
                 [path_corr](ROOT::VecOps::RVec<double> ht0) {
    ROOT::VecOps::RVec<double> result;
    for (const auto& v : ht0) result.push_back(v - path_corr);
    return result;
  }, {"hittime_0"});

  df = df.Define("hittime_1_corrected", 
                 [path_corr](ROOT::VecOps::RVec<double> ht1) {
    ROOT::VecOps::RVec<double> result;
    for (const auto& v : ht1) result.push_back(v - path_corr);
    return result;
  }, {"hittime_1"});

  df = df.Define("tof_1_corrected_front_veto", 
                 [path_corr](ROOT::VecOps::RVec<double> pd0, ROOT::VecOps::RVec<double> pd1, 
                            ROOT::VecOps::RVec<double> tof1) {
    ROOT::VecOps::RVec<double> result;
    for (size_t j = 0; j < pd1.size(); ++j) {
      bool has_front_hit = false;
      for (size_t i = 0; i < pd0.size(); ++i) {
        if (pd0[i] == pd1[j]) {
          has_front_hit = true;
          break;
        }
      }
      if (!has_front_hit) {
        result.push_back(tof1[j] - path_corr);
      }
    }
    return result;
  }, {"paddle_0", "paddle_1", "tof_1"});

  df = df.Define("hittime_1_corrected_front_veto", 
                 [path_corr](ROOT::VecOps::RVec<double> pd0, ROOT::VecOps::RVec<double> pd1, 
                            ROOT::VecOps::RVec<double> ht1) {
    ROOT::VecOps::RVec<double> result;
    for (size_t j = 0; j < pd1.size(); ++j) {
      bool has_front_hit = false;
      for (size_t i = 0; i < pd0.size(); ++i) {
        if (pd0[i] == pd1[j]) {
          has_front_hit = true;
          break;
        }
      }
      if (!has_front_hit) {
        result.push_back(ht1[j] - path_corr);
      }
    }
    return result;
  }, {"paddle_0", "paddle_1", "hittime_1"});

  // ----------------------------------------------------------------
  using RH1 = ROOT::RDF::RResultPtr<::TH1D>;
  using RH2 = ROOT::RDF::RResultPtr<::TH2D>;
  std::map<std::string, std::vector<std::vector<RH1>>> hpad; // [var][plane][paddle]
  std::map<std::string, std::vector<RH1>> hsum;              // [var][plane]
  std::map<std::string, RH1> htot;                           // [var]

  // Punch-through histograms
  std::map<std::string, std::vector<RH1>> hpt_pad; // [type][paddle] - type: "punchthrough"
  std::map<std::string, std::vector<RH1>> hpt_pad_01; // punchthrough for planes 0-1
  std::map<std::string, std::vector<RH1>> hpt_pad_23; // punchthrough for planes 2-3

  // Front veto histograms
  std::map<std::string, RH1> h_hittime_1_veto;  // hittime for back with no front hit
  std::map<std::string, RH1> h_tof_1_veto;      // tof for back with no front hit

  // Corrected tof and hittime histograms
  std::map<std::string, RH1> h_tof_0_corr;      // corrected tof side 0
  std::map<std::string, RH1> h_tof_1_corr;      // corrected tof side 1
  std::map<std::string, RH1> h_hittime_0_corr;  // corrected hittime side 0
  std::map<std::string, RH1> h_hittime_1_corr;  // corrected hittime side 1
  std::map<std::string, RH1> h_tof_1_corr_veto; // corrected tof side 1 with front veto
  std::map<std::string, RH1> h_hittime_1_corr_veto; // corrected hittime side 1 with front veto

  // 2D histograms for edep vs punch-through
  std::map<std::string, std::vector<RH2>> h_edep_vs_pt; // [side_edeptype][paddle]

  // 2D histograms for tof vs edep
  std::map<std::string, std::vector<RH2>> h_tof_vs_edep; // [side_edeptype][paddle]

  for (const auto &v : vars) {
    hpad[v.name].resize(N_PLANES);
    hsum[v.name].resize(N_PLANES);
    for (int p = 0; p < N_PLANES; ++p)
      hpad[v.name][p].resize(N_PADDLES);
  }

  for (const auto &v : vars) {
    for (int plane = 0; plane < N_PLANES; ++plane) {
      for (int paddle = 0; paddle < N_PADDLES; ++paddle) {
        const std::string col = v.name + "_p" + std::to_string(plane) + "_b" + std::to_string(paddle);
        const std::string title =
            v.name + " plane " + std::to_string(plane) + " paddle " + std::to_string(paddle) + ";" + v.name + ";Counts";
        hpad[v.name][plane][paddle] = df.Histo1D({col.c_str(), title.c_str(), v.nbins, v.xmin, v.xmax}, col);
      }
      const std::string col   = v.name + "_p" + std::to_string(plane) + "_sum";
      const std::string title = v.name + " plane " + std::to_string(plane) + " (all paddles);" + v.name + ";Counts";
      hsum[v.name][plane]     = df.Histo1D({col.c_str(), title.c_str(), v.nbins, v.xmin, v.xmax}, col);
    }
    const std::string col   = v.name + "_total";
    const std::string title = v.name + " all planes;" + v.name + ";Counts";
    htot[v.name]            = df.Histo1D({col.c_str(), title.c_str(), v.nbins, v.xmin, v.xmax}, col);
  }

  // Book punch-through histograms
  hpt_pad["punchthrough"].resize(N_PADDLES);
  for (int paddle = 0; paddle < N_PADDLES; ++paddle) {
    const std::string col   = "punchthrough_b" + std::to_string(paddle);
    const std::string title = "Punch-through time paddle " + std::to_string(paddle) + ";Punch-through time (ns);Counts";
    hpt_pad["punchthrough"][paddle] =
        df.Histo1D({col.c_str(), title.c_str(), NBINS_PUNCHTHROUGH, XMIN_PUNCHTHROUGH, XMAX_PUNCHTHROUGH}, col);
  }

  // Book punch-through 01 histograms (planes 0-1)
  hpt_pad_01["punchthrough_01"].resize(N_PADDLES);
  for (int paddle = 0; paddle < N_PADDLES; ++paddle) {
    const std::string col   = "punchthrough_01_b" + std::to_string(paddle);
    const std::string title = "Punch-through planes 0-1 paddle " + std::to_string(paddle) + ";Punch-through time (ns);Counts";
    hpt_pad_01["punchthrough_01"][paddle] =
        df.Histo1D({col.c_str(), title.c_str(), NBINS_PUNCHTHROUGH, XMIN_PUNCHTHROUGH, XMAX_PUNCHTHROUGH}, col);
  }

  // Book punch-through 23 histograms (planes 2-3)
  hpt_pad_23["punchthrough_23"].resize(N_PADDLES);
  for (int paddle = 0; paddle < N_PADDLES; ++paddle) {
    const std::string col   = "punchthrough_23_b" + std::to_string(paddle);
    const std::string title = "Punch-through planes 2-3 paddle " + std::to_string(paddle) + ";Punch-through time (ns);Counts";
    hpt_pad_23["punchthrough_23"][paddle] =
        df.Histo1D({col.c_str(), title.c_str(), NBINS_PUNCHTHROUGH, XMIN_PUNCHTHROUGH, XMAX_PUNCHTHROUGH}, col);
  }

  // Book front veto histograms
  {
    const std::string col_ht = "hittime_1_front_veto";
    const std::string title_ht = "hittime back (no front hit);hittime (ns);Counts";
    h_hittime_1_veto["hittime_1_veto"] = df.Histo1D({col_ht.c_str(), title_ht.c_str(), 
                                                      NBINS_HITTIME, XMIN_HITTIME, XMAX_HITTIME}, col_ht);

    const std::string col_tof = "tof_1_front_veto";
    const std::string title_tof = "tof back (no front hit);tof (ns);Counts";
    h_tof_1_veto["tof_1_veto"] = df.Histo1D({col_tof.c_str(), title_tof.c_str(), 
                                             NBINS_TOF, XMIN_TOF, XMAX_TOF}, col_tof);
  }

  // Book corrected tof and hittime histograms
  {
    const std::string col_tof0 = "tof_0_corrected";
    const std::string title_tof0 = "tof (photon corrected) side 0;tof (ns);Counts";
    h_tof_0_corr["tof_0_corr"] = df.Histo1D({col_tof0.c_str(), title_tof0.c_str(), 
                                             NBINS_TOF, XMIN_TOF, XMAX_TOF}, col_tof0);

    const std::string col_tof1 = "tof_1_corrected";
    const std::string title_tof1 = "tof (photon corrected) side 1;tof (ns);Counts";
    h_tof_1_corr["tof_1_corr"] = df.Histo1D({col_tof1.c_str(), title_tof1.c_str(), 
                                             NBINS_TOF, XMIN_TOF, XMAX_TOF}, col_tof1);

    const std::string col_ht0 = "hittime_0_corrected";
    const std::string title_ht0 = "hittime (photon corrected) side 0;hittime (ns);Counts";
    h_hittime_0_corr["hittime_0_corr"] = df.Histo1D({col_ht0.c_str(), title_ht0.c_str(), 
                                                      NBINS_HITTIME, XMIN_HITTIME, XMAX_HITTIME}, col_ht0);

    const std::string col_ht1 = "hittime_1_corrected";
    const std::string title_ht1 = "hittime (photon corrected) side 1;hittime (ns);Counts";
    h_hittime_1_corr["hittime_1_corr"] = df.Histo1D({col_ht1.c_str(), title_ht1.c_str(), 
                                                      NBINS_HITTIME, XMIN_HITTIME, XMAX_HITTIME}, col_ht1);

    const std::string col_tof1_veto = "tof_1_corrected_front_veto";
    const std::string title_tof1_veto = "tof (photon corrected) back no front;tof (ns);Counts";
    h_tof_1_corr_veto["tof_1_corr_veto"] = df.Histo1D({col_tof1_veto.c_str(), title_tof1_veto.c_str(), 
                                                        NBINS_TOF, XMIN_TOF, XMAX_TOF}, col_tof1_veto);

    const std::string col_ht1_veto = "hittime_1_corrected_front_veto";
    const std::string title_ht1_veto = "hittime (photon corrected) back no front;hittime (ns);Counts";
    h_hittime_1_corr_veto["hittime_1_corr_veto"] = df.Histo1D({col_ht1_veto.c_str(), title_ht1_veto.c_str(), 
                                                                 NBINS_HITTIME, XMIN_HITTIME, XMAX_HITTIME}, col_ht1_veto);
  }

  // Book 2D histograms for edep vs punch-through
  h_edep_vs_pt["edep_0"].resize(N_PADDLES);
  h_edep_vs_pt["edep_1"].resize(N_PADDLES);
  h_edep_vs_pt["edep_amp_0"].resize(N_PADDLES);
  h_edep_vs_pt["edep_amp_1"].resize(N_PADDLES);

  for (int paddle = 0; paddle < N_PADDLES; ++paddle) {
    // edep_0 vs punchthrough
    {
      const std::string x_col = "pt_0_b" + std::to_string(paddle);
      const std::string y_col = "edep_0_b_pt" + std::to_string(paddle);
      const std::string title =
          "edep vs punch-through (side 0, paddle " + std::to_string(paddle) + ");Punch-through time (ns);edep (MeV)";
      h_edep_vs_pt["edep_0"][paddle] =
          df.Histo2D({("h_edep_0_vs_pt_" + std::to_string(paddle)).c_str(), title.c_str(), NBINS_PUNCHTHROUGH,
                      XMIN_PUNCHTHROUGH, XMAX_PUNCHTHROUGH, NBINS_EDEP, XMIN_EDEP, XMAX_EDEP},
                     x_col, y_col);
    }

    // edep_1 vs punchthrough
    {
      const std::string x_col = "pt_1_b" + std::to_string(paddle);
      const std::string y_col = "edep_1_b_pt" + std::to_string(paddle);
      const std::string title =
          "edep vs punch-through (side 1, paddle " + std::to_string(paddle) + ");Punch-through time (ns);edep (MeV)";
      h_edep_vs_pt["edep_1"][paddle] =
          df.Histo2D({("h_edep_1_vs_pt_" + std::to_string(paddle)).c_str(), title.c_str(), NBINS_PUNCHTHROUGH,
                      XMIN_PUNCHTHROUGH, XMAX_PUNCHTHROUGH, NBINS_EDEP, XMIN_EDEP, XMAX_EDEP},
                     x_col, y_col);
    }

    // edep_amp_0 vs punchthrough
    {
      const std::string x_col = "ptamp_0_b" + std::to_string(paddle);
      const std::string y_col = "edepamp_0_b_pt" + std::to_string(paddle);
      const std::string title =
          "edep_amp vs punch-through (side 0, paddle " + std::to_string(paddle) + ");Punch-through time (ns);edep_amp";
      h_edep_vs_pt["edep_amp_0"][paddle] =
          df.Histo2D({("h_edepamp_0_vs_pt_" + std::to_string(paddle)).c_str(), title.c_str(), NBINS_PUNCHTHROUGH,
                      XMIN_PUNCHTHROUGH, XMAX_PUNCHTHROUGH, NBINS_EDEP_AMP, XMIN_EDEP_AMP, XMAX_EDEP_AMP},
                     x_col, y_col);
    }

    // edep_amp_1 vs punchthrough
    {
      const std::string x_col = "ptamp_1_b" + std::to_string(paddle);
      const std::string y_col = "edepamp_1_b_pt" + std::to_string(paddle);
      const std::string title =
          "edep_amp vs punch-through (side 1, paddle " + std::to_string(paddle) + ");Punch-through time (ns);edep_amp";
      h_edep_vs_pt["edep_amp_1"][paddle] =
          df.Histo2D({("h_edepamp_1_vs_pt_" + std::to_string(paddle)).c_str(), title.c_str(), NBINS_PUNCHTHROUGH,
                      XMIN_PUNCHTHROUGH, XMAX_PUNCHTHROUGH, NBINS_EDEP_AMP, XMIN_EDEP_AMP, XMAX_EDEP_AMP},
                     x_col, y_col);
    }
  }

  // Book 2D histograms for tof vs edep
  h_tof_vs_edep["edep_0"].resize(N_PADDLES);
  h_tof_vs_edep["edep_1"].resize(N_PADDLES);
  h_tof_vs_edep["edep_amp_0"].resize(N_PADDLES);
  h_tof_vs_edep["edep_amp_1"].resize(N_PADDLES);

  for (int paddle = 0; paddle < N_PADDLES; ++paddle) {
    // tof vs edep_0
    {
      const std::string x_col = "tof_0_b" + std::to_string(paddle);
      const std::string y_col = "edep_0_b_tof" + std::to_string(paddle);
      const std::string title = "edep vs tof (side 0, paddle " + std::to_string(paddle) + ");tof (ns);edep (MeV)";
      h_tof_vs_edep["edep_0"][paddle] = df.Histo2D({("h_tof_edep_0_" + std::to_string(paddle)).c_str(), title.c_str(),
                                                    NBINS_TOF, XMIN_TOF, XMAX_TOF, NBINS_EDEP, XMIN_EDEP, XMAX_EDEP},
                                                   x_col, y_col);
    }

    // tof vs edep_1
    {
      const std::string x_col = "tof_1_b" + std::to_string(paddle);
      const std::string y_col = "edep_1_b_tof" + std::to_string(paddle);
      const std::string title = "edep vs tof (side 1, paddle " + std::to_string(paddle) + ");tof (ns);edep (MeV)";
      h_tof_vs_edep["edep_1"][paddle] = df.Histo2D({("h_tof_edep_1_" + std::to_string(paddle)).c_str(), title.c_str(),
                                                    NBINS_TOF, XMIN_TOF, XMAX_TOF, NBINS_EDEP, XMIN_EDEP, XMAX_EDEP},
                                                   x_col, y_col);
    }

    // tof vs edep_amp_0
    {
      const std::string x_col = "tofamp_0_b" + std::to_string(paddle);
      const std::string y_col = "edepamp_0_b_tof" + std::to_string(paddle);
      const std::string title = "edep_amp vs tof (side 0, paddle " + std::to_string(paddle) + ");tof (ns);edep_amp";
      h_tof_vs_edep["edep_amp_0"][paddle] =
          df.Histo2D({("h_tofamp_edep_0_" + std::to_string(paddle)).c_str(), title.c_str(), NBINS_TOF, XMIN_TOF,
                      XMAX_TOF, NBINS_EDEP_AMP, XMIN_EDEP_AMP, XMAX_EDEP_AMP},
                     x_col, y_col);
    }

    // tof vs edep_amp_1
    {
      const std::string x_col = "tofamp_1_b" + std::to_string(paddle);
      const std::string y_col = "edepamp_1_b_tof" + std::to_string(paddle);
      const std::string title = "edep_amp vs tof (side 1, paddle " + std::to_string(paddle) + ");tof (ns);edep_amp";
      h_tof_vs_edep["edep_amp_1"][paddle] =
          df.Histo2D({("h_tofamp_edep_1_" + std::to_string(paddle)).c_str(), title.c_str(), NBINS_TOF, XMIN_TOF,
                      XMAX_TOF, NBINS_EDEP_AMP, XMIN_EDEP_AMP, XMAX_EDEP_AMP},
                     x_col, y_col);
    }
  }

  // ----------------------------------------------------------------
  // 5. Trigger the event loop once, then write canvases.
  // ----------------------------------------------------------------
  std::cout << "[lad_tof_fast] Running event loop (IMT enabled)...\n";
  (void)htot[vars.front().name]->GetEntries(); // force materialization
#ifndef LAD_HAS_RDF_PROGRESSBAR
  std::fprintf(stderr, "\n"); // close the manual progress line
#endif

  TFile fout(out_file, "RECREATE");
  if (fout.IsZombie()) {
    std::cerr << "[lad_tof_fast] ERROR: cannot open output '" << out_file << "'\n";
    return;
  }

  for (const auto &v : vars) {
    TDirectory *dir = fout.mkdir(v.name.c_str());
    dir->cd();

    // 5 per-plane canvases, each with 11 paddle pads
    for (int plane = 0; plane < N_PLANES; ++plane) {
      const std::string cname  = "c_" + v.name + "_plane" + std::to_string(plane);
      const std::string ctitle = v.name + " plane " + std::to_string(plane);
      TCanvas *c               = new TCanvas(cname.c_str(), ctitle.c_str(), 1600, 1000);
      c->Divide(4, 3); // 12 pads; we use the first 11
      for (int paddle = 0; paddle < N_PADDLES; ++paddle) {
        c->cd(paddle + 1);
        hpad[v.name][plane][paddle]->DrawCopy();
      }
      c->Write();
      delete c;
    }

    // summary canvas: 5 plane sums + grand total
    const std::string cname  = "c_" + v.name + "_summary";
    const std::string ctitle = v.name + " summary (plane sums + total)";
    TCanvas *c               = new TCanvas(cname.c_str(), ctitle.c_str(), 1600, 1000);
    c->Divide(3, 2);
    for (int plane = 0; plane < N_PLANES; ++plane) {
      c->cd(plane + 1);
      hsum[v.name][plane]->DrawCopy();
    }
    c->cd(6);
    htot[v.name]->DrawCopy();
    c->Write();
    delete c;

    fout.cd();
  }

  // Write punch-through histograms
  {
    TDirectory *dir = fout.mkdir("punchthrough");
    dir->cd();

    TCanvas *c = new TCanvas("c_punchthrough_summary", "Punch-through time summary", 1600, 1000);
    c->Divide(4, 3);
    for (int paddle = 0; paddle < N_PADDLES; ++paddle) {
      c->cd(paddle + 1);
      hpt_pad["punchthrough"][paddle]->DrawCopy();
    }
    c->Write();
    delete c;

    fout.cd();
  }

  // Write edep vs punch-through 2D histograms
  {
    TDirectory *dir = fout.mkdir("edep_vs_punchthrough");
    dir->cd();

    // Side 0 - edep
    {
      TCanvas *c = new TCanvas("c_edep_0_vs_pt", "edep vs punch-through (side 0)", 1600, 1000);
      c->Divide(4, 3);
      for (int paddle = 0; paddle < N_PADDLES; ++paddle) {
        c->cd(paddle + 1);
        h_edep_vs_pt["edep_0"][paddle]->DrawCopy("COLZ");
      }
      c->Write();
      delete c;
    }

    // Side 1 - edep
    {
      TCanvas *c = new TCanvas("c_edep_1_vs_pt", "edep vs punch-through (side 1)", 1600, 1000);
      c->Divide(4, 3);
      for (int paddle = 0; paddle < N_PADDLES; ++paddle) {
        c->cd(paddle + 1);
        h_edep_vs_pt["edep_1"][paddle]->DrawCopy("COLZ");
      }
      c->Write();
      delete c;
    }

    // Side 0 - edep_amp
    {
      TCanvas *c = new TCanvas("c_edepamp_0_vs_pt", "edep_amp vs punch-through (side 0)", 1600, 1000);
      c->Divide(4, 3);
      for (int paddle = 0; paddle < N_PADDLES; ++paddle) {
        c->cd(paddle + 1);
        h_edep_vs_pt["edep_amp_0"][paddle]->DrawCopy("COLZ");
      }
      c->Write();
      delete c;
    }

    // Side 1 - edep_amp
    {
      TCanvas *c = new TCanvas("c_edepamp_1_vs_pt", "edep_amp vs punch-through (side 1)", 1600, 1000);
      c->Divide(4, 3);
      for (int paddle = 0; paddle < N_PADDLES; ++paddle) {
        c->cd(paddle + 1);
        h_edep_vs_pt["edep_amp_1"][paddle]->DrawCopy("COLZ");
      }
      c->Write();
      delete c;
    }

    fout.cd();
  }

  // Write tof vs edep 2D histograms
  {
    TDirectory *dir = fout.mkdir("tof_vs_edep");
    dir->cd();

    // Side 0 - edep
    {
      TCanvas *c = new TCanvas("c_tof_vs_edep_0", "tof vs edep (side 0)", 1600, 1000);
      c->Divide(4, 3);
      for (int paddle = 0; paddle < N_PADDLES; ++paddle) {
        c->cd(paddle + 1);
        h_tof_vs_edep["edep_0"][paddle]->DrawCopy("COLZ");
      }
      c->Write();
      delete c;
    }

    // Side 1 - edep
    {
      TCanvas *c = new TCanvas("c_tof_vs_edep_1", "tof vs edep (side 1)", 1600, 1000);
      c->Divide(4, 3);
      for (int paddle = 0; paddle < N_PADDLES; ++paddle) {
        c->cd(paddle + 1);
        h_tof_vs_edep["edep_1"][paddle]->DrawCopy("COLZ");
      }
      c->Write();
      delete c;
    }

    // Side 0 - edep_amp
    {
      TCanvas *c = new TCanvas("c_tof_vs_edepamp_0", "tof vs edep_amp (side 0)", 1600, 1000);
      c->Divide(4, 3);
      for (int paddle = 0; paddle < N_PADDLES; ++paddle) {
        c->cd(paddle + 1);
        h_tof_vs_edep["edep_amp_0"][paddle]->DrawCopy("COLZ");
      }
      c->Write();
      delete c;
    }

    // Side 1 - edep_amp
    {
      TCanvas *c = new TCanvas("c_tof_vs_edepamp_1", "tof vs edep_amp (side 1)", 1600, 1000);
      c->Divide(4, 3);
      for (int paddle = 0; paddle < N_PADDLES; ++paddle) {
        c->cd(paddle + 1);
        h_tof_vs_edep["edep_amp_1"][paddle]->DrawCopy("COLZ");
      }
      c->Write();
      delete c;
    }

    fout.cd();
  }

  // Write punch-through 01 and 23 canvases
  {
    TDirectory *dir = fout.mkdir("punchthrough_planes");
    dir->cd();

    TCanvas *c_01 = new TCanvas("c_punchthrough_01", "Punch-through planes 0-1", 1600, 1000);
    c_01->Divide(4, 3);
    for (int paddle = 0; paddle < N_PADDLES; ++paddle) {
      c_01->cd(paddle + 1);
      hpt_pad_01["punchthrough_01"][paddle]->DrawCopy();
    }
    c_01->Write();
    delete c_01;

    TCanvas *c_23 = new TCanvas("c_punchthrough_23", "Punch-through planes 2-3", 1600, 1000);
    c_23->Divide(4, 3);
    for (int paddle = 0; paddle < N_PADDLES; ++paddle) {
      c_23->cd(paddle + 1);
      hpt_pad_23["punchthrough_23"][paddle]->DrawCopy();
    }
    c_23->Write();
    delete c_23;

    fout.cd();
  }

  // Write front veto histograms
  {
    TDirectory *dir = fout.mkdir("front_veto");
    dir->cd();

    TCanvas *c = new TCanvas("c_front_veto", "Front veto (back with no front hit)", 1200, 600);
    c->Divide(2, 1);
    c->cd(1);
    h_hittime_1_veto["hittime_1_veto"]->DrawCopy();
    c->cd(2);
    h_tof_1_veto["tof_1_veto"]->DrawCopy();
    c->Write();
    delete c;

    fout.cd();
  }

  // Write corrected tof and hittime histograms
  {
    TDirectory *dir = fout.mkdir("corrected");
    dir->cd();

    TCanvas *c_tof = new TCanvas("c_tof_corrected", "TOF with photon path correction", 1200, 600);
    c_tof->Divide(2, 1);
    c_tof->cd(1);
    h_tof_0_corr["tof_0_corr"]->DrawCopy();
    c_tof->cd(2);
    h_tof_1_corr["tof_1_corr"]->DrawCopy();
    c_tof->Write();
    delete c_tof;

    TCanvas *c_ht = new TCanvas("c_hittime_corrected", "Hittime with photon path correction", 1200, 600);
    c_ht->Divide(2, 1);
    c_ht->cd(1);
    h_hittime_0_corr["hittime_0_corr"]->DrawCopy();
    c_ht->cd(2);
    h_hittime_1_corr["hittime_1_corr"]->DrawCopy();
    c_ht->Write();
    delete c_ht;

    TCanvas *c_veto = new TCanvas("c_corrected_veto", "Corrected with front veto", 1200, 600);
    c_veto->Divide(2, 1);
    c_veto->cd(1);
    h_tof_1_corr_veto["tof_1_corr_veto"]->DrawCopy();
    c_veto->cd(2);
    h_hittime_1_corr_veto["hittime_1_corr_veto"]->DrawCopy();
    c_veto->Write();
    delete c_veto;

    fout.cd();
  }

  fout.Close();
  std::cout << "[lad_tof_fast] Done. Wrote " << out_file << "\n";
}