// make_RF_timing_plots.C
//
// RF-timing peak analysis, one output ROOT file for a whole list of runs.
//
// For every ROOT file listed in a .dat file this macro builds the
// path-length-corrected, RF-corrected time-of-flight distribution (analogous to
// the "P_c_tof_corr_summary" canvas produced by lad_tof/lad_tof_fast.C, but
// using the goodhit_hit_tof_rfcorr_* branches instead of goodhit_hit_tof_*),
// then:
//   * writes a per-file summary canvas (planes 000,001,100,101,200 + all-plane
//     total) with a Gaussian fit drawn on the peak of each pad,
//   * finds the tallest bin inside the coincidence-peak search band
//     [10, 30] ns and fits a Gaussian over [peak-1.5, peak+1.5] ns,
//   * accumulates the fitted peak CENTER (mean) and peak WIDTH (sigma), with
//     fit errors, versus run number.
//
// The run number is the first 5-digit token in the file name (e.g.
// LAD_COIN_22614_1_1_-1.root -> 22614).  Several files may share a run number;
// each contributes its own point (so a run can appear multiple times on x).
//
// Both spectrometers (P and H) are analysed.  Hits are taken with no tracking
// cut (the tracking-independent "all_hits" category of lad_tof_fast).
//
// Output ROOT file layout:
//   per_file/<stem>/{P,H}_c_tof_corr_summary   -- summary canvas per input file
//   per_file/<stem>/<spec>_tofcorr_<label>     -- the fitted histograms
//   summary/{P,H}_c_peak_center                -- center-vs-run, 6 pads
//   summary/{P,H}_c_peak_width                 -- width-vs-run,  6 pads
//
// Usage:
//   root -l -b -q 'make_RF_timing_plots.C("files.dat","out.root")'
//   root -l -b -q 'make_RF_timing_plots.C("files.dat","out.root",8)'  // 8 threads
// The 3rd arg caps the RDataFrame worker threads (default 4); <=0 uses all cores.
// Each input file is processed with its own RDataFrame under implicit MT, so the
// per-file event loop is parallelised across cores while memory stays bounded to
// one file at a time.

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TCanvas.h>
#include <TDirectory.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include <cstdio>

#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <regex>
#include <string>
#include <vector>

// =====================================================================
// Geometry / binning constants (kept identical to lad_tof_fast.C so the
// distribution is directly comparable to P_c_tof_corr_summary)
// =====================================================================
const int    N_PLANES = 5, N_PADDLES = 11, N_SPECS = 2;
const int    N_TARGETS = N_PLANES + 1;               // 5 planes + all-plane total
const double hodo_radii[N_PLANES] = {615., 655.6, 523., 563.6, 615.}; // cm
const char *const plane_names[N_PLANES] = {"000", "001", "100", "101", "200"};
const char specs[N_SPECS] = {'P', 'H'};

const int    NBINS_TCORR = 900;
const double XMIN_TCORR = -150., XMAX_TCORR = 300.;
const double FIT_HALF_WINDOW = 1.5;   // Gaussian fit range = peak +/- this (ns)
const double MIN_ENTRIES_FOR_FIT = 50.;
// The RF-corrected coincidence peak sits at ~20 ns in tof_rfcorr - L/c; search
// for the maximum only inside this band so the fit cannot latch onto a taller
// noise bin elsewhere in the (mostly flat) background.
const double PEAK_SEARCH_MIN = 10.;
const double PEAK_SEARCH_MAX = 30.;

const char *DEFAULT_DAT_FILE =
  "/Users/lucasehinger/Library/CloudStorage/OneDrive-Personal/Documents/MIT/"
  "Research/LAD/LAD_coding/lad_replay/CALIBRATION/lad_tof/files/all_C3_test.dat";
const char *DEFAULT_OUT_FILE = "RF_timing_plots.root";

// label for target index t: 0..4 -> plane name, 5 -> "total"
static std::string target_label(int t) {
  return (t < N_PLANES) ? std::string(plane_names[t]) : std::string("total");
}

// result of a single Gaussian peak fit
struct PeakFit {
  bool   ok = false;
  double center = 0., center_err = 0.;   // Gaussian mean      (par 1)
  double height = 0., height_err = 0.;   // Gaussian amplitude (par 0)
  double width = 0., width_err = 0.;     // Gaussian sigma     (par 2)
};

// =====================================================================
// Extract the 5-digit run number from a file path (first such token).
// Returns -1 if none found.
// =====================================================================
static int parse_run_number(const std::string &path) {
  std::string base = gSystem->BaseName(path.c_str());
  std::smatch m;
  static const std::regex re("([0-9]{5})");
  if (std::regex_search(base, m, re)) return std::stoi(m[1].str());
  return -1;
}

// file-name stem (basename without .root) for the per-file output directory
static std::string file_stem(const std::string &path) {
  std::string base = gSystem->BaseName(path.c_str());
  size_t dot = base.rfind(".root");
  if (dot != std::string::npos) base = base.substr(0, dot);
  return base;
}

// =====================================================================
// Analyse one file: build the 12 corrected-tof histograms (6 per spec),
// fit the peaks, write the per-file summary canvases, and return the fits.
// =====================================================================
static bool analyze_file(const std::string &path, TDirectory *perfile,
                         const std::string &stem,
                         std::array<std::array<PeakFit, N_TARGETS>, N_SPECS> &fits) {
  using RVd = ROOT::VecOps::RVec<double>;

  // Guard against an unreadable file / missing tree.
  {
    TFile ftest(path.c_str(), "READ");
    if (ftest.IsZombie() || !ftest.Get("T")) {
      std::cerr << "[make_RF_timing_plots] cannot open tree T in " << path << " -- skipping\n";
      return false;
    }
  }

  ROOT::RDataFrame rdf("T", path);
  ROOT::RDF::RNode df = rdf;

  // --- column definitions: per-plane corrected-tof sums + all-plane total ---
  for (int is = 0; is < N_SPECS; ++is) {
    const std::string sp(1, specs[is]);
    const std::string pfx = sp + ".ladhod.goodhit_";
    std::vector<std::string> plane_cols;
    for (int pl = 0; pl < N_PLANES; ++pl) {
      const std::string side = (pl % 2 == 0) ? "0" : "1";
      const double R = hodo_radii[pl];
      const double plv = static_cast<double>(pl);
      const bool excl = (pl == 2 || pl == 3);   // planes 100/101 drop paddles 1 and 9
      const std::string col = sp + "_tofcorr_p" + std::to_string(pl) + "_sum";
      df = df.Define(col,
        [R, plv, excl](const RVd &pln, const RVd &pdl, const RVd &yp, const RVd &tf) {
          RVd r;
          for (size_t i = 0; i < pln.size(); ++i) {
            if (pln[i] != plv) continue;
            if (excl && (pdl[i] == 1. || pdl[i] == 9.)) continue;
            const double dx  = 22. * (pdl[i] - 6.);
            const double p2d = std::sqrt(yp[i] * yp[i] + dx * dx);
            r.push_back(tf[i] - std::sqrt(p2d * p2d + R * R) / 100. / 0.3);
          }
          return r;
        },
        {pfx + "plane_" + side, pfx + "paddle_" + side,
         pfx + "hit_ypos_" + side, pfx + "hit_tof_rfcorr_" + side});
      plane_cols.push_back(col);
    }
    df = df.Define(sp + "_tofcorr_total",
      [](const RVd &a, const RVd &b, const RVd &c, const RVd &d, const RVd &e) {
        auto r = ROOT::VecOps::Concatenate(a, b);
        r = ROOT::VecOps::Concatenate(r, c);
        r = ROOT::VecOps::Concatenate(r, d);
        return ROOT::VecOps::Concatenate(r, e);
      },
      plane_cols);
  }

  // --- book histograms (lazy) ---
  std::array<std::array<ROOT::RDF::RResultPtr<TH1D>, N_TARGETS>, N_SPECS> hbk;
  for (int is = 0; is < N_SPECS; ++is) {
    const std::string sp(1, specs[is]);
    for (int t = 0; t < N_TARGETS; ++t) {
      const std::string col = (t < N_PLANES)
        ? sp + "_tofcorr_p" + std::to_string(t) + "_sum"
        : sp + "_tofcorr_total";
      const std::string hn = sp + "_tofcorr_" + target_label(t) + "_" + stem;
      const std::string ht = sp + " tof corr " + target_label(t) +
                             " (" + stem + ");tof_{rfcorr}-L/c (ns);Counts";
      hbk[is][t] = df.Histo1D({hn.c_str(), ht.c_str(), NBINS_TCORR, XMIN_TCORR, XMAX_TCORR}, col);
    }
  }

  // Trigger the (multi-threaded) event loop once for this file.
  (void)hbk[0][0]->GetEntries();

  // --- fit peaks + draw per-file summary canvases ---
  TDirectory *fdir = perfile->mkdir(stem.c_str());
  if (!fdir) fdir = perfile;   // fall back if the stem already exists
  fdir->cd();

  for (int is = 0; is < N_SPECS; ++is) {
    const std::string sp(1, specs[is]);
    std::array<TH1D *, N_TARGETS> hc{};
    std::array<TF1 *, N_TARGETS> fc{};
    std::array<TPaveText *, N_TARGETS> pt{};

    for (int t = 0; t < N_TARGETS; ++t) {
      TH1D *h = static_cast<TH1D *>(hbk[is][t]->Clone());
      h->SetDirectory(nullptr);
      hc[t] = h;

      PeakFit pf;
      if (h->GetEntries() >= MIN_ENTRIES_FOR_FIT && h->Integral() > 0.) {
        // Find the tallest bin only within the search band, then restore the
        // full axis range for drawing.
        h->GetXaxis()->SetRangeUser(PEAK_SEARCH_MIN, PEAK_SEARCH_MAX);
        const int    mbin  = h->GetMaximumBin();
        const double xpeak = h->GetBinCenter(mbin);
        const double ypeak = h->GetBinContent(mbin);
        h->GetXaxis()->SetRange(0, 0);       // reset zoom -> full range
        TF1 *g = new TF1((sp + "_gaus_" + target_label(t) + "_" + stem).c_str(), "gaus",
                         xpeak - FIT_HALF_WINDOW, xpeak + FIT_HALF_WINDOW);
        g->SetParameters(ypeak, xpeak, 0.5);
        g->SetParLimits(2, 1e-3, 20.);       // keep sigma positive/bounded
        g->SetLineColor(kRed);
        const int status = h->Fit(g, "RQ");  // R: use g's range; Q: quiet
        fc[t] = g;
        if (status == 0) {
          pf.ok         = true;
          pf.center     = g->GetParameter(1);
          pf.center_err = g->GetParError(1);
          pf.height     = g->GetParameter(0);
          pf.height_err = g->GetParError(0);
          pf.width      = g->GetParameter(2);
          pf.width_err  = g->GetParError(2);
        }
      }
      fits[is][t] = pf;

      // Stamp the fit results onto the pad so they travel with the file
      // (independent of the viewer's gStyle OptFit setting).
      TPaveText *box = new TPaveText(0.13, 0.72, 0.45, 0.89, "NDC");
      box->SetFillColor(0);
      box->SetBorderSize(1);
      box->SetTextAlign(12);
      box->SetTextSize(0.035);
      char line[128];
      if (pf.ok) {
        std::snprintf(line, sizeof line, "peak @ %.2f ns", fc[t]->GetParameter(1)); box->AddText(line);
        std::snprintf(line, sizeof line, "height = %.1f #pm %.1f", pf.height, pf.height_err); box->AddText(line);
        std::snprintf(line, sizeof line, "#sigma = %.3f #pm %.3f ns", pf.width, pf.width_err); box->AddText(line);
      } else {
        box->AddText("fit failed / too few entries");
      }
      pt[t] = box;
    }

    // Summary canvas: 3x2 -> planes 000,001,100,101,200 then all-plane total.
    TCanvas *cs = new TCanvas((sp + "_c_tof_corr_summary").c_str(),
                              (sp + " tof corr summary (" + stem + ")").c_str(), 1600, 1000);
    cs->Divide(3, 2);
    for (int t = 0; t < N_TARGETS; ++t) {
      cs->cd(t + 1);
      hc[t]->Draw();
      if (fc[t]) fc[t]->Draw("SAME");
      pt[t]->Draw();
    }
    fdir->cd();
    cs->Write();
    for (int t = 0; t < N_TARGETS; ++t) hc[t]->Write();  // keep the fitted histos too
    delete cs;
    for (int t = 0; t < N_TARGETS; ++t) { delete pt[t]; delete fc[t]; delete hc[t]; }
  }
  return true;
}

// =====================================================================
// Main
// =====================================================================
void make_RF_timing_plots(const char *dat_file = DEFAULT_DAT_FILE,
                          const char *out_file = DEFAULT_OUT_FILE,
                          int nthreads = 4) {
  gROOT->SetBatch(kTRUE);
  TH1::AddDirectory(kFALSE);
  gStyle->SetOptFit(1111);
  gStyle->SetOptStat(1110);

  if (nthreads > 0) { ROOT::EnableImplicitMT(nthreads);
    std::cout << "[make_RF_timing_plots] implicit MT: " << nthreads << " threads\n"; }
  else              { ROOT::EnableImplicitMT();
    std::cout << "[make_RF_timing_plots] implicit MT: all cores\n"; }

  // ---- read the file list ----
  std::vector<std::string> files;
  {
    std::ifstream fin(dat_file);
    if (!fin.is_open()) { std::cerr << "cannot open list " << dat_file << "\n"; return; }
    std::string ln;
    while (std::getline(fin, ln)) {
      size_t a = ln.find_first_not_of(" \t\r\n");
      if (a == std::string::npos) continue;
      std::string p = ln.substr(a, ln.find_last_not_of(" \t\r\n") - a + 1);
      if (p.empty() || p[0] == '#') continue;
      files.push_back(p);
    }
  }
  if (files.empty()) { std::cerr << "no files in " << dat_file << "\n"; return; }
  std::cout << "[make_RF_timing_plots] " << files.size() << " file(s) to process\n";

  TFile fout(out_file, "RECREATE");
  if (fout.IsZombie()) { std::cerr << "cannot open output " << out_file << "\n"; return; }
  TDirectory *perfile = fout.mkdir("per_file");

  // per [spec][target] accumulators of (run, value, error): peak center & width
  std::array<std::array<std::vector<double>, N_TARGETS>, N_SPECS> gx, gyc, gec, gyw, gew;

  // ---- process each file ----
  for (size_t i = 0; i < files.size(); ++i) {
    const std::string &path = files[i];
    const int run = parse_run_number(path);
    const std::string stem = file_stem(path);
    std::cout << "[make_RF_timing_plots] (" << (i + 1) << "/" << files.size()
              << ") run " << run << "  " << stem << std::endl;
    if (run < 0)
      std::cerr << "[make_RF_timing_plots] WARNING: no 5-digit run number in " << stem << "\n";

    std::array<std::array<PeakFit, N_TARGETS>, N_SPECS> fits;
    if (!analyze_file(path, perfile, stem, fits)) continue;

    for (int is = 0; is < N_SPECS; ++is)
      for (int t = 0; t < N_TARGETS; ++t) {
        const PeakFit &pf = fits[is][t];
        if (!pf.ok || run < 0) continue;
        gx[is][t].push_back(run);
        gyc[is][t].push_back(pf.center); gec[is][t].push_back(pf.center_err);
        gyw[is][t].push_back(pf.width);  gew[is][t].push_back(pf.width_err);
      }
  }

  // ---- build peak-center / peak-width vs run-number graphs ----
  TDirectory *sdir = fout.mkdir("summary");

  auto make_graph = [](const std::vector<double> &x, const std::vector<double> &y,
                       const std::vector<double> &ey, const std::string &name,
                       const std::string &title) -> TGraphErrors * {
    const int n = static_cast<int>(x.size());
    std::vector<double> ex(n, 0.);
    TGraphErrors *g = (n > 0)
      ? new TGraphErrors(n, const_cast<double *>(x.data()), const_cast<double *>(y.data()),
                         ex.data(), const_cast<double *>(ey.data()))
      : new TGraphErrors();
    g->SetName(name.c_str());
    g->SetTitle(title.c_str());
    g->SetMarkerStyle(20);
    g->SetMarkerSize(0.9);
    g->SetMarkerColor(kBlue + 1);
    g->SetLineColor(kBlue + 1);
    return g;
  };

  for (int is = 0; is < N_SPECS; ++is) {
    const std::string sp(1, specs[is]);

    TCanvas *cc = new TCanvas((sp + "_c_peak_center").c_str(),
                              (sp + " peak center vs run").c_str(), 1600, 1000);
    TCanvas *cw = new TCanvas((sp + "_c_peak_width").c_str(),
                              (sp + " peak width vs run").c_str(), 1600, 1000);
    cc->Divide(3, 2);
    cw->Divide(3, 2);

    std::vector<TGraphErrors *> owned;   // graphs live on the canvas pads until written
    for (int t = 0; t < N_TARGETS; ++t) {
      const std::string lab = target_label(t);
      TGraphErrors *gc = make_graph(gx[is][t], gyc[is][t], gec[is][t],
        sp + "_peakcenter_" + lab, sp + " peak center " + lab + ";run number;peak center (ns)");
      TGraphErrors *gw = make_graph(gx[is][t], gyw[is][t], gew[is][t],
        sp + "_peakwidth_" + lab, sp + " peak width " + lab + ";run number;peak width #sigma (ns)");

      cc->cd(t + 1); if (gc->GetN() > 0) gc->Draw("AP");
      cw->cd(t + 1); if (gw->GetN() > 0) gw->Draw("AP");
      owned.push_back(gc);
      owned.push_back(gw);
    }
    // Only the combined canvases are stored (individual TGraphs are not written).
    sdir->cd();
    cc->Write();
    cw->Write();
    delete cc;
    delete cw;
    for (auto *g : owned) delete g;
  }

  fout.Close();
  std::cout << "[make_RF_timing_plots] Done. Wrote " << out_file << "\n";
}
