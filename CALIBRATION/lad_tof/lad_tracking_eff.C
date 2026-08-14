// lad_tracking_eff.C
// -------------------------------------------------------------------------
// Tracking-efficiency macro. Derived from lad_tof_fast.C but produces ONLY the
// proton-ID plots (<spec>/proton_id/<variant>/...), for every tracking variant
// present in the data -- now including the new 1D-cluster projective variants.
//
// For speed, none of the tof / ypos / hittime / edep / punchthrough / front_veto
// / corrected / bgsub plot families from lad_tof_fast.C are reproduced; only the
// proton-tagged corrected-tof histograms and their track/total efficiency ratios.
//
// Tracking variants (per-hit chiSquare branches, X = P or H):
//   standard          -> X.ladhod.goodhit_chiSquare
//   xz                -> X.ladhod.goodhit_chiSquare_xz
//   noTrackVertex     -> X.ladhod.goodhit_chiSquare_noTrackVertex
//   noTrackVertex_xz  -> X.ladhod.goodhit_chiSquare_noTrackVertex_xz
//   1D_xz_GEM0/1/both -> X.ladhod.goodhit_chiSquare_1D_xz_GEM0 / _GEM1 / _GEMboth
//   1D_y_GEM0/1/both  -> X.ladhod.goodhit_chiSquare_1D_y_GEM0  / _GEM1 / _GEMboth
//   1D_GEM0/1/both    -> X.ladhod.goodhit_chiSquare_1D_GEM0    / _GEM1 / _GEMboth
// Variants absent from the input are skipped automatically. The 1D variants use
// -1 as the "no track" sentinel, so a hit counts as tracked when its chiSquare is
// in [0, CHI_CUT_1D); the 2D variants (large sentinel) use chiSquare < CHI_CUT_2D.
//
// Usage:
//   root -l -b -q 'lad_tracking_eff.C("input.dat","out.root")'
//   root -l -b -q 'lad_tracking_eff.C("input.dat","out.root",8)'  // 8 MT threads
// -------------------------------------------------------------------------

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <RVersion.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TDirectory.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMath.h>
#include <TNamed.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TVirtualPad.h>
#include <cmath>

#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 30, 0)
#if __has_include(<ROOT/RDFHelpers.hxx>)
#include <ROOT/RDFHelpers.hxx>
#define LAD_HAS_RDF_PROGRESSBAR 1
#elif __has_include(<ROOT/RDF/RDFHelpers.hxx>)
#include <ROOT/RDF/RDFHelpers.hxx>
#define LAD_HAS_RDF_PROGRESSBAR 1
#endif
#endif

#include <array>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// =====================================================================
// Constants
// =====================================================================
const int NBINS_TCORR = 900;
const double XMIN_TCORR = -150., XMAX_TCORR = 300.;

const int N_PLANES = 5, N_PADDLES = 11, N_SPECS = 2;
// N_TRACKS is the compile-time capacity for the per-variant histogram arrays.
// Runtime count (ntracks) is <= this. 4 legacy + 9 1D variants = 13.
const int N_TRACKS = 13;
const int PROTON_REBIN = 15; // rebin before sideband-subtracted efficiency
// Two-sided sidebands (tof-L/c ns) for the _c_proton_tof track/total ratio.
const double SB_LO1 = -150., SB_HI1 = -100., SB_LO2 = 125., SB_HI2 = 175.;

// "Has track" chiSquare window per family. 2D variants use a large "no track"
// sentinel, so any chiSquare < CHI_CUT_2D is a track. 1D variants use -1 as the
// sentinel, so a track is chiSquare in [0, CHI_CUT_1D).
const double CHI_CUT_2D = 100.0;
const double CHI_CUT_1D = 100.0;

// Chi-square cut variants. The stored chiSquare branches are RAW chi-square
// (sum of weighted squared residuals; see THcLADKine::FitTrack / FitProj), NOT
// chi-square per DOF -- so a fixed cut is comparatively looser for fits with
// more points (e.g. a 1D GEMboth fit over 3-4 points, or a combined xz+y fit)
// than for a GEM0/GEM1 fit over 2-3 points. We emit one output folder per cut
// (chi2cut_<value>) so the sensitivity to the cut can be seen: the nominal cut
// and its half and double. The window is [chi_lo, chi_hi * scale).
const double CHI_CUT_BASE = 100.0; // nominal cut (== CHI_CUT_2D / CHI_CUT_1D)
const int N_CUTS = 3;
const std::array<double, N_CUTS> CHI_CUT_SCALES = {0.5, 1.0, 2.0};

// edep vs inter-plane tof (t_back - t_front) proton-ID 2x2 panels
const int NBINS_DT = 150;
const double XMIN_DT = -10., XMAX_DT = 30.; // t_back - t_front (ns)
const int NBINS_E2 = 200;
const double EMIN_F = 0., EMAX_F = 400.; // front-plane edep (MeV)
const double EMIN_B = 0., EMAX_B = 150.; // back-plane  edep (MeV)

const double hodo_radii[N_PLANES] = {615., 655.6, 523., 563.6, 615.}; // cm
const char *const plane_names[N_PLANES] = {"000", "001", "100", "101", "200"};
const std::array<char, N_SPECS> specs = {'P', 'H'};

const char *DEFAULT_DAT_FILE = "../files/run-lists/all_C3_runlist_SHMS_13p5.dat";
const char *DEFAULT_OUT_FILE = "files/tracking_eff/tracking_eff_C3_SHMS_13p5_v2_PH.root";

// Sideband-subtract a histogram: mean bin content in [xmin, xmin+sideband_ns]
// is subtracted as a flat background from all bins.
TH1D *flat_bgsub(const TH1D *h, double sideband_ns = 50.) {
  TH1D *out = (TH1D *)h->Clone((std::string(h->GetName()) + "_sb").c_str());
  out->SetTitle((std::string(h->GetTitle()) + " (sb-sub)").c_str());
  double bg_end = h->GetXaxis()->GetXmin() + sideband_ns;
  int n_sb = 0;
  double sum_sb = 0.;
  for (int b = 1; b <= h->GetNbinsX(); ++b) {
    if (h->GetBinCenter(b) >= bg_end)
      break;
    sum_sb += h->GetBinContent(b);
    ++n_sb;
  }
  double bg = (n_sb > 0) ? sum_sb / n_sb : 0.;
  for (int b = 1; b <= h->GetNbinsX(); ++b)
    out->SetBinContent(b, h->GetBinContent(b) - bg);
  return out;
}

// Two-sided sideband subtraction: flat background = mean bin content over the
// bins whose center falls in [lo1,hi1) or [lo2,hi2), subtracted from all bins.
TH1D *flat_bgsub2(const TH1D *h, double lo1, double hi1, double lo2, double hi2) {
  TH1D *out = (TH1D *)h->Clone((std::string(h->GetName()) + "_sb2").c_str());
  out->SetTitle((std::string(h->GetTitle()) + " (2-sb-sub)").c_str());
  int n = 0;
  double sum = 0.;
  for (int b = 1; b <= h->GetNbinsX(); ++b) {
    double x = h->GetBinCenter(b);
    if ((x >= lo1 && x < hi1) || (x >= lo2 && x < hi2)) {
      sum += h->GetBinContent(b);
      ++n;
    }
  }
  double bg = (n > 0) ? sum / n : 0.;
  for (int b = 1; b <= h->GetNbinsX(); ++b)
    out->SetBinContent(b, h->GetBinContent(b) - bg);
  return out;
}

// A tracking variant: 'dir' is its output sub-directory; 'tsuf' is the suffix
// appended to "goodhit_chiSquare"; a hit is tracked when chiSquare is in
// [chi_lo, chi_hi).
struct Track {
  std::string dir, tsuf;
  double chi_lo, chi_hi;
};

// cache_file: optional path to a histogram cache. When non-empty, the first run
// fills the RDataFrame histograms and writes them here; subsequent runs whose
// configuration signature matches load the histograms from it and SKIP the
// event loop entirely (so edits to the plotting/fitting in section 6 re-run in
// seconds). The signature invalidates the cache automatically if the binning,
// chi-square cut scales, tracking-variant list, or the run list changes.
void lad_tracking_eff(const char *dat_file = DEFAULT_DAT_FILE, const char *out_file = DEFAULT_OUT_FILE,
                      int nthreads = 100, const char *cache_file = "") {

  const auto t_start = std::chrono::steady_clock::now();
  gROOT->SetBatch(kTRUE);
  TH1::AddDirectory(kFALSE);
  if (nthreads > 0) {
    ROOT::EnableImplicitMT(nthreads);
    std::cout << "[lad_tracking_eff] implicit MT: " << nthreads << " threads\n";
  } else {
    ROOT::EnableImplicitMT();
    std::cout << "[lad_tracking_eff] implicit MT: all cores\n";
  }

  // ---------------------------------------------------------------
  // 1. TChain
  // ---------------------------------------------------------------
  TChain chain("T");
  std::string datlist; // accumulated kept run-file paths (folded into the cache signature)
  {
    std::ifstream fin(dat_file);
    if (!fin.is_open()) {
      std::cerr << "cannot open " << dat_file << "\n";
      return;
    }
    std::string ln;
    while (std::getline(fin, ln)) {
      size_t a = ln.find_first_not_of(" \t\r\n");
      if (a == std::string::npos)
        continue;
      std::string p = ln.substr(a, ln.find_last_not_of(" \t\r\n") - a + 1);
      if (p.empty() || p[0] == '#')
        continue;
      chain.Add(p.c_str());
      datlist += p + "\n";
    }
  }
  std::cout << "[lad_tracking_eff] entries: " << chain.GetEntries() << "\n";
  if (!chain.GetEntries()) {
    std::cerr << "empty chain\n";
    return;
  }

  // ---------------------------------------------------------------
  // 1b. Detect which tracking variants are present (need both spectrometers so
  //     every alias resolves). Standard tracking is always assumed present.
  // ---------------------------------------------------------------
  chain.LoadTree(0);
  auto has_branch = [&chain](const std::string &n) { return chain.GetBranch(n.c_str()) != nullptr; };
  std::vector<Track> tracks;
  tracks.push_back({"standard", "", -1e30, CHI_CUT_2D});
  {
    const std::array<Track, 12> variants = {{
        {"xz", "_xz", -1e30, CHI_CUT_2D},
        {"noTrackVertex", "_noTrackVertex", -1e30, CHI_CUT_2D},
        {"noTrackVertex_xz", "_noTrackVertex_xz", -1e30, CHI_CUT_2D},
        {"1D_xz_GEM0", "_1D_xz_GEM0", 0., CHI_CUT_1D},
        {"1D_xz_GEM1", "_1D_xz_GEM1", 0., CHI_CUT_1D},
        {"1D_xz_GEMboth", "_1D_xz_GEMboth", 0., CHI_CUT_1D},
        {"1D_y_GEM0", "_1D_y_GEM0", 0., CHI_CUT_1D},
        {"1D_y_GEM1", "_1D_y_GEM1", 0., CHI_CUT_1D},
        {"1D_y_GEMboth", "_1D_y_GEMboth", 0., CHI_CUT_1D},
        {"1D_GEM0", "_1D_GEM0", 0., CHI_CUT_1D},
        {"1D_GEM1", "_1D_GEM1", 0., CHI_CUT_1D},
        {"1D_GEMboth", "_1D_GEMboth", 0., CHI_CUT_1D},
    }};
    for (const auto &v : variants) {
      bool ok = true;
      for (int is = 0; is < N_SPECS; ++is)
        if (!has_branch(std::string(1, specs[is]) + ".ladhod.goodhit_chiSquare" + v.tsuf))
          ok = false;
      if (ok)
        tracks.push_back(v);
      else
        std::cout << "[lad_tracking_eff] tracking variant '" << v.dir << "' absent in data; skipping\n";
    }
  }
  const int ntracks = (int)tracks.size();
  std::cout << "[lad_tracking_eff] tracking variants (" << ntracks << "):";
  for (const auto &t : tracks)
    std::cout << " " << t.dir;
  std::cout << "\n";

  // ---------------------------------------------------------------
  // 1c. Histogram cache decision. Build a configuration signature; if a cache
  //     file exists with a matching signature we load histograms from it and
  //     skip the (expensive) event loop. Bump CACHE_VERSION whenever the set of
  //     booked histograms changes so old caches are rejected.
  // ---------------------------------------------------------------
  const char *CACHE_VERSION = "v1";
  std::string sig = std::string("lad_tracking_eff;") + CACHE_VERSION + ";";
  sig += "tof=" + std::to_string(NBINS_TCORR) + "," + std::to_string(XMIN_TCORR) + "," + std::to_string(XMAX_TCORR) +
         ";dt=" + std::to_string(NBINS_DT) + "," + std::to_string(XMIN_DT) + "," + std::to_string(XMAX_DT) +
         ";e2=" + std::to_string(NBINS_E2) + "," + std::to_string(EMIN_F) + "," + std::to_string(EMAX_F) + "," +
         std::to_string(EMIN_B) + "," + std::to_string(EMAX_B) + ";cuts=";
  for (int ic = 0; ic < N_CUTS; ++ic)
    sig += std::to_string(CHI_CUT_SCALES[ic]) + ",";
  sig += ";vars=";
  for (const auto &t : tracks)
    sig += t.dir + "|" + std::to_string(t.chi_lo) + "|" + std::to_string(t.chi_hi) + ",";
  sig += ";runlist=" + std::to_string((unsigned long long)std::hash<std::string>{}(datlist));

  const bool cache_on = (cache_file && cache_file[0] != '\0');
  bool load = false;
  TFile *fcache = nullptr;
  if (cache_on && !gSystem->AccessPathName(cache_file)) { // AccessPathName == false means the file exists
    fcache = TFile::Open(cache_file, "READ");
    if (fcache && !fcache->IsZombie()) {
      TNamed *s = dynamic_cast<TNamed *>(fcache->Get("signature"));
      if (s && sig == std::string(s->GetTitle()))
        load = true;
    }
    if (!load && fcache) {
      fcache->Close();
      delete fcache;
      fcache = nullptr;
    }
  }
  if (cache_on)
    std::cout << "[lad_tracking_eff] cache " << (load ? "HIT -> loading histograms, skipping event loop" : "MISS")
              << ": " << cache_file << "\n";

  // ---------------------------------------------------------------
  // 2. RDataFrame + aliases (only the columns the proton-ID plots need)
  // ---------------------------------------------------------------
  ROOT::RDataFrame rdf(chain);
  ROOT::RDF::RNode df = rdf;

  // Booking helpers shared by the fill and load paths. In fill mode BK1/BK2 book
  // an RDataFrame action and remember (histogram slot, result) so the slot can
  // be resolved to a raw pointer after the single event loop; in load mode they
  // read the histogram straight from the cache file by name.
  struct HBind1 {
    TH1D **slot;
    ROOT::RDF::RResultPtr<TH1D> res;
  };
  struct HBind2 {
    TH2D **slot;
    ROOT::RDF::RResultPtr<TH2D> res;
  };
  std::vector<HBind1> bind1;
  std::vector<HBind2> bind2;
  auto BK1 = [&](TH1D *&slot, const std::string &nm, const std::string &tt) {
    if (load) {
      slot = dynamic_cast<TH1D *>(fcache->Get(nm.c_str()));
      return;
    }
    bind1.push_back({&slot, df.Histo1D({nm.c_str(), tt.c_str(), NBINS_TCORR, XMIN_TCORR, XMAX_TCORR}, nm)});
  };
  auto BK2 = [&](TH2D *&slot, const std::string &xcol, const std::string &ycol, const std::string &nm,
                 const std::string &tt, double ymin, double ymax) {
    if (load) {
      slot = dynamic_cast<TH2D *>(fcache->Get(nm.c_str()));
      return;
    }
    bind2.push_back(
        {&slot, df.Histo2D({nm.c_str(), tt.c_str(), NBINS_DT, XMIN_DT, XMAX_DT, NBINS_E2, ymin, ymax}, xcol, ycol)});
  };

  if (!load) { // sections 2-3 (aliases + column definitions) are only needed to fill
  ROOT::RDF::RResultPtr<ULong64_t> progress_count;
#ifdef LAD_HAS_RDF_PROGRESSBAR
  ROOT::RDF::Experimental::AddProgressBar(df);
#else
  {
    const ULong64_t total = chain.GetEntries();
    const ULong64_t step = std::max<ULong64_t>(1ULL, total / 200ULL);
    auto t0 = std::make_shared<std::chrono::steady_clock::time_point>(std::chrono::steady_clock::now());
    progress_count = df.Count();
    progress_count.OnPartialResult(step, [total, t0](ULong64_t n) {
      double f = total ? double(n) / double(total) : 0.;
      double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - *t0).count();
      std::fprintf(stderr, "\r[lad_tracking_eff] %6.2f%% (%llu/%llu) %.1fs eta %.1fs   ", 100. * f,
                   (unsigned long long)n, (unsigned long long)total, s, f > 0 ? s * (1. / f - 1.) : 0.);
      std::fflush(stderr);
    });
  }
#endif

  for (int is = 0; is < N_SPECS; ++is) {
    const std::string sp(1, specs[is]);
    const std::string pfx = sp + ".ladhod.goodhit_";
    // proton-ID uses side-1 hits only
    df = df.Alias(sp + "_plane_1", pfx + "plane_1");
    df = df.Alias(sp + "_paddle_1", pfx + "paddle_1");
    df = df.Alias(sp + "_tof_1", pfx + "hit_tof_1");
    df = df.Alias(sp + "_ypos_1", pfx + "hit_ypos_1");
    df = df.Alias(sp + "_isProton_1", pfx + "isProton_1");
    // front/back plane times and calibrated edep (for the edep-vs-tof panels)
    df = df.Alias(sp + "_hittime_0", pfx + "hittime_0");
    df = df.Alias(sp + "_hittime_1", pfx + "hittime_1");
    df = df.Alias(sp + "_edepMeV_0", pfx + "hitedep_MeV_0");
    df = df.Alias(sp + "_edepMeV_1", pfx + "hitedep_MeV_1");
    for (const auto &tk : tracks)
      df = df.Alias(sp + "_chiSquare" + tk.tsuf, pfx + "chiSquare" + tk.tsuf);
  }

  // ---------------------------------------------------------------
  // 3. Proton-tagged corrected-tof column definitions (planes 001/101 only).
  //    proton-only (no track) is tracking-independent; proton+track is built
  //    once per variant with that variant's [chi_lo, chi_hi) track window.
  // ---------------------------------------------------------------
  using RVd = ROOT::VecOps::RVec<double>;
  for (int is = 0; is < N_SPECS; ++is) {
    const std::string sp(1, specs[is]);

    auto mk_proton = [&](const std::string &col, bool req_track, const std::string &chicol, double clo, double chi_hi) {
      df = df.Define(col,
                     [req_track, clo, chi_hi](const RVd &pl1, const RVd &pd1, const RVd &y1, const RVd &t1,
                                              const RVd &ip1, const RVd &chi) {
                       RVd r;
                       for (size_t i = 0; i < pl1.size(); ++i) {
                         if (ip1[i] != 1.)
                           continue;
                         if (req_track && !(chi[i] >= clo && chi[i] < chi_hi))
                           continue;
                         int pi = (int)std::round(pl1[i]);
                         if (pi != 1 && pi != 3)
                           continue; // planes 001 and 101 only
                         double R = hodo_radii[pi];
                         double dx = 22. * (pd1[i] - 6.);
                         double p2d = std::sqrt(y1[i] * y1[i] + dx * dx);
                         r.push_back(t1[i] - std::sqrt(p2d * p2d + R * R) / 100. / 0.3);
                       }
                       return r;
                     },
                     {sp + "_plane_1", sp + "_paddle_1", sp + "_ypos_1", sp + "_tof_1", sp + "_isProton_1", chicol});
    };
    mk_proton(sp + "_tof_corr_proton", false, sp + "_chiSquare", 0., 0.);
    for (int ic = 0; ic < N_CUTS; ++ic)
      for (const auto &tk : tracks)
        mk_proton(sp + "_tof_corr_proton_track" + tk.tsuf + "_cut" + std::to_string(ic), true,
                  sp + "_chiSquare" + tk.tsuf, tk.chi_lo, tk.chi_hi * CHI_CUT_SCALES[ic]);

    for (int pp = 0; pp < 2; ++pp) {
      const int pi_c = (pp == 0) ? 1 : 3;
      const double R_c = hodo_radii[pi_c];
      const bool excl_c = (pi_c == 3); // plane 101 excludes paddles 1,9 in the sum

      for (int paddle = 0; paddle < N_PADDLES; ++paddle) {
        const double pv_c = static_cast<double>(paddle);
        const std::string ps = std::to_string(paddle);
        auto mk_p = [&](const std::string &col, bool req_track, const std::string &chicol, double clo, double chi_hi) {
          df = df.Define(col,
                         [pi_c, R_c, pv_c, req_track, clo, chi_hi](const RVd &pl1, const RVd &pd1, const RVd &y1,
                                                                   const RVd &t1, const RVd &ip1, const RVd &chi) {
                           RVd r;
                           for (size_t i = 0; i < pl1.size(); ++i) {
                             if (ip1[i] != 1.)
                               continue;
                             if ((int)std::round(pl1[i]) != pi_c)
                               continue;
                             if (pd1[i] != pv_c)
                               continue;
                             if (req_track && !(chi[i] >= clo && chi[i] < chi_hi))
                               continue;
                             double dx = 22. * (pd1[i] - 6.);
                             double p2d = std::sqrt(y1[i] * y1[i] + dx * dx);
                             r.push_back(t1[i] - std::sqrt(p2d * p2d + R_c * R_c) / 100. / 0.3);
                           }
                           return r;
                         },
                         {sp + "_plane_1", sp + "_paddle_1", sp + "_ypos_1", sp + "_tof_1", sp + "_isProton_1", chicol});
        };
        mk_p(sp + "_tof_corr_proton_p" + std::to_string(pi_c) + "_b" + ps, false, sp + "_chiSquare", 0., 0.);
        for (int ic = 0; ic < N_CUTS; ++ic)
          for (const auto &tk : tracks)
            mk_p(sp + "_tof_corr_proton_track_p" + std::to_string(pi_c) + "_b" + ps + tk.tsuf + "_cut" +
                     std::to_string(ic),
                 true, sp + "_chiSquare" + tk.tsuf, tk.chi_lo, tk.chi_hi * CHI_CUT_SCALES[ic]);
      }

      auto mk_sum = [&](const std::string &col, bool req_track, const std::string &chicol, double clo, double chi_hi) {
        df = df.Define(col,
                       [pi_c, R_c, excl_c, req_track, clo, chi_hi](const RVd &pl1, const RVd &pd1, const RVd &y1,
                                                                   const RVd &t1, const RVd &ip1, const RVd &chi) {
                         RVd r;
                         for (size_t i = 0; i < pl1.size(); ++i) {
                           if (ip1[i] != 1.)
                             continue;
                           if ((int)std::round(pl1[i]) != pi_c)
                             continue;
                           if (excl_c && (pd1[i] == 1. || pd1[i] == 9.))
                             continue;
                           if (req_track && !(chi[i] >= clo && chi[i] < chi_hi))
                             continue;
                           double dx = 22. * (pd1[i] - 6.);
                           double p2d = std::sqrt(y1[i] * y1[i] + dx * dx);
                           r.push_back(t1[i] - std::sqrt(p2d * p2d + R_c * R_c) / 100. / 0.3);
                         }
                         return r;
                       },
                       {sp + "_plane_1", sp + "_paddle_1", sp + "_ypos_1", sp + "_tof_1", sp + "_isProton_1", chicol});
      };
      mk_sum(sp + "_tof_corr_proton_p" + std::to_string(pi_c) + "_sum", false, sp + "_chiSquare", 0., 0.);
      for (int ic = 0; ic < N_CUTS; ++ic)
        for (const auto &tk : tracks)
          mk_sum(sp + "_tof_corr_proton_track_p" + std::to_string(pi_c) + "_sum" + tk.tsuf + "_cut" +
                     std::to_string(ic),
                 true, sp + "_chiSquare" + tk.tsuf, tk.chi_lo, tk.chi_hi * CHI_CUT_SCALES[ic]);
    }

    // edep-vs-(inter-plane tof) columns for the proton-ID 2x2 panels. dt =
    // t_back - t_front; front/back edep in MeV. Selection = proton cut
    // (isProton_1==1), optionally AND the variant's [chi_lo, chi_hi) track window.
    auto mk_pid = [&](const std::string &tag, bool req_track, const std::string &chicol, double clo, double chi_hi) {
      auto keep = [req_track, clo, chi_hi](double ip, double chi) {
        return ip == 1. && (!req_track || (chi >= clo && chi < chi_hi));
      };
      df = df.Define(sp + "_pid_dt_" + tag,
                     [keep](const RVd &t0, const RVd &t1, const RVd &ip1, const RVd &chi) {
                       RVd r;
                       for (size_t i = 0; i < ip1.size(); ++i)
                         if (keep(ip1[i], chi[i]))
                           r.push_back(t1[i] - t0[i]);
                       return r;
                     },
                     {sp + "_hittime_0", sp + "_hittime_1", sp + "_isProton_1", chicol});
      df = df.Define(sp + "_pid_ef_" + tag,
                     [keep](const RVd &e0, const RVd &ip1, const RVd &chi) {
                       RVd r;
                       for (size_t i = 0; i < ip1.size(); ++i)
                         if (keep(ip1[i], chi[i]))
                           r.push_back(e0[i]);
                       return r;
                     },
                     {sp + "_edepMeV_0", sp + "_isProton_1", chicol});
      df = df.Define(sp + "_pid_eb_" + tag,
                     [keep](const RVd &e1, const RVd &ip1, const RVd &chi) {
                       RVd r;
                       for (size_t i = 0; i < ip1.size(); ++i)
                         if (keep(ip1[i], chi[i]))
                           r.push_back(e1[i]);
                       return r;
                     },
                     {sp + "_edepMeV_1", sp + "_isProton_1", chicol});
    };
    mk_pid("pro", false, sp + "_chiSquare", 0., 0.);
    for (int ic = 0; ic < N_CUTS; ++ic)
      for (const auto &tk : tracks)
        mk_pid("trk" + tk.tsuf + "_cut" + std::to_string(ic), true, sp + "_chiSquare" + tk.tsuf, tk.chi_lo,
               tk.chi_hi * CHI_CUT_SCALES[ic]);
  }
  } // end if(!load): RDF alias/define setup

  // ---------------------------------------------------------------
  // 4. Book proton histograms (BK1/BK2 either book RDF actions or, when loading
  //    from cache, fetch the histograms by name). Raw pointers are used so both
  //    paths share one type; fill-path pointers are resolved after the loop.
  // ---------------------------------------------------------------
  std::array<TH1D *, N_SPECS> h_proton_tof{};
  // Track histograms carry a leading [N_CUTS] chi-square-cut index.
  std::array<std::array<std::array<TH1D *, N_TRACKS>, N_CUTS>, N_SPECS> h_proton_track_tof{};
  std::array<std::array<std::vector<TH1D *>, 2>, N_SPECS> h_proton_pad;
  std::array<std::array<std::array<std::array<std::vector<TH1D *>, 2>, N_TRACKS>, N_CUTS>, N_SPECS> h_proton_track_pad;
  std::array<std::array<TH1D *, 2>, N_SPECS> h_proton_psum{};
  std::array<std::array<std::array<std::array<TH1D *, 2>, N_TRACKS>, N_CUTS>, N_SPECS> h_proton_track_psum{};
  // edep-vs-tof 2D panels: proton-only (variant-independent) and proton+track.
  std::array<TH2D *, N_SPECS> h_pid_ef_pro{}, h_pid_eb_pro{};                            // front/back, proton cut
  std::array<std::array<std::array<TH2D *, N_TRACKS>, N_CUTS>, N_SPECS> h_pid_ef_trk{}, h_pid_eb_trk{}; // proton+track

  for (int is = 0; is < N_SPECS; ++is) {
    const std::string sp(1, specs[is]);
    BK1(h_proton_tof[is], sp + "_tof_corr_proton", sp + " tof corr proton;tof-L/c(ns);Counts");
    for (int ic = 0; ic < N_CUTS; ++ic)
      for (int it = 0; it < ntracks; ++it) {
        const std::string &td = tracks[it].dir;
        const std::string &tu = tracks[it].tsuf;
        BK1(h_proton_track_tof[is][ic][it], sp + "_tof_corr_proton_track" + tu + "_cut" + std::to_string(ic),
            sp + " tof corr proton+track [" + td + "];tof-L/c(ns);Counts");
      }
    for (int pp = 0; pp < 2; ++pp) {
      const int pi = (pp == 0) ? 1 : 3;
      const std::string pn = plane_names[pi];
      h_proton_pad[is][pp].resize(N_PADDLES);
      for (int ic = 0; ic < N_CUTS; ++ic)
        for (int it = 0; it < ntracks; ++it)
          h_proton_track_pad[is][ic][it][pp].resize(N_PADDLES);
      for (int pa = 0; pa < N_PADDLES; ++pa) {
        const std::string ps = std::to_string(pa);
        const std::string cp = sp + "_tof_corr_proton_p" + std::to_string(pi) + "_b" + ps;
        BK1(h_proton_pad[is][pp][pa], cp, sp + " tof corr proton " + pn + " pd" + ps + ";tof-L/c(ns);Counts");
        for (int ic = 0; ic < N_CUTS; ++ic)
          for (int it = 0; it < ntracks; ++it) {
            const std::string &td = tracks[it].dir;
            const std::string &tu = tracks[it].tsuf;
            const std::string ct =
                sp + "_tof_corr_proton_track_p" + std::to_string(pi) + "_b" + ps + tu + "_cut" + std::to_string(ic);
            BK1(h_proton_track_pad[is][ic][it][pp][pa], ct,
                sp + " tof corr proton+track [" + td + "] " + pn + " pd" + ps + ";tof-L/c(ns);Counts");
          }
      }
      const std::string cs2 = sp + "_tof_corr_proton_p" + std::to_string(pi) + "_sum";
      BK1(h_proton_psum[is][pp], cs2, sp + " tof corr proton " + pn + " sum;tof-L/c(ns);Counts");
      for (int ic = 0; ic < N_CUTS; ++ic)
        for (int it = 0; it < ntracks; ++it) {
          const std::string &td = tracks[it].dir;
          const std::string &tu = tracks[it].tsuf;
          const std::string ts2 =
              sp + "_tof_corr_proton_track_p" + std::to_string(pi) + "_sum" + tu + "_cut" + std::to_string(ic);
          BK1(h_proton_track_psum[is][ic][it][pp], ts2,
              sp + " tof corr proton+track [" + td + "] " + pn + " sum;tof-L/c(ns);Counts");
        }
    }

    // edep-vs-tof 2D panels
    BK2(h_pid_ef_pro[is], sp + "_pid_dt_pro", sp + "_pid_ef_pro", sp + "_pid_front_pro",
        sp + " front edep vs tof (proton);t_{back}-t_{front} (ns);front edep (MeV)", EMIN_F, EMAX_F);
    BK2(h_pid_eb_pro[is], sp + "_pid_dt_pro", sp + "_pid_eb_pro", sp + "_pid_back_pro",
        sp + " back edep vs tof (proton);t_{back}-t_{front} (ns);back edep (MeV)", EMIN_B, EMAX_B);
    for (int ic = 0; ic < N_CUTS; ++ic)
      for (int it = 0; it < ntracks; ++it) {
        const std::string &td = tracks[it].dir;
        const std::string &tu = tracks[it].tsuf;
        const std::string cc = "_cut" + std::to_string(ic);
        BK2(h_pid_ef_trk[is][ic][it], sp + "_pid_dt_trk" + tu + cc, sp + "_pid_ef_trk" + tu + cc,
            sp + "_pid_front_trk" + tu + cc,
            sp + " front edep vs tof (proton+track [" + td + "]);t_{back}-t_{front} (ns);front edep (MeV)", EMIN_F,
            EMAX_F);
        BK2(h_pid_eb_trk[is][ic][it], sp + "_pid_dt_trk" + tu + cc, sp + "_pid_eb_trk" + tu + cc,
            sp + "_pid_back_trk" + tu + cc,
            sp + " back edep vs tof (proton+track [" + td + "]);t_{back}-t_{front} (ns);back edep (MeV)", EMIN_B,
            EMAX_B);
      }
  }

  // ---------------------------------------------------------------
  // 5. Trigger the event loop (fill mode) and resolve histogram pointers, or
  //    just report that histograms came from the cache (load mode).
  // ---------------------------------------------------------------
  if (!load) {
    std::cout << "[lad_tracking_eff] Running event loop...\n";
    // The first GetPtr() runs the single event loop that fills every booked action.
    for (auto &b : bind1)
      *b.slot = b.res.GetPtr();
    for (auto &b : bind2)
      *b.slot = b.res.GetPtr();
#ifndef LAD_HAS_RDF_PROGRESSBAR
    std::fprintf(stderr, "\n");
#endif
    if (cache_on) {
      TFile *wf = TFile::Open(cache_file, "RECREATE");
      if (wf && !wf->IsZombie()) {
        wf->cd();
        TNamed signm("signature", sig.c_str());
        signm.Write();
        for (auto &b : bind1)
          if (*b.slot)
            (*b.slot)->Write((*b.slot)->GetName(), TObject::kOverwrite);
        for (auto &b : bind2)
          if (*b.slot)
            (*b.slot)->Write((*b.slot)->GetName(), TObject::kOverwrite);
        wf->Close();
        delete wf;
        std::cout << "[lad_tracking_eff] wrote histogram cache: " << cache_file << "\n";
      } else {
        std::cerr << "[lad_tracking_eff] warning: could not write cache " << cache_file << "\n";
      }
    }
  } else {
    std::cout << "[lad_tracking_eff] Histograms loaded from cache; event loop skipped.\n";
  }

  // ---------------------------------------------------------------
  // 6. Write proton_id plots (one sub-directory per tracking variant)
  // ---------------------------------------------------------------
  TFile fout(out_file, "RECREATE");
  if (fout.IsZombie()) {
    std::cerr << "cannot open output\n";
    return;
  }
  auto wc = [](TCanvas *c) {
    c->Write();
    delete c;
  };
  auto mkpath = [](TDirectory *base, const std::string &path) -> TDirectory * {
    TDirectory *dcur = base;
    std::string seg;
    std::istringstream ss(path);
    while (std::getline(ss, seg, '/')) {
      if (seg.empty())
        continue;
      TDirectory *nd = dcur->GetDirectory(seg.c_str());
      if (!nd)
        nd = dcur->mkdir(seg.c_str());
      dcur = nd;
    }
    return dcur ? dcur : base;
  };

  // Mean and stdev of a histogram's bin contents over one or more x-ranges
  // [lo,hi). Truly-empty bins (content == 0 and error == 0) are skipped so a
  // region that extends past the data does not drag the average toward 0. The
  // stdev (spread of the region's bin values) is returned as the error.
  auto region_stats = [](const TH1D *h, const std::vector<std::array<double, 2>> &ranges, double &mean, double &err) {
    double s = 0., s2 = 0.;
    int n = 0;
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
      double x = h->GetBinCenter(b);
      bool in = false;
      for (const auto &r : ranges)
        if (x >= r[0] && x < r[1]) {
          in = true;
          break;
        }
      if (!in)
        continue;
      double v = h->GetBinContent(b), e = h->GetBinError(b);
      if (v == 0.0 && e == 0.0)
        continue; // empty bin, outside the filled data
      s += v;
      s2 += v * v;
      ++n;
    }
    if (n == 0) {
      mean = 0.;
      err = 0.;
      return;
    }
    mean = s / n;
    double var = (n > 1) ? (s2 - n * mean * mean) / (n - 1) : 0.0;
    err = (var > 0.) ? std::sqrt(var) : 0.0;
  };

  // Event-weighted version of region_stats: each bin's ratio value is weighted
  // by the number of events in that bin (numerator + denominator counts,
  // wnum + wden). Returns the weighted mean and the weighted stdev of the
  // values (as the error). Truly-empty ratio bins and zero-weight bins skipped.
  auto region_wstats = [](const TH1D *h, const TH1D *wnum, const TH1D *wden,
                          const std::vector<std::array<double, 2>> &ranges, double &mean, double &err) {
    double sw = 0., swv = 0., swv2 = 0.;
    int n = 0;
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
      double x = h->GetBinCenter(b);
      bool in = false;
      for (const auto &r : ranges)
        if (x >= r[0] && x < r[1]) {
          in = true;
          break;
        }
      if (!in)
        continue;
      double v = h->GetBinContent(b), e = h->GetBinError(b);
      if (v == 0.0 && e == 0.0)
        continue; // empty ratio bin, outside the filled data
      double w = wnum->GetBinContent(b) + wden->GetBinContent(b);
      if (w <= 0.)
        continue;
      sw += w;
      swv += w * v;
      swv2 += w * v * v;
      ++n;
    }
    if (sw <= 0.) {
      mean = 0.;
      err = 0.;
      return;
    }
    mean = swv / sw;
    double var = (n > 1) ? (swv2 / sw - mean * mean) * n / (n - 1) : 0.0;
    err = (var > 0.) ? std::sqrt(var) : 0.0;
  };

  for (int is = 0; is < N_SPECS; ++is) {
    const std::string sp(1, specs[is]);
    TDirectory *sdir = fout.mkdir(sp.c_str());

    TDirectory *d = sdir->mkdir("proton_id");
    // Per-1D-variant summary accumulators, one series per chi-square cut. Filled
    // inside the cut loop and overlaid on a single plot per metric afterwards.
    std::vector<std::string> sum_names; // variant names (cut-independent)
    std::array<std::vector<double>, N_CUTS> sum_eff, sum_eff_err; // peak-region avg of bg-sub ratio
    std::array<std::vector<double>, N_CUTS> sum_bkg, sum_bkg_err; // sideband avg of raw track/total
    // Extra summary metrics from the pad-2 (proton+track) fit. All background
    // quantities are scaled to the signal's effective width so they are directly
    // comparable to the signal event count (see where they are computed below).
    std::array<std::vector<double>, N_CUTS> sum_sig;   // signal: peak gaussian integral (events)
    std::array<std::vector<double>, N_CUTS> sum_oob;   // out-of-time bg: sideband integral, signal-width norm
    std::array<std::vector<double>, N_CUTS> sum_itb;   // in-time  bg: trapezoid integral, signal-width norm
    std::array<std::vector<double>, N_CUTS> sum_soob;  // signal / out-of-time bg
    std::array<std::vector<double>, N_CUTS> sum_sitb;  // signal / in-time  bg
    std::array<std::vector<double>, N_CUTS> sum_itoob; // in-time bg / out-of-time bg
    std::array<std::vector<double>, N_CUTS> sum_zero;  // shared zero-error series for the metrics above
    // No-track-cut reference from the pad-1 (proton total, top-left) fit: the
    // same signal / out-of-time-bg / ratio, but with NO chi-square or track cut.
    // Overlaid on the signal, out-of-time-bg, and signal/oob summary plots.
    std::array<std::vector<double>, N_CUTS> sum_sig_nc;  // no-cut signal (peak integral)
    std::array<std::vector<double>, N_CUTS> sum_oob_nc;  // no-cut out-of-time bg
    std::array<std::vector<double>, N_CUTS> sum_soob_nc; // no-cut signal / out-of-time bg
    // One folder per chi-square cut: proton_id/chi2cut_<value>/<variant>/...
    for (int ic = 0; ic < N_CUTS; ++ic) {
      const double cutval = CHI_CUT_SCALES[ic] * CHI_CUT_BASE;
      const std::string cc = "_cut" + std::to_string(ic);
      const std::string cutstr = std::to_string((int)std::lround(cutval));
      TDirectory *dc = d->mkdir(("chi2cut_" + cutstr).c_str());
      for (int it = 0; it < ntracks; ++it) {
        const std::string &tu = tracks[it].tsuf;
        const bool is1D = tracks[it].dir.rfind("1D_", 0) == 0;
        double eff_m = 0., eff_e = 0., bkg_m = 0., bkg_e = 0.;
        TDirectory *td = mkpath(dc, tracks[it].dir);
        td->cd();

      TCanvas *c = new TCanvas((sp + "_c_proton_tof").c_str(),
                               (sp + " proton tof [" + tracks[it].dir + ", chi2<" + cutstr + "]").c_str(), 1400, 1000);
      c->Divide(2, 2);
      // Mark regions with HATCHED bands rather than a solid fill. A hatch
      // pattern (diagonal lines, fill style 33xx) is inherently see-through, so
      // the histogram shows through the gaps regardless of paint order -- this
      // sidesteps JSROOT (the VSCode ROOT File Viewer), which composites the
      // histogram into the frame layer and paints TBoxes on top no matter the
      // primitive order, and which renders alpha-transparent solid fills as
      // outline-only. We still draw axes first, then the boxes, then the
      // histogram content on top, so native ROOT layers it cleanly too. Each
      // band carries its own {x1, x2, color, fillstyle}. TBoxes are leaked so
      // they survive until the canvas is written.
      auto draw_shaded = [](TH1D *h, const std::vector<std::array<double, 4>> &bands) {
        h->DrawCopy(); // full draw: frame + title + content (first = main painter)
        gPad->Update();
        double y1 = gPad->GetUymin(), y2 = gPad->GetUymax();
        for (const auto &bnd : bands) {
          TBox *b = new TBox(bnd[0], y1, bnd[1], y2);
          b->SetFillColor((int)bnd[2]);
          b->SetFillStyle((int)bnd[3]);
          b->SetLineColor((int)bnd[2]);
          b->SetLineWidth(1);
          b->Draw();
        }
        h->DrawCopy("same"); // histogram content back on top of the boxes
        gPad->RedrawAxis();
      };
      const int kSB = kGray + 2;  // sideband hatch color
      const int kPk = kRed;       // peak hatch color
      const int fsSB = 3004;      // sideband fill: diagonal hatch  ///
      const int fsPk = 3005;      // peak fill:     anti-diagonal hatch  \\\
      // Fits for pads 1 & 2, over (-150,150). Pad 2 (proton+track) = flat line +
      // a fixed-corner trapezoid (corners at x = -75,-25,50,125) + a gaussian;
      // pad 1 (proton total) = flat line + gaussian whose width is FIXED to the
      // pad-2 gaussian width ("same width"). Both display the gaussian integral
      // expressed as a number of events (curve area / bin width).
      auto trapgaus = [](double *xx, double *p) -> double {
        const double x = xx[0];
        const double flat = p[0];
        const double amp = p[1] - p[0]; // trapezoid plateau height above the flat line
        double trap;
        if (x < -75.)
          trap = 0.;
        else if (x < -25.)
          trap = amp * (x + 75.) / 50.; // rise  -75 -> -25
        else if (x < 50.)
          trap = amp; // plateau -25 -> 50
        else if (x < 125.)
          trap = amp * (125. - x) / 75.; // fall   50 -> 125
        else
          trap = 0.;
        const double g = p[2] * std::exp(-0.5 * std::pow((x - p[3]) / p[4], 2));
        return flat + trap + g;
      };
      auto flatgaus = [](double *xx, double *p) -> double {
        const double x = xx[0];
        return p[0] + p[1] * std::exp(-0.5 * std::pow((x - p[2]) / p[3], 2));
      };
      // Draw the number of events under a fitted gaussian (curve area / bin
      // width) in the current pad.
      auto draw_gaus_integral = [](TH1D *h, double amp, double sigma) {
        const double binw = h->GetXaxis()->GetBinWidth(1);
        const double area = amp * std::fabs(sigma) * std::sqrt(TMath::TwoPi());
        const double nev = (binw > 0.) ? area / binw : 0.;
        TLatex *tx = new TLatex();
        tx->SetNDC();
        tx->SetTextSize(0.035);
        tx->SetTextColor(kBlue + 2);
        tx->DrawLatex(0.14, 0.84, Form("Gaus integral: %.0f", nev));
      };
      // Subtract a fit-derived background from a (possibly rebinned) histogram.
      // fbg is the non-gaussian part of the fit (counts per ORIGINAL bin), so
      // the background counts in a wider rebinned bin = integral of fbg over the
      // bin / original bin width.
      auto subtract_fit_bg = [](TH1D *h, TF1 *fbg, double orig_binw) {
        for (int b = 1; b <= h->GetNbinsX(); ++b) {
          double xlo = h->GetXaxis()->GetBinLowEdge(b), xhi = h->GetXaxis()->GetBinUpEdge(b);
          double bg = (orig_binw > 0.) ? fbg->Integral(xlo, xhi) / orig_binw : 0.;
          h->SetBinContent(b, h->GetBinContent(b) - bg);
        }
      };

      // Fit pad-2 histogram FIRST so its gaussian width can be reused for pad 1.
      TH1D *ht2 = (TH1D *)h_proton_track_tof[is][ic][it]->Clone((sp + "_proton_track_tof_p2" + tu).c_str());
      TF1 *f2 = new TF1((sp + "_fit_trapgaus" + tu + cc).c_str(), trapgaus, -150., 150., 5);
      f2->SetParNames("flat", "trap_top", "gaus_h", "gaus_mean", "gaus_sigma");
      f2->SetParameters(50., 70., 100., 41., 5.);
      f2->FixParameter(3, 41.); // gaussian center fixed at 41 ns
      f2->SetLineColor(kGreen + 2);
      f2->SetNpx(600);
      // N = don't attach the function to the histogram (so deleting the clone
      // won't delete it); 0 = don't draw now. We Draw() the TF1 explicitly as a
      // standalone "same" primitive after the histogram, and leak it so it
      // survives until the canvas is written.
      ht2->Fit(f2, "RQN0");
      const double sig2 = f2->GetParameter(4);

      // pad 1: proton total (full binning), sidebands + peak region shaded.
      // Fit = flat + gaussian, gaussian width fixed to the pad-2 gaussian width.
      c->cd(1);
      TH1D *hp1 = (TH1D *)h_proton_tof[is]->Clone((sp + "_proton_tof_p1" + tu).c_str());
      TF1 *f1 = new TF1((sp + "_fit_flatgaus" + tu + cc).c_str(), flatgaus, -150., 150., 4);
      f1->SetParNames("flat", "gaus_h", "gaus_mean", "gaus_sigma");
      f1->SetParameters(50., 100., 41., sig2);
      f1->FixParameter(2, 41.);  // gaussian center fixed at 41 ns
      f1->FixParameter(3, sig2); // same width as the pad-2 gaussian
      f1->SetLineColor(kGreen + 2);
      f1->SetNpx(600);
      hp1->Fit(f1, "RQN0"); // N = don't attach to hist; 0 = don't draw; drawn explicitly below
      draw_shaded(hp1, {{SB_LO1, SB_HI1, (double)kSB, (double)fsSB},
                        {SB_LO2, SB_HI2, (double)kSB, (double)fsSB},
                        {30., 50., (double)kPk, (double)fsPk}});
      f1->Draw("same"); // fit curve on top, as its own canvas primitive
      draw_gaus_integral(hp1, f1->GetParameter(1), f1->GetParameter(3));
      delete hp1;
      // pad 2: proton+track (full binning), sidebands + peak region shaded.
      c->cd(2);
      draw_shaded(ht2, {{SB_LO1, SB_HI1, (double)kSB, (double)fsSB},
                        {SB_LO2, SB_HI2, (double)kSB, (double)fsSB},
                        {30., 50., (double)kPk, (double)fsPk}});
      f2->Draw("same"); // fit curve on top, as its own canvas primitive
      draw_gaus_integral(ht2, f2->GetParameter(2), f2->GetParameter(4));
      delete ht2;

      // ----- extra summary metrics from the pad-2 (proton+track) fit -----
      // Signal        = full gaussian integral (events under the peak),
      //                 amp * sigma * sqrt(2pi) / binw. Its "effective width"
      //                 (the width of a rectangle of height amp with the same
      //                 area) is w_sig = sigma * sqrt(2pi); the backgrounds are
      //                 rescaled to this width so all three are event counts in
      //                 a signal-width window.
      // Out-of-time bg = proton+track counts in the two sidebands (total width
      //                 100 ns), rescaled to the signal width.
      // In-time    bg = integral of the trapezoid component ABOVE the flat line
      //                 (the non-random in-time excess), from its 2nd corner
      //                 (x=-25) to x=170 (width 195 ns), rescaled to the signal
      //                 width. The flat baseline is the out-of-time level and is
      //                 excluded here.
      TH1D *hpt = h_proton_track_tof[is][ic][it];
      const double p2binw = hpt->GetXaxis()->GetBinWidth(1);
      const double g_amp = f2->GetParameter(2), g_sig = std::fabs(f2->GetParameter(4));
      const double w_sig = g_sig * std::sqrt(TMath::TwoPi());
      const double sig_ev = (p2binw > 0.) ? g_amp * w_sig / p2binw : 0.;
      const double w_sb = (SB_HI1 - SB_LO1) + (SB_HI2 - SB_LO2);
      const double oob_raw = hpt->Integral(hpt->FindBin(SB_LO1 + 1e-6), hpt->FindBin(SB_HI1 - 1e-6)) +
                             hpt->Integral(hpt->FindBin(SB_LO2 + 1e-6), hpt->FindBin(SB_HI2 - 1e-6));
      const double oob_ev = (w_sb > 0.) ? oob_raw * w_sig / w_sb : 0.;
      const double trap_amp = f2->GetParameter(1) - f2->GetParameter(0);
      TF1 ftrap((sp + "_trap_only" + tu + cc).c_str(), trapgaus, -150., 150., 5);
      ftrap.SetParameters(0., trap_amp, 0., f2->GetParameter(3), f2->GetParameter(4)); // flat=0, gaus off
      const double w_trap = 170. - (-25.);
      const double itb_raw = (p2binw > 0.) ? ftrap.Integral(-25., 170.) / p2binw : 0.;
      const double itb_ev = (w_trap > 0.) ? itb_raw * w_sig / w_trap : 0.;
      const double soob = (oob_ev > 0.) ? sig_ev / oob_ev : 0.;
      const double sitb = (itb_ev > 0.) ? sig_ev / itb_ev : 0.;
      const double itoob = (oob_ev > 0.) ? itb_ev / oob_ev : 0.;
      // No-cut reference from the pad-1 (proton total) fit f1: same signal and
      // out-of-time-bg definitions, on the histogram with no track/chi cut.
      // (No trapezoid in that fit, so there is no in-time-background analogue.)
      TH1D *hpp = h_proton_tof[is];
      const double p1binw = hpp->GetXaxis()->GetBinWidth(1);
      const double g_amp_nc = f1->GetParameter(1), g_sig_nc = std::fabs(f1->GetParameter(3));
      const double w_sig_nc = g_sig_nc * std::sqrt(TMath::TwoPi());
      const double sig_nc = (p1binw > 0.) ? g_amp_nc * w_sig_nc / p1binw : 0.;
      const double oob_raw_nc = hpp->Integral(hpp->FindBin(SB_LO1 + 1e-6), hpp->FindBin(SB_HI1 - 1e-6)) +
                                hpp->Integral(hpp->FindBin(SB_LO2 + 1e-6), hpp->FindBin(SB_HI2 - 1e-6));
      const double oob_nc = (w_sb > 0.) ? oob_raw_nc * w_sig_nc / w_sb : 0.;
      const double soob_nc = (oob_nc > 0.) ? sig_nc / oob_nc : 0.;
      // pad 3: track/total. Both histograms are rebinned by 5 BEFORE the divide,
      // so the Rebin scale factor cancels in the ratio (no extra /5 needed). The
      // background subtracted from each is now the NON-GAUSSIAN part of the pad
      // 1 & 2 fits (track: flat+trapezoid; total: flat), so only the gaussian
      // peak survives. Fixed to x in [-50,150], y in [0,3]; peak region shaded.
      c->cd(3);
      const double orig_binw = h_proton_track_tof[is][ic][it]->GetXaxis()->GetBinWidth(1);
      // background models = the fits with the gaussian height set to zero.
      TF1 *fbg_trk = new TF1((sp + "_bg_trk" + tu + cc).c_str(), trapgaus, -150., 150., 5);
      fbg_trk->SetParameters(f2->GetParameter(0), f2->GetParameter(1), 0., f2->GetParameter(3), f2->GetParameter(4));
      TF1 *fbg_tot = new TF1((sp + "_bg_tot" + tu + cc).c_str(), flatgaus, -150., 150., 4);
      fbg_tot->SetParameters(f1->GetParameter(0), 0., f1->GetParameter(2), f1->GetParameter(3));
      TH1D *ht_rb5r = (TH1D *)h_proton_track_tof[is][ic][it]->Clone((sp + "_proton_track_ratio_trk_rb5" + tu).c_str());
      ht_rb5r->Rebin(5);
      TH1D *hp_rb5r = (TH1D *)h_proton_tof[is]->Clone((sp + "_proton_track_ratio_tot_rb5" + tu).c_str());
      hp_rb5r->Rebin(5);
      // subtract the fit background into copies (keep ht_rb5r/hp_rb5r as raw
      // event counts for the weighting below).
      TH1D *ht_sb2 = (TH1D *)ht_rb5r->Clone((sp + "_proton_track_ratio_trk_fitbg" + tu).c_str());
      TH1D *hp_sb2 = (TH1D *)hp_rb5r->Clone((sp + "_proton_track_ratio_tot_fitbg" + tu).c_str());
      subtract_fit_bg(ht_sb2, fbg_trk, orig_binw);
      subtract_fit_bg(hp_sb2, fbg_tot, orig_binw);
      TH1D *hratio = (TH1D *)ht_sb2->Clone((sp + "_proton_track_ratio" + tu).c_str());
      hratio->SetTitle((sp + " proton (track-fitbg)/(total-fitbg);tof-L/c(ns);ratio").c_str());
      hratio->Divide(hp_sb2);
      // event-weighted peak-region efficiency: weight each ratio bin by its
      // numerator+denominator event counts (rebinned, pre-bg-sub).
      region_wstats(hratio, ht_rb5r, hp_rb5r, {{30., 50.}}, eff_m, eff_e);
      hratio->GetXaxis()->SetRangeUser(-50., 150.);
      hratio->GetYaxis()->SetRangeUser(0., 3);
      draw_shaded(hratio, {{30., 50., (double)kPk, (double)fsPk}});
      delete hratio;
      delete ht_sb2;
      delete hp_sb2;
      delete ht_rb5r;
      delete hp_rb5r;
      delete fbg_trk;
      delete fbg_tot;
      // pad 4: raw track/total ratio (NO background subtraction, and NOT
      // rebinned -- full binning), same axis ranges as pad 3; peak region shaded.
      c->cd(4);
      TH1D *ht_rawn = (TH1D *)h_proton_track_tof[is][ic][it]->Clone((sp + "_proton_track_ratio_raw_trk" + tu).c_str());
      TH1D *hp_rawn = (TH1D *)h_proton_tof[is]->Clone((sp + "_proton_track_ratio_raw_tot" + tu).c_str());
      TH1D *hratio_raw = (TH1D *)ht_rawn->Clone((sp + "_proton_track_ratio_raw" + tu).c_str());
      hratio_raw->SetTitle((sp + " proton track/total (no bg-sub);tof-L/c(ns);ratio").c_str());
      hratio_raw->Divide(hp_rawn);
      // sideband avg of the raw ratio = fraction of random background surviving
      // the track cut ("background reduction fraction").
      region_stats(hratio_raw, {{SB_LO1, SB_HI1}, {SB_LO2, SB_HI2}}, bkg_m, bkg_e);
      // x range left at the full default (like pads 1 & 2); y range fixed to [0,1].
      hratio_raw->GetYaxis()->SetRangeUser(0., 1.);
      draw_shaded(hratio_raw, {{SB_LO1, SB_HI1, (double)kSB, (double)fsSB},
                               {SB_LO2, SB_HI2, (double)kSB, (double)fsSB},
                               {30., 50., (double)kPk, (double)fsPk}});
      delete hratio_raw;
      delete ht_rawn;
      delete hp_rawn;
      // Record this variant's summary points (1D variants only). Names are the
      // same for every cut, so store them once (on the first cut).
      if (is1D) {
        if (ic == 0)
          sum_names.push_back(tracks[it].dir);
        sum_eff[ic].push_back(eff_m);
        sum_eff_err[ic].push_back(eff_e);
        sum_bkg[ic].push_back(bkg_m);
        sum_bkg_err[ic].push_back(bkg_e);
        sum_sig[ic].push_back(sig_ev);
        sum_oob[ic].push_back(oob_ev);
        sum_itb[ic].push_back(itb_ev);
        sum_soob[ic].push_back(soob);
        sum_sitb[ic].push_back(sitb);
        sum_itoob[ic].push_back(itoob);
        sum_sig_nc[ic].push_back(sig_nc);
        sum_oob_nc[ic].push_back(oob_nc);
        sum_soob_nc[ic].push_back(soob_nc);
        sum_zero[ic].push_back(0.);
      }
      wc(c);

      // Sideband-subtracted efficiency: (track-bg)/(total-bg)
      TCanvas *csb = new TCanvas((sp + "_c_proton_tof_sb").c_str(), (sp + " proton tof sb-sub").c_str(), 1800, 600);
      csb->Divide(3, 1);
      TH1D *hp_rb = (TH1D *)h_proton_tof[is]->Clone((sp + "_proton_rb" + tu).c_str());
      hp_rb->Rebin(PROTON_REBIN);
      TH1D *ht_rb = (TH1D *)h_proton_track_tof[is][ic][it]->Clone((sp + "_proton_track_rb" + tu).c_str());
      ht_rb->Rebin(PROTON_REBIN);
      TH1D *hpb = flat_bgsub(hp_rb);
      TH1D *htb = flat_bgsub(ht_rb);
      delete hp_rb;
      delete ht_rb;
      TH1D *hratio_sb = (TH1D *)htb->Clone((sp + "_proton_track_ratio_sb" + tu).c_str());
      hratio_sb->SetTitle((sp + " proton (track-bg)/(total-bg);tof-L/c(ns);ratio").c_str());
      hratio_sb->Divide(hpb);
      csb->cd(1);
      hpb->DrawCopy();
      csb->cd(2);
      htb->DrawCopy();
      csb->cd(3);
      hratio_sb->DrawCopy();
      wc(csb);
      delete hpb;
      delete htb;
      delete hratio_sb;

      // Per-paddle canvases for planes 001 and 101
      for (int pp = 0; pp < 2; ++pp) {
        const int pi = (pp == 0) ? 1 : 3;
        const std::string pn = plane_names[pi];

        TCanvas *cp =
            new TCanvas((sp + "_c_proton_" + pn).c_str(), (sp + " proton " + pn + " per paddle").c_str(), 1600, 1000);
        cp->Divide(4, 3);
        for (int pa = 0; pa < N_PADDLES; ++pa) {
          cp->cd(pa + 1);
          h_proton_pad[is][pp][pa]->DrawCopy();
        }
        wc(cp);

        TCanvas *ct = new TCanvas((sp + "_c_proton_track_" + pn).c_str(),
                                  (sp + " proton+track " + pn + " per paddle").c_str(), 1600, 1000);
        ct->Divide(4, 3);
        for (int pa = 0; pa < N_PADDLES; ++pa) {
          ct->cd(pa + 1);
          h_proton_track_pad[is][ic][it][pp][pa]->DrawCopy();
        }
        wc(ct);

        TCanvas *cr = new TCanvas((sp + "_c_proton_ratio_" + pn).c_str(),
                                  (sp + " proton ratio " + pn + " per paddle").c_str(), 1600, 1000);
        cr->Divide(4, 3);
        for (int pa = 0; pa < N_PADDLES; ++pa) {
          cr->cd(pa + 1);
          TH1D *hr = (TH1D *)h_proton_track_pad[is][ic][it][pp][pa]->Clone(
              (sp + "_proton_ratio_p" + std::to_string(pi) + "_b" + std::to_string(pa) + tu).c_str());
          hr->Divide(h_proton_pad[is][pp][pa]);
          hr->DrawCopy();
          delete hr;
        }
        wc(cr);
      }

      // Summary canvases: plane 001, plane 101, total
      TCanvas *csp2 = new TCanvas((sp + "_c_proton_sum").c_str(), (sp + " proton plane sums").c_str(), 1800, 600);
      csp2->Divide(3, 1);
      csp2->cd(1);
      h_proton_psum[is][0]->DrawCopy();
      csp2->cd(2);
      h_proton_psum[is][1]->DrawCopy();
      csp2->cd(3);
      h_proton_tof[is]->DrawCopy();
      wc(csp2);

      TCanvas *cst2 =
          new TCanvas((sp + "_c_proton_track_sum").c_str(), (sp + " proton+track plane sums").c_str(), 1800, 600);
      cst2->Divide(3, 1);
      cst2->cd(1);
      h_proton_track_psum[is][ic][it][0]->DrawCopy();
      cst2->cd(2);
      h_proton_track_psum[is][ic][it][1]->DrawCopy();
      cst2->cd(3);
      h_proton_track_tof[is][ic][it]->DrawCopy();
      wc(cst2);

      TCanvas *csr2 =
          new TCanvas((sp + "_c_proton_ratio_sum").c_str(), (sp + " proton ratio plane sums").c_str(), 1800, 600);
      csr2->Divide(3, 1);
      for (int pp = 0; pp < 2; ++pp) {
        csr2->cd(pp + 1);
        TH1D *hr = (TH1D *)h_proton_track_psum[is][ic][it][pp]->Clone(
            (sp + "_proton_ratio_sum_" + std::to_string(pp) + tu).c_str());
        hr->Divide(h_proton_psum[is][pp]);
        hr->DrawCopy();
        delete hr;
      }
      csr2->cd(3);
      {
        TH1D *hr = (TH1D *)h_proton_track_tof[is][ic][it]->Clone((sp + "_proton_ratio_total" + tu).c_str());
        hr->Divide(h_proton_tof[is]);
        hr->DrawCopy();
        delete hr;
      }
      wc(csr2);

      // edep vs inter-plane tof, 2x2: rows = front (top) / back (bottom) plane
      // edep; cols = proton cut (left) / proton+track cut (right).
      TCanvas *cpid = new TCanvas((sp + "_c_proton_edep_tof").c_str(),
                                  (sp + " edep vs tof [" + tracks[it].dir + ", chi2<" + cutstr + "]").c_str(), 1400,
                                  1200);
      cpid->Divide(2, 2);
      cpid->cd(1);
      h_pid_ef_pro[is]->DrawCopy("COLZ");
      cpid->cd(2);
      h_pid_ef_trk[is][ic][it]->DrawCopy("COLZ");
      cpid->cd(3);
      h_pid_eb_pro[is]->DrawCopy("COLZ");
      cpid->cd(4);
      h_pid_eb_trk[is][ic][it]->DrawCopy("COLZ");
      wc(cpid);
      } // end tracking-variant loop
    } // end chi-square-cut loop

    // -------------------------------------------------------------
    // Combined summary: one point per 1D tracking variant, with ALL chi-square
    // cuts overlaid on a single plot per metric (proton_id/summary/).
    //   efficiency        = peak-region [30,50] event-weighted average of the
    //                       bg-subtracted (track-bg)/(total-bg) ratio (pad 3),
    //                       weighted by numerator+denominator counts per bin.
    //   bkg reduction frac = sideband average of the raw track/total ratio
    //                       (pad 4, no bg subtraction) = fraction of random
    //                       background surviving the track cut.
    // The (weighted) stdev of the region's bin values is the error bar.
    if (!sum_names.empty()) {
      TDirectory *dsum = d->mkdir("summary");
      dsum->cd();
      const int nb = (int)sum_names.size();
      const int cutColor[3] = {kRed + 1, kBlack, kBlue + 1};
      const int cutMark[3] = {20, 21, 22};
      // Nominal cut index (scale 1.0) -- used as the single representative for
      // the optional "no cut" reference series (the proton-total plot has no
      // chi-square cut, so it needs only one representative width).
      int inom = 0;
      for (int ic = 0; ic < N_CUTS; ++ic)
        if (std::fabs(CHI_CUT_SCALES[ic] - 1.0) < 1e-9)
          inom = ic;
      // Overlay every cut's series on one canvas. The per-cut TH1s are leaked
      // (like the TBoxes/TF1s elsewhere) so they persist as canvas primitives
      // the legend can reference until the canvas is written. If ncval is given
      // (the no-track-cut reference from pad 1, proton total), one representative
      // series (nominal cut) is overlaid with open markers + a dashed line.
      auto draw_summary = [&](const std::string &cname, const std::string &ytitle, const std::string &ttl,
                              const std::array<std::vector<double>, N_CUTS> &val,
                              const std::array<std::vector<double>, N_CUTS> &err,
                              const std::array<std::vector<double>, N_CUTS> *ncval = nullptr) {
        double vmax = 0.;
        for (int ic = 0; ic < N_CUTS; ++ic)
          for (int i = 0; i < nb; ++i)
            vmax = std::max(vmax, val[ic][i] + err[ic][i]);
        if (ncval)
          for (int i = 0; i < nb; ++i)
            vmax = std::max(vmax, (*ncval)[inom][i]);
        TCanvas *cs = new TCanvas((sp + cname).c_str(), (sp + ttl).c_str(), 1000, 700);
        cs->SetBottomMargin(0.22);
        cs->SetGridy();
        TLegend *leg = new TLegend(0.74, 0.77, 0.98, 0.93);
        leg->SetHeader("chi-square cut");
        for (int ic = 0; ic < N_CUTS; ++ic) {
          const std::string cs2 = std::to_string((int)std::lround(CHI_CUT_SCALES[ic] * CHI_CUT_BASE));
          TH1D *hs = new TH1D((sp + cname + "_cut" + std::to_string(ic)).c_str(),
                              (sp + ttl + ";;" + ytitle).c_str(), nb, 0, nb);
          for (int i = 0; i < nb; ++i) {
            hs->SetBinContent(i + 1, val[ic][i]);
            hs->SetBinError(i + 1, err[ic][i]);
            hs->GetXaxis()->SetBinLabel(i + 1, sum_names[i].c_str());
          }
          hs->SetStats(0);
          hs->SetMarkerStyle(cutMark[ic % 3]);
          hs->SetMarkerSize(1.2);
          hs->SetMarkerColor(cutColor[ic % 3]);
          hs->SetLineColor(cutColor[ic % 3]);
          hs->SetLineWidth(2);
          hs->GetXaxis()->CenterLabels(kTRUE); // tick/label centered on the bin (= point) center
          hs->GetXaxis()->LabelsOption("v");   // vertical labels for the long names
          hs->GetYaxis()->SetRangeUser(0., vmax > 0. ? vmax * 1.15 : 1.);
          hs->Draw(ic == 0 ? "PE1" : "PE1 SAME"); // leaked on purpose (legend + canvas primitive)
          leg->AddEntry(hs, ("chi2 < " + cs2).c_str(), "pl");
        }
        if (ncval) {
          // No-track-cut reference (proton total, pad 1). One representative
          // series (nominal cut width), open markers + dashed line, leaked.
          TH1D *hnc = new TH1D((sp + cname + "_nocut").c_str(), (sp + ttl + ";;" + ytitle).c_str(), nb, 0, nb);
          for (int i = 0; i < nb; ++i) {
            hnc->SetBinContent(i + 1, (*ncval)[inom][i]);
            hnc->GetXaxis()->SetBinLabel(i + 1, sum_names[i].c_str());
          }
          hnc->SetStats(0);
          hnc->SetMarkerStyle(24); // open circle
          hnc->SetMarkerSize(1.3);
          hnc->SetMarkerColor(kGreen + 3);
          hnc->SetLineColor(kGreen + 3);
          hnc->SetLineWidth(2);
          hnc->SetLineStyle(2); // dashed
          hnc->GetXaxis()->CenterLabels(kTRUE);
          hnc->GetXaxis()->LabelsOption("v");
          hnc->GetYaxis()->SetRangeUser(0., vmax > 0. ? vmax * 1.15 : 1.);
          hnc->Draw("PL SAME"); // leaked on purpose
          leg->AddEntry(hnc, "no cut (proton total)", "pl");
        }
        leg->Draw();
        wc(cs);
      };
      draw_summary("_c_summary_eff", "<(track-bg)/(total-bg)> over [30,50]", " 1D tracking efficiency summary",
                   sum_eff, sum_eff_err);
      draw_summary("_c_summary_bkgfrac", "<track/total> over sidebands", " 1D background reduction summary", sum_bkg,
                   sum_bkg_err);
      // Signal / background metrics (all backgrounds normalized to the signal width).
      draw_summary("_c_summary_signal", "peak gaussian integral (events)", " 1D signal (peak) summary", sum_sig,
                   sum_zero, &sum_sig_nc);
      draw_summary("_c_summary_oobkg", "sideband integral, signal-width norm (events)",
                   " 1D out-of-time background summary", sum_oob, sum_zero, &sum_oob_nc);
      draw_summary("_c_summary_itbkg", "trapezoid [-25,170] integral, signal-width norm (events)",
                   " 1D in-time background summary", sum_itb, sum_zero);
      draw_summary("_c_summary_s_over_oob", "signal / out-of-time bkg", " 1D signal / out-of-time bkg summary",
                   sum_soob, sum_zero, &sum_soob_nc);
      draw_summary("_c_summary_s_over_itb", "signal / in-time bkg", " 1D signal / in-time bkg summary", sum_sitb,
                   sum_zero);
      draw_summary("_c_summary_it_over_oob", "in-time bkg / out-of-time bkg", " 1D in-time / out-of-time bkg summary",
                   sum_itoob, sum_zero);
    }
    fout.cd();
  } // end spec loop

  fout.Close();
  // In load mode the histograms live in the cache file; close it only now that
  // all plotting (which reads/clones them) is done.
  if (fcache) {
    fcache->Close();
    delete fcache;
    fcache = nullptr;
  }
  const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count();
  const int mins = (int)(elapsed / 60.);
  const double secs = elapsed - 60. * mins;
  std::cout << "[lad_tracking_eff] Done. Wrote " << out_file << "\n";
  std::printf("[lad_tracking_eff] Total time: %.1f s (%dm %.1fs)\n", elapsed, mins, secs);
}
