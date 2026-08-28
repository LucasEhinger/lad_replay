// lad_hodo_eff.C
// -------------------------------------------------------------------------
// Hodoscope-distribution macro. A sibling of lad_tracking_eff.C: it reads the
// SAME input files and uses the SAME proton selection (back-paddle proton cut:
// isProton_1 == 1 on the punchthrough planes 001/101), but produces a DIFFERENT
// set of plots -- it does NOT reproduce the tracking-eff plot families (GEM,
// cluster-ADC, punchthrough, edep-tof, summary, ...). It produces, per
// spectrometer / chi-square cut / tracking variant / hodoscope plane:
//
//   * y-position-on-the-bar distribution, one canvas per paddle (0..10).
//   * hodoscope number (paddle #) distribution.
//   * back-plane energy deposition.
//   * front-back inter-plane time difference (t_back - t_front).
//
// Each of those canvases has 6 panels -- the 4 tof regions (all / OOT / IT /
// peak) plus the two background-subtracted variants IT-OOT and peak-IT-OOT --
// and on every panel three overlaid marker histograms: all hits, hits WITH a
// track, and hits WITHOUT a track (= all - with-track). The background
// subtraction reuses lad_tracking_eff.C's scheme exactly (tof-window flat and
// trapezoid scale factors from the pad-2 proton+track fit).
//
// It ALSO reproduces the _c_proton_tof canvas (per spectrometer / cut / variant).
//
// NEW vs lad_hodo_dist.C: a 4-stage tracking-EFFICIENCY breakdown ("funnel").
// For each spectrometer / chi-square cut it plots, per tracking variant, the
// fraction of proton-cut hits surviving each successive gate:
//   (1) proton-cut hits                         [denominator]
//   (2) ... in an event with a reconstructed vertex (react.ok != 0)
//   (3) ... whose variant chiSquare is a valid fit (>= 0 and < 1e29)
//   (4) ... whose chiSquare falls in the [chi_lo, chi_hi) "has-track" window
// The gap between successive stages localizes where a variant loses hits: (1->2)
// the vertex gate (a hard gate for the 1D and vertex-using 2D variants; only
// informational for the noTrackVertex variants), (2->3) no GEM cluster / failed
// fit, (3->4) the chiSquare window cut. Stage 4 / stage 1 == the plotted
// with-track/all efficiency for the 1D variants. One canvas per (spec, cut):
//   <spec>/chi2cut_<val>/<spec>_c_eff_funnel
//
// Tracking variants and the chi-square-cut folders are the same as
// lad_tracking_eff.C. Directory layout:
//   <spec>/chi2cut_<val>/_c_eff_funnel                (NEW: efficiency funnel)
//   <spec>/chi2cut_<val>/<variant>/_c_proton_tof
//   <spec>/chi2cut_<val>/<variant>/plane_<001|101>/<the per-quantity canvases>
//
// Usage:
//   root -l -b -q 'lad_hodo_eff.C("input.dat","out.root")'
//   root -l -b -q 'lad_hodo_eff.C("input.dat","out.root",8)'          // 8 threads
//   root -l -b -q 'lad_hodo_eff.C("input.dat","out.root",8,"cache.root")' // + cache
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
#include <TH3D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMath.h>
#include <TPad.h>
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

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// =====================================================================
// Constants (shared with lad_tracking_eff.C where relevant)
// =====================================================================
const int NBINS_TCORR = 650; // 0.5 ns bins over [-150, 175] (proton_tof)
const double XMIN_TCORR = -150., XMAX_TCORR = 175.;

const int N_PADDLES = 11, N_SPECS = 2;
const int N_TRACKS = 13; // 4 legacy + 9 1D variants (compile-time capacity)
const int PROTON_REBIN = 10;

// Two-sided sidebands (tof-L/c ns) used by the proton_tof ratio pads.
const double SB_LO1 = -150., SB_HI1 = -100., SB_LO2 = 125., SB_HI2 = 175.;

// "Has track" chiSquare window per family (2D: any chi < cut; 1D: chi in [0,cut)).
const double CHI_CUT_2D = 100.0;
const double CHI_CUT_1D = 100.0;
const double CHI_CUT_BASE = 100.0;
const int N_CUTS = 3;
const std::array<double, N_CUTS> CHI_CUT_SCALES = {0.5, 1.0, 2.0};

// Coarse tof axis carried on every quantity histogram so the OOT/IT/peak regions
// can be projected out at plot time. 5-ns bins over [-150,175]: every region
// boundary (-150,-100,-25,30,50,125,175) lands on a bin edge.
const int TOF2_NBINS = 65;

// Quantity (X) axis binnings.
const int YP_NB = 20;                        // y-on-bar: 24 cm bins (4x coarser)
const double YP_LO = -240., YP_HI = 240.;    // bars span ~+/-210 cm
const int HN_NB = N_PADDLES;                  // hodoscope number: one bin/paddle
const double HN_LO = -0.5, HN_HI = N_PADDLES - 0.5;
const int ED_NB = 38;                         // back-plane edep (MeV), ~4x coarser (~3.95 MeV bins)
const double ED_LO = 0., ED_HI = 150.;
const int DT_NB = 40;                         // t_back - t_front (ns): 0.25 ns bins over [0,10]
const double DT_LO = 0., DT_HI = 10.;

// tof-L/c regions: all / out-of-time sideband / in-time sideband / peak.
const int N_GREG = 4;
const char *const GREG_NAME[N_GREG] = {"all", "oot", "it", "peak"};
const std::vector<std::vector<std::array<double, 2>>> GREG_INT = {
    {{-1e9, 1e9}},                        // all
    {{SB_LO1, SB_HI1}, {SB_LO2, SB_HI2}}, // oot: [-150,-100] u [125,175]
    {{-25., 30.}, {50., 125.}},           // it:  [-25,30] u [50,125]
    {{30., 50.}}};                        // peak: [30,50]

const double hodo_radii[5] = {615., 655.6, 523., 563.6, 615.}; // cm, by plane index
const char *const plane_names[5] = {"000", "001", "100", "101", "200"};
const std::array<char, N_SPECS> specs = {'P', 'H'};
const int PLANE_IDX[2] = {1, 3}; // the two punchthrough planes: 001 (idx 1), 101 (idx 3)

const char *DEFAULT_DAT_FILE = "../files/run-lists/all_C3_runlist_SHMS_13p5.dat";
const char *DEFAULT_OUT_FILE = "files/hodo_eff/hodo_eff_C3_SHMS_13p5_v1_PH.root";

// ----- fit models reused from lad_tracking_eff.C for the proton_tof canvas -----
// flat + fixed-corner trapezoid (corners -75,-25,50,125) + gaussian.
static double trapgaus(double *xx, double *p) {
  const double x = xx[0];
  const double flat = p[0];
  const double amp = p[1] - p[0];
  double trap;
  if (x < -75.)
    trap = 0.;
  else if (x < -25.)
    trap = amp * (x + 75.) / 50.;
  else if (x < 50.)
    trap = amp;
  else if (x < 125.)
    trap = amp * (125. - x) / 75.;
  else
    trap = 0.;
  const double g = p[2] * std::exp(-0.5 * std::pow((x - p[3]) / p[4], 2));
  return flat + trap + g;
}
static double flatgaus(double *xx, double *p) {
  const double x = xx[0];
  return p[0] + p[1] * std::exp(-0.5 * std::pow((x - p[2]) / p[3], 2));
}

// Track-cut OOT (accidental/flat) suppression factor, from the proton corrected-tof
// spectrum. Per tof bin b, f(b) = (all(b) - h(b)) / (all(b) - oot), where all = the
// no-track (proton-cut) spectrum (hAll), h = the with-track spectrum (hTrk), and oot
// = the with-track OOT-window mean level. f is the fraction of the flat/accidental
// background that survives the track requirement, so wherever the OOT template is
// subtracted from tof window W we use OOT * f_W instead of OOT. f_W is formed as a
// ratio of window sums over region 'reg' (GREG_INT index: 2 = IT, 3 = peak):
// f_W = (sum_W all - sum_W h) / (sum_W all - oot * N_W). The result is clamped to
// [0,1] (the track cut can only suppress accidentals) and falls back to 1.0 if the
// denominator is degenerate or the inputs are missing (unmodified subtraction).
static double oot_scale_f(const TH1 *hAll, const TH1 *hTrk, int reg) {
  if (!hAll || !hTrk)
    return 1.0;
  const TAxis *ax = hAll->GetXaxis();
  double sH = 0.;
  int nO = 0;
  for (const auto &iv : GREG_INT[1]) { // OOT window -> with-track flat level
    int b1 = ax->FindBin(iv[0] + 1e-6), b2 = ax->FindBin(iv[1] - 1e-6);
    for (int b = b1; b <= b2; ++b) {
      sH += hTrk->GetBinContent(b);
      ++nO;
    }
  }
  const double oot = (nO > 0) ? sH / nO : 0.;
  // ratio of window sums: f_W = (sum all - sum h) / (sum all - oot * N_W)
  double sumA = 0., sumH = 0.;
  int nW = 0;
  for (const auto &iv : GREG_INT[reg]) {
    int b1 = ax->FindBin(iv[0] + 1e-6), b2 = ax->FindBin(iv[1] - 1e-6);
    for (int b = b1; b <= b2; ++b) {
      sumA += hAll->GetBinContent(b);
      sumH += hTrk->GetBinContent(b);
      ++nW;
    }
  }
  const double den = sumA - oot * nW;
  if (std::fabs(den) < 1e-9)
    return 1.0; // sum all -> oot * N_W over the window: undefined
  double f = (sumA - sumH) / den;
  if (f < 0.)
    f = 0.;
  if (f > 1.)
    f = 1.;
  return f;
}

// A tracking variant: 'dir' output sub-directory; 'tsuf' suffix appended to
// "goodhit_chiSquare"; a hit is tracked when chiSquare is in [chi_lo, chi_hi).
struct Track {
  std::string dir, tsuf;
  double chi_lo, chi_hi;
};

void lad_hodo_eff(const char *dat_file = DEFAULT_DAT_FILE, const char *out_file = DEFAULT_OUT_FILE, int nthreads = 100,
                   const char *cache_file = "") {

  const auto t_start = std::chrono::steady_clock::now();
  gROOT->SetBatch(kTRUE);
  TH1::AddDirectory(kFALSE);
  if (nthreads > 0) {
    ROOT::EnableImplicitMT(nthreads);
    std::cout << "[lad_hodo_eff] implicit MT: " << nthreads << " threads\n";
  } else {
    ROOT::EnableImplicitMT();
    std::cout << "[lad_hodo_eff] implicit MT: all cores\n";
  }

  // ---------------------------------------------------------------
  // 1. TChain (same runlist mechanism as lad_tracking_eff.C)
  // ---------------------------------------------------------------
  TChain chain("T");
  std::string datlist;
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
  std::cout << "[lad_hodo_eff] entries: " << chain.GetEntries() << "\n";
  if (!chain.GetEntries()) {
    std::cerr << "empty chain\n";
    return;
  }

  // ---------------------------------------------------------------
  // 1b. Detect which tracking variants are present.
  // ---------------------------------------------------------------
  chain.LoadTree(0);
  auto has_branch = [&chain](const std::string &n) { return chain.GetBranch(n.c_str()) != nullptr; };
  std::vector<Track> tracks;
  tracks.push_back({"standard", "", -1e30, CHI_CUT_2D});
  {
    const std::array<Track, 12> variants = {{
        {"x", "_xz", -1e30, CHI_CUT_2D},
        {"noTrackVertex", "_noTrackVertex", -1e30, CHI_CUT_2D},
        {"noTrackVertex_x", "_noTrackVertex_xz", -1e30, CHI_CUT_2D},
        {"1D_x_GEM0", "_1D_xz_GEM0", 0., CHI_CUT_1D},
        {"1D_x_GEM1", "_1D_xz_GEM1", 0., CHI_CUT_1D},
        {"1D_x_GEMboth", "_1D_xz_GEMboth", 0., CHI_CUT_1D},
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
        std::cout << "[lad_hodo_eff] tracking variant '" << v.dir << "' absent in data; skipping\n";
    }
  }
  const int ntracks = (int)tracks.size();
  std::cout << "[lad_hodo_eff] tracking variants (" << ntracks << "):";
  for (const auto &t : tracks)
    std::cout << " " << t.dir;
  std::cout << "\n";

  // Vertex flag for the efficiency-funnel's stage-2 gate. react.ok (from the
  // spectrometer THcReactionPoint) tracks THcLADKine's fVertexModule->HasVertex()
  // -- the same condition that gates the 1D tracking. If it is absent the vertex
  // stage is collapsed (treated as always-passing) so the funnel still builds.
  bool has_react[N_SPECS];
  for (int is = 0; is < N_SPECS; ++is) {
    has_react[is] = has_branch(std::string(1, specs[is]) + ".react.ok");
    if (!has_react[is])
      std::cout << "[lad_hodo_eff] " << specs[is]
                << ".react.ok absent; funnel vertex stage collapsed (stage 2 = stage 1)\n";
  }

  // ---------------------------------------------------------------
  // 1c. Histogram cache decision (opt-in; identical scheme to lad_tracking_eff).
  // ---------------------------------------------------------------
  const char *CACHE_VERSION = "he_v1"; // he_v1: lad_hodo_dist hd_v2 + 4-stage efficiency funnel
  std::string sig = std::string("lad_hodo_eff;") + CACHE_VERSION + ";funnel=4;";
  sig += "tof=" + std::to_string(NBINS_TCORR) + "," + std::to_string(XMIN_TCORR) + "," + std::to_string(XMAX_TCORR) +
         ";tof2=" + std::to_string(TOF2_NBINS) + ";yp=" + std::to_string(YP_NB) + "," + std::to_string(YP_LO) + "," +
         std::to_string(YP_HI) + ";hn=" + std::to_string(HN_NB) + ";ed=" + std::to_string(ED_NB) + "," +
         std::to_string(ED_LO) + "," + std::to_string(ED_HI) + ";dt=" + std::to_string(DT_NB) + "," +
         std::to_string(DT_LO) + "," + std::to_string(DT_HI) + ";cuts=";
  for (int ic = 0; ic < N_CUTS; ++ic)
    sig += std::to_string(CHI_CUT_SCALES[ic]) + ",";
  sig += ";vars=";
  for (const auto &t : tracks)
    sig += t.dir + "|" + std::to_string(t.chi_lo) + "|" + std::to_string(t.chi_hi) + ",";
  sig += ";runlist=" + std::to_string((unsigned long long)std::hash<std::string>{}(datlist));

  const bool cache_on = (cache_file && cache_file[0] != '\0');
  bool load = false;
  TFile *fcache = nullptr;
  if (cache_on && !gSystem->AccessPathName(cache_file)) {
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
    std::cout << "[lad_hodo_eff] cache " << (load ? "HIT -> loading histograms, skipping event loop" : "MISS") << ": "
              << cache_file << "\n";

  // ---------------------------------------------------------------
  // 2. RDataFrame + cache-aware booking helpers
  // ---------------------------------------------------------------
  ROOT::RDataFrame rdf(chain);
  ROOT::RDF::RNode df = rdf;

  struct HBind1 {
    TH1D **slot;
    ROOT::RDF::RResultPtr<TH1D> res;
  };
  struct HBind2 {
    TH2D **slot;
    ROOT::RDF::RResultPtr<TH2D> res;
  };
  struct HBind3 {
    TH3D **slot;
    ROOT::RDF::RResultPtr<TH3D> res;
  };
  std::vector<HBind1> bind1;
  std::vector<HBind2> bind2;
  std::vector<HBind3> bind3;

  // Proton-tof 1D (name == fill column, like lad_tracking_eff.C).
  auto BKtof = [&](TH1D *&slot, const std::string &col, const std::string &nm, const std::string &tt) {
    if (load) {
      slot = dynamic_cast<TH1D *>(fcache->Get(nm.c_str()));
      return;
    }
    bind1.push_back({&slot, df.Histo1D({nm.c_str(), tt.c_str(), NBINS_TCORR, XMIN_TCORR, XMAX_TCORR}, col)});
  };
  // 4-bin efficiency-funnel TH1D (fill column pushes the stage numbers 1..4 a hit
  // reaches, so bin k counts hits surviving to stage k).
  auto BKfun = [&](TH1D *&slot, const std::string &col, const std::string &nm, const std::string &tt) {
    if (load) {
      slot = dynamic_cast<TH1D *>(fcache->Get(nm.c_str()));
      return;
    }
    bind1.push_back({&slot, df.Histo1D({nm.c_str(), tt.c_str(), 4, 0.5, 4.5}, col)});
  };
  // TH3D(paddle x ypos x tof).
  auto BK3 = [&](TH3D *&slot, const std::string &xcol, const std::string &ycol, const std::string &zcol,
                 const std::string &nm, const std::string &tt) {
    if (load) {
      slot = dynamic_cast<TH3D *>(fcache->Get(nm.c_str()));
      return;
    }
    bind3.push_back({&slot, df.Histo3D({nm.c_str(), tt.c_str(), HN_NB, HN_LO, HN_HI, YP_NB, YP_LO, YP_HI, TOF2_NBINS,
                                        XMIN_TCORR, XMAX_TCORR},
                                       xcol, ycol, zcol)});
  };
  // TH2D(quantity x tof) with a caller-supplied x binning.
  auto BK2Q = [&](TH2D *&slot, const std::string &xcol, const std::string &zcol, const std::string &nm,
                  const std::string &tt, int nx, double xlo, double xhi) {
    if (load) {
      slot = dynamic_cast<TH2D *>(fcache->Get(nm.c_str()));
      return;
    }
    bind2.push_back(
        {&slot, df.Histo2D({nm.c_str(), tt.c_str(), nx, xlo, xhi, TOF2_NBINS, XMIN_TCORR, XMAX_TCORR}, xcol, zcol)});
  };

  using RVd = ROOT::VecOps::RVec<double>;

  // Histogram storage.
  TH1D *h_proton_tof[N_SPECS] = {nullptr};
  TH1D *h_proton_track_tof[N_SPECS][N_CUTS][N_TRACKS] = {{{nullptr}}};
  // per (spec, plane) "all" (no-track) and per (spec, cut, variant, plane) "with track".
  TH3D *h_all_yp3[N_SPECS][2] = {{nullptr}};
  TH2D *h_all_ed2[N_SPECS][2] = {{nullptr}};
  TH2D *h_all_dt2[N_SPECS][2] = {{nullptr}};
  TH3D *h_trk_yp3[N_SPECS][N_CUTS][N_TRACKS][2] = {{{{nullptr}}}};
  TH2D *h_trk_ed2[N_SPECS][N_CUTS][N_TRACKS][2] = {{{{nullptr}}}};
  TH2D *h_trk_dt2[N_SPECS][N_CUTS][N_TRACKS][2] = {{{{nullptr}}}};
  // Efficiency funnel: a 4-bin TH1D per (spec, cut, variant). Bin k = # proton-cut
  // hits surviving to stage k (1=proton, 2=+vertex, 3=+valid fit, 4=+in-window).
  TH1D *h_funnel[N_SPECS][N_CUTS][N_TRACKS] = {{{nullptr}}};

#ifdef LAD_HAS_RDF_PROGRESSBAR
  if (!load)
    ROOT::RDF::Experimental::AddProgressBar(df);
#endif

  for (int is = 0; is < N_SPECS; ++is) {
    const std::string sp(1, specs[is]);
    const std::string pfx = sp + ".ladhod.goodhit_";
    if (!load) {
      df = df.Alias(sp + "_plane_1", pfx + "plane_1");
      df = df.Alias(sp + "_paddle_1", pfx + "paddle_1");
      df = df.Alias(sp + "_tof_1", pfx + "hit_tof_1");
      df = df.Alias(sp + "_ypos_1", pfx + "hit_ypos_1");
      df = df.Alias(sp + "_isProton_1", pfx + "isProton_1");
      df = df.Alias(sp + "_hittime_0", pfx + "hittime_0");
      df = df.Alias(sp + "_hittime_1", pfx + "hittime_1");
      df = df.Alias(sp + "_edepMeV_1", pfx + "hitedep_MeV_1");
      for (const auto &tk : tracks)
        df = df.Alias(sp + "_chiSquare" + tk.tsuf, pfx + "chiSquare" + tk.tsuf);
      // Vertex flag (per event) for the funnel's stage-2 gate.
      if (has_react[is])
        df = df.Alias(sp + "_vtxok", sp + ".react.ok");
      else
        df = df.Define(sp + "_vtxok", "1.0"); // no react.ok -> vertex stage collapses
    }

    // ---- proton-tagged corrected-tof columns (planes 001 & 101 combined) ----
    auto mk_proton = [&](const std::string &col, bool req_track, const std::string &chicol, double clo, double chi_hi) {
      if (load)
        return;
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
                           continue;
                         double R = hodo_radii[pi];
                         double dx = 22. * (pd1[i] - 6.);
                         double p2d = std::sqrt(y1[i] * y1[i] + dx * dx);
                         r.push_back(t1[i] - std::sqrt(p2d * p2d + R * R) / 100. / 0.3);
                       }
                       return r;
                     },
                     {sp + "_plane_1", sp + "_paddle_1", sp + "_ypos_1", sp + "_tof_1", sp + "_isProton_1", chicol});
    };
    mk_proton(sp + "_ptof", false, sp + "_chiSquare", 0., 0.);
    BKtof(h_proton_tof[is], sp + "_ptof", sp + "_ptof", (sp + " proton tof;tof-L/c(ns);Counts"));
    for (int ic = 0; ic < N_CUTS; ++ic)
      for (int it = 0; it < ntracks; ++it) {
        const std::string col = sp + "_pttof" + tracks[it].tsuf + "_cut" + std::to_string(ic);
        mk_proton(col, true, sp + "_chiSquare" + tracks[it].tsuf, tracks[it].chi_lo,
                  tracks[it].chi_hi * CHI_CUT_SCALES[ic]);
        BKtof(h_proton_track_tof[is][ic][it], col, col, (sp + " proton+track tof;tof-L/c(ns);Counts"));
      }

    // ---- efficiency-funnel column: for each proton-cut hit, push the stage
    //      numbers 1..4 it reaches, so the booked 4-bin TH1D's bin k = # hits
    //      surviving to stage k (cumulative, strictly nested). ----
    auto mk_funnel = [&](const std::string &col, const std::string &chicol, double clo, double chi_hi) {
      if (load)
        return;
      df = df.Define(col,
                     [clo, chi_hi](const RVd &pl1, const RVd &ip1, const RVd &chi, double vtxok) {
                       RVd r;
                       const bool vtx = (vtxok != 0.0);
                       for (size_t i = 0; i < pl1.size(); ++i) {
                         if (ip1[i] != 1.)
                           continue;
                         int pi = (int)std::round(pl1[i]);
                         if (pi != 1 && pi != 3)
                           continue;
                         r.push_back(1.); // stage 1: proton-cut hit
                         if (!vtx)
                           continue;
                         r.push_back(2.); // stage 2: + event has a vertex
                         const double c = (i < chi.size()) ? chi[i] : 1e30;
                         if (!(c >= 0. && c < 1e29))
                           continue;
                         r.push_back(3.); // stage 3: + a valid variant fit
                         if (c >= clo && c < chi_hi)
                           r.push_back(4.); // stage 4: + chiSquare in the has-track window
                       }
                       return r;
                     },
                     {sp + "_plane_1", sp + "_isProton_1", chicol, sp + "_vtxok"});
    };
    for (int ic = 0; ic < N_CUTS; ++ic)
      for (int it = 0; it < ntracks; ++it) {
        const std::string col = sp + "_funnel" + tracks[it].tsuf + "_cut" + std::to_string(ic);
        mk_funnel(col, sp + "_chiSquare" + tracks[it].tsuf, tracks[it].chi_lo,
                  tracks[it].chi_hi * CHI_CUT_SCALES[ic]);
        BKfun(h_funnel[is][ic][it], col, col, (sp + " eff funnel [" + tracks[it].dir + "];stage;hits"));
      }

    // ---- per-plane packed {tof, paddle, ypos, edep_back, dt} columns ----
    // One packed column per selection; five aligned unpacked RVec columns feed
    // the TH3D(paddle,ypos,tof) and the two TH2D(quantity,tof).
    auto mk_pack = [&](const std::string &col, int pi, bool req_track, const std::string &chicol, double clo,
                       double chi_hi) {
      if (load)
        return;
      const double R = hodo_radii[pi];
      df = df.Define(col,
                     [pi, R, req_track, clo, chi_hi](const RVd &pl1, const RVd &pd1, const RVd &yp, const RVd &tf,
                                                     const RVd &ip1, const RVd &ht0, const RVd &ht1, const RVd &eb,
                                                     const RVd &chi) {
                       RVd r;
                       for (size_t i = 0; i < pl1.size(); ++i) {
                         if (ip1[i] != 1.)
                           continue;
                         if ((int)std::round(pl1[i]) != pi)
                           continue;
                         if (req_track && !(chi[i] >= clo && chi[i] < chi_hi))
                           continue;
                         double dx = 22. * (pd1[i] - 6.), p2d = std::sqrt(yp[i] * yp[i] + dx * dx);
                         double tofc = tf[i] - std::sqrt(p2d * p2d + R * R) / 100. / 0.3;
                         r.push_back(tofc);
                         r.push_back(pd1[i]);
                         r.push_back(yp[i]);
                         r.push_back(i < eb.size() ? eb[i] : 0.);
                         r.push_back((i < ht1.size() && i < ht0.size()) ? (ht1[i] - ht0[i]) : 0.);
                       }
                       return r;
                     },
                     {sp + "_plane_1", sp + "_paddle_1", sp + "_ypos_1", sp + "_tof_1", sp + "_isProton_1",
                      sp + "_hittime_0", sp + "_hittime_1", sp + "_edepMeV_1", chicol});
    };
    auto unpack = [&](const std::string &pk) {
      if (load)
        return;
      df = df.Define(pk + "_pad", [](const RVd &p) { RVd r; for (size_t h = 0; h + 4 < p.size(); h += 5) r.push_back(p[h + 1]); return r; }, {pk});
      df = df.Define(pk + "_yp", [](const RVd &p) { RVd r; for (size_t h = 0; h + 4 < p.size(); h += 5) r.push_back(p[h + 2]); return r; }, {pk});
      df = df.Define(pk + "_tof", [](const RVd &p) { RVd r; for (size_t h = 0; h + 4 < p.size(); h += 5) r.push_back(p[h + 0]); return r; }, {pk});
      df = df.Define(pk + "_ed", [](const RVd &p) { RVd r; for (size_t h = 0; h + 4 < p.size(); h += 5) r.push_back(p[h + 3]); return r; }, {pk});
      df = df.Define(pk + "_dt", [](const RVd &p) { RVd r; for (size_t h = 0; h + 4 < p.size(); h += 5) r.push_back(p[h + 4]); return r; }, {pk});
    };

    for (int pl = 0; pl < 2; ++pl) {
      const int pi = PLANE_IDX[pl];
      const std::string ptag = "_pl" + std::to_string(pi);
      // "all" (no track): variant/cut independent.
      const std::string pka = sp + "_hpack_all" + ptag;
      mk_pack(pka, pi, false, sp + "_chiSquare", 0., 0.);
      unpack(pka);
      BK3(h_all_yp3[is][pl], pka + "_pad", pka + "_yp", pka + "_tof", sp + "_yp3_all" + ptag,
          (sp + " ypos-vs-tof (all) " + plane_names[pi]));
      BK2Q(h_all_ed2[is][pl], pka + "_ed", pka + "_tof", sp + "_ed2_all" + ptag,
           (sp + " edep-vs-tof (all) " + plane_names[pi]), ED_NB, ED_LO, ED_HI);
      BK2Q(h_all_dt2[is][pl], pka + "_dt", pka + "_tof", sp + "_dt2_all" + ptag,
           (sp + " dt-vs-tof (all) " + plane_names[pi]), DT_NB, DT_LO, DT_HI);
      // "with track": one per cut and variant.
      for (int ic = 0; ic < N_CUTS; ++ic)
        for (int it = 0; it < ntracks; ++it) {
          const std::string tag = tracks[it].dir + "_cut" + std::to_string(ic) + ptag;
          const std::string pk = sp + "_hpack_" + tag;
          mk_pack(pk, pi, true, sp + "_chiSquare" + tracks[it].tsuf, tracks[it].chi_lo,
                  tracks[it].chi_hi * CHI_CUT_SCALES[ic]);
          unpack(pk);
          BK3(h_trk_yp3[is][ic][it][pl], pk + "_pad", pk + "_yp", pk + "_tof", sp + "_yp3_" + tag, (sp + " ypos-vs-tof"));
          BK2Q(h_trk_ed2[is][ic][it][pl], pk + "_ed", pk + "_tof", sp + "_ed2_" + tag, (sp + " edep-vs-tof"), ED_NB,
               ED_LO, ED_HI);
          BK2Q(h_trk_dt2[is][ic][it][pl], pk + "_dt", pk + "_tof", sp + "_dt2_" + tag, (sp + " dt-vs-tof"), DT_NB, DT_LO,
               DT_HI);
        }
    }
  }

  // ---------------------------------------------------------------
  // 3. Trigger the single event loop and resolve slots (fill mode only).
  // ---------------------------------------------------------------
  std::vector<TH1 *> cache_list;
  if (!load) {
    std::cout << "[lad_hodo_eff] booked " << bind1.size() << " TH1D + " << bind2.size() << " TH2D + " << bind3.size()
              << " TH3D; running event loop...\n";
    for (auto &b : bind1) {
      *b.slot = b.res.GetPtr();
      cache_list.push_back(*b.slot);
    }
    for (auto &b : bind2) {
      *b.slot = b.res.GetPtr();
      cache_list.push_back(*b.slot);
    }
    for (auto &b : bind3) {
      *b.slot = b.res.GetPtr();
      cache_list.push_back(*b.slot);
    }
    std::fprintf(stderr, "\n");
    if (cache_on) {
      TFile fc(cache_file, "RECREATE");
      if (!fc.IsZombie()) {
        fc.cd();
        for (TH1 *h : cache_list)
          if (h)
            h->Write();
        TNamed sg("signature", sig.c_str());
        sg.Write();
        fc.Close();
        std::cout << "[lad_hodo_eff] wrote cache: " << cache_file << "\n";
      }
    }
  }

  // ---------------------------------------------------------------
  // 4. Plotting
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
  // shaded-band overlay + gaussian-integral label, reused for _c_proton_tof.
  auto draw_shaded = [](TH1D *h, const std::vector<std::array<double, 4>> &bands) {
    h->DrawCopy();
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
    h->DrawCopy("same");
    gPad->RedrawAxis();
  };
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
  auto subtract_fit_bg = [](TH1D *h, TF1 *fbg, double orig_binw) {
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
      double xlo = h->GetXaxis()->GetBinLowEdge(b), xhi = h->GetXaxis()->GetBinUpEdge(b);
      double bg = (orig_binw > 0.) ? fbg->Integral(xlo, xhi) / orig_binw : 0.;
      h->SetBinContent(b, h->GetBinContent(b) - bg);
    }
  };
  const int kSB = kGray + 2, kPk = kRed, fsSB = 3004, fsPk = 3005;

  // Unique-name counter for temporary derived histograms.
  long uqn = 0;
  auto uq = [&uqn]() { return std::string("hdtmp_") + std::to_string(uqn++); };

  // Project the tof windows of region `reg` from a (quantity x tof) TH2 into a
  // 1D quantity histogram (X = quantity, summed over the region's tof bins).
  auto regionX = [&](TH2 *h2, int reg) -> TH1D * {
    const int ny = h2->GetNbinsY();
    if (reg == 0)
      return (TH1D *)h2->ProjectionX(uq().c_str(), 1, ny);
    TH1D *out = nullptr;
    for (const auto &iv : GREG_INT[reg]) {
      int b1 = h2->GetYaxis()->FindBin(iv[0] + 1e-6);
      int b2 = h2->GetYaxis()->FindBin(iv[1] - 1e-6);
      if (b1 < 1)
        b1 = 1;
      if (b2 > ny)
        b2 = ny;
      TH1D *p = (TH1D *)h2->ProjectionX(uq().c_str(), b1, b2);
      if (!out)
        out = p;
      else {
        out->Add(p);
        delete p;
      }
    }
    return out;
  };
  // (ypos x tof) slice of a TH3(paddle,ypos,tof) at a fixed paddle bin.
  auto paddleSlice = [&](TH3 *h3, int pbin) -> TH2D * {
    const int ny = h3->GetNbinsY(), nz = h3->GetNbinsZ();
    TH2D *out = new TH2D(uq().c_str(), "", ny, h3->GetYaxis()->GetXmin(), h3->GetYaxis()->GetXmax(), nz,
                         h3->GetZaxis()->GetXmin(), h3->GetZaxis()->GetXmax());
    for (int iy = 1; iy <= ny; ++iy)
      for (int iz = 1; iz <= nz; ++iz)
        out->SetBinContent(iy, iz, h3->GetBinContent(pbin, iy, iz));
    return out;
  };
  // (paddle x tof) of a TH3(paddle,ypos,tof), summed over ypos.
  auto paddleProj = [&](TH3 *h3) -> TH2D * {
    const int nx = h3->GetNbinsX(), ny = h3->GetNbinsY(), nz = h3->GetNbinsZ();
    TH2D *out = new TH2D(uq().c_str(), "", nx, h3->GetXaxis()->GetXmin(), h3->GetXaxis()->GetXmax(), nz,
                         h3->GetZaxis()->GetXmin(), h3->GetZaxis()->GetXmax());
    for (int ix = 1; ix <= nx; ++ix)
      for (int iz = 1; iz <= nz; ++iz) {
        double s = 0.;
        for (int iy = 1; iy <= ny; ++iy)
          s += h3->GetBinContent(ix, iy, iz);
        out->SetBinContent(ix, iz, s);
      }
    return out;
  };

  // Build one 6-panel canvas for a quantity: base2D[3] = (quantity x tof) for the
  // all / with-track / no-track selections. Panels = all/OOT/IT/peak/IT-OOT/
  // peak-IT-OOT; every panel overlays the three selections as marker histograms.
  auto draw_qcanvas = [&](const std::string &cname, const std::string &ctitle, TH2 *base_all, TH2 *base_trk,
                          TH2 *base_ntr, double sflat_it, double sflat_pk, double strap, const std::string &xtitle) {
    if (!base_all || !base_trk || !base_ntr)
      return;
    TH2 *base[3] = {base_all, base_trk, base_ntr};
    const char *selname[3] = {"all", "with track", "no track"};
    const int selcol[3] = {kBlack, kRed + 1, kAzure + 2};
    const int selmk[3] = {20, 21, 24};
    const char *rlab[6] = {"all tof", "OOT", "IT", "peak", "IT-OOT", "peak-IT-OOT"};
    // reg[s][0..5]: all, oot, it, peak, IT-OOT, peak-IT-OOT for selection s.
    TH1D *reg[3][6];
    for (int s = 0; s < 3; ++s) {
      TH1D *hall = regionX(base[s], 0);
      TH1D *hoot = regionX(base[s], 1);
      TH1D *hit = regionX(base[s], 2);
      TH1D *hpk = regionX(base[s], 3);
      TH1D *hio = (TH1D *)hit->Clone(uq().c_str());
      hio->Add(hoot, -sflat_it);
      TH1D *hpb = (TH1D *)hpk->Clone(uq().c_str());
      hpb->Add(hoot, -sflat_pk);
      {
        TH1D *t = (TH1D *)hit->Clone(uq().c_str());
        t->Add(hoot, -sflat_it);
        hpb->Add(t, -strap);
        delete t;
      }
      reg[s][0] = hall;
      reg[s][1] = hoot;
      reg[s][2] = hit;
      reg[s][3] = hpk;
      reg[s][4] = hio;
      reg[s][5] = hpb;
    }
    TCanvas *c = new TCanvas(cname.c_str(), ctitle.c_str(), 1800, 1200);
    c->Divide(3, 2);
    std::vector<TH1D *> rtrash; // ratio clones, deleted after the canvas is written
    for (int p = 0; p < 6; ++p) {
      c->cd(p + 1);
      // Split this panel: upper pad = the three overlaid distributions, lower pad
      // = the with-track/all and no-track/all ratios (y limited to [0,1]).
      TPad *pMain = new TPad((cname + "_m" + std::to_string(p)).c_str(), "", 0., 0.30, 1., 1.);
      TPad *pRat = new TPad((cname + "_r" + std::to_string(p)).c_str(), "", 0., 0.0, 1., 0.30);
      pMain->SetBottomMargin(0.03);
      pRat->SetTopMargin(0.04);
      pRat->SetBottomMargin(0.30);
      pRat->SetGridy();
      pMain->Draw();
      pRat->Draw();
      // ---- upper: all / with-track / no-track overlaid ----
      pMain->cd();
      double ymax = 0., ymin = 0.;
      for (int s = 0; s < 3; ++s) {
        ymax = std::max(ymax, reg[s][p]->GetMaximum());
        ymin = std::min(ymin, reg[s][p]->GetMinimum());
      }
      TLegend *lg = new TLegend(0.66, 0.72, 0.98, 0.92);
      lg->SetFillStyle(0);
      for (int s = 0; s < 3; ++s) {
        TH1D *h = reg[s][p];
        h->SetStats(0);
        h->SetLineColor(selcol[s]);
        h->SetMarkerColor(selcol[s]);
        h->SetMarkerStyle(selmk[s]);
        h->SetMarkerSize(0.7);
        h->SetTitle((std::string(rlab[p]) + ";;# hits").c_str()); // x labels live on the ratio pad
        h->GetYaxis()->SetRangeUser(ymin < 0. ? 1.15 * ymin : 0., ymax > 0. ? 1.15 * ymax : 1.);
        h->GetXaxis()->SetLabelSize(0.);
        TH1 *drawn = h->DrawCopy(s == 0 ? "E1" : "E1 SAME");
        lg->AddEntry(drawn, selname[s], "lp");
      }
      lg->Draw();
      // ---- lower: ratios to "all" (with-track/all, no-track/all) ----
      pRat->cd();
      TH1D *rr[2];
      rr[0] = (TH1D *)reg[1][p]->Clone(uq().c_str()); // with-track / all
      rr[0]->Divide(reg[0][p]);
      rr[1] = (TH1D *)reg[2][p]->Clone(uq().c_str()); // no-track / all
      rr[1]->Divide(reg[0][p]);
      const int rcol[2] = {selcol[1], selcol[2]};
      const int rmk[2] = {selmk[1], selmk[2]};
      for (int k = 0; k < 2; ++k) {
        TH1D *r = rr[k];
        rtrash.push_back(r);
        r->SetStats(0);
        r->SetTitle((";" + xtitle + ";ratio /all").c_str());
        r->SetLineColor(rcol[k]);
        r->SetMarkerColor(rcol[k]);
        r->SetMarkerStyle(rmk[k]);
        r->SetMarkerSize(0.6);
        r->GetYaxis()->SetRangeUser(0., 1.);
        r->GetYaxis()->SetNdivisions(505);
        r->GetYaxis()->SetTitleSize(0.11);
        r->GetYaxis()->SetTitleOffset(0.42);
        r->GetYaxis()->SetLabelSize(0.10);
        r->GetXaxis()->SetTitleSize(0.11);
        r->GetXaxis()->SetTitleOffset(1.05);
        r->GetXaxis()->SetLabelSize(0.10);
        r->DrawCopy(k == 0 ? "E1" : "E1 SAME");
      }
    }
    wc(c);
    for (auto *h : rtrash)
      delete h;
    for (int s = 0; s < 3; ++s)
      for (int p = 0; p < 6; ++p)
        delete reg[s][p];
  };

  for (int is = 0; is < N_SPECS; ++is) {
    const std::string sp(1, specs[is]);
    TDirectory *sdir = fout.mkdir(sp.c_str());
    for (int ic = 0; ic < N_CUTS; ++ic) {
      const double cutval = CHI_CUT_SCALES[ic] * CHI_CUT_BASE;
      const std::string cutstr = std::to_string((int)std::lround(cutval));
      TDirectory *dc = sdir->mkdir(("chi2cut_" + cutstr).c_str());
      for (int it = 0; it < ntracks; ++it) {
        const std::string &tu = tracks[it].tsuf;
        const std::string cc = "_cut" + std::to_string(ic);
        TDirectory *td = mkpath(dc, tracks[it].dir);
        td->cd();

        // ---- fit the proton+track tof (pad-2 model): needed for BOTH the
        //      _c_proton_tof canvas and the background-subtraction scale factors.
        TH1D *hpt = h_proton_track_tof[is][ic][it];
        TF1 *f2 = new TF1((sp + "_fit_trapgaus" + tu + cc).c_str(), trapgaus, -150., 157., 5);
        f2->SetParNames("flat", "trap_top", "gaus_h", "gaus_mean", "gaus_sigma");
        f2->SetParameters(50., 70., 100., 41., 5.);
        f2->FixParameter(3, 41.);
        f2->SetLineColor(kGreen + 2);
        f2->SetNpx(600);
        if (hpt)
          hpt->Fit(f2, "RQN0");
        const double sig2 = f2->GetParameter(4);

        // ---- _c_proton_tof (2x2): pad1 proton total, pad2 proton+track, pad3
        //      (bg-subtracted) ratio, pad4 raw ratio. Reproduced from
        //      lad_tracking_eff.C.
        {
          TCanvas *c = new TCanvas((sp + "_c_proton_tof").c_str(),
                                   (sp + " proton tof [" + tracks[it].dir + ", chi2<" + cutstr + "]").c_str(), 1400,
                                   1000);
          c->Divide(2, 2);
          // pad 1: proton total, gaussian width fixed to the pad-2 width.
          c->cd(1);
          TH1D *hp1 = (TH1D *)h_proton_tof[is]->Clone((sp + "_ptof_p1" + tu + cc).c_str());
          TF1 *f1 = new TF1((sp + "_fit_flatgaus" + tu + cc).c_str(), flatgaus, -150., 157., 4);
          f1->SetParNames("flat", "gaus_h", "gaus_mean", "gaus_sigma");
          f1->SetParameters(50., 100., 41., sig2);
          f1->FixParameter(2, 41.);
          f1->FixParameter(3, sig2);
          f1->SetLineColor(kGreen + 2);
          f1->SetNpx(600);
          hp1->Fit(f1, "RQN0");
          draw_shaded(hp1, {{SB_LO1, SB_HI1, (double)kSB, (double)fsSB},
                            {SB_LO2, SB_HI2, (double)kSB, (double)fsSB},
                            {30., 50., (double)kPk, (double)fsPk}});
          f1->Draw("same");
          draw_gaus_integral(hp1, f1->GetParameter(1), f1->GetParameter(3));
          delete hp1;
          // pad 2: proton+track.
          c->cd(2);
          TH1D *ht2 = (TH1D *)hpt->Clone((sp + "_pttof_p2" + tu + cc).c_str());
          draw_shaded(ht2, {{SB_LO1, SB_HI1, (double)kSB, (double)fsSB},
                            {SB_LO2, SB_HI2, (double)kSB, (double)fsSB},
                            {30., 50., (double)kPk, (double)fsPk}});
          f2->Draw("same");
          draw_gaus_integral(ht2, f2->GetParameter(2), f2->GetParameter(4));
          delete ht2;
          // pad 3: (track - fitbg) / (total - fitbg), rebinned by 5.
          c->cd(3);
          const double orig_binw = hpt->GetXaxis()->GetBinWidth(1);
          TF1 *fbg_trk = new TF1((sp + "_bg_trk" + tu + cc).c_str(), trapgaus, -150., 157., 5);
          fbg_trk->SetParameters(f2->GetParameter(0), f2->GetParameter(1), 0., f2->GetParameter(3), f2->GetParameter(4));
          TF1 *fbg_tot = new TF1((sp + "_bg_tot" + tu + cc).c_str(), flatgaus, -150., 157., 4);
          fbg_tot->SetParameters(f1->GetParameter(0), 0., f1->GetParameter(2), f1->GetParameter(3));
          TH1D *ht_rb5 = (TH1D *)hpt->Clone((sp + "_ratio_trk_rb5" + tu + cc).c_str());
          ht_rb5->Rebin(5);
          TH1D *hp_rb5 = (TH1D *)h_proton_tof[is]->Clone((sp + "_ratio_tot_rb5" + tu + cc).c_str());
          hp_rb5->Rebin(5);
          subtract_fit_bg(ht_rb5, fbg_trk, orig_binw);
          subtract_fit_bg(hp_rb5, fbg_tot, orig_binw);
          TH1D *hratio = (TH1D *)ht_rb5->Clone((sp + "_ratio" + tu + cc).c_str());
          hratio->SetTitle((sp + " proton (track-fitbg)/(total-fitbg);tof-L/c(ns);ratio").c_str());
          hratio->Divide(hp_rb5);
          hratio->GetXaxis()->SetRangeUser(-50., 175.);
          hratio->GetYaxis()->SetRangeUser(0., 3.);
          draw_shaded(hratio, {{30., 50., (double)kPk, (double)fsPk}});
          delete hratio;
          delete ht_rb5;
          delete hp_rb5;
          delete fbg_trk;
          delete fbg_tot;
          // pad 4: raw track/total ratio (no bg-sub, full binning).
          c->cd(4);
          TH1D *ht_raw = (TH1D *)hpt->Clone((sp + "_ratio_raw_trk" + tu + cc).c_str());
          TH1D *hp_raw = (TH1D *)h_proton_tof[is]->Clone((sp + "_ratio_raw_tot" + tu + cc).c_str());
          TH1D *hratio_raw = (TH1D *)ht_raw->Clone((sp + "_ratio_raw" + tu + cc).c_str());
          hratio_raw->SetTitle((sp + " proton track/total (no bg-sub);tof-L/c(ns);ratio").c_str());
          hratio_raw->Divide(hp_raw);
          hratio_raw->GetYaxis()->SetRangeUser(0., 1.);
          draw_shaded(hratio_raw, {{SB_LO1, SB_HI1, (double)kSB, (double)fsSB},
                                   {SB_LO2, SB_HI2, (double)kSB, (double)fsSB},
                                   {30., 50., (double)kPk, (double)fsPk}});
          delete hratio_raw;
          delete ht_raw;
          delete hp_raw;
          wc(c);
        }

        // ---- background-subtraction scale factors from the pad-2 fit ----
        auto regW = [](int r) { double w = 0.; for (const auto &iv : GREG_INT[r]) w += iv[1] - iv[0]; return w; };
        const double wOOT = regW(1), wIT = regW(2), wPK = regW(3);
        // Track-cut OOT suppression from the proton tof spectrum (no-track vs with-track),
        // folded into the flat scale factors so every OOT subtraction becomes OOT * <f>.
        const double fIT = oot_scale_f(h_proton_tof[is], h_proton_track_tof[is][ic][it], 2);
        const double fPK = oot_scale_f(h_proton_tof[is], h_proton_track_tof[is][ic][it], 3);
        const double sflat_it = ((wOOT > 0.) ? wIT / wOOT : 0.) * fIT;
        const double sflat_pk = ((wOOT > 0.) ? wPK / wOOT : 0.) * fPK;
        TF1 ftr((sp + "_htrap" + tu + cc).c_str(), trapgaus, -150., 157., 5);
        ftr.SetParameters(0., f2->GetParameter(1) - f2->GetParameter(0), 0., f2->GetParameter(3), f2->GetParameter(4));
        auto Itr = [&](int r) { double s = 0.; for (const auto &iv : GREG_INT[r]) s += ftr.Integral(iv[0], iv[1]); return s; };
        const double itrIT = Itr(2), itrPK = Itr(3);
        const double strap = (itrIT != 0.) ? itrPK / itrIT : 0.;
        delete f2;

        // ---- per-plane quantity canvases ----
        for (int pl = 0; pl < 2; ++pl) {
          const int pi = PLANE_IDX[pl];
          const std::string pn = plane_names[pi];
          TDirectory *pdir = mkpath(td, "plane_" + pn);
          pdir->cd();
          const std::string tinfo = " [" + tracks[it].dir + ", chi2<" + cutstr + ", " + pn + "]";

          // no-track = all - with-track, for each quantity histogram.
          TH3 *y3[3] = {h_all_yp3[is][pl], h_trk_yp3[is][ic][it][pl], nullptr};
          TH2 *e2[3] = {h_all_ed2[is][pl], h_trk_ed2[is][ic][it][pl], nullptr};
          TH2 *d2[3] = {h_all_dt2[is][pl], h_trk_dt2[is][ic][it][pl], nullptr};
          if (!y3[0] || !y3[1] || !e2[0] || !e2[1] || !d2[0] || !d2[1])
            continue;
          y3[2] = (TH3 *)y3[0]->Clone(uq().c_str());
          y3[2]->Add(y3[1], -1.);
          e2[2] = (TH2 *)e2[0]->Clone(uq().c_str());
          e2[2]->Add(e2[1], -1.);
          d2[2] = (TH2 *)d2[0]->Clone(uq().c_str());
          d2[2]->Add(d2[1], -1.);

          // y-position on the bar, one canvas per paddle.
          for (int pa = 0; pa < N_PADDLES; ++pa) {
            TH2 *b[3];
            for (int s = 0; s < 3; ++s)
              b[s] = paddleSlice(y3[s], pa + 1);
            draw_qcanvas(sp + "_c_ypos_paddle" + std::to_string(pa),
                         (sp + " y-on-bar paddle " + std::to_string(pa) + tinfo).c_str(), b[0], b[1], b[2], sflat_it,
                         sflat_pk, strap, "y on bar (cm)");
            for (int s = 0; s < 3; ++s)
              delete b[s];
          }
          // hodoscope number (paddle #).
          {
            TH2 *b[3];
            for (int s = 0; s < 3; ++s)
              b[s] = paddleProj(y3[s]);
            draw_qcanvas(sp + "_c_hodo_num", (sp + " hodo number" + tinfo).c_str(), b[0], b[1], b[2], sflat_it, sflat_pk,
                         strap, "paddle #");
            for (int s = 0; s < 3; ++s)
              delete b[s];
          }
          // back-plane energy deposition.
          draw_qcanvas(sp + "_c_edep_back", (sp + " back edep" + tinfo).c_str(), e2[0], e2[1], e2[2], sflat_it, sflat_pk,
                       strap, "back edep (MeV)");
          // front-back inter-plane time difference.
          draw_qcanvas(sp + "_c_dt_frontback", (sp + " t_back - t_front" + tinfo).c_str(), d2[0], d2[1], d2[2], sflat_it,
                       sflat_pk, strap, "t_back - t_front (ns)");

          delete y3[2];
          delete e2[2];
          delete d2[2];
        }
      }

      // ================= efficiency funnel (this spec, this cut) =============
      // One canvas: x = tracking variant, four overlaid marker series giving the
      // fraction of proton-cut hits that survive each successive gate. The gap
      // between series localizes the loss: 1->2 vertex gate, 2->3 no fit/cluster,
      // 3->4 the chiSquare window. Series 4 / series 1 is the with-track/all
      // efficiency (exact for the 1D variants).
      {
        dc->cd();
        const std::string cc = "_cut" + std::to_string(ic);
        const int nv = ntracks;
        const char *stlab[4] = {"1: proton-cut", "2: + has vertex", "3: + valid fit", "4: + in window (eff)"};
        const int stcol[4] = {kGray + 2, kOrange + 7, kAzure + 1, kRed + 1};
        const int stmk[4] = {24, 25, 26, 20};
        TH1D *hs[4] = {nullptr, nullptr, nullptr, nullptr};
        double n_proton = 0.;
        for (int st = 0; st < 4; ++st) {
          hs[st] = new TH1D((sp + "_funnelfrac_s" + std::to_string(st + 1) + cc).c_str(),
                            (sp + " tracking-efficiency funnel [chi2<" + cutstr +
                             "];tracking variant;fraction of proton-cut hits")
                                .c_str(),
                            nv, 0.5, nv + 0.5);
          hs[st]->SetDirectory(nullptr);
          hs[st]->SetStats(0);
          hs[st]->SetMarkerStyle(stmk[st]);
          hs[st]->SetMarkerColor(stcol[st]);
          hs[st]->SetLineColor(stcol[st]);
          hs[st]->SetMarkerSize(1.3);
          hs[st]->SetLineWidth(2);
          for (int it = 0; it < nv; ++it) {
            hs[st]->GetXaxis()->SetBinLabel(it + 1, tracks[it].dir.c_str());
            TH1D *f = h_funnel[is][ic][it];
            const double s1 = f ? f->GetBinContent(1) : 0.;
            const double sk = f ? f->GetBinContent(st + 1) : 0.;
            if (st == 0)
              n_proton = std::max(n_proton, s1);
            const double frac = (s1 > 0.) ? sk / s1 : 0.;
            hs[st]->SetBinContent(it + 1, frac);
            if (s1 > 0.)
              hs[st]->SetBinError(it + 1, std::sqrt(std::max(0., frac * (1. - frac) / s1)));
          }
        }
        TCanvas *cf = new TCanvas((sp + "_c_eff_funnel").c_str(),
                                  (sp + " tracking-efficiency funnel [chi2<" + cutstr + "]").c_str(), 1600, 800);
        cf->SetGridy();
        cf->SetBottomMargin(0.24); // room for vertical variant labels
        hs[0]->GetYaxis()->SetRangeUser(0., 1.05);
        hs[0]->GetXaxis()->LabelsOption("v");
        for (int st = 0; st < 4; ++st)
          hs[st]->DrawCopy(st == 0 ? "PL E1" : "PL E1 SAME");
        TLegend *lg = new TLegend(0.68, 0.70, 0.99, 0.93);
        lg->SetHeader(Form("N proton-cut = %.0f", n_proton));
        lg->SetFillStyle(0);
        for (int st = 0; st < 4; ++st)
          lg->AddEntry(hs[st], stlab[st], "pl");
        lg->Draw();
        wc(cf);
        // compact numeric table to the log (plot is the primary artifact)
        std::printf("[lad_hodo_eff] %s chi2<%s funnel  (frac of %.0f proton-cut hits)\n", sp.c_str(), cutstr.c_str(),
                    n_proton);
        std::printf("    %-18s %8s %8s %8s\n", "variant", "vtx", "fit", "window");
        for (int it = 0; it < nv; ++it) {
          TH1D *f = h_funnel[is][ic][it];
          if (!f)
            continue;
          const double s1 = f->GetBinContent(1);
          std::printf("    %-18s %8.3f %8.3f %8.3f\n", tracks[it].dir.c_str(), s1 > 0 ? f->GetBinContent(2) / s1 : 0.,
                      s1 > 0 ? f->GetBinContent(3) / s1 : 0., s1 > 0 ? f->GetBinContent(4) / s1 : 0.);
        }
        for (int st = 0; st < 4; ++st)
          delete hs[st];
      }
    }
  }

  fout.Close();
  if (fcache) {
    fcache->Close();
    delete fcache;
    fcache = nullptr;
  }
  const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count();
  const int mins = (int)(elapsed / 60.);
  std::cout << "[lad_hodo_eff] Done. Wrote " << out_file << "\n";
  std::printf("[lad_hodo_eff] Total time: %.1f s (%dm %.1fs)\n", elapsed, mins, elapsed - 60. * mins);
}
