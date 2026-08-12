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
#include <TCanvas.h>
#include <TChain.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TROOT.h>
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
const double SB_LO1 = -150., SB_HI1 = -100., SB_LO2 = 100., SB_HI2 = 150.;

// "Has track" chiSquare window per family. 2D variants use a large "no track"
// sentinel, so any chiSquare < CHI_CUT_2D is a track. 1D variants use -1 as the
// sentinel, so a track is chiSquare in [0, CHI_CUT_1D).
const double CHI_CUT_2D = 100.0;
const double CHI_CUT_1D = 100.0;

// edep vs inter-plane tof (t_back - t_front) proton-ID 2x2 panels
const int NBINS_DT = 150;
const double XMIN_DT = -10., XMAX_DT = 30.; // t_back - t_front (ns)
const int NBINS_E2 = 200;
const double EMIN_F = 0., EMAX_F = 400.; // front-plane edep (MeV)
const double EMIN_B = 0., EMAX_B = 150.; // back-plane  edep (MeV)

const double hodo_radii[N_PLANES] = {615., 655.6, 523., 563.6, 615.}; // cm
const char *const plane_names[N_PLANES] = {"000", "001", "100", "101", "200"};
const std::array<char, N_SPECS> specs = {'P', 'H'};

const char *DEFAULT_DAT_FILE = "../files/run-lists/all_C3_runlist_22745-23590.dat";
const char *DEFAULT_OUT_FILE = "files/root_fast/tracking_eff_C3_22745-23590_PH.root";

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

void lad_tracking_eff(const char *dat_file = DEFAULT_DAT_FILE, const char *out_file = DEFAULT_OUT_FILE,
                      int nthreads = 100) {

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
  // 2. RDataFrame + aliases (only the columns the proton-ID plots need)
  // ---------------------------------------------------------------
  ROOT::RDataFrame rdf(chain);
  ROOT::RDF::RNode df = rdf;

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
    for (const auto &tk : tracks)
      mk_proton(sp + "_tof_corr_proton_track" + tk.tsuf, true, sp + "_chiSquare" + tk.tsuf, tk.chi_lo, tk.chi_hi);

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
        for (const auto &tk : tracks)
          mk_p(sp + "_tof_corr_proton_track_p" + std::to_string(pi_c) + "_b" + ps + tk.tsuf, true,
               sp + "_chiSquare" + tk.tsuf, tk.chi_lo, tk.chi_hi);
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
      for (const auto &tk : tracks)
        mk_sum(sp + "_tof_corr_proton_track_p" + std::to_string(pi_c) + "_sum" + tk.tsuf, true,
               sp + "_chiSquare" + tk.tsuf, tk.chi_lo, tk.chi_hi);
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
    for (const auto &tk : tracks)
      mk_pid("trk" + tk.tsuf, true, sp + "_chiSquare" + tk.tsuf, tk.chi_lo, tk.chi_hi);
  }

  // ---------------------------------------------------------------
  // 4. Book proton histograms
  // ---------------------------------------------------------------
  using RH1 = ROOT::RDF::RResultPtr<::TH1D>;
  std::array<RH1, N_SPECS> h_proton_tof;
  std::array<std::array<RH1, N_TRACKS>, N_SPECS> h_proton_track_tof;
  std::array<std::array<std::vector<RH1>, 2>, N_SPECS> h_proton_pad;
  std::array<std::array<std::array<std::vector<RH1>, 2>, N_TRACKS>, N_SPECS> h_proton_track_pad;
  std::array<std::array<RH1, 2>, N_SPECS> h_proton_psum;
  std::array<std::array<std::array<RH1, 2>, N_TRACKS>, N_SPECS> h_proton_track_psum;
  // edep-vs-tof 2D panels: proton-only (variant-independent) and proton+track.
  using RH2 = ROOT::RDF::RResultPtr<::TH2D>;
  std::array<RH2, N_SPECS> h_pid_ef_pro, h_pid_eb_pro;                 // front/back, proton cut
  std::array<std::array<RH2, N_TRACKS>, N_SPECS> h_pid_ef_trk, h_pid_eb_trk; // front/back, proton+track

  for (int is = 0; is < N_SPECS; ++is) {
    const std::string sp(1, specs[is]);
    auto bkp = [&](const std::string &col, const std::string &ttl) -> RH1 {
      return df.Histo1D({col.c_str(), ttl.c_str(), NBINS_TCORR, XMIN_TCORR, XMAX_TCORR}, col);
    };
    h_proton_tof[is] = bkp(sp + "_tof_corr_proton", sp + " tof corr proton;tof-L/c(ns);Counts");
    for (int it = 0; it < ntracks; ++it) {
      const std::string &td = tracks[it].dir;
      const std::string &tu = tracks[it].tsuf;
      h_proton_track_tof[is][it] =
          bkp(sp + "_tof_corr_proton_track" + tu, sp + " tof corr proton+track [" + td + "];tof-L/c(ns);Counts");
    }
    for (int pp = 0; pp < 2; ++pp) {
      const int pi = (pp == 0) ? 1 : 3;
      const std::string pn = plane_names[pi];
      h_proton_pad[is][pp].resize(N_PADDLES);
      for (int it = 0; it < ntracks; ++it)
        h_proton_track_pad[is][it][pp].resize(N_PADDLES);
      for (int pa = 0; pa < N_PADDLES; ++pa) {
        const std::string ps = std::to_string(pa);
        const std::string cp = sp + "_tof_corr_proton_p" + std::to_string(pi) + "_b" + ps;
        h_proton_pad[is][pp][pa] = bkp(cp, sp + " tof corr proton " + pn + " pd" + ps + ";tof-L/c(ns);Counts");
        for (int it = 0; it < ntracks; ++it) {
          const std::string &td = tracks[it].dir;
          const std::string &tu = tracks[it].tsuf;
          const std::string ct = sp + "_tof_corr_proton_track_p" + std::to_string(pi) + "_b" + ps + tu;
          h_proton_track_pad[is][it][pp][pa] =
              bkp(ct, sp + " tof corr proton+track [" + td + "] " + pn + " pd" + ps + ";tof-L/c(ns);Counts");
        }
      }
      const std::string cs2 = sp + "_tof_corr_proton_p" + std::to_string(pi) + "_sum";
      h_proton_psum[is][pp] = bkp(cs2, sp + " tof corr proton " + pn + " sum;tof-L/c(ns);Counts");
      for (int it = 0; it < ntracks; ++it) {
        const std::string &td = tracks[it].dir;
        const std::string &tu = tracks[it].tsuf;
        const std::string ts2 = sp + "_tof_corr_proton_track_p" + std::to_string(pi) + "_sum" + tu;
        h_proton_track_psum[is][it][pp] =
            bkp(ts2, sp + " tof corr proton+track [" + td + "] " + pn + " sum;tof-L/c(ns);Counts");
      }
    }

    // edep-vs-tof 2D panels
    auto bk2 = [&](const std::string &xcol, const std::string &ycol, const std::string &name, const std::string &ttl,
                   double ymin, double ymax) -> RH2 {
      return df.Histo2D({name.c_str(), ttl.c_str(), NBINS_DT, XMIN_DT, XMAX_DT, NBINS_E2, ymin, ymax}, xcol, ycol);
    };
    h_pid_ef_pro[is] = bk2(sp + "_pid_dt_pro", sp + "_pid_ef_pro", sp + "_pid_front_pro",
                           sp + " front edep vs tof (proton);t_{back}-t_{front} (ns);front edep (MeV)", EMIN_F, EMAX_F);
    h_pid_eb_pro[is] = bk2(sp + "_pid_dt_pro", sp + "_pid_eb_pro", sp + "_pid_back_pro",
                           sp + " back edep vs tof (proton);t_{back}-t_{front} (ns);back edep (MeV)", EMIN_B, EMAX_B);
    for (int it = 0; it < ntracks; ++it) {
      const std::string &td = tracks[it].dir;
      const std::string &tu = tracks[it].tsuf;
      h_pid_ef_trk[is][it] =
          bk2(sp + "_pid_dt_trk" + tu, sp + "_pid_ef_trk" + tu, sp + "_pid_front_trk" + tu,
              sp + " front edep vs tof (proton+track [" + td + "]);t_{back}-t_{front} (ns);front edep (MeV)", EMIN_F,
              EMAX_F);
      h_pid_eb_trk[is][it] =
          bk2(sp + "_pid_dt_trk" + tu, sp + "_pid_eb_trk" + tu, sp + "_pid_back_trk" + tu,
              sp + " back edep vs tof (proton+track [" + td + "]);t_{back}-t_{front} (ns);back edep (MeV)", EMIN_B,
              EMAX_B);
    }
  }

  // ---------------------------------------------------------------
  // 5. Trigger the event loop
  // ---------------------------------------------------------------
  std::cout << "[lad_tracking_eff] Running event loop...\n";
  (void)h_proton_tof[0]->GetEntries();
#ifndef LAD_HAS_RDF_PROGRESSBAR
  std::fprintf(stderr, "\n");
#endif

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

  for (int is = 0; is < N_SPECS; ++is) {
    const std::string sp(1, specs[is]);
    TDirectory *sdir = fout.mkdir(sp.c_str());

    TDirectory *d = sdir->mkdir("proton_id");
    for (int it = 0; it < ntracks; ++it) {
      const std::string &tu = tracks[it].tsuf;
      TDirectory *td = mkpath(d, tracks[it].dir);
      td->cd();

      TCanvas *c = new TCanvas((sp + "_c_proton_tof").c_str(), (sp + " proton tof").c_str(), 1800, 600);
      c->Divide(3, 1);
      c->cd(1);
      h_proton_tof[is]->DrawCopy();
      c->cd(2);
      h_proton_track_tof[is][it]->DrawCopy();
      c->cd(3);
      // track/total after subtracting a flat two-sideband background from each.
      TH1D *ht_sb2 = flat_bgsub2(h_proton_track_tof[is][it].GetPtr(), SB_LO1, SB_HI1, SB_LO2, SB_HI2);
      TH1D *hp_sb2 = flat_bgsub2(h_proton_tof[is].GetPtr(), SB_LO1, SB_HI1, SB_LO2, SB_HI2);
      TH1D *hratio = (TH1D *)ht_sb2->Clone((sp + "_proton_track_ratio" + tu).c_str());
      hratio->SetTitle((sp + " proton (track-bg)/(total-bg);tof-L/c(ns);ratio").c_str());
      hratio->Divide(hp_sb2);
      hratio->DrawCopy();
      delete hratio;
      delete ht_sb2;
      delete hp_sb2;
      wc(c);

      // Sideband-subtracted efficiency: (track-bg)/(total-bg)
      TCanvas *csb = new TCanvas((sp + "_c_proton_tof_sb").c_str(), (sp + " proton tof sb-sub").c_str(), 1800, 600);
      csb->Divide(3, 1);
      TH1D *hp_rb = (TH1D *)h_proton_tof[is]->Clone((sp + "_proton_rb" + tu).c_str());
      hp_rb->Rebin(PROTON_REBIN);
      TH1D *ht_rb = (TH1D *)h_proton_track_tof[is][it]->Clone((sp + "_proton_track_rb" + tu).c_str());
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
          h_proton_track_pad[is][it][pp][pa]->DrawCopy();
        }
        wc(ct);

        TCanvas *cr = new TCanvas((sp + "_c_proton_ratio_" + pn).c_str(),
                                  (sp + " proton ratio " + pn + " per paddle").c_str(), 1600, 1000);
        cr->Divide(4, 3);
        for (int pa = 0; pa < N_PADDLES; ++pa) {
          cr->cd(pa + 1);
          TH1D *hr = (TH1D *)h_proton_track_pad[is][it][pp][pa]->Clone(
              (sp + "_proton_ratio_p" + std::to_string(pi) + "_b" + std::to_string(pa) + tu).c_str());
          hr->Divide(h_proton_pad[is][pp][pa].GetPtr());
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
      h_proton_track_psum[is][it][0]->DrawCopy();
      cst2->cd(2);
      h_proton_track_psum[is][it][1]->DrawCopy();
      cst2->cd(3);
      h_proton_track_tof[is][it]->DrawCopy();
      wc(cst2);

      TCanvas *csr2 =
          new TCanvas((sp + "_c_proton_ratio_sum").c_str(), (sp + " proton ratio plane sums").c_str(), 1800, 600);
      csr2->Divide(3, 1);
      for (int pp = 0; pp < 2; ++pp) {
        csr2->cd(pp + 1);
        TH1D *hr = (TH1D *)h_proton_track_psum[is][it][pp]->Clone(
            (sp + "_proton_ratio_sum_" + std::to_string(pp) + tu).c_str());
        hr->Divide(h_proton_psum[is][pp].GetPtr());
        hr->DrawCopy();
        delete hr;
      }
      csr2->cd(3);
      {
        TH1D *hr = (TH1D *)h_proton_track_tof[is][it]->Clone((sp + "_proton_ratio_total" + tu).c_str());
        hr->Divide(h_proton_tof[is].GetPtr());
        hr->DrawCopy();
        delete hr;
      }
      wc(csr2);

      // edep vs inter-plane tof, 2x2: rows = front (top) / back (bottom) plane
      // edep; cols = proton cut (left) / proton+track cut (right).
      TCanvas *cpid = new TCanvas((sp + "_c_proton_edep_tof").c_str(),
                                  (sp + " edep vs tof [" + tracks[it].dir + "]").c_str(), 1400, 1200);
      cpid->Divide(2, 2);
      cpid->cd(1);
      h_pid_ef_pro[is]->DrawCopy("COLZ");
      cpid->cd(2);
      h_pid_ef_trk[is][it]->DrawCopy("COLZ");
      cpid->cd(3);
      h_pid_eb_pro[is]->DrawCopy("COLZ");
      cpid->cd(4);
      h_pid_eb_trk[is][it]->DrawCopy("COLZ");
      wc(cpid);
    } // end tracking-variant loop
    fout.cd();
  } // end spec loop

  fout.Close();
  std::cout << "[lad_tracking_eff] Done. Wrote " << out_file << "\n";
}
