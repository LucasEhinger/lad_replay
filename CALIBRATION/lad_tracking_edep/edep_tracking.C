// edep_tracking.C
// -------------------------------------------------------------------------
// GEM cluster-ADC ("edep") tracking discriminant for the LAD hodoscope.
//
// lad_tof/lad_tracking_eff.C tags a hodoscope good hit as "tracked" when a fitted
// track's chi-square is small.  This macro replaces that fit entirely: it draws the
// straight line from the spectrometer target vertex (<sp>.react.x/y/z) to the
// hodoscope good-hit position, intersects it with each GEM plane, and sums the ADC
// of every 1D cluster that lands close to the crossing point, weighted by how close.
// The per-DOF weighted sum
//
//     adc_sum / (number of (layer, axis) slots the variant uses)
//
// then replaces chi-square as the discriminant: a hit is "tracked" when the per-DOF
// sum is ABOVE the cut (the opposite sense to a chi-square cut).
//
// Tracking variants -- the same nine 1D combinations lad_tracking_eff.C uses:
//   1D_x_GEM0 / _GEM1 / _GEMboth      V (x) strips, front / back / both   (dof 1,1,2)
//   1D_y_GEM0 / _GEM1 / _GEMboth      U (y) strips, front / back / both   (dof 1,1,2)
//   1D_GEM0   / _GEM1 / _GEMboth      both axes,   front / back / both    (dof 2,2,4)
// No replay branch is needed for these -- the variants are built here from the raw
// clust.* arrays, so this macro works on any replay that wrote <sp>.gem.clust.* and
// <sp>.react.*, with or without ldo_1Dcluster_tracking.
//
// RESIDUAL CONVENTION.  Residuals are taken ALONG THE STRIP MEASUREMENT AXES, i.e.
// exactly the quantity THcLADKine::FitProj fits:
//     r = axisHat . (P_crossing - O_plane)  -  axisHat . (clusterLabPoint - O_plane)
// x uses V-hat (clust.axis == 1), y uses U-hat (clust.axis == 0).  This matters in x:
// V-hat is 52 deg from lab x (V-hat_x = 0.61), so a window along V is 1.64x the same
// number quoted in lab x.  It does not matter in y (U-hat_y = 0.9998).
//
// WEIGHTING.  The residual distribution is NOT the same shape in the two axes:
//   x  the 22 cm paddle gives a TOP-HAT of half-width 11 cm at the hodoscope, scaled
//      to the GEM by the lever arm, then smeared by the vertex-z resolution.  So a
//      flat core followed by a Gaussian shoulder is the correct shape.
//   y  the along-bar position comes from timing (sigma ~ 10 cm), so the projected
//      distribution is a PURE Gaussian with no flat core.
// Projecting both through the real geometry (see the table in WEIGHT MODE 1 below)
// gives sigma_y > sigma_x for every combination -- y is the WORSE coordinate, despite
// 10 cm < 22 cm, because the bar's RMS is 22/sqrt(12) = 6.35 cm, not 22 cm.
//
// Four weighting schemes are selectable with WEIGHT_MODE:
//   0  PER_LAYER   (default) geometry-derived, one (a, sigma) pair per GEM layer per
//                  axis, averaged over hodoscope plane.  Removes the lever-arm bias
//                  that would otherwise make GEM1 variants look worse than GEM0 for
//                  purely geometric reasons.  Gaussian shoulder, truncated at 3 sigma.
//   1  PER_LAYER_PLANE  same, but the full (layer x hodo plane) table.
//   2  FIXED_WIDE  one linear-taper window for everything: x 2.0 -> 3.5 cm,
//                  y 1.0 -> 4.0 cm (along the strip axes).
//   3  FIXED_ORIG  x 2.0 -> 3.0 cm, y 1.0 -> 2.0 cm.
// Modes 2 and 3 use a LINEAR taper, not 1/(r-a)^2: that form diverges as r -> a+ and
// equals 1 only at the outer edge, i.e. it up-weights hits just outside the core.
//
// Output (per spectrometer P / H):
//   <sp>/residuals/                     signed residual per GEM layer and axis, split
//                                       by tof region -- use these to check the widths
//                                       assumed by the weighting against the data
//   <sp>/adc_sum/<variant>/             adc_sum-per-dof spectra (cut lines drawn)
//   <sp>/proton_id/adccut_<v>/<variant> corrected-tof total / tracked / ratios
//   <sp>/summary/                       peak efficiency and background vs variant
//
// Usage:
//   root -l -b -q 'edep_tracking.C("input.dat","out.root")'
//   root -l -b -q 'edep_tracking.C("input.dat","out.root",8)'          // 8 MT threads
//   root -l -b -q 'edep_tracking.C("input.dat","out.root",8,"cache.root")'
// -------------------------------------------------------------------------

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <RVersion.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TNamed.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TString.h>
#include <TSystem.h>
#include <TVirtualPad.h>
#include <cmath>
#include <functional>

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
// Constants
// =====================================================================
const int N_SPECS = 2;
const std::array<char, N_SPECS> specs = {'P', 'H'};

// --- corrected-tof axis (identical to lad_tracking_eff.C so the plots overlay) ---
const int NBINS_TCORR = 650; // 0.5 ns bins over [-150, 175]
const double XMIN_TCORR = -150., XMAX_TCORR = 175.;
const int PROTON_REBIN = 10; // rebin to 5 ns before the sideband-subtracted ratio
const double SB_LO1 = -150., SB_HI1 = -100., SB_LO2 = 125., SB_HI2 = 175.;
const std::array<double, 2> PEAK_WIN = {30., 50.};

// tof regions the spectra are split by (a hit is in-region if its corrected tof
// falls in any of that region's intervals). Same definition as lad_tracking_eff.C.
const int N_GREG = 4;
const char *const GREG_NAME[N_GREG] = {"all", "oot", "it", "peak"};
const std::vector<std::vector<std::array<double, 2>>> GREG_INT = {
    {{-1e9, 1e9}},                        // all
    {{SB_LO1, SB_HI1}, {SB_LO2, SB_HI2}}, // oot: out-of-time sidebands
    {{-25., 30.}, {50., 125.}},           // it:  in-time sidebands
    {{30., 50.}}};                        // peak
// Hodoscope plane groups: plane 001 (index 1) and plane 101 (index 3).
const int N_GGRP = 2;
const char *const GGRP_NAME[N_GGRP] = {"000_001", "100_101"};
const char *const GGRP_TITLE[N_GGRP] = {"001", "101"};
const int GGRP_COLOR[N_GGRP] = {kRed + 1, kAzure + 2};

// --- hodoscope geometry (PARAM/LAD/HODO/lhodo_geom.param) ---
// Lab position of a hit, as THcLADHodoscope::GetHitPositionLab builds it:
//   (PosCenter(paddle), ypos, zpos) rotated about y by theta.
const int N_PLANES = 5;
const char *const plane_names[N_PLANES] = {"000", "001", "100", "101", "200"};
const double HODO_ZPOS[N_PLANES] = {618.8125006, 658.5516251, 526.1277344, 568.079082, 614.5200935};
const double HODO_THETA[N_PLANES] = {2.612452591, 2.615772334, 2.212325822, 2.213860777, 1.814590772};
// Paddle offsets. Each plane has 11 paddles of 22 cm; a hit is placed at the paddle
// CENTRE, so the offsets are the 11 values 110, 88, ... -88, -110 cm (the
// lladhodo_<plane>_center table in lhodo_geom.param), i.e. +/-110 about the middle
// paddle -- NOT 0 .. 220 measured from an edge.
// The index is 0-BASED: THcLADHodoscope stores goodhit paddles as
// GetPaddleNumber() - 1 (THcLADHodoscope.cxx:593/601/634), so goodhit_paddle_* runs
// 0..10, with -1 meaning "no hit". Confirmed against data (range -1 .. 10).
const int HODO_NPADDLE = 11;
const double HODO_POSCENTER_0 = 110.0;  // centre of paddle 0 (cm)
const double HODO_POSCENTER_DX = -22.0; // centre-to-centre step per paddle (cm)
// Radii used ONLY for the tof path-length correction (the transverse offset comes
// from HODO_POSCENTER_* above, so the tof path and the tracking line now use the
// same paddle geometry).
const double hodo_radii[N_PLANES] = {615., 655.6, 523., 563.6, 615.};

// --- GEM plane geometry, lab frame -------------------------------------
// Computed from PARAM/LAD/GEM/lgem_align_shms13p5.param + lgem_geom.param with the
// same construction as THcLADGEMModule::ComputeLabGeometry:
//   O = R_gem * position,  axisHat = (R_gem * R_module * (cos ang, sin ang, 0)).Unit(),
//   n = (R_gem * R_module * zhat).Unit(),   uangle = 90 deg, vangle = 180 deg.
// *** The module position/angle/rotation inputs come from the per-SHMS-setting
// *** alignment file and genuinely differ between settings (GEM0: 84.95 cm /
// *** 126.585 deg at SHMS 13.5, vs 85.91 cm / 127.008 deg at SHMS 17), so these
// *** constants are SHMS 13.5 only (runs >= 23426). The setting-independent part
// *** (uangle/vangle, pitch, nstrips) comes from lgem_geom.param. For the SHMS 17
// *** run list, regenerate O/N/V/U from lgem_align_shms17.param.
const int N_GEMLAY = 2;
struct GemPlane {
  double O[3], N[3], V[3], U[3]; // origin, normal, x-strip (V) axis, y-strip (U) axis
};
const GemPlane GEM[N_GEMLAY] = {
    // GEM0 (front, |O| = 84.96 cm)
    {{68.241528, 1.377531, -50.586652},
     {0.792455, -0.020872, -0.609573},
     {0.609729, 0.001398, 0.792609},
     {0.015691, 0.999781, -0.013835}},
    // GEM1 (back, |O| = 104.11 cm)
    {{83.655189, 1.308164, -61.953673},
     {0.794123, -0.024063, -0.607280},
     {0.607471, 0.000795, 0.794341},
     {0.018632, 0.999710, -0.015249}}};
// clust.axis values (THcLADKine::Do1DClusterTracking is authoritative): the x-z
// projection is measured by the V strips, the y projection by the U strips.
const int CLUST_AXIS_X = 1; // LADGEM::kVaxis
const int CLUST_AXIS_Y = 0; // LADGEM::kUaxis

// --- hit weighting ------------------------------------------------------
// WEIGHT_MODE selects the scheme (see the file header).
const int WM_PER_LAYER = 0, WM_PER_LAYER_PLANE = 1, WM_FIXED_WIDE = 2, WM_FIXED_ORIG = 3;
const int WEIGHT_MODE = WM_PER_LAYER;

// Logical axis index used by the weight tables and the slot numbering.
const int AX_X = 0, AX_Y = 1;

// Geometry-derived widths, in cm along the strip measurement axis. Obtained by
// propagating each hodoscope measurement through the exact geometry above:
//   x flat core = 11 cm (half a paddle) x d(residual_V)/d(hodo offset)
//   x sigma     = 1 cm  (vertex z resolution) x d(residual_V)/d(vertex z)
//   y sigma     = 10 cm (along-bar timing resolution) x d(residual_U)/d(hodo y)
// Full (layer, hodo plane) table, used by WM_PER_LAYER_PLANE:
//        combo          x core   x sigma   y core   y sigma
//   GEM0 / plane 001     1.66      0.47       0       1.39
//   GEM0 / plane 101     1.64      0.68       0       1.50
//   GEM1 / plane 001     2.04      0.45       0       1.71
//   GEM1 / plane 101     2.02      0.65       0       1.83
// Indexed [axis][layer][plane group].
const double W_CORE_LP[2][N_GEMLAY][N_GGRP] = {{{1.657, 1.645}, {2.035, 2.015}}, {{0.0, 0.0}, {0.0, 0.0}}};
const double W_SIG_LP[2][N_GEMLAY][N_GGRP] = {{{0.467, 0.681}, {0.450, 0.654}}, {{1.394, 1.495}, {1.710, 1.832}}};
// Per-layer table (the plane-averaged version), used by WM_PER_LAYER. Indexed [axis][layer].
const double W_CORE_L[2][N_GEMLAY] = {{1.651, 2.025}, {0.0, 0.0}};
const double W_SIG_L[2][N_GEMLAY] = {{0.574, 0.552}, {1.445, 1.771}};
// Multiple scattering between target and GEM (target foil + chamber window + ~1 m of
// air, ~4 mrad over ~90 cm), added in quadrature to the geometric sigma. Set to 0 to
// use pure geometry.
const double SIGMA_MS = 0.35; // cm
// Gaussian shoulder truncation, in units of sigma (weight there is 0.011).
const double W_NSIGMA = 3.0;
// Fixed windows for WM_FIXED_WIDE / WM_FIXED_ORIG: [mode][axis] flat to r1, linear to r2.
const double W_FIX_R1[2][2] = {{2.0, 1.0}, {2.0, 1.0}};
const double W_FIX_R2[2][2] = {{3.5, 4.0}, {3.0, 2.0}};

// --- adc_sum-per-dof cut ------------------------------------------------
// A hit is "tracked" when adc_sum/dof >= cut (note the sense is inverted relative to
// a chi-square cut). Cluster ADC sums span ~[400, 1e5] with a median near 3k, so the
// nominal cut sits just below the median single-cluster amplitude; the three cuts are
// the nominal, half and double, mirroring the three chi-square cuts of the tracking
// macro. Check the <sp>/adc_sum/ spectra before trusting the nominal value.
const double ADC_CUT_BASE = 2000.0;
const int N_CUTS = 3;
const std::array<double, N_CUTS> ADC_CUT_SCALES = {0.5, 1.0, 2.0};

// adc_sum-per-dof spectrum binning
const int ADC_NBINS = 250;
const double ADC_LO = 0., ADC_HI = 50000.;
// signed-residual spectrum binning (cm along the strip axis)
const int RES_NBINS = 200;
const double RES_LO = -10., RES_HI = 10.;

// Vertex validity: THaVertexModule::VertexClear sets the vertex to kBig (1e38) when
// there is no golden track, so a finite-range test is the validity test.
const double VTX_ZMAX = 100.0; // cm

// --- tracking variants ---------------------------------------------------
// Slot numbering, shared by the packed columns and the weight tables:
//   0 = GEM0 x (layer 0, V)   1 = GEM1 x (layer 1, V)
//   2 = GEM0 y (layer 0, U)   3 = GEM1 y (layer 1, U)
const int N_SLOT = 4;
const int SLOT_LAYER[N_SLOT] = {0, 1, 0, 1};
const int SLOT_AXIS[N_SLOT] = {AX_X, AX_X, AX_Y, AX_Y};
const char *const SLOT_NAME[N_SLOT] = {"GEM0_x", "GEM1_x", "GEM0_y", "GEM1_y"};

struct Variant {
  const char *dir;
  int nslot;
  int slot[N_SLOT];
};
const int N_VAR = 9;
const Variant VARIANTS[N_VAR] = {
    {"1D_x_GEM0", 1, {0, 0, 0, 0}},    {"1D_x_GEM1", 1, {1, 0, 0, 0}},    {"1D_x_GEMboth", 2, {0, 1, 0, 0}},
    {"1D_y_GEM0", 1, {2, 0, 0, 0}},    {"1D_y_GEM1", 1, {3, 0, 0, 0}},    {"1D_y_GEMboth", 2, {2, 3, 0, 0}},
    {"1D_GEM0", 2, {0, 2, 0, 0}},      {"1D_GEM1", 2, {1, 3, 0, 0}},      {"1D_GEMboth", 4, {0, 1, 2, 3}}};

const char *DEFAULT_DAT_FILE = "../files/run-lists/all_C3_runlist_SHMS_13p5.dat";
const char *DEFAULT_OUT_FILE = "files/edep_tracking/edep_tracking_C3_SHMS_13p5_PH.root";

// =====================================================================
// Helpers
// =====================================================================

// Weight of a cluster whose residual along its strip axis is r (cm), on GEM layer
// `lay`, logical axis `ax` (AX_X / AX_Y), for a hodoscope hit in plane group `pg`.
// Returns 0 outside the window, 1 inside the flat core.
static double hit_weight(double r, int lay, int ax, int pg) {
  const double a = std::fabs(r);
  if (WEIGHT_MODE == WM_FIXED_WIDE || WEIGHT_MODE == WM_FIXED_ORIG) {
    const int m = (WEIGHT_MODE == WM_FIXED_WIDE) ? 0 : 1;
    const double r1 = W_FIX_R1[m][ax], r2 = W_FIX_R2[m][ax];
    if (a <= r1)
      return 1.0;
    if (a >= r2)
      return 0.0;
    return (r2 - a) / (r2 - r1);
  }
  double core, sig;
  if (WEIGHT_MODE == WM_PER_LAYER_PLANE) {
    core = W_CORE_LP[ax][lay][pg];
    sig = W_SIG_LP[ax][lay][pg];
  } else { // WM_PER_LAYER
    core = W_CORE_L[ax][lay];
    sig = W_SIG_L[ax][lay];
  }
  sig = std::sqrt(sig * sig + SIGMA_MS * SIGMA_MS);
  if (a <= core)
    return 1.0;
  const double d = a - core;
  if (d > W_NSIGMA * sig)
    return 0.0;
  return std::exp(-0.5 * d * d / (sig * sig));
}

// Outer edge of the weight window (cm) -- beyond this the weight is exactly 0.
static double weight_reach(int lay, int ax, int pg) {
  if (WEIGHT_MODE == WM_FIXED_WIDE || WEIGHT_MODE == WM_FIXED_ORIG)
    return W_FIX_R2[(WEIGHT_MODE == WM_FIXED_WIDE) ? 0 : 1][ax];
  const double core = (WEIGHT_MODE == WM_PER_LAYER_PLANE) ? W_CORE_LP[ax][lay][pg] : W_CORE_L[ax][lay];
  double sig = (WEIGHT_MODE == WM_PER_LAYER_PLANE) ? W_SIG_LP[ax][lay][pg] : W_SIG_L[ax][lay];
  sig = std::sqrt(sig * sig + SIGMA_MS * SIGMA_MS);
  return core + W_NSIGMA * sig;
}

// Transverse offset (cm) of the CENTRE of a 0-based paddle, as
// THcLADHodoscopePlane::GetPosCenter returns it: 110, 88, ... -88, -110.
static double hodo_paddle_center(double paddle) {
  return HODO_POSCENTER_0 + HODO_POSCENTER_DX * paddle;
}

// True for a real paddle index (0 .. 10). goodhit_paddle_* uses -1 for "no hit".
static bool hodo_paddle_ok(double paddle) {
  const int p = (int)std::round(paddle);
  return p >= 0 && p < HODO_NPADDLE;
}

// Path length (cm) from the target to a hodoscope hit, used for the tof correction.
// The transverse offset is the paddle CENTRE, so this matches the geometry the
// tracking line uses. NOTE: lad_tof/lad_tracking_eff.C instead uses
// 22*(paddle - 6), which is the correct expression only for a 1-BASED paddle index;
// applied to the 0-based goodhit_paddle_* it spans -132 .. +88 cm rather than
// -110 .. +110, i.e. it is off by one paddle and asymmetric. The resulting tof
// difference is at most ~0.15 ns (well under the 0.5 ns bin width), so the axis here
// is still directly comparable to that macro's.
static double hodo_path_len(int plane, double paddle, double ypos) {
  const double dx = hodo_paddle_center(paddle);
  const double p2d = std::sqrt(ypos * ypos + dx * dx);
  return std::sqrt(p2d * p2d + hodo_radii[plane] * hodo_radii[plane]);
}

// Lab position of a hodoscope hit, matching THcLADHodoscope::GetHitPositionLab.
static void hodo_lab(int plane, double paddle, double ypos, double out[3]) {
  const double off = hodo_paddle_center(paddle);
  const double z = HODO_ZPOS[plane], th = HODO_THETA[plane];
  const double c = std::cos(th), s = std::sin(th);
  out[0] = off * c + z * s; // RotateY
  out[1] = ypos;
  out[2] = -off * s + z * c;
}

// Intersection of the ray vtx + t*dir (t > 0) with GEM plane `lay`. Returns false if
// the ray is parallel to the plane or crosses it behind the vertex.
static bool gem_cross(const double vtx[3], const double dir[3], int lay, double P[3]) {
  const GemPlane &g = GEM[lay];
  const double nd = g.N[0] * dir[0] + g.N[1] * dir[1] + g.N[2] * dir[2];
  if (std::fabs(nd) < 1e-6)
    return false;
  const double num =
      g.N[0] * (g.O[0] - vtx[0]) + g.N[1] * (g.O[1] - vtx[1]) + g.N[2] * (g.O[2] - vtx[2]);
  const double t = num / nd;
  if (t <= 0.)
    return false;
  for (int k = 0; k < 3; ++k)
    P[k] = vtx[k] + t * dir[k];
  return true;
}

// Coordinate of a lab point along a GEM plane's strip axis, measured from the plane
// origin -- the same "meas" THcLADGEMModule stores for a cluster.
static double along_axis(const double P[3], int lay, int ax) {
  const GemPlane &g = GEM[lay];
  const double *h = (ax == AX_X) ? g.V : g.U;
  return h[0] * (P[0] - g.O[0]) + h[1] * (P[1] - g.O[1]) + h[2] * (P[2] - g.O[2]);
}

// Two-sided sideband subtraction: flat background = mean bin content over the bins
// whose center falls in [lo1,hi1) or [lo2,hi2), subtracted from all bins.
static TH1D *flat_bgsub2(const TH1D *h, double lo1, double hi1, double lo2, double hi2) {
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

// Event-weighted mean (and weighted stdev) of a ratio histogram's bin values over one
// or more x-ranges, each bin weighted by numerator+denominator counts.
static void region_wstats(const TH1D *h, const TH1D *wnum, const TH1D *wden,
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
      continue;
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
}

// =====================================================================
// Main
// =====================================================================
void edep_tracking(const char *dat_file = DEFAULT_DAT_FILE, const char *out_file = DEFAULT_OUT_FILE,
                   int nthreads = 100, const char *cache_file = "") {

  const auto t_start = std::chrono::steady_clock::now();
  gROOT->SetBatch(kTRUE);
  TH1::AddDirectory(kFALSE);
  gStyle->SetOptStat(0);
  if (nthreads > 0) {
    ROOT::EnableImplicitMT(nthreads);
    std::cout << "[edep_tracking] implicit MT: " << nthreads << " threads\n";
  } else {
    ROOT::EnableImplicitMT();
    std::cout << "[edep_tracking] implicit MT: all cores\n";
  }
  {
    const char *wmname[4] = {"PER_LAYER", "PER_LAYER_PLANE", "FIXED_WIDE", "FIXED_ORIG"};
    std::cout << "[edep_tracking] weight mode: " << wmname[WEIGHT_MODE] << "; sigma_MS = " << SIGMA_MS << " cm\n";
    for (int lay = 0; lay < N_GEMLAY; ++lay)
      std::cout << "[edep_tracking]   GEM" << lay << " reach: x " << weight_reach(lay, AX_X, 0) << " cm, y "
                << weight_reach(lay, AX_Y, 0) << " cm\n";
  }

  // ---------------------------------------------------------------
  // 1. TChain
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
  std::cout << "[edep_tracking] entries: " << chain.GetEntries() << "\n";
  if (!chain.GetEntries()) {
    std::cerr << "empty chain\n";
    return;
  }

  // ---------------------------------------------------------------
  // 1b. Per-spectrometer branch availability. Everything here is built from the
  //     raw cluster arrays plus the hodoscope good hits and the vertex, so no
  //     tracking branch (chiSquare, trackid, clidx) is required.
  // ---------------------------------------------------------------
  chain.LoadTree(0);
  auto has_branch = [&chain](const std::string &n) { return chain.GetBranch(n.c_str()) != nullptr; };
  std::array<bool, N_SPECS> spec_ok{};
  for (int is = 0; is < N_SPECS; ++is) {
    const std::string sp(1, specs[is]);
    const std::string pfx = sp + ".ladhod.goodhit_";
    bool ok = has_branch(pfx + "plane_1") && has_branch(pfx + "paddle_1") && has_branch(pfx + "hit_ypos_1") &&
              has_branch(pfx + "hit_tof_1") && has_branch(pfx + "isProton_1");
    ok = ok && has_branch(sp + ".react.x") && has_branch(sp + ".react.y") && has_branch(sp + ".react.z");
    ok = ok && has_branch(sp + ".gem.clust.layer") && has_branch(sp + ".gem.clust.axis") &&
         has_branch(sp + ".gem.clust.adc") && has_branch(sp + ".gem.clust.labx") &&
         has_branch(sp + ".gem.clust.laby") && has_branch(sp + ".gem.clust.labz");
    spec_ok[is] = ok;
    std::cout << "[edep_tracking] spectrometer " << sp << ": " << (ok ? "enabled" : "missing branches, skipped")
              << "\n";
  }
  if (!spec_ok[0] && !spec_ok[1]) {
    std::cerr << "[edep_tracking] no spectrometer has the required branches\n";
    return;
  }

  // ---------------------------------------------------------------
  // 1c. Histogram cache. Same mechanism as lad_tracking_eff.C: when cache_file is
  //     given and its stored signature matches, load the histograms and skip the
  //     event loop entirely.
  // ---------------------------------------------------------------
  const char *CACHE_VERSION = "v1";
  std::string sig = std::string("edep_tracking;") + CACHE_VERSION + ";";
  sig += "tof=" + std::to_string(NBINS_TCORR) + "," + std::to_string(XMIN_TCORR) + "," + std::to_string(XMAX_TCORR) +
         ";adc=" + std::to_string(ADC_NBINS) + "," + std::to_string(ADC_LO) + "," + std::to_string(ADC_HI) +
         ";res=" + std::to_string(RES_NBINS) + "," + std::to_string(RES_LO) + "," + std::to_string(RES_HI) +
         ";wm=" + std::to_string(WEIGHT_MODE) + "," + std::to_string(SIGMA_MS) + "," + std::to_string(W_NSIGMA) +
         ";cut=" + std::to_string(ADC_CUT_BASE) + ",";
  for (int ic = 0; ic < N_CUTS; ++ic)
    sig += std::to_string(ADC_CUT_SCALES[ic]) + ",";
  sig += ";spec=" + std::to_string((int)spec_ok[0]) + std::to_string((int)spec_ok[1]);
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
    std::cout << "[edep_tracking] cache " << (load ? "HIT -> loading histograms, skipping event loop" : "MISS") << ": "
              << cache_file << "\n";

  // ---------------------------------------------------------------
  // 2. Histogram slots + cache-aware booking helpers
  // ---------------------------------------------------------------
  TH1D *h_tof_tot[N_SPECS] = {};                 // all proton hits (no vertex requirement)
  TH1D *h_tof_den[N_SPECS] = {};                 // proton hits WITH a valid vertex -- efficiency denominator
  TH1D *h_tof_trk[N_SPECS][N_CUTS][N_VAR] = {};  // + adc_sum/dof >= cut
  TH1D *h_res[N_SPECS][N_SLOT][N_GREG] = {};     // signed residual along the strip axis
  TH1D *h_adc[N_SPECS][N_VAR][N_GREG][N_GGRP] = {};

  ROOT::RDataFrame rdf(chain);
  ROOT::RDF::RNode df = rdf;

  struct HBind1 {
    TH1D **slot;
    ROOT::RDF::RResultPtr<TH1D> res;
  };
  std::vector<HBind1> bind1;
  auto BK = [&](TH1D *&slot, const std::string &col, const std::string &nm, const std::string &tt, int nb, double lo,
                double hi) {
    if (load) {
      slot = dynamic_cast<TH1D *>(fcache->Get(nm.c_str()));
      return;
    }
    bind1.push_back({&slot, df.Histo1D({nm.c_str(), tt.c_str(), nb, lo, hi}, col)});
  };

  using RVd = ROOT::VecOps::RVec<double>;

  if (!load) {
#ifdef LAD_HAS_RDF_PROGRESSBAR
    ROOT::RDF::Experimental::AddProgressBar(df);
#else
    ROOT::RDF::RResultPtr<ULong64_t> progress_count;
    {
      const ULong64_t total = chain.GetEntries();
      const ULong64_t step = std::max<ULong64_t>(1ULL, total / 200ULL);
      auto t0 = std::make_shared<std::chrono::steady_clock::time_point>(std::chrono::steady_clock::now());
      progress_count = df.Count();
      progress_count.OnPartialResult(step, [total, t0](ULong64_t n) {
        double f = total ? double(n) / double(total) : 0.;
        double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - *t0).count();
        std::fprintf(stderr, "\r[edep_tracking] %6.2f%% (%llu/%llu) %.1fs eta %.1fs   ", 100. * f,
                     (unsigned long long)n, (unsigned long long)total, s, f > 0 ? s * (1. / f - 1.) : 0.);
        std::fflush(stderr);
      });
    }
#endif

    for (int is = 0; is < N_SPECS; ++is) {
      if (!spec_ok[is])
        continue;
      const std::string sp(1, specs[is]);
      const std::string pfx = sp + ".ladhod.goodhit_";
      df = df.Alias(sp + "_plane_1", pfx + "plane_1");
      df = df.Alias(sp + "_paddle_1", pfx + "paddle_1");
      df = df.Alias(sp + "_ypos_1", pfx + "hit_ypos_1");
      df = df.Alias(sp + "_tof_1", pfx + "hit_tof_1");
      df = df.Alias(sp + "_isProton_1", pfx + "isProton_1");

      // -------------------------------------------------------------
      // 3. The two packed columns. Both loop over the proton good hits on planes
      //    001/101, build the vertex -> hodoscope line, cross it with each GEM
      //    plane and compare every cluster to the crossing point.
      //
      //    <sp>_epack  stride 2 + N_VAR per hit:
      //        [0] corrected tof, [1] plane group (0 = 001, 1 = 101),
      //        [2 + iv] adc_sum / dof for variant iv, or -1 when any slot the
      //                 variant needs has no cluster inside its window.
      //    <sp>_rpack  stride 4 per (hit, nearby cluster):
      //        [0] corrected tof, [1] slot index, [2] signed residual, [3] ADC.
      // -------------------------------------------------------------
      df = df.Define(
          sp + "_epack",
          [](const RVd &pl1, const RVd &pd1, const RVd &yp, const RVd &t1, const RVd &ip1, double vx, double vy,
             double vz, const RVd &clay, const RVd &cax, const RVd &cadc, const RVd &clx, const RVd &cly,
             const RVd &clz) {
            RVd r;
            if (!(std::fabs(vz) < VTX_ZMAX && std::fabs(vx) < VTX_ZMAX && std::fabs(vy) < VTX_ZMAX))
              return r; // no golden-track vertex -> no line to draw
            const double vtx[3] = {vx, vy, vz};
            for (size_t i = 0; i < pl1.size(); ++i) {
              if (ip1[i] != 1.)
                continue;
              const int pi = (int)std::round(pl1[i]);
              if (pi != 1 && pi != 3)
                continue; // planes 001 and 101 only
              if (!hodo_paddle_ok(pd1[i]))
                continue;
              const int pg = (pi == 1) ? 0 : 1;
              const double tofc = t1[i] - hodo_path_len(pi, pd1[i], yp[i]) / 100. / 0.3;
              // vertex -> hodoscope line
              double H[3];
              hodo_lab(pi, pd1[i], yp[i], H);
              double d[3] = {H[0] - vtx[0], H[1] - vtx[1], H[2] - vtx[2]};
              const double dm = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
              if (dm < 1e-6)
                continue;
              for (int k = 0; k < 3; ++k)
                d[k] /= dm;
              // predicted strip coordinate at each GEM plane
              double proj[N_SLOT];
              bool pok[N_SLOT] = {false, false, false, false};
              for (int lay = 0; lay < N_GEMLAY; ++lay) {
                double P[3];
                if (!gem_cross(vtx, d, lay, P))
                  continue;
                for (int s = 0; s < N_SLOT; ++s)
                  if (SLOT_LAYER[s] == lay) {
                    proj[s] = along_axis(P, lay, SLOT_AXIS[s]);
                    pok[s] = true;
                  }
              }
              // weighted ADC sum per slot
              double slot_sum[N_SLOT] = {0., 0., 0., 0.};
              bool slot_hit[N_SLOT] = {false, false, false, false};
              for (size_t j = 0; j < clay.size(); ++j) {
                const int lay = (int)std::llround(clay[j]);
                if (lay < 0 || lay >= N_GEMLAY)
                  continue;
                const int cx = (int)std::llround(cax[j]);
                const int ax = (cx == CLUST_AXIS_X) ? AX_X : ((cx == CLUST_AXIS_Y) ? AX_Y : -1);
                if (ax < 0)
                  continue;
                int s = -1;
                for (int k = 0; k < N_SLOT; ++k)
                  if (SLOT_LAYER[k] == lay && SLOT_AXIS[k] == ax)
                    s = k;
                if (s < 0 || !pok[s])
                  continue;
                const double lab[3] = {clx[j], cly[j], clz[j]};
                const double meas = along_axis(lab, lay, ax);
                const double w = hit_weight(proj[s] - meas, lay, ax, pg);
                if (w <= 0.)
                  continue;
                slot_sum[s] += w * ((j < cadc.size()) ? cadc[j] : 0.);
                slot_hit[s] = true;
              }
              r.push_back(tofc);
              r.push_back((double)pg);
              for (int iv = 0; iv < N_VAR; ++iv) {
                const Variant &V = VARIANTS[iv];
                double s = 0.;
                bool ok = true;
                for (int k = 0; k < V.nslot; ++k) {
                  if (!slot_hit[V.slot[k]]) {
                    ok = false;
                    break;
                  }
                  s += slot_sum[V.slot[k]];
                }
                // per-DOF: divide by the number of (layer, axis) slots the variant uses
                r.push_back(ok ? s / (double)V.nslot : -1.);
              }
            }
            return r;
          },
          {sp + "_plane_1", sp + "_paddle_1", sp + "_ypos_1", sp + "_tof_1", sp + "_isProton_1", sp + ".react.x",
           sp + ".react.y", sp + ".react.z", sp + ".gem.clust.layer", sp + ".gem.clust.axis", sp + ".gem.clust.adc",
           sp + ".gem.clust.labx", sp + ".gem.clust.laby", sp + ".gem.clust.labz"});

      df = df.Define(
          sp + "_rpack",
          [](const RVd &pl1, const RVd &pd1, const RVd &yp, const RVd &t1, const RVd &ip1, double vx, double vy,
             double vz, const RVd &clay, const RVd &cax, const RVd &cadc, const RVd &clx, const RVd &cly,
             const RVd &clz) {
            RVd r;
            if (!(std::fabs(vz) < VTX_ZMAX && std::fabs(vx) < VTX_ZMAX && std::fabs(vy) < VTX_ZMAX))
              return r;
            const double vtx[3] = {vx, vy, vz};
            for (size_t i = 0; i < pl1.size(); ++i) {
              if (ip1[i] != 1.)
                continue;
              const int pi = (int)std::round(pl1[i]);
              if (pi != 1 && pi != 3)
                continue;
              if (!hodo_paddle_ok(pd1[i]))
                continue;
              const double tofc = t1[i] - hodo_path_len(pi, pd1[i], yp[i]) / 100. / 0.3;
              double H[3];
              hodo_lab(pi, pd1[i], yp[i], H);
              double d[3] = {H[0] - vtx[0], H[1] - vtx[1], H[2] - vtx[2]};
              const double dm = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
              if (dm < 1e-6)
                continue;
              for (int k = 0; k < 3; ++k)
                d[k] /= dm;
              double proj[N_SLOT];
              bool pok[N_SLOT] = {false, false, false, false};
              for (int lay = 0; lay < N_GEMLAY; ++lay) {
                double P[3];
                if (!gem_cross(vtx, d, lay, P))
                  continue;
                for (int s = 0; s < N_SLOT; ++s)
                  if (SLOT_LAYER[s] == lay) {
                    proj[s] = along_axis(P, lay, SLOT_AXIS[s]);
                    pok[s] = true;
                  }
              }
              for (size_t j = 0; j < clay.size(); ++j) {
                const int lay = (int)std::llround(clay[j]);
                if (lay < 0 || lay >= N_GEMLAY)
                  continue;
                const int cx = (int)std::llround(cax[j]);
                const int ax = (cx == CLUST_AXIS_X) ? AX_X : ((cx == CLUST_AXIS_Y) ? AX_Y : -1);
                if (ax < 0)
                  continue;
                int s = -1;
                for (int k = 0; k < N_SLOT; ++k)
                  if (SLOT_LAYER[k] == lay && SLOT_AXIS[k] == ax)
                    s = k;
                if (s < 0 || !pok[s])
                  continue;
                const double lab[3] = {clx[j], cly[j], clz[j]};
                const double res = proj[s] - along_axis(lab, lay, ax);
                if (res <= RES_LO || res >= RES_HI)
                  continue;
                r.push_back(tofc);
                r.push_back((double)s);
                r.push_back(res);
                r.push_back((j < cadc.size()) ? cadc[j] : 0.);
              }
            }
            return r;
          },
          {sp + "_plane_1", sp + "_paddle_1", sp + "_ypos_1", sp + "_tof_1", sp + "_isProton_1", sp + ".react.x",
           sp + ".react.y", sp + ".react.z", sp + ".gem.clust.layer", sp + ".gem.clust.axis", sp + ".gem.clust.adc",
           sp + ".gem.clust.labx", sp + ".gem.clust.laby", sp + ".gem.clust.labz"});

      // ------------------------------------------------------------- unpack
      // Denominator with NO vertex requirement: rebuilt from the hodoscope branches
      // alone so the two denominators can be compared.
      df = df.Define(sp + "_tof_all",
                     [](const RVd &pl1, const RVd &pd1, const RVd &yp, const RVd &t1, const RVd &ip1) {
                       RVd r;
                       for (size_t i = 0; i < pl1.size(); ++i) {
                         if (ip1[i] != 1.)
                           continue;
                         const int pi = (int)std::round(pl1[i]);
                         if (pi != 1 && pi != 3)
                           continue;
                         if (!hodo_paddle_ok(pd1[i]))
                           continue;
                         r.push_back(t1[i] - hodo_path_len(pi, pd1[i], yp[i]) / 100. / 0.3);
                       }
                       return r;
                     },
                     {sp + "_plane_1", sp + "_paddle_1", sp + "_ypos_1", sp + "_tof_1", sp + "_isProton_1"});
      // Denominator WITH a valid vertex (stride 2 + N_VAR, offset 0 of each hit).
      df = df.Define(sp + "_tof_den",
                     [](const RVd &p) {
                       RVd r;
                       for (size_t h = 0; h + 2 + N_VAR <= p.size(); h += 2 + N_VAR)
                         r.push_back(p[h]);
                       return r;
                     },
                     {sp + "_epack"});
      // Numerators: one per (cut, variant).
      for (int ic = 0; ic < N_CUTS; ++ic) {
        const double cut = ADC_CUT_BASE * ADC_CUT_SCALES[ic];
        for (int iv = 0; iv < N_VAR; ++iv) {
          df = df.Define(sp + "_tof_trk_c" + std::to_string(ic) + "_v" + std::to_string(iv),
                         [cut, iv](const RVd &p) {
                           RVd r;
                           for (size_t h = 0; h + 2 + N_VAR <= p.size(); h += 2 + N_VAR) {
                             const double a = p[h + 2 + iv];
                             if (a >= 0. && a >= cut)
                               r.push_back(p[h]);
                           }
                           return r;
                         },
                         {sp + "_epack"});
        }
      }
      // adc_sum/dof spectra, per variant x tof region x hodoscope plane group.
      for (int iv = 0; iv < N_VAR; ++iv)
        for (int rg = 0; rg < N_GREG; ++rg) {
          const auto ivals = GREG_INT[rg];
          for (int gp = 0; gp < N_GGRP; ++gp) {
            df = df.Define(sp + "_adc_v" + std::to_string(iv) + "_r" + std::to_string(rg) + "_g" +
                               std::to_string(gp),
                           [ivals, gp, iv](const RVd &p) {
                             RVd r;
                             for (size_t h = 0; h + 2 + N_VAR <= p.size(); h += 2 + N_VAR) {
                               if ((int)std::round(p[h + 1]) != gp)
                                 continue;
                               bool in = false;
                               for (const auto &v : ivals)
                                 if (p[h] >= v[0] && p[h] < v[1])
                                   in = true;
                               if (!in)
                                 continue;
                               const double a = p[h + 2 + iv];
                               if (a >= 0.)
                                 r.push_back(a);
                             }
                             return r;
                           },
                           {sp + "_epack"});
          }
        }
      // residual spectra, per slot x tof region
      for (int s = 0; s < N_SLOT; ++s)
        for (int rg = 0; rg < N_GREG; ++rg) {
          const auto ivals = GREG_INT[rg];
          df = df.Define(sp + "_res_s" + std::to_string(s) + "_r" + std::to_string(rg),
                         [ivals, s](const RVd &p) {
                           RVd r;
                           for (size_t h = 0; h + 4 <= p.size(); h += 4) {
                             if ((int)std::round(p[h + 1]) != s)
                               continue;
                             bool in = false;
                             for (const auto &v : ivals)
                               if (p[h] >= v[0] && p[h] < v[1])
                                 in = true;
                             if (in)
                               r.push_back(p[h + 2]);
                           }
                           return r;
                         },
                         {sp + "_rpack"});
        }
    }
  } // !load

  // ---------------------------------------------------------------
  // 4. Book every histogram (fill or load)
  // ---------------------------------------------------------------
  for (int is = 0; is < N_SPECS; ++is) {
    if (!spec_ok[is])
      continue;
    const std::string sp(1, specs[is]);
    BK(h_tof_tot[is], sp + "_tof_all", sp + "_tof_all", sp + " proton, all;tof-L/c (ns);Counts", NBINS_TCORR,
       XMIN_TCORR, XMAX_TCORR);
    BK(h_tof_den[is], sp + "_tof_den", sp + "_tof_den", sp + " proton, vertex ok;tof-L/c (ns);Counts", NBINS_TCORR,
       XMIN_TCORR, XMAX_TCORR);
    for (int ic = 0; ic < N_CUTS; ++ic)
      for (int iv = 0; iv < N_VAR; ++iv) {
        const std::string col = sp + "_tof_trk_c" + std::to_string(ic) + "_v" + std::to_string(iv);
        BK(h_tof_trk[is][ic][iv], col, col,
           sp + " proton + " + VARIANTS[iv].dir + ";tof-L/c (ns);Counts", NBINS_TCORR, XMIN_TCORR, XMAX_TCORR);
      }
    for (int iv = 0; iv < N_VAR; ++iv)
      for (int rg = 0; rg < N_GREG; ++rg)
        for (int gp = 0; gp < N_GGRP; ++gp) {
          const std::string col =
              sp + "_adc_v" + std::to_string(iv) + "_r" + std::to_string(rg) + "_g" + std::to_string(gp);
          BK(h_adc[is][iv][rg][gp], col, col,
             sp + " " + VARIANTS[iv].dir + " " + GREG_NAME[rg] + " " + GGRP_TITLE[gp] +
                 ";adc_sum / dof;Counts",
             ADC_NBINS, ADC_LO, ADC_HI);
        }
    for (int s = 0; s < N_SLOT; ++s)
      for (int rg = 0; rg < N_GREG; ++rg) {
        const std::string col = sp + "_res_s" + std::to_string(s) + "_r" + std::to_string(rg);
        BK(h_res[is][s][rg], col, col,
           sp + " " + SLOT_NAME[s] + " " + GREG_NAME[rg] + ";residual along strip axis (cm);Clusters", RES_NBINS,
           RES_LO, RES_HI);
      }
  }

  // ---------------------------------------------------------------
  // 5. Event loop
  // ---------------------------------------------------------------
  if (!load) {
    std::cout << "[edep_tracking] Running event loop (" << bind1.size() << " histograms)...\n";
    for (auto &b : bind1)
      *b.slot = (TH1D *)b.res.GetPtr()->Clone();
    std::fprintf(stderr, "\n");
    if (cache_on) {
      TFile fc(cache_file, "RECREATE");
      if (!fc.IsZombie()) {
        TNamed sg("signature", sig.c_str());
        sg.Write();
        for (auto &b : bind1)
          if (*b.slot)
            (*b.slot)->Write();
        fc.Close();
        std::cout << "[edep_tracking] wrote histogram cache: " << cache_file << "\n";
      } else {
        std::cerr << "[edep_tracking] warning: could not write cache " << cache_file << "\n";
      }
    }
  } else {
    std::cout << "[edep_tracking] Histograms loaded from cache; event loop skipped.\n";
  }

  // ---------------------------------------------------------------
  // 6. Output
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
  // Shade the peak window on the current pad.
  auto draw_peak_lines = [](TH1D *h) {
    double y0 = h->GetMinimum(), y1 = h->GetMaximum();
    for (double x : {PEAK_WIN[0], PEAK_WIN[1]}) {
      TLine *l = new TLine(x, y0, x, y1);
      l->SetLineColor(kGreen + 2);
      l->SetLineStyle(2);
      l->Draw();
    }
  };

  for (int is = 0; is < N_SPECS; ++is) {
    if (!spec_ok[is])
      continue;
    const std::string sp(1, specs[is]);
    TDirectory *sdir = fout.mkdir(sp.c_str());

    // ---- residual diagnostics -------------------------------------------
    // One canvas per logical axis, 4 pads (tof regions), GEM0 and GEM1 overlaid.
    // These are the plots to compare against the widths the weighting assumes:
    // the peak sits on a flat combinatoric background from unrelated clusters.
    {
      TDirectory *d = mkpath(sdir, "residuals");
      d->cd();
      for (int ax = 0; ax < 2; ++ax) {
        TCanvas *c = new TCanvas((sp + "_c_res_" + (ax == AX_X ? "x" : "y")).c_str(),
                                 (sp + " residual along " + (ax == AX_X ? "V (x)" : "U (y)") + " strips").c_str(),
                                 1400, 1000);
        c->Divide(2, 2);
        for (int rg = 0; rg < N_GREG; ++rg) {
          c->cd(rg + 1);
          TLegend *lg = new TLegend(0.62, 0.72, 0.89, 0.89);
          lg->SetBorderSize(0);
          bool first = true;
          for (int lay = 0; lay < N_GEMLAY; ++lay) {
            int s = -1;
            for (int k = 0; k < N_SLOT; ++k)
              if (SLOT_LAYER[k] == lay && SLOT_AXIS[k] == ax)
                s = k;
            TH1D *h = h_res[is][s][rg];
            if (!h)
              continue;
            h->SetLineColor(lay == 0 ? kRed + 1 : kAzure + 2);
            h->SetLineWidth(2);
            h->SetTitle((sp + " " + (ax == AX_X ? "x (V)" : "y (U)") + " residual, " + GREG_NAME[rg] +
                         ";residual (cm);Clusters")
                            .c_str());
            h->Draw(first ? "HIST" : "HIST SAME");
            lg->AddEntry(h, (std::string("GEM") + std::to_string(lay) + " (window +/-" +
                             TString::Format("%.2f", weight_reach(lay, ax, 0)).Data() + " cm)")
                                .c_str(),
                         "l");
            first = false;
          }
          lg->Draw();
        }
        wc(c);
      }
      for (int s = 0; s < N_SLOT; ++s)
        for (int rg = 0; rg < N_GREG; ++rg)
          if (h_res[is][s][rg])
            h_res[is][s][rg]->Write();
    }

    // ---- adc_sum/dof spectra --------------------------------------------
    for (int iv = 0; iv < N_VAR; ++iv) {
      TDirectory *d = mkpath(sdir, std::string("adc_sum/") + VARIANTS[iv].dir);
      d->cd();
      TCanvas *c = new TCanvas((sp + "_c_adc_" + VARIANTS[iv].dir).c_str(),
                               (sp + " " + VARIANTS[iv].dir + " adc_sum/dof (dof = " +
                                std::to_string(VARIANTS[iv].nslot) + ")")
                                   .c_str(),
                               1400, 1000);
      c->Divide(2, 2);
      for (int rg = 0; rg < N_GREG; ++rg) {
        c->cd(rg + 1);
        gPad->SetLogy();
        double ymax = 1.;
        for (int gp = 0; gp < N_GGRP; ++gp)
          if (h_adc[is][iv][rg][gp])
            ymax = std::max(ymax, h_adc[is][iv][rg][gp]->GetMaximum());
        bool first = true;
        TLegend *lg = new TLegend(0.60, 0.70, 0.89, 0.89);
        lg->SetBorderSize(0);
        for (int gp = 0; gp < N_GGRP; ++gp) {
          TH1D *h = h_adc[is][iv][rg][gp];
          if (!h)
            continue;
          h->SetLineColor(GGRP_COLOR[gp]);
          h->SetLineWidth(2);
          h->SetTitle((sp + " " + VARIANTS[iv].dir + " " + GREG_NAME[rg] + ";adc_sum / dof;Counts").c_str());
          h->GetYaxis()->SetRangeUser(0.5, ymax * 2.);
          h->Draw(first ? "HIST" : "HIST SAME");
          lg->AddEntry(h, (std::string("plane ") + GGRP_TITLE[gp]).c_str(), "l");
          first = false;
        }
        for (int ic = 0; ic < N_CUTS; ++ic) {
          const double x = ADC_CUT_BASE * ADC_CUT_SCALES[ic];
          TLine *l = new TLine(x, 0.5, x, ymax * 2.);
          l->SetLineColor(kGreen + 2);
          l->SetLineStyle(ic == 1 ? 1 : 2);
          l->Draw();
        }
        lg->Draw();
      }
      wc(c);
      for (int rg = 0; rg < N_GREG; ++rg)
        for (int gp = 0; gp < N_GGRP; ++gp)
          if (h_adc[is][iv][rg][gp])
            h_adc[is][iv][rg][gp]->Write();
    }

    // ---- efficiency ------------------------------------------------------
    // The two denominators are variant- and cut-independent, so they live once at
    // the top of proton_id/ rather than being copied into every variant folder.
    {
      TDirectory *d = mkpath(sdir, "proton_id");
      d->cd();
      if (h_tof_tot[is])
        h_tof_tot[is]->Write();
      if (h_tof_den[is])
        h_tof_den[is]->Write();
    }
    // Per-variant summary series, one per cut, filled inside the cut loop.
    std::vector<std::string> sum_names;
    std::array<std::vector<double>, N_CUTS> sum_eff, sum_eff_err; // peak avg of the bg-sub ratio
    std::array<std::vector<double>, N_CUTS> sum_bkg, sum_bkg_err; // sideband avg of the raw ratio

    for (int ic = 0; ic < N_CUTS; ++ic) {
      const double cutval = ADC_CUT_BASE * ADC_CUT_SCALES[ic];
      char cbuf[64];
      std::snprintf(cbuf, sizeof(cbuf), "%g", cutval);
      for (int iv = 0; iv < N_VAR; ++iv) {
        TDirectory *d =
            mkpath(sdir, std::string("proton_id/adccut_") + cbuf + "/" + VARIANTS[iv].dir);
        d->cd();
        TH1D *hden = h_tof_den[is];
        TH1D *htrk = h_tof_trk[is][ic][iv];
        if (!hden || !htrk)
          continue;

        TCanvas *c = new TCanvas((sp + "_c_proton_tof").c_str(),
                                 (sp + " " + VARIANTS[iv].dir + "  adc_sum/dof >= " + cbuf).c_str(), 1400, 1000);
        c->Divide(2, 2);

        // pad 1: denominator (proton hits with a valid vertex), all-hits overlaid
        c->cd(1);
        hden->SetLineColor(kBlack);
        hden->SetLineWidth(2);
        hden->SetTitle((sp + " proton total;tof-L/c (ns);Counts").c_str());
        hden->Draw("HIST");
        if (h_tof_tot[is]) {
          h_tof_tot[is]->SetLineColor(kGray + 2);
          h_tof_tot[is]->SetLineStyle(2);
          h_tof_tot[is]->Draw("HIST SAME");
        }
        draw_peak_lines(hden);
        {
          TLegend *lg = new TLegend(0.55, 0.72, 0.89, 0.89);
          lg->SetBorderSize(0);
          lg->AddEntry(hden, "vertex ok (denominator)", "l");
          if (h_tof_tot[is])
            lg->AddEntry(h_tof_tot[is], "all proton hits", "l");
          lg->Draw();
        }

        // pad 2: numerator
        c->cd(2);
        htrk->SetLineColor(kBlue + 1);
        htrk->SetLineWidth(2);
        htrk->SetTitle((sp + " proton + " + VARIANTS[iv].dir + " (adc/dof >= " + cbuf + ");tof-L/c (ns);Counts")
                           .c_str());
        htrk->Draw("HIST");
        draw_peak_lines(htrk);

        // pad 3: sideband-subtracted track/total ratio (rebinned to 5 ns)
        c->cd(3);
        TH1D *ht_rb = (TH1D *)htrk->Clone((sp + "_rb_trk").c_str());
        TH1D *hd_rb = (TH1D *)hden->Clone((sp + "_rb_tot").c_str());
        ht_rb->Rebin(PROTON_REBIN);
        hd_rb->Rebin(PROTON_REBIN);
        TH1D *ht_sb = flat_bgsub2(ht_rb, SB_LO1, SB_HI1, SB_LO2, SB_HI2);
        TH1D *hd_sb = flat_bgsub2(hd_rb, SB_LO1, SB_HI1, SB_LO2, SB_HI2);
        TH1D *hratio = (TH1D *)ht_sb->Clone((sp + "_ratio_" + VARIANTS[iv].dir).c_str());
        hratio->SetTitle((sp + " (track-bg)/(total-bg);tof-L/c (ns);ratio").c_str());
        hratio->Divide(hd_sb);
        double eff_m = 0., eff_e = 0.;
        region_wstats(hratio, ht_rb, hd_rb, {{PEAK_WIN[0], PEAK_WIN[1]}}, eff_m, eff_e);
        hratio->GetXaxis()->SetRangeUser(-50., 175.);
        hratio->GetYaxis()->SetRangeUser(0., 1.5);
        hratio->SetLineColor(kBlue + 1);
        hratio->SetLineWidth(2);
        hratio->Draw("HIST");
        draw_peak_lines(hratio);
        {
          TLatex tx;
          tx.SetNDC();
          tx.SetTextSize(0.045);
          tx.DrawLatex(0.15, 0.85, TString::Format("peak eff = %.3f #pm %.3f", eff_m, eff_e));
        }

        // pad 4: raw ratio (no background subtraction). Its sideband average is the
        // fraction of random background surviving the cut.
        c->cd(4);
        TH1D *ht_rw = (TH1D *)ht_rb->Clone((sp + "_raw_trk").c_str());
        TH1D *hd_rw = (TH1D *)hd_rb->Clone((sp + "_raw_tot").c_str());
        TH1D *hraw = (TH1D *)ht_rw->Clone((sp + "_ratio_raw_" + VARIANTS[iv].dir).c_str());
        hraw->SetTitle((sp + " track/total, no bg-sub;tof-L/c (ns);ratio").c_str());
        hraw->Divide(hd_rw);
        double bkg_m = 0., bkg_e = 0.;
        region_wstats(hraw, ht_rw, hd_rw, {{SB_LO1, SB_HI1}, {SB_LO2, SB_HI2}}, bkg_m, bkg_e);
        hraw->GetXaxis()->SetRangeUser(-150., 175.);
        hraw->GetYaxis()->SetRangeUser(0., 1.5);
        hraw->SetLineColor(kMagenta + 2);
        hraw->SetLineWidth(2);
        hraw->Draw("HIST");
        {
          TLatex tx;
          tx.SetNDC();
          tx.SetTextSize(0.045);
          tx.DrawLatex(0.15, 0.85, TString::Format("sideband = %.3f #pm %.3f", bkg_m, bkg_e));
        }

        htrk->Write();
        hratio->Write();
        hraw->Write();
        wc(c);

        if (ic == 0)
          sum_names.push_back(VARIANTS[iv].dir);
        sum_eff[ic].push_back(eff_m);
        sum_eff_err[ic].push_back(eff_e);
        sum_bkg[ic].push_back(bkg_m);
        sum_bkg_err[ic].push_back(bkg_e);

        delete ht_rb;
        delete hd_rb;
        delete ht_sb;
        delete hd_sb;
        delete hratio;
        delete ht_rw;
        delete hd_rw;
        delete hraw;
      }
    }

    // ---- summary: efficiency / background vs variant, one series per cut ----
    {
      TDirectory *d = mkpath(sdir, "summary");
      d->cd();
      const int cols[N_CUTS] = {kGreen + 2, kBlue + 1, kRed + 1};
      auto mk_summary = [&](const char *nm, const char *ttl, const char *yttl,
                            const std::array<std::vector<double>, N_CUTS> &vals,
                            const std::array<std::vector<double>, N_CUTS> &errs, double ymax) {
        if (sum_names.empty())
          return;
        const int n = (int)sum_names.size();
        TCanvas *c = new TCanvas((sp + "_c_" + nm).c_str(), (sp + " " + ttl).c_str(), 1400, 800);
        c->SetGridy();
        c->SetBottomMargin(0.22);
        TH1D *frame = new TH1D((sp + "_frame_" + nm).c_str(), (sp + " " + ttl + ";;" + yttl).c_str(), n, 0., n);
        for (int i = 0; i < n; ++i)
          frame->GetXaxis()->SetBinLabel(i + 1, sum_names[i].c_str());
        frame->GetXaxis()->LabelsOption("v");
        frame->GetYaxis()->SetRangeUser(0., ymax);
        frame->Draw();
        TLegend *lg = new TLegend(0.70, 0.74, 0.89, 0.89);
        lg->SetBorderSize(0);
        for (int ic = 0; ic < N_CUTS; ++ic) {
          if ((int)vals[ic].size() != n)
            continue;
          auto *g = new TGraphErrors(n);
          for (int i = 0; i < n; ++i) {
            g->SetPoint(i, i + 0.5, vals[ic][i]);
            g->SetPointError(i, 0., errs[ic][i]);
          }
          g->SetMarkerStyle(20 + ic);
          g->SetMarkerColor(cols[ic]);
          g->SetLineColor(cols[ic]);
          g->SetMarkerSize(1.3);
          g->Draw("P SAME");
          lg->AddEntry(g, TString::Format("adc/dof >= %g", ADC_CUT_BASE * ADC_CUT_SCALES[ic]), "p");
        }
        lg->Draw();
        wc(c);
        delete frame;
      };
      mk_summary("eff_vs_variant", "peak-region tracking efficiency", "efficiency", sum_eff, sum_eff_err, 1.2);
      mk_summary("bkg_vs_variant", "sideband survival (random background)", "surviving fraction", sum_bkg,
                 sum_bkg_err, 1.2);
    }
  }

  fout.cd();
  TNamed sg("signature", sig.c_str());
  sg.Write();
  fout.Close();
  if (fcache) {
    fcache->Close();
    delete fcache;
  }
  const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count();
  std::cout << "[edep_tracking] wrote " << out_file << " in " << secs << " s\n";
}
