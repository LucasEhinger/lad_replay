// Tracking-efficiency diagnostic skim: one row per SHMS (P) proton hit, with the
// GEM geometry seen by that hit, so the efficiency can be studied against any
// quantity without re-reading the replay files.
//
//   root -l -b -q 'skim_trk.C+("files.txt","out.root")'
//
// Selection (same as lad_tracking_eff.C): P.react.ok != 0, |P.react.z| < 20,
// goodhit isProton_1 == 1, plane_1 in {1 (001), 3 (101)}, tof - L/c in [-150, 175].
//
// Per hit and per GEM slot s (0 = x0 [GEM0 V], 1 = y0 [GEM0 U], 2 = x1 [GEM1 V],
// 3 = y1 [GEM1 U]):
//   pred[s]  where the straight line vertex -> back hodo point crosses the module,
//            along that slot's strip axis (cm; the road centre in Do1DClusterTracking)
//   res[s]   nearest cluster - pred (cm), 999 if the slot has no cluster
//   res2[s]  second-nearest
//   adcn/tn/nsn[s]  ADC sum, time, nstrip of that nearest cluster
//   ncl[s]   clusters in the slot;  n10[s] clusters within 10 cm
//   rw[s]    winning (goodhit_trk1D_clidx) cluster - pred, 999 if none
// GEM plane geometry (origin, strip axes) is the 13.5 deg alignment, fitted
// exactly from clust.lab{x,y,z} = origin + pos * axisHat.
#include <ROOT/RVec.hxx>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <TVector3.h>
#include <TMath.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace sk {
// slot -> (layer, axis): axis 1 = V (x), 0 = U (y)
const int SL_LAYER[4] = {0, 0, 1, 1};
const int SL_AXIS[4] = {1, 0, 1, 0};
const TVector3 ORIG[2] = {TVector3(68.2415, 1.3775, -50.5867), TVector3(83.6552, 1.3082, -61.9537)};
const TVector3 AXU[2] = {TVector3(0.015691, 0.999781, -0.013835), TVector3(0.018632, 0.99971, -0.015249)};
const TVector3 AXV[2] = {TVector3(0.609729, 0.001398, 0.792609), TVector3(0.607471, 0.000795, 0.794341)};
// hodo planes 000 001 100 101 200
const double HZ[5] = {618.8125006, 658.5516251, 526.1277344, 568.079082, 614.5200935};
const double HTH[5] = {2.612452591, 2.615772334, 2.212325822, 2.213860777, 1.814590772};
const double HR[5] = {615., 655.6, 523., 563.6, 615.}; // radii used for tof - L/c in the macros
const char *VARS[9] = {"1D_xz_GEM0", "1D_xz_GEM1", "1D_xz_GEMboth", "1D_y_GEM0",   "1D_y_GEM1",
                       "1D_y_GEMboth", "1D_GEM0",  "1D_GEM1",      "1D_GEMboth"};

TVector3 hodo_lab(int plane, int paddle, double ypos) {
  TVector3 h(110. - 22. * paddle, ypos, HZ[plane]);
  h.RotateY(HTH[plane]);
  return h;
}
} // namespace sk

// allhits = 1: keep every back-plane (001/101) good hit, not only isProton_1 == 1, and skip the
// per-cluster vectors (proton-ID study).
void skim_trk(const char *list, const char *out, Long64_t maxev = -1, int allhits = 0) {
  using namespace sk;
  TChain ch("T");
  std::ifstream in(list);
  std::string fn;
  while (in >> fn)
    ch.Add(fn.c_str());
  ch.SetBranchStatus("*", 0);
  for (const char *b : {"P.react.*", "P.ladhod.goodhit_*", "Ndata.P.ladhod.goodhit_*", "P.gem.clust.*",
                        "Ndata.P.gem.clust.*", "g.evnum", "g.runnum", "g.evtyp", "g.evtime",
                        "T.shms.pTRIG1_tdcTime", "T.shms.hTRIG1_tdcTime", "T.shms.shmsTrigLAD_tdcTime",
                        "T.shms.hmsTrigLAD_tdcTime"})
    ch.SetBranchStatus(b, 1);

  TTreeReader r(&ch);
  if (maxev > 0)
    r.SetEntriesRange(0, maxev);
  TTreeReaderValue<double> ok(r, "P.react.ok"), vx(r, "P.react.x"), vy(r, "P.react.y"), vz(r, "P.react.z");
  TTreeReaderValue<double> evn(r, "g.evnum"), runn(r, "g.runnum"), evtyp(r, "g.evtyp"), evtime(r, "g.evtime");
  TTreeReaderValue<double> tp1(r, "T.shms.pTRIG1_tdcTime"), th1(r, "T.shms.hTRIG1_tdcTime"),
      tsl(r, "T.shms.shmsTrigLAD_tdcTime"), thl(r, "T.shms.hmsTrigLAD_tdcTime");
  auto A = [&](const char *n) { return TTreeReaderArray<double>(r, Form("P.ladhod.goodhit_%s", n)); };
  TTreeReaderArray<double> pl0 = A("plane_0"), pl1 = A("plane_1"), pd0 = A("paddle_0"), pd1 = A("paddle_1"),
                           ip1 = A("isProton_1"), ip0 = A("isProton_0"), yp0 = A("hit_ypos_0"),
                           yp1 = A("hit_ypos_1"), tof1 = A("hit_tof_1"), ht0 = A("hittime_0"), ht1 = A("hittime_1"),
                           e0 = A("hitedep_MeV_0"), e1 = A("hitedep_MeV_1");
  std::vector<TTreeReaderArray<double>> chi;
  chi.reserve(9);
  for (auto v : VARS)
    chi.emplace_back(r, Form("P.ladhod.goodhit_chiSquare_%s", v));
  TTreeReaderArray<double> cx0 = A("trk1D_clidx_x0"), cy0 = A("trk1D_clidx_y0"), cx1 = A("trk1D_clidx_x1"),
                           cy1 = A("trk1D_clidx_y1");
  TTreeReaderArray<double> clpos(r, "P.gem.clust.pos"), cllay(r, "P.gem.clust.layer"), clax(r, "P.gem.clust.axis"),
      cladc(r, "P.gem.clust.adc"), cltime(r, "P.gem.clust.time"), clidx(r, "P.gem.clust.index"),
      clns(r, "P.gem.clust.nstrip"), clmpd(r, "P.gem.clust.mpd"), clapv(r, "P.gem.clust.apv");

  TFile fo(out, "RECREATE");
  TTree t("H", "P proton hits");
  Int_t isp0, isp1, run, etyp, mpdn[6], apvmask[6], nprot, ngood, ipl0, ipd0, ipl1, ipd1, ncl[4], n10[4];
  Long64_t evt;
  Float_t dtprev, dtprev1, trg[4], tofc, y0, y1, ed0, ed1, dt, x_v, y_v, z_v, c[9], pred[4], res[4], res2[4], adcn[4], tn[4], nsn[4], rw[4];
  t.Branch("run", &run);
  t.Branch("evt", &evt);
  t.Branch("evtyp", &etyp);
  t.Branch("dtprev", &dtprev);   // g.evtime ticks since the previous event (any type)
  t.Branch("dtprev1", &dtprev1); // ticks since the previous evtyp-1 (SHMS, GEM read) event
  t.Branch("mpdn", mpdn, "mpdn[6]/I");       // clusters per MPD in the event
  t.Branch("apvmask", apvmask, "apvmask[6]/I"); // bit a set if APV a of that MPD has a cluster
  t.Branch("trg", trg, "trg[4]/F"); // pTRIG1, hTRIG1, shmsTrigLAD, hmsTrigLAD tdc times
  t.Branch("nprot", &nprot);
  t.Branch("isp0", &isp0);
  t.Branch("isp1", &isp1);
  t.Branch("ngood", &ngood);
  t.Branch("pl0", &ipl0);
  t.Branch("pd0", &ipd0);
  t.Branch("pl1", &ipl1);
  t.Branch("pd1", &ipd1);
  t.Branch("tofc", &tofc);
  t.Branch("y0", &y0);
  t.Branch("y1", &y1);
  t.Branch("ed0", &ed0);
  t.Branch("ed1", &ed1);
  t.Branch("dt", &dt);
  t.Branch("vx", &x_v);
  t.Branch("vy", &y_v);
  t.Branch("vz", &z_v);
  t.Branch("chi", c, "chi[9]/F");
  t.Branch("pred", pred, "pred[4]/F");
  t.Branch("res", res, "res[4]/F");
  t.Branch("res2", res2, "res2[4]/F");
  t.Branch("adcn", adcn, "adcn[4]/F");
  t.Branch("tn", tn, "tn[4]/F");
  t.Branch("nsn", nsn, "nsn[4]/F");
  t.Branch("rw", rw, "rw[4]/F");
  t.Branch("ncl", ncl, "ncl[4]/I");
  t.Branch("n10", n10, "n10[4]/I");

  std::vector<int> cs; // slot of each cluster
  std::vector<float> cr, ct, ca;
  t.Branch("cs", &cs);
  t.Branch("cr", &cr); // cluster pos - this hit's prediction (cm)
  t.Branch("ct", &ct); // cluster time
  t.Branch("ca", &ca); // cluster ADC sum
  Long64_t nev = 0;
  double last_t = -1, last_t1 = -1;
  while (r.Next()) {
    const double now = *evtime;
    dtprev = last_t > 0 ? now - last_t : -1;
    dtprev1 = last_t1 > 0 ? now - last_t1 : -1;
    last_t = now;
    if ((int)*evtyp == 1)
      last_t1 = now;
    if (++nev % 200000 == 0)
      std::cout << nev << " events, " << t.GetEntries() << " hits" << std::endl;
    if (*ok == 0 || std::fabs(*vz) >= 20.)
      continue;
    const size_t nh = pl1.GetSize();
    std::vector<size_t> sel;
    for (size_t i = 0; i < nh; ++i) {
      int p = (int)std::round(pl1[i]);
      if ((allhits || ip1[i] == 1.) && (p == 1 || p == 3))
        sel.push_back(i);
    }
    if (sel.empty())
      continue;
    const TVector3 vtx(*vx, *vy, *vz);
    const size_t nc = clpos.GetSize();
    for (int m = 0; m < 6; ++m)
      mpdn[m] = apvmask[m] = 0;
    for (size_t k = 0; k < nc; ++k) {
      int m = (int)clmpd[k], a = (int)clapv[k];
      if (m >= 0 && m < 6) {
        ++mpdn[m];
        if (a >= 0 && a < 31)
          apvmask[m] |= (1 << a);
      }
    }
    for (size_t i : sel) {
      int p = (int)std::round(pl1[i]);
      double dx = 110. - 22. * pd1[i];
      double tc = tof1[i] - std::sqrt(yp1[i] * yp1[i] + dx * dx + HR[p] * HR[p]) / 30.;
      if (!(tc >= -150. && tc <= 175.))
        continue;
      run = (int)*runn;
      evt = (Long64_t)*evn;
      etyp = (int)*evtyp;
      trg[0] = *tp1;
      trg[1] = *th1;
      trg[2] = *tsl;
      trg[3] = *thl;
      nprot = (int)sel.size();
      ngood = (int)nh;
      isp0 = (int)ip0[i];
      isp1 = (int)ip1[i];
      ipl0 = (int)std::round(pl0[i]);
      ipd0 = (int)std::round(pd0[i]);
      ipl1 = p;
      ipd1 = (int)std::round(pd1[i]);
      tofc = tc;
      y0 = yp0[i];
      y1 = yp1[i];
      ed0 = e0[i];
      ed1 = e1[i];
      dt = ht1[i] - ht0[i];
      x_v = *vx;
      y_v = *vy;
      z_v = *vz;
      for (int k = 0; k < 9; ++k)
        c[k] = chi[k][i];
      const TVector3 hb = hodo_lab(p, ipd1, yp1[i]);
      const TVector3 d = hb - vtx;
      const double win[4] = {cx0[i], cy0[i], cx1[i], cy1[i]};
      cs.clear();
      cr.clear();
      ct.clear();
      ca.clear();
      for (int s = 0; s < 4; ++s) {
        const int L = SL_LAYER[s], AX = SL_AXIS[s];
        const TVector3 nrm = AXU[L].Cross(AXV[L]).Unit();
        const TVector3 &ah = AX == 1 ? AXV[L] : AXU[L];
        const TVector3 P = vtx + (nrm.Dot(ORIG[L] - vtx) / nrm.Dot(d)) * d;
        pred[s] = ah.Dot(P - ORIG[L]);
        res[s] = res2[s] = rw[s] = 999.;
        adcn[s] = tn[s] = nsn[s] = -1.;
        ncl[s] = n10[s] = 0;
        int best = -1;
        for (size_t k = 0; k < nc; ++k) {
          if ((int)cllay[k] != L || (int)clax[k] != AX)
            continue;
          ++ncl[s];
          double rr = clpos[k] - pred[s];
          if (!allhits) {
          cs.push_back(s);
          cr.push_back(rr);
          ct.push_back(cltime[k]);
          ca.push_back(cladc[k]);
          }
          if (std::fabs(rr) < 10.)
            ++n10[s];
          if (win[s] >= 0 && (int)clidx[k] == (int)win[s])
            rw[s] = rr;
          if (std::fabs(rr) < std::fabs(res[s])) {
            res2[s] = res[s];
            res[s] = rr;
            best = (int)k;
          } else if (std::fabs(rr) < std::fabs(res2[s]))
            res2[s] = rr;
        }
        if (best >= 0) {
          adcn[s] = cladc[best];
          tn[s] = cltime[best];
          nsn[s] = clns[best];
        }
      }
      t.Fill();
    }
  }
  std::cout << "done: " << nev << " events, " << t.GetEntries() << " hits" << std::endl;
  t.Write();
}
