#include <TFile.h>
#include <TH2D.h>
#include <TMinuit.h>
#include <TTree.h>
#include <TVector3.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

struct fit_range {
  double min;
  double max;
  int nsteps;

  fit_range(double min_val, double max_val, int steps) : min(min_val), max(max_val), nsteps(steps) {}
};

struct fit_param {
  double gem0_dx;
  double gem0_dy;
  double gem1_dx;
  double gem1_dy;
  double gem0_theta;
  double gem1_theta;
  double gem0_z;
  double gem1_z;
};

struct good_track {
  double trk_x[2];
  double trk_y[2];
  TVector3 p_hit;
  TVector3 p_react;
};

const fit_range gem0_x_range(-10.0, 10.0, 100);
// const fit_range gem0_y_range(-50.0, 50.0, 10);
const fit_range gem0_z_range(65, 90, 10);
const fit_range gem1_x_range(-10.0, 10.0, 100);
// const fit_range gem1_y_range(-50.0, 50.0, 10);
const fit_range gem1_z_range(85, 115.0, 10);
const fit_range gem_theta_range(127.0, 127.01, 1);
int run_number = 300000;

const int MAX_DATA  = 10000;
const int maxTracks = 30;
const int nPlanes   = 5;

const double d0_cut               = 30.0;
const double delta_trans_cut[2]   = {-50, 50.0};
const double delta_long_cut[2]    = {-50.0, 50.0};
const double plane_theta[nPlanes] = {150.0, 150.0, 127.0, 127.0, 104.0}; // Angle in degrees
const double plane_r[nPlanes]     = {615.0, 655.6, 523.0, 563.6, 615.0}; // Radius of the second point
std::vector<good_track> goodTracks;

// Should ultimately get these from the param file
// const double gem_theta = 127.0;            // Angle in degrees
// const double gem_phi   = 0.0;              // Angle in degrees
// const double gem_r[2]  = {77.571, 95.571}; // Radius of the GEM's

const int nFixedz = 3;
// const double fixedz[nFixedz] = {0.0}; // Fixed z positions for the planes
const double fixedz[nFixedz] = {-10.0, 0.0, 10.0}; // Fixed z positions for the planes
// const double z_range[nFixedz + 1] = {-15.0, -5.0, 5.0, 15.0}; // Tolerance range for the fixed z positions (low,
// high)
const double z_range[nFixedz + 1] = {-20.0, -5, 5.0, 15.0}; // Tolerance range for the fixed z positions (low,
// high)
// const double z_range[nFixedz + 1] = {-20.0, 20.0}; // Tolerance range for the fixed z positions (low, high)

const bool use_projz = true; // Fix z position to GEM projz

const double dx_min = -30.0;
const double dx_max = 30.0;
const int dx_NBINS  = 60;
const double dy_min = -30.0;
const double dy_max = 30.0;
const int dy_NBINS  = 60;

TVector3 LinePlaneIntersection(const TVector3 &p1, const TVector3 &p2, const TVector3 &p3, const TVector3 &l1,
                               const TVector3 &l2) {
  // Define the plane normal
  TVector3 planeNormal = (p2 - p1).Cross(p3 - p1);
  planeNormal          = planeNormal.Unit();

  // Line direction
  TVector3 lineDir = l2 - l1;

  // Check if the line is parallel to the plane
  double denom = planeNormal.Dot(lineDir);
  if (fabs(denom) < 1e-6) {
    std::cerr << "The line is parallel to the plane and does not intersect." << std::endl;
    return TVector3(0, 0, 0); // Return a zero vector if the line is parallel to the plane
  }

  // Calculate the intersection point
  double t              = planeNormal.Dot(p1 - l1) / denom;
  TVector3 intersection = l1 + t * lineDir;

  return intersection;
}

TVector3 GetHodoHitPosition(const int paddle, const int plane) {

  double dTrans = (5 - paddle) * 22.0;
  // Define the plane using three points
  TVector3 p_hit(plane_r[plane] * cos((plane_theta[plane] - 90) * TMath::DegToRad()), 0,
                 -plane_r[plane] * sin((plane_theta[plane] - 90) * TMath::DegToRad()));
  p_hit = p_hit + TVector3(-dTrans * cos((180 - plane_theta[plane]) * TMath::DegToRad()), 0,
                           -dTrans * sin((180 - plane_theta[plane]) * TMath::DegToRad()));

  return p_hit;
}

TVector3 GetGEMHitPosition(const double x_loc, const double y_loc, const double dx, const double dy, const double r,
                           const double gem_theta) {
  // Define the plane using three points
  TVector3 p_hit(r * cos((gem_theta - 90) * TMath::DegToRad()), 0, -r * sin((gem_theta - 90) * TMath::DegToRad()));
  p_hit = p_hit + TVector3(-(x_loc + dx) * cos((180 - gem_theta) * TMath::DegToRad()), dy,
                           -(x_loc + dx) * sin((180 - gem_theta) * TMath::DegToRad()));

  return p_hit;
}

double GetProjZ(const TVector3 &gem_hit0, const TVector3 &gem_hit1, const int x_targ) {

  // Calculate the direction vector of the line
  TVector3 direction = gem_hit1 - gem_hit0;

  // Check if the direction vector is valid
  if (fabs(direction.X()) < 1e-6) {
    std::cerr << "Warning: The line is parallel to the x-axis and cannot be projected. Setting proj_z to 100."
              << std::endl;
    return 100.0;
  }

  // Calculate the parameter t for the line equation
  double t = (x_targ - gem_hit0.X()) / direction.X();

  // Calculate the z-coordinate at the projected point
  double proj_z = gem_hit0.Z() + t * direction.Z();

  return proj_z;
}

double CostFunction(const std::vector<good_track> &goodTracks, const fit_param &fitParams) {

  double cost = 0;
  int nPoints = 0;
  TVector3 p1[2], p2[2], p3[2];
  p1[0] = TVector3(fitParams.gem0_z * cos((fitParams.gem0_theta - 90) * TMath::DegToRad()), 0,
                   -fitParams.gem0_z * sin((fitParams.gem0_theta - 90) * TMath::DegToRad()));
  p1[1] = TVector3(fitParams.gem1_z * cos((fitParams.gem1_theta - 90) * TMath::DegToRad()), 0,
                   -fitParams.gem1_z * sin((fitParams.gem1_theta - 90) * TMath::DegToRad()));
  p2[0] = p1[0] + TVector3(-10.0 * cos((180 - fitParams.gem0_theta) * TMath::DegToRad()), fitParams.gem0_dy,
                           -10.0 * sin((180 - fitParams.gem0_theta) * TMath::DegToRad()));
  p2[1] = p1[1] + TVector3(-10.0 * cos((180 - fitParams.gem1_theta) * TMath::DegToRad()), fitParams.gem1_dy,
                           -10.0 * sin((180 - fitParams.gem1_theta) * TMath::DegToRad()));
  p3[0] = p1[0] + TVector3(0, 10.0, 0); // Arbitrary y value of 10.0
  p3[1] = p1[1] + TVector3(0, 10.0, 0); // Arbitrary y value of 10.0

  for (const auto &track : goodTracks) {
    TVector3 p_gem_hit0 = GetGEMHitPosition(track.trk_x[0], track.trk_y[0], fitParams.gem0_dx, fitParams.gem0_dy,
                                            fitParams.gem0_z, fitParams.gem0_theta);
    TVector3 p_gem_hit1 = GetGEMHitPosition(track.trk_x[1], track.trk_y[1], fitParams.gem1_dx, fitParams.gem1_dy,
                                            fitParams.gem1_z, fitParams.gem1_theta);
    TVector3 p_hit      = track.p_hit;
    TVector3 p_react    = track.p_react;
    double proj_z       = GetProjZ(p_gem_hit0, p_gem_hit1, p_react.X());
    double min_dist_sq  = std::numeric_limits<double>::max();
    for (int k = 0; k < nFixedz; ++k) {
      double dist_sq = (proj_z - fixedz[k]) * (proj_z - fixedz[k]);
      if (dist_sq < min_dist_sq) {
        min_dist_sq = dist_sq;
      }
    }
    if (min_dist_sq < 500) {
      nPoints++;
      cost += min_dist_sq*5;
    }

    // TVector3 p_gem0_intersection = LinePlaneIntersection(p1[0], p2[0], p3[0], track.p_hit, track.p_react);
    // TVector3 p_gem1_intersection = LinePlaneIntersection(p1[1], p2[1], p3[1], track.p_hit, track.p_react);
    // double xz_separation[2];
    // xz_separation[0] =
    //     sqrt(pow(p_gem_hit0.X() - p_gem0_intersection.X(), 2) + pow(p_gem_hit0.Z() - p_gem0_intersection.Z(), 2));
    // xz_separation[1] =
    //     sqrt(pow(p_gem_hit1.X() - p_gem1_intersection.X(), 2) + pow(p_gem_hit1.Z() - p_gem1_intersection.Z(), 2));
    // double x_diff       = p_gem_hit0.X() - p_gem0_intersection.X();
    // double y_separation = p_gem_hit0.Y() - p_gem0_intersection.Y();
    // double xz_separation_sq_sum =
    //     xz_separation[0] * xz_separation[0] + xz_separation[1] * xz_separation[1];
    // if (xz_separation_sq_sum < 500) {
    //   cost += (xz_separation_sq_sum*10);
    // }

    // Calculate the line between the two GEM hits
    TVector3 lineDir = p_gem_hit1 - p_gem_hit0;
    TVector3 linePoint = p_gem_hit0;

    // Vector from the line point to p_hit
    TVector3 vecToPoint = p_hit - linePoint;

    // Project vecToPoint onto the line direction
    double t = vecToPoint.Dot(lineDir) / lineDir.Mag2();

    // Closest point on the line to p_hit
    TVector3 closestPoint = linePoint + t * lineDir;

    // Calculate the distance of closest approach in the x-z plane
    double dca_xz = sqrt(pow(closestPoint.X() - p_hit.X(), 2) + pow(closestPoint.Z() - p_hit.Z(), 2));

    // Add the distance to the cost function if it is within a reasonable range
    if (dca_xz < 500.0) { // Example threshold
      // cost += dca_xz * dca_xz;
    }
    // cost += xz_separation[0] * xz_separation[0] + xz_separation[1] * xz_separation[1];
  }
  cost /= nPoints;
  return cost;
}

void iterate_fit() {
  // Define the plane using three points
  // TVector3 p1[2], p2[2], p3[2];
  // for (int i = 0; i < 2; ++i) {
  //   double theta = gem_theta;
  //   p1[i]        = TVector3(gem_r[i] * cos((theta - 90) * TMath::DegToRad()), 0,
  //                           -gem_r[i] * sin((theta - 90) * TMath::DegToRad()));
  //   p2[i]        = p1[i] + TVector3(-gem_r[i] * cos((180 - theta) * TMath::DegToRad()), 0,
  //                                   -gem_r[i] * sin((180 - theta) * TMath::DegToRad()));
  //   p3[i]        = p1[i] + TVector3(0, 10.0, 0); // Arbitrary y value of 10.0
  // }

  TString fileName = "/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/"
                     //  "LAD_COIN_22282_-1_inverted.root";
                     //  "LAD_COIN_22282_-1_500trks_good_timing.root";
                     "LAD_COIN_22282_300005.root";

  TString outputFileName = "gem_pos_fixing_22282_300006.root";
  // Open the ROOT file
  TFile *file = TFile::Open(fileName);
  if (!file || file->IsZombie()) {
    std::cerr << "Error: Cannot open the ROOT file!" << std::endl;
    return;
  }

  // Get the TTree
  TTree *T = dynamic_cast<TTree *>(file->Get("T"));
  if (!T) {
    std::cerr << "Error: Cannot find the TTree named 'T'!" << std::endl;
    file->Close();
    return;
  }
  // Define arrays to hold the data
  Double_t trk_d0[MAX_DATA];
  Double_t trk_d0_good[MAX_DATA];
  Double_t trk_projz[MAX_DATA];
  Double_t trk_x[2][MAX_DATA], trk_y[2][MAX_DATA], trk_z[2][MAX_DATA];
  Double_t trk_x_local[2][MAX_DATA], trk_y_local[2][MAX_DATA];
  Double_t kin_trackID_0[MAX_DATA], kin_trackID_1[MAX_DATA];
  Double_t kin_plane_0[MAX_DATA], kin_plane_1[MAX_DATA];
  Double_t kin_paddle_0[MAX_DATA], kin_paddle_1[MAX_DATA];
  Double_t kin_hittime_0[MAX_DATA], kin_hittime_1[MAX_DATA];
  Double_t kin_hittheta_0[MAX_DATA], kin_hittheta_1[MAX_DATA];
  Double_t kin_hitphi_0[MAX_DATA], kin_hitphi_1[MAX_DATA];
  Double_t kin_hitedep_0[MAX_DATA], kin_hitedep_1[MAX_DATA];
  Double_t kin_deltapostrans_0[MAX_DATA], kin_deltapostrans_1[MAX_DATA];
  Double_t kin_deltaposlong_0[MAX_DATA], kin_deltaposlong_1[MAX_DATA];
  Int_t nTracks, nGoodHits;
  Double_t vertex_x, vertex_y, vertex_z;

  T->SetBranchAddress("Ndata.H.gem.trk.d0", &nTracks);
  T->SetBranchAddress("Ndata.H.ladkin.goodhit_trackid_0", &nGoodHits);
  T->SetBranchAddress("H.gem.trk.d0", &trk_d0);
  T->SetBranchAddress("H.gem.trk.d0_good", &trk_d0_good);
  T->SetBranchAddress("H.gem.trk.projz", &trk_projz);
  T->SetBranchAddress("H.gem.trk.x1", &trk_x[0]);
  T->SetBranchAddress("H.gem.trk.y1", &trk_y[0]);
  T->SetBranchAddress("H.gem.trk.z1", &trk_z[0]);
  T->SetBranchAddress("H.gem.trk.x2", &trk_x[1]);
  T->SetBranchAddress("H.gem.trk.y2", &trk_y[1]);
  T->SetBranchAddress("H.gem.trk.z2", &trk_z[1]);
  T->SetBranchAddress("H.gem.trk.x1_local", &trk_x_local[0]);
  T->SetBranchAddress("H.gem.trk.y1_local", &trk_y_local[0]);
  T->SetBranchAddress("H.gem.trk.x2_local", &trk_x_local[1]);
  T->SetBranchAddress("H.gem.trk.y2_local", &trk_y_local[1]);
  T->SetBranchAddress("H.ladkin.goodhit_trackid_0", &kin_trackID_0);
  T->SetBranchAddress("H.ladkin.goodhit_trackid_1", &kin_trackID_1);
  T->SetBranchAddress("H.ladkin.goodhit_plane_0", &kin_plane_0);
  T->SetBranchAddress("H.ladkin.goodhit_plane_1", &kin_plane_1);
  T->SetBranchAddress("H.ladkin.goodhit_paddle_0", &kin_paddle_0);
  T->SetBranchAddress("H.ladkin.goodhit_paddle_1", &kin_paddle_1);
  T->SetBranchAddress("H.ladkin.goodhit_hittime_0", &kin_hittime_0);
  T->SetBranchAddress("H.ladkin.goodhit_hittime_1", &kin_hittime_1);
  T->SetBranchAddress("H.ladkin.goodhit_hittheta_0", &kin_hittheta_0);
  T->SetBranchAddress("H.ladkin.goodhit_hittheta_1", &kin_hittheta_1);
  T->SetBranchAddress("H.ladkin.goodhit_hitphi_0", &kin_hitphi_0);
  T->SetBranchAddress("H.ladkin.goodhit_hitphi_1", &kin_hitphi_1);
  T->SetBranchAddress("H.ladkin.goodhit_hitedep_0", &kin_hitedep_0);
  T->SetBranchAddress("H.ladkin.goodhit_hitedep_1", &kin_hitedep_1);
  T->SetBranchAddress("H.ladkin.goodhit_deltapostrans_0", &kin_deltapostrans_0);
  T->SetBranchAddress("H.ladkin.goodhit_deltapostrans_1", &kin_deltapostrans_1);
  T->SetBranchAddress("H.ladkin.goodhit_deltaposlong_0", &kin_deltaposlong_0);
  T->SetBranchAddress("H.ladkin.goodhit_deltaposlong_1", &kin_deltaposlong_1);
  T->SetBranchAddress("H.react.x", &vertex_x);
  T->SetBranchAddress("H.react.y", &vertex_y);
  T->SetBranchAddress("H.react.z", &vertex_z);

  ///////////////////////////////////////////////////////////////////

  // Loop through the tree entries
  Long64_t nEntries = T->GetEntries();
  // nEntries          = 10000; // For testing purposes, limit to 1000 entries
  for (Long64_t i = 0; i < nEntries; ++i) {
    T->GetEntry(i);
    // Skip processing if there are too many tracks
    if (nTracks > maxTracks) {
      continue;
    }

    // Create a vector to label tracks as good or not
    std::vector<bool> goodTrack(nTracks, false);
    // Create vectors to store hodoscope paddle, plane, and y position
    std::vector<int> hodoPaddle(nTracks, -1);
    std::vector<int> hodoPlane(nTracks, -1);
    std::vector<double> hodoYPos(nTracks, 0.0);
    // Loop through all hodo hits
    for (int j = 0; j < nGoodHits; ++j) {
      // Get the track ID associated with the hodo hit
      int trackID = static_cast<int>(kin_trackID_0[j]);

      // Ensure the track ID is within bounds
      if (trackID < 0 || trackID >= nTracks) {
        continue;
      }

      // Apply the cuts
      if (fabs(trk_d0[trackID]) < d0_cut && // d0 cut
          kin_deltapostrans_0[j] > delta_trans_cut[0] &&
          kin_deltapostrans_0[j] < delta_trans_cut[1] && // delta_pos_trans cut
          kin_deltaposlong_0[j] > delta_long_cut[0] &&
          kin_deltaposlong_0[j] < delta_long_cut[1]) { // delta_pos_long cut
        goodTrack[trackID]  = true;
        hodoPaddle[trackID] = static_cast<int>(kin_paddle_0[j]);
        hodoPlane[trackID]  = static_cast<int>(kin_plane_0[j]);
        // TODO: Need to get Y position somehow
      }
    }

    // Loop through all tracks
    for (int j = 0; j < nTracks; ++j) {
      // Fill the histograms based on the goodTrack vector
      if (!goodTrack[j]) {
        continue;
      }
      TVector3 p_hit = GetHodoHitPosition(hodoPaddle[j], hodoPlane[j]);
      TVector3 p_react(vertex_x, vertex_y, vertex_z);
      bool good_vertez = true;
      if (nFixedz <= 0) {
        if (vertex_x > 1000 || vertex_y > 1000 || vertex_z > 1000) {
          good_vertez = false;
        }
        if (vertex_x > 1000)
          p_react.SetX(0);
        if (vertex_y > 1000)
          p_react.SetY(0);
        if (vertex_z > 1000)
          p_react.SetZ(0);
      } else {
        good_vertez = false;
        for (int k = 0; k < nFixedz; ++k) {
          if (use_projz) {
            if (trk_projz[j] > z_range[k] && trk_projz[j] < z_range[k + 1]) {
              good_vertez = true;
              p_react.SetZ(fixedz[k]);
              break;
            }
          } else {
            if (vertex_z > z_range[k] && vertex_z < z_range[k + 1]) {
              good_vertez = true;
              p_react.SetZ(fixedz[k]);
              break;
            }
          }
        }
      }

      // Fill the goodTracks vector
      good_track track;
      track.trk_x[0] = trk_x_local[0][j];
      track.trk_y[0] = trk_y_local[0][j];
      track.trk_x[1] = trk_x_local[1][j];
      track.trk_y[1] = trk_y_local[1][j];
      track.p_hit    = p_hit;
      track.p_react  = TVector3(vertex_x, vertex_y, vertex_z);
      goodTracks.push_back(track);
    }
    // Print the status as a percentage
    if (i % (nEntries / 100) == 0) {
      std::cout << "\rProcessing: " << int(i * 100.0 / nEntries) << "% completed." << std::flush;
    }
  }
  // End Event Loop
  std::cout << "\nProcessing completed." << std::endl;

  // Define the initial fit parameters
  fit_param bestFitParams;
  double minCost = std::numeric_limits<double>::max();

  auto MinimizeCostFunction = [](int &npar, double *grad, double &fval, double *par, int flag) {
    fit_param currentParams = {par[0], 0.0, par[1], 0.0, par[2], par[2], par[3], par[4]};
    fval                    = CostFunction(goodTracks, currentParams);
  };

  TMinuit minuit(5); // 5 parameters to fit
  minuit.SetFCN(MinimizeCostFunction);

  // Set parameter names, initial values, and ranges
  minuit.DefineParameter(0, "gem0_dx", 0.0, 0.1, gem0_x_range.min, gem0_x_range.max);
  minuit.DefineParameter(1, "gem1_dx", 0.0, 0.1, gem1_x_range.min, gem1_x_range.max);
  minuit.DefineParameter(2, "gem_theta", 127, 0.01, gem_theta_range.min, gem_theta_range.max);
  minuit.DefineParameter(3, "gem0_z", 77.571, 0.1, gem0_z_range.min, gem0_z_range.max);
  minuit.DefineParameter(4, "gem1_z", 95.571, 0.1, gem1_z_range.min, gem1_z_range.max);

  // Perform the minimization
  minuit.Migrad();

  // Retrieve the best fit parameters
  double par[5], err[5];
  for (int i = 0; i < 5; ++i) {
    minuit.GetParameter(i, par[i], err[i]);
  }

  // Print the best fit parameters and their errors
  std::cout << "Best Fit Parameters and Errors:" << std::endl;
  std::cout << "GEM0 dx: " << par[0] << " ± " << err[0] << std::endl;
  std::cout << "GEM1 dx: " << par[1] << " ± " << err[1] << std::endl;
  std::cout << "GEM theta: " << par[2] << " ± " << err[2] << std::endl;
  std::cout << "GEM0 z: " << par[3] << " ± " << err[3] << std::endl;
  std::cout << "GEM1 z: " << par[4] << " ± " << err[4] << std::endl;
  bestFitParams.gem0_dx    = par[0];
  bestFitParams.gem1_dx    = par[1];
  bestFitParams.gem0_theta = par[2];
  bestFitParams.gem1_theta = par[2];
  bestFitParams.gem0_z     = par[3];
  bestFitParams.gem1_z     = par[4];
  minCost                  = CostFunction(goodTracks, bestFitParams);

  int nSteps =
      gem0_x_range.nsteps * gem0_z_range.nsteps * gem1_x_range.nsteps * gem1_z_range.nsteps * gem_theta_range.nsteps;
  int currentStep = 0;

  // // Iterate over the parameter ranges
  // for (double gem0_dx = gem0_x_range.min; gem0_dx <= gem0_x_range.max;
  //      gem0_dx += (gem0_x_range.max - gem0_x_range.min) / gem0_x_range.nsteps) {
  //   for (double gem1_dx = gem1_x_range.min; gem1_dx <= gem1_x_range.max;
  //        gem1_dx += (gem1_x_range.max - gem1_x_range.min) / gem1_x_range.nsteps) {
  //     for (double gem0_z = gem0_z_range.min; gem0_z <= gem0_z_range.max;
  //          gem0_z += (gem0_z_range.max - gem0_z_range.min) / gem0_z_range.nsteps) {
  //       for (double gem1_z = gem1_z_range.min; gem1_z <= gem1_z_range.max;
  //            gem1_z += (gem1_z_range.max - gem1_z_range.min) / gem1_z_range.nsteps) {
  //         for (double gem_theta = gem_theta_range.min; gem_theta <= gem_theta_range.max;
  //              gem_theta += (gem_theta_range.max - gem_theta_range.min) / gem_theta_range.nsteps) {
  //           // Set the current fit parameters
  //           fit_param currentParams = {gem0_dx, 0.0, gem1_dx, 0.0, gem_theta, gem_theta, gem0_z, gem1_z};

  //           // Calculate the cost function
  //           double cost = CostFunction(goodTracks, currentParams);

  //           currentStep++;
  //           // Print the progress
  //           if (currentStep % (nSteps / 100) == 0) {
  //             std::cout << "\rProcessing: " << int(currentStep * 100.0 / nSteps) << "% completed." << std::flush;
  //           }
  //           // Update the best fit parameters if a lower cost is found
  //           if (cost < minCost) {
  //             minCost       = cost;
  //             bestFitParams = currentParams;
  //           }
  //         }
  //       }
  //     }
  //   }
  // }

  // // // Output the best fit parameters
  // std::cout << "Best Fit Parameters:" << std::endl;
  // std::cout << "GEM0 dx: " << bestFitParams.gem0_dx << ", GEM0 z: " << bestFitParams.gem0_z << std::endl;
  // std::cout << "GEM1 dx: " << bestFitParams.gem1_dx << ", GEM1 z: " << bestFitParams.gem1_z << std::endl;
  // std::cout << "GEM theta: " << bestFitParams.gem0_theta << std::endl;
  // std::cout << "Minimum Cost: " << minCost << std::endl;

  // Close the file
  file->Close();

  return;
}