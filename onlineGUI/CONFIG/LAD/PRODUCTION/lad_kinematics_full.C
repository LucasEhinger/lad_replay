void lad_kinematics_full(const char* spec = "P", const char* treename = "T") {
    TTree* tree = (TTree*)gDirectory->Get(treename);
    if (!tree) {
        std::cerr << "Tree '" << treename << "' not found!" << std::endl;
        return;
    }

    TString prefix = spec;
    prefix.ToUpper();

    // Dynamic branch names
    TString kin_qx     = Form("%s.kin.q_x", prefix.Data());
    TString kin_qy     = Form("%s.kin.q_y", prefix.Data());
    TString kin_qz     = Form("%s.kin.q_z", prefix.Data());
    TString kin_Q2     = Form("%s.kin.Q2", prefix.Data());
    TString kin_omega  = Form("%s.kin.omega", prefix.Data());
    TString react_z = Form("%s.react.z", prefix.Data());
    TString lad_beta   = Form("%s.ladhod.goodhit_beta", prefix.Data());
    TString lad_theta  = Form("%s.ladhod.goodhit_hit_theta", prefix.Data());
    TString lad_phi    = Form("%s.ladhod.goodhit_hit_phi", prefix.Data());
    TString lad_delta_long = Form("%s.ladhod.goodhit_delta_pos_long", prefix.Data());

    // Check if all required branches exist
    std::vector<TString> requiredBranches = {
        kin_qx, kin_qy, kin_qz, kin_Q2, kin_omega, react_z,
        lad_beta, lad_theta, lad_phi, lad_delta_long
    };
    for (auto& bname : requiredBranches) {
        if (!tree->GetBranch(bname)) {
            std::cerr << "Branch '" << bname << "' not found. Skipping macro." << std::endl;
            return;
        }
    }

    // Variables
    double qx, qy, qz, Q2, omega, v_e;
    double beta, theta_deg, phi_deg, v_p;

    tree->SetBranchAddress(kin_qx, &qx);
    tree->SetBranchAddress(kin_qy, &qy);
    tree->SetBranchAddress(kin_qz, &qz);
    tree->SetBranchAddress(kin_Q2, &Q2);
    tree->SetBranchAddress(kin_omega, &omega);
    tree->SetBranchAddress(react_z, &v_e);
    tree->SetBranchAddress(lad_beta, &beta);
    tree->SetBranchAddress(lad_theta, &theta_deg);
    tree->SetBranchAddress(lad_phi, &phi_deg);
    tree->SetBranchAddress(lad_delta_long, &v_p);

    const double Mp = 0.938;
    const double Md = 1.8756;

    // Histograms
    TH1D* h_xprime            = new TH1D("h_xprime", "x';x';Counts", 100, 0, 2);
    TH1D* h_alpha             = new TH1D("h_alpha", "#alpha;#alpha;Counts", 100, 0, 3);
    TH1D* h_pperp             = new TH1D("h_pperp", "p_{#perp};p_{#perp} [GeV];Counts", 100, 0, 1);
    TH1D* h_theta_p           = new TH1D("h_theta_p", "#theta_{p};#theta_{p} [deg];Counts", 100, 0, 180);
    TH1D* h_pp                = new TH1D("h_pp", "p_{p};p_{p} [GeV];Counts", 100, 0, 1.5);
    TH1D* h_theta_pq          = new TH1D("h_theta_pq", "#theta_{pq};#theta_{pq} [deg];Counts", 100, 0, 180);
    TH2D* h2_theta_pq_vs_pp   = new TH2D("h2_theta_pq_vs_pp", "#theta_{pq} vs p_{p};p_{p} [GeV];#theta_{pq} [deg]", 100, 0, 1.5, 100, 0, 180);
    TH1D* h_dv                = new TH1D("h_dv", "#Delta v = v_{e} - v_{p};#Delta v;Counts", 100, -2, 2);

    // Optional: not drawn
    //TH1D* h_ve                = new TH1D("h_ve", "v_{e};v_{e};Counts", 100, 0, 10);
    //TH1D* h_vp                = new TH1D("h_vp", "v_{p};v_{p};Counts", 100, 0, 10);
    TH2D* h2_xprime_vs_alpha  = new TH2D("h2_xprime_vs_alpha", "x' vs #alpha;#alpha;x'", 100, 0, 3, 100, 0, 2);
    TH2D* h2_theta_p_vs_pp    = new TH2D("h2_theta_p_vs_pp", "#theta_{p} vs p_{p};p_{p} [GeV];#theta_{p} [deg]", 100, 0, 1.5, 100, 0, 180);
    TH2D* h2_ve_vs_vp         = new TH2D("h2_ve_vs_vp", "v_{e} vs v_{p};v_{p};v_{e}", 100, 0, 3, 100, 0, 3);

    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);
        double theta = TMath::DegToRad() * theta_deg;
        double phi   = TMath::DegToRad() * phi_deg;

        double gamma = 1.0 / sqrt(1 - beta * beta);
        double Er = gamma * Mp;
        double pr = sqrt(Er * Er - Mp * Mp);

        double px = pr * sin(theta) * cos(phi);
        double py = pr * sin(theta) * sin(phi);
        double pz = pr * cos(theta);

        TVector3 p_recoil(px, py, pz);
        TVector3 q(qx, qy, qz);
        double qmag = q.Mag();

        double alpha = (Er - p_recoil.Dot(q.Unit())) / Mp;
        double xprime = Q2 / (2.0 * ((Md - Er) * omega - p_recoil.Dot(q)));
        double pperp = (p_recoil.Cross(q.Unit())).Mag();
        double theta_pq = p_recoil.Angle(q) * TMath::RadToDeg();
        double ve = v_e;
        double vp = v_p;
        double dv = ve - vp;

        //double W2 = Mp * Mp + 2 * Mp * omega - Q2;
        //double W = (W2 > 0) ? sqrt(W2) : 0;

        // if (pperp > 0.3) continue;
        // if (fabs(dv) > 0.2) continue;
        // if (W < 1.8) continue;

        h_xprime->Fill(xprime);
        h_alpha->Fill(alpha);
        h_pperp->Fill(pperp);
        h_theta_p->Fill(theta_deg);
        h_pp->Fill(pr);
        h_theta_pq->Fill(theta_pq);
        h2_theta_pq_vs_pp->Fill(pr, theta_pq);
        h_dv->Fill(dv);
        h2_xprime_vs_alpha->Fill(alpha, xprime);
        h2_theta_p_vs_pp->Fill(pr, theta_deg);
        h2_ve_vs_vp->Fill(vp, ve);
        // h_ve->Fill(ve);
        // h_vp->Fill(vp);
    }

    TCanvas* c = (TCanvas*)gPad->GetCanvas();
    if (!c) {
        std::cerr << "Error: No GUI canvas detected. Exiting draw." << std::endl;
        return;
    }

    c->cd(1);  h_xprime->Draw();
    c->cd(2);  h_alpha->Draw();
    c->cd(3);  h_pperp->Draw();
    c->cd(4);  h_theta_p->Draw();
    c->cd(5);  h_pp->Draw();
    c->cd(6);  h_theta_pq->Draw();
    c->cd(7);  h2_theta_pq_vs_pp->Draw("COLZ");
    //c->cd(8);  h_ve->Draw();         // optional
    // c->cd(9);  h_vp->Draw();         // optional
    c->cd(8);  h_dv->Draw();         // optional
    c->cd(9);  h2_xprime_vs_alpha->Draw("COLZ");
    c->cd(10);  h2_theta_p_vs_pp->Draw("COLZ");
    c->cd(11); h2_ve_vs_vp->Draw("COLZ");  // optional
}
