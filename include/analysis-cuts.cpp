#include "analysis-cuts.h"
#include "TCut.h"
#include "TMath.h"
#include "TString.h"
#include "analysis-constants.h"
#include "analysis-binning.h"

// Trigger cuts
const bool trigger_lines = true;

// Jet cuts
const double jet_eta_min = 2.5;
const double jet_eta_max = 4.0;
const double jet_pt_min  = 15;

// Jet ID cuts
const double mpf_max = 0.8; 
const double cpf_min = 0.1;
const double mpt_min = 1.2;
const double nPVtrks_min = 1.5;

// Topological cuts
const double deltaphi_z_jet_min  = 7*TMath::Pi()/8.;
const double jet_radius = 0.5;

// Z boson cuts
const double muon_pt_min  = 20.;
const double lhcb_eta_min = 2;
const double lhcb_eta_max = 4.5;
const double dimuon_mass_min = 60;
const double dimuon_mass_max = 120;
const double muon_trackprob_min = 0.001;

// Track cuts
const double track_chi2ndf_max     = 3;
const double track_p_min           = 10;
const double track_p_max           = 1000;
const double track_pt_min          = 0.25;
const double track_probnnghost_max = 0.5;

// Jet pt dependent cuts
const double weight_pt_cut[]  = {0.12,0.10,0.085};
const double min_efficiency[] = {0.14,0.10,0.05};
const double max_relerror[]   = {0.50,0.35,0.45};

TCut pair_matching_cut = "R_L_truth!=-999";

TCut eec_cut  = Form("weight*(jet_pt>%f&&jet_pt<%f)",jet_pt_min_nom,jet_pt_max);
TCut pair_cut = Form("jet_pt>%f&&jet_pt<%f",jet_pt_min_nom,jet_pt_max);

TCut pair_jet_pt_cut[] = {
	Form("jet_pt>%f&&jet_pt<%f", jet_pt_binning[0], jet_pt_binning[1]),
	Form("jet_pt>%f&&jet_pt<%f", jet_pt_binning[1], jet_pt_binning[2]),
	Form("jet_pt>%f&&jet_pt<%f", jet_pt_binning[2], jet_pt_binning[3])
};

TCut pair_signal_cut   = Form("TMath::Abs(R_L_truth-R_L)<%f",rl_resolution);
TCut single_signal_cut = "key_match==1";

bool apply_jet_cuts(double jet_eta, double jet_pt)
{
        if (jet_eta < jet_eta_min || jet_eta > jet_eta_max) 
                return false;
        
        if (jet_pt < unfolding_jet_pt_binning[0])
                return false;

        return true;    
}

bool apply_jet_id_cuts(double mpt, double nPVtrk, double cpf, double mpf)
{
        if (mpt < mpt_min)
                return false;
        
        if (nPVtrk < nPVtrks_min)
                return false;
        
        if (cpf < cpf_min)
                return false;
        
        if (mpf > mpf_max)
                return false;

        return true;    
}

bool apply_muon_cuts(double deltaR_mu_jet, double mu_pt, double mu_eta)
{
        if (deltaR_mu_jet < jet_radius) 
                return false; 
        
        if (mu_pt < muon_pt_min)
                return false;

        if (mu_eta < lhcb_eta_min || mu_eta > lhcb_eta_max) 
                return false;

        return true;
}

bool apply_zboson_cuts(double deltaphi_zboson_jet, double zboson_mass)
{
    if (deltaphi_zboson_jet < deltaphi_z_jet_min) 
            return false;

    if (zboson_mass < dimuon_mass_min || zboson_mass > dimuon_mass_max) 
            return false;

    return true;
}

bool apply_chargedtrack_cuts(double charge, double p, double pt, double chi2ndf, double probnnghost, double eta)
{
        if (charge == 0)
                return false;

        if (p < track_p_min || p > track_p_max)
                return false;

        if (pt < track_pt_min)
                return false;

        if (chi2ndf > track_chi2ndf_max)
                return false;

        if (probnnghost > track_probnnghost_max)
                return false;

        if (eta < lhcb_eta_min || eta > lhcb_eta_max)
                return false;

        return true;
}

bool apply_chargedtrack_cuts(double charge, double p, double pt, double chi2ndf, double probnnghost, double eta, double deltaR_h_jet)
{
        if (charge == 0)
                return false;

        if (p < track_p_min || p > track_p_max)
                return false;

        if (pt < track_pt_min)
                return false;

        if (chi2ndf > track_chi2ndf_max)
                return false;

        if (probnnghost > track_probnnghost_max)
                return false;

        if (eta < lhcb_eta_min || eta > lhcb_eta_max)
                return false;

        if (deltaR_h_jet > jet_radius)
                return false;

        return true;
}

bool apply_chargedtrack_momentum_cuts(double charge, double p, double pt, double eta)
{
        if (charge == 0)
                return false;

        if (p < track_p_min || p > track_p_max)
                return false;

        if (pt < track_pt_min)
                return false;

        if (eta < lhcb_eta_min || eta > lhcb_eta_max)
                return false;

        return true;
}

bool apply_chargedtrack_momentum_cuts(double charge, double p, double pt, double eta, double deltaR_h_jet)
{
        if (charge == 0)
                return false;

        if (p < track_p_min || p > track_p_max)
                return false;

        if (pt < track_pt_min)
                return false;

        if (eta < lhcb_eta_min || eta > lhcb_eta_max)
                return false;

        if (deltaR_h_jet > jet_radius)
                return false;

        return true;
}
