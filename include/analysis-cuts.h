#ifndef ANALYSIS_CUTS_H
#define ANALYSIS_CUTS_H

#include "TCut.h"
#include "TString.h"
#include "analysis-constants.h"

// _______________________________ Z Tagged Jet Cuts _______________________________ //

// Trigger cuts
const bool trigger_lines = true;

// Jet cuts
const double jet_eta_min = 2.5;
const double jet_eta_max = 4.0;
const double jet_pt_min  = 15;

// Topological cuts
const double deltaphi_z_jet_min  = 7*TMath::Pi()/8.;
const double jet_radius = 0.5;

// Z boson cuts
const double muon_pt_min  = 20.;
const double muon_eta_min = 2;
const double muon_eta_max = 4.5;
const double muonmuon_mass_min = 60;
const double muonmuon_mass_max = 120;
const double muon_trackprob_min = 0.001;

// Track cuts
const double track_chi2ndf_max     = 3;
const double track_p_min           = 4;
const double track_p_max           = 1000;
const double track_pt_min          = 0.25;
const double track_probnnghost_max = 0.5;

// _______________________________ Nominal Analysis Cuts _______________________________ //

TCut e2c_cut  = Form("weight*(jet_pt>%f&&jet_pt<%f)",jet_pt_min_nom,jet_pt_max);
TCut pair_cut = Form("jet_pt>%f&&jet_pt<%f",jet_pt_min_nom,jet_pt_max);

// Cuts as function of jet pt
TCut e2c_jetpt_cut[] = 
{Form("weight*(jet_pt>%f&&jet_pt<%f)",jet_pt_binning[0],jet_pt_binning[1]),
 Form("weight*(jet_pt>%f&&jet_pt<%f)",jet_pt_binning[1],jet_pt_binning[2]),
 Form("weight*(jet_pt>%f&&jet_pt<%f)",jet_pt_binning[2],jet_pt_binning[3])};

TCut pair_jetpt_cut[] = 
{Form("jet_pt>%f&&jet_pt<%f",jet_pt_binning[0],jet_pt_binning[1]),
 Form("jet_pt>%f&&jet_pt<%f",jet_pt_binning[1],jet_pt_binning[2]),
 Form("jet_pt>%f&&jet_pt<%f",jet_pt_binning[2],jet_pt_binning[3])};

// _______________________________ Purity Analysis Cuts _______________________________ //
TCut pair_signal_cut   = Form("TMath::Abs(R_L_truth-R_L)<%f&&jet_pt>%f&&jet_pt<%f",R_L_res,jet_pt_min_nom,jet_pt_max);
TCut pair_pairbg_cut   = Form("(h1truth_y==-999&&h2truth_y==-999)&&jet_pt>%f&&jet_pt<%f",jet_pt_min_nom,jet_pt_max);
TCut pair_singlebg_cut = Form("((h1truth_y==-999&&h2truth_y!=-999)||(h1truth_y!=-999&&h2truth_y==-999))&&jet_pt>%f&&jet_pt<%f",jet_pt_min_nom,jet_pt_max);
TCut single_signal_cut = Form("key_match==1&&jet_pt>%f&&jet_pt<%f",jet_pt_min_nom,jet_pt_max); // This was designed for single particle tuples

TCut e2c_signal_cut   = Form("weight*(TMath::Abs(R_L_truth-R_L)<%f&&jet_pt>%f&&jet_pt<%f)",R_L_res,jet_pt_min_nom,jet_pt_max);
TCut e2c_pairbg_cut   = Form("weight*((h1truth_y==-999&&h2truth_y==-999)&&jet_pt>%f&&jet_pt<%f)",jet_pt_min_nom,jet_pt_max);
TCut e2c_singlebg_cut = Form("weight*(((h1truth_y==-999&&h2truth_y!=-999)||(h1truth_y!=-999&&h2truth_y==-999))&&jet_pt>%f&&jet_pt<%f)",jet_pt_min_nom,jet_pt_max);

TCut purity_corr_singletrack         = Form("purity*(purity_relerror<%f&&jet_pt>%f&&jet_pt<%f)",corr_rel_error,jet_pt_min_nom,jet_pt_max);
TCut efficiency_corr_singletrack     = Form("(1./efficiency)*(efficiency_relerror<%f&&jet_pt>%f&&jet_pt<%f)",corr_rel_error,jet_pt_min_nom,jet_pt_max);
TCut full_corr_singletrack           = Form("purity*(1./efficiency)*(efficiency_relerror<%f&&purity_relerror<%f&&jet_pt>%f&&jet_pt<%f)",corr_rel_error,corr_rel_error,jet_pt_min_nom,jet_pt_max);

TCut e2c_purity_corr_singletrack     = Form("weight*purity*(purity_relerror<%f&&jet_pt>%f&&jet_pt<%f)",corr_rel_error,jet_pt_min_nom,jet_pt_max);
TCut e2c_efficiency_corr_singletrack = Form("weight*(1./efficiency)*(efficiency_relerror<%f&&jet_pt>%f&&jet_pt<%f)",corr_rel_error,jet_pt_min_nom,jet_pt_max);
TCut e2c_full_corr_singletrack       = Form("weight*purity*(1./efficiency)*(efficiency_relerror<%f&&purity_relerror<%f&&jet_pt>%f&&jet_pt<%f)",corr_rel_error,corr_rel_error,jet_pt_min_nom,jet_pt_max);

TCut pair_jetpt_signal_cut[] = 
{
Form("TMath::Abs(R_L_truth-R_L)<%f&&jet_pt>%f&&jet_pt<%f",R_L_res,jet_pt_binning[0],jet_pt_binning[1]),
Form("TMath::Abs(R_L_truth-R_L)<%f&&jet_pt>%f&&jet_pt<%f",R_L_res,jet_pt_binning[1],jet_pt_binning[2]),
Form("TMath::Abs(R_L_truth-R_L)<%f&&jet_pt>%f&&jet_pt<%f",R_L_res,jet_pt_binning[2],jet_pt_binning[3])
};

TCut pair_jetpt_pairbg_cut[] = 
{
Form("(h1truth_y==-999&&h2truth_y==-999)&&jet_pt>%f&&jet_pt<%f",jet_pt_binning[0],jet_pt_binning[1]),
Form("(h1truth_y==-999&&h2truth_y==-999)&&jet_pt>%f&&jet_pt<%f",jet_pt_binning[1],jet_pt_binning[2]),
Form("(h1truth_y==-999&&h2truth_y==-999)&&jet_pt>%f&&jet_pt<%f",jet_pt_binning[2],jet_pt_binning[3])
};

TCut pair_jetpt_singlebg_cut[] = 
{
Form("((h1truth_y==-999&&h2truth_y!=-999)||(h1truth_y!=-999&&h2truth_y==-999))&&jet_pt>%f&&jet_pt<%f",jet_pt_binning[0],jet_pt_binning[1]),
Form("((h1truth_y==-999&&h2truth_y!=-999)||(h1truth_y!=-999&&h2truth_y==-999))&&jet_pt>%f&&jet_pt<%f",jet_pt_binning[1],jet_pt_binning[2]),
Form("((h1truth_y==-999&&h2truth_y!=-999)||(h1truth_y!=-999&&h2truth_y==-999))&&jet_pt>%f&&jet_pt<%f",jet_pt_binning[2],jet_pt_binning[3])
};

// FUNCTIONS TO APPLY CUTS
bool apply_jet_cuts(double jet_eta, double jet_pt)
{
    if(jet_eta<jet_eta_min||jet_eta>jet_eta_max) return false;
    if(jet_pt<em_jetptunfolding_binning[0]) return false;

    return true;    
}

bool apply_muon_cuts(double deltaR_mu_jet, double mu_pt, double mu_eta)
{
    if(deltaR_mu_jet<jet_radius) return false; 
    if(mu_pt<muon_pt_min) return false;
    if(mu_eta<muon_eta_min||mu_eta>muon_eta_max) return false;

    return true;
}

bool apply_zboson_cuts(double deltaphi_zboson_jet, double zboson_mass)
{
    if(deltaphi_zboson_jet<deltaphi_z_jet_min) return false;
    if(zboson_mass<muonmuon_mass_min||zboson_mass>muonmuon_mass_max) return false;

    return true;
}

bool apply_chargedtrack_cuts(double charge, double p, double pt, double chi2ndf, double probnnghost, double deltaR_h_jet)
{
    if(charge==0) return false;
    if(p<track_p_min||p>track_p_max) return false;
    if(pt<track_pt_min) return false;
    if(chi2ndf>track_chi2ndf_max) return false;
    if(probnnghost>track_probnnghost_max) return false;
    if(deltaR_h_jet>jet_radius) return false;

    return true;
}

bool apply_chargedtrack_momentum_cuts(double charge, double p, double pt, double deltaR_h_jet)
{
    if(charge==0) return false;
    if(p<track_p_min||p>track_p_max) return false;
    if(pt<track_pt_min) return false;
    if(deltaR_h_jet>jet_radius) return false;

    return true;
}

#endif