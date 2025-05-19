#ifndef ANALYSIS_CUTS_H
#define ANALYSIS_CUTS_H

#include "TCut.h"
#include "TString.h"
#include "analysis-constants.h"
#include "analysis-binning.h"

// _______________________________ Z Tagged Jet Cuts _______________________________ //

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

TCut e2c_jetpt_cut_weightpt[] = 
{Form("weight_pt*(jet_pt>%f&&jet_pt<%f)",jet_pt_binning[0],jet_pt_binning[1]),
 Form("weight_pt*(jet_pt>%f&&jet_pt<%f)",jet_pt_binning[1],jet_pt_binning[2]),
 Form("weight_pt*(jet_pt>%f&&jet_pt<%f)",jet_pt_binning[2],jet_pt_binning[3])};

TCut e2c_zpt_cut_weightpt[] = 
{Form("weight_pt*(z_pt>%f&&z_pt<%f&&TMath::Abs(deltaphi_z_h1)>TMath::Pi()/2.&&TMath::Abs(deltaphi_z_h2)>TMath::Pi()/2.)",z_pt_binning[0],z_pt_binning[1]),
 Form("weight_pt*(z_pt>%f&&z_pt<%f&&TMath::Abs(deltaphi_z_h1)>TMath::Pi()/2.&&TMath::Abs(deltaphi_z_h2)>TMath::Pi()/2.)",z_pt_binning[1],z_pt_binning[2]),
 Form("weight_pt*(z_pt>%f&&z_pt<%f&&TMath::Abs(deltaphi_z_h1)>TMath::Pi()/2.&&TMath::Abs(deltaphi_z_h2)>TMath::Pi()/2.)",z_pt_binning[2],z_pt_binning[3])};

TCut pair_jetpt_cut[] = 
{Form("jet_pt>%f&&jet_pt<%f",jet_pt_binning[0],jet_pt_binning[1]),
 Form("jet_pt>%f&&jet_pt<%f",jet_pt_binning[1],jet_pt_binning[2]),
 Form("jet_pt>%f&&jet_pt<%f",jet_pt_binning[2],jet_pt_binning[3])};

TCut pair_zpt_cut[] = 
{Form("z_pt>%f&&z_pt<%f",z_pt_binning[0],z_pt_binning[1]),
 Form("z_pt>%f&&z_pt<%f",z_pt_binning[1],z_pt_binning[2]),
 Form("z_pt>%f&&z_pt<%f",z_pt_binning[2],z_pt_binning[3])};

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

TCut pair_purity_corr_singletrack_weightpt = Form("purity*(purity_relerror<%f)",corr_rel_error);

TString muons_eff = "1./(mum_eff_id*mup_eff_id*mum_eff_trk*mup_eff_trk*(mum_eff_trg+mup_eff_trg-mum_eff_trg*mup_eff_trg))";

// TCut jet_full_corr[] =
// {
// Form("jet_purity*(1./jet_efficiency)*(jet_pt>%f&&jet_pt<%f)",jet_pt_binning[0],jet_pt_binning[1]),
// Form("jet_purity*(1./jet_efficiency)*(jet_pt>%f&&jet_pt<%f)",jet_pt_binning[1],jet_pt_binning[2]),
// Form("jet_purity*(1./jet_efficiency)*(jet_pt>%f&&jet_pt<%f)",jet_pt_binning[2],jet_pt_binning[3])
// };

TCut jet_full_corr[] =
{
Form("jet_purity*(1./jet_efficiency)*"+muons_eff+"*(jet_pt>%f&&jet_pt<%f)",jet_pt_binning[0],jet_pt_binning[1]),
Form("jet_purity*(1./jet_efficiency)*"+muons_eff+"*(jet_pt>%f&&jet_pt<%f)",jet_pt_binning[1],jet_pt_binning[2]),
Form("jet_purity*(1./jet_efficiency)*"+muons_eff+"*(jet_pt>%f&&jet_pt<%f)",jet_pt_binning[2],jet_pt_binning[3])
};

TCut e2c_jetpt_full_corr_singletrack[] =
{
Form("weight*purity*(1./efficiency)*(efficiency_relerror<%f&&purity_relerror<%f&&jet_pt>%f&&jet_pt<%f)",corr_rel_error,corr_rel_error,jet_pt_binning[0],jet_pt_binning[1]),
Form("weight*purity*(1./efficiency)*(efficiency_relerror<%f&&purity_relerror<%f&&jet_pt>%f&&jet_pt<%f)",corr_rel_error,corr_rel_error,jet_pt_binning[1],jet_pt_binning[2]),
Form("weight*purity*(1./efficiency)*(efficiency_relerror<%f&&purity_relerror<%f&&jet_pt>%f&&jet_pt<%f)",corr_rel_error,corr_rel_error,jet_pt_binning[2],jet_pt_binning[3])
};

TCut e2c_jetpt_full_corr_singletrack_weightpt[] =
{
Form("weight_pt*purity*(1./efficiency)*(efficiency_relerror<%f&&purity_relerror<%f&&jet_pt>%f&&jet_pt<%f)",corr_rel_error,corr_rel_error,jet_pt_binning[0],jet_pt_binning[1]),
Form("weight_pt*purity*(1./efficiency)*(efficiency_relerror<%f&&purity_relerror<%f&&jet_pt>%f&&jet_pt<%f)",corr_rel_error,corr_rel_error,jet_pt_binning[1],jet_pt_binning[2]),
Form("weight_pt*purity*(1./efficiency)*(efficiency_relerror<%f&&purity_relerror<%f&&jet_pt>%f&&jet_pt<%f)",corr_rel_error,corr_rel_error,jet_pt_binning[2],jet_pt_binning[3])
};

TCut e2c_jetpt_purity_corr_singletrack_weightpt[] =
{
Form("weight_pt*purity*(purity_relerror<%f&&jet_pt>%f&&jet_pt<%f)",corr_rel_error,jet_pt_binning[0],jet_pt_binning[1]),
Form("weight_pt*purity*(purity_relerror<%f&&jet_pt>%f&&jet_pt<%f)",corr_rel_error,jet_pt_binning[1],jet_pt_binning[2]),
Form("weight_pt*purity*(purity_relerror<%f&&jet_pt>%f&&jet_pt<%f)",corr_rel_error,jet_pt_binning[2],jet_pt_binning[3])
};

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

TCut jet_signal_cut = Form("jet_pt>%f&&jet_pt<%f&&jet_pt_truth!=999",jet_pt_min_nom,jet_pt_max);

// ANALYSIS VARIATIONS CUTS
TCut e2c_jetpt_cut_weightpt_pp[] = 
{Form("weight_pt*(h1_charge==3&&h2_charge==3&&jet_pt>%f&&jet_pt<%f)",jet_pt_binning[0],jet_pt_binning[1]),
 Form("weight_pt*(h1_charge==3&&h2_charge==3&&jet_pt>%f&&jet_pt<%f)",jet_pt_binning[1],jet_pt_binning[2]),
 Form("weight_pt*(h1_charge==3&&h2_charge==3&&jet_pt>%f&&jet_pt<%f)",jet_pt_binning[2],jet_pt_binning[3])};

TCut e2c_jetpt_cut_weightpt_mm[] = 
{Form("weight_pt*(h1_charge==-3&&h2_charge==-3&&jet_pt>%f&&jet_pt<%f)",jet_pt_binning[0],jet_pt_binning[1]),
 Form("weight_pt*(h1_charge==-3&&h2_charge==-3&&jet_pt>%f&&jet_pt<%f)",jet_pt_binning[1],jet_pt_binning[2]),
 Form("weight_pt*(h1_charge==-3&&h2_charge==-3&&jet_pt>%f&&jet_pt<%f)",jet_pt_binning[2],jet_pt_binning[3])};

TCut e2c_jetpt_cut_weightpt_pm[] = 
{Form("weight_pt*(h1_charge*h2_charge<0&&jet_pt>%f&&jet_pt<%f)",jet_pt_binning[0],jet_pt_binning[1]),
 Form("weight_pt*(h1_charge*h2_charge<0&&jet_pt>%f&&jet_pt<%f)",jet_pt_binning[1],jet_pt_binning[2]),
 Form("weight_pt*(h1_charge*h2_charge<0&&jet_pt>%f&&jet_pt<%f)",jet_pt_binning[2],jet_pt_binning[3])};

TCut e2c_jetpt_cut_weightpt_kaon[] = 
{Form("weight_pt*((TMath::Abs(h1_pid)==321||TMath::Abs(h2_pid)==321)&&jet_pt>%f&&jet_pt<%f)",jet_pt_binning[0],jet_pt_binning[1]),
 Form("weight_pt*((TMath::Abs(h1_pid)==321||TMath::Abs(h2_pid)==321)&&jet_pt>%f&&jet_pt<%f)",jet_pt_binning[1],jet_pt_binning[2]),
 Form("weight_pt*((TMath::Abs(h1_pid)==321||TMath::Abs(h2_pid)==321)&&jet_pt>%f&&jet_pt<%f)",jet_pt_binning[2],jet_pt_binning[3])};

TCut e2c_jetpt_cut_weightpt_nokaon[] = 
{Form("weight_pt*((TMath::Abs(h1_pid)!=321&&TMath::Abs(h2_pid)!=321)&&jet_pt>%f&&jet_pt<%f)",jet_pt_binning[0],jet_pt_binning[1]),
 Form("weight_pt*((TMath::Abs(h1_pid)!=321&&TMath::Abs(h2_pid)!=321)&&jet_pt>%f&&jet_pt<%f)",jet_pt_binning[1],jet_pt_binning[2]),
 Form("weight_pt*((TMath::Abs(h1_pid)!=321&&TMath::Abs(h2_pid)!=321)&&jet_pt>%f&&jet_pt<%f)",jet_pt_binning[2],jet_pt_binning[3])};

// FUNCTIONS TO APPLY CUTS
bool apply_jet_cuts(double jet_eta, double jet_pt)
{
    if(jet_eta<jet_eta_min||jet_eta>jet_eta_max) return false;
    if(jet_pt<unfolding_jetpt_binning[0]) return false;

    return true;    
}

bool apply_jet_id_cuts(double mpt, double nPVtrk, double cpf, double mpf)
{
    if(mpt < mpt_min) return false;
    if(nPVtrk < nPVtrks_min) return false;
    if(cpf < cpf_min) return false;
    if(mpf > mpf_max) return false;

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

bool apply_chargedtrack_cuts(double charge, double p, double pt, double chi2ndf, double probnnghost, double eta)
{
    if(charge==0) return false;
    if(p<track_p_min||p>track_p_max) return false;
    if(pt<track_pt_min) return false;
    if(chi2ndf>track_chi2ndf_max) return false;
    if(probnnghost>track_probnnghost_max) return false;
    if(eta<muon_eta_min||eta>muon_eta_max) return false;

    return true;
}

bool apply_chargedtrack_momentum_cuts(double charge, double p, double pt, double eta)
{
    if(charge==0) return false;
    if(p<track_p_min||p>track_p_max) return false;
    if(pt<track_pt_min) return false;
    if(eta<muon_eta_min||eta>muon_eta_max) return false;
    
    return true;
}

#endif