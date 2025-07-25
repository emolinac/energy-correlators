#ifndef NAMES_H
#define NAMES_H

#include<map>

// Name of the ntuples
std::string name_ntuple_unfold          = "ntuple_unfold";
std::string name_ntuple_data            = "ntuple_hadron";
std::string name_ntuple_mcjetmatch      = "ntuple_mcjetmatch";
std::string name_ntuple_mcreco          = "ntuple_mcreco";
std::string name_ntuple_mc              = "ntuple_mc";
std::string name_ntuple_purity          = "ntuple_purity";
std::string name_ntuple_unfolding       = "ntuple_unfolding";
std::string name_ntuple_correction_mc   = "ntuple_correction_mc";
std::string name_ntuple_correction_reco = "ntuple_correction_reco";
std::string name_ntuple_mc_jet          = "ntuple_mc_jet";
std::string name_ntuple_mcreco_jet      = "ntuple_mcreco_jet";


// Names of the files
std::string namef_ntuple_e2c                = "ntuple_e2c.root";
std::string namef_ntuple_e2c_corr           = "ntuple_corre2c.root";

std::string namef_ntuple_e2c_paircorr       = "ntuple_paircorre2c.root";
std::string namef_ntuple_e2c_paircorr_ct    = "ntuple_paircorre2c_ct.root";
std::string namef_ntuple_e2c_paircorr_jes   = "ntuple_paircorre2c_jes.root";
std::string namef_ntuple_e2c_paircorr_jer   = "ntuple_paircorre2c_jer.root";
std::string namef_ntuple_e2c_paircorr_prior = "ntuple_paircorre2c_prior.root";

std::string namef_ntuple_e2c_hadroncorrections   = "ntuple_e2c_hadroncorrections.root";

std::string namef_ntuple_e2c_paircorrections       = "ntuple_e2c_paircorrections.root";
std::string namef_ntuple_e2c_paircorrections_ct    = "ntuple_e2c_paircorrections_ct.root";
std::string namef_ntuple_e2c_paircorrections_jes   = "ntuple_e2c_paircorrections_jes.root";
std::string namef_ntuple_e2c_paircorrections_jer   = "ntuple_e2c_paircorrections_jer.root";

std::string namef_ntuple_mc_e2c    = "ntuple_mc_e2c.root";
std::string namef_ntuple_mc_at_e2c = "ntuple_mc_at_e2c.root";
std::string namef_ntuple_hadron    = "ntuple_hadron.root";

std::string namef_histos_corr_e2c             = "histos_corr_e2c.root";
std::string namef_histos_paircorr_e2c         = "histos_paircorr_e2c.root";
std::string namef_histos_paircorr_e2c_ct      = "histos_paircorr_e2c_ct.root";
std::string namef_histos_paircorr_e2c_jes     = "histos_paircorr_e2c_jes.root";
std::string namef_histos_paircorr_e2c_prior   = "histos_paircorr_e2c_prior.root";
std::string namef_histos_paircorr_e2c_shapect = "histos_paircorr_e2c_shapect.root";
std::string namef_histos_paircorr_e2c_jer     = "histos_paircorr_e2c_jer.root";

std::string namef_histos_devfromnom_e2c                  = "histos_devfromnom_e2c.root";
std::string namef_histos_paircorr_devfromnom_e2c_prior   = "histos_paircorr_devfromnom_e2c_prior.root";
std::string namef_histos_paircorr_devfromnom_e2c_shapect = "histos_paircorr_devfromnom_e2c_shapect.root";
std::string namef_histos_paircorr_devfromnom_e2c_jes     = "histos_paircorr_devfromnom_e2c_jes.root";
std::string namef_histos_paircorr_devfromnom_e2c_jer     = "histos_paircorr_devfromnom_e2c_jer.root";
std::string namef_histos_paircorr_devfromnom_e2c_ct      = "histos_paircorr_devfromnom_e2c_ct.root";

// About systematics
std::string available_systematics[] = {
                "corr-paradigm", "ct", "jer", "jes", "prior", "shape-ct"
};

std::map<std::string, std::string> systematic_name  = {
        {available_systematics[0],"Correction Paradigm"},
        {available_systematics[1],"Closure Test"},
        {available_systematics[2],"JER"},
        {available_systematics[3],"JES"},
        {available_systematics[4],"Prior Variation"},
        {available_systematics[5],"Shape Closure Test"}
};

std::map<std::string, std::string> systematic_namef = {
        {available_systematics[0],namef_histos_corr_e2c},
        {available_systematics[1],namef_histos_paircorr_e2c_ct},
        {available_systematics[2],namef_histos_paircorr_e2c_jer},
        {available_systematics[3],namef_histos_paircorr_e2c_jes},
        {available_systematics[4],namef_histos_paircorr_e2c_prior},
        {available_systematics[5],namef_histos_paircorr_e2c_shapect}
};

std::map<std::string, std::string> devfromnom_namef = {
        {available_systematics[0],namef_histos_devfromnom_e2c},
        {available_systematics[1],namef_histos_paircorr_devfromnom_e2c_ct},
        {available_systematics[2],namef_histos_paircorr_devfromnom_e2c_jer},
        {available_systematics[3],namef_histos_paircorr_devfromnom_e2c_jes},
        {available_systematics[4],namef_histos_paircorr_devfromnom_e2c_prior},
        {available_systematics[5],namef_histos_paircorr_devfromnom_e2c_shapect}
};

std::map<std::string, std::string> systematic_errtype = {
        {available_systematics[0],"uniform"},
        {available_systematics[1],"normal"},
        {available_systematics[2],"uniform"},
        {available_systematics[3],"uniform"},
        {available_systematics[4],"uniform"},
        {available_systematics[5],"uniform"}
};

//--------------------------------------------------------------------------------------//
const int Nvars_unfold = 26; 

const char* ntuple_unfold_vars = "weight:R_L:h1_eta:h2_eta:h1_y:h2_y:h1_p:h2_p:h1_pt:h2_pt:jet_eta:mum_pt:mum_eta:mup_pt:mup_eta:jet_pt:jet_pt_truth:reco_passed:deltaR_h1:deltaR_h2:h1_y_truth:h2_y_truth:R_L_truth:weight_truth:weight_pt_truth:weight_pt";

//--------------------------------------------------------------------------------------//
const int Nvars_mcreco          = 21; 
const int Nvars_mc              = 21; 
const int Nvars_data            = 17; 
const int Nvars_paircorrdata    = 26; 

const char* ntuple_mcreco_vars       = "weight:R_L:h1_eta:h2_eta:h1_y:h2_y:h1_charge:h2_charge:h1_p:h2_p:h1_pt:h2_pt:jet_pt:jet_eta:weight_pt:mum_pt:mum_eta:mup_pt:mup_eta:h1_pid:h2_pid";
const char* ntuple_data_vars         = "weight:R_L:h1_eta:h2_eta:h1_y:h2_y:h1_p:h2_p:h1_pt:h2_pt:jet_pt:jet_eta:mum_pt:mum_eta:mup_pt:mup_eta:jet_e";
const char* ntuple_paircorrdata_vars = "weight:efficiency:purity:efficiency_relerror:purity_relerror:R_L:h1_eta:h2_eta:h1_y:h2_y:h1_p:h2_p:h1_pt:h2_pt:jet_pt:jet_eta:weight_pt:jet_e:h1_e:h2_e:year:n_reco_ok:n_reco:n_truth_ok:n_truth:eq_charge";
const char* ntuple_mc_vars           = "weight:R_L:h1_eta:h2_eta:h1_y:h2_y:h1_charge:h2_charge:h1_p:h2_p:h1_pt:h2_pt:jet_pt:jet_eta:weight_pt:mum_pt:mum_eta:mup_pt:mup_eta:h1_pid:h2_pid";

//--------------------------------------------------------------------------------------//
const int   Nvars_corrections_mc         = 14;
const int   Nvars_corrections_mcreco     = 26;
const char* ntuple_corrections_mc_vars   = "weight:R_L:h1_eta:h2_eta:h1_y:h2_y:h1_p:h2_p:h1_pt:h2_pt:jet_eta:jet_pt:weight_pt:eq_charge";
const char* ntuple_corrections_reco_vars = "weight:R_L:h1_eta:h2_eta:h1_y:h2_y:h1_p:h2_p:h1_pt:h2_pt:jet_eta:h1_p_truth:h2_p_truth:h1_pt_truth:h2_pt_truth:jet_pt:jet_pt_truth:deltaR_h1:deltaR_h2:h1_y_truth:h2_y_truth:R_L_truth:weight_truth:weight_pt_truth:weight_pt:eq_charge";

//--------------------------------------------------------------------------------------//
// Deprecated vars
const int   Nvars_hadroncorrections_mc         = 8;
const char* ntuple_hadroncorrections_mc_vars   = "h_eta:h_y:h_p:h_pt:jet_pt:jet_eta:jet_e:deltaR_jet_h";
const int   Nvars_hadroncorrections_reco       = 14;
const char* ntuple_hadroncorrections_reco_vars = "h_eta:h_y:h_p:h_pt:jet_pt:jet_eta:jet_e:jet_e_truth:jet_pt_truth:h_y_truth:h_eta_truth:h_p_truth:deltaR_jet_h:key_match";

// const int   Nvars_purity       = 16;
// const char* ntuple_purity_vars = "h_eta:h_y:h_p:h_pt:jet_pt:jet_eta:mum_pt:mum_eta:mup_pt:mup_eta:jet_e:jet_e_truth:h_y_truth:h_eta_truth:deltaR_jet_h:key_match";

const int   Nvars_corrdata       = 29; 
const char* ntuple_corrdata_vars = "weight:efficiency:purity:efficiency_relerror:purity_relerror:R_L:h1_eta:h2_eta:h1_y:h2_y:h1_p:h2_p:h1_pt:h2_pt:jet_pt:jet_eta:weight_pt:jet_e:h1_e:h2_e:year:n_h1_reco_ok:n_h1_reco:n_h1_truth_ok:n_h1_truth:n_h2_reco_ok:n_h2_reco:n_h2_truth_ok:n_h2_truth";

const int Nvars_mc_at             = 4;
const char* ntuple_mc_at_vars     = "weight:weight_pt:R_L:jet_pt";
const char* ntuple_mcreco_at_vars = "weight:weight_pt:R_L:jet_pt";



// JET NTUPLES SPECS
std::string namef_ntuple_jet_corr       = "ntuple_corrjet.root";
std::string namef_ntuple_jet_purity     = "ntuple_jet_purity.root";
std::string namef_ntuple_jet_efficiency = "ntuple_jet_efficiency.root";
std::string namef_ntuple_jes_jer        = "ntuple_jes_jer.root";

std::string name_ntuple_corrjet       = "ntuple_jet";
std::string name_ntuple_jetpurity     = "ntuple_jetpurity";
std::string name_ntuple_jetefficiency = "ntuple_jetefficiency";
std::string name_ntuple_jes_data      = "ntuple_jes_data";
std::string name_ntuple_jes_reco      = "ntuple_jes_reco";

//--------------------------------------------------------------------------------------//
const int Nvars_corrjet         = 18;
const int Nvars_corrjet_nominal = 19;
const int Nvars_jetpurity       = 8;
const int Nvars_jetefficiency   = 7;
const int Nvars_jet_minimal     = 2;

const char* ntuple_jet_vars           = "jet_pt:jet_e:jet_ndtr:jet_efficiency:jet_purity:jet_efficiency_error:jet_purity_error:mup_eff_id:mup_eff_trk:mup_eff_trg:mum_eff_id:mum_eff_trk:mum_eff_trg:year:jet_eta:z_pt:z_eta:z_y";
const char* ntuple_jet_vars_nominal   = "jet_pt:jet_e:jet_ndtr:jet_efficiency:jet_purity:jet_efficiency_error:jet_purity_error:mup_eff_id:mup_eff_trk:mup_eff_trg:mum_eff_id:mum_eff_trk:mum_eff_trg:year:jet_eta:z_pt:z_eta:z_y:jet_ndtr";
const char* ntuple_jetpurity_vars     = "jet_pt:jet_e:jet_ndtr:jet_pt_truth:jet_e_truth:jet_ndtr_truth:deltaR_matchedjets:jet_eta";
const char* ntuple_jetefficiency_vars = "jet_pt_truth:jet_e_truth:jet_ndtr_truth:jet_pt:jet_e:jet_ndtr:jet_eta_truth";
const char* ntuple_jetminimal_vars    = "jet_pt:jet_eta";
//--------------------------------------------------------------------------------------//
const int   Nvars_jes            = 5;
const char* ntuple_jec_reco_vars = "z_pt:jet_pt:jet_eta:jet_jec_cor:jet_jec_err";
const char* ntuple_jec_data_vars = "z_pt:jet_pt:jet_eta:jet_jec_cor:jet_jec_err";
const char* ntuple_jec_mc_vars   = "z_pt:jet_pt:jet_eta:jet_jec_cor:jet_jec_err";

#endif