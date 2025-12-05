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
std::string namef_ntuple_reco2truth_match         = "ntuple_reco2truth_match.root";
std::string namef_ntuple_reco2truth_match_ct      = "ntuple_reco2truth_match_ct.root";
std::string namef_ntuple_reco2truth_match_jes_jer = "ntuple_reco2truth_match_jes_jer.root";

std::string namef_ntuple_truth2reco_match     = "ntuple_truth2reco_match.root";
std::string namef_ntuple_truth2reco_match_ct  = "ntuple_truth2reco_match_ct.root";

std::string namef_3dpaircorr_rl_jetpt_weightpt_histos         = "histos_3dpaircorr_rl_jetpt_weightpt_eec.root";
std::string namef_3dpaircorr_rl_jetpt_weightpt_histos_jes_jer = "histos_3dpaircorr_rl_jetpt_weightpt_eec_jes_jer.root";
std::string namef_3dpaircorr_rl_jetpt_weightpt_histos_ct      = "histos_3dpaircorr_rl_jetpt_weightpt_eec_ct.root";
std::string namef_3dpaircorr_rl_jetpt_weightpt_histos_muon    = "histos_3dpaircorr_rl_jetpt_weightpt_eec_muon.root";

std::string namef_ntuple_eec      = "ntuple_eec.root";
std::string namef_ntuple_eec_corr = "ntuple_correec.root";

std::string namef_ntuple_mc_eec    = "ntuple_mc_eec.root";
std::string namef_ntuple_mc_at_eec = "ntuple_mc_at_eec.root";
std::string namef_ntuple_hadron    = "ntuple_hadron.root";

// About names and options
std::map<std::string, std::string> namef_reco_corrections = {
        {"--get-nominal", namef_ntuple_reco2truth_match},
        {"--get-jes-jer", namef_ntuple_reco2truth_match_jes_jer},
        {"--get-prior"  , namef_ntuple_reco2truth_match},
        {"--get-regpar" , namef_ntuple_reco2truth_match},
        {"--get-muon"   , namef_ntuple_reco2truth_match}
};

std::map<std::string, std::string> namef_all_corrections = {
        {"--get-nominal", namef_3dpaircorr_rl_jetpt_weightpt_histos},
        {"--get-jes-jer", namef_3dpaircorr_rl_jetpt_weightpt_histos_jes_jer},
        {"--get-prior"  , namef_3dpaircorr_rl_jetpt_weightpt_histos},
        {"--get-regpar" , namef_3dpaircorr_rl_jetpt_weightpt_histos},
        {"--get-muon"   , namef_3dpaircorr_rl_jetpt_weightpt_histos_muon}
};

// About systematics
std::string available_systematics[] = {
        "ct-stat", "ct-shape", "jes-jer", "prior", "regpar", "muon"
};

std::map<std::string, std::string> systematic_name  = {
        {available_systematics[0],"Statistical closure test"},
        {available_systematics[1],"Shape closure test"},
        {available_systematics[2],"JES-JER"},
        {available_systematics[3],"Prior variation"},
        {available_systematics[4],"Regularization parameter"},
        {available_systematics[5],"Muon eff"}
};

std::map<std::string, std::string> systematic_namef = {
        {available_systematics[0],"histos_eec_3dcorr_rl_jetpt_weightpt_niter4_niterjet4_statct_niterct10.root"},
        {available_systematics[1],"histos_eec_3dcorr_rl_jetpt_weightpt_niter4_niterjet4_shapect.root"},
        {available_systematics[2],"histos_eec_3dcorr_rl_jetpt_weightpt_niter4_niterjet4--get-jes-jer.root"},
        {available_systematics[3],"histos_eec_3dcorr_rl_jetpt_weightpt_niter4_niterjet4--get-prior.root"},
        {available_systematics[4],"histos_eec_3dcorr_rl_jetpt_weightpt_niterjet4--get-regpar.root"},
        {available_systematics[5],"histos_eec_3dcorr_rl_jetpt_weightpt_niter4_niterjet4--get-muon.root"}
};

std::map<std::string, std::string> systematic_errtype = {
        {available_systematics[0],"normal"},
        {available_systematics[1],"normal"},
        {available_systematics[2],"normal"},
        {available_systematics[3],"normal"},
        {available_systematics[4],"normal"},
        {available_systematics[5],"normal"}
};

std::string namef_ntuple_jes_jer = "ntuple_jes_jer.root";

std::string name_ntuple_corrjet              = "ntuple_jet";
std::string name_ntuple_jet_reco2truth_match = "ntuple_jetpurity";
std::string name_ntuple_jet_truth2reco_match = "ntuple_jetefficiency";
std::string name_ntuple_jes_data             = "ntuple_jes_data";
std::string name_ntuple_jes_reco             = "ntuple_jes_reco";


//--------------------------------------------------------------------------------------//
const int Nvars_mcreco          = 21; 
const int Nvars_mc              = 21; 
const int Nvars_mc_match        = 22; 
const int Nvars_data            = 17; 
const int Nvars_paircorrdata    = 26; 

const char* ntuple_mcreco_vars       = "weight:R_L:h1_eta:h2_eta:h1_y:h2_y:h1_charge:h2_charge:h1_p:h2_p:h1_pt:h2_pt:jet_pt:jet_eta:weight_pt:mum_pt:mum_eta:mup_pt:mup_eta:h1_pid:h2_pid";
const char* ntuple_data_vars         = "event_weight:R_L:h1_eta:h2_eta:h1_y:h2_y:h1_p:h2_p:h1_pt:h2_pt:jet_pt:jet_eta:mum_pt:mum_eta:mup_pt:mup_eta:jet_e";
const char* ntuple_paircorrdata_vars = "event_weight:efficiency:purity:efficiency_relerror:purity_relerror:R_L:h1_eta:h2_eta:h1_y:h2_y:h1_p:h2_p:h1_pt:h2_pt:jet_pt:jet_eta:weight_pt:jet_e:h1_e:h2_e:year:n_reco_ok:n_reco:n_truth_ok:n_truth:eq_charge";
const char* ntuple_mc_vars           = "eq_charge:R_L:h1_eta:h2_eta:h1_y:h2_y:h1_charge:h2_charge:h1_p:h2_p:h1_pt:h2_pt:jet_pt:jet_eta:weight_pt:mum_pt:mum_eta:mup_pt:mup_eta:h1_pid:h2_pid";
const char* ntuple_mc_match_vars     = "eq_charge:R_L:h1_eta:h2_eta:h1_y:h2_y:h1_charge:h2_charge:h1_p:h2_p:h1_pt:h2_pt:jet_pt:jet_eta:weight_pt:mum_pt:mum_eta:mup_pt:mup_eta:h1_pid:h2_pid:weight_pt_reco";

const int Nvars_pairandsinglecorrdata = 28;
const char* ntuple_pairandsinglecorrdata_vars = "event_weight:efficiency:purity:h1_efficiency:h2_efficiency:h1_purity:h2_purity:R_L:h1_eta:h2_eta:h1_y:h2_y:h1_p:h2_p:h1_pt:h2_pt:jet_pt:jet_eta:weight_pt:jet_e:h1_e:h2_e:year:n_reco_ok:n_reco:n_truth_ok:n_truth:eq_charge";

//--------------------------------------------------------------------------------------//
const int   Nvars_corrections_mc         = 13;
const int   Nvars_corrections_mcreco     = 25;
const int   Nvars_corrections_mcreco_plus= 26;
const char* ntuple_corrections_mc_vars   = "eq_charge:R_L:h1_eta:h2_eta:h1_y:h2_y:h1_p:h2_p:h1_pt:h2_pt:jet_eta:jet_pt:weight_pt";
const char* ntuple_corrections_reco_vars = "eq_charge:R_L:h1_eta:h2_eta:h1_y:h2_y:h1_p:h2_p:h1_pt:h2_pt:jet_eta:h1_p_truth:h2_p_truth:h1_pt_truth:h2_pt_truth:jet_pt:jet_pt_truth:deltaR_h1:deltaR_h2:h1_eta_truth:h2_eta_truth:R_L_truth:weight_truth:weight_pt_truth:weight_pt";
const char* ntuple_corrections_reco_vars_plus = "eq_charge:R_L:h1_eta:h2_eta:h1_y:h2_y:h1_p:h2_p:h1_pt:h2_pt:jet_eta:h1_p_truth:h2_p_truth:h1_pt_truth:h2_pt_truth:jet_pt:jet_pt_truth:deltaR_h1:deltaR_h2:h1_eta_truth:h2_eta_truth:R_L_truth:weight_truth:weight_pt_truth:weight_pt:jet_eta_truth";

//--------------------------------------------------------------------------------------//
const int   Nvars_hadroncorrections_mc         = 8;
const int   Nvars_hadroncorrections_reco       = 14;
const char* ntuple_hadroncorrections_mc_vars   = "h_eta:h_y:h_p:h_pt:jet_pt:jet_eta:jet_e:deltaR_jet_h";
const char* ntuple_hadroncorrections_reco_vars = "h_eta:h_y:h_p:h_pt:jet_pt:jet_eta:jet_e:jet_e_truth:jet_pt_truth:h_y_truth:h_eta_truth:h_p_truth:deltaR_jet_h:key_match";

const int   Nvars_corrdata       = 30; 
const char* ntuple_corrdata_vars = "event_weight:efficiency:purity:efficiency_relerror:purity_relerror:R_L:h1_eta:h2_eta:h1_y:h2_y:h1_p:h2_p:h1_pt:h2_pt:jet_pt:jet_eta:weight_pt:jet_e:h1_e:h2_e:year:n_h1_reco_ok:n_h1_reco:n_h1_truth_ok:n_h1_truth:n_h2_reco_ok:n_h2_reco:n_h2_truth_ok:n_h2_truth:eq_charge";

//--------------------------------------------------------------------------------------//
const int Nvars_corrjet         = 18;
const int Nvars_jetpurity       = 8;
const int Nvars_jetefficiency   = 7;
const int Nvars_jet_minimal     = 2;

const char* ntuple_jet_vars           = "jet_pt:jet_e:jet_ndtr:jet_efficiency:jet_purity:jet_efficiency_error:jet_purity_error:mup_eff_id:mup_eff_trk:mup_eff_trg:mum_eff_id:mum_eff_trk:mum_eff_trg:year:jet_eta:z_pt:z_eta:z_y";
const char* ntuple_jetpurity_vars     = "jet_pt:jet_e:jet_ndtr:jet_pt_truth:jet_e_truth:jet_ndtr_truth:deltaR_matchedjets:jet_eta";
const char* ntuple_jetefficiency_vars = "jet_pt_truth:jet_e_truth:jet_ndtr_truth:jet_pt:jet_e:jet_ndtr:jet_eta_truth";
const char* ntuple_jetminimal_vars    = "jet_pt:jet_eta";
//--------------------------------------------------------------------------------------//
const int   Nvars_jes            = 5;
const char* ntuple_jec_reco_vars = "z_pt:jet_pt:jet_eta:jet_jec_cor:jet_jec_err";
const char* ntuple_jec_data_vars = "z_pt:jet_pt:jet_eta:jet_jec_cor:jet_jec_err";
const char* ntuple_jec_mc_vars   = "z_pt:jet_pt:jet_eta:jet_jec_cor:jet_jec_err";

#endif