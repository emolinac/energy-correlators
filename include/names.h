#ifndef NAMES_H
#define NAMES_H

// JET NTUPLES SPECS
std::string namef_ntuple_jet_purity     = "ntuple_jet_purity.root";
std::string namef_ntuple_jet_efficiency = "ntuple_jet_efficiency.root";

std::string name_ntuple_jetpurity     = "ntuple_jetpurity";
std::string name_ntuple_jetefficiency = "ntuple_jetefficiency";

const int Nvars_jetpurity     = 6;
const int Nvars_jetefficiency = 6;
const char* ntuple_jetpurity_vars  = "jet_pt:jet_e:jet_ndtr:jet_pt_truth:jet_e_truth:jet_ndtr_truth";
const char* ntuple_jetefficiency_vars  = "jet_pt_truth:jet_e_truth:jet_ndtr_truth:jet_pt:jet_e:jet_ndtr";

// HADRONIC NTUPLES SPECS

// Names of the files
std::string namef_ntuple_e2c            = "ntuple_e2c.root";
std::string namef_ntuple_e2c_corr       = "ntuple_corre2c.root";
std::string namef_ntuple_e2c_purity     = "ntuple_e2c_purity.root";
std::string namef_ntuple_e2c_pairpurity = "ntuple_e2c_pairpurity.root";
std::string namef_ntuple_e2c_unfolding  = "ntuple_e2c_unfolding.root";
std::string namef_ntuple_e2c_efficiency = "ntuple_e2c_efficiency.root";
std::string namef_ntuple_jetptpurity    = "ntuple_jetpurity.root";

// Name of the ntuples
std::string name_ntuple_unfold          = "ntuple_unfold";
std::string name_ntuple_data            = "ntuple_data";
std::string name_ntuple_mcjetmatch      = "ntuple_mcjetmatch";
std::string name_ntuple_mcreco          = "ntuple_mcreco";
std::string name_ntuple_mc              = "ntuple_mc";
std::string name_ntuple_purity          = "ntuple_purity";
std::string name_ntuple_unfolding       = "ntuple_unfolding";
std::string name_ntuple_efficiency_mc   = "ntuple_efficiency_mc";
std::string name_ntuple_efficiency_reco = "ntuple_efficiency_reco";

// TNTuples variables
const int Nvars_unfold          = 18; 
const int Nvars_mcreco          = 22; 
const int Nvars_data            = 23; 
const int Nvars_corrdata        = 23; 
const int Nvars_mc              = 22; 
const int Nvars_purity          = 23;
const int Nvars_efficiency_mc   = 14;
const int Nvars_efficiency_reco = 21;
const int Nvars_pairpurity      = 32;

const char* ntuple_mcreco_vars          = "weight:R_L:h1_eta:h2_eta:h1_y:h2_y:h1_phi:h2_phi:h1_p:h2_p:h1_pt:h2_pt:jet_pt:jet_eta:jet_phi:deltaphi_z_jet:R_L_mum_jet:mum_pt:mum_eta:R_L_mup_jet:mup_pt:mup_eta";
const char* ntuple_data_vars            = "weight:R_L:h1_eta:h2_eta:h1_y:h2_y:h1_phi:h2_phi:h1_p:h2_p:h1_pt:h2_pt:jet_pt:jet_eta:jet_phi:deltaphi_z_jet:R_L_mum_jet:mum_pt:mum_eta:R_L_mup_jet:mup_pt:mup_eta:jet_e";
const char* ntuple_corrdata_vars        = "weight:efficiency:purity:efficiency_relerror:purity_relerror:R_L:h1_eta:h2_eta:h1_y:h2_y:h1_phi:h2_phi:h1_p:h2_p:h1_pt:h2_pt:jet_pt:jet_eta:jet_phi:deltaphi_z_jet:jet_e:h1_e:h2_e";
const char* ntuple_mc_vars              = "weight:R_L:h1_eta:h2_eta:h1_y:h2_y:h1_phi:h2_phi:h1_p:h2_p:h1_pt:h2_pt:jet_pt:jet_eta:jet_phi:deltaphi_z_jet:R_L_mum_jet:mum_pt:mum_eta:R_L_mup_jet:mup_pt:mup_eta";

const char* ntuple_purity_vars          = "h_eta:h_y:h_phi:h_p:h_pt:jet_pt:jet_eta:deltaphi_z_jet_root:deltaphi_z_jet_em:R_L_mum_jet:mum_pt:mum_eta:R_L_mup_jet:mup_pt:mup_eta:recojet_e:truthjet_e:ndtr:h_y_truth:h_eta_truth:h_phi_truth:R_jet_h:key_match";
const char* ntuple_pairpurity_vars      = "weight:R_L:h1_eta:h2_eta:h1_y:h2_y:h1_phi:h2_phi:h1_p:h2_p:h1_pt:h2_pt:jet_eta:jet_phi:deltaphi_z_jet:R_L_mum_jet:mum_pt:mum_eta:R_L_mup_jet:mup_pt:mup_eta:jet_pt:jet_pt_truth:ndtr:deltaR_h1:deltaR_h2:h1truth_y:h2truth_y:h1truth_phi:h2truth_phi:R_L_truth:weight_truth";

const char* ntuple_efficiency_reco_vars = "h_eta:h_y:h_phi:h_p:h_pt:jet_pt:jet_eta:deltaphi_z_jet_root:deltaphi_z_jet_em:R_L_mum_jet:R_L_mup_jet:recojet_e:truthjet_e:jet_pt_truth:ndtr:h_y_truth:h_eta_truth:h_phi_truth:h_p_truth:R_jet_h:key_match";
const char* ntuple_efficiency_mc_vars   = "h_eta:h_y:h_phi:h_p:h_pt:jet_pt:jet_eta:deltaphi_z_jet_root:deltaphi_z_jet_em:R_L_mum_jet:R_L_mup_jet:jet_e:ndtr:R_jet_h";
#endif