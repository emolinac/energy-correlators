#ifndef NAMES_H
#define NAMES_H

// HADRONIC NTUPLES SPECS

// Names of the files
std::string namef_ntuple_e2c            = "ntuple_e2c.root";
std::string namef_ntuple_e2c_corr       = "ntuple_corre2c.root";
std::string namef_ntuple_e2c_purity     = "ntuple_e2c_purity.root";
std::string namef_ntuple_e2c_pairpurity = "ntuple_e2c_pairpurity.root";
std::string namef_ntuple_e2c_unfolding  = "ntuple_e2c_unfolding.root";
std::string namef_ntuple_e2c_efficiency = "ntuple_e2c_efficiency.root";
std::string namef_ntuple_jetptpurity    = "ntuple_jetpurity.root";
std::string namef_ntuple_mc_e2c         = "ntuple_mc_e2c.root";
std::string namef_ntuple_mc_at_e2c      = "ntuple_mc_at_e2c.root";

// Name of the ntuples
std::string name_ntuple_unfold          = "ntuple_unfold";
std::string name_ntuple_data            = "ntuple_hadron";
std::string name_ntuple_mcjetmatch      = "ntuple_mcjetmatch";
std::string name_ntuple_mcreco          = "ntuple_mcreco";
std::string name_ntuple_mc              = "ntuple_mc";
std::string name_ntuple_purity          = "ntuple_purity";
std::string name_ntuple_unfolding       = "ntuple_unfolding";
std::string name_ntuple_efficiency_mc   = "ntuple_efficiency_mc";
std::string name_ntuple_efficiency_reco = "ntuple_efficiency_reco";
std::string name_ntuple_mc_jet          = "ntuple_mc_jet";
std::string name_ntuple_mcreco_jet      = "ntuple_mcreco_jet";

// TNTuples variables
const int Nvars_unfold          = 34; 
const int Nvars_mcreco          = 24; 
const int Nvars_data            = 23; 
const int Nvars_corrdata        = 24; 
const int Nvars_mc              = 24; 
const int Nvars_purity          = 23;
const int Nvars_efficiency_mc   = 14;
const int Nvars_efficiency_reco = 21;
const int Nvars_pairpurity      = 34;
const int Nvars_mc_at           = 4;

const char* ntuple_mc_at_vars     = "weight:weight_pt:R_L:jet_pt";
const char* ntuple_mcreco_at_vars = "weight:weight_pt:R_L:jet_pt";

const char* ntuple_mcreco_vars          = "weight:R_L:h1_eta:h2_eta:h1_y:h2_y:h1_charge:h2_charge:h1_p:h2_p:h1_pt:h2_pt:jet_pt:jet_eta:weight_pt:deltaphi_z_jet:R_L_mum_jet:mum_pt:mum_eta:R_L_mup_jet:mup_pt:mup_eta:h1_pid:h2_pid";
const char* ntuple_data_vars            = "weight:R_L:h1_eta:h2_eta:h1_y:h2_y:h1_phi:h2_phi:h1_p:h2_p:h1_pt:h2_pt:jet_pt:jet_eta:jet_phi:deltaphi_z_jet:R_L_mum_jet:mum_pt:mum_eta:R_L_mup_jet:mup_pt:mup_eta:jet_e";
const char* ntuple_corrdata_vars        = "weight:efficiency:purity:efficiency_relerror:purity_relerror:R_L:h1_eta:h2_eta:h1_y:h2_y:h1_phi:h2_phi:h1_p:h2_p:h1_pt:h2_pt:jet_pt:jet_eta:jet_phi:weight_pt:jet_e:h1_e:h2_e:year";
const char* ntuple_mc_vars              = "weight:R_L:h1_eta:h2_eta:h1_y:h2_y:h1_charge:h2_charge:h1_p:h2_p:h1_pt:h2_pt:jet_pt:jet_eta:weight_pt:deltaphi_z_jet:R_L_mum_jet:mum_pt:mum_eta:R_L_mup_jet:mup_pt:mup_eta:h1_pid:h2_pid";

const char* ntuple_purity_vars          = "h_eta:h_y:h_phi:h_p:h_pt:jet_pt:jet_eta:deltaphi_z_jet_root:deltaphi_z_jet_em:R_L_mum_jet:mum_pt:mum_eta:R_L_mup_jet:mup_pt:mup_eta:recojet_e:truthjet_e:ndtr:h_y_truth:h_eta_truth:h_phi_truth:R_jet_h:key_match";
const char* ntuple_pairpurity_vars      = "weight:R_L:h1_eta:h2_eta:h1_y:h2_y:h1_phi:h2_phi:h1_p:h2_p:h1_pt:h2_pt:jet_eta:jet_phi:deltaphi_z_jet:R_L_mum_jet:mum_pt:mum_eta:R_L_mup_jet:mup_pt:mup_eta:jet_pt:jet_pt_truth:ndtr:deltaR_h1:deltaR_h2:h1truth_y:h2truth_y:h1truth_phi:h2truth_phi:R_L_truth:weight_truth:weight_pt_truth:weight_pt";
const char* ntuple_unfold_vars          = "weight:R_L:h1_eta:h2_eta:h1_y:h2_y:h1_phi:h2_phi:h1_p:h2_p:h1_pt:h2_pt:jet_eta:jet_phi:deltaphi_z_jet:R_L_mum_jet:mum_pt:mum_eta:R_L_mup_jet:mup_pt:mup_eta:jet_pt:jet_pt_truth:reco_passed:deltaR_h1:deltaR_h2:h1truth_y:h2truth_y:h1truth_phi:h2truth_phi:R_L_truth:weight_truth:weight_pt_truth:weight_pt";

const char* ntuple_efficiency_reco_vars = "h_eta:h_y:h_phi:h_p:h_pt:jet_pt:jet_eta:deltaphi_z_jet_root:deltaphi_z_jet_em:R_L_mum_jet:R_L_mup_jet:recojet_e:truthjet_e:jet_pt_truth:ndtr:h_y_truth:h_eta_truth:h_phi_truth:h_p_truth:R_jet_h:key_match";
const char* ntuple_efficiency_mc_vars   = "h_eta:h_y:h_phi:h_p:h_pt:jet_pt:jet_eta:deltaphi_z_jet_root:deltaphi_z_jet_em:R_L_mum_jet:R_L_mup_jet:jet_e:ndtr:R_jet_h";


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
std::string name_ntuple_jes_mc        = "ntuple_jes_mc";
std::string name_ntuple_jer           = "ntuple_jer";

const int Nvars_corrjet       = 14;
const int Nvars_jetpurity     = 7;
const int Nvars_jetefficiency = 6;
const int Nvars_jes_reco      = 3;
const int Nvars_jes           = 2;

const char* ntuple_jet_vars           = "jet_pt:jet_e:jet_ndtr:jet_efficiency:jet_purity:jet_efficiency_error:jet_purity_error:mup_eff_id:mup_eff_trk:mup_eff_trg:mum_eff_id:mum_eff_trk:mum_eff_trg:year";
const char* ntuple_jetpurity_vars     = "jet_pt:jet_e:jet_ndtr:jet_pt_truth:jet_e_truth:jet_ndtr_truth:deltaR_matchedjets";
const char* ntuple_jetefficiency_vars = "jet_pt_truth:jet_e_truth:jet_ndtr_truth:jet_pt:jet_e:jet_ndtr";
const char* ntuple_jes_reco_vars      = "jet_pt:z_pt:jet_jes_cor";
const char* ntuple_jes_mc_vars        = "jet_pt:z_pt";
const char* ntuple_jes_data_vars      = "jet_pt:z_pt";

// MUON NTUPLE SPECS

std::string namef_ntuple_muon = "ntuple_muons.root";
std::string name_ntuple_muon  = "ntuple_muons";

const int Nvars_muons = 10;

const char* ntuple_muons_vars = "mup_eff_id:mup_eff_trk:mup_eff_trg:mum_eff_id:mum_eff_trk:mum_eff_trg";

#endif