#ifndef NAMES_H
#define NAMES_H

// Names of the files
std::string namef_ntuple_e2c            = "ntuple_e2c.root";
std::string namef_ntuple_e2c_corr       = "ntuple_corre2c.root";
std::string namef_ntuple_e2c_purity     = "ntuple_e2c_purity.root";
std::string namef_ntuple_e2c_efficiency = "ntuple_e2c_efficiency.root";
std::string namef_ntuple_jetptpurity    = "ntuple_jetpurity.root";

// Name of the ntuples
std::string name_ntuple_unfold          = "ntuple_unfold";
std::string name_ntuple_data            = "ntuple_data";
std::string name_ntuple_mcjetmatch      = "ntuple_mcjetmatch";
std::string name_ntuple_mcreco          = "ntuple_mcreco";
std::string name_ntuple_mc              = "ntuple_mc";
std::string name_ntuple_purity          = "ntuple_purity";
std::string name_ntuple_jetpurity       = "ntuple_jetpurity";
std::string name_ntuple_efficiency_mc   = "ntuple_efficiency_mc";
std::string name_ntuple_efficiency_reco = "ntuple_efficiency_reco";

// TNTuples variables
const int Nvars_unfold          = 18; 
const int Nvars_mcreco          = 22; 
const int Nvars_data            = 22; 
const int Nvars_corrdata        = 26; 
const int Nvars_mc              = 22; 
const int Nvars_purity          = 23;
const int Nvars_efficiency_mc   = 15;
const int Nvars_efficiency_reco = 21;

const char* ntuple_mcreco_vars          = "weight:R_L:h1_eta:h2_eta:h1_y:h2_y:h1_phi:h2_phi:h1_p:h2_p:h1_pt:h2_pt:jet_pt:jet_eta:jet_phi:deltaphi_z_jet:R_L_mum_jet:mum_pt:mum_eta:R_L_mup_jet:mup_pt:mup_eta";
const char* ntuple_data_vars            = "weight:R_L:h1_eta:h2_eta:h1_y:h2_y:h1_phi:h2_phi:h1_p:h2_p:h1_pt:h2_pt:jet_pt:jet_eta:jet_phi:deltaphi_z_jet:R_L_mum_jet:mum_pt:mum_eta:R_L_mup_jet:mup_pt:mup_eta";
const char* ntuple_corrdata_vars        = "weight:efficiency:purity:efficiency_relerror:purity_relerror:R_L:h1_eta:h2_eta:h1_y:h2_y:h1_phi:h2_phi:h1_p:h2_p:h1_pt:h2_pt:jet_pt:jet_eta:jet_phi:deltaphi_z_jet:R_L_mum_jet:mum_pt:mum_eta:R_L_mup_jet:mup_pt:mup_eta";
const char* ntuple_mc_vars              = "weight:R_L:h1_eta:h2_eta:h1_y:h2_y:h1_phi:h2_phi:h1_p:h2_p:h1_pt:h2_pt:jet_pt:jet_eta:jet_phi:deltaphi_z_jet:R_L_mum_jet:mum_pt:mum_eta:R_L_mup_jet:mup_pt:mup_eta";
const char* ntuple_purity_vars          = "h_eta:h_y:h_phi:h_p:h_pt:jet_pt:jet_eta:deltaphi_z_jet_root:deltaphi_z_jet_em:R_L_mum_jet:mum_pt:mum_eta:R_L_mup_jet:mup_pt:mup_eta:recojet_e:truthjet_e:ndtr:htruth_y:htruth_eta:htruth_phi:R_jet_h:key_match";
const char* ntuple_efficiency_reco_vars = "h_eta:h_y:h_phi:h_p:h_pt:jet_pt:jet_eta:deltaphi_z_jet_root:deltaphi_z_jet_em:R_L_mum_jet:mum_pt:mum_eta:R_L_mup_jet:mup_pt:mup_eta:recojet_e:truthjet_e:ndtr:htruth_y:htruth_eta:htruth_phi";
const char* ntuple_efficiency_mc_vars   = "h_eta:h_y:h_phi:h_p:h_pt:jet_pt:jet_eta:deltaphi_z_jet_em:R_L_mum_jet:mum_pt:mum_eta:R_L_mup_jet:mup_pt:mup_eta";
#endif