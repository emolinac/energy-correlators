#ifndef NAMES_H
#define NAMES_H

// Names of the files
std::string namef_ntuple_e2c            = "ntuple_e2c.root";
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
const int Nvars_mc              = 22; 
const int Nvars_purity          = 32;
const int Nvars_efficiency_mc   = 22;
const int Nvars_efficiency_reco = 32;

const char* ntuple_mcreco_vars          = "weight:R_L:h1_eta:h2_eta:h1_y:h2_y:h1_phi:h2_phi:h1_p:h2_p:h1_pt:h2_pt:jet_pt:jet_eta:jet_phi:deltaphi_z_jet:R_L_mum_jet:mum_pt:mum_eta:R_L_mup_jet:mup_pt:mup_eta";
const char* ntuple_data_vars            = "weight:R_L:h1_eta:h2_eta:h1_y:h2_y:h1_phi:h2_phi:h1_p:h2_p:h1_pt:h2_pt:jet_pt:jet_eta:jet_phi:deltaphi_z_jet:R_L_mum_jet:mum_pt:mum_eta:R_L_mup_jet:mup_pt:mup_eta";
const char* ntuple_mc_vars              = "weight:R_L:h1_eta:h2_eta:h1_y:h2_y:h1_phi:h2_phi:h1_p:h2_p:h1_pt:h2_pt:jet_pt:jet_eta:jet_phi:deltaphi_z_jet:R_L_mum_jet:mum_pt:mum_eta:R_L_mup_jet:mup_pt:mup_eta";
const char* ntuple_purity_vars          = "weight:R_L:h1_eta:h2_eta:h1_y:h2_y:h1_phi:h2_phi:h1_p:h2_p:h1_pt:h2_pt:jet_pt:jet_eta:jet_phi:deltaphi_z_jet:R_L_mum_jet:mum_pt:mum_eta:R_L_mup_jet:mup_pt:mup_eta:recojet_e:truthjet_e:ndtr:deltaR_h1:deltaR_h2:h1truth_y:h2truth_y:h1truth_phi:h2truth_phi:R_L_truth";
const char* ntuple_efficiency_reco_vars = "weight:R_L:h1_eta:h2_eta:h1_y:h2_y:h1_phi:h2_phi:h1_p:h2_p:h1_pt:h2_pt:jet_pt:jet_eta:jet_phi:deltaphi_z_jet:R_L_mum_jet:mum_pt:mum_eta:R_L_mup_jet:mup_pt:mup_eta:recojet_e:truthjet_e:ndtr:deltaR_h1:deltaR_h2:h1truth_y:h2truth_y:h1truth_phi:h2truth_phi:R_L_truth";
const char* ntuple_efficiency_mc_vars   = "weight:R_L:h1_eta:h2_eta:h1_y:h2_y:h1_phi:h2_phi:h1_p:h2_p:h1_pt:h2_pt:jet_pt:jet_eta:jet_phi:deltaphi_z_jet:R_L_mum_jet:mum_pt:mum_eta:R_L_mup_jet:mup_pt:mup_eta";
#endif