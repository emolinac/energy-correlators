#ifndef NAMES_H
#define NAMES_H

// Names of the files
std::string namef_ntuple_e2c        = "ntuple_e2c.root";
std::string namef_ntuple_e2c_purity = "ntuple_e2c_purity.root";

// Name of the ntuples
std::string name_ntuple_unfold     = "ntuple_unfold";
std::string name_ntuple_data       = "ntuple_data";
std::string name_ntuple_mcjetmatch = "ntuple_mcjetmatch";
std::string name_ntuple_mcreco     = "ntuple_mcreco";
std::string name_ntuple_mc         = "ntuple_mc";
std::string name_ntuple_purity     = "ntuple_purity";

// TNTuples variables
const int Nvars_unfold     = 18; 
const int Nvars_mcjetmatch = 23; 
const int Nvars_mcreco     = 44; 
const int Nvars_data       = 44; 
const int Nvars_mc         = 36; 
const int Nvars_purity     = 46;

const char* ntuple_mcjetmatch_vars = "weight:R_L:h1_id:h2_id:h1_eta:h2_eta:h1_y:h2_y:h1_phi:h2_phi:h1_motherid:h2_motherid:h1_topmotherid:h2_topmotherid:h1_p:h2_p:h1_pt:h2_pt:h1_pz:h2_pz:jet_pt:jet_eta:jet_phi";
const char* ntuple_mcreco_vars     = "weight:R_L:h1_id:h2_id:h1_eta:h2_eta:h1_y:h2_y:h1_phi:h2_phi:h1_chi2:h2_chi2:h1_ndf:h2_ndf:h1_probnnghost:h2_probnnghost:h1_p:h2_p:h1_pt:h2_pt:h1_pz:h2_pz:jet_pt:jet_eta:jet_phi:z0_phi:mum_phi:mum_pt:mum_eta:mum_px:mum_py:mum_pz:mum_pe:mum_m:mum_probchi2:mup_phi:mup_pt:mup_eta:mup_px:mup_py:mup_pz:mup_pe:mup_m:mup_probchi2";
const char* ntuple_data_vars       = "weight:R_L:h1_id:h2_id:h1_eta:h2_eta:h1_y:h2_y:h1_phi:h2_phi:h1_chi2:h2_chi2:h1_ndf:h2_ndf:h1_probnnghost:h2_probnnghost:h1_p:h2_p:h1_pt:h2_pt:h1_pz:h2_pz:jet_pt:jet_eta:jet_phi:z0_phi:mum_phi:mum_pt:mum_eta:mum_px:mum_py:mum_pz:mum_pe:mum_m:mum_probchi2:mup_phi:mup_pt:mup_eta:mup_px:mup_py:mup_pz:mup_pe:mup_m:mup_probchi2";
const char* ntuple_mc_vars         = "weight:R_L:h1_id:h2_id:h1_eta:h2_eta:h1_y:h2_y:h1_phi:h2_phi:h1_p:h2_p:h1_pt:h2_pt:h1_pz:h2_pz:jet_pt:jet_eta:jet_phi:z0_phi:mum_phi:mum_pt:mum_eta:mum_px:mum_py:mum_pz:mum_pe:mum_m:mup_phi:mup_pt:mup_eta:mup_px:mup_py:mup_pz:mup_pe:mup_m";
const char* ntuple_purity_vars     = "weight:R_L:h1_id:h2_id:h1_eta:h2_eta:h1_y:h2_y:h1_phi:h2_phi:h1_chi2:h2_chi2:h1_ndf:h2_ndf:h1_probnnghost:h2_probnnghost:h1_p:h2_p:h1_pt:h2_pt:h1_pz:h2_pz:jet_pt:jet_eta:jet_phi:z0_phi:mum_phi:mum_pt:mum_eta:mum_px:mum_py:mum_pz:mum_pe:mum_m:mum_probchi2:mup_phi:mup_pt:mup_eta:mup_px:mup_py:mup_pz:mup_pe:mup_m:mup_probchi2:signal_h1:signal_h2";

#endif