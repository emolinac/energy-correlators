#ifndef NAMES_H
#define NAMES_H

// Names of the files
std::string namef_ntuple_dihadron    = "ntuple_dihadron.root";
std::string namef_ntuple_effmuons    = "ntuple_effmuons.root";
std::string namef_ntuple_purity      = "ntuple_purity.root";
std::string namef_ntuple_unfold      = "ntuple_unfold.root";
std::string namef_ntuple_resolution  = "ntuple_resolution.root";
std::string namef_ntuple_decays      = "ntuple_decays.root";
std::string namef_ntuple_datadecays  = "ntuple_datadecays.root";
std::string namef_ntuple_invmass     = "ntuple_invmass.root";
std::string namef_ntuple_invmasstrio = "ntuple_invmasstrio.root";

// Name of the ntuples
std::string name_ntuple_data       = "ntuple_data";
std::string name_ntuple_mc         = "ntuple_mc";
std::string name_ntuple_mcreco     = "ntuple_mcreco";
std::string name_ntuple_purity     = "ntuple_purity";
std::string name_ntuple_unfold     = "ntuple_unfold";
std::string name_ntuple_resolution = "ntuple_resolution";
std::string name_ntuple_decays     = "ntuple_decays";
std::string name_ntuple_datadecays = "ntuple_datadecays";
std::string name_ntuple_invmass    = "ntuple_invmass";
std::string name_ntuple_invmasstrio= "ntuple_invmasstrio";
std::string name_ntuple_muons      = "ntuple_muons";
std::string name_ntuple_reco_muons = "ntuple_reco_muons";

// TNTuples variables
const int Nvars_mc         = 32; 
const int Nvars_mcreco     = 40; 
const int Nvars_data       = 41; 
const int Nvars_purity     = 39;
const int Nvars_unfold     = 41;
const int Nvars_resolution = 44;
const int Nvars_decays     = 20; 
const int Nvars_datadecays = 41; 
const int Nvars_invmass    = 41; 
const int Nvars_invmasstrio= 42; 
const int Nvars_muons      = 34; 

const char* ntuple_mc_vars          = "eq_charge:lh_id:nlh_id:lh_p:nlh_p:lh_pt:nlh_pt:lh_z:nlh_z:dh_kt:lh_pz:nlh_pz:jet_pt:jet_eta:jet_phi:z0_phi:mum_phi:mum_pt:mum_eta:mum_px:mum_py:mum_pz:mum_pe:mum_m:mup_phi:mup_pt:mup_eta:mup_px:mup_py:mup_pz:mup_pe:mup_m";
const char* ntuple_mcreco_vars      = "eq_charge:lh_id:nlh_id:lh_chi2:nlh_chi2:lh_ndf:nlh_ndf:lh_probnnghost:nlh_probnnghost:lh_p:nlh_p:lh_pt:nlh_pt:lh_z:nlh_z:dh_kt:lh_pz:nlh_pz:jet_pt:jet_eta:jet_phi:z0_phi:mum_phi:mum_pt:mum_eta:mum_px:mum_py:mum_pz:mum_pe:mum_m:mum_probchi2:mup_phi:mup_pt:mup_eta:mup_px:mup_py:mup_pz:mup_pe:mup_m:mup_probchi2";
const char* ntuple_data_vars        = "eq_charge:dh_m:lh_id:nlh_id:lh_chi2:nlh_chi2:lh_ndf:nlh_ndf:lh_probnnghost:nlh_probnnghost:lh_p:nlh_p:lh_pt:nlh_pt:lh_z:nlh_z:dh_kt:lh_pz:nlh_pz:jet_pt:jet_eta:jet_phi:z0_phi:mum_phi:mum_pt:mum_eta:mum_px:mum_py:mum_pz:mum_pe:mum_m:mum_probchi2:mup_phi:mup_pt:mup_eta:mup_px:mup_py:mup_pz:mup_pe:mup_m:mup_probchi2";
const char* ntuple_purity_vars      = "eq_charge:signal:lh_chi2:nlh_chi2:lh_ndf:nlh_ndf:lh_probnnghost:nlh_probnnghost:lh_p:nlh_p:lh_pt:nlh_pt:lh_z:nlh_z:dh_kt:lh_pz:nlh_pz:jet_pt:jet_eta:jet_phi:z0_phi:mum_phi:mum_pt:mum_eta:mum_px:mum_py:mum_pz:mum_pe:mum_m:mum_probchi2:mup_phi:mup_pt:mup_eta:mup_px:mup_py:mup_pz:mup_pe:mup_m:mup_probchi2";
const char* ntuple_unfold_vars      = "eq_charge:lh_chi2:nlh_chi2:lh_ndf:nlh_ndf:lh_probnnghost:nlh_probnnghost:lh_p:nlh_p:lh_pt:nlh_pt:nlh_z_mcreco:nlh_z_mc:dh_kt_mcreco:dh_kt_mc:lh_pz:nlh_pz:jet_pt_mcreco:jet_pt_mc:jet_eta:jet_phi:z0_phi:mum_phi:mum_pt:mum_eta:mum_px:mum_py:mum_pz:mum_pe:mum_m:mum_probchi2:mup_phi:mup_pt:mup_eta:mup_px:mup_py:mup_pz:mup_pe:mup_m:mup_probchi2:signal";
const char* ntuple_resolution_vars  = "eq_charge:lh_chi2:nlh_chi2:lh_ndf:nlh_ndf:lh_probnnghost:nlh_probnnghost:nlh_z_mcreco:nlh_z_mc:dh_kt_mcreco:dh_kt_mc:dh_m_mcreco:dh_m_mc:lh_m_mcreco:lh_m_mc:nlh_m_mcreco:nlh_m_mc:jet_pt_mcreco:jet_pt_mc:jet_eta:jet_phi:z0_phi:mum_phi:mum_pt:mum_eta:mum_px:mum_py:mum_pz:mum_pe:mum_m:mum_probchi2:mup_phi:mup_pt:mup_eta:mup_px:mup_py:mup_pz:mup_pe:mup_m:mup_probchi2:lh_p:nlh_p:lh_pt:nlh_pt";
const char* ntuple_decays_vars      = "eq_charge:prob:lh_id:nlh_id:lh_motherid:nlh_motherid:lh_topmotherid:nlh_topmotherid:lh_p:nlh_p:lh_pt:nlh_pt:lh_z:nlh_z:dh_kt:lh_pz:nlh_pz:jet_pt:jet_eta:jet_phi";
const char* ntuple_invmass_vars     = "eq_charge:pair_m:dhcomp_id:h_id:lh_chi2:nlh_chi2:lh_ndf:nlh_ndf:lh_probnnghost:nlh_probnnghost:lh_p:nlh_p:lh_pt:nlh_pt:lh_z:nlh_z:dh_kt:lh_pz:nlh_pz:jet_pt:jet_eta:jet_phi:z0_phi:mum_phi:mum_pt:mum_eta:mum_px:mum_py:mum_pz:mum_pe:mum_m:mum_probchi2:mup_phi:mup_pt:mup_eta:mup_px:mup_py:mup_pz:mup_pe:mup_m:mup_probchi2";
const char* ntuple_datadecays_vars  = "eq_charge:prob:ha_id:hb_id:lh_chi2:nlh_chi2:lh_ndf:nlh_ndf:lh_probnnghost:nlh_probnnghost:lh_p:nlh_p:lh_pt:nlh_pt:lh_z:nlh_z:dh_kt:lh_pz:nlh_pz:jet_pt:jet_eta:jet_phi:z0_phi:mum_phi:mum_pt:mum_eta:mum_px:mum_py:mum_pz:mum_pe:mum_m:mum_probchi2:mup_phi:mup_pt:mup_eta:mup_px:mup_py:mup_pz:mup_pe:mup_m:mup_probchi2";
const char* ntuple_invmasstrio_vars = "eq_charge:trio_m:ha_id:hb_id:hc_id:lh_chi2:nlh_chi2:lh_ndf:nlh_ndf:lh_probnnghost:nlh_probnnghost:lh_p:nlh_p:lh_pt:nlh_pt:lh_z:nlh_z:dh_kt:lh_pz:nlh_pz:jet_pt:jet_eta:jet_phi:z0_phi:mum_phi:mum_pt:mum_eta:mum_px:mum_py:mum_pz:mum_pe:mum_m:mum_probchi2:mup_phi:mup_pt:mup_eta:mup_px:mup_py:mup_pz:mup_pe:mup_m:mup_probchi2";
const char* ntuple_muons_vars       = "eq_charge:dh_kt:lh_pz:nlh_pz:jet_pt:jet_eta:jet_phi:z0_phi:mum_phi:mum_pt:mum_eta:mum_px:mum_py:mum_pz:mum_pe:mum_m:mum_probchi2:mum_ismuon:mum_probnnmuon:mum_ipchi2:mum_pvchi2ndf:mup_phi:mup_pt:mup_eta:mup_px:mup_py:mup_pz:mup_pe:mup_m:mup_probchi2:mup_ismuon:mup_probnnmuon:mup_ipchi2:mup_pvchi2ndf";
#endif