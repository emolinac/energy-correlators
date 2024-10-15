#ifndef ANALYSIS_CONSTANTS_H
#define ANALYSIS_CONSTANTS_H

#include "TCut.h"
#include "TString.h"

// PID 
const int p_id     = 2212;
const int pbar_id  = -p_id;
const int pip_id   = 211;
const int pim_id   = -pip_id;
const int kp_id    = 321;
const int km_id    = -kp_id;
const int gamma_id = 22;

// Masses (GeV)
const double rho_mass      = 0.77526;  // PDG 2023
const double omega_mass    = 0.78266; // PDG 2023
const double eta_mass      = 0.547862; // PDG 2023
const double etaprime_mass = 0.95778;  // PDG 2023
const double kaonp_mass    = 0.493677; // PDG 2023
const double kaonm_mass    = 0.493677; // PDG 2023
const double kaon_mass     = 0.497611; // PDG 2023
const double pi_mass       = 0.134977;
const double phi_mass      = 1.019455;
const double mass_res      = 0.008; // Mass resolution parameter (see src-resolution)

// Define matching parameter
const double R_match_max = 0.009;
const double R_match_min = 0.; // THIS VALUE CANT BE MOVED!

// Define binning
const double R_L_res  = 0.02;
const int Nbin_R_L    = 15;
const int Nbin_jet_pt = 5;

const double R_L_min    = 0.0099;
const double R_L_max    = 1.2;
const double jet_pt_min = 20; 
const double jet_pt_max = 150;
const double eta_min    = 2.5;
const double eta_max    = 4.;

// Binning
//const double jet_pt_binning[] = {jet_pt_min, 30., 50., jet_pt_max};
//const double jet_pt_binning[] = {jet_pt_min, 26.868, 38.404, jet_pt_max};
const double jet_pt_binning[] = {jet_pt_min, 25.9735, 33.1105, 43.0815, 59.6175, jet_pt_max};

// Visual constants
const double std_marker_size  = 1.3;
const int    std_marker_style = 8;
const int    std_line_width   = 3;

const int std_marker_style_jet_pt[]  = {24,25,26,32,27};
const int corr_marker_style_jet_pt[] = {20,21,22,23,33};

const int std_marker_color_jet_pt[]  = {879,433,417,632,802}; // darker
const int corr_marker_color_jet_pt[] = {880,432,416,633,797}; // lighter

#endif