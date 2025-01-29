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
const double R_L_res  = 0.0455;
// const double R_L_res  = 3;
const int Nbin_R_L    = 15;
const int Nbin_jet_pt = 3;
const int Nbin_weight = 10;

const double corr_rel_error = 0.1;
const int ndim_corr = 50;

const double R_L_min        = 0.0099;
const double R_L_max        = 0.49;
const double jet_pt_min_nom = 20; 
const double jet_pt_max     = 100;
const double eta_min        = 2.;
const double eta_max        = 4.5;

// Binning
const double jet_pt_binning[] = {jet_pt_min_nom, 30., 50., jet_pt_max};
const double weight_binning[] = {0, 0.000366, 0.000686, 0.001082, 0.00159, 0.002286, 0.003274, 0.004818, 0.007506, 0.013538, 0.4};

// Visual constants
const double std_marker_size  = 1.0;
const int    std_marker_style = 8;
const int    std_line_width   = 3;

const int std_marker_style_jet_pt[]  = {24,25,26,32,27};
const int corr_marker_style_jet_pt[] = {20,21,22,23,33};

const int std_marker_color_jet_pt[]  = {879,433,417,632,802}; // darker
const int corr_marker_color_jet_pt[] = {880,432,416,633,797}; // lighter

const double em_rl_binning[] = {R_L_min, 0.0419067, 0.0739133, 0.10592, 0.137927, 0.169933, 0.20194,
                                0.233947, 0.265953, 0.29796, 0.329967, 0.361973, 0.39398, 0.425987, 
                                0.457993, R_L_max};
const double em_rlunfolding_binning[] = {R_L_min-0.0005,R_L_min, 0.0419067, 0.0739133, 0.10592, 0.137927, 0.169933,
                                         0.20194, 0.233947, 0.265953, 0.29796, 0.329967, 0.361973, 
                                         0.39398, 0.425987, 0.457993, R_L_max, R_L_max + 0.04};

// SL binning
const double sl_p_binning[]   = {4., 10., 16., 22., 28., 40., 52., 76., 100., 150., 250., 500., 1000.};
const double ic_p_binning[]   = {4, 5, 7.5, 10, 12.5, 15, 20, 25, 30, 40, 50, 75, 100, 150, 200, 250, 300, 500, 1000};
const double sl_eta_binning[] = {2.,2.25,2.5,2.75,3.,3.25,3.5,3.75,4.,4.25,4.5};
const int sl_eta_nbins = 10;
const int sl_p_nbins   = 12;
const int ic_p_nbins   = 18;

#endif