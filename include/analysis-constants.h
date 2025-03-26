#ifndef ANALYSIS_CONSTANTS_H
#define ANALYSIS_CONSTANTS_H

#include "TCut.h"
#include "TString.h"
#include "TColor.h"

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

// Correction-related constants
const double corr_rel_error = 0.15;
// const double corr_rel_error = 0.2;
const int ndim_corr = 50;

// Limits of variables
const double R_L_absmin     = 0.001; // Lowest value achievable for MC. Used as limit for the underflow in the unfolding
const double R_L_min        = 0.01;
// const double R_L_max        = 0.5;
// const double R_L_max        = 1.0;
const double R_L_max        = 3;
const double jet_pt_min_nom = 20; 
const double jet_pt_max     = 100;
const double jet_e_min      = 100; 
const double jet_e_max      = 4000;
const double eta_min        = 2.;
const double eta_max        = 4.5;
const double weight_max     = 0.2;
const double weight_min     = 0.00001;

const double R_L_res  = 0.0455;

// Visual constants
const double std_marker_size  = 1.0;
const int    std_marker_style = 8;
const int    std_line_width   = 3;

const int std_marker_style_jet_pt[]  = {24,25,26,32,27};
const int corr_marker_style_jet_pt[] = {20,21,22,23,33};

const int std_marker_color_jet_pt[]  = {868,797,618,633,820,418,810,616,600,1}; // darker
const int corr_marker_color_jet_pt[] = {868,797,618,633,820,418,810,616,600,1}; // lighter

// Binning
const int Nbin_R_L    = 10;
const int Nbin_jet_pt = 3;
const int Nbin_weight = 10;
const int Nbin_jet_e  = 3;
const int Nbin_z_pt   = 3;

const double z_pt_binning[]   = {10,15,30,150};
const double jet_pt_binning[] = {jet_pt_min_nom, 30., 50., jet_pt_max};
const double weight_binning[] = {1e-05, 0.000370882, 0.000689666, 0.00108385, 0.00159442, 0.00228839, 0.00327794, 0.00482046, 0.00751032, 0.0135402, 0.2};
const double jet_e_binning[]  = {jet_e_min,350,560,jet_e_max};
const double rl_binning[]     = {R_L_min, 0.0176894, 0.0312913, 0.0553524, 0.0979148, 0.173205, 0.306389, 0.541982, 0.958732, 1.69593, R_L_max};

const double unfolding_jetpt_binning[] = {15,20,30,50,100,150};
const double unfolding_rl_binning[]    = {R_L_absmin,R_L_min, 0.0176894, 0.0312913, 0.0553524, 0.0979148, 0.173205, 0.306389, 0.541982, 0.958732, 1.69593, R_L_max, R_L_max + 0.04};

const int em_jetptcorrection_nbins = 7;
const double corrections_jetpt_binning[] = {20,22.5,25,30,40,50,75,100};

// Alternative binnings
const double sl_p_binning[]   = {4., 10., 16., 22., 28., 40., 52., 76., 100., 150., 250., 500., 1000.};
const double ic_p_binning[]   = {4, 5, 7.5, 10, 12.5, 15, 20, 25, 30, 40, 50, 75, 100, 150, 200, 250, 300, 500, 1000};
const double sl_eta_binning[] = {2.,2.25,2.5,2.75,3.,3.25,3.5,3.75,4.,4.25,4.5};
const int sl_eta_nbins = 10;
const int sl_p_nbins   = 12;
const int ic_p_nbins   = 18;

#endif