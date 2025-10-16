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
// const double corr_rel_error = 0.05;
// const double corr_rel_error = 0.2; // nominal
const double corr_rel_error = 0.5;
const int ndim_corr = 50;

const double rl_resolution = 0.015;

const int nominal_niter = 4;

const int smearing_constant = 30;

// Weight mapping
const double slope = 1.0382;
const double intercept = 0.000147376;

const double slope_jetpt[]     = {1.08357, 1.08947, 1.1131};
const double intercept_jetpt[] = {0.000161732, -4.37428e-06, 5.46739e-06};

// Visual constants
const double std_marker_size  = 1.0;
const int    std_marker_style = 8;
const int    std_line_width   = 3;

const int std_marker_style_jet_pt[]  = {24,25,26,32,27};
const int corr_marker_style_jet_pt[] = {20,21,22,23,33};

const int std_marker_color_jet_pt[]  = {868,797,618,633,820,418,810,616,600,1};
const int corr_marker_color_jet_pt[] = {868,797,618,633,820,418,810,616,600,1};

const double text_size_correction_plots = 0.015;

#endif