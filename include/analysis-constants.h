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

// Define binning
const int Nbin_X_L    = 25;
const int Nbin_jet_pt = 3;

const double X_L_min    = 0.005;
const double X_L_max    = 4.;
const double jet_pt_min = 20; 
const double jet_pt_max = 100;
const double eta_min    = 2.5;
const double eta_max    = 4.;

// Binning
const double jet_pt_limits[] = {jet_pt_min, 25.644, 35.532, jet_pt_max};

// Visual constants
const double std_marker_size  = 1.3;
const int    std_marker_style = 8;
const int    std_line_width   = 3;

#endif