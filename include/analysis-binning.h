#ifndef ANALYSIS_BINNING_H
#define ANALYSIS_BINNING_H

// Limits of variables
const double R_L_absmin     = 0.001; // Lowest value achievable for MC. Used as limit for the underflow in the unfolding
const double R_L_absmax     = 1.;
const double R_L_min        = 0.008;
const double R_L_max        = 0.8;
const double jet_pt_min_nom = 20; 
const double jet_pt_max     = 100;
const double jet_e_min      = 100; 
const double jet_e_max      = 4000;
const double eta_min        = 2.;
const double eta_max        = 4.5;
const double weight_max     = 0.2;
const double weight_min     = 0.00001;

const double R_L_min_at     = R_L_min;
const double R_L_max_at     = TMath::Pi();

const double R_L_res  = 0.015;

// Binning
// const int Nbin_R_L    = 10; // Golden binning
const int Nbin_R_L    = 20;
const int Nbin_jet_pt = 3;
const int Nbin_weight = 10;
const int Nbin_jet_e  = 3;
const int Nbin_z_pt   = 3;

const double z_pt_binning[]   = {15,20,30,75};
const double jet_pt_binning[] = {jet_pt_min_nom, 30., 50., jet_pt_max};
const double weight_binning[] = {1e-05, 0.000370882, 0.000689666, 0.00108385, 0.00159442, 0.00228839, 0.00327794, 0.00482046, 0.00751032, 0.0135402, 0.2};
const double jet_e_binning[]  = {jet_e_min,350,560,jet_e_max};
// const double rl_binning[]     = {R_L_min, 0.0903675, 0.116661, 0.150606, 0.194427, 0.250998, 0.32403, 0.418311, 0.540025, 0.697153, R_L_max}; // Golden binning
const double rl_binning[]     = {R_L_min, 0.044, 0.08, 0.116, 0.152, 0.188, 0.224, 0.26, 0.296, 0.332, 0.368, 0.404, 0.44, 0.476, 0.512, 0.548, 0.584, 0.62, 0.656, 0.692, 0.728, 0.764, R_L_max};
const double rl_binning_at[]  = {R_L_min_at, 0.102401, 0.149799, 0.219136, 0.320567, 0.468947, 0.686008, 1.00354, 1.46805, 2.14756, R_L_max_at};

const int Nbin_jetpt_corrections = 7;
const double corrections_jetpt_binning[] = {20,22.5,25,30,40,50,75,100};

const int Nbin_R_L_unfolding    = Nbin_R_L + 2;
const int Nbin_jet_pt_unfolding = Nbin_jet_pt+2;
const double unfolding_jetpt_binning[] = {15,20,30,50,100,150};
// const double unfolding_jetpt_binning[] = {15,20,22.5,25,30,40,50,75,100,150};
const double unfolding_rl_binning[]    = {R_L_absmin,R_L_min, 0.044, 0.08, 0.116, 0.152, 0.188, 0.224, 0.26, 0.296, 0.332, 0.368, 0.404, 0.44, 0.476, 0.512, 0.548, 0.584, 0.62, 0.656, 0.692, 0.728, 0.764, R_L_max, R_L_absmax};

// Alternative binnings
const double sl_p_binning[]   = {4., 10., 16., 22., 28., 40., 52., 76., 100., 150., 250., 500., 1000.};
const double ic_p_binning[]   = {4, 5, 7.5, 10, 12.5, 15, 20, 25, 30, 40, 50, 75, 100, 150, 200, 250, 300, 500, 1000};
const double sl_eta_binning[] = {2.,2.25,2.5,2.75,3.,3.25,3.5,3.75,4.,4.25,4.5};
const int sl_eta_nbins = 10;
const int sl_p_nbins   = 12;
const int ic_p_nbins   = 18;

#endif