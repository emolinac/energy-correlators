#ifndef ANALYSIS_BINNING_H
#define ANALYSIS_BINNING_H

// Limits of variables
const double rl_absmin     = 0.001; // Lowest value achievable for MC. Used as limit for the underflow in the unfolding
const double rl_absmax     = 1.;
const double rl_min        = 0.008;
const double rl_max        = 0.5;
const double rl_logmin     = 0.02;
// const double rl_logmax     = 0.8; // nominal
const double rl_logmax     = 0.5; // nominal
const double jet_pt_min_nom = 20; 
const double jet_pt_max     = 100;
const double jet_e_min      = 100; 
const double jet_e_max      = 4000;
const double eta_min        = 2.;
const double eta_max        = 4.5;
const double weight_max     = 0.1;
const double weight_min     = 0.0001;
const double weight_absmax  = 0.2;
const double weight_absmin  = 1E-5;//5E-6
const double tau_max        = 60;
const double tau_min        = 0.5;

const double rl_min_at     = rl_min;
const double rl_max_at     = TMath::Pi();

// Binning
const int nbin_jet_pt            = 3;
const int nbin_weight            = 10;
const int nbin_weight_unfolding  = nbin_weight+2;
const int nbin_jet_e             = 3;
const int nbin_z_pt              = 3;
const int nbin_jet_pt_corrections = 3;
const int nbin_jet_pt_unfolding  = nbin_jet_pt+2;

const double z_pt_binning[]   = {15,20,30,75};
const double jet_pt_binning[] = {jet_pt_min_nom, 30., 50., jet_pt_max};
const double weight_binning[]            = {0.0001, 0.000396231, 0.000714927, 0.00111206, 0.00163222, 0.0023315, 0.00332664, 0.00485622, 0.00753926, 0.0136004, 0.1};
const double weight_unfoldingbinning[]   = {weight_absmin,0.0001, 0.000396231, 0.000714927, 0.00111206, 0.00163222, 0.0023315, 0.00332664, 0.00485622, 0.00753926, 0.0136004, 0.1,weight_absmax};
const double jet_e_binning[]             = {jet_e_min,350,560,jet_e_max};
const double jet_pt_corrections_binning[] = {20,30,50,100};
const double unfolding_jet_pt_binning[]   = {15,20,30,50,100,150};

// Angular distance binning
const int nbin_rl_nominal           = 15; 
const int nbin_rl_nominal_unfolding = nbin_rl_nominal + 2;

const double rl_chargedeec_binning[]              = {rl_min, 0.0408, 0.0736, 0.1064, 0.1392, 0.172, 0.2048, 0.2376, 0.2704, 0.3032, 0.336, 0.3688, 0.4016, 0.4344, 0.4672, rl_max};
const double unfolding_rl_chargedeec_binning[]    = {rl_absmin,rl_min, 0.0408, 0.0736, 0.1064, 0.1392, 0.172, 0.2048, 0.2376, 0.2704, 0.3032, 0.336, 0.3688, 0.4016, 0.4344, 0.4672, rl_max, rl_absmax};

const double rl_nominal_binning[]           = {rl_logmin, 0.0247871, 0.0307201, 0.0380731, 0.0471861, 0.0584804, 0.072478, 0.089826, 0.111326, 0.137973, 0.170998, 0.211927, 0.262653, 0.32552, 0.403435, rl_logmax};
const double unfolding_rl_nominal_binning[] = {rl_absmin,rl_logmin, 0.0247871, 0.0307201, 0.0380731, 0.0471861, 0.0584804, 0.072478, 0.089826, 0.111326, 0.137973, 0.170998, 0.211927, 0.262653, 0.32552, 0.403435, rl_logmax, rl_absmax};

const int nbin_tau_logbin          = nbin_rl_nominal;
const double tau_nominal_binning[] = {tau_min, 0.68799, 0.94666, 1.30259, 1.79233, 2.46621, 3.39346, 4.66933, 6.4249, 8.84054, 12.1644, 16.738, 23.0311, 31.6904, 43.6053, tau_max};

// Alternative binnings
const int nbin_rl_altlogbin           = 20; 
const int nbin_rl_altlogbin_unfolding = nbin_rl_nominal + 2;

const double rl_altlogbinning[]           = {rl_logmin, 0.0234924, 0.0275946, 0.0324131, 0.0380731, 0.0447214, 0.0525306, 0.0617034, 0.072478, 0.085134, 0.1, 0.117462, 0.137973, 0.162066, 0.190365, 0.223607, 0.262653, 0.308517, 0.36239, 0.42567, rl_logmax};
const double unfolding_rl_altlogbinning[] = {rl_absmin,rl_logmin, 0.0234924, 0.0275946, 0.0324131, 0.0380731, 0.0447214, 0.0525306, 0.0617034, 0.072478, 0.085134, 0.1, 0.117462, 0.137973, 0.162066, 0.190365, 0.223607, 0.262653, 0.308517, 0.36239, 0.42567, rl_logmax, rl_absmax};

const double sl_p_binning[]   = {4., 10., 16., 22., 28., 40., 52., 76., 100., 150., 250., 500., 1000.};
const double ic_p_binning[]   = {4, 5, 7.5, 10, 12.5, 15, 20, 25, 30, 40, 50, 75, 100, 150, 200, 250, 300, 500, 1000};
const double sl_eta_binning[] = {2.,2.25,2.5,2.75,3.,3.25,3.5,3.75,4.,4.25,4.5};
const int sl_eta_nbins = 10;
const int sl_p_nbins   = 12;
const int ic_p_nbins   = 18;

#endif