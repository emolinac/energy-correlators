#ifndef ANALYSIS_BINNING_H
#define ANALYSIS_BINNING_H

// #include "analysis-cuts.h"

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
const double weight_absmax  = 0.35;
const double weight_absmin  = 2E-6;

// // SQRT WEIGHT BINNING
// const double weight_absmax  = 0.6;
// const double weight_absmin  = 1E-3;
// // SQRT WEIGHT BINNING

const double ptprod_absmax  = 2000.;
const double ptprod_absmin  = 0.05;
const double ptprod_resolution = 1.;

const double tau_max        = 60;
const double tau_min        = 0.5;

const double rl_min_at     = rl_min;
const double rl_max_at     = TMath::Pi();

// Binning
const int nbin_jet_pt = 3;
const int nbin_weight = 30;
const int nbin_jet_e  = 3;
const int nbin_z_pt   = 3;
const int nbin_h_pt   = 20;

const int nbin_jet_pt_corrections = 3;
const int nbin_jet_pt_unfolding  = nbin_jet_pt+2;

const double z_pt_binning[]   = {15,20,30,75};

// const double weight_binning[] = {weight_absmin, 3.65748e-06, 6.68858e-06, 1.22317e-05, 2.23685e-05, 4.09062e-05, 7.48069e-05, 0.000136802, 0.000250176, 0.000457506, 0.00083666, 0.00153003, 0.00279803, 0.00511687, 0.00935743, 0.0171123, 0.031294, 0.0572285, 0.104656, 0.191389, weight_absmax};

// const double weight_binning[] = {weight_absmin, 3.24154e-06, 5.25379e-06, 8.51519e-06, 1.38012e-05, 2.23685e-05, 3.62543e-05, 5.87598e-05, 9.52362e-05, 0.000154356, 0.000250176, 0.000405477, 0.000657186, 0.00106515, 0.00172636, 0.00279803, 0.00453497, 0.00735014, 0.0119129, 0.0193081, 0.031294, 0.0507203, 0.082206, 0.133237, 0.215947, weight_absmax};

const double weight_binning[] = {weight_absmin, 2.99087e-06, 4.47266e-06, 6.68858e-06, 1.00023e-05, 1.49579e-05, 2.23685e-05, 3.34507e-05, 5.00234e-05, 7.48069e-05, 0.000111869, 0.000167293, 0.000250176, 0.000374122, 0.000559476, 0.00083666, 0.00125117, 0.00187105, 0.00279803, 0.00418428, 0.00625732, 0.00935743, 0.0139934, 0.0209263, 0.031294, 0.0467981, 0.0699836, 0.104656, 0.156506, 0.234045, weight_absmax};

// // SQRT WEIGHT BINNING
// const double weight_binning[] = {weight_absmin, 0.00137692, 0.0018959, 0.00261049, 0.00359443, 0.00494923, 0.00681468, 0.00938324, 0.0129199, 0.0177897, 0.0244949, 0.0337274, 0.0464398, 0.0639438, 0.0880452, 0.121231, 0.166925, 0.229842, 0.316473, 0.435756, weight_absmax};
// // SQRT WEIGHT BINNING


// // LOG BINNING PTPROD
// const int nbin_ptprod = 20;
// const double ptprod_binning[] = {ptprod_absmin, 0.0849323, 0.14427, 0.245064, 0.416277, 0.707107, 1.20112, 2.04029, 3.46572, 5.88704, 10, 16.9865, 28.854, 49.0127, 83.2553, 141.421, 240.225, 408.057, 693.145, 1177.41, ptprod_absmax};
// // LOG BINNING PTPROD

const int nbin_ptprod = 5;
const double ptprod_binning[] = {ptprod_absmin, 1.01998, 2.27994, 4.81988, 13.7597, ptprod_absmax};

const double h_pt_binning[] = {0.25, 0.32583, 0.424662, 0.553471, 0.72135, 0.940151, 1.22532, 1.59698, 2.08138, 2.71271, 3.53553, 4.60794, 6.00562, 7.82726, 10.2014, 13.2957, 17.3286, 22.5848, 29.4352, 38.3635, 50};

const double jet_pt_binning[]           = {20,30,50,100};
// const double unfolding_jet_pt_binning[] = {15,20,30,50,100,150};
const double unfolding_jet_pt_binning[] = {12.5,15,20,30,50,100};

// Angular distance binning
const int nbin_rl_nominal           = 15; 
const int nbin_chargedeec_nominal   = 15;
const int nbin_rl_nominal_unfolding = nbin_rl_nominal + 2;

const double rl_chargedeec_binning[]              = {rl_logmin, 0.052, 0.084, 0.116, 0.148, 0.18, 0.212, 0.244, 0.276, 0.308, 0.34, 0.372, 0.404, 0.436, 0.468, rl_logmax};
const double unfolding_rl_chargedeec_binning[]    = {rl_absmin,rl_logmin, 0.052, 0.084, 0.116, 0.148, 0.18, 0.212, 0.244, 0.276, 0.308, 0.34, 0.372, 0.404, 0.436, 0.468, rl_logmax, rl_absmax};

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