#ifndef ANALYSIS_CUTS_H_H
#define ANALYSIS_CUTS_H_H

bool apply_jet_cuts(double jet_eta, double jet_pt);

bool apply_muon_cuts(double deltaR_mu_jet, double mu_pt, double mu_eta);

bool apply_zboson_cuts(double deltaphi_zboson_jet, double zboson_mass);

bool apply_chargedtrack_cuts(double charge, double p, double pt, double chi2ndf, double probnnghost, double deltaR_h_jet);

bool apply_chargedtrack_momentum_cuts(double charge, double p, double pt, double deltaR_h_jet);

#endif