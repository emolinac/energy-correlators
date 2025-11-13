#ifndef ANALYSIS_CUTS_H
#define ANALYSIS_CUTS_H

#include "TCut.h"
#include "TMath.h"
#include "TString.h"
#include "analysis-constants.h"
#include "analysis-binning.h"

extern TCut pair_matching_cut;
extern TCut eec_cut;
extern TCut pair_cut;
extern TCut pair_jet_pt_cut[];
extern TCut pair_signal_cut;
extern TCut single_signal_cut;

bool apply_jet_cuts(double jet_eta, double jet_pt);

bool apply_jet_id_cuts(double mpt, double nPVtrk, double cpf, double mpf);

bool apply_muon_cuts(double deltaR_mu_jet, double mu_pt, double mu_eta);

bool apply_zboson_cuts(double deltaphi_zboson_jet, double zboson_mass);

bool apply_chargedtrack_cuts(double charge, double p, double pt, double chi2ndf, double probnnghost, double eta);

bool apply_chargedtrack_cuts(double charge, double p, double pt, double chi2ndf, double probnnghost, double eta, double deltaR_h_jet);

bool apply_chargedtrack_momentum_cuts(double charge, double p, double pt, double eta);

bool apply_chargedtrack_momentum_cuts(double charge, double p, double pt, double eta, double deltaR_h_jet);

#endif