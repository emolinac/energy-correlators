#ifndef ANALYSIS_CUTS_H
#define ANALYSIS_CUTS_H

#include "TCut.h"
#include "TMath.h"
#include "TString.h"
#include "analysis-constants.h"
#include "analysis-binning.h"

// Misc analysis cuts
extern TCut pair_matching_cut;

// Nominal Analysis Cuts
extern TCut eec_cut;
extern TCut pair_cut;

// Cuts as function of jet pt
extern TCut eec_jet_pt_cut[];
extern TCut eec_eqcharge_jet_pt_cut[];
extern TCut eec_neqcharge_jet_pt_cut[];
extern TCut eec_zpt_cut_weightpt[];
extern TCut pair_jet_pt_cut[];
extern TCut pair_zpt_cut[];

// _______________________________ Purity Analysis Cuts _______________________________ //
extern TCut pair_signal_cut;
extern TCut single_signal_cut;

extern TCut purity_corr_singletrack;
extern TCut efficiency_corr_singletrack;
extern TCut full_corr_singletrack;
extern TCut eec_purity_corr_singletrack;
extern TCut eec_efficiency_corr_singletrack;
extern TCut eec_full_corr_singletrack;

extern TCut pair_purity_corr_singletrack_weightpt;

extern TString muons_eff;

extern TCut jet_full_corr[];
extern TCut pair_jet_pt_signal_cut[];

// ANALYSIS VARIATIONS CUTS
extern TCut eec_jet_pt_cut_pp[];
extern TCut eec_jet_pt_cut_mm[];
extern TCut eec_jet_pt_cut_pm[];
extern TCut eec_jet_pt_cut_kaon[];
extern TCut eec_jet_pt_cut_nokaon[];


bool apply_jet_cuts(double jet_eta, double jet_pt);

bool apply_jet_id_cuts(double mpt, double nPVtrk, double cpf, double mpf);

bool apply_muon_cuts(double deltaR_mu_jet, double mu_pt, double mu_eta);

bool apply_zboson_cuts(double deltaphi_zboson_jet, double zboson_mass);

bool apply_chargedtrack_cuts(double charge, double p, double pt, double chi2ndf, double probnnghost, double eta);

bool apply_chargedtrack_cuts(double charge, double p, double pt, double chi2ndf, double probnnghost, double eta, double deltaR_h_jet);

bool apply_chargedtrack_momentum_cuts(double charge, double p, double pt, double eta);

bool apply_chargedtrack_momentum_cuts(double charge, double p, double pt, double eta, double deltaR_h_jet);

#endif