#ifndef ANALYSIS_CUTS_H
#define ANALYSIS_CUTS_H

#include "TCut.h"
#include "TString.h"

// Jet cuts
TCut jet_eta_cut = "jet_eta>2.5&&jet_eta<4.";
TCut jet_pt_cut  = Form("jet_pt>%f",jet_pt_min);
TCut jet_cuts    = jet_eta_cut+jet_pt_cut;

// Track cuts
TCut chi2ndf_cut  = "h1_chi2/h1_ndf<3&&h2_chi2/h2_ndf<3";
TCut p_cut        = "h1_p>4&&h1_p<1000&&h2_p>4&&h2_p<1000";
TCut pt_cut       = "h1_pt>0.250&&h2_pt>0.250";
TCut pnnghost_cut = "h1_probnnghost<0.5&&h2_probnnghost<0.5";

TCut trackmc_cuts = p_cut+pt_cut;
TCut track_cuts   = chi2ndf_cut+p_cut+pt_cut+pnnghost_cut;

// Topological cuts
TCut phi_zjet_cut     = "TMath::Abs(z0_phi-jet_phi)>7*TMath::Pi()/8.";
TCut phi_mumjet_cut   = "TMath::Abs(jet_phi-mum_phi)>0.4";
TCut phi_mupjet_cut   = "TMath::Abs(jet_phi-mup_phi)>0.4";
TCut topological_cuts = phi_zjet_cut+phi_mumjet_cut+phi_mupjet_cut;

// Z boson cuts
TCut mu_pt_cut        = "mum_pt>20.&&mup_pt>20.";
TCut mu_eta_cut       = "mum_eta>2&&mum_eta<4.5&&mup_eta>2&&mup_eta<4.5";
TCut mum_mup_mass_cut = "sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))>60 && \
                         sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))<120";
TCut mu_trackprob_cut = "mum_probchi2>0.001&&mup_probchi2>0.001";

TCut Zboson_cuts   = mu_pt_cut+mu_eta_cut+mum_mup_mass_cut+mu_trackprob_cut;
TCut Zbosonmc_cuts = mu_pt_cut+mu_eta_cut+mum_mup_mass_cut;

TCut e2c_nominal_cut = Form("weight*(jet_eta>2.5&&jet_eta<4.&&jet_pt>%f&&h1_chi2/h1_ndf<3&&h2_chi2/h2_ndf<3&&h1_p>4&&h1_p<1000&&h2_p>4&&h2_p<1000&&h1_pt>0.250&&h2_pt>0.250&&\
                             h1_probnnghost<0.5&&h2_probnnghost<0.5&&TMath::Abs(z0_phi-jet_phi)>7*TMath::Pi()/8.&&TMath::Abs(jet_phi-mum_phi)>0.4&&TMath::Abs(jet_phi-mup_phi)>0.4&&\
                             mum_pt>20.&&mup_pt>20.&&mum_eta>2&&mum_eta<4.5&&mup_eta>2&&mup_eta<4.5&&sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))>60&&\
                             sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))<120&&mum_probchi2>0.001&&mup_probchi2>0.001)",jet_pt_min);


// All analysis cuts
TCut diffsign_cut_mc     = jet_cuts+trackmc_cuts+topological_cuts+Zbosonmc_cuts;
TCut diffsign_cut_mcreco = jet_cuts+track_cuts  +topological_cuts+Zboson_cuts;
TCut diffsign_cut_data   = jet_cuts+track_cuts  +topological_cuts+Zboson_cuts;

#endif