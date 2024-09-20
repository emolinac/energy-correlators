#ifndef ANALYSIS_CUTS_H
#define ANALYSIS_CUTS_H

#include "TCut.h"
#include "TString.h"
#include "analysis-constants.h"

// _______________________________ Z Tagged Jet Cuts _______________________________ //

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

// _______________________________ Nominal Analysis Cuts _______________________________ //

// MC cut
TCut e2c_mc_cut = Form("weight*(jet_eta>2.5&&jet_eta<4.&&jet_pt>%f&&jet_pt<%f&&h1_p>4&&h1_p<1000&&h2_p>4&&h2_p<1000&&h1_pt>0.250&&h2_pt>0.250&&\
                  TMath::Abs(z0_phi-jet_phi)>7*TMath::Pi()/8.&&TMath::Abs(jet_phi-mum_phi)>0.4&&TMath::Abs(jet_phi-mup_phi)>0.4&&\
                  mum_pt>20.&&mup_pt>20.&&mum_eta>2&&mum_eta<4.5&&mup_eta>2&&mup_eta<4.5&&\
                  sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))>60&&\
                  sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))<120)",jet_pt_min,jet_pt_max);

TCut pair_mc_cut = Form("jet_eta>2.5&&jet_eta<4.&&jet_pt>%f&&jet_pt<%f&&h1_p>4&&h1_p<1000&&h2_p>4&&h2_p<1000&&h1_pt>0.250&&h2_pt>0.250&&\
              TMath::Abs(z0_phi-jet_phi)>7*TMath::Pi()/8.&&TMath::Abs(jet_phi-mum_phi)>0.4&&TMath::Abs(jet_phi-mup_phi)>0.4&&\
              mum_pt>20.&&mup_pt>20.&&mum_eta>2&&mum_eta<4.5&&mup_eta>2&&mup_eta<4.5&&\
              sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))>60&&\
              sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))<120",jet_pt_min,jet_pt_max);

// Data and MCReco cut
TCut e2c_data_cut = Form("weight*(jet_eta>2.5&&jet_eta<4.&&jet_pt>%f&&jet_pt<%f&&h1_chi2/h1_ndf<3&&h2_chi2/h2_ndf<3&&h1_p>4&&h1_p<1000&&h2_p>4&&h2_p<1000&&h1_pt>0.250&&h2_pt>0.250&&\
                    h1_probnnghost<0.5&&h2_probnnghost<0.5&&TMath::Abs(z0_phi-jet_phi)>7*TMath::Pi()/8.&&TMath::Abs(jet_phi-mum_phi)>0.4&&TMath::Abs(jet_phi-mup_phi)>0.4&&\
                    mum_pt>20.&&mup_pt>20.&&mum_eta>2&&mum_eta<4.5&&mup_eta>2&&mup_eta<4.5&&sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))>60&&\
                    sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))<120&&mum_probchi2>0.001&&mup_probchi2>0.001)",jet_pt_min,jet_pt_max);

TCut pair_data_cut = Form("jet_eta>2.5&&jet_eta<4.&&jet_pt>%f&&jet_pt<%f&&h1_chi2/h1_ndf<3&&h2_chi2/h2_ndf<3&&h1_p>4&&h1_p<1000&&h2_p>4&&h2_p<1000&&h1_pt>0.250&&h2_pt>0.250&&\
                           h1_probnnghost<0.5&&h2_probnnghost<0.5&&TMath::Abs(z0_phi-jet_phi)>7*TMath::Pi()/8.&&TMath::Abs(jet_phi-mum_phi)>0.4&&TMath::Abs(jet_phi-mup_phi)>0.4&&\
                           mum_pt>20.&&mup_pt>20.&&mum_eta>2&&mum_eta<4.5&&mup_eta>2&&mup_eta<4.5&&sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))>60&&\
                           sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))<120&&mum_probchi2>0.001&&mup_probchi2>0.001",jet_pt_min,jet_pt_max);

// Cuts as function of jet pt
TCut e2c_mc_jetpt_cut[] = {Form("weight*(jet_eta>2.5&&jet_eta<4.&&jet_pt>%f&&jet_pt<%f&&h1_p>4&&h1_p<1000&&h2_p>4&&h2_p<1000&&h1_pt>0.250&&h2_pt>0.250&&\
                           TMath::Abs(z0_phi-jet_phi)>7*TMath::Pi()/8.&&TMath::Abs(jet_phi-mum_phi)>0.4&&TMath::Abs(jet_phi-mup_phi)>0.4&&\
                           mum_pt>20.&&mup_pt>20.&&mum_eta>2&&mum_eta<4.5&&mup_eta>2&&mup_eta<4.5&&\
                           sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))>60&&\
                           sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))<120)",jet_pt_binning[0],jet_pt_binning[1]),
                           Form("weight*(jet_eta>2.5&&jet_eta<4.&&jet_pt>%f&&jet_pt<%f&&h1_p>4&&h1_p<1000&&h2_p>4&&h2_p<1000&&h1_pt>0.250&&h2_pt>0.250&&\
                           TMath::Abs(z0_phi-jet_phi)>7*TMath::Pi()/8.&&TMath::Abs(jet_phi-mum_phi)>0.4&&TMath::Abs(jet_phi-mup_phi)>0.4&&\
                           mum_pt>20.&&mup_pt>20.&&mum_eta>2&&mum_eta<4.5&&mup_eta>2&&mup_eta<4.5&&\
                           sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))>60&&\
                           sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))<120)",jet_pt_binning[1],jet_pt_binning[2]),
                           Form("weight*(jet_eta>2.5&&jet_eta<4.&&jet_pt>%f&&jet_pt<%f&&h1_p>4&&h1_p<1000&&h2_p>4&&h2_p<1000&&h1_pt>0.250&&h2_pt>0.250&&\
                           TMath::Abs(z0_phi-jet_phi)>7*TMath::Pi()/8.&&TMath::Abs(jet_phi-mum_phi)>0.4&&TMath::Abs(jet_phi-mup_phi)>0.4&&\
                           mum_pt>20.&&mup_pt>20.&&mum_eta>2&&mum_eta<4.5&&mup_eta>2&&mup_eta<4.5&&\
                           sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))>60&&\
                           sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))<120)",jet_pt_binning[2],jet_pt_binning[3]),
                           Form("weight*(jet_eta>2.5&&jet_eta<4.&&jet_pt>%f&&jet_pt<%f&&h1_p>4&&h1_p<1000&&h2_p>4&&h2_p<1000&&h1_pt>0.250&&h2_pt>0.250&&\
                           TMath::Abs(z0_phi-jet_phi)>7*TMath::Pi()/8.&&TMath::Abs(jet_phi-mum_phi)>0.4&&TMath::Abs(jet_phi-mup_phi)>0.4&&\
                           mum_pt>20.&&mup_pt>20.&&mum_eta>2&&mum_eta<4.5&&mup_eta>2&&mup_eta<4.5&&\
                           sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))>60&&\
                           sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))<120)",jet_pt_binning[3],jet_pt_binning[4]),
                           Form("weight*(jet_eta>2.5&&jet_eta<4.&&jet_pt>%f&&jet_pt<%f&&h1_p>4&&h1_p<1000&&h2_p>4&&h2_p<1000&&h1_pt>0.250&&h2_pt>0.250&&\
                           TMath::Abs(z0_phi-jet_phi)>7*TMath::Pi()/8.&&TMath::Abs(jet_phi-mum_phi)>0.4&&TMath::Abs(jet_phi-mup_phi)>0.4&&\
                           mum_pt>20.&&mup_pt>20.&&mum_eta>2&&mum_eta<4.5&&mup_eta>2&&mup_eta<4.5&&\
                           sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))>60&&\
                           sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))<120)",jet_pt_binning[4],jet_pt_binning[5])};


TCut e2c_jetpt_cut[] = {Form("weight*(jet_eta>2.5&&jet_eta<4.&&jet_pt>%f&&jet_pt<%f&&h1_chi2/h1_ndf<3&&h2_chi2/h2_ndf<3&&h1_p>4&&h1_p<1000&&h2_p>4&&h2_p<1000&&h1_pt>0.250&&h2_pt>0.250&&\
                        h1_probnnghost<0.5&&h2_probnnghost<0.5&&TMath::Abs(z0_phi-jet_phi)>7*TMath::Pi()/8.&&TMath::Abs(jet_phi-mum_phi)>0.4&&TMath::Abs(jet_phi-mup_phi)>0.4&&\
                        mum_pt>20.&&mup_pt>20.&&mum_eta>2&&mum_eta<4.5&&mup_eta>2&&mup_eta<4.5&&sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))>60&&\
                        sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))<120&&mum_probchi2>0.001&&mup_probchi2>0.001)",jet_pt_binning[0],jet_pt_binning[1]),
                        Form("weight*(jet_eta>2.5&&jet_eta<4.&&jet_pt>%f&&jet_pt<%f&&h1_chi2/h1_ndf<3&&h2_chi2/h2_ndf<3&&h1_p>4&&h1_p<1000&&h2_p>4&&h2_p<1000&&h1_pt>0.250&&h2_pt>0.250&&\
                        h1_probnnghost<0.5&&h2_probnnghost<0.5&&TMath::Abs(z0_phi-jet_phi)>7*TMath::Pi()/8.&&TMath::Abs(jet_phi-mum_phi)>0.4&&TMath::Abs(jet_phi-mup_phi)>0.4&&\
                        mum_pt>20.&&mup_pt>20.&&mum_eta>2&&mum_eta<4.5&&mup_eta>2&&mup_eta<4.5&&sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))>60&&\
                        sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))<120&&mum_probchi2>0.001&&mup_probchi2>0.001)",jet_pt_binning[1],jet_pt_binning[2]),
                        Form("weight*(jet_eta>2.5&&jet_eta<4.&&jet_pt>%f&&jet_pt<%f&&h1_chi2/h1_ndf<3&&h2_chi2/h2_ndf<3&&h1_p>4&&h1_p<1000&&h2_p>4&&h2_p<1000&&h1_pt>0.250&&h2_pt>0.250&&\
                        h1_probnnghost<0.5&&h2_probnnghost<0.5&&TMath::Abs(z0_phi-jet_phi)>7*TMath::Pi()/8.&&TMath::Abs(jet_phi-mum_phi)>0.4&&TMath::Abs(jet_phi-mup_phi)>0.4&&\
                        mum_pt>20.&&mup_pt>20.&&mum_eta>2&&mum_eta<4.5&&mup_eta>2&&mup_eta<4.5&&sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))>60&&\
                        sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))<120&&mum_probchi2>0.001&&mup_probchi2>0.001)",jet_pt_binning[2],jet_pt_binning[3]),
                        Form("weight*(jet_eta>2.5&&jet_eta<4.&&jet_pt>%f&&jet_pt<%f&&h1_chi2/h1_ndf<3&&h2_chi2/h2_ndf<3&&h1_p>4&&h1_p<1000&&h2_p>4&&h2_p<1000&&h1_pt>0.250&&h2_pt>0.250&&\
                        h1_probnnghost<0.5&&h2_probnnghost<0.5&&TMath::Abs(z0_phi-jet_phi)>7*TMath::Pi()/8.&&TMath::Abs(jet_phi-mum_phi)>0.4&&TMath::Abs(jet_phi-mup_phi)>0.4&&\
                        mum_pt>20.&&mup_pt>20.&&mum_eta>2&&mum_eta<4.5&&mup_eta>2&&mup_eta<4.5&&sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))>60&&\
                        sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))<120&&mum_probchi2>0.001&&mup_probchi2>0.001)",jet_pt_binning[3],jet_pt_binning[4]),
                        Form("weight*(jet_eta>2.5&&jet_eta<4.&&jet_pt>%f&&jet_pt<%f&&h1_chi2/h1_ndf<3&&h2_chi2/h2_ndf<3&&h1_p>4&&h1_p<1000&&h2_p>4&&h2_p<1000&&h1_pt>0.250&&h2_pt>0.250&&\
                        h1_probnnghost<0.5&&h2_probnnghost<0.5&&TMath::Abs(z0_phi-jet_phi)>7*TMath::Pi()/8.&&TMath::Abs(jet_phi-mum_phi)>0.4&&TMath::Abs(jet_phi-mup_phi)>0.4&&\
                        mum_pt>20.&&mup_pt>20.&&mum_eta>2&&mum_eta<4.5&&mup_eta>2&&mup_eta<4.5&&sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))>60&&\
                        sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))<120&&mum_probchi2>0.001&&mup_probchi2>0.001)",jet_pt_binning[4],jet_pt_binning[5])};

TCut pair_data_jetpt_cut[] = {Form("jet_eta>2.5&&jet_eta<4.&&jet_pt>%f&&jet_pt<%f&&h1_chi2/h1_ndf<3&&h2_chi2/h2_ndf<3&&h1_p>4&&h1_p<1000&&h2_p>4&&h2_p<1000&&h1_pt>0.250&&h2_pt>0.250&&\
                         h1_probnnghost<0.5&&h2_probnnghost<0.5&&TMath::Abs(z0_phi-jet_phi)>7*TMath::Pi()/8.&&TMath::Abs(jet_phi-mum_phi)>0.4&&TMath::Abs(jet_phi-mup_phi)>0.4&&\
                         mum_pt>20.&&mup_pt>20.&&mum_eta>2&&mum_eta<4.5&&mup_eta>2&&mup_eta<4.5&&sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))>60&&\
                         sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))<120&&mum_probchi2>0.001&&mup_probchi2>0.001",jet_pt_binning[0],jet_pt_binning[1]),
                         Form("jet_eta>2.5&&jet_eta<4.&&jet_pt>%f&&jet_pt<%f&&h1_chi2/h1_ndf<3&&h2_chi2/h2_ndf<3&&h1_p>4&&h1_p<1000&&h2_p>4&&h2_p<1000&&h1_pt>0.250&&h2_pt>0.250&&\
                         h1_probnnghost<0.5&&h2_probnnghost<0.5&&TMath::Abs(z0_phi-jet_phi)>7*TMath::Pi()/8.&&TMath::Abs(jet_phi-mum_phi)>0.4&&TMath::Abs(jet_phi-mup_phi)>0.4&&\
                         mum_pt>20.&&mup_pt>20.&&mum_eta>2&&mum_eta<4.5&&mup_eta>2&&mup_eta<4.5&&sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))>60&&\
                         sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))<120&&mum_probchi2>0.001&&mup_probchi2>0.001",jet_pt_binning[1],jet_pt_binning[2]),
                         Form("jet_eta>2.5&&jet_eta<4.&&jet_pt>%f&&jet_pt<%f&&h1_chi2/h1_ndf<3&&h2_chi2/h2_ndf<3&&h1_p>4&&h1_p<1000&&h2_p>4&&h2_p<1000&&h1_pt>0.250&&h2_pt>0.250&&\
                         h1_probnnghost<0.5&&h2_probnnghost<0.5&&TMath::Abs(z0_phi-jet_phi)>7*TMath::Pi()/8.&&TMath::Abs(jet_phi-mum_phi)>0.4&&TMath::Abs(jet_phi-mup_phi)>0.4&&\
                         mum_pt>20.&&mup_pt>20.&&mum_eta>2&&mum_eta<4.5&&mup_eta>2&&mup_eta<4.5&&sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))>60&&\
                         sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))<120&&mum_probchi2>0.001&&mup_probchi2>0.001",jet_pt_binning[2],jet_pt_binning[3]),
                         Form("jet_eta>2.5&&jet_eta<4.&&jet_pt>%f&&jet_pt<%f&&h1_chi2/h1_ndf<3&&h2_chi2/h2_ndf<3&&h1_p>4&&h1_p<1000&&h2_p>4&&h2_p<1000&&h1_pt>0.250&&h2_pt>0.250&&\
                         h1_probnnghost<0.5&&h2_probnnghost<0.5&&TMath::Abs(z0_phi-jet_phi)>7*TMath::Pi()/8.&&TMath::Abs(jet_phi-mum_phi)>0.4&&TMath::Abs(jet_phi-mup_phi)>0.4&&\
                         mum_pt>20.&&mup_pt>20.&&mum_eta>2&&mum_eta<4.5&&mup_eta>2&&mup_eta<4.5&&sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))>60&&\
                         sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))<120&&mum_probchi2>0.001&&mup_probchi2>0.001",jet_pt_binning[3],jet_pt_binning[4]),
                         Form("jet_eta>2.5&&jet_eta<4.&&jet_pt>%f&&jet_pt<%f&&h1_chi2/h1_ndf<3&&h2_chi2/h2_ndf<3&&h1_p>4&&h1_p<1000&&h2_p>4&&h2_p<1000&&h1_pt>0.250&&h2_pt>0.250&&\
                         h1_probnnghost<0.5&&h2_probnnghost<0.5&&TMath::Abs(z0_phi-jet_phi)>7*TMath::Pi()/8.&&TMath::Abs(jet_phi-mum_phi)>0.4&&TMath::Abs(jet_phi-mup_phi)>0.4&&\
                         mum_pt>20.&&mup_pt>20.&&mum_eta>2&&mum_eta<4.5&&mup_eta>2&&mup_eta<4.5&&sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))>60&&\
                         sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))<120&&mum_probchi2>0.001&&mup_probchi2>0.001",jet_pt_binning[4],jet_pt_binning[5])};

// _______________________________ Purity Analysis Cuts _______________________________ //
TCut e2c_signal_cut = Form("weight*(jet_eta>2.5&&jet_eta<4.&&jet_pt>%f&&jet_pt<%f&&h1_chi2/h1_ndf<3&&h2_chi2/h2_ndf<3&&h1_p>4&&h1_p<1000&&h2_p>4&&h2_p<1000&&h1_pt>0.250&&h2_pt>0.250&&\
                      h1_probnnghost<0.5&&h2_probnnghost<0.5&&TMath::Abs(z0_phi-jet_phi)>7*TMath::Pi()/8.&&TMath::Abs(jet_phi-mum_phi)>0.4&&TMath::Abs(jet_phi-mup_phi)>0.4&&\
                      mum_pt>20.&&mup_pt>20.&&mum_eta>2&&mum_eta<4.5&&mup_eta>2&&mup_eta<4.5&&sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))>60&&\
                      sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))<120&&mum_probchi2>0.001&&mup_probchi2>0.001)",jet_pt_min,jet_pt_max);

TCut pair_signal_cut = Form("(R_h1>%f&&R_h1<%f)&&(R_h2>%f&&R_h2<%f)&&jet_eta>2.5&&jet_eta<4.&&jet_pt>%f&&jet_pt<%f&&h1_chi2/h1_ndf<3&&h2_chi2/h2_ndf<3&&h1_p>4&&h1_p<1000&&h2_p>4&&h2_p<1000&&h1_pt>0.250&&h2_pt>0.250&&\
                             h1_probnnghost<0.5&&h2_probnnghost<0.5&&TMath::Abs(z0_phi-jet_phi)>7*TMath::Pi()/8.&&TMath::Abs(jet_phi-mum_phi)>0.4&&TMath::Abs(jet_phi-mup_phi)>0.4&&\
                             mum_pt>20.&&mup_pt>20.&&mum_eta>2&&mum_eta<4.5&&mup_eta>2&&mup_eta<4.5&&sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))>60&&\
                             sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))<120&&mum_probchi2>0.001&&mup_probchi2>0.001",
                             R_match_min,R_match_max,R_match_min,R_match_max,jet_pt_min,jet_pt_max);

TCut e2c_jetpt_signal_cut[] = {Form("weight*(jet_eta>2.5&&jet_eta<4.&&jet_pt>%f&&jet_pt<%f&&h1_chi2/h1_ndf<3&&h2_chi2/h2_ndf<3&&h1_p>4&&h1_p<1000&&h2_p>4&&h2_p<1000&&h1_pt>0.250&&h2_pt>0.250&&\
                               h1_probnnghost<0.5&&h2_probnnghost<0.5&&TMath::Abs(z0_phi-jet_phi)>7*TMath::Pi()/8.&&TMath::Abs(jet_phi-mum_phi)>0.4&&TMath::Abs(jet_phi-mup_phi)>0.4&&\
                               mum_pt>20.&&mup_pt>20.&&mum_eta>2&&mum_eta<4.5&&mup_eta>2&&mup_eta<4.5&&sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))>60&&\
                               sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))<120&&mum_probchi2>0.001&&mup_probchi2>0.001)"
                               ,jet_pt_binning[0],jet_pt_binning[1]),
                               Form("weight*(jet_eta>2.5&&jet_eta<4.&&jet_pt>%f&&jet_pt<%f&&h1_chi2/h1_ndf<3&&h2_chi2/h2_ndf<3&&h1_p>4&&h1_p<1000&&h2_p>4&&h2_p<1000&&h1_pt>0.250&&h2_pt>0.250&&\
                               h1_probnnghost<0.5&&h2_probnnghost<0.5&&TMath::Abs(z0_phi-jet_phi)>7*TMath::Pi()/8.&&TMath::Abs(jet_phi-mum_phi)>0.4&&TMath::Abs(jet_phi-mup_phi)>0.4&&\
                               mum_pt>20.&&mup_pt>20.&&mum_eta>2&&mum_eta<4.5&&mup_eta>2&&mup_eta<4.5&&sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))>60&&\
                               sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))<120&&mum_probchi2>0.001&&mup_probchi2>0.001)"
                               ,jet_pt_binning[1],jet_pt_binning[2]),
                               Form("weight*(jet_eta>2.5&&jet_eta<4.&&jet_pt>%f&&jet_pt<%f&&h1_chi2/h1_ndf<3&&h2_chi2/h2_ndf<3&&h1_p>4&&h1_p<1000&&h2_p>4&&h2_p<1000&&h1_pt>0.250&&h2_pt>0.250&&\
                               h1_probnnghost<0.5&&h2_probnnghost<0.5&&TMath::Abs(z0_phi-jet_phi)>7*TMath::Pi()/8.&&TMath::Abs(jet_phi-mum_phi)>0.4&&TMath::Abs(jet_phi-mup_phi)>0.4&&\
                               mum_pt>20.&&mup_pt>20.&&mum_eta>2&&mum_eta<4.5&&mup_eta>2&&mup_eta<4.5&&sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))>60&&\
                               sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))<120&&mum_probchi2>0.001&&mup_probchi2>0.001)"
                               ,jet_pt_binning[2],jet_pt_binning[3]),
                               Form("weight*(jet_eta>2.5&&jet_eta<4.&&jet_pt>%f&&jet_pt<%f&&h1_chi2/h1_ndf<3&&h2_chi2/h2_ndf<3&&h1_p>4&&h1_p<1000&&h2_p>4&&h2_p<1000&&h1_pt>0.250&&h2_pt>0.250&&\
                               h1_probnnghost<0.5&&h2_probnnghost<0.5&&TMath::Abs(z0_phi-jet_phi)>7*TMath::Pi()/8.&&TMath::Abs(jet_phi-mum_phi)>0.4&&TMath::Abs(jet_phi-mup_phi)>0.4&&\
                               mum_pt>20.&&mup_pt>20.&&mum_eta>2&&mum_eta<4.5&&mup_eta>2&&mup_eta<4.5&&sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))>60&&\
                               sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))<120&&mum_probchi2>0.001&&mup_probchi2>0.001)"
                               ,jet_pt_binning[3],jet_pt_binning[4]),
                               Form("weight*(jet_eta>2.5&&jet_eta<4.&&jet_pt>%f&&jet_pt<%f&&h1_chi2/h1_ndf<3&&h2_chi2/h2_ndf<3&&h1_p>4&&h1_p<1000&&h2_p>4&&h2_p<1000&&h1_pt>0.250&&h2_pt>0.250&&\
                               h1_probnnghost<0.5&&h2_probnnghost<0.5&&TMath::Abs(z0_phi-jet_phi)>7*TMath::Pi()/8.&&TMath::Abs(jet_phi-mum_phi)>0.4&&TMath::Abs(jet_phi-mup_phi)>0.4&&\
                               mum_pt>20.&&mup_pt>20.&&mum_eta>2&&mum_eta<4.5&&mup_eta>2&&mup_eta<4.5&&sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))>60&&\
                               sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))<120&&mum_probchi2>0.001&&mup_probchi2>0.001)"
                               ,jet_pt_binning[4],jet_pt_binning[5])};

TCut pair_jetpt_signal_cut[] = {Form("(R_h1>%f&&R_h1<%f)&&(R_h2>%f&&R_h2<%f)&&jet_eta>2.5&&jet_eta<4.&&jet_pt>%f&&jet_pt<%f&&h1_chi2/h1_ndf<3&&h2_chi2/h2_ndf<3&&h1_p>4&&h1_p<1000&&h2_p>4&&h2_p<1000&&h1_pt>0.250&&h2_pt>0.250&&\
                           h1_probnnghost<0.5&&h2_probnnghost<0.5&&TMath::Abs(z0_phi-jet_phi)>7*TMath::Pi()/8.&&TMath::Abs(jet_phi-mum_phi)>0.4&&TMath::Abs(jet_phi-mup_phi)>0.4&&\
                           mum_pt>20.&&mup_pt>20.&&mum_eta>2&&mum_eta<4.5&&mup_eta>2&&mup_eta<4.5&&sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))>60&&\
                           sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))<120&&mum_probchi2>0.001&&mup_probchi2>0.001"
                           ,R_match_min,R_match_max,R_match_min,R_match_max,jet_pt_binning[0],jet_pt_binning[1]),
                           Form("(R_h1>%f&&R_h1<%f)&&(R_h2>%f&&R_h2<%f)&&jet_eta>2.5&&jet_eta<4.&&jet_pt>%f&&jet_pt<%f&&h1_chi2/h1_ndf<3&&h2_chi2/h2_ndf<3&&h1_p>4&&h1_p<1000&&h2_p>4&&h2_p<1000&&h1_pt>0.250&&h2_pt>0.250&&\
                           h1_probnnghost<0.5&&h2_probnnghost<0.5&&TMath::Abs(z0_phi-jet_phi)>7*TMath::Pi()/8.&&TMath::Abs(jet_phi-mum_phi)>0.4&&TMath::Abs(jet_phi-mup_phi)>0.4&&\
                           mum_pt>20.&&mup_pt>20.&&mum_eta>2&&mum_eta<4.5&&mup_eta>2&&mup_eta<4.5&&sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))>60&&\
                           sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))<120&&mum_probchi2>0.001&&mup_probchi2>0.001"
                           ,R_match_min,R_match_max,R_match_min,R_match_max,jet_pt_binning[1],jet_pt_binning[2]),
                           Form("(R_h1>%f&&R_h1<%f)&&(R_h2>%f&&R_h2<%f)&&jet_eta>2.5&&jet_eta<4.&&jet_pt>%f&&jet_pt<%f&&h1_chi2/h1_ndf<3&&h2_chi2/h2_ndf<3&&h1_p>4&&h1_p<1000&&h2_p>4&&h2_p<1000&&h1_pt>0.250&&h2_pt>0.250&&\
                           h1_probnnghost<0.5&&h2_probnnghost<0.5&&TMath::Abs(z0_phi-jet_phi)>7*TMath::Pi()/8.&&TMath::Abs(jet_phi-mum_phi)>0.4&&TMath::Abs(jet_phi-mup_phi)>0.4&&\
                           mum_pt>20.&&mup_pt>20.&&mum_eta>2&&mum_eta<4.5&&mup_eta>2&&mup_eta<4.5&&sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))>60&&\
                           sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))<120&&mum_probchi2>0.001&&mup_probchi2>0.001"
                           ,R_match_min,R_match_max,R_match_min,R_match_max,jet_pt_binning[2],jet_pt_binning[3]),
                           Form("(R_h1>%f&&R_h1<%f)&&(R_h2>%f&&R_h2<%f)&&jet_eta>2.5&&jet_eta<4.&&jet_pt>%f&&jet_pt<%f&&h1_chi2/h1_ndf<3&&h2_chi2/h2_ndf<3&&h1_p>4&&h1_p<1000&&h2_p>4&&h2_p<1000&&h1_pt>0.250&&h2_pt>0.250&&\
                           h1_probnnghost<0.5&&h2_probnnghost<0.5&&TMath::Abs(z0_phi-jet_phi)>7*TMath::Pi()/8.&&TMath::Abs(jet_phi-mum_phi)>0.4&&TMath::Abs(jet_phi-mup_phi)>0.4&&\
                           mum_pt>20.&&mup_pt>20.&&mum_eta>2&&mum_eta<4.5&&mup_eta>2&&mup_eta<4.5&&sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))>60&&\
                           sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))<120&&mum_probchi2>0.001&&mup_probchi2>0.001"
                           ,R_match_min,R_match_max,R_match_min,R_match_max,jet_pt_binning[3],jet_pt_binning[4]),
                           Form("(R_h1>%f&&R_h1<%f)&&(R_h2>%f&&R_h2<%f)&&jet_eta>2.5&&jet_eta<4.&&jet_pt>%f&&jet_pt<%f&&h1_chi2/h1_ndf<3&&h2_chi2/h2_ndf<3&&h1_p>4&&h1_p<1000&&h2_p>4&&h2_p<1000&&h1_pt>0.250&&h2_pt>0.250&&\
                           h1_probnnghost<0.5&&h2_probnnghost<0.5&&TMath::Abs(z0_phi-jet_phi)>7*TMath::Pi()/8.&&TMath::Abs(jet_phi-mum_phi)>0.4&&TMath::Abs(jet_phi-mup_phi)>0.4&&\
                           mum_pt>20.&&mup_pt>20.&&mum_eta>2&&mum_eta<4.5&&mup_eta>2&&mup_eta<4.5&&sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))>60&&\
                           sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))<120&&mum_probchi2>0.001&&mup_probchi2>0.001"
                           ,R_match_min,R_match_max,R_match_min,R_match_max,jet_pt_binning[4],jet_pt_binning[5])};

#endif