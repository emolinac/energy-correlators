#include <iostream>
#include "TZJetsMC.h"
#include "TZJetsMC.C"
#include "TZJetsMCReco.h"
#include "TZJetsMCReco.C"
#include "TZJets2016Data.h"
#include "TZJets2016Data.C"
#include "TZJets2017Data.h"
#include "TZJets2017Data.C"
#include "TZJets2018Data.h"
#include "TZJets2018Data.C"
#include "TROOT.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TH3.h"
#include "analysis-constants.h"
#include "analysis-binning.h"
#include "analysis-cuts.h"
#include "analysis-functions.h"
#include "directories.h"
#include "names.h"

int main()
{
  // Open correction files
  TFile* fpurity         = new TFile((output_folder + namef_ntuple_e2c_pairpurity).c_str());
  TFile* fefficiency     = new TFile((output_folder + namef_ntuple_e2c_pairefficiency).c_str());
  
  TFile* fpurity_jet     = new TFile((output_folder + namef_ntuple_jet_purity).c_str());
  TFile* fefficiency_jet = new TFile((output_folder + namef_ntuple_jet_efficiency).c_str());
  
  TFile* fefficiency_muon_2016_id  = new TFile((muons_folder + "IDEff_Data_2016.root").c_str());
  TFile* fefficiency_muon_2016_trk = new TFile((muons_folder + "TRKEff_Data_2016.root").c_str());
  TFile* fefficiency_muon_2016_trg = new TFile((muons_folder + "TRGEff_Data_2016.root").c_str());
  TFile* fefficiency_muon_2017_id  = new TFile((muons_folder + "IDEff_Data_2017.root").c_str());
  TFile* fefficiency_muon_2017_trk = new TFile((muons_folder + "TRKEff_Data_2017.root").c_str());
  TFile* fefficiency_muon_2017_trg = new TFile((muons_folder + "TRGEff_Data_2017.root").c_str());
  TFile* fefficiency_muon_2018_id  = new TFile((muons_folder + "IDEff_Data_2018.root").c_str());
  TFile* fefficiency_muon_2018_trk = new TFile((muons_folder + "TRKEff_Data_2018.root").c_str());
  TFile* fefficiency_muon_2018_trg = new TFile((muons_folder + "TRGEff_Data_2018.root").c_str());
  
  // Create output file
  TFile* fout = new TFile((output_folder + namef_ntuple_e2c_paircorr).c_str(),"RECREATE");
  
  // Declare the TTrees to be used to build the ntuples
  TZJets2016Data* datatree_2016 = new TZJets2016Data();
  TZJets2017Data* datatree_2017 = new TZJets2017Data();
  TZJets2018Data* datatree_2018 = new TZJets2018Data();
  
  // Create Ntuples
  TNtuple* ntuple_purity          = (TNtuple*) fpurity->Get((name_ntuple_purity.c_str()));
  TNtuple* ntuple_efficiency_mc   = (TNtuple*) fefficiency->Get((name_ntuple_correction_mc.c_str()));
  TNtuple* ntuple_efficiency_reco = (TNtuple*) fefficiency->Get((name_ntuple_correction_reco.c_str()));
  TNtuple* ntuple_purity_jet      = (TNtuple*) fpurity_jet->Get((name_ntuple_jetpurity.c_str()));
  TNtuple* ntuple_efficiency_jet  = (TNtuple*) fefficiency_jet->Get((name_ntuple_jetefficiency.c_str()));
  TNtuple* ntuple_data            = new TNtuple(name_ntuple_data.c_str(),"All Data",ntuple_paircorrdata_vars); 
  TNtuple* ntuple_corrjet         = new TNtuple(name_ntuple_corrjet.c_str(),"All Data",ntuple_jet_vars); 
  ntuple_data->SetAutoSave(0);
  ntuple_corrjet->SetAutoSave(0);

  // Muon corrections
  TH2D* h2_muon_2016_ideff_data  = (TH2D*) fefficiency_muon_2016_id->Get("Hist_ALL_2016_ETA_PT_Eff");
  TH2D* h2_muon_2016_trkeff_data = (TH2D*) fefficiency_muon_2016_trk->Get("Hist_ALL_2016_ETA_PT_Eff");
  TH2D* h2_muon_2016_trgeff_data = (TH2D*) fefficiency_muon_2016_trg->Get("Hist_ALL_2016_ETA_PT_Eff");
  TH2D* h2_muon_2017_ideff_data  = (TH2D*) fefficiency_muon_2017_id->Get("Hist_ALL_2017_ETA_PT_Eff");
  TH2D* h2_muon_2017_trkeff_data = (TH2D*) fefficiency_muon_2017_trk->Get("Hist_ALL_2017_ETA_PT_Eff");
  TH2D* h2_muon_2017_trgeff_data = (TH2D*) fefficiency_muon_2017_trg->Get("Hist_ALL_2017_ETA_PT_Eff");
  TH2D* h2_muon_2018_ideff_data  = (TH2D*) fefficiency_muon_2018_id->Get("Hist_ALL_2018_ETA_PT_Eff");
  TH2D* h2_muon_2018_trkeff_data = (TH2D*) fefficiency_muon_2018_trk->Get("Hist_ALL_2018_ETA_PT_Eff");
  TH2D* h2_muon_2018_trgeff_data = (TH2D*) fefficiency_muon_2018_trg->Get("Hist_ALL_2018_ETA_PT_Eff");

  // Jet corrections
  TH1F* hnum_pur_jet = new TH1F("hnum_pur_jet", "", Nbin_jetpt_corrections, corrections_jetpt_binning);
  TH1F* hden_pur_jet = new TH1F("hden_pur_jet", "", Nbin_jetpt_corrections, corrections_jetpt_binning);
  TH1F* hpurity_jet  = new TH1F("hpurity_jet" , "", Nbin_jetpt_corrections, corrections_jetpt_binning);
  hnum_pur_jet->Sumw2();
  hden_pur_jet->Sumw2();

  TH1F* hnum_eff_jet    = new TH1F("hnum_eff_jet"   , "", Nbin_jetpt_corrections, corrections_jetpt_binning);
  TH1F* hden_eff_jet    = new TH1F("hden_eff_jet"   , "", Nbin_jetpt_corrections, corrections_jetpt_binning);
  TH1F* hefficiency_jet = new TH1F("hefficiency_jet", "", Nbin_jetpt_corrections, corrections_jetpt_binning);
  hnum_eff_jet->Sumw2();
  hden_eff_jet->Sumw2();

  ntuple_purity_jet->Project("hnum_pur_jet", "jet_pt", "jet_pt_truth!=-999");
  ntuple_purity_jet->Project("hden_pur_jet", "jet_pt");
  ntuple_efficiency_jet->Project("hnum_eff_jet", "jet_pt_truth", "jet_pt!=-999");
  ntuple_efficiency_jet->Project("hden_eff_jet", "jet_pt_truth");

  hpurity_jet->Divide(hnum_pur_jet, hden_pur_jet, 1, 1, "B");
  hefficiency_jet->Divide(hnum_eff_jet, hden_eff_jet, 1, 1, "B");

  // Hadron corrections
  TH2F* hnum_pur    = new TH2F("hnum_pur"   , "", Nbin_R_L_nominal, rl_nominal_binning, Nbin_jetpt_corrections, corrections_jetpt_binning);
  TH2F* hden_pur    = new TH2F("hden_pur"   , "", Nbin_R_L_nominal, rl_nominal_binning, Nbin_jetpt_corrections, corrections_jetpt_binning);
  TH2F* hpurity     = new TH2F("hpurity"    , "", Nbin_R_L_nominal, rl_nominal_binning, Nbin_jetpt_corrections, corrections_jetpt_binning);
  TH2F* hnum_eff    = new TH2F("hnum_eff"   , "", Nbin_R_L_nominal, rl_nominal_binning, Nbin_jetpt_corrections, corrections_jetpt_binning);
  TH2F* hden_eff    = new TH2F("hden_eff"   , "", Nbin_R_L_nominal, rl_nominal_binning, Nbin_jetpt_corrections, corrections_jetpt_binning);
  TH2F* hefficiency = new TH2F("hefficiency", "", Nbin_R_L_nominal, rl_nominal_binning, Nbin_jetpt_corrections, corrections_jetpt_binning);
  
  hnum_pur->Sumw2();
  hden_pur->Sumw2();
  hnum_eff->Sumw2();
  hden_eff->Sumw2();
  
  ntuple_purity->Project("hnum_pur", "jet_pt:R_L", pair_matching_cut);
  ntuple_purity->Project("hden_pur", "jet_pt:R_L");
  ntuple_efficiency_reco->Project("hnum_eff", "jet_pt_truth:R_L_truth", pair_matching_cut);
  ntuple_efficiency_mc->Project("hden_eff", "jet_pt:R_L");

  hpurity->Divide(hnum_pur, hden_pur, 1, 1, "B");
  hefficiency->Divide(hnum_eff, hden_eff, 1, 1, "B");

  // DELETE LATER
  // DELETE LATER
  TH2F* hnum_pur_eqcharge    = new TH2F("hnum_pur_eqcharge"   , "", Nbin_R_L_nominal, rl_nominal_binning, Nbin_jetpt_corrections, corrections_jetpt_binning);
  TH2F* hden_pur_eqcharge    = new TH2F("hden_pur_eqcharge"   , "", Nbin_R_L_nominal, rl_nominal_binning, Nbin_jetpt_corrections, corrections_jetpt_binning);
  TH2F* hpurity_eqcharge     = new TH2F("hpurity_eqcharge"    , "", Nbin_R_L_nominal, rl_nominal_binning, Nbin_jetpt_corrections, corrections_jetpt_binning);
  TH2F* hnum_eff_eqcharge    = new TH2F("hnum_eff_eqcharge"   , "", Nbin_R_L_nominal, rl_nominal_binning, Nbin_jetpt_corrections, corrections_jetpt_binning);
  TH2F* hden_eff_eqcharge    = new TH2F("hden_eff_eqcharge"   , "", Nbin_R_L_nominal, rl_nominal_binning, Nbin_jetpt_corrections, corrections_jetpt_binning);
  TH2F* hefficiency_eqcharge = new TH2F("hefficiency_eqcharge", "", Nbin_R_L_nominal, rl_nominal_binning, Nbin_jetpt_corrections, corrections_jetpt_binning);
  
  TH2F* hnum_pur_neqcharge    = new TH2F("hnum_pur_neqcharge"   , "", Nbin_R_L_nominal, rl_nominal_binning, Nbin_jetpt_corrections, corrections_jetpt_binning);
  TH2F* hden_pur_neqcharge    = new TH2F("hden_pur_neqcharge"   , "", Nbin_R_L_nominal, rl_nominal_binning, Nbin_jetpt_corrections, corrections_jetpt_binning);
  TH2F* hpurity_neqcharge     = new TH2F("hpurity_neqcharge"    , "", Nbin_R_L_nominal, rl_nominal_binning, Nbin_jetpt_corrections, corrections_jetpt_binning);
  TH2F* hnum_eff_neqcharge    = new TH2F("hnum_eff_neqcharge"   , "", Nbin_R_L_nominal, rl_nominal_binning, Nbin_jetpt_corrections, corrections_jetpt_binning);
  TH2F* hden_eff_neqcharge    = new TH2F("hden_eff_neqcharge"   , "", Nbin_R_L_nominal, rl_nominal_binning, Nbin_jetpt_corrections, corrections_jetpt_binning);
  TH2F* hefficiency_neqcharge = new TH2F("hefficiency_neqcharge", "", Nbin_R_L_nominal, rl_nominal_binning, Nbin_jetpt_corrections, corrections_jetpt_binning);
  
  ntuple_purity->Project("hnum_pur_eqcharge", "jet_pt:R_L",pair_matching_cut + "eq_charge>1");
  ntuple_purity->Project("hden_pur_eqcharge", "jet_pt:R_L","eq_charge>1");
  ntuple_efficiency_reco->Project("hnum_eff_eqcharge", "jet_pt_truth:R_L_truth", pair_matching_cut + "eq_charge>1");
  ntuple_efficiency_mc->Project("hden_eff_eqcharge", "jet_pt:R_L","eq_charge>1");

  ntuple_purity->Project("hnum_pur_neqcharge", "jet_pt:R_L",pair_matching_cut + "eq_charge<-1");
  ntuple_purity->Project("hden_pur_neqcharge", "jet_pt:R_L","eq_charge<-1");
  ntuple_efficiency_reco->Project("hnum_eff_neqcharge", "jet_pt_truth:R_L_truth", pair_matching_cut + "eq_charge<-1");
  ntuple_efficiency_mc->Project("hden_eff_neqcharge", "jet_pt:R_L","eq_charge<-1");

  hpurity_eqcharge->Divide(hnum_pur_eqcharge, hden_pur_eqcharge, 1, 1, "B");
  hefficiency_eqcharge->Divide(hnum_eff_eqcharge, hden_eff_eqcharge, 1, 1, "B");

  hpurity_neqcharge->Divide(hnum_pur_neqcharge, hden_pur_neqcharge, 1, 1, "B");
  hefficiency_neqcharge->Divide(hnum_eff_neqcharge, hden_eff_neqcharge, 1, 1, "B");

  TCanvas* c = new TCanvas("c", "", 1920, 1080);
  c->Draw();

  gStyle->SetOptStat("");
  gStyle->SetPaintTextFormat("4.2f");
  hpurity->Draw("col text");
  hpurity->SetTitle("Purity Correction;R_{L};p^{jet}_{T}(GeV)");
  hpurity->GetXaxis()->SetRangeUser(R_L_min, R_L_max);
  hpurity->GetYaxis()->SetRangeUser(jet_pt_binning[0], jet_pt_binning[3]);
  gPad->SetLogx(1);
  gPad->SetLogy(1);
  c->Print("../src-analysis/plots/pair_purity_correction_logbin.pdf");

  hefficiency->Draw("col text");
  hefficiency->SetTitle("Efficiency Correction;R_{L};p^{jet}_{T}(GeV)");
  hefficiency->GetXaxis()->SetRangeUser(R_L_min, R_L_max);
  hefficiency->GetYaxis()->SetRangeUser(jet_pt_binning[0], jet_pt_binning[3]);
  gPad->SetLogx(1);
  gPad->SetLogy(1);
  c->Print("../src-analysis/plots/pair_efficiency_correction_logbin.pdf");

  // DELETE LATER!!!!
  hpurity_eqcharge->Draw("col text");
  hpurity_eqcharge->SetTitle("Purity Correction;R_{L};p^{jet}_{T}(GeV)");
  hpurity_eqcharge->GetXaxis()->SetRangeUser(R_L_min, R_L_max);
  hpurity_eqcharge->GetYaxis()->SetRangeUser(jet_pt_binning[0], jet_pt_binning[3]);
  gPad->SetLogx(1);
  gPad->SetLogy(1);
  c->Print("../src-analysis/plots/pair_purity_correction_eqcharge_logbin.pdf");

  hefficiency_eqcharge->Draw("col text");
  hefficiency_eqcharge->SetTitle("Efficiency Correction;R_{L};p^{jet}_{T}(GeV)");
  hefficiency_eqcharge->GetXaxis()->SetRangeUser(R_L_min, R_L_max);
  hefficiency_eqcharge->GetYaxis()->SetRangeUser(jet_pt_binning[0], jet_pt_binning[3]);
  gPad->SetLogx(1);
  gPad->SetLogy(1);
  c->Print("../src-analysis/plots/pair_efficiency_correction_eqcharge_logbin.pdf");

  hpurity_neqcharge->Draw("col text");
  hpurity_neqcharge->SetTitle("Purity Correction;R_{L};p^{jet}_{T}(GeV)");
  hpurity_neqcharge->GetXaxis()->SetRangeUser(R_L_min, R_L_max);
  hpurity_neqcharge->GetYaxis()->SetRangeUser(jet_pt_binning[0], jet_pt_binning[3]);
  gPad->SetLogx(1);
  gPad->SetLogy(1);
  c->Print("../src-analysis/plots/pair_purity_correction_neqcharge_logbin.pdf");

  hefficiency_neqcharge->Draw("col text");
  hefficiency_neqcharge->SetTitle("Efficiency Correction;R_{L};p^{jet}_{T}(GeV)");
  hefficiency_neqcharge->GetXaxis()->SetRangeUser(R_L_min, R_L_max);
  hefficiency_neqcharge->GetYaxis()->SetRangeUser(jet_pt_binning[0], jet_pt_binning[3]);
  gPad->SetLogx(1);
  gPad->SetLogy(1);
  c->Print("../src-analysis/plots/pair_efficiency_correction_neqcharge_logbin.pdf");
  
  // Create necessary 4vectors
  TLorentzVector* Jet_4vector = new TLorentzVector();
  TLorentzVector* Z0_4vector  = new TLorentzVector();
  TLorentzVector* mum_4vector = new TLorentzVector();
  TLorentzVector* mup_4vector = new TLorentzVector();
  TLorentzVector* h1_4vector  = new TLorentzVector();
  TLorentzVector* h2_4vector  = new TLorentzVector();
  int eventNum;
  unsigned long long last_eventNum = 0;
  int events = 0;
  bool maxjetpT_found = false;
    
  // Define array carrying the variables
  float vars[Nvars_paircorrdata];
  float vars_jet[Nvars_corrjet];

  // Fill the data TNtuple
  std::cout<<"Working with 2016 data."<<std::endl;
  for (int evt = 0 ; evt < datatree_2016->fChain->GetEntries() ; evt++)
  {
    // Access entry of tree
    datatree_2016->GetEntry(evt);
    if (evt%10000 == 0)
    {
      double percentage = 100.*evt/datatree_2016->fChain->GetEntries();
      std::cout<<"\r"<<percentage<<"\% jets processed."<< std::flush;
    }
    // Access entry of tree
    datatree_2016->GetEntry(evt);

    if (evt != 0)
    {
      if (last_eventNum == datatree_2016->eventNumber) continue;
    }

    // Apply PV cut
    if (datatree_2016->nPV != 1) continue;

    // Apply trigger cut
    bool mum_trigger = (datatree_2016->mum_L0MuonEWDecision_TOS == 1 && datatree_2016->mum_Hlt1SingleMuonHighPTDecision_TOS == 1 && datatree_2016->mum_Hlt2EWSingleMuonVHighPtDecision_TOS == 1);
    bool mup_trigger = (datatree_2016->mup_L0MuonEWDecision_TOS == 1 && datatree_2016->mup_Hlt1SingleMuonHighPTDecision_TOS == 1 && datatree_2016->mup_Hlt2EWSingleMuonVHighPtDecision_TOS == 1);

    if (!mum_trigger && !mup_trigger) continue;
    
    // Set Jet-associated 4 vectors and apply cuts
    Jet_4vector->SetPxPyPzE(datatree_2016->Jet_PX/1000.,datatree_2016->Jet_PY/1000.,datatree_2016->Jet_PZ/1000.,datatree_2016->Jet_PE/1000.);
    if (!apply_jet_cuts(Jet_4vector->Eta(), Jet_4vector->Pt())) continue;
    
    mum_4vector->SetPxPyPzE(datatree_2016->mum_PX/1000.,datatree_2016->mum_PY/1000.,datatree_2016->mum_PZ/1000.,datatree_2016->mum_PE/1000.);
    if (!apply_muon_cuts(Jet_4vector->DeltaR(*mum_4vector), mum_4vector->Pt(), mum_4vector->Eta())) continue;
    
    mup_4vector->SetPxPyPzE(datatree_2016->mup_PX/1000.,datatree_2016->mup_PY/1000.,datatree_2016->mup_PZ/1000.,datatree_2016->mup_PE/1000.);
    if (!apply_muon_cuts(Jet_4vector->DeltaR(*mup_4vector), mup_4vector->Pt(), mup_4vector->Eta())) continue;
    
    Z0_4vector->SetPxPyPzE(mup_4vector->Px()+mum_4vector->Px(),mup_4vector->Py()+mum_4vector->Py(),mup_4vector->Pz()+mum_4vector->Pz(),mup_4vector->E() +mum_4vector->E());
    if (!apply_zboson_cuts(TMath::Abs(Jet_4vector->DeltaPhi(*Z0_4vector)), Z0_4vector->M())) continue;

    double mup_pt  = (mup_4vector->Pt() >= 70.) ? 69. : mup_4vector->Pt();
    double mum_pt  = (mum_4vector->Pt() >= 70.) ? 69. : mum_4vector->Pt();
    double mup_eta = mup_4vector->Eta();
    double mum_eta = mum_4vector->Eta();

    double mup_eff_id  = h2_muon_2016_ideff_data->GetBinContent(h2_muon_2016_ideff_data->FindBin(mup_eta, mup_pt));
    double mup_eff_trk = h2_muon_2016_trkeff_data->GetBinContent(h2_muon_2016_trkeff_data->FindBin(mup_eta, mup_pt));
    double mup_eff_trg = h2_muon_2016_trgeff_data->GetBinContent(h2_muon_2016_trgeff_data->FindBin(mup_eta, mup_pt));

    double mum_eff_id  = h2_muon_2016_ideff_data->GetBinContent(h2_muon_2016_ideff_data->FindBin(mum_eta, mum_pt));
    double mum_eff_trk = h2_muon_2016_trkeff_data->GetBinContent(h2_muon_2016_trkeff_data->FindBin(mum_eta, mum_pt));
    double mum_eff_trg = h2_muon_2016_trgeff_data->GetBinContent(h2_muon_2016_trgeff_data->FindBin(mum_eta, mum_pt));

    double jet_efficiency = hefficiency_jet->GetBinContent(hefficiency_jet->FindBin(Jet_4vector->Pt()));
    double jet_purity     = hpurity_jet->GetBinContent(hpurity_jet->FindBin(Jet_4vector->Pt()));
    
    double jet_efficiency_error = hefficiency_jet->GetBinError(hefficiency_jet->FindBin(Jet_4vector->Pt()));
    double jet_purity_error     = hpurity_jet->GetBinError(hpurity_jet->FindBin(Jet_4vector->Pt()));
    
    vars_jet[0]  = Jet_4vector->Pt();
    vars_jet[1]  = Jet_4vector->E();
    vars_jet[2]  = datatree_2016->Jet_NDtr;
    vars_jet[3]  = jet_efficiency;
    vars_jet[4]  = jet_purity;
    vars_jet[5]  = jet_efficiency_error;
    vars_jet[6]  = jet_purity_error;
    vars_jet[7]  = mup_eff_id;
    vars_jet[8]  = mup_eff_trk;
    vars_jet[9]  = mup_eff_trg;
    vars_jet[10] = mum_eff_id;
    vars_jet[11] = mum_eff_trk;
    vars_jet[12] = mum_eff_trg;
    vars_jet[13] = 2016;
    vars_jet[14] = Jet_4vector->Eta();
    vars_jet[15] = Z0_4vector->Pt();
    vars_jet[16] = Z0_4vector->Eta();
    vars_jet[17] = Z0_4vector->Rapidity();
    
    ntuple_corrjet->Fill(vars_jet);

    // Loop over hadron 1
    for (int h1_index = 0 ; h1_index < datatree_2016->Jet_NDtr ; h1_index++)
    {
      // Skip non-hadronic particles
      if (datatree_2016->Jet_Dtr_IsMeson[h1_index] != 1 && datatree_2016->Jet_Dtr_IsBaryon[h1_index] != 1) continue;

      h1_4vector->SetPxPyPzE(datatree_2016->Jet_Dtr_PX[h1_index]/1000.,datatree_2016->Jet_Dtr_PY[h1_index]/1000.,datatree_2016->Jet_Dtr_PZ[h1_index]/1000.,datatree_2016->Jet_Dtr_E[h1_index]/1000.);
      if (!apply_chargedtrack_cuts(datatree_2016->Jet_Dtr_ThreeCharge[h1_index],
                                  h1_4vector->P(),
                                  h1_4vector->Pt(),
                                  datatree_2016->Jet_Dtr_TrackChi2[h1_index]/datatree_2016->Jet_Dtr_TrackNDF[h1_index],
                                  datatree_2016->Jet_Dtr_ProbNNghost[h1_index],
                                  h1_4vector->Eta())) continue;

      // Loop over hadron 2
      for (int h2_index = h1_index+1 ; h2_index < datatree_2016->Jet_NDtr ; h2_index++)
      {
        // Skip non-hadronic particles
        if (datatree_2016->Jet_Dtr_IsMeson[h2_index] != 1 && datatree_2016->Jet_Dtr_IsBaryon[h2_index] != 1) continue;

        h2_4vector->SetPxPyPzE(datatree_2016->Jet_Dtr_PX[h2_index]/1000.,datatree_2016->Jet_Dtr_PY[h2_index]/1000.,datatree_2016->Jet_Dtr_PZ[h2_index]/1000.,datatree_2016->Jet_Dtr_E[h2_index]/1000.);
        if (!apply_chargedtrack_cuts(datatree_2016->Jet_Dtr_ThreeCharge[h2_index],
                                    h2_4vector->P(),
                                    h2_4vector->Pt(),
                                    datatree_2016->Jet_Dtr_TrackChi2[h2_index]/datatree_2016->Jet_Dtr_TrackNDF[h2_index],
                                    datatree_2016->Jet_Dtr_ProbNNghost[h2_index],
                                    h2_4vector->Eta())) continue;

        double R_L = h1_4vector->DeltaR(*h2_4vector);

        double purity           = hpurity->GetBinContent(hpurity->FindBin(R_L, Jet_4vector->Pt()));        
        double efficiency       = hefficiency->GetBinContent(hefficiency->FindBin(R_L, Jet_4vector->Pt()));
        double purity_error     = hpurity->GetBinError(hpurity->FindBin(R_L, Jet_4vector->Pt()));        
        double efficiency_error = hefficiency->GetBinError(hefficiency->FindBin(R_L, Jet_4vector->Pt()));

        double nreco_ok  = hnum_pur->GetBinContent(hnum_pur->FindBin(R_L, Jet_4vector->Pt()));
        double nreco     = hden_pur->GetBinContent(hden_pur->FindBin(R_L, Jet_4vector->Pt()));
        double ntruth_ok = hnum_eff->GetBinContent(hnum_eff->FindBin(R_L, Jet_4vector->Pt()));
        double ntruth    = hden_eff->GetBinContent(hden_eff->FindBin(R_L, Jet_4vector->Pt()));
        
        vars[0 ] = weight(h1_4vector->E(), h2_4vector->E(), Jet_4vector->E());
        vars[1 ] = efficiency;
        vars[2 ] = purity;
        vars[3 ] = efficiency_error/efficiency;
        vars[4 ] = purity_error/purity;
        vars[5 ] = h1_4vector->DeltaR(*h2_4vector);
        vars[6 ] = h1_4vector->Eta();
        vars[7 ] = h2_4vector->Eta();
        vars[8 ] = h1_4vector->Rapidity();
        vars[9 ] = h2_4vector->Rapidity();
        vars[10] = h1_4vector->P();
        vars[11] = h2_4vector->P();
        vars[12] = h1_4vector->Pt();
        vars[13] = h2_4vector->Pt();
        vars[14] = Jet_4vector->Pt();
        vars[15] = Jet_4vector->Eta();
        vars[16] = weight(h1_4vector->Pt(), h2_4vector->Pt(), Jet_4vector->Pt());
        vars[17] = Jet_4vector->E();
        vars[18] = h1_4vector->E();
        vars[19] = h2_4vector->E();
        vars[20] = 2016;
        vars[21] = nreco_ok;
        vars[22] = nreco;
        vars[23] = ntruth_ok;
        vars[24] = ntruth;
        vars[25] = datatree_2016->Jet_Dtr_ThreeCharge[h1_index]*datatree_2016->Jet_Dtr_ThreeCharge[h2_index];
        
        // Fill the TNtuple
        ntuple_data->Fill(vars);
      }
    }

    last_eventNum = datatree_2016->eventNumber;
  }

  std::cout<<"Working with 2017 data."<<std::endl;
  for (int evt = 0 ; evt < datatree_2017->fChain->GetEntries() ; evt++)
  {
    // Access entry of tree
    datatree_2017->GetEntry(evt);
    if (evt%10000 == 0)
    {
      double percentage = 100.*evt/datatree_2017->fChain->GetEntries();
      std::cout<<"\r"<<percentage<<"\% jets processed."<< std::flush;
    }
    // Access entry of tree
    datatree_2017->GetEntry(evt);

    if (evt != 0)
    {
      if (last_eventNum == datatree_2017->eventNumber) continue;
    }

    // Apply PV cut
    if (datatree_2017->nPV != 1) continue;

    // Apply trigger cut
    bool mum_trigger = (datatree_2017->mum_L0MuonEWDecision_TOS == 1 && datatree_2017->mum_Hlt1SingleMuonHighPTDecision_TOS == 1 && datatree_2017->mum_Hlt2EWSingleMuonVHighPtDecision_TOS == 1);
    bool mup_trigger = (datatree_2017->mup_L0MuonEWDecision_TOS == 1 && datatree_2017->mup_Hlt1SingleMuonHighPTDecision_TOS == 1 && datatree_2017->mup_Hlt2EWSingleMuonVHighPtDecision_TOS == 1);

    if (!mum_trigger && !mup_trigger) continue;
    
    // Set Jet-associated 4 vectors and apply cuts
    Jet_4vector->SetPxPyPzE(datatree_2017->Jet_PX/1000.,datatree_2017->Jet_PY/1000.,datatree_2017->Jet_PZ/1000.,datatree_2017->Jet_PE/1000.);
    if (!apply_jet_cuts(Jet_4vector->Eta(), Jet_4vector->Pt())) continue;
    
    mum_4vector->SetPxPyPzE(datatree_2017->mum_PX/1000.,datatree_2017->mum_PY/1000.,datatree_2017->mum_PZ/1000.,datatree_2017->mum_PE/1000.);
    if (!apply_muon_cuts(Jet_4vector->DeltaR(*mum_4vector), mum_4vector->Pt(), mum_4vector->Eta())) continue;
    
    mup_4vector->SetPxPyPzE(datatree_2017->mup_PX/1000.,datatree_2017->mup_PY/1000.,datatree_2017->mup_PZ/1000.,datatree_2017->mup_PE/1000.);
    if (!apply_muon_cuts(Jet_4vector->DeltaR(*mup_4vector), mup_4vector->Pt(), mup_4vector->Eta())) continue;
    
    Z0_4vector->SetPxPyPzE(mup_4vector->Px()+mum_4vector->Px(),mup_4vector->Py()+mum_4vector->Py(),mup_4vector->Pz()+mum_4vector->Pz(),mup_4vector->E() +mum_4vector->E());
    if (!apply_zboson_cuts(TMath::Abs(Jet_4vector->DeltaPhi(*Z0_4vector)), Z0_4vector->M())) continue;

    double mup_pt  = (mup_4vector->Pt() >= 70.) ? 69. : mup_4vector->Pt();
    double mum_pt  = (mum_4vector->Pt() >= 70.) ? 69. : mum_4vector->Pt();
    double mup_eta = mup_4vector->Eta();
    double mum_eta = mum_4vector->Eta();

    double mup_eff_id  = h2_muon_2017_ideff_data->GetBinContent(h2_muon_2017_ideff_data->FindBin(mup_eta, mup_pt));
    double mup_eff_trk = h2_muon_2017_trkeff_data->GetBinContent(h2_muon_2017_trkeff_data->FindBin(mup_eta, mup_pt));
    double mup_eff_trg = h2_muon_2017_trgeff_data->GetBinContent(h2_muon_2017_trgeff_data->FindBin(mup_eta, mup_pt));

    double mum_eff_id  = h2_muon_2017_ideff_data->GetBinContent(h2_muon_2017_ideff_data->FindBin(mum_eta, mum_pt));
    double mum_eff_trk = h2_muon_2017_trkeff_data->GetBinContent(h2_muon_2017_trkeff_data->FindBin(mum_eta, mum_pt));
    double mum_eff_trg = h2_muon_2017_trgeff_data->GetBinContent(h2_muon_2017_trgeff_data->FindBin(mum_eta, mum_pt));

    double jet_efficiency = hefficiency_jet->GetBinContent(hefficiency_jet->FindBin(Jet_4vector->Pt()));
    double jet_purity     = hpurity_jet->GetBinContent(hpurity_jet->FindBin(Jet_4vector->Pt()));
    
    double jet_efficiency_error = hefficiency_jet->GetBinError(hefficiency_jet->FindBin(Jet_4vector->Pt()));
    double jet_purity_error     = hpurity_jet->GetBinError(hpurity_jet->FindBin(Jet_4vector->Pt()));
    
    vars_jet[0]  = Jet_4vector->Pt();
    vars_jet[1]  = Jet_4vector->E();
    vars_jet[2]  = datatree_2017->Jet_NDtr;
    vars_jet[3]  = jet_efficiency;
    vars_jet[4]  = jet_purity;
    vars_jet[5]  = jet_efficiency_error;
    vars_jet[6]  = jet_purity_error;
    vars_jet[7]  = mup_eff_id;
    vars_jet[8]  = mup_eff_trk;
    vars_jet[9]  = mup_eff_trg;
    vars_jet[10] = mum_eff_id;
    vars_jet[11] = mum_eff_trk;
    vars_jet[12] = mum_eff_trg;
    vars_jet[13] = 2017;
    vars_jet[14] = Jet_4vector->Eta();
    vars_jet[15] = Z0_4vector->Pt();
    vars_jet[16] = Z0_4vector->Eta();
    vars_jet[17] = Z0_4vector->Rapidity();
    
    ntuple_corrjet->Fill(vars_jet);

    // Loop over hadron 1
    for (int h1_index = 0 ; h1_index < datatree_2017->Jet_NDtr ; h1_index++)
    {
      // Skip non-hadronic particles
      if (datatree_2017->Jet_Dtr_IsMeson[h1_index] != 1 && datatree_2017->Jet_Dtr_IsBaryon[h1_index] != 1) continue;

      h1_4vector->SetPxPyPzE(datatree_2017->Jet_Dtr_PX[h1_index]/1000.,datatree_2017->Jet_Dtr_PY[h1_index]/1000.,datatree_2017->Jet_Dtr_PZ[h1_index]/1000.,datatree_2017->Jet_Dtr_E[h1_index]/1000.);
      if (!apply_chargedtrack_cuts(datatree_2017->Jet_Dtr_ThreeCharge[h1_index],
                                  h1_4vector->P(),
                                  h1_4vector->Pt(),
                                  datatree_2017->Jet_Dtr_TrackChi2[h1_index]/datatree_2017->Jet_Dtr_TrackNDF[h1_index],
                                  datatree_2017->Jet_Dtr_ProbNNghost[h1_index],
                                  h1_4vector->Eta())) continue;

      // Loop over hadron 2
      for (int h2_index = h1_index+1 ; h2_index < datatree_2017->Jet_NDtr ; h2_index++)
      {
        // Skip non-hadronic particles
        if (datatree_2017->Jet_Dtr_IsMeson[h2_index] != 1 && datatree_2017->Jet_Dtr_IsBaryon[h2_index] != 1) continue;

        h2_4vector->SetPxPyPzE(datatree_2017->Jet_Dtr_PX[h2_index]/1000.,datatree_2017->Jet_Dtr_PY[h2_index]/1000.,datatree_2017->Jet_Dtr_PZ[h2_index]/1000.,datatree_2017->Jet_Dtr_E[h2_index]/1000.);
        if (!apply_chargedtrack_cuts(datatree_2017->Jet_Dtr_ThreeCharge[h2_index],
                                    h2_4vector->P(),
                                    h2_4vector->Pt(),
                                    datatree_2017->Jet_Dtr_TrackChi2[h2_index]/datatree_2017->Jet_Dtr_TrackNDF[h2_index],
                                    datatree_2017->Jet_Dtr_ProbNNghost[h2_index],
                                    h2_4vector->Eta())) continue;

        double R_L = h1_4vector->DeltaR(*h2_4vector);

        double purity           = hpurity->GetBinContent(hpurity->FindBin(R_L, Jet_4vector->Pt()));        
        double efficiency       = hefficiency->GetBinContent(hefficiency->FindBin(R_L, Jet_4vector->Pt()));
        double purity_error     = hpurity->GetBinError(hpurity->FindBin(R_L, Jet_4vector->Pt()));        
        double efficiency_error = hefficiency->GetBinError(hefficiency->FindBin(R_L, Jet_4vector->Pt()));

        double nreco_ok  = hnum_pur->GetBinContent(hnum_pur->FindBin(R_L, Jet_4vector->Pt()));
        double nreco     = hden_pur->GetBinContent(hden_pur->FindBin(R_L, Jet_4vector->Pt()));
        double ntruth_ok = hnum_eff->GetBinContent(hnum_eff->FindBin(R_L, Jet_4vector->Pt()));
        double ntruth    = hden_eff->GetBinContent(hden_eff->FindBin(R_L, Jet_4vector->Pt()));
        
        vars[0 ] = weight(h1_4vector->E(), h2_4vector->E(), Jet_4vector->E());
        vars[1 ] = efficiency;
        vars[2 ] = purity;
        vars[3 ] = efficiency_error/efficiency;
        vars[4 ] = purity_error/purity;
        vars[5 ] = h1_4vector->DeltaR(*h2_4vector);
        vars[6 ] = h1_4vector->Eta();
        vars[7 ] = h2_4vector->Eta();
        vars[8 ] = h1_4vector->Rapidity();
        vars[9 ] = h2_4vector->Rapidity();
        vars[10] = h1_4vector->P();
        vars[11] = h2_4vector->P();
        vars[12] = h1_4vector->Pt();
        vars[13] = h2_4vector->Pt();
        vars[14] = Jet_4vector->Pt();
        vars[15] = Jet_4vector->Eta();
        vars[16] = weight(h1_4vector->Pt(), h2_4vector->Pt(), Jet_4vector->Pt());
        vars[17] = Jet_4vector->E();
        vars[18] = h1_4vector->E();
        vars[19] = h2_4vector->E();
        vars[20] = 2017;
        vars[21] = nreco_ok;
        vars[22] = nreco;
        vars[23] = ntruth_ok;
        vars[24] = ntruth;
        vars[25] = datatree_2017->Jet_Dtr_ThreeCharge[h1_index]*datatree_2017->Jet_Dtr_ThreeCharge[h2_index];
        
        // Fill the TNtuple
        ntuple_data->Fill(vars);
      }
    }

    last_eventNum = datatree_2017->eventNumber;
  }

  std::cout<<"Working with 2018 data."<<std::endl;
  for (int evt = 0 ; evt < datatree_2018->fChain->GetEntries() ; evt++)
  {
    // Access entry of tree
    datatree_2018->GetEntry(evt);
    if (evt%10000 == 0)
    {
      double percentage = 100.*evt/datatree_2018->fChain->GetEntries();
      std::cout<<"\r"<<percentage<<"\% jets processed."<< std::flush;
    }
    // Access entry of tree
    datatree_2018->GetEntry(evt);

    if (evt != 0)
    {
      if (last_eventNum == datatree_2018->eventNumber) continue;
    }

    // Apply PV cut
    if (datatree_2018->nPV != 1) continue;

    // Apply trigger cut
    bool mum_trigger = (datatree_2018->mum_L0MuonEWDecision_TOS == 1 && datatree_2018->mum_Hlt1SingleMuonHighPTDecision_TOS == 1 && datatree_2018->mum_Hlt2EWSingleMuonVHighPtDecision_TOS == 1);
    bool mup_trigger = (datatree_2018->mup_L0MuonEWDecision_TOS == 1 && datatree_2018->mup_Hlt1SingleMuonHighPTDecision_TOS == 1 && datatree_2018->mup_Hlt2EWSingleMuonVHighPtDecision_TOS == 1);

    if (!mum_trigger && !mup_trigger) continue;
    
    // Set Jet-associated 4 vectors and apply cuts
    Jet_4vector->SetPxPyPzE(datatree_2018->Jet_PX/1000.,datatree_2018->Jet_PY/1000.,datatree_2018->Jet_PZ/1000.,datatree_2018->Jet_PE/1000.);
    if (!apply_jet_cuts(Jet_4vector->Eta(), Jet_4vector->Pt())) continue;
    
    mum_4vector->SetPxPyPzE(datatree_2018->mum_PX/1000.,datatree_2018->mum_PY/1000.,datatree_2018->mum_PZ/1000.,datatree_2018->mum_PE/1000.);
    if (!apply_muon_cuts(Jet_4vector->DeltaR(*mum_4vector), mum_4vector->Pt(), mum_4vector->Eta())) continue;
    
    mup_4vector->SetPxPyPzE(datatree_2018->mup_PX/1000.,datatree_2018->mup_PY/1000.,datatree_2018->mup_PZ/1000.,datatree_2018->mup_PE/1000.);
    if (!apply_muon_cuts(Jet_4vector->DeltaR(*mup_4vector), mup_4vector->Pt(), mup_4vector->Eta())) continue;
    
    Z0_4vector->SetPxPyPzE(mup_4vector->Px()+mum_4vector->Px(),mup_4vector->Py()+mum_4vector->Py(),mup_4vector->Pz()+mum_4vector->Pz(),mup_4vector->E() +mum_4vector->E());
    if (!apply_zboson_cuts(TMath::Abs(Jet_4vector->DeltaPhi(*Z0_4vector)), Z0_4vector->M())) continue;

    double mup_pt  = (mup_4vector->Pt() >= 70.) ? 69. : mup_4vector->Pt();
    double mum_pt  = (mum_4vector->Pt() >= 70.) ? 69. : mum_4vector->Pt();
    double mup_eta = mup_4vector->Eta();
    double mum_eta = mum_4vector->Eta();

    double mup_eff_id  = h2_muon_2018_ideff_data->GetBinContent(h2_muon_2018_ideff_data->FindBin(mup_eta, mup_pt));
    double mup_eff_trk = h2_muon_2018_trkeff_data->GetBinContent(h2_muon_2018_trkeff_data->FindBin(mup_eta, mup_pt));
    double mup_eff_trg = h2_muon_2018_trgeff_data->GetBinContent(h2_muon_2018_trgeff_data->FindBin(mup_eta, mup_pt));

    double mum_eff_id  = h2_muon_2018_ideff_data->GetBinContent(h2_muon_2018_ideff_data->FindBin(mum_eta, mum_pt));
    double mum_eff_trk = h2_muon_2018_trkeff_data->GetBinContent(h2_muon_2018_trkeff_data->FindBin(mum_eta, mum_pt));
    double mum_eff_trg = h2_muon_2018_trgeff_data->GetBinContent(h2_muon_2018_trgeff_data->FindBin(mum_eta, mum_pt));

    double jet_efficiency = hefficiency_jet->GetBinContent(hefficiency_jet->FindBin(Jet_4vector->Pt()));
    double jet_purity     = hpurity_jet->GetBinContent(hpurity_jet->FindBin(Jet_4vector->Pt()));
    
    double jet_efficiency_error = hefficiency_jet->GetBinError(hefficiency_jet->FindBin(Jet_4vector->Pt()));
    double jet_purity_error     = hpurity_jet->GetBinError(hpurity_jet->FindBin(Jet_4vector->Pt()));
    
    vars_jet[0]  = Jet_4vector->Pt();
    vars_jet[1]  = Jet_4vector->E();
    vars_jet[2]  = datatree_2018->Jet_NDtr;
    vars_jet[3]  = jet_efficiency;
    vars_jet[4]  = jet_purity;
    vars_jet[5]  = jet_efficiency_error;
    vars_jet[6]  = jet_purity_error;
    vars_jet[7]  = mup_eff_id;
    vars_jet[8]  = mup_eff_trk;
    vars_jet[9]  = mup_eff_trg;
    vars_jet[10] = mum_eff_id;
    vars_jet[11] = mum_eff_trk;
    vars_jet[12] = mum_eff_trg;
    vars_jet[13] = 2018;
    vars_jet[14] = Jet_4vector->Eta();
    vars_jet[15] = Z0_4vector->Pt();
    vars_jet[16] = Z0_4vector->Eta();
    vars_jet[17] = Z0_4vector->Rapidity();
    
    ntuple_corrjet->Fill(vars_jet);

    // Loop over hadron 1
    for (int h1_index = 0 ; h1_index < datatree_2018->Jet_NDtr ; h1_index++)
    {
      // Skip non-hadronic particles
      if (datatree_2018->Jet_Dtr_IsMeson[h1_index] != 1 && datatree_2018->Jet_Dtr_IsBaryon[h1_index] != 1) continue;

      h1_4vector->SetPxPyPzE(datatree_2018->Jet_Dtr_PX[h1_index]/1000.,datatree_2018->Jet_Dtr_PY[h1_index]/1000.,datatree_2018->Jet_Dtr_PZ[h1_index]/1000.,datatree_2018->Jet_Dtr_E[h1_index]/1000.);
      if (!apply_chargedtrack_cuts(datatree_2018->Jet_Dtr_ThreeCharge[h1_index],
                                  h1_4vector->P(),
                                  h1_4vector->Pt(),
                                  datatree_2018->Jet_Dtr_TrackChi2[h1_index]/datatree_2018->Jet_Dtr_TrackNDF[h1_index],
                                  datatree_2018->Jet_Dtr_ProbNNghost[h1_index],
                                  h1_4vector->Eta())) continue;

      // Loop over hadron 2
      for (int h2_index = h1_index+1 ; h2_index < datatree_2018->Jet_NDtr ; h2_index++)
      {
        // Skip non-hadronic particles
        if (datatree_2018->Jet_Dtr_IsMeson[h2_index] != 1 && datatree_2018->Jet_Dtr_IsBaryon[h2_index] != 1) continue;

        h2_4vector->SetPxPyPzE(datatree_2018->Jet_Dtr_PX[h2_index]/1000.,datatree_2018->Jet_Dtr_PY[h2_index]/1000.,datatree_2018->Jet_Dtr_PZ[h2_index]/1000.,datatree_2018->Jet_Dtr_E[h2_index]/1000.);
        if (!apply_chargedtrack_cuts(datatree_2018->Jet_Dtr_ThreeCharge[h2_index],
                                    h2_4vector->P(),
                                    h2_4vector->Pt(),
                                    datatree_2018->Jet_Dtr_TrackChi2[h2_index]/datatree_2018->Jet_Dtr_TrackNDF[h2_index],
                                    datatree_2018->Jet_Dtr_ProbNNghost[h2_index],
                                    h2_4vector->Eta())) continue;

        double R_L = h1_4vector->DeltaR(*h2_4vector);

        double purity           = hpurity->GetBinContent(hpurity->FindBin(R_L, Jet_4vector->Pt()));        
        double efficiency       = hefficiency->GetBinContent(hefficiency->FindBin(R_L, Jet_4vector->Pt()));
        double purity_error     = hpurity->GetBinError(hpurity->FindBin(R_L, Jet_4vector->Pt()));        
        double efficiency_error = hefficiency->GetBinError(hefficiency->FindBin(R_L, Jet_4vector->Pt()));

        double nreco_ok  = hnum_pur->GetBinContent(hnum_pur->FindBin(R_L, Jet_4vector->Pt()));
        double nreco     = hden_pur->GetBinContent(hden_pur->FindBin(R_L, Jet_4vector->Pt()));
        double ntruth_ok = hnum_eff->GetBinContent(hnum_eff->FindBin(R_L, Jet_4vector->Pt()));
        double ntruth    = hden_eff->GetBinContent(hden_eff->FindBin(R_L, Jet_4vector->Pt()));
        
        vars[0 ] = weight(h1_4vector->E(), h2_4vector->E(), Jet_4vector->E());
        vars[1 ] = efficiency;
        vars[2 ] = purity;
        vars[3 ] = efficiency_error/efficiency;
        vars[4 ] = purity_error/purity;
        vars[5 ] = h1_4vector->DeltaR(*h2_4vector);
        vars[6 ] = h1_4vector->Eta();
        vars[7 ] = h2_4vector->Eta();
        vars[8 ] = h1_4vector->Rapidity();
        vars[9 ] = h2_4vector->Rapidity();
        vars[10] = h1_4vector->P();
        vars[11] = h2_4vector->P();
        vars[12] = h1_4vector->Pt();
        vars[13] = h2_4vector->Pt();
        vars[14] = Jet_4vector->Pt();
        vars[15] = Jet_4vector->Eta();
        vars[16] = weight(h1_4vector->Pt(), h2_4vector->Pt(), Jet_4vector->Pt());
        vars[17] = Jet_4vector->E();
        vars[18] = h1_4vector->E();
        vars[19] = h2_4vector->E();
        vars[20] = 2018;
        vars[21] = nreco_ok;
        vars[22] = nreco;
        vars[23] = ntruth_ok;
        vars[24] = ntruth;
        vars[25] = datatree_2018->Jet_Dtr_ThreeCharge[h1_index]*datatree_2018->Jet_Dtr_ThreeCharge[h2_index];
        
        // Fill the TNtuple
        ntuple_data->Fill(vars);
      }
    }

    last_eventNum = datatree_2018->eventNumber;
  }

  fout->cd();
  ntuple_data->Write();
  ntuple_corrjet->Write();
  fout->Close();
  
  std::cout<<std::endl;

  return 0;
}
