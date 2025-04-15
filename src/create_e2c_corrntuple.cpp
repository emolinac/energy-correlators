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
#include "analysis-cuts.h"
#include "analysis-functions.h"
#include "directories.h"
#include "names.h"

int main()
{
  // Open correction files
  TFile* fpurity         = new TFile((output_folder+namef_ntuple_e2c_purity).c_str());
  TFile* fefficiency     = new TFile((output_folder+namef_ntuple_e2c_efficiency).c_str());
  
  TFile* fpurity_jet     = new TFile((output_folder+namef_ntuple_jet_purity).c_str());
  TFile* fefficiency_jet = new TFile((output_folder+namef_ntuple_jet_efficiency).c_str());
  
  TFile* fefficiency_muon_2016_id  = new TFile((muons_folder+"IDEff_Data_2016.root").c_str());
  TFile* fefficiency_muon_2016_trk = new TFile((muons_folder+"TRKEff_Data_2016.root").c_str());
  TFile* fefficiency_muon_2016_trg = new TFile((muons_folder+"TRGEff_Data_2016.root").c_str());
  TFile* fefficiency_muon_2017_id  = new TFile((muons_folder+"IDEff_Data_2017.root").c_str());
  TFile* fefficiency_muon_2017_trk = new TFile((muons_folder+"TRKEff_Data_2017.root").c_str());
  TFile* fefficiency_muon_2017_trg = new TFile((muons_folder+"TRGEff_Data_2017.root").c_str());
  TFile* fefficiency_muon_2018_id  = new TFile((muons_folder+"IDEff_Data_2018.root").c_str());
  TFile* fefficiency_muon_2018_trk = new TFile((muons_folder+"TRKEff_Data_2018.root").c_str());
  TFile* fefficiency_muon_2018_trg = new TFile((muons_folder+"TRGEff_Data_2018.root").c_str());
  
  // Create output file
  TFile* fout = new TFile((output_folder+namef_ntuple_e2c_corr).c_str(),"RECREATE");
  
  // Declare the TTrees to be used to build the ntuples
  TZJets2016Data* datatree_2016 = new TZJets2016Data();
  TZJets2017Data* datatree_2017 = new TZJets2017Data();
  TZJets2018Data* datatree_2018 = new TZJets2018Data();
  
  // Create Ntuples
  TNtuple* ntuple_purity          = (TNtuple*) fpurity->Get((name_ntuple_purity.c_str()));
  TNtuple* ntuple_efficiency_mc   = (TNtuple*) fefficiency->Get((name_ntuple_efficiency_mc.c_str()));
  TNtuple* ntuple_efficiency_reco = (TNtuple*) fefficiency->Get((name_ntuple_efficiency_reco.c_str()));
  TNtuple* ntuple_purity_jet      = (TNtuple*) fpurity_jet->Get((name_ntuple_jetpurity.c_str()));
  TNtuple* ntuple_efficiency_jet  = (TNtuple*) fefficiency_jet->Get((name_ntuple_jetefficiency.c_str()));
  TNtuple* ntuple_data            = new TNtuple(name_ntuple_data.c_str(),"All Data",ntuple_corrdata_vars); 
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
  TH1F* hsigp_jet   = new TH1F("hsigp_jet"  ,"",5,unfolding_jetpt_binning);
  TH1F* hallp_jet   = new TH1F("hallp_jet"  ,"",5,unfolding_jetpt_binning);
  TH1F* hpurity_jet = new TH1F("hpurity_jet","",5,unfolding_jetpt_binning);
  hsigp_jet->Sumw2();
  hallp_jet->Sumw2();

  TH1F* hsigeff_jet     = new TH1F("hsigeff_jet"    ,"",5,unfolding_jetpt_binning);
  TH1F* halleff_jet     = new TH1F("halleff_jet"    ,"",5,unfolding_jetpt_binning);
  TH1F* hefficiency_jet = new TH1F("hefficiency_jet","",5,unfolding_jetpt_binning);
  hsigeff_jet->Sumw2();
  halleff_jet->Sumw2();

  ntuple_purity_jet->Project("hsigp_jet","jet_pt","jet_pt_truth!=-999");
  ntuple_purity_jet->Project("hallp_jet","jet_pt");
  ntuple_efficiency_jet->Project("hsigeff_jet","jet_pt_truth","jet_pt!=-999");
  ntuple_efficiency_jet->Project("halleff_jet","jet_pt_truth");

  hpurity_jet->Divide(hsigp_jet,hallp_jet,1,1,"B");
  hefficiency_jet->Divide(hsigeff_jet,halleff_jet,1,1,"B");

  // Hadron corrections
  TH3F* hsigp   = new TH3F("hsigp"  ,"",ic_p_nbins,ic_p_binning,sl_eta_nbins,sl_eta_binning,5,unfolding_jetpt_binning);
  TH3F* hallp   = new TH3F("hallp"  ,"",ic_p_nbins,ic_p_binning,sl_eta_nbins,sl_eta_binning,5,unfolding_jetpt_binning);
  TH3F* hpurity = new TH3F("hpurity","",ic_p_nbins,ic_p_binning,sl_eta_nbins,sl_eta_binning,5,unfolding_jetpt_binning);
  hsigp->Sumw2();
  hallp->Sumw2();
  
  TH3F* hsigeff     = new TH3F("hsigeff"    ,"",ic_p_nbins,ic_p_binning,sl_eta_nbins,sl_eta_binning,5,unfolding_jetpt_binning);
  TH3F* halleff     = new TH3F("halleff"    ,"",ic_p_nbins,ic_p_binning,sl_eta_nbins,sl_eta_binning,5,unfolding_jetpt_binning);
  TH3F* hefficiency = new TH3F("hefficiency","",ic_p_nbins,ic_p_binning,sl_eta_nbins,sl_eta_binning,5,unfolding_jetpt_binning);
  hsigeff->Sumw2();
  halleff->Sumw2();
  
  ntuple_purity->Project("hsigp","jet_pt:h_eta:h_p",single_signal_cut);
  ntuple_purity->Project("hallp","jet_pt:h_eta:h_p",pair_cut         );
  ntuple_efficiency_reco->Project("hsigeff","jet_pt_truth:h_eta_truth:h_p_truth",single_signal_cut);
  ntuple_efficiency_mc->Project("halleff","jet_pt:h_eta:h_p",pair_cut);

  hpurity->Divide(hsigp,hallp,1,1,"B");
  hefficiency->Divide(hsigeff,halleff,1,1,"B");

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
  float vars[Nvars_corrdata];
  float vars_jet[Nvars_corrjet];

  // Fill the data TNtuple
  std::cout<<"Working with 2016 data."<<std::endl;
  for(int evt = 0 ; evt < datatree_2016->fChain->GetEntries() ; evt++)
  {
    // Access entry of tree
    datatree_2016->GetEntry(evt);
    if(evt%10000==0)
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
    if(datatree_2016->nPV!=1) continue;

    // Apply trigger cut
    bool mum_trigger = (datatree_2016->mum_L0MuonEWDecision_TOS==1&&datatree_2016->mum_Hlt1SingleMuonHighPTDecision_TOS==1&&datatree_2016->mum_Hlt2EWSingleMuonVHighPtDecision_TOS==1);
    bool mup_trigger = (datatree_2016->mup_L0MuonEWDecision_TOS==1&&datatree_2016->mup_Hlt1SingleMuonHighPTDecision_TOS==1&&datatree_2016->mup_Hlt2EWSingleMuonVHighPtDecision_TOS==1);

    if(!mum_trigger&&!mup_trigger) continue;
    
    // Set Jet-associated 4 vectors and apply cuts
    Jet_4vector->SetPxPyPzE(datatree_2016->Jet_PX/1000.,datatree_2016->Jet_PY/1000.,datatree_2016->Jet_PZ/1000.,datatree_2016->Jet_PE/1000.);
    if(!apply_jet_cuts(Jet_4vector->Eta(),Jet_4vector->Pt())) continue;
    
    mum_4vector->SetPxPyPzE(datatree_2016->mum_PX/1000.,datatree_2016->mum_PY/1000.,datatree_2016->mum_PZ/1000.,datatree_2016->mum_PE/1000.);
    if(!apply_muon_cuts(Jet_4vector->DeltaR(*mum_4vector),mum_4vector->Pt(),mum_4vector->Eta())) continue;
    
    mup_4vector->SetPxPyPzE(datatree_2016->mup_PX/1000.,datatree_2016->mup_PY/1000.,datatree_2016->mup_PZ/1000.,datatree_2016->mup_PE/1000.);
    if(!apply_muon_cuts(Jet_4vector->DeltaR(*mup_4vector),mup_4vector->Pt(),mup_4vector->Eta())) continue;
    
    Z0_4vector->SetPxPyPzE(mup_4vector->Px()+mum_4vector->Px(),mup_4vector->Py()+mum_4vector->Py(),mup_4vector->Pz()+mum_4vector->Pz(),mup_4vector->E() +mum_4vector->E());
    if(!apply_zboson_cuts(TMath::Abs(Jet_4vector->DeltaPhi(*Z0_4vector)),Z0_4vector->M())) continue;

    double mup_eff_id  = h2_muon_2016_ideff_data->GetBinContent(h2_muon_2016_ideff_data->FindBin(mup_4vector->Eta(),mup_4vector->Pt()));
    double mup_eff_trk = h2_muon_2016_trkeff_data->GetBinContent(h2_muon_2016_trkeff_data->FindBin(mup_4vector->Eta(),mup_4vector->Pt()));
    double mup_eff_trg = h2_muon_2016_trgeff_data->GetBinContent(h2_muon_2016_trgeff_data->FindBin(mup_4vector->Eta(),mup_4vector->Pt()));

    double mum_eff_id  = h2_muon_2016_ideff_data->GetBinContent(h2_muon_2016_ideff_data->FindBin(mum_4vector->Eta(),mum_4vector->Pt()));
    double mum_eff_trk = h2_muon_2016_trkeff_data->GetBinContent(h2_muon_2016_trkeff_data->FindBin(mum_4vector->Eta(),mum_4vector->Pt()));
    double mum_eff_trg = h2_muon_2016_trgeff_data->GetBinContent(h2_muon_2016_trgeff_data->FindBin(mum_4vector->Eta(),mum_4vector->Pt()));

    double jet_efficiency = hefficiency_jet->GetBinContent(hefficiency_jet->FindBin(datatree_2016->Jet_PT/1000.));
    double jet_purity     = hpurity_jet->GetBinContent(hpurity_jet->FindBin(datatree_2016->Jet_PT/1000.));
    
    double jet_efficiency_error = hefficiency_jet->GetBinError(hefficiency_jet->FindBin(datatree_2016->Jet_PT/1000.));
    double jet_purity_error     = hpurity_jet->GetBinError(hpurity_jet->FindBin(datatree_2016->Jet_PT/1000.));
    
    vars_jet[0]  = datatree_2016->Jet_PT/1000.;
    vars_jet[1]  = datatree_2016->Jet_PE/1000.;
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
    for(int h1_index = 0 ; h1_index < datatree_2016->Jet_NDtr ; h1_index++)
    {
      // Skip non-hadronic particles
      if(datatree_2016->Jet_Dtr_IsMeson[h1_index]!=1&&datatree_2016->Jet_Dtr_IsBaryon[h1_index]!=1) continue;

      h1_4vector->SetPxPyPzE(datatree_2016->Jet_Dtr_PX[h1_index]/1000.,datatree_2016->Jet_Dtr_PY[h1_index]/1000.,datatree_2016->Jet_Dtr_PZ[h1_index]/1000.,datatree_2016->Jet_Dtr_E[h1_index]/1000.);
      if(!apply_chargedtrack_cuts(datatree_2016->Jet_Dtr_ThreeCharge[h1_index],
                                  datatree_2016->Jet_Dtr_P[h1_index]/1000.,
                                  datatree_2016->Jet_Dtr_PT[h1_index]/1000.,
                                  datatree_2016->Jet_Dtr_TrackChi2[h1_index]/datatree_2016->Jet_Dtr_TrackNDF[h1_index],
                                  datatree_2016->Jet_Dtr_ProbNNghost[h1_index],
                                  Jet_4vector->DeltaR(*h1_4vector))) continue;

      double h1_purity     = hpurity->GetBinContent(hpurity->FindBin(datatree_2016->Jet_Dtr_P[h1_index]/1000.,datatree_2016->Jet_Dtr_ETA[h1_index],datatree_2016->Jet_PT/1000.));
      double h1_efficiency = hefficiency->GetBinContent(hefficiency->FindBin(datatree_2016->Jet_Dtr_P[h1_index]/1000.,datatree_2016->Jet_Dtr_ETA[h1_index],datatree_2016->Jet_PT/1000.));
      if(h1_purity>1.||h1_efficiency>1.) continue;

      // Loop over hadron 2
      for(int h2_index = h1_index+1 ; h2_index < datatree_2016->Jet_NDtr ; h2_index++)
      {
        // Skip non-hadronic particles
        if(datatree_2016->Jet_Dtr_IsMeson[h2_index]!=1&&datatree_2016->Jet_Dtr_IsBaryon[h2_index]!=1) continue;

        h2_4vector->SetPxPyPzE(datatree_2016->Jet_Dtr_PX[h2_index]/1000.,datatree_2016->Jet_Dtr_PY[h2_index]/1000.,datatree_2016->Jet_Dtr_PZ[h2_index]/1000.,datatree_2016->Jet_Dtr_E[h2_index]/1000.);
        if(!apply_chargedtrack_cuts(datatree_2016->Jet_Dtr_ThreeCharge[h2_index],
                                    datatree_2016->Jet_Dtr_P[h2_index]/1000.,
                                    datatree_2016->Jet_Dtr_PT[h2_index]/1000.,
                                    datatree_2016->Jet_Dtr_TrackChi2[h2_index]/datatree_2016->Jet_Dtr_TrackNDF[h2_index],
                                    datatree_2016->Jet_Dtr_ProbNNghost[h2_index],
                                    Jet_4vector->DeltaR(*h2_4vector))) continue;

        double h2_purity     = hpurity->GetBinContent(hpurity->FindBin(datatree_2016->Jet_Dtr_P[h2_index]/1000.,datatree_2016->Jet_Dtr_ETA[h2_index],datatree_2016->Jet_PT/1000.));
        double h2_efficiency = hefficiency->GetBinContent(hefficiency->FindBin(datatree_2016->Jet_Dtr_P[h2_index]/1000.,datatree_2016->Jet_Dtr_ETA[h2_index],datatree_2016->Jet_PT/1000.));
        if(h2_purity>1.||h2_efficiency>1.) continue;

        double purity_correction = (h1_purity)*(h2_purity);

        double h1_purity_err = hpurity->GetBinError(hpurity->FindBin(datatree_2016->Jet_Dtr_P[h1_index]/1000.,datatree_2016->Jet_Dtr_ETA[h1_index],datatree_2016->Jet_PT/1000.));
        double h2_purity_err = hpurity->GetBinError(hpurity->FindBin(datatree_2016->Jet_Dtr_P[h2_index]/1000.,datatree_2016->Jet_Dtr_ETA[h2_index],datatree_2016->Jet_PT/1000.));
        double purity_error = sqrt(pow((h1_purity)*(h2_purity_err),2) + pow((h1_purity_err)*(h2_purity),2));

        double efficiency_correction = (h1_efficiency)*(h2_efficiency);

        double h1_efficiency_err = hefficiency->GetBinError(hefficiency->FindBin(datatree_2016->Jet_Dtr_P[h1_index]/1000.,datatree_2016->Jet_Dtr_ETA[h1_index],datatree_2016->Jet_PT/1000.));
        double h2_efficiency_err = hefficiency->GetBinError(hefficiency->FindBin(datatree_2016->Jet_Dtr_P[h2_index]/1000.,datatree_2016->Jet_Dtr_ETA[h2_index],datatree_2016->Jet_PT/1000.));
        double efficiency_error = sqrt(pow((h1_efficiency)*(h2_efficiency_err),2) + pow((h1_efficiency_err)*(h2_efficiency),2));

        vars[0 ] = weight(datatree_2016->Jet_Dtr_E[h1_index], datatree_2016->Jet_Dtr_E[h2_index], datatree_2016->Jet_PE);
        vars[1 ] = efficiency_correction;
        vars[2 ] = purity_correction;
        vars[3 ] = efficiency_error/efficiency_correction;
        vars[4 ] = purity_error/purity_correction;
        vars[5 ] = h1_4vector->DeltaR(*h2_4vector);
        vars[6 ] = datatree_2016->Jet_Dtr_ETA[h1_index];
        vars[7 ] = datatree_2016->Jet_Dtr_ETA[h2_index];
        vars[8 ] = h1_4vector->Rapidity();
        vars[9 ] = h2_4vector->Rapidity();
        vars[10] = datatree_2016->Jet_Dtr_P[h1_index]/1000.;
        vars[11] = datatree_2016->Jet_Dtr_P[h2_index]/1000.;
        vars[12] = datatree_2016->Jet_Dtr_PT[h1_index]/1000.;
        vars[13] = datatree_2016->Jet_Dtr_PT[h2_index]/1000.;
        vars[14] = datatree_2016->Jet_PT/1000.;
        vars[15] = Jet_4vector->Eta();
        vars[16] = weight(datatree_2016->Jet_Dtr_PT[h1_index], datatree_2016->Jet_Dtr_PT[h2_index], datatree_2016->Jet_PT);
        vars[17] = datatree_2016->Jet_PE/1000.;
        vars[18] = datatree_2016->Jet_Dtr_E[h1_index]/1000.;
        vars[19] = datatree_2016->Jet_Dtr_E[h2_index]/1000.;
        vars[20] = 2016;
        
        // Fill the TNtuple
        ntuple_data->Fill(vars);
      }
    }

    last_eventNum = datatree_2016->eventNumber;
  }

  std::cout<<"Working with 2017 data."<<std::endl;
  for(int evt = 0 ; evt < datatree_2017->fChain->GetEntries() ; evt++)
  {
    // Access entry of tree
    datatree_2017->GetEntry(evt);
    if(evt%10000==0)
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
    if(datatree_2017->nPV!=1) continue;

    // Apply trigger cut
    bool mum_trigger = (datatree_2017->mum_L0MuonEWDecision_TOS==1&&datatree_2017->mum_Hlt1SingleMuonHighPTDecision_TOS==1&&datatree_2017->mum_Hlt2EWSingleMuonVHighPtDecision_TOS==1);
    bool mup_trigger = (datatree_2017->mup_L0MuonEWDecision_TOS==1&&datatree_2017->mup_Hlt1SingleMuonHighPTDecision_TOS==1&&datatree_2017->mup_Hlt2EWSingleMuonVHighPtDecision_TOS==1);

    if(!mum_trigger&&!mup_trigger) continue;
    
    // Set Jet-associated 4 vectors and apply cuts
    Jet_4vector->SetPxPyPzE(datatree_2017->Jet_PX/1000.,datatree_2017->Jet_PY/1000.,datatree_2017->Jet_PZ/1000.,datatree_2017->Jet_PE/1000.);
    if(!apply_jet_cuts(Jet_4vector->Eta(),Jet_4vector->Pt())) continue;
    
    mum_4vector->SetPxPyPzE(datatree_2017->mum_PX/1000.,datatree_2017->mum_PY/1000.,datatree_2017->mum_PZ/1000.,datatree_2017->mum_PE/1000.);
    if(!apply_muon_cuts(Jet_4vector->DeltaR(*mum_4vector),mum_4vector->Pt(),mum_4vector->Eta())) continue;
    
    mup_4vector->SetPxPyPzE(datatree_2017->mup_PX/1000.,datatree_2017->mup_PY/1000.,datatree_2017->mup_PZ/1000.,datatree_2017->mup_PE/1000.);
    if(!apply_muon_cuts(Jet_4vector->DeltaR(*mup_4vector),mup_4vector->Pt(),mup_4vector->Eta())) continue;
    
    Z0_4vector->SetPxPyPzE(mup_4vector->Px()+mum_4vector->Px(),mup_4vector->Py()+mum_4vector->Py(),mup_4vector->Pz()+mum_4vector->Pz(),mup_4vector->E() +mum_4vector->E());
    if(!apply_zboson_cuts(TMath::Abs(Jet_4vector->DeltaPhi(*Z0_4vector)),Z0_4vector->M())) continue;

    double mup_eff_id  = h2_muon_2017_ideff_data->GetBinContent(h2_muon_2017_ideff_data->FindBin(mup_4vector->Eta(),mup_4vector->Pt()));
    double mup_eff_trk = h2_muon_2017_trkeff_data->GetBinContent(h2_muon_2017_trkeff_data->FindBin(mup_4vector->Eta(),mup_4vector->Pt()));
    double mup_eff_trg = h2_muon_2017_trgeff_data->GetBinContent(h2_muon_2017_trgeff_data->FindBin(mup_4vector->Eta(),mup_4vector->Pt()));

    double mum_eff_id  = h2_muon_2017_ideff_data->GetBinContent(h2_muon_2017_ideff_data->FindBin(mum_4vector->Eta(),mum_4vector->Pt()));
    double mum_eff_trk = h2_muon_2017_trkeff_data->GetBinContent(h2_muon_2017_trkeff_data->FindBin(mum_4vector->Eta(),mum_4vector->Pt()));
    double mum_eff_trg = h2_muon_2017_trgeff_data->GetBinContent(h2_muon_2017_trgeff_data->FindBin(mum_4vector->Eta(),mum_4vector->Pt()));

    double jet_efficiency = hefficiency_jet->GetBinContent(hefficiency_jet->FindBin(datatree_2017->Jet_PT/1000.));
    double jet_purity     = hpurity_jet->GetBinContent(hpurity_jet->FindBin(datatree_2017->Jet_PT/1000.));
    
    double jet_efficiency_error = hefficiency_jet->GetBinError(hefficiency_jet->FindBin(datatree_2017->Jet_PT/1000.));
    double jet_purity_error     = hpurity_jet->GetBinError(hpurity_jet->FindBin(datatree_2017->Jet_PT/1000.));
    
    vars_jet[0]  = datatree_2017->Jet_PT/1000.;
    vars_jet[1]  = datatree_2017->Jet_PE/1000.;
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
    for(int h1_index = 0 ; h1_index < datatree_2017->Jet_NDtr ; h1_index++)
    {
      // Skip non-hadronic particles
      if(datatree_2017->Jet_Dtr_IsMeson[h1_index]!=1&&datatree_2017->Jet_Dtr_IsBaryon[h1_index]!=1) continue;

      h1_4vector->SetPxPyPzE(datatree_2017->Jet_Dtr_PX[h1_index]/1000.,datatree_2017->Jet_Dtr_PY[h1_index]/1000.,datatree_2017->Jet_Dtr_PZ[h1_index]/1000.,datatree_2017->Jet_Dtr_E[h1_index]/1000.);
      if(!apply_chargedtrack_cuts(datatree_2017->Jet_Dtr_ThreeCharge[h1_index],
                                  datatree_2017->Jet_Dtr_P[h1_index]/1000.,
                                  datatree_2017->Jet_Dtr_PT[h1_index]/1000.,
                                  datatree_2017->Jet_Dtr_TrackChi2[h1_index]/datatree_2017->Jet_Dtr_TrackNDF[h1_index],
                                  datatree_2017->Jet_Dtr_ProbNNghost[h1_index],
                                  Jet_4vector->DeltaR(*h1_4vector))) continue;

      double h1_purity     = hpurity->GetBinContent(hpurity->FindBin(datatree_2017->Jet_Dtr_P[h1_index]/1000.,datatree_2017->Jet_Dtr_ETA[h1_index],datatree_2017->Jet_PT/1000.));
      double h1_efficiency = hefficiency->GetBinContent(hefficiency->FindBin(datatree_2017->Jet_Dtr_P[h1_index]/1000.,datatree_2017->Jet_Dtr_ETA[h1_index],datatree_2017->Jet_PT/1000.));
      if(h1_purity>1.||h1_efficiency>1.) continue;

      // Loop over hadron 2
      for(int h2_index = h1_index+1 ; h2_index < datatree_2017->Jet_NDtr ; h2_index++)
      {
        // Skip non-hadronic particles
        if(datatree_2017->Jet_Dtr_IsMeson[h2_index]!=1&&datatree_2017->Jet_Dtr_IsBaryon[h2_index]!=1) continue;

        h2_4vector->SetPxPyPzE(datatree_2017->Jet_Dtr_PX[h2_index]/1000.,datatree_2017->Jet_Dtr_PY[h2_index]/1000.,datatree_2017->Jet_Dtr_PZ[h2_index]/1000.,datatree_2017->Jet_Dtr_E[h2_index]/1000.);
        if(!apply_chargedtrack_cuts(datatree_2017->Jet_Dtr_ThreeCharge[h2_index],
                                    datatree_2017->Jet_Dtr_P[h2_index]/1000.,
                                    datatree_2017->Jet_Dtr_PT[h2_index]/1000.,
                                    datatree_2017->Jet_Dtr_TrackChi2[h2_index]/datatree_2017->Jet_Dtr_TrackNDF[h2_index],
                                    datatree_2017->Jet_Dtr_ProbNNghost[h2_index],
                                    Jet_4vector->DeltaR(*h2_4vector))) continue;

        double h2_purity     = hpurity->GetBinContent(hpurity->FindBin(datatree_2017->Jet_Dtr_P[h2_index]/1000.,datatree_2017->Jet_Dtr_ETA[h2_index],datatree_2017->Jet_PT/1000.));
        double h2_efficiency = hefficiency->GetBinContent(hefficiency->FindBin(datatree_2017->Jet_Dtr_P[h2_index]/1000.,datatree_2017->Jet_Dtr_ETA[h2_index],datatree_2017->Jet_PT/1000.));
        if(h2_purity>1.||h2_efficiency>1.) continue;

        double purity_correction = (h1_purity)*(h2_purity);

        double h1_purity_err = hpurity->GetBinError(hpurity->FindBin(datatree_2017->Jet_Dtr_P[h1_index]/1000.,datatree_2017->Jet_Dtr_ETA[h1_index],datatree_2017->Jet_PT/1000.));
        double h2_purity_err = hpurity->GetBinError(hpurity->FindBin(datatree_2017->Jet_Dtr_P[h2_index]/1000.,datatree_2017->Jet_Dtr_ETA[h2_index],datatree_2017->Jet_PT/1000.));
        double purity_error = sqrt(pow((h1_purity)*(h2_purity_err),2) + pow((h1_purity_err)*(h2_purity),2));

        double efficiency_correction = (h1_efficiency)*(h2_efficiency);

        double h1_efficiency_err = hefficiency->GetBinError(hefficiency->FindBin(datatree_2017->Jet_Dtr_P[h1_index]/1000.,datatree_2017->Jet_Dtr_ETA[h1_index],datatree_2017->Jet_PT/1000.));
        double h2_efficiency_err = hefficiency->GetBinError(hefficiency->FindBin(datatree_2017->Jet_Dtr_P[h2_index]/1000.,datatree_2017->Jet_Dtr_ETA[h2_index],datatree_2017->Jet_PT/1000.));
        double efficiency_error = sqrt(pow((h1_efficiency)*(h2_efficiency_err),2) + pow((h1_efficiency_err)*(h2_efficiency),2));

        vars[0 ] = weight(datatree_2017->Jet_Dtr_E[h1_index], datatree_2017->Jet_Dtr_E[h2_index], datatree_2017->Jet_PE);
        vars[1 ] = efficiency_correction;
        vars[2 ] = purity_correction;
        vars[3 ] = efficiency_error/efficiency_correction;
        vars[4 ] = purity_error/purity_correction;
        vars[5 ] = h1_4vector->DeltaR(*h2_4vector);
        vars[6 ] = datatree_2017->Jet_Dtr_ETA[h1_index];
        vars[7 ] = datatree_2017->Jet_Dtr_ETA[h2_index];
        vars[8 ] = h1_4vector->Rapidity();
        vars[9 ] = h2_4vector->Rapidity();
        vars[10] = datatree_2017->Jet_Dtr_P[h1_index]/1000.;
        vars[11] = datatree_2017->Jet_Dtr_P[h2_index]/1000.;
        vars[12] = datatree_2017->Jet_Dtr_PT[h1_index]/1000.;
        vars[13] = datatree_2017->Jet_Dtr_PT[h2_index]/1000.;
        vars[14] = datatree_2017->Jet_PT/1000.;
        vars[15] = Jet_4vector->Eta();
        vars[16] = weight(datatree_2017->Jet_Dtr_PT[h1_index], datatree_2017->Jet_Dtr_PT[h2_index], datatree_2017->Jet_PT);
        vars[17] = datatree_2017->Jet_PE/1000.;
        vars[18] = datatree_2017->Jet_Dtr_E[h1_index]/1000.;
        vars[19] = datatree_2017->Jet_Dtr_E[h2_index]/1000.;
        vars[20] = 2017;
        
        // Fill the TNtuple
        ntuple_data->Fill(vars);
      }
    }

    last_eventNum = datatree_2017->eventNumber;
  }

  std::cout<<"Working with 2018 data."<<std::endl;
  for(int evt = 0 ; evt < datatree_2018->fChain->GetEntries() ; evt++)
  {
    // Access entry of tree
    datatree_2018->GetEntry(evt);
    if(evt%10000==0)
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
    if(datatree_2018->nPV!=1) continue;

    // Apply trigger cut
    bool mum_trigger = (datatree_2018->mum_L0MuonEWDecision_TOS==1&&datatree_2018->mum_Hlt1SingleMuonHighPTDecision_TOS==1&&datatree_2018->mum_Hlt2EWSingleMuonVHighPtDecision_TOS==1);
    bool mup_trigger = (datatree_2018->mup_L0MuonEWDecision_TOS==1&&datatree_2018->mup_Hlt1SingleMuonHighPTDecision_TOS==1&&datatree_2018->mup_Hlt2EWSingleMuonVHighPtDecision_TOS==1);

    if(!mum_trigger&&!mup_trigger) continue;
    
    // Set Jet-associated 4 vectors and apply cuts
    Jet_4vector->SetPxPyPzE(datatree_2018->Jet_PX/1000.,datatree_2018->Jet_PY/1000.,datatree_2018->Jet_PZ/1000.,datatree_2018->Jet_PE/1000.);
    if(!apply_jet_cuts(Jet_4vector->Eta(),Jet_4vector->Pt())) continue;
    
    mum_4vector->SetPxPyPzE(datatree_2018->mum_PX/1000.,datatree_2018->mum_PY/1000.,datatree_2018->mum_PZ/1000.,datatree_2018->mum_PE/1000.);
    if(!apply_muon_cuts(Jet_4vector->DeltaR(*mum_4vector),mum_4vector->Pt(),mum_4vector->Eta())) continue;
    
    mup_4vector->SetPxPyPzE(datatree_2018->mup_PX/1000.,datatree_2018->mup_PY/1000.,datatree_2018->mup_PZ/1000.,datatree_2018->mup_PE/1000.);
    if(!apply_muon_cuts(Jet_4vector->DeltaR(*mup_4vector),mup_4vector->Pt(),mup_4vector->Eta())) continue;
    
    Z0_4vector->SetPxPyPzE(mup_4vector->Px()+mum_4vector->Px(),mup_4vector->Py()+mum_4vector->Py(),mup_4vector->Pz()+mum_4vector->Pz(),mup_4vector->E() +mum_4vector->E());
    if(!apply_zboson_cuts(TMath::Abs(Jet_4vector->DeltaPhi(*Z0_4vector)),Z0_4vector->M())) continue;

    double mup_eff_id  = h2_muon_2018_ideff_data->GetBinContent(h2_muon_2018_ideff_data->FindBin(mup_4vector->Eta(),mup_4vector->Pt()));
    double mup_eff_trk = h2_muon_2018_trkeff_data->GetBinContent(h2_muon_2018_trkeff_data->FindBin(mup_4vector->Eta(),mup_4vector->Pt()));
    double mup_eff_trg = h2_muon_2018_trgeff_data->GetBinContent(h2_muon_2018_trgeff_data->FindBin(mup_4vector->Eta(),mup_4vector->Pt()));

    double mum_eff_id  = h2_muon_2018_ideff_data->GetBinContent(h2_muon_2018_ideff_data->FindBin(mum_4vector->Eta(),mum_4vector->Pt()));
    double mum_eff_trk = h2_muon_2018_trkeff_data->GetBinContent(h2_muon_2018_trkeff_data->FindBin(mum_4vector->Eta(),mum_4vector->Pt()));
    double mum_eff_trg = h2_muon_2018_trgeff_data->GetBinContent(h2_muon_2018_trgeff_data->FindBin(mum_4vector->Eta(),mum_4vector->Pt()));

    double jet_efficiency = hefficiency_jet->GetBinContent(hefficiency_jet->FindBin(datatree_2018->Jet_PT/1000.));
    double jet_purity     = hpurity_jet->GetBinContent(hpurity_jet->FindBin(datatree_2018->Jet_PT/1000.));
    
    double jet_efficiency_error = hefficiency_jet->GetBinError(hefficiency_jet->FindBin(datatree_2018->Jet_PT/1000.));
    double jet_purity_error     = hpurity_jet->GetBinError(hpurity_jet->FindBin(datatree_2018->Jet_PT/1000.));
    
    vars_jet[0]  = datatree_2018->Jet_PT/1000.;
    vars_jet[1]  = datatree_2018->Jet_PE/1000.;
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
    for(int h1_index = 0 ; h1_index < datatree_2018->Jet_NDtr ; h1_index++)
    {
      // Skip non-hadronic particles
      if(datatree_2018->Jet_Dtr_IsMeson[h1_index]!=1&&datatree_2018->Jet_Dtr_IsBaryon[h1_index]!=1) continue;

      h1_4vector->SetPxPyPzE(datatree_2018->Jet_Dtr_PX[h1_index]/1000.,datatree_2018->Jet_Dtr_PY[h1_index]/1000.,datatree_2018->Jet_Dtr_PZ[h1_index]/1000.,datatree_2018->Jet_Dtr_E[h1_index]/1000.);
      if(!apply_chargedtrack_cuts(datatree_2018->Jet_Dtr_ThreeCharge[h1_index],
                                  datatree_2018->Jet_Dtr_P[h1_index]/1000.,
                                  datatree_2018->Jet_Dtr_PT[h1_index]/1000.,
                                  datatree_2018->Jet_Dtr_TrackChi2[h1_index]/datatree_2018->Jet_Dtr_TrackNDF[h1_index],
                                  datatree_2018->Jet_Dtr_ProbNNghost[h1_index],
                                  Jet_4vector->DeltaR(*h1_4vector))) continue;

      double h1_purity     = hpurity->GetBinContent(hpurity->FindBin(datatree_2018->Jet_Dtr_P[h1_index]/1000.,datatree_2018->Jet_Dtr_ETA[h1_index],datatree_2018->Jet_PT/1000.));
      double h1_efficiency = hefficiency->GetBinContent(hefficiency->FindBin(datatree_2018->Jet_Dtr_P[h1_index]/1000.,datatree_2018->Jet_Dtr_ETA[h1_index],datatree_2018->Jet_PT/1000.));
      if(h1_purity>1.||h1_efficiency>1.) continue;

      // Loop over hadron 2
      for(int h2_index = h1_index+1 ; h2_index < datatree_2018->Jet_NDtr ; h2_index++)
      {
        // Skip non-hadronic particles
        if(datatree_2018->Jet_Dtr_IsMeson[h2_index]!=1&&datatree_2018->Jet_Dtr_IsBaryon[h2_index]!=1) continue;

        h2_4vector->SetPxPyPzE(datatree_2018->Jet_Dtr_PX[h2_index]/1000.,datatree_2018->Jet_Dtr_PY[h2_index]/1000.,datatree_2018->Jet_Dtr_PZ[h2_index]/1000.,datatree_2018->Jet_Dtr_E[h2_index]/1000.);
        if(!apply_chargedtrack_cuts(datatree_2018->Jet_Dtr_ThreeCharge[h2_index],
                                    datatree_2018->Jet_Dtr_P[h2_index]/1000.,
                                    datatree_2018->Jet_Dtr_PT[h2_index]/1000.,
                                    datatree_2018->Jet_Dtr_TrackChi2[h2_index]/datatree_2018->Jet_Dtr_TrackNDF[h2_index],
                                    datatree_2018->Jet_Dtr_ProbNNghost[h2_index],
                                    Jet_4vector->DeltaR(*h2_4vector))) continue;

        double h2_purity     = hpurity->GetBinContent(hpurity->FindBin(datatree_2018->Jet_Dtr_P[h2_index]/1000.,datatree_2018->Jet_Dtr_ETA[h2_index],datatree_2018->Jet_PT/1000.));
        double h2_efficiency = hefficiency->GetBinContent(hefficiency->FindBin(datatree_2018->Jet_Dtr_P[h2_index]/1000.,datatree_2018->Jet_Dtr_ETA[h2_index],datatree_2018->Jet_PT/1000.));
        if(h2_purity>1.||h2_efficiency>1.) continue;

        double purity_correction = (h1_purity)*(h2_purity);

        double h1_purity_err = hpurity->GetBinError(hpurity->FindBin(datatree_2018->Jet_Dtr_P[h1_index]/1000.,datatree_2018->Jet_Dtr_ETA[h1_index],datatree_2018->Jet_PT/1000.));
        double h2_purity_err = hpurity->GetBinError(hpurity->FindBin(datatree_2018->Jet_Dtr_P[h2_index]/1000.,datatree_2018->Jet_Dtr_ETA[h2_index],datatree_2018->Jet_PT/1000.));
        double purity_error = sqrt(pow((h1_purity)*(h2_purity_err),2) + pow((h1_purity_err)*(h2_purity),2));

        double efficiency_correction = (h1_efficiency)*(h2_efficiency);

        double h1_efficiency_err = hefficiency->GetBinError(hefficiency->FindBin(datatree_2018->Jet_Dtr_P[h1_index]/1000.,datatree_2018->Jet_Dtr_ETA[h1_index],datatree_2018->Jet_PT/1000.));
        double h2_efficiency_err = hefficiency->GetBinError(hefficiency->FindBin(datatree_2018->Jet_Dtr_P[h2_index]/1000.,datatree_2018->Jet_Dtr_ETA[h2_index],datatree_2018->Jet_PT/1000.));
        double efficiency_error = sqrt(pow((h1_efficiency)*(h2_efficiency_err),2) + pow((h1_efficiency_err)*(h2_efficiency),2));

        vars[0 ] = weight(datatree_2018->Jet_Dtr_E[h1_index], datatree_2018->Jet_Dtr_E[h2_index], datatree_2018->Jet_PE);
        vars[1 ] = efficiency_correction;
        vars[2 ] = purity_correction;
        vars[3 ] = efficiency_error/efficiency_correction;
        vars[4 ] = purity_error/purity_correction;
        vars[5 ] = h1_4vector->DeltaR(*h2_4vector);
        vars[6 ] = datatree_2018->Jet_Dtr_ETA[h1_index];
        vars[7 ] = datatree_2018->Jet_Dtr_ETA[h2_index];
        vars[8 ] = h1_4vector->Rapidity();
        vars[9 ] = h2_4vector->Rapidity();
        vars[10] = datatree_2018->Jet_Dtr_P[h1_index]/1000.;
        vars[11] = datatree_2018->Jet_Dtr_P[h2_index]/1000.;
        vars[12] = datatree_2018->Jet_Dtr_PT[h1_index]/1000.;
        vars[13] = datatree_2018->Jet_Dtr_PT[h2_index]/1000.;
        vars[14] = datatree_2018->Jet_PT/1000.;
        vars[15] = Jet_4vector->Eta();
        vars[16] = weight(datatree_2018->Jet_Dtr_PT[h1_index], datatree_2018->Jet_Dtr_PT[h2_index], datatree_2018->Jet_PT);
        vars[17] = datatree_2018->Jet_PE/1000.;
        vars[18] = datatree_2018->Jet_Dtr_E[h1_index]/1000.;
        vars[19] = datatree_2018->Jet_Dtr_E[h2_index]/1000.;
        vars[20] = 2018;
        
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
