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
#include "analysis-constants.h"
#include "analysis-cuts.h"
#include "analysis-functions.h"
#include "directories.h"
#include "names.h"

int main()
{
  // Create output file
  TFile* fout = new TFile((output_folder+namef_ntuple_hadron).c_str(),"RECREATE");
  
  // Declare the TTrees to be used to build the ntuples
  TZJetsMC* mctree              = new TZJetsMC();
  TZJetsMCReco* mcrecotree      = new TZJetsMCReco();
  TZJets2016Data* datatree_2016 = new TZJets2016Data();
  TZJets2017Data* datatree_2017 = new TZJets2017Data();
  TZJets2018Data* datatree_2018 = new TZJets2018Data();
  
  // Create Ntuples
  TNtuple* ntuple_mc     = new TNtuple(name_ntuple_mc.c_str()    ,"MC Sim"  ,"h_p:h_pt:h_eta:h_phi");
  TNtuple* ntuple_mcreco = new TNtuple(name_ntuple_mcreco.c_str(),"Reco Sim","h_p:h_pt:h_eta:h_phi");
  TNtuple* ntuple_data   = new TNtuple(name_ntuple_data.c_str()  ,"Data"    ,"h_p:h_pt:h_eta:h_phi");
  ntuple_mc->SetAutoSave(0);
  ntuple_mcreco->SetAutoSave(0);
  ntuple_data->SetAutoSave(0);
  
  TLorentzVector* Jet_4vector = new TLorentzVector();
  TLorentzVector* Z0_4vector  = new TLorentzVector();
  TLorentzVector* mum_4vector = new TLorentzVector();
  TLorentzVector* mup_4vector = new TLorentzVector();
  TLorentzVector* h_4vector  = new TLorentzVector();
  
  int eventNum;
  unsigned long long last_eventNum = 0;
  int events = 0;
  bool maxjetpT_found = false;
  
  // Fill the MC TNtuple
  float vars[4];
  for(int evt = 0 ; evt < mctree->fChain->GetEntries() ; evt++)
  {
    // Access entry of tree
    mctree->GetEntry(evt);

    if(evt%10000==0)
    {
      double percentage = 100.*evt/mctree->fChain->GetEntries();
      std::cout<<"\r"<<percentage<<"\% jets processed."<< std::flush;
    }
    
    if (evt != 0)
    {
      if (mctree->eventNumber != last_eventNum) maxjetpT_found = false;
      if (last_eventNum == mctree->eventNumber) continue;
    }
    
    // Apply PV cut
    if(mctree->nPVs!=1) continue;

    // Set Jet-associated 4 vectors and apply cuts
    Jet_4vector->SetPxPyPzE(mctree->MCJet_PX/1000.,mctree->MCJet_PY/1000.,mctree->MCJet_PZ/1000.,mctree->MCJet_PE/1000.);
    if(!apply_jet_cuts(Jet_4vector->Eta(),Jet_4vector->Pt())) continue;
    
    mum_4vector->SetPxPyPzE(mctree->MCJet_truth_mum_PX/1000.,mctree->MCJet_truth_mum_PY/1000.,mctree->MCJet_truth_mum_PZ/1000.,mctree->MCJet_truth_mum_PE/1000.);
    if(!apply_muon_cuts(Jet_4vector->DeltaR(*mum_4vector),mum_4vector->Pt(),mum_4vector->Eta())) continue;
    
    mup_4vector->SetPxPyPzE(mctree->MCJet_truth_mup_PX/1000.,mctree->MCJet_truth_mup_PY/1000.,mctree->MCJet_truth_mup_PZ/1000.,mctree->MCJet_truth_mup_PE/1000.);
    if(!apply_muon_cuts(Jet_4vector->DeltaR(*mup_4vector),mup_4vector->Pt(),mup_4vector->Eta())) continue;
    
    Z0_4vector->SetPxPyPzE(mup_4vector->Px()+mum_4vector->Px(),mup_4vector->Py()+mum_4vector->Py(),mup_4vector->Pz()+mum_4vector->Pz(),mup_4vector->E() +mum_4vector->E());
    if(!apply_zboson_cuts(TMath::Abs(Jet_4vector->DeltaPhi(*Z0_4vector)),Z0_4vector->M())) continue;

    // ntuple_mc_jet->Fill(Jet_4vector->Pt(),Jet_4vector->Eta(),Z0_4vector->Pt(),Z0_4vector->Eta(),Z0_4vector->Rapidity());
    
    for(int h1_index = 0 ; h1_index < mctree->MCJet_Dtr_nmcdtrs ; h1_index++)
    {
      // Skip non-hadronic particles
      if(mctree->MCJet_Dtr_IsMeson[h1_index]!=1&&mctree->MCJet_Dtr_IsBaryon[h1_index]!=1) continue;

      h_4vector->SetPxPyPzE(mctree->MCJet_Dtr_PX[h1_index]/1000.,
                             mctree->MCJet_Dtr_PY[h1_index]/1000.,
                             mctree->MCJet_Dtr_PZ[h1_index]/1000., 
                             mctree->MCJet_Dtr_E[h1_index]/1000.);

      if(!apply_chargedtrack_momentum_cuts(mctree->MCJet_Dtr_ThreeCharge[h1_index],
                                           h_4vector->P(),
                                           h_4vector->Pt(),
                                           Jet_4vector->DeltaR(*h_4vector))) continue;
      
      vars[0] = h_4vector->P();
      vars[1] = h_4vector->Pt();
      vars[2] = h_4vector->Eta();
      vars[3] = h_4vector->Phi();

      ntuple_mc->Fill(vars);
    }

    last_eventNum = mctree->eventNumber;
  }

  // Fill the MCReco TNtuple
  for(int evt = 0 ; evt < mcrecotree->fChain->GetEntries() ; evt++)
  {
    // Access entry of tree
    mcrecotree->GetEntry(evt);

    if(evt%10000==0)
    {
      double percentage = 100.*evt/mcrecotree->fChain->GetEntries();
      std::cout<<"\r"<<percentage<<"\% jets processed."<< std::flush;
    }
    
    if (evt != 0){if (last_eventNum == mcrecotree->eventNumber) continue;}

    // -999 means there is not matched jet
    if(mcrecotree->Jet_mcjet_nmcdtrs==-999) continue;

    // Apply PV cut
    if(mcrecotree->nPV!=1) continue;

    // Apply trigger cut
    bool mum_trigger = (mcrecotree->mum_L0MuonEWDecision_TOS==1&&mcrecotree->mum_Hlt1SingleMuonHighPTDecision_TOS==1&&mcrecotree->mum_Hlt2EWSingleMuonVHighPtDecision_TOS==1);
    bool mup_trigger = (mcrecotree->mup_L0MuonEWDecision_TOS==1&&mcrecotree->mup_Hlt1SingleMuonHighPTDecision_TOS==1&&mcrecotree->mup_Hlt2EWSingleMuonVHighPtDecision_TOS==1);

    if(!mum_trigger&&!mup_trigger) continue;

    // Set Jet-associated 4 vectors and apply cuts
    Jet_4vector->SetPxPyPzE(mcrecotree->Jet_PX/1000.,mcrecotree->Jet_PY/1000.,mcrecotree->Jet_PZ/1000.,mcrecotree->Jet_PE/1000.);
    if(!apply_jet_cuts(Jet_4vector->Eta(),Jet_4vector->Pt())) continue;
    
    mum_4vector->SetPxPyPzE(mcrecotree->mum_PX/1000.,mcrecotree->mum_PY/1000.,mcrecotree->mum_PZ/1000.,mcrecotree->mum_PE/1000.);
    if(!apply_muon_cuts(Jet_4vector->DeltaR(*mum_4vector),mum_4vector->Pt(),mum_4vector->Eta())) continue;
    
    mup_4vector->SetPxPyPzE(mcrecotree->mup_PX/1000.,mcrecotree->mup_PY/1000.,mcrecotree->mup_PZ/1000.,mcrecotree->mup_PE/1000.);
    if(!apply_muon_cuts(Jet_4vector->DeltaR(*mup_4vector),mup_4vector->Pt(),mup_4vector->Eta())) continue;
    
    Z0_4vector->SetPxPyPzE(mup_4vector->Px()+mum_4vector->Px(),mup_4vector->Py()+mum_4vector->Py(),mup_4vector->Pz()+mum_4vector->Pz(),mup_4vector->E() +mum_4vector->E());
    if(!apply_zboson_cuts(TMath::Abs(Jet_4vector->DeltaPhi(*Z0_4vector)),Z0_4vector->M())) continue;

    // ntuple_mcreco_jet->Fill(Jet_4vector->Pt(),Jet_4vector->Eta(),Z0_4vector->Pt(),Z0_4vector->Eta(),Z0_4vector->Rapidity());
            
    // Loop over hadron 1
    for(int h1_index = 0 ; h1_index < mcrecotree->Jet_NDtr ; h1_index++)
    {
      // Skip non-hadronic particles
      if(mcrecotree->Jet_Dtr_IsMeson[h1_index]!=1&&mcrecotree->Jet_Dtr_IsBaryon[h1_index]!=1) continue;

      h_4vector->SetPxPyPzE(mcrecotree->Jet_Dtr_PX[h1_index]/1000.,mcrecotree->Jet_Dtr_PY[h1_index]/1000.,mcrecotree->Jet_Dtr_PZ[h1_index]/1000.,mcrecotree->Jet_Dtr_E[h1_index]/1000.);
      if(!apply_chargedtrack_cuts(mcrecotree->Jet_Dtr_ThreeCharge[h1_index],
                                  h_4vector->P(),
                                  h_4vector->Pt(),
                                  mcrecotree->Jet_Dtr_TrackChi2[h1_index]/mcrecotree->Jet_Dtr_TrackNDF[h1_index],
                                  mcrecotree->Jet_Dtr_ProbNNghost[h1_index],
                                  Jet_4vector->DeltaR(*h_4vector))) continue;

      vars[0] = h_4vector->P();
      vars[1] = h_4vector->Pt();
      vars[2] = h_4vector->Eta();
      vars[3] = h_4vector->Phi();

      ntuple_mcreco->Fill(vars);
    }

    last_eventNum = mcrecotree->eventNumber;
  }

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

    // Loop over hadron 1
    for(int h1_index = 0 ; h1_index < datatree_2016->Jet_NDtr ; h1_index++)
    {
      // Skip non-hadronic particles
      if(datatree_2016->Jet_Dtr_IsMeson[h1_index]!=1&&datatree_2016->Jet_Dtr_IsBaryon[h1_index]!=1) continue;

      h_4vector->SetPxPyPzE(datatree_2016->Jet_Dtr_PX[h1_index]/1000.,datatree_2016->Jet_Dtr_PY[h1_index]/1000.,datatree_2016->Jet_Dtr_PZ[h1_index]/1000.,datatree_2016->Jet_Dtr_E[h1_index]/1000.);
      if(!apply_chargedtrack_cuts(datatree_2016->Jet_Dtr_ThreeCharge[h1_index],
                                  h_4vector->P(),
                                  h_4vector->Pt(),
                                  datatree_2016->Jet_Dtr_TrackChi2[h1_index]/datatree_2016->Jet_Dtr_TrackNDF[h1_index],
                                  datatree_2016->Jet_Dtr_ProbNNghost[h1_index],
                                  Jet_4vector->DeltaR(*h_4vector))) continue;

      vars[0] = h_4vector->P();
      vars[1] = h_4vector->Pt();
      vars[2] = h_4vector->Eta();
      vars[3] = h_4vector->Phi();

      ntuple_data->Fill(vars);
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

    // Loop over hadron 1
    for(int h1_index = 0 ; h1_index < datatree_2017->Jet_NDtr ; h1_index++)
    {
      // Skip non-hadronic particles
      if(datatree_2017->Jet_Dtr_IsMeson[h1_index]!=1&&datatree_2017->Jet_Dtr_IsBaryon[h1_index]!=1) continue;

      h_4vector->SetPxPyPzE(datatree_2017->Jet_Dtr_PX[h1_index]/1000.,datatree_2017->Jet_Dtr_PY[h1_index]/1000.,datatree_2017->Jet_Dtr_PZ[h1_index]/1000.,datatree_2017->Jet_Dtr_E[h1_index]/1000.);
      if(!apply_chargedtrack_cuts(datatree_2017->Jet_Dtr_ThreeCharge[h1_index],
                                  h_4vector->P(),
                                  h_4vector->Pt(),
                                  datatree_2017->Jet_Dtr_TrackChi2[h1_index]/datatree_2017->Jet_Dtr_TrackNDF[h1_index],
                                  datatree_2017->Jet_Dtr_ProbNNghost[h1_index],
                                  Jet_4vector->DeltaR(*h_4vector))) continue;

      vars[0] = h_4vector->P();
      vars[1] = h_4vector->Pt();
      vars[2] = h_4vector->Eta();
      vars[3] = h_4vector->Phi();

      ntuple_data->Fill(vars);
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

    // Loop over hadron 1
    for(int h1_index = 0 ; h1_index < datatree_2018->Jet_NDtr ; h1_index++)
    {
      // Skip non-hadronic particles
      if(datatree_2018->Jet_Dtr_IsMeson[h1_index]!=1&&datatree_2018->Jet_Dtr_IsBaryon[h1_index]!=1) continue;

      h_4vector->SetPxPyPzE(datatree_2018->Jet_Dtr_PX[h1_index]/1000.,datatree_2018->Jet_Dtr_PY[h1_index]/1000.,datatree_2018->Jet_Dtr_PZ[h1_index]/1000.,datatree_2018->Jet_Dtr_E[h1_index]/1000.);
      if(!apply_chargedtrack_cuts(datatree_2018->Jet_Dtr_ThreeCharge[h1_index],
                                  h_4vector->P(),
                                  h_4vector->Pt(),
                                  datatree_2018->Jet_Dtr_TrackChi2[h1_index]/datatree_2018->Jet_Dtr_TrackNDF[h1_index],
                                  datatree_2018->Jet_Dtr_ProbNNghost[h1_index],
                                  Jet_4vector->DeltaR(*h_4vector))) continue;

      vars[0] = h_4vector->P();
      vars[1] = h_4vector->Pt();
      vars[2] = h_4vector->Eta();
      vars[3] = h_4vector->Phi();

      ntuple_data->Fill(vars);
    }

    last_eventNum = datatree_2018->eventNumber;
  }

  fout->cd();
  ntuple_mc->Write();
  ntuple_mcreco->Write();
  ntuple_data->Write();
  fout->Close();

  std::cout<<std::endl;
  
  return 0;
}

