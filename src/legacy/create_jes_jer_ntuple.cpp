#include <iostream>
#include "TZJetsMC.h"
#include "TZJetsMC.C"
#include "TZJetsMCReco.h"
#include "TZJetsMCReco.C"
#include "TZJetsData.h"
#include "TZJetsData.C"
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
  TFile* fout = new TFile((output_folder+namef_ntuple_jes_jer).c_str(),"RECREATE");
  
  // Declare the TTrees to be used to build the ntuples
  TZJetsMC*     mctree     = new TZJetsMC();
  TZJetsMCReco* mcrecotree = new TZJetsMCReco();
  TZJetsData*   datatree   = new TZJetsData();

  // Create Ntuples
  TNtuple* ntuple_jes_data = new TNtuple(name_ntuple_jes_data.c_str(),"",ntuple_jes_data_vars);
  TNtuple* ntuple_jes_reco = new TNtuple(name_ntuple_jes_reco.c_str(),"",ntuple_jes_reco_vars);
  TNtuple* ntuple_jes_mc   = new TNtuple(name_ntuple_jes_mc.c_str()  ,"",ntuple_jes_mc_vars);
  ntuple_jes_data->SetAutoSave(0);
  ntuple_jes_reco->SetAutoSave(0);
  ntuple_jes_mc->SetAutoSave(0);
  
  // Create necessary 4vectors
  TLorentzVector* Jet_4vector   = new TLorentzVector();
  TLorentzVector* Z0_4vector    = new TLorentzVector();
  TLorentzVector* mum_4vector   = new TLorentzVector();
  TLorentzVector* mup_4vector   = new TLorentzVector();
  
  TLorentzVector* true_Jet_4vector   = new TLorentzVector();
  TLorentzVector* true_Z0_4vector    = new TLorentzVector();
  TLorentzVector* true_mum_4vector   = new TLorentzVector();
  TLorentzVector* true_mup_4vector   = new TLorentzVector();

  int eventNum;
  unsigned long long last_eventNum = 0;
  int events = 0;
  bool maxjetpT_found = false;
  
  // Define array carrying the variables
  float vars_jes_data[Nvars_jes];
  float vars_jes_reco[Nvars_jes_reco];
  float vars_jes_mc[Nvars_jes];
  
  // Fill the reco ntuple
  for(int evt = 0 ; evt < mcrecotree->fChain->GetEntries() ; evt++)
  {
    if(evt%10000==0)
    {
      double percentage = 100*evt/mcrecotree->fChain->GetEntries();
      std::cout<<"\r"<<percentage<<"\% jets processed."<< std::flush;
    }
    
    // Access entry of tree
    mcrecotree->GetEntry(evt);

    if (evt != 0)
    {
      if (mcrecotree->eventNumber != last_eventNum) maxjetpT_found = false;
      if (last_eventNum == mcrecotree->eventNumber) continue;
    }

    last_eventNum = mcrecotree->eventNumber;
    if (maxjetpT_found) continue;

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
    //if(mcrecotree->mum_TRACK_PCHI2<muon_trackprob_min) continue;

    mup_4vector->SetPxPyPzE(mcrecotree->mup_PX/1000.,mcrecotree->mup_PY/1000.,mcrecotree->mup_PZ/1000.,mcrecotree->mup_PE/1000.);
    if(!apply_muon_cuts(Jet_4vector->DeltaR(*mup_4vector),mup_4vector->Pt(),mup_4vector->Eta())) continue;
    //if(mcrecotree->mup_TRACK_PCHI2<muon_trackprob_min) continue;

    Z0_4vector->SetPxPyPzE(mup_4vector->Px()+mum_4vector->Px(),mup_4vector->Py()+mum_4vector->Py(),mup_4vector->Pz()+mum_4vector->Pz(),mup_4vector->E() +mum_4vector->E());
    if(!apply_zboson_cuts(TMath::Abs(Jet_4vector->DeltaPhi(*Z0_4vector)),Z0_4vector->M())) continue;

    double ncandidates = mcrecotree->totCandidates;
    double subleading_jet_pt = 0;
    double leading_jet_pt    = mcrecotree->Jet_PT/1000.;
    if(ncandidates>1)
    {
      for(int cand = evt+1 ; cand < (evt+ncandidates) ; cand++)
      {
        mcrecotree->GetEntry(cand);
        
        Jet_4vector->SetPxPyPzE(mcrecotree->Jet_PX/1000.,mcrecotree->Jet_PY/1000.,mcrecotree->Jet_PZ/1000.,mcrecotree->Jet_PE/1000.);
        if(!apply_jet_cuts(Jet_4vector->Eta(),Jet_4vector->Pt())) continue;
        if(Jet_4vector->Pt() > subleading_jet_pt) subleading_jet_pt = Jet_4vector->Pt();
      }
    }
      
    if(subleading_jet_pt>0.25*leading_jet_pt) continue;

    mcrecotree->GetEntry(evt);
    if(mcrecotree->Jet_PT/1000.!=leading_jet_pt) {std::cout<<"UPS"<<std::endl;return 0;}
    vars_jes_reco[0] = mcrecotree->Jet_PT/1000.;
    vars_jes_reco[1] = Z0_4vector->Pt();
    vars_jes_reco[2] = mcrecotree->Jet_JEC_Cor;

    // Fill the TNtuple
    ntuple_jes_reco->Fill(vars_jes_reco); 
  } 

  std::cout<<"Reco done"<<std::endl;

  // // Fill the mc tuple
  // for(int evt = 0 ; evt < mctree->fChain->GetEntries() ; evt++)
  // {
  //   if(evt%10000==0)
  //   {
  //     double percentage = 100*evt/mctree->fChain->GetEntries();
  //     std::cout<<"\r"<<percentage<<"\% jets processed."<< std::flush;
  //   }
    
  //   // Access entry of tree
  //   mctree->GetEntry(evt);

  //   if (evt != 0)
  //   {
  //     if (mctree->eventNumber != last_eventNum) maxjetpT_found = false;
  //     if (last_eventNum == mctree->eventNumber) continue;
  //   }

  //   last_eventNum = mctree->eventNumber;
  //   if (maxjetpT_found) continue;

  //   // Apply PV cut
  //   if(mctree->nPVs!=1) continue;
    
  //   // Set Jet-associated 4 vectors and apply cuts
  //   Jet_4vector->SetPxPyPzE(mctree->MCJet_PX/1000.,mctree->MCJet_PY/1000.,mctree->MCJet_PZ/1000.,mctree->MCJet_PE/1000.);
  //   if(!apply_jet_cuts(Jet_4vector->Eta(),Jet_4vector->Pt())) continue;
    
  //   mum_4vector->SetPxPyPzE(mctree->MCJet_truth_mum_PX/1000.,mctree->MCJet_truth_mum_PY/1000.,mctree->MCJet_truth_mum_PZ/1000.,mctree->MCJet_truth_mum_PE/1000.);
  //   if(!apply_muon_cuts(Jet_4vector->DeltaR(*mum_4vector),mum_4vector->Pt(),mum_4vector->Eta())) continue;
    
  //   mup_4vector->SetPxPyPzE(mctree->MCJet_truth_mup_PX/1000.,mctree->MCJet_truth_mup_PY/1000.,mctree->MCJet_truth_mup_PZ/1000.,mctree->MCJet_truth_mup_PE/1000.);
  //   if(!apply_muon_cuts(Jet_4vector->DeltaR(*mup_4vector),mup_4vector->Pt(),mup_4vector->Eta())) continue;
    
  //   Z0_4vector->SetPxPyPzE(mup_4vector->Px()+mum_4vector->Px(),mup_4vector->Py()+mum_4vector->Py(),mup_4vector->Pz()+mum_4vector->Pz(),mup_4vector->E() +mum_4vector->E());
  //   if(!apply_zboson_cuts(TMath::Abs(Jet_4vector->DeltaPhi(*Z0_4vector)),Z0_4vector->M())) continue;
    
  //   double ncandidates = mctree->totCandidates;
  //   double subleading_jet_pt = 0;
  //   double leading_jet_pt    = mctree->MCJet_PT/1000.;

  //   std::cout<<"----------------------------------------------------------------"<<std::endl;
  //   std::cout<<"Working on event "<<mctree->eventNumber<<std::endl;
  //   std::cout<<"Leading jet pt = "<<leading_jet_pt<<std::endl;
  //   if(ncandidates>1)
  //   {
  //     std::cout<<"  Checking the other candidates..."<<std::endl;
  //     for(int cand = evt+1 ; cand < (evt+ncandidates) ; cand++)
  //     {
  //       mctree->GetEntry(cand);
        
  //       Jet_4vector->SetPxPyPzE(mctree->MCJet_PX/1000.,mctree->MCJet_PY/1000.,mctree->MCJet_PZ/1000.,mctree->MCJet_PE/1000.);
  //       if(!apply_jet_cuts(Jet_4vector->Eta(),Jet_4vector->Pt())) continue;
  //       if(Jet_4vector->Pt() > subleading_jet_pt) {std::cout<<"   There is subleading"<<std::endl; subleading_jet_pt = Jet_4vector->Pt();}
  //     }
  //   }

  //   std::cout<<"Leading pt    = "<<leading_jet_pt<<std::endl;
  //   std::cout<<"Subleading pt = "<<subleading_jet_pt<<std::endl;
  //   if(subleading_jet_pt>0.25*leading_jet_pt) continue;
    
  //   std::cout<<"Saving entry"<<std::endl;

  //   mctree->GetEntry(evt);
  //   if(mctree->MCJet_PT/1000.!=leading_jet_pt) {std::cout<<"UPS"<<std::endl;return 0;}
  //   vars_jes_mc[0] = mctree->MCJet_PT/1000.;
  //   vars_jes_mc[1] = Z0_4vector->Pt();

  //   // Fill the TNtuple
  //   ntuple_jes_mc->Fill(vars_jes_mc); 
  // } 

  // std::cout<<"MC done"<<std::endl;
  
  for(int evt = 0 ; evt < datatree->fChain->GetEntries() ; evt++)
  {
    // Access entry of tree
    datatree->GetEntry(evt);
    if(evt%10000==0)
    {
      double percentage = 100.*evt/datatree->fChain->GetEntries();
      std::cout<<"\r"<<percentage<<"\% jets processed."<< std::flush;
    }
    // Access entry of tree
    datatree->GetEntry(evt);

    if (evt != 0)
    {
      if (datatree->eventNumber != last_eventNum) maxjetpT_found = false;
      if (last_eventNum == datatree->eventNumber) continue;
    }

    last_eventNum = datatree->eventNumber;
    if (maxjetpT_found) continue;

    // Apply PV cut
    if(datatree->nPV!=1) continue;

    // Apply trigger cut
    bool mum_trigger = (datatree->mum_L0MuonEWDecision_TOS==1&&datatree->mum_Hlt1SingleMuonHighPTDecision_TOS==1&&datatree->mum_Hlt2EWSingleMuonVHighPtDecision_TOS==1);
    bool mup_trigger = (datatree->mup_L0MuonEWDecision_TOS==1&&datatree->mup_Hlt1SingleMuonHighPTDecision_TOS==1&&datatree->mup_Hlt2EWSingleMuonVHighPtDecision_TOS==1);

    if(!mum_trigger&&!mup_trigger) continue;
    
    // Set Jet-associated 4 vectors and apply cuts
    Jet_4vector->SetPxPyPzE(datatree->Jet_PX/1000.,datatree->Jet_PY/1000.,datatree->Jet_PZ/1000.,datatree->Jet_PE/1000.);
    if(!apply_jet_cuts(Jet_4vector->Eta(),Jet_4vector->Pt())) continue;
    
    mum_4vector->SetPxPyPzE(datatree->mum_PX/1000.,datatree->mum_PY/1000.,datatree->mum_PZ/1000.,datatree->mum_PE/1000.);
    if(!apply_muon_cuts(Jet_4vector->DeltaR(*mum_4vector),mum_4vector->Pt(),mum_4vector->Eta())) continue;
    //if(datatree->mum_TRACK_PCHI2<muon_trackprob_min) continue;
    
    mup_4vector->SetPxPyPzE(datatree->mup_PX/1000.,datatree->mup_PY/1000.,datatree->mup_PZ/1000.,datatree->mup_PE/1000.);
    if(!apply_muon_cuts(Jet_4vector->DeltaR(*mup_4vector),mup_4vector->Pt(),mup_4vector->Eta())) continue;
    //if(datatree->mup_TRACK_PCHI2<muon_trackprob_min) continue;
    
    Z0_4vector->SetPxPyPzE(mup_4vector->Px()+mum_4vector->Px(),mup_4vector->Py()+mum_4vector->Py(),mup_4vector->Pz()+mum_4vector->Pz(),mup_4vector->E() +mum_4vector->E());
    if(!apply_zboson_cuts(TMath::Abs(Jet_4vector->DeltaPhi(*Z0_4vector)),Z0_4vector->M())) continue;

    double ncandidates = datatree->totCandidates;
    double subleading_jet_pt = 0;
    double leading_jet_pt    = datatree->Jet_PT/1000.;
    if(ncandidates>1)
    {
      for(int cand = evt+1 ; cand < (evt+ncandidates) ; cand++)
      {
        datatree->GetEntry(cand);
        
        Jet_4vector->SetPxPyPzE(datatree->Jet_PX/1000.,datatree->Jet_PY/1000.,datatree->Jet_PZ/1000.,datatree->Jet_PE/1000.);
        if(!apply_jet_cuts(Jet_4vector->Eta(),Jet_4vector->Pt())) continue;
        if(Jet_4vector->Pt() > subleading_jet_pt) subleading_jet_pt = Jet_4vector->Pt();
      }
    }
      
    if(subleading_jet_pt>0.25*leading_jet_pt) continue;
    
    datatree->GetEntry(evt);
    if(datatree->Jet_PT/1000.!=leading_jet_pt) {std::cout<<"QTF"<<std::endl;return 0;}
    vars_jes_data[0] = datatree->Jet_PT/1000.;
    vars_jes_data[1] = Z0_4vector->Pt();

    // Fill the TNtuple
    ntuple_jes_data->Fill(vars_jes_data); 
  }

  std::cout<<"Data done"<<std::endl;

  fout->cd();
  ntuple_jes_data->Write();
  ntuple_jes_reco->Write();
  ntuple_jes_mc->Write();
  fout->Close();
  
  return 0;
}
