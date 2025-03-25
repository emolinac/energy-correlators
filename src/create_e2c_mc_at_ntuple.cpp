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
  TFile* fout = new TFile((output_folder+namef_ntuple_mc_at_e2c).c_str(),"RECREATE");
  
  // Declare the TTrees to be used to build the ntuples
  TZJetsMC* mctree = new TZJetsMC();
  TZJetsMCReco* mcrecotree = new TZJetsMCReco();
  
  // Create Ntuples
  TNtuple* ntuple_mcreco     = new TNtuple(name_ntuple_mcreco.c_str()    ,"Reco Sim" ,ntuple_mcreco_at_vars);
  TNtuple* ntuple_mcreco_jet = new TNtuple(name_ntuple_mcreco_jet.c_str(),"Reco Sim","jet_pt");
  ntuple_mcreco->SetAutoSave(0);
  ntuple_mcreco_jet->SetAutoSave(0);

  // Create Ntuples
  TNtuple* ntuple_mc     = new TNtuple(name_ntuple_mc.c_str()    ,"MC Sim",ntuple_mc_at_vars);
  TNtuple* ntuple_mc_jet = new TNtuple(name_ntuple_mc_jet.c_str(),"MC Sim","jet_pt");
  ntuple_mc->SetAutoSave(0);
  ntuple_mc_jet->SetAutoSave(0);
  
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
  
  // Fill the MC TNtuple
  float vars_mc[Nvars_mc_at];
  float vars_mc_jet[1];
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

    vars_mc_jet[0] = mctree->MCJet_PT/1000.;
    ntuple_mc_jet->Fill(vars_mc_jet);
    
    for(int h1_index = 0 ; h1_index < mctree->MCJet_AssociatedTracks_n_Tracks ; h1_index++)
    {
      // Skip non-hadronic particles
      if(mctree->MCJet_Dtr_IsMeson[h1_index]!=1&&mctree->MCJet_Dtr_IsBaryon[h1_index]!=1) continue;

      h1_4vector->SetPxPyPzE(mctree->MCJet_Dtr_PX[h1_index]/1000.,
                             mctree->MCJet_Dtr_PY[h1_index]/1000.,
                             mctree->MCJet_Dtr_PZ[h1_index]/1000., 
                             mctree->MCJet_Dtr_E[h1_index]/1000.);

      if(!apply_chargedtrack_momentum_cuts(mctree->MCJet_Dtr_ThreeCharge[h1_index],
                                           mctree->MCJet_Dtr_P[h1_index]/1000.,
                                           mctree->MCJet_Dtr_PT[h1_index]/1000.,
                                           Jet_4vector->DeltaR(*h1_4vector))) continue;

      for(int h2_index = h1_index+1 ; h2_index < mctree->MCJet_Dtr_nmcdtrs ; h2_index++)
      {
          // Skip non-hadronic particles
          if(mctree->MCJet_Dtr_IsMeson[h2_index]!=1&&mctree->MCJet_Dtr_IsBaryon[h2_index]!=1) continue;

          h2_4vector->SetPxPyPzE(mctree->MCJet_Dtr_PX[h2_index]/1000.,
                                 mctree->MCJet_Dtr_PY[h2_index]/1000.,
                                 mctree->MCJet_Dtr_PZ[h2_index]/1000., 
                                 mctree->MCJet_Dtr_E[h2_index]/1000.);

          if(!apply_chargedtrack_momentum_cuts(mctree->MCJet_Dtr_ThreeCharge[h2_index],
                                               mctree->MCJet_Dtr_P[h2_index]/1000.,
                                               mctree->MCJet_Dtr_PT[h2_index]/1000.,
                                               Jet_4vector->DeltaR(*h2_4vector))) continue;


          double h1_y = rapidity(mctree->MCJet_Dtr_E[h1_index],mctree->MCJet_Dtr_PZ[h1_index]); 
          double h2_y = rapidity(mctree->MCJet_Dtr_E[h2_index],mctree->MCJet_Dtr_PZ[h2_index]);

          vars_mc[0]  = weight(mctree->MCJet_Dtr_E[h1_index]/1000.,mctree->MCJet_Dtr_E[h2_index]/1000.,mctree->MCJet_PE/1000.);
          vars_mc[1]  = R_L(h1_y, h2_y, mctree->MCJet_Dtr_PHI[h1_index], mctree->MCJet_Dtr_PHI[h2_index]);
          vars_mc[2]  = mctree->MCJet_Dtr_ETA[h1_index];
          vars_mc[3]  = mctree->MCJet_Dtr_ETA[h2_index];
          vars_mc[4]  = h1_y;
          vars_mc[5]  = h2_y;
          vars_mc[6]  = mctree->MCJet_Dtr_ThreeCharge[h1_index];
          vars_mc[7]  = mctree->MCJet_Dtr_ThreeCharge[h2_index];
          vars_mc[8]  = mctree->MCJet_Dtr_P[h1_index]/1000.;
          vars_mc[9]  = mctree->MCJet_Dtr_P[h2_index]/1000.;
          vars_mc[10] = mctree->MCJet_Dtr_PT[h1_index]/1000.;
          vars_mc[11] = mctree->MCJet_Dtr_PT[h2_index]/1000.; 
          vars_mc[12] = mctree->MCJet_PT/1000.;
          vars_mc[13] = mctree->MCJet_ETA;
          vars_mc[14] = weight(mctree->MCJet_Dtr_PT[h1_index], mctree->MCJet_Dtr_PT[h2_index], mctree->MCJet_PT);
          vars_mc[15] = delta_phi(mctree->MCJet_PHI,Z0_4vector->Phi());
          vars_mc[16] = R_L(Jet_4vector->Rapidity(),mum_4vector->Rapidity(),mctree->MCJet_PHI,mctree->MCJet_truth_mum_PHI);
          vars_mc[17] = mctree->MCJet_truth_mum_PT/1000.;
          vars_mc[18] = mctree->MCJet_truth_mum_ETA;
          vars_mc[19] = R_L(Jet_4vector->Rapidity(),mup_4vector->Rapidity(),mctree->MCJet_PHI,mctree->MCJet_truth_mup_PHI);
          vars_mc[20] = mctree->MCJet_truth_mup_PT/1000.;
          vars_mc[21] = mctree->MCJet_truth_mup_ETA;
          vars_mc[22] = mctree->MCJet_Dtr_ID[h1_index];
          vars_mc[23] = mctree->MCJet_Dtr_ID[h2_index];
            
          // Fill the TNtuple
          ntuple_mc->Fill(vars_mc);        
      }   
    }

    last_eventNum = mctree->eventNumber;
  }

  float vars[Nvars_mcreco];
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

    ntuple_mcreco_jet->Fill(Jet_4vector->Pt());
            
    // Loop over hadron 1
    for(int h1_index = 0 ; h1_index < mcrecotree->Jet_NDtr ; h1_index++)
    {
        // Skip non-hadronic particles
        if(mcrecotree->Jet_Dtr_IsMeson[h1_index]!=1&&mcrecotree->Jet_Dtr_IsBaryon[h1_index]!=1) continue;

        h1_4vector->SetPxPyPzE(mcrecotree->Jet_Dtr_PX[h1_index]/1000.,mcrecotree->Jet_Dtr_PY[h1_index]/1000.,mcrecotree->Jet_Dtr_PZ[h1_index]/1000.,mcrecotree->Jet_Dtr_E[h1_index]/1000.);
        if(!apply_chargedtrack_cuts(mcrecotree->Jet_Dtr_ThreeCharge[h1_index],
                                    mcrecotree->Jet_Dtr_P[h1_index]/1000.,
                                    mcrecotree->Jet_Dtr_PT[h1_index]/1000.,
                                    mcrecotree->Jet_Dtr_TrackChi2[h1_index]/mcrecotree->Jet_Dtr_TrackNDF[h1_index],
                                    mcrecotree->Jet_Dtr_ProbNNghost[h1_index],
                                    Jet_4vector->DeltaR(*h1_4vector))) continue;

        // Loop over hadron 2
        for(int h2_index = h1_index+1 ; h2_index < mcrecotree->Jet_NDtr ; h2_index++)
        {
            // Skip non-hadronic particles
            if(mcrecotree->Jet_Dtr_IsMeson[h2_index]!=1&&mcrecotree->Jet_Dtr_IsBaryon[h2_index]!=1) continue;

            h2_4vector->SetPxPyPzE(mcrecotree->Jet_Dtr_PX[h2_index]/1000.,mcrecotree->Jet_Dtr_PY[h2_index]/1000.,mcrecotree->Jet_Dtr_PZ[h2_index]/1000.,mcrecotree->Jet_Dtr_E[h2_index]/1000.);
            if(!apply_chargedtrack_cuts(mcrecotree->Jet_Dtr_ThreeCharge[h2_index],
                                        mcrecotree->Jet_Dtr_P[h2_index]/1000.,
                                        mcrecotree->Jet_Dtr_PT[h2_index]/1000.,
                                        mcrecotree->Jet_Dtr_TrackChi2[h2_index]/mcrecotree->Jet_Dtr_TrackNDF[h2_index],
                                        mcrecotree->Jet_Dtr_ProbNNghost[h2_index],
                                        Jet_4vector->DeltaR(*h2_4vector))) continue;

            double h1_y = rapidity(mcrecotree->Jet_Dtr_E[h1_index],mcrecotree->Jet_Dtr_PZ[h1_index]); 
            double h2_y = rapidity(mcrecotree->Jet_Dtr_E[h2_index],mcrecotree->Jet_Dtr_PZ[h2_index]);

            // If all good, fille Ntuple
            vars[0]  = weight(mcrecotree->Jet_Dtr_E[h1_index], mcrecotree->Jet_Dtr_E[h2_index], mcrecotree->Jet_PE);
            vars[1]  = R_L(h1_y, h2_y, mcrecotree->Jet_Dtr_PHI[h1_index], mcrecotree->Jet_Dtr_PHI[h2_index]);
            vars[2]  = mcrecotree->Jet_Dtr_ETA[h1_index];
            vars[3]  = mcrecotree->Jet_Dtr_ETA[h2_index];
            vars[4]  = h1_y;
            vars[5]  = h2_y;
            vars[6]  = mcrecotree->Jet_Dtr_ThreeCharge[h1_index];
            vars[7]  = mcrecotree->Jet_Dtr_ThreeCharge[h2_index];
            vars[8]  = mcrecotree->Jet_Dtr_P[h1_index]/1000.;
            vars[9]  = mcrecotree->Jet_Dtr_P[h2_index]/1000.;
            vars[10] = mcrecotree->Jet_Dtr_PT[h1_index]/1000.;
            vars[11] = mcrecotree->Jet_Dtr_PT[h2_index]/1000.;
            vars[12] = mcrecotree->Jet_PT/1000.;
            vars[13] = Jet_4vector->Eta();
            vars[14] = weight(mcrecotree->Jet_Dtr_PT[h1_index], mcrecotree->Jet_Dtr_PT[h2_index], mcrecotree->Jet_PT);
            vars[15] = delta_phi(Jet_4vector->Phi(),Z0_4vector->Phi());
            vars[16] = R_L(Jet_4vector->Rapidity(),mum_4vector->Rapidity(),Jet_4vector->Phi(),mum_4vector->Phi());
            vars[17] = mum_4vector->Pt();
            vars[18] = mum_4vector->Eta();
            vars[19] = R_L(Jet_4vector->Rapidity(),mup_4vector->Rapidity(),Jet_4vector->Phi(),mup_4vector->Phi());
            vars[20] = mup_4vector->Pt();
            vars[21] = mup_4vector->Eta();
            vars[22] = mcrecotree->Jet_Dtr_ID[h1_index];
            vars[23] = mcrecotree->Jet_Dtr_ID[h2_index];
            
            // Fill the TNtuple
            ntuple_mcreco->Fill(vars);
        }
    }

    last_eventNum = mcrecotree->eventNumber;
  }

  fout->cd();
  ntuple_mc->Write();
  ntuple_mc_jet->Write();
  ntuple_mcreco->Write();
  ntuple_mcreco_jet->Write();
  fout->Close();

  std::cout<<std::endl;
  
  return 0;
}

