#include <iostream>
#include "TZJetsMC.h"
#include "TZJetsMC.C"
#include "TZJetsMCReco.h"
#include "TZJetsMCReco.C"
#include "TROOT.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "analysis-constants.h"
#include "analysis-binning.h"
#include "analysis-cuts.h"
#include "analysis-functions.h"
#include "directories.h"
#include "names.h"

int main()
{
  // Create output file
  TFile* fout = new TFile((output_folder + namef_ntuple_e2c_pairpurity).c_str(),"RECREATE");
  
  // Declare the TTrees to be used to build the ntuples
  TZJetsMCReco* mcrecotree = new TZJetsMCReco();

  // Create Ntuples
  TNtuple* ntuple_jet_match = new TNtuple(name_ntuple_purity.c_str(),"Reco jet&dtr matched 2 MC",ntuple_corrections_reco_vars); 
  
  ntuple_jet_match->SetAutoSave(0);
  
  // Create necessary 4vectors
  TLorentzVector* Jet_4vector   = new TLorentzVector();
  TLorentzVector* Z0_4vector    = new TLorentzVector();
  TLorentzVector* mum_4vector   = new TLorentzVector();
  TLorentzVector* mup_4vector   = new TLorentzVector();
  TLorentzVector* h1_4vector     = new TLorentzVector();
  TLorentzVector* h2_4vector     = new TLorentzVector();

  TLorentzVector* true_Jet_4vector = new TLorentzVector();
  TLorentzVector* true_Z0_4vector  = new TLorentzVector();
  TLorentzVector* true_mum_4vector = new TLorentzVector();
  TLorentzVector* true_mup_4vector = new TLorentzVector();
  TLorentzVector* true_h1_4vector  = new TLorentzVector();
  TLorentzVector* true_h2_4vector  = new TLorentzVector();
  
  int eventNum;
  unsigned long long last_eventNum = 0;
  int events = 0;
  bool maxjetpT_found = false;
  
  // Define array carrying the variables
  float vars[Nvars_corrections_mcreco];

  // Fill the matched jets Ntuple
  for (int evt = 0 ; evt < mcrecotree->fChain->GetEntries() ; evt++)
  {
    if (evt%10000 == 0)
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

    // -999 means there is not matched jet
    if (mcrecotree->Jet_mcjet_nmcdtrs == -999) continue;

    // Apply PV cut
    if (mcrecotree->nPV != 1) continue;

    // Apply trigger cut
    bool mum_trigger = (mcrecotree->mum_L0MuonEWDecision_TOS == 1 && mcrecotree->mum_Hlt1SingleMuonHighPTDecision_TOS == 1 && mcrecotree->mum_Hlt2EWSingleMuonVHighPtDecision_TOS == 1);
    bool mup_trigger = (mcrecotree->mup_L0MuonEWDecision_TOS == 1 && mcrecotree->mup_Hlt1SingleMuonHighPTDecision_TOS == 1 && mcrecotree->mup_Hlt2EWSingleMuonVHighPtDecision_TOS == 1);

    if (!mum_trigger && !mup_trigger) continue;

    // Set Jet-associated 4 vectors and apply cuts
    Jet_4vector->SetPxPyPzE(mcrecotree->Jet_PX/1000.,mcrecotree->Jet_PY/1000.,mcrecotree->Jet_PZ/1000.,mcrecotree->Jet_PE/1000.);
    if (!apply_jet_cuts(Jet_4vector->Eta(), Jet_4vector->Pt())) continue;
    
    true_Jet_4vector->SetPxPyPzE(mcrecotree->Jet_mcjet_PX/1000.,mcrecotree->Jet_mcjet_PY/1000.,mcrecotree->Jet_mcjet_PZ/1000.,mcrecotree->Jet_mcjet_PE/1000.);
    if (!apply_jet_cuts(true_Jet_4vector->Eta(),true_Jet_4vector->Pt())) continue;
    
    mum_4vector->SetPxPyPzE(mcrecotree->mum_PX/1000.,mcrecotree->mum_PY/1000.,mcrecotree->mum_PZ/1000.,mcrecotree->mum_PE/1000.);
    if (!apply_muon_cuts(Jet_4vector->DeltaR(*mum_4vector), mum_4vector->Pt(), mum_4vector->Eta())) continue;
    //if (mcrecotree->mum_TRACK_PCHI2<muon_trackprob_min) continue;

    mup_4vector->SetPxPyPzE(mcrecotree->mup_PX/1000.,mcrecotree->mup_PY/1000.,mcrecotree->mup_PZ/1000.,mcrecotree->mup_PE/1000.);
    if (!apply_muon_cuts(Jet_4vector->DeltaR(*mup_4vector), mup_4vector->Pt(), mup_4vector->Eta())) continue;
    //if (mcrecotree->mup_TRACK_PCHI2<muon_trackprob_min) continue;

    true_mum_4vector->SetPxPyPzE(mcrecotree->Jet_mcjet_mum_PX/1000.,mcrecotree->Jet_mcjet_mum_PY/1000.,mcrecotree->Jet_mcjet_mum_PZ/1000.,mcrecotree->Jet_mcjet_mum_PE/1000.);
    if (!apply_muon_cuts(true_Jet_4vector->DeltaR(*true_mum_4vector),true_mum_4vector->Pt(),true_mum_4vector->Eta())) continue;
    //if (mcrecotree->mum_TRACK_PCHI2<muon_trackprob_min) continue;

    true_mup_4vector->SetPxPyPzE(mcrecotree->Jet_mcjet_mup_PX/1000.,mcrecotree->Jet_mcjet_mup_PY/1000.,mcrecotree->Jet_mcjet_mup_PZ/1000.,mcrecotree->Jet_mcjet_mup_PE/1000.);
    if (!apply_muon_cuts(true_Jet_4vector->DeltaR(*true_mup_4vector),true_mup_4vector->Pt(),true_mup_4vector->Eta())) continue;
    
    Z0_4vector->SetPxPyPzE(mup_4vector->Px()+mum_4vector->Px(),mup_4vector->Py()+mum_4vector->Py(),mup_4vector->Pz()+mum_4vector->Pz(),mup_4vector->E() +mum_4vector->E());
    if (!apply_zboson_cuts(TMath::Abs(Jet_4vector->DeltaPhi(*Z0_4vector)), Z0_4vector->M())) continue;
    
    true_Z0_4vector->SetPxPyPzE(true_mup_4vector->Px()+true_mum_4vector->Px(),true_mup_4vector->Py()+true_mum_4vector->Py(),true_mup_4vector->Pz()+true_mum_4vector->Pz(),true_mup_4vector->E() +true_mum_4vector->E());
    if (!apply_zboson_cuts(TMath::Abs(true_Jet_4vector->DeltaPhi(*true_Z0_4vector)),true_Z0_4vector->M())) continue;
    
    // Loop over hadron 2
    for (int h1_index = 0 ; h1_index < mcrecotree->Jet_NDtr ; h1_index++)
    {
      // Skip non-hadronic particles
      if (mcrecotree->Jet_Dtr_IsMeson[h1_index] != 1 && mcrecotree->Jet_Dtr_IsBaryon[h1_index] != 1) continue;

      h1_4vector->SetPxPyPzE(mcrecotree->Jet_Dtr_PX[h1_index]/1000.,mcrecotree->Jet_Dtr_PY[h1_index]/1000.,mcrecotree->Jet_Dtr_PZ[h1_index]/1000.,mcrecotree->Jet_Dtr_E[h1_index]/1000.);
      if (!apply_chargedtrack_cuts(mcrecotree->Jet_Dtr_ThreeCharge[h1_index],
                                  h1_4vector->P(),
                                  h1_4vector->Pt(),
                                  mcrecotree->Jet_Dtr_TrackChi2[h1_index]/mcrecotree->Jet_Dtr_TrackNDF[h1_index],
                                  mcrecotree->Jet_Dtr_ProbNNghost[h1_index],
                                  h1_4vector->Eta())) continue;

      int key1_match = 0;
      if (mcrecotree->Jet_Dtr_TRUE_ETA[h1_index] != -999)
      {
        key1_match++;

       true_h1_4vector->SetPxPyPzE(mcrecotree->Jet_Dtr_TRUE_PX[h1_index]/1000.,
                                   mcrecotree->Jet_Dtr_TRUE_PY[h1_index]/1000.,
                                   mcrecotree->Jet_Dtr_TRUE_PZ[h1_index]/1000.,
                                   mcrecotree->Jet_Dtr_TRUE_E[h1_index]/1000.);
       
       if (!apply_chargedtrack_momentum_cuts(mcrecotree->Jet_Dtr_TRUE_ThreeCharge[h1_index],
                                            true_h1_4vector->P(),
                                            true_h1_4vector->Pt(),
                                            true_h1_4vector->Eta())) key1_match = 0;
      } 

      for (int h2_index = h1_index+1 ; h2_index < mcrecotree->Jet_NDtr ; h2_index++)
      {
        // Skip non-hadronic particles
        if (mcrecotree->Jet_Dtr_IsMeson[h2_index] != 1 && mcrecotree->Jet_Dtr_IsBaryon[h2_index] != 1) continue;

        h2_4vector->SetPxPyPzE(mcrecotree->Jet_Dtr_PX[h2_index]/1000.,mcrecotree->Jet_Dtr_PY[h2_index]/1000.,mcrecotree->Jet_Dtr_PZ[h2_index]/1000.,mcrecotree->Jet_Dtr_E[h2_index]/1000.);
        if (!apply_chargedtrack_cuts(mcrecotree->Jet_Dtr_ThreeCharge[h2_index],
                                    h2_4vector->P(),
                                    h2_4vector->Pt(),
                                    mcrecotree->Jet_Dtr_TrackChi2[h2_index]/mcrecotree->Jet_Dtr_TrackNDF[h2_index],
                                    mcrecotree->Jet_Dtr_ProbNNghost[h2_index],
                                    h2_4vector->Eta())) continue;

        int key2_match = 0;
        if (mcrecotree->Jet_Dtr_TRUE_ETA[h2_index] != -999)
        {
          key2_match++;

          true_h2_4vector->SetPxPyPzE(mcrecotree->Jet_Dtr_TRUE_PX[h2_index]/1000.,
                                      mcrecotree->Jet_Dtr_TRUE_PY[h2_index]/1000.,
                                      mcrecotree->Jet_Dtr_TRUE_PZ[h2_index]/1000.,
                                      mcrecotree->Jet_Dtr_TRUE_E[h2_index]/1000.);
          
          if (!apply_chargedtrack_momentum_cuts(mcrecotree->Jet_Dtr_TRUE_ThreeCharge[h2_index],
                                               true_h2_4vector->P(),
                                               true_h2_4vector->Pt(),
                                               true_h2_4vector->Eta())) key2_match = 0;
        } 
      
        vars[0]  = weight(h1_4vector->E(), h2_4vector->E(), Jet_4vector->E());
        vars[1]  = h1_4vector->DeltaR(*h2_4vector);
        vars[2]  = h1_4vector->Eta();
        vars[3]  = h2_4vector->Eta();
        vars[4]  = h1_4vector->Rapidity();
        vars[5]  = h2_4vector->Rapidity();
        vars[6]  = h1_4vector->P();
        vars[7]  = h2_4vector->P();
        vars[8]  = h1_4vector->Pt();
        vars[9]  = h2_4vector->Pt();
        vars[10] = Jet_4vector->Eta();
        vars[11] = (key1_match==0) ? -999 : true_h1_4vector->P();
        vars[12] = (key2_match==0) ? -999 : true_h2_4vector->P();
        vars[13] = (key1_match==0) ? -999 : true_h1_4vector->Pt();
        vars[14] = (key2_match==0) ? -999 : true_h2_4vector->Pt();
        vars[15] = Jet_4vector->Pt();
        vars[16] = true_Jet_4vector->Pt();
        vars[17] = (key1_match==0) ? -999 : true_h1_4vector->DeltaR(*h1_4vector);
        vars[18] = (key2_match==0) ? -999 : true_h2_4vector->DeltaR(*h2_4vector);
        vars[19] = (key1_match==0) ? -999 : true_h1_4vector->Rapidity();
        vars[20] = (key2_match==0) ? -999 : true_h2_4vector->Rapidity();
        vars[21] = (key1_match==0||key2_match==0) ? -999 : true_h1_4vector->DeltaR(*true_h2_4vector);

        double weight_truth = weight(true_h1_4vector->E(), true_h2_4vector->E(), true_Jet_4vector->E());
        double weight_pt_truth = weight(true_h1_4vector->Pt(), true_h2_4vector->Pt(), true_Jet_4vector->Pt());
        vars[22] = (key1_match==0||key2_match==0) ? -999 : weight_truth;
        vars[23] = (key1_match==0||key2_match==0) ? -999 : weight_pt_truth;
        vars[24] = weight(h1_4vector->Pt(), h2_4vector->Pt(), Jet_4vector->Pt());

        // Fill the TNtuple
        ntuple_jet_match->Fill(vars);
      }
    }

    last_eventNum = mcrecotree->eventNumber;
  }

  fout->cd();
  ntuple_jet_match->Write();
  fout->Close();

  std::cout<<std::endl;

  return 0;
}
