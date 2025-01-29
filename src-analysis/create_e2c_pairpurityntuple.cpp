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
    TFile* fout = new TFile((output_folder+namef_ntuple_e2c_pairpurity).c_str(),"RECREATE");
    
    // Declare the TTrees to be used to build the ntuples
    TZJetsMCReco* mcrecotree = new TZJetsMCReco();

    // Create Ntuples
    TNtuple* ntuple_jet_match = new TNtuple(name_ntuple_purity.c_str(),"Reco jet&dtr matched 2 MC",ntuple_pairpurity_vars); 
    
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
    float vars[Nvars_pairpurity];

    // Fill the matched jets Ntuple
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
        
        true_Jet_4vector->SetPxPyPzE(mcrecotree->Jet_mcjet_PX/1000.,mcrecotree->Jet_mcjet_PY/1000.,mcrecotree->Jet_mcjet_PZ/1000.,mcrecotree->Jet_mcjet_PE/1000.);
        if(!apply_jet_cuts(true_Jet_4vector->Eta(),true_Jet_4vector->Pt())) continue;
        
        mum_4vector->SetPxPyPzE(mcrecotree->mum_PX/1000.,mcrecotree->mum_PY/1000.,mcrecotree->mum_PZ/1000.,mcrecotree->mum_PE/1000.);
        if(!apply_muon_cuts(Jet_4vector->DeltaR(*mum_4vector),mum_4vector->Pt(),mum_4vector->Eta())) continue;
        //if(mcrecotree->mum_TRACK_PCHI2<muon_trackprob_min) continue;

        mup_4vector->SetPxPyPzE(mcrecotree->mup_PX/1000.,mcrecotree->mup_PY/1000.,mcrecotree->mup_PZ/1000.,mcrecotree->mup_PE/1000.);
        if(!apply_muon_cuts(Jet_4vector->DeltaR(*mup_4vector),mup_4vector->Pt(),mup_4vector->Eta())) continue;
        //if(mcrecotree->mup_TRACK_PCHI2<muon_trackprob_min) continue;

        true_mum_4vector->SetPxPyPzE(mcrecotree->Jet_mcjet_mum_PX/1000.,mcrecotree->Jet_mcjet_mum_PY/1000.,mcrecotree->Jet_mcjet_mum_PZ/1000.,mcrecotree->Jet_mcjet_mum_PE/1000.);
        if(!apply_muon_cuts(true_Jet_4vector->DeltaR(*true_mum_4vector),true_mum_4vector->Pt(),true_mum_4vector->Eta())) continue;
        //if(mcrecotree->mum_TRACK_PCHI2<muon_trackprob_min) continue;

        true_mup_4vector->SetPxPyPzE(mcrecotree->Jet_mcjet_mup_PX/1000.,mcrecotree->Jet_mcjet_mup_PY/1000.,mcrecotree->Jet_mcjet_mup_PZ/1000.,mcrecotree->Jet_mcjet_mup_PE/1000.);
        if(!apply_muon_cuts(true_Jet_4vector->DeltaR(*true_mup_4vector),true_mup_4vector->Pt(),true_mup_4vector->Eta())) continue;
        
        Z0_4vector->SetPxPyPzE(mup_4vector->Px()+mum_4vector->Px(),mup_4vector->Py()+mum_4vector->Py(),mup_4vector->Pz()+mum_4vector->Pz(),mup_4vector->E() +mum_4vector->E());
        if(!apply_zboson_cuts(TMath::Abs(Jet_4vector->DeltaPhi(*Z0_4vector)),Z0_4vector->M())) continue;
        
        true_Z0_4vector->SetPxPyPzE(true_mup_4vector->Px()+true_mum_4vector->Px(),true_mup_4vector->Py()+true_mum_4vector->Py(),true_mup_4vector->Pz()+true_mum_4vector->Pz(),true_mup_4vector->E() +true_mum_4vector->E());
        if(!apply_zboson_cuts(TMath::Abs(true_Jet_4vector->DeltaPhi(*true_Z0_4vector)),true_Z0_4vector->M())) continue;
        
        // Loop over hadron 2
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

            int key1_match = 0;
            if(mcrecotree->Jet_Dtr_TRUE_ETA[h1_index]!=-999)
            {
//                for(int h_subindex = 0 ; h_subindex < mcrecotree->Jet_NDtr ; h_subindex++)
//                {   
//                    if(mcrecotree->Jet_Dtr_TRUE_KEY[h1_index]==mcrecotree->Jet_mcjet_dtrKeys[h_subindex]) {key_match++; key_match1_index = h_subindex;}
//                }
                key1_match++;

                true_h1_4vector->SetPxPyPzE(mcrecotree->Jet_Dtr_TRUE_PX[h1_index]/1000.,
                                            mcrecotree->Jet_Dtr_TRUE_PY[h1_index]/1000.,
                                            mcrecotree->Jet_Dtr_TRUE_PZ[h1_index]/1000.,
                                            mcrecotree->Jet_Dtr_TRUE_E[h1_index]/1000.);
                
                if(!apply_chargedtrack_momentum_cuts(mcrecotree->Jet_Dtr_TRUE_ID[h1_index],
                                                     mcrecotree->Jet_Dtr_TRUE_P[h1_index]/1000.,
                                                     mcrecotree->Jet_Dtr_TRUE_PT[h1_index],
                                                     true_Jet_4vector->DeltaR(*true_h1_4vector))) key1_match = 0;
            } 

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

              int key2_match = 0;
              if(mcrecotree->Jet_Dtr_TRUE_ETA[h2_index]!=-999)
              {
//                 for(int h_subindex = 0 ; h_subindex < mcrecotree->Jet_NDtr ; h_subindex++)
//                 {   
//                     if(mcrecotree->Jet_Dtr_TRUE_KEY[h2_index]==mcrecotree->Jet_mcjet_dtrKeys[h_subindex]) {key_match++; key_match2_index = h_subindex;}
//                 }
                key2_match++;

                true_h2_4vector->SetPxPyPzE(mcrecotree->Jet_Dtr_TRUE_PX[h2_index]/1000.,
                                            mcrecotree->Jet_Dtr_TRUE_PY[h2_index]/1000.,
                                            mcrecotree->Jet_Dtr_TRUE_PZ[h2_index]/1000.,
                                            mcrecotree->Jet_Dtr_TRUE_E[h2_index]/1000.);
                
                if(!apply_chargedtrack_momentum_cuts(mcrecotree->Jet_Dtr_TRUE_ID[h2_index],
                                                     mcrecotree->Jet_Dtr_TRUE_P[h2_index]/1000.,
                                                     mcrecotree->Jet_Dtr_TRUE_PT[h2_index],
                                                     true_Jet_4vector->DeltaR(*true_h2_4vector))) key2_match = 0;
              } 
            
              double h1_y = rapidity(mcrecotree->Jet_Dtr_E[h1_index],mcrecotree->Jet_Dtr_PZ[h1_index]);                
              double h2_y = rapidity(mcrecotree->Jet_Dtr_E[h2_index],mcrecotree->Jet_Dtr_PZ[h2_index]);

              // If all good, fill Ntuple
              vars[0]  = weight(mcrecotree->Jet_Dtr_E[h1_index]/1000., mcrecotree->Jet_Dtr_E[h2_index]/1000., mcrecotree->Jet_PE/1000.);
              vars[1]  = R_L(h1_y, h2_y, mcrecotree->Jet_Dtr_PHI[h1_index], mcrecotree->Jet_Dtr_PHI[h2_index]);
              vars[2]  = mcrecotree->Jet_Dtr_ETA[h1_index];
              vars[3]  = mcrecotree->Jet_Dtr_ETA[h2_index];
              vars[4]  = h1_y;
              vars[5]  = h2_y;
              vars[6]  = mcrecotree->Jet_Dtr_PHI[h1_index];
              vars[7]  = mcrecotree->Jet_Dtr_PHI[h2_index];
              vars[8]  = mcrecotree->Jet_Dtr_P[h1_index]/1000.;
              vars[9]  = mcrecotree->Jet_Dtr_P[h2_index]/1000.;
              vars[10] = mcrecotree->Jet_Dtr_PT[h1_index]/1000.;
              vars[11] = mcrecotree->Jet_Dtr_PT[h2_index]/1000.;
              vars[12] = Jet_4vector->Eta();
              vars[13] = Jet_4vector->DeltaPhi(*Z0_4vector);//Jet_4vector->Phi();
              vars[14] = delta_phi(Jet_4vector->Phi(),Z0_4vector->Phi());
              vars[15] = Jet_4vector->DeltaR(*mum_4vector, 1);
              vars[16] = mum_4vector->Pt();
              vars[17] = mum_4vector->Eta();
              vars[18] = Jet_4vector->DeltaR(*mup_4vector, 1);
              vars[19] = mup_4vector->Pt();
              vars[20] = mup_4vector->Eta();
              vars[21] = mcrecotree->Jet_PT/1000.;
              vars[22] = mcrecotree->Jet_mcjet_PT/1000.;
              vars[23] = mcrecotree->Jet_mcjet_nmcdtrs;

              double matchedmc_y1   = rapidity(mcrecotree->Jet_Dtr_TRUE_E[h1_index],mcrecotree->Jet_Dtr_TRUE_PZ[h1_index]);
              double matchedmc_y2   = rapidity(mcrecotree->Jet_Dtr_TRUE_E[h2_index],mcrecotree->Jet_Dtr_TRUE_PZ[h2_index]);
              double matchedmc_phi1 = mcrecotree->Jet_Dtr_TRUE_PHI[h1_index];
              double matchedmc_phi2 = mcrecotree->Jet_Dtr_TRUE_PHI[h2_index];

              vars[24] = (key1_match==0) ? -999 : R_L(h1_y, matchedmc_y1, mcrecotree->Jet_Dtr_PHI[h1_index], matchedmc_phi1);
              vars[25] = (key2_match==0) ? -999 : R_L(h2_y, matchedmc_y2, mcrecotree->Jet_Dtr_PHI[h2_index], matchedmc_phi2);
              vars[26] = (key1_match==0) ? -999 : matchedmc_y1;
              vars[27] = (key2_match==0) ? -999 : matchedmc_y2;
              vars[28] = (key1_match==0) ? -999 : matchedmc_phi1;
              vars[29] = (key2_match==0) ? -999 : matchedmc_phi2;
              vars[30] = (key1_match==0||key2_match==0) ? -999 : R_L(matchedmc_y1,matchedmc_y2,matchedmc_phi1,matchedmc_phi2);

              double weight_truth = weight(mcrecotree->Jet_Dtr_TRUE_E[h1_index],mcrecotree->Jet_Dtr_TRUE_E[h2_index],mcrecotree->Jet_mcjet_PE);
              vars[31] = (key1_match==0||key2_match==0) ? -999 : weight_truth;

              // Fill the TNtuple
              ntuple_jet_match->Fill(vars);
            }
        }
    }

    fout->cd();
    ntuple_jet_match->Write();
    fout->Close();

    return 0;
}
