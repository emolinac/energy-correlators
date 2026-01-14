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
#include "TRandom3.h"
#include "analysis-constants.h"
#include "analysis-binning.h"
#include "analysis-cuts.h"
#include "directories.h"
#include "names.h"
#include "syst-jes-jer.h"
#include "utils.h"

int main()
{
        // Create output file
        TFile* fout = new TFile((output_folder + namef_ntuple_reco2truth_singlehadron_match).c_str(),"RECREATE");
        
        // Declare the TTrees to be used to build the ntuples
        TZJetsMCReco* mcrecotree = new TZJetsMCReco();

        // Create Ntuples
        TNtuple* ntuple_reco = new TNtuple(name_ntuple_correction_reco.c_str(),"",ntuple_hadroncorrections_reco_vars); 
        TNtuple* ntuple_mc   = new TNtuple(name_ntuple_correction_mc.c_str()  ,"",ntuple_hadroncorrections_mc_vars); 
        
        ntuple_reco->SetAutoSave(0);
        ntuple_mc->SetAutoSave(0);
        
        // Create necessary 4vectors
        TLorentzVector* Jet_4vector = new TLorentzVector();
        TLorentzVector* Z0_4vector  = new TLorentzVector();
        TLorentzVector* mum_4vector = new TLorentzVector();
        TLorentzVector* mup_4vector = new TLorentzVector();
        TLorentzVector* h_4vector   = new TLorentzVector();

        TLorentzVector* true_Jet_4vector = new TLorentzVector();
        TLorentzVector* true_Z0_4vector  = new TLorentzVector();
        TLorentzVector* true_mum_4vector = new TLorentzVector();
        TLorentzVector* true_mup_4vector = new TLorentzVector();
        TLorentzVector* true_h_4vector   = new TLorentzVector();
        
        int eventNum;
        unsigned long long last_eventNum = 0;
        int events = 0;
        bool maxjetpT_found = false;
        
        // Define array carrying the variables
        float vars_reco[Nvars_hadroncorrections_reco];
        float vars_mc[Nvars_hadroncorrections_mc];

        // Fill the matched jets Ntuple
        for (int evt = 0 ; evt < mcrecotree->fChain->GetEntries() ; evt++) {
                if (evt%10000 == 0) {
                        double percentage = 100*evt/mcrecotree->fChain->GetEntries();
                        std::cout<<"\r"<<percentage<<"\% jets processed."<< std::flush;
                }

                // Access entry of tree
                mcrecotree->GetEntry(evt);

                if (evt != 0)
                        if (last_eventNum == mcrecotree->eventNumber) 
                                continue;

                // -There must be a matched truth-level jet
                if (mcrecotree->Jet_mcjet_nmcdtrs == -999) 
                        continue;

                // Apply PV cut
                if (mcrecotree->nPV != 1) 
                        continue;

                // Apply trigger cut
                bool mum_trigger = (mcrecotree->mum_L0MuonEWDecision_TOS == 1 && 
                                    mcrecotree->mum_Hlt1SingleMuonHighPTDecision_TOS == 1 && 
                                    mcrecotree->mum_Hlt2EWSingleMuonVHighPtDecision_TOS == 1);

                bool mup_trigger = (mcrecotree->mup_L0MuonEWDecision_TOS == 1 && 
                                    mcrecotree->mup_Hlt1SingleMuonHighPTDecision_TOS == 1 && 
                                    mcrecotree->mup_Hlt2EWSingleMuonVHighPtDecision_TOS == 1);                                    

                if (!mum_trigger && !mup_trigger) 
                        continue;

                // Apply cuts to MCReco jets
                Jet_4vector->SetPxPyPzE(mcrecotree->Jet_PX/1000.,
                                        mcrecotree->Jet_PY/1000.,
                                        mcrecotree->Jet_PZ/1000.,
                                        mcrecotree->Jet_PE/1000.);
                
                mum_4vector->SetPxPyPzE(mcrecotree->mum_PX/1000.,
                                        mcrecotree->mum_PY/1000.,
                                        mcrecotree->mum_PZ/1000.,
                                        mcrecotree->mum_PE/1000.);
                
                mup_4vector->SetPxPyPzE(mcrecotree->mup_PX/1000.,
                                        mcrecotree->mup_PY/1000.,
                                        mcrecotree->mup_PZ/1000.,
                                        mcrecotree->mup_PE/1000.);
                
                Z0_4vector->SetPxPyPzE(mup_4vector->Px()+mum_4vector->Px(),
                                       mup_4vector->Py()+mum_4vector->Py(),
                                       mup_4vector->Pz()+mum_4vector->Pz(),
                                       mup_4vector->E() +mum_4vector->E());

                if (!apply_jet_cuts(Jet_4vector->Rapidity(), Jet_4vector->Pt())) 
                        continue;
                
                if (!apply_muon_cuts(Jet_4vector->DeltaR(*mum_4vector, true), mum_4vector->Pt(), mum_4vector->Eta())) 
                        continue;
                
                if (!apply_muon_cuts(Jet_4vector->DeltaR(*mup_4vector, true), mup_4vector->Pt(), mup_4vector->Eta())) 
                        continue;
                
                if (!apply_zboson_cuts(TMath::Abs(Jet_4vector->DeltaPhi(*Z0_4vector)), Z0_4vector->M())) 
                        continue;
                
                // Apply cuts to MC jets
                true_Jet_4vector->SetPxPyPzE(mcrecotree->Jet_mcjet_PX/1000.,
                                             mcrecotree->Jet_mcjet_PY/1000.,
                                             mcrecotree->Jet_mcjet_PZ/1000.,
                                             mcrecotree->Jet_mcjet_PE/1000.);

                if (!apply_jet_cuts(true_Jet_4vector->Rapidity(),true_Jet_4vector->Pt())) 
                        continue;
                
                // Loop over reco
                for (int h_index = 0 ; h_index < mcrecotree->Jet_NDtr ; h_index++) {
                        if (abs(mcrecotree->Jet_Dtr_ID[h_index]) < 100) 
                                        continue;

                        h_4vector->SetPxPyPzE(mcrecotree->Jet_Dtr_PX[h_index]/1000.,
                                              mcrecotree->Jet_Dtr_PY[h_index]/1000.,
                                              mcrecotree->Jet_Dtr_PZ[h_index]/1000.,
                                              mcrecotree->Jet_Dtr_E[h_index]/1000.);

                        if (!apply_chargedtrack_cuts(mcrecotree->Jet_Dtr_ThreeCharge[h_index],
                                                     h_4vector->P(),
                                                     h_4vector->Pt(),
                                                     mcrecotree->Jet_Dtr_TrackChi2[h_index]/mcrecotree->Jet_Dtr_TrackNDF[h_index],
                                                     mcrecotree->Jet_Dtr_ProbNNghost[h_index],
                                                     h_4vector->Eta())) 
                                continue;

                        int key_match = 0;
                        int key1_sim  = 0;
                        if (mcrecotree->Jet_Dtr_TRUE_ETA[h_index] != -999 && std::abs(mcrecotree->Jet_Dtr_TRUE_ID[h_index]) > 100) {
                                key_match++;

                                for(int i = 0 ; i < mcrecotree->Jet_mcjet_nmcdtrs ; i++) {
                                        if(mcrecotree->Jet_Dtr_TRUE_KEY[h_index] == mcrecotree->Jet_mcjet_dtrKeys[i]) {
                                                key1_sim++;
                                                break;
                                        }
                                }

                                // Verifies if the matched particle is in the matched jet
                                if (key1_sim == 0) 
                                        key_match = 0;

                                true_h_4vector->SetPxPyPzE(mcrecotree->Jet_Dtr_TRUE_PX[h_index]/1000.,
                                                           mcrecotree->Jet_Dtr_TRUE_PY[h_index]/1000.,
                                                           mcrecotree->Jet_Dtr_TRUE_PZ[h_index]/1000.,
                                                           mcrecotree->Jet_Dtr_TRUE_E[h_index]/1000.);
                                
                                if (!apply_chargedtrack_momentum_cuts(mcrecotree->Jet_Dtr_TRUE_ThreeCharge[h_index],
                                                                      true_h_4vector->P(),
                                                                      true_h_4vector->Pt(),
                                                                      true_h_4vector->Eta())) 
                                        key_match = 0;
                        } 

                        // If all good, fill Ntuple
                        vars_reco[0]  = h_4vector->Eta();
                        vars_reco[1]  = h_4vector->Rapidity();
                        vars_reco[2]  = h_4vector->P();
                        vars_reco[3]  = h_4vector->Pt();
                        vars_reco[4]  = Jet_4vector->Pt();
                        vars_reco[5]  = Jet_4vector->Eta();
                        vars_reco[6]  = Jet_4vector->E();
                        vars_reco[7]  = true_Jet_4vector->E();
                        vars_reco[8]  = true_Jet_4vector->Pt();
                        vars_reco[9]  = (key_match == 0) ? -999 : true_h_4vector->Rapidity();
                        vars_reco[10] = (key_match == 0) ? -999 : true_h_4vector->Eta();
                        vars_reco[11] = (key_match == 0) ? -999 : true_h_4vector->P();
                        vars_reco[12] = Jet_4vector->DeltaR(*h_4vector, true);            
                        vars_reco[13] = key_match;
                        vars_reco[14] = std::sqrt(std::pow(h_4vector->Eta(),2) + std::pow(h_4vector->Phi(),2));
                        vars_reco[15] = (key_match == 0) ? -999 : std::sqrt(std::pow(true_h_4vector->Eta(),2) + std::pow(true_h_4vector->Phi(),2));


                        // Fill the TNtuple
                        ntuple_reco->Fill(vars_reco);
                }

                // Loop over mc
                for (int h_index = 0 ; h_index < mcrecotree->Jet_mcjet_nmcdtrs ; h_index++) {
                        if (abs(mcrecotree->Jet_Dtr_ID[h_index]) < 100) 
                                        continue;
                        
                        true_h_4vector->SetPxPyPzE(mcrecotree->Jet_mcjet_dtrPX[h_index]/1000.,
                                                   mcrecotree->Jet_mcjet_dtrPY[h_index]/1000.,
                                                   mcrecotree->Jet_mcjet_dtrPZ[h_index]/1000., 
                                                   mcrecotree->Jet_mcjet_dtrE[h_index]/1000.);

                        if (!apply_chargedtrack_momentum_cuts(mcrecotree->Jet_mcjet_dtrThreeCharge[h_index],
                                                              true_h_4vector->P(),
                                                              true_h_4vector->Pt(),
                                                              true_h_4vector->Eta())) 
                                continue;

                        // If all good, fill Ntuple
                        vars_mc[0] = true_h_4vector->Eta();
                        vars_mc[1] = true_h_4vector->Rapidity();
                        vars_mc[2] = true_h_4vector->P();
                        vars_mc[3] = true_h_4vector->Pt();
                        vars_mc[4] = true_Jet_4vector->Pt();
                        vars_mc[5] = true_Jet_4vector->Eta();
                        vars_mc[6] = true_Jet_4vector->E();
                        vars_mc[7] = true_Jet_4vector->DeltaR(*true_h_4vector);
                        
                        // Fill the TNtuple
                        ntuple_mc->Fill(vars_mc);
                }

                last_eventNum = mcrecotree->eventNumber;
        }

        fout->cd();
        ntuple_reco->Write();
        ntuple_mc->Write();
        fout->Close();

        std::cout<<std::endl;

        return 0;
}
