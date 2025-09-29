#include <iostream>
#include "TZJetsMCCTCorr.h"
#include "TZJetsMCCTCorr.C"
#include "TZJetsMCRecoCTCorr.h"
#include "TZJetsMCRecoCTCorr.C"
#include "TROOT.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "analysis-constants.h"
#include "analysis-binning.h"
#include "analysis-cuts.h"
#include "analysis-functions.h"
#include "directories.h"
#include "names.h"

int main()
{
        // Create output file
        TFile* fout = new TFile((output_folder + namef_ntuple_eec_paircorrections_nojetmatch_ct).c_str(),"RECREATE");
        
        // Declare the TTrees to be used to build the ntuples
        TZJetsMCRecoCTCorr* mcrecotree = new TZJetsMCRecoCTCorr();
        TZJetsMCCTCorr* mctree         = new TZJetsMCCTCorr();

        // Create Ntuples
        TNtuple* ntuple_reco = new TNtuple(name_ntuple_correction_reco.c_str(),"",ntuple_corrections_reco_vars); 
        TNtuple* ntuple_mc   = new TNtuple(name_ntuple_correction_mc.c_str()  ,"",ntuple_corrections_mc_vars); 
        
        ntuple_reco->SetAutoSave(0);
        ntuple_mc->SetAutoSave(0);

        TNtuple* ntuple_reco_jet = new TNtuple(name_ntuple_jetpurity.c_str()    ,"",ntuple_jetpurity_vars); 
        TNtuple* ntuple_mc_jet   = new TNtuple(name_ntuple_jetefficiency.c_str(),"",ntuple_jetefficiency_vars); 
        
        ntuple_mc_jet->SetAutoSave(0);
        ntuple_reco_jet->SetAutoSave(0);
        
        TRandom3* rndm = new TRandom3(0);

        // Create necessary 4vectors
        TLorentzVector* Jet_4vector   = new TLorentzVector();
        TLorentzVector* Z0_4vector    = new TLorentzVector();
        TLorentzVector* mum_4vector   = new TLorentzVector();
        TLorentzVector* mup_4vector   = new TLorentzVector();
        TLorentzVector* h1_4vector    = new TLorentzVector();
        TLorentzVector* h2_4vector    = new TLorentzVector();
        TLorentzVector* h_4vector     = new TLorentzVector();

        TLorentzVector* true_Jet_4vector = new TLorentzVector();
        TLorentzVector* true_Z0_4vector  = new TLorentzVector();
        TLorentzVector* true_mum_4vector = new TLorentzVector();
        TLorentzVector* true_mup_4vector = new TLorentzVector();
        TLorentzVector* true_h1_4vector  = new TLorentzVector();
        TLorentzVector* true_h2_4vector  = new TLorentzVector();
        
        unsigned long long last_eventNum = 0;
        bool maxjetpT_found = false;
        
        // Define array carrying the variables
        float vars_reco[Nvars_corrections_mcreco];
        float vars_mc[Nvars_corrections_mc];

        float vars_reco_jet[Nvars_jetpurity];
        float vars_mc_jet[Nvars_jetefficiency];

        std::cout<<"Processing events ..."<<std::endl;
        for (int evt = 0 ; evt < mcrecotree->fChain->GetEntries() ; evt++) {
                if (evt%10000 == 0) {
                        double percentage = 100*evt/mcrecotree->fChain->GetEntries();
                        std::cout<<"\r"<<percentage<<"\%"<< std::flush;
                }

                // Access entry of tree
                mcrecotree->GetEntry(evt);

                if (evt != 0)
                        if (last_eventNum == mcrecotree->eventNumber) 
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

                // Set Jet-associated 4 vectors and apply cuts
                Jet_4vector->SetPxPyPzE(mcrecotree->Jet_PX/1000.,
                                        mcrecotree->Jet_PY/1000.,
                                        mcrecotree->Jet_PZ/1000.,
                                        mcrecotree->Jet_PE/1000.);
                
                if (!apply_jet_cuts(Jet_4vector->Eta(), Jet_4vector->Pt())) 
                        continue;
                
                mum_4vector->SetPxPyPzE(mcrecotree->mum_PX/1000.,
                                        mcrecotree->mum_PY/1000.,
                                        mcrecotree->mum_PZ/1000.,
                                        mcrecotree->mum_PE/1000.);
                
                if (!apply_muon_cuts(Jet_4vector->DeltaR(*mum_4vector), mum_4vector->Pt(), mum_4vector->Eta())) 
                        continue;
                
                mup_4vector->SetPxPyPzE(mcrecotree->mup_PX/1000.,
                                        mcrecotree->mup_PY/1000.,
                                        mcrecotree->mup_PZ/1000.,
                                        mcrecotree->mup_PE/1000.);

                if (!apply_muon_cuts(Jet_4vector->DeltaR(*mup_4vector), mup_4vector->Pt(), mup_4vector->Eta())) 
                        continue;
                
                Z0_4vector->SetPxPyPzE(mup_4vector->Px()+mum_4vector->Px(),
                                       mup_4vector->Py()+mum_4vector->Py(),
                                       mup_4vector->Pz()+mum_4vector->Pz(),
                                       mup_4vector->E() +mum_4vector->E());

                if (!apply_zboson_cuts(TMath::Abs(Jet_4vector->DeltaPhi(*Z0_4vector)), Z0_4vector->M())) 
                        continue;

                bool truth_passed = false;
                
                if (mcrecotree->Jet_mcjet_nmcdtrs!=-999) {
                        true_Jet_4vector->SetPxPyPzE(mcrecotree->Jet_mcjet_PX/1000.,
                                                     mcrecotree->Jet_mcjet_PY/1000.,
                                                     mcrecotree->Jet_mcjet_PZ/1000.,
                                                     mcrecotree->Jet_mcjet_PE/1000.);

                        true_mum_4vector->SetPxPyPzE(mcrecotree->mum_TRUEP_X/1000.,
                                                     mcrecotree->mum_TRUEP_Y/1000.,
                                                     mcrecotree->mum_TRUEP_Z/1000.,
                                                     mcrecotree->mum_TRUEP_E/1000.);

                        true_mup_4vector->SetPxPyPzE(mcrecotree->mup_TRUEP_X/1000.,
                                                     mcrecotree->mup_TRUEP_Y/1000.,
                                                     mcrecotree->mup_TRUEP_Z/1000.,
                                                     mcrecotree->mup_TRUEP_E/1000.);

                        true_Z0_4vector->SetPxPyPzE(true_mup_4vector->Px()+true_mum_4vector->Px(),
                                                    true_mup_4vector->Py()+true_mum_4vector->Py(),
                                                    true_mup_4vector->Pz()+true_mum_4vector->Pz(),
                                                    true_mup_4vector->E() +true_mum_4vector->E());

                        // if (apply_jet_cuts(true_Jet_4vector->Eta(),true_Jet_4vector->Pt())&&\
                        //     apply_muon_cuts(true_Jet_4vector->DeltaR(*true_mum_4vector),true_mum_4vector->Pt(),true_mum_4vector->Eta())&&\
                        //     apply_muon_cuts(true_Jet_4vector->DeltaR(*true_mup_4vector),true_mup_4vector->Pt(),true_mup_4vector->Eta())&&\
                        //     apply_zboson_cuts(TMath::Abs(true_Jet_4vector->DeltaPhi(*true_Z0_4vector)),true_Z0_4vector->M())) 
                        //         truth_passed = true;

                        if (apply_jet_cuts(true_Jet_4vector->Eta(),true_Jet_4vector->Pt())) 
                                truth_passed = true;
                }

                
                for (int h1_index = 0 ; h1_index < mcrecotree->Jet_NDtr ; h1_index++) {
                        if (mcrecotree->Jet_Dtr_IsMeson[h1_index] != 1 && mcrecotree->Jet_Dtr_IsBaryon[h1_index] != 1) 
                                continue;

                        h1_4vector->SetPxPyPzE(mcrecotree->Jet_Dtr_PX[h1_index]/1000.,
                                               mcrecotree->Jet_Dtr_PY[h1_index]/1000.,
                                               mcrecotree->Jet_Dtr_PZ[h1_index]/1000.,
                                               mcrecotree->Jet_Dtr_E[h1_index]/1000.);

                        if (!apply_chargedtrack_cuts(mcrecotree->Jet_Dtr_ThreeCharge[h1_index],
                                                     h1_4vector->P(),
                                                     h1_4vector->Pt(),
                                                     mcrecotree->Jet_Dtr_TrackChi2[h1_index]/mcrecotree->Jet_Dtr_TrackNDF[h1_index],
                                                     mcrecotree->Jet_Dtr_ProbNNghost[h1_index],
                                                     h1_4vector->Eta())) 
                                continue;

                        int key1_match = 0;
                        if (mcrecotree->Jet_Dtr_TRUE_ETA[h1_index] != -999) {
                                key1_match++;

                                true_h1_4vector->SetPxPyPzE(mcrecotree->Jet_Dtr_TRUE_PX[h1_index]/1000.,
                                                            mcrecotree->Jet_Dtr_TRUE_PY[h1_index]/1000.,
                                                            mcrecotree->Jet_Dtr_TRUE_PZ[h1_index]/1000.,
                                                            mcrecotree->Jet_Dtr_TRUE_E[h1_index]/1000.);
                                
                                if (!apply_chargedtrack_momentum_cuts(mcrecotree->Jet_Dtr_TRUE_ThreeCharge[h1_index],
                                                                      true_h1_4vector->P(),
                                                                      true_h1_4vector->Pt(),
                                                                      true_h1_4vector->Eta())) 
                                        key1_match = 0;
                        } 

                        for (int h2_index = h1_index+1 ; h2_index < mcrecotree->Jet_NDtr ; h2_index++) {
                                if (mcrecotree->Jet_Dtr_IsMeson[h2_index] != 1 && mcrecotree->Jet_Dtr_IsBaryon[h2_index] != 1) 
                                        continue;

                                h2_4vector->SetPxPyPzE(mcrecotree->Jet_Dtr_PX[h2_index]/1000.,
                                                       mcrecotree->Jet_Dtr_PY[h2_index]/1000.,
                                                       mcrecotree->Jet_Dtr_PZ[h2_index]/1000.,
                                                       mcrecotree->Jet_Dtr_E[h2_index]/1000.);

                                if (!apply_chargedtrack_cuts(mcrecotree->Jet_Dtr_ThreeCharge[h2_index],
                                                             h2_4vector->P(),
                                                             h2_4vector->Pt(),
                                                             mcrecotree->Jet_Dtr_TrackChi2[h2_index]/mcrecotree->Jet_Dtr_TrackNDF[h2_index],
                                                             mcrecotree->Jet_Dtr_ProbNNghost[h2_index],
                                                             h2_4vector->Eta())) 
                                        continue;

                                int key2_match = 0;
                                if (mcrecotree->Jet_Dtr_TRUE_ETA[h2_index] != -999) {
                                        key2_match++;

                                        true_h2_4vector->SetPxPyPzE(mcrecotree->Jet_Dtr_TRUE_PX[h2_index]/1000.,
                                                                    mcrecotree->Jet_Dtr_TRUE_PY[h2_index]/1000.,
                                                                    mcrecotree->Jet_Dtr_TRUE_PZ[h2_index]/1000.,
                                                                    mcrecotree->Jet_Dtr_TRUE_E[h2_index]/1000.);
                                        
                                        if (!apply_chargedtrack_momentum_cuts(mcrecotree->Jet_Dtr_TRUE_ThreeCharge[h2_index],
                                                                              true_h2_4vector->P(),
                                                                              true_h2_4vector->Pt(),
                                                                              true_h2_4vector->Eta())) 
                                                key2_match = 0;
                                } 
                        
                                vars_reco[0]  = mcrecotree->Jet_Dtr_ThreeCharge[h1_index]*mcrecotree->Jet_Dtr_ThreeCharge[h2_index];
                                vars_reco[1]  = h1_4vector->DeltaR(*h2_4vector);
                                vars_reco[2]  = h1_4vector->Eta();
                                vars_reco[3]  = h2_4vector->Eta();
                                vars_reco[4]  = h1_4vector->Rapidity();
                                vars_reco[5]  = h2_4vector->Rapidity();
                                vars_reco[6]  = h1_4vector->P();
                                vars_reco[7]  = h2_4vector->P();
                                vars_reco[8]  = h1_4vector->Pt();
                                vars_reco[9]  = h2_4vector->Pt();
                                vars_reco[10] = Jet_4vector->Eta();
                                vars_reco[11] = (key1_match==0) ? -999 : true_h1_4vector->P();
                                vars_reco[12] = (key2_match==0) ? -999 : true_h2_4vector->P();
                                vars_reco[13] = (key1_match==0) ? -999 : true_h1_4vector->Pt();
                                vars_reco[14] = (key2_match==0) ? -999 : true_h2_4vector->Pt();
                                vars_reco[15] = Jet_4vector->Pt();
                                vars_reco[16] = true_Jet_4vector->Pt();
                                vars_reco[17] = (key1_match==0) ? -999 : true_h1_4vector->DeltaR(*h1_4vector);
                                vars_reco[18] = (key2_match==0) ? -999 : true_h2_4vector->DeltaR(*h2_4vector);
                                vars_reco[19] = (key1_match==0) ? -999 : true_h1_4vector->Eta();
                                vars_reco[20] = (key2_match==0) ? -999 : true_h2_4vector->Eta();
                                vars_reco[21] = (key1_match==0||key2_match==0) ? -999 : true_h1_4vector->DeltaR(*true_h2_4vector);

                                double weight_truth = weight(true_h1_4vector->E(), true_h2_4vector->E(), true_Jet_4vector->E());
                                double weight_pt_truth = weight(true_h1_4vector->Pt(), true_h2_4vector->Pt(), true_Jet_4vector->Pt());
                                vars_reco[22] = (key1_match==0||key2_match==0) ? -999 : weight_truth;
                                vars_reco[23] = (key1_match==0||key2_match==0) ? -999 : weight_pt_truth;
                                vars_reco[24] = weight(h1_4vector->Pt(), h2_4vector->Pt(), Jet_4vector->Pt());

                                // Fill the TNtuple
                                ntuple_reco->Fill(vars_reco);
                        }
                }

                vars_reco_jet[0] = Jet_4vector->Pt();
                vars_reco_jet[1] = Jet_4vector->E();
                vars_reco_jet[2] = mcrecotree->Jet_NDtr;
                vars_reco_jet[3] = (truth_passed) ? true_Jet_4vector->Pt() : -999 ;
                vars_reco_jet[4] = (truth_passed) ? true_Jet_4vector->E()  : -999 ;
                vars_reco_jet[5] = (truth_passed) ? mcrecotree->Jet_mcjet_nmcdtrs : -999 ;
                vars_reco_jet[6] = (truth_passed) ? Jet_4vector->DeltaR(*true_Jet_4vector) : -999;
                vars_reco_jet[7] = Jet_4vector->Eta();

                ntuple_reco_jet->Fill(vars_reco_jet); 

                last_eventNum = mcrecotree->eventNumber;
        }

        for (int evt = 0 ; evt < mctree->fChain->GetEntries() ; evt++)
        {
                mctree->GetEntry(evt);

                if (evt%10000 == 0) {
                        double percentage = 100*evt/mctree->fChain->GetEntries();
                        std::cout<<"\r"<<percentage<<"\% jets processed."<< std::flush;
                }

                if (evt != 0)
                        if (last_eventNum == mctree->eventNumber) 
                                continue;
                
                if (mctree->nPVs!=1) 
                        continue;

                true_Jet_4vector->SetPxPyPzE(mctree->MCJet_PX/1000.,
                                             mctree->MCJet_PY/1000.,
                                             mctree->MCJet_PZ/1000.,
                                             mctree->MCJet_PE/1000.);

                if (!apply_jet_cuts(true_Jet_4vector->Eta(),true_Jet_4vector->Pt())) 
                        continue;

                true_mum_4vector->SetPxPyPzE(mctree->MCJet_truth_mum_PX/1000.,
                                             mctree->MCJet_truth_mum_PY/1000.,
                                             mctree->MCJet_truth_mum_PZ/1000.,
                                             mctree->MCJet_truth_mum_PE/1000.);

                if (!apply_muon_cuts(true_Jet_4vector->DeltaR(*true_mum_4vector),true_mum_4vector->Pt(),true_mum_4vector->Eta())) 
                        continue;
                
                true_mup_4vector->SetPxPyPzE(mctree->MCJet_truth_mup_PX/1000.,
                                             mctree->MCJet_truth_mup_PY/1000.,
                                             mctree->MCJet_truth_mup_PZ/1000.,
                                             mctree->MCJet_truth_mup_PE/1000.);

                if (!apply_muon_cuts(true_Jet_4vector->DeltaR(*true_mup_4vector),true_mup_4vector->Pt(),true_mup_4vector->Eta())) 
                        continue;
                
                true_Z0_4vector->SetPxPyPzE(true_mup_4vector->Px()+true_mum_4vector->Px(),
                                            true_mup_4vector->Py()+true_mum_4vector->Py(),
                                            true_mup_4vector->Pz()+true_mum_4vector->Pz(),
                                            true_mup_4vector->E() +true_mum_4vector->E());

                if (!apply_zboson_cuts(TMath::Abs(true_Jet_4vector->DeltaPhi(*true_Z0_4vector)),true_Z0_4vector->M())) 
                        continue;

                bool reco_passed = false;
                
                if (mctree->MCJet_recojet_PX!=-999) {
                        double mum_energy = sqrt(pow(mctree->MCJet_truth_match_mum_PX,2) + 
                                                 pow(mctree->MCJet_truth_match_mum_PY,2) +
                                                 pow(mctree->MCJet_truth_match_mum_PZ,2));

                        double mup_energy = sqrt(pow(mctree->MCJet_truth_match_mup_PX,2) + 
                                                 pow(mctree->MCJet_truth_match_mup_PY,2) +
                                                 pow(mctree->MCJet_truth_match_mup_PZ,2));
                        
                        Jet_4vector->SetPxPyPzE(mctree->MCJet_recojet_PX/1000.,
                                                mctree->MCJet_recojet_PY/1000.,
                                                mctree->MCJet_recojet_PZ/1000.,
                                                mctree->MCJet_recojet_PE/1000.);

                        mum_4vector->SetPxPyPzE(mctree->MCJet_truth_match_mum_PX/1000.,
                                                mctree->MCJet_truth_match_mum_PY/1000.,
                                                mctree->MCJet_truth_match_mum_PZ/1000.,
                                                mum_energy/1000.);

                        mup_4vector->SetPxPyPzE(mctree->MCJet_truth_match_mup_PX/1000.,
                                                mctree->MCJet_truth_match_mup_PY/1000.,
                                                mctree->MCJet_truth_match_mup_PZ/1000.,
                                                mup_energy/1000.);

                        Z0_4vector->SetPxPyPzE(mup_4vector->Px()+mum_4vector->Px(),
                                               mup_4vector->Py()+mum_4vector->Py(),
                                               mup_4vector->Pz()+mum_4vector->Pz(),
                                               mup_4vector->E() +mum_4vector->E());

                        // if (apply_jet_cuts(Jet_4vector->Eta(), Jet_4vector->Pt())&&\
                        //     apply_muon_cuts(Jet_4vector->DeltaR(*mum_4vector), mum_4vector->Pt(), mum_4vector->Eta())&&\
                        //     apply_muon_cuts(Jet_4vector->DeltaR(*mup_4vector), mup_4vector->Pt(), mup_4vector->Eta())&&\
                        //     apply_zboson_cuts(TMath::Abs(Jet_4vector->DeltaPhi(*Z0_4vector)), Z0_4vector->M())) 
                        //         reco_passed = true;

                        if (apply_jet_cuts(Jet_4vector->Eta(), Jet_4vector->Pt())) 
                                reco_passed = true;
                }
                

                for (int h1_index = 0 ; h1_index < mctree->MCJet_Dtr_nmcdtrs ; h1_index++) {
                        if (mctree->MCJet_Dtr_IsMeson[h1_index] != 1 && mctree->MCJet_Dtr_IsBaryon[h1_index] != 1) 
                                continue;

                        h1_4vector->SetPxPyPzE(mctree->MCJet_Dtr_PX[h1_index]/1000.,
                                               mctree->MCJet_Dtr_PY[h1_index]/1000.,
                                               mctree->MCJet_Dtr_PZ[h1_index]/1000.,
                                               mctree->MCJet_Dtr_E[h1_index]/1000.);
                        
                        if (!apply_chargedtrack_momentum_cuts(mctree->MCJet_Dtr_ThreeCharge[h1_index],
                                                              h1_4vector->P(),
                                                              h1_4vector->Pt(),
                                                              h1_4vector->Eta())) 
                                continue;

                        for (int h2_index = h1_index+1 ; h2_index < mctree->MCJet_Dtr_nmcdtrs ; h2_index++) {
                                if (mctree->MCJet_Dtr_IsMeson[h2_index] != 1 && mctree->MCJet_Dtr_IsBaryon[h2_index] != 1) 
                                        continue;

                                h2_4vector->SetPxPyPzE(mctree->MCJet_Dtr_PX[h2_index]/1000.,
                                                       mctree->MCJet_Dtr_PY[h2_index]/1000.,
                                                       mctree->MCJet_Dtr_PZ[h2_index]/1000.,
                                                       mctree->MCJet_Dtr_E[h2_index]/1000.);

                                if (!apply_chargedtrack_momentum_cuts(mctree->MCJet_Dtr_ThreeCharge[h2_index],
                                                                      h2_4vector->P(),
                                                                      h2_4vector->Pt(),
                                                                      h2_4vector->Eta())) 
                                        continue;

                                vars_mc[0]  = mctree->MCJet_Dtr_ThreeCharge[h1_index]*mctree->MCJet_Dtr_ThreeCharge[h2_index];
                                vars_mc[1]  = h1_4vector->DeltaR(*h2_4vector);
                                vars_mc[2]  = h1_4vector->Eta();
                                vars_mc[3]  = h2_4vector->Eta();
                                vars_mc[4]  = h1_4vector->Rapidity();
                                vars_mc[5]  = h2_4vector->Rapidity();
                                vars_mc[6]  = h1_4vector->P();
                                vars_mc[7]  = h2_4vector->P();
                                vars_mc[8]  = h1_4vector->Pt();
                                vars_mc[9]  = h2_4vector->Pt();
                                vars_mc[10] = Jet_4vector->Eta();
                                vars_mc[11] = Jet_4vector->Pt();
                                vars_mc[12] = weight(h1_4vector->Pt(), h2_4vector->Pt(), Jet_4vector->Pt());

                                // Fill the TNtuple
                                ntuple_mc->Fill(vars_mc);
                        }
                }
                
                vars_mc_jet[0] = true_Jet_4vector->Pt();
                vars_mc_jet[1] = true_Jet_4vector->E();
                vars_mc_jet[2] = mctree->MCJet_Dtr_nmcdtrs;
                vars_mc_jet[3] = (reco_passed) ? Jet_4vector->Pt() : -999;
                vars_mc_jet[4] = (reco_passed) ? Jet_4vector->E()  : -999;
                vars_mc_jet[5] = (reco_passed) ? mctree->MCJet_recojet_nrecodtrs : -999;
                vars_mc_jet[6] = true_Jet_4vector->Eta();

                ntuple_mc_jet->Fill(vars_mc_jet);
                
                last_eventNum = mctree->eventNumber;
        }

        fout->cd();
        ntuple_reco->Write();
        ntuple_mc->Write();
        ntuple_reco_jet->Write();
        ntuple_mc_jet->Write();
        fout->Close();

        std::cout<<std::endl;

        return 0;
}
