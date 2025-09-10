#include <iostream>
#include "TZJetsMC.h"
#include "TZJetsMC.C"
#include "TZJetsMCRecoCTCorr.h"
#include "TZJetsMCRecoCTCorr.C"
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
        TFile* fout = new TFile((output_folder + namef_ntuple_jet_purity_ct).c_str(),"RECREATE");
        
        // Declare the TTrees to be used to build the ntuples
        TZJetsMCRecoCTCorr* mcrecotree = new TZJetsMCRecoCTCorr();

        // Create Ntuples
        TNtuple* ntuple_jet_match = new TNtuple(name_ntuple_jetpurity.c_str(),"Reco jet matched 2 MC",ntuple_jetpurity_vars); 
        
        ntuple_jet_match->SetAutoSave(0);
        
        // Create necessary 4vectors
        TLorentzVector* Jet_4vector   = new TLorentzVector();
        TLorentzVector* Z0_4vector    = new TLorentzVector();
        TLorentzVector* mum_4vector   = new TLorentzVector();
        TLorentzVector* mup_4vector   = new TLorentzVector();
        TLorentzVector* h_4vector     = new TLorentzVector();
        
        TLorentzVector* true_Jet_4vector = new TLorentzVector();
        TLorentzVector* true_Z0_4vector  = new TLorentzVector();
        TLorentzVector* true_mum_4vector = new TLorentzVector();
        TLorentzVector* true_mup_4vector = new TLorentzVector();
        
        int eventNum;
        unsigned long long last_eventNum = 0;
        int events = 0;
        bool maxjetpT_found = false;
        
        // Define array carrying the variables
        float vars[Nvars_jetpurity];

        // Fill the matched jets Ntuple
        for (int evt = 0 ; evt < mcrecotree->fChain->GetEntries() ; evt++)
        {
                // Access entry of tree
                mcrecotree->GetEntry(evt);

                if (evt%10000 == 0) {
                        double percentage = 100*evt/mcrecotree->fChain->GetEntries();
                        std::cout<<"\r"<<percentage<<"\% jets processed."<< std::flush;
                }
                
                if (evt != 0)
                        if (last_eventNum == mcrecotree->eventNumber) 
                                continue;

                if (maxjetpT_found) 
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

                double ndtrs_mcrecojet = 0; 
                
                for (int h_index = 0 ; h_index < mcrecotree->Jet_NDtr ; h_index++) {
                        if (mcrecotree->Jet_Dtr_IsMeson[h_index] != 1 && mcrecotree->Jet_Dtr_IsBaryon[h_index] != 1) 
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

                        ndtrs_mcrecojet++;
                }

                if (ndtrs_mcrecojet < 2)
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

                        int ndtrsmatched_mcrecojet = 0;
                        for (int h_index = 0 ; h_index < mcrecotree->Jet_NDtr ; h_index++) {
                                if (mcrecotree->Jet_Dtr_TRUE_IsMeson[h_index] != 1 && mcrecotree->Jet_Dtr_TRUE_IsBaryon[h_index] != 1) 
                                        continue;

                                h_4vector->SetPxPyPzE(mcrecotree->Jet_Dtr_TRUE_PX[h_index]/1000.,
                                                      mcrecotree->Jet_Dtr_TRUE_PY[h_index]/1000.,
                                                      mcrecotree->Jet_Dtr_TRUE_PZ[h_index]/1000.,
                                                      mcrecotree->Jet_Dtr_TRUE_E[h_index]/1000.);

                                if (!apply_chargedtrack_momentum_cuts(mcrecotree->Jet_Dtr_TRUE_ThreeCharge[h_index],
                                                                      h_4vector->P(),
                                                                      h_4vector->Pt(),
                                                                      h_4vector->Eta())) 
                                        continue;

                                ndtrsmatched_mcrecojet++;
                        }

                        // if (apply_jet_cuts(true_Jet_4vector->Eta(),true_Jet_4vector->Pt())&&\
                        //     apply_muon_cuts(true_Jet_4vector->DeltaR(*true_mum_4vector),true_mum_4vector->Pt(),true_mum_4vector->Eta())&&\
                        //     apply_muon_cuts(true_Jet_4vector->DeltaR(*true_mup_4vector),true_mup_4vector->Pt(),true_mup_4vector->Eta())&&\
                        //     apply_zboson_cuts(TMath::Abs(true_Jet_4vector->DeltaPhi(*true_Z0_4vector)),true_Z0_4vector->M())&&(ndtrsmatched_mcrecojet > 1)) 
                        //         truth_passed = true;

                        if (apply_jet_cuts(true_Jet_4vector->Eta(),true_Jet_4vector->Pt())) 
                                truth_passed = true;
                }
                        
                vars[0] = Jet_4vector->Pt();
                vars[1] = Jet_4vector->E();
                vars[2] = mcrecotree->Jet_NDtr;
                vars[3] = (truth_passed) ? true_Jet_4vector->Pt() : -999 ;
                vars[4] = (truth_passed) ? true_Jet_4vector->E()  : -999 ;
                vars[5] = (truth_passed) ? mcrecotree->Jet_mcjet_nmcdtrs : -999 ;
                vars[6] = (truth_passed) ? Jet_4vector->DeltaR(*true_Jet_4vector) : -999;
                vars[7] = Jet_4vector->Eta();

                last_eventNum = mcrecotree->eventNumber;
                
                // Fill the TNtuple
                ntuple_jet_match->Fill(vars); 
        } 

        fout->cd();
        ntuple_jet_match->Write();
        fout->Close();

        std::cout<<std::endl;

        return 0;
}
