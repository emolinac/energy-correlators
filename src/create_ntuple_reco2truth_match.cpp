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
#include "analysis-functions.h"
#include "directories.h"
#include "names.h"
#include "syst-jes-jer.h"

int main(int argc, char* argv[])
{
        bool get_nominal = false;
        bool get_jes     = false;
        bool get_jer     = false;

        if (argc < 2 || std::string(argv[1]) == "--help") {
                std::cout<<"You have to pass one argument to this code. Possible options are:"<<std::endl;
                std::cout<<"--get-nominal : Gets the nominal pair corrections."<<std::endl;
                std::cout<<"--get-jes     : Gets the pair corrections with the JES variation."<<std::endl;
                std::cout<<"--get-jer     : Gets the pair corrections with the JER variation."<<std::endl;

                return 0;
        }

        if (std::string(argv[1]) == "--get-nominal")
                get_nominal = true;
        else if (std::string(argv[1]) == "--get-jes")
                get_jes = true;
        else if (std::string(argv[1]) == "--get-jer")
                get_jer = true;
        else {
                std::cout<<"No valid options provided"<<std::endl;
                return 0;
        }
        
        if ((get_jes && syst_jes_array[0] == -999)||(get_jer && syst_jer_array[0] == -999)) {
                std::cout<<"No values were found for custom JEC. Set the values related to the JEC:"<<std::endl;
                std::cout<<"1 - In the bin folder execute: ./create_jes_jer_ntuple"<<std::endl;
                std::cout<<"2 - In the src-analysis-jets folder execute: root -b -q macro_print_jes_chisquare.cpp"<<std::endl;
                std::cout<<"3 - Copy the indicated output in include/syst-jes-jer.h"<<std::endl;
                std::cout<<"4 - In the src-analysis-jets folder execute: root -b -q macro_print_jer_chisquare.cpp"<<std::endl;
                std::cout<<"5 - Copy the indicated output in include/syst-jes-jer.h"<<std::endl;

                return 0;
        }

        // Create output file
        TFile* fout = new TFile((output_folder + namef_reco_corrections[argv[1]]).c_str(),"RECREATE");
        
        // Declare the TTrees to be used to build the ntuples
        TZJetsMCReco* mcrecotree = new TZJetsMCReco();

        // Create Ntuples
        TNtuple* ntuple_reco = new TNtuple(name_ntuple_correction_reco.c_str(),"",ntuple_corrections_reco_vars); 
        TNtuple* ntuple_mc   = new TNtuple(name_ntuple_correction_mc.c_str()  ,"",ntuple_corrections_mc_vars); 

        TNtuple* ntuple_reco_jet = new TNtuple(name_ntuple_jet_reco2truth_match.c_str(),"",ntuple_jetpurity_vars);
        
        ntuple_reco->SetAutoSave(0);
        ntuple_mc->SetAutoSave(0);
        ntuple_reco_jet->SetAutoSave(0);
        
        // Create necessary 4vectors
        TLorentzVector* Jet_4vector   = new TLorentzVector();
        TLorentzVector* Z0_4vector    = new TLorentzVector();
        TLorentzVector* mum_4vector   = new TLorentzVector();
        TLorentzVector* mup_4vector   = new TLorentzVector();
        TLorentzVector* h1_4vector    = new TLorentzVector();
        TLorentzVector* h2_4vector    = new TLorentzVector();
        TLorentzVector* h_4vector     = new TLorentzVector();

        TLorentzVector* true_Jet_4vector = new TLorentzVector();
        TLorentzVector* true_h1_4vector  = new TLorentzVector();
        TLorentzVector* true_h2_4vector  = new TLorentzVector();
        
        unsigned long long last_eventNum = 0;
        bool maxjetpT_found = false;
        
        // Define array carrying the variables
        float vars_reco[Nvars_corrections_mcreco];
        float vars_mc[Nvars_corrections_mc];
        float vars_reco_jet[Nvars_jetpurity];

        TRandom3* rndm = new TRandom3(0);

        std::cout<<"Processing events ..."<<std::endl;
        for (int evt = 0 ; evt < mcrecotree->fChain->GetEntries() ; evt++) {
                if (evt%10000 == 0) {
                        double percentage = 100*evt/mcrecotree->fChain->GetEntries();
                        std::cout<<"\r"<<percentage<<"\%"<< std::flush;
                }

                mcrecotree->GetEntry(evt);

                if (evt != 0)
                        if (last_eventNum == mcrecotree->eventNumber) 
                                continue;

                if (mcrecotree->nPV != 1)
                        continue;

                bool mum_trigger = (mcrecotree->mum_L0MuonEWDecision_TOS == 1 && 
                                    mcrecotree->mum_Hlt1SingleMuonHighPTDecision_TOS == 1 && 
                                    mcrecotree->mum_Hlt2EWSingleMuonVHighPtDecision_TOS == 1);

                bool mup_trigger = (mcrecotree->mup_L0MuonEWDecision_TOS == 1 && 
                                    mcrecotree->mup_Hlt1SingleMuonHighPTDecision_TOS == 1 && 
                                    mcrecotree->mup_Hlt2EWSingleMuonVHighPtDecision_TOS == 1);                                    

                if (!mum_trigger && !mup_trigger) 
                        continue;

                Jet_4vector->SetPxPyPzE(mcrecotree->Jet_PX/1000.,
                                        mcrecotree->Jet_PY/1000.,
                                        mcrecotree->Jet_PZ/1000.,
                                        mcrecotree->Jet_PE/1000.);

                if (!apply_jet_cuts(Jet_4vector->Eta(), Jet_4vector->Pt())) 
                        continue;
                
                if (get_jes) {
                        double beta_star = 1.;
                        
                        for (int jet_pt_bin = 0 ; jet_pt_bin < nbin_jet_pt_unfolding ; jet_pt_bin++)
                                if (Jet_4vector->Pt()>unfolding_jet_pt_binning[jet_pt_bin]&&Jet_4vector->Pt()<unfolding_jet_pt_binning[jet_pt_bin + 1])
                                        beta_star = syst_jes_array[jet_pt_bin];

                        if (beta_star == 1)
                                continue;

                        double new_jes_cor_effect = std::abs(1. - beta_star);

                        if (rndm->Integer(2))
                                beta_star = 1 + new_jes_cor_effect;
                        else
                                beta_star = 1 - new_jes_cor_effect;

                        Jet_4vector->SetPxPyPzE(beta_star*mcrecotree->Jet_PX/1000.,
                                                beta_star*mcrecotree->Jet_PY/1000.,
                                                beta_star*mcrecotree->Jet_PZ/1000.,
                                                beta_star*mcrecotree->Jet_PE/1000.);
                } else if (get_jer) {
                        double alpha_star = 1.;

                        for (int jet_pt_bin = 0 ; jet_pt_bin < nbin_jet_pt_unfolding ; jet_pt_bin++)
                                if (Jet_4vector->Pt()>unfolding_jet_pt_binning[jet_pt_bin]&&Jet_4vector->Pt()<unfolding_jet_pt_binning[jet_pt_bin + 1]) 
                                        alpha_star = syst_jer_array[jet_pt_bin];
                        
                        if (alpha_star == 1)
                                continue;
                        
                        double smearing_factor = rndm->Gaus(1., alpha_star);
                        
                        Jet_4vector->SetPxPyPzE(smearing_factor*mcrecotree->Jet_PX/1000.,
                                                smearing_factor*mcrecotree->Jet_PY/1000.,
                                                smearing_factor*mcrecotree->Jet_PZ/1000.,
                                                smearing_factor*mcrecotree->Jet_PE/1000.);
                }

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
                
                bool truth_jet_passed = false;
                
                if (mcrecotree->Jet_mcjet_nmcdtrs!=-999) {
                        true_Jet_4vector->SetPxPyPzE(mcrecotree->Jet_mcjet_PX/1000.,
                                                     mcrecotree->Jet_mcjet_PY/1000.,
                                                     mcrecotree->Jet_mcjet_PZ/1000.,
                                                     mcrecotree->Jet_mcjet_PE/1000.);

                        if (apply_jet_cuts(true_Jet_4vector->Eta(),true_Jet_4vector->Pt())) 
                                truth_jet_passed = true;
                }

                vars_reco_jet[0] = Jet_4vector->Pt();
                vars_reco_jet[1] = Jet_4vector->E();
                vars_reco_jet[2] = mcrecotree->Jet_NDtr;
                vars_reco_jet[3] = (truth_jet_passed) ? true_Jet_4vector->Pt() : -999 ;
                vars_reco_jet[4] = (truth_jet_passed) ? true_Jet_4vector->E()  : -999 ;
                vars_reco_jet[5] = (truth_jet_passed) ? mcrecotree->Jet_mcjet_nmcdtrs : -999 ;
                vars_reco_jet[6] = (truth_jet_passed) ? Jet_4vector->DeltaR(*true_Jet_4vector) : -999;
                vars_reco_jet[7] = Jet_4vector->Eta();
                

                double jet_ndtr_mc   = 0;
                double jet_ndtr_reco = 0;

                if (truth_jet_passed) {
                        for (int h1_index = 0 ; h1_index < mcrecotree->Jet_NDtr ; h1_index++) {
                                if (abs(mcrecotree->Jet_Dtr_ID[h1_index]) < 100) 
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
                                int key1_sim   = 0;
                                if (mcrecotree->Jet_Dtr_TRUE_ETA[h1_index] != -999 && std::abs(mcrecotree->Jet_Dtr_TRUE_ID[h1_index]) > 100) {
                                        key1_match++;

                                        for(int i = 0 ; i < mcrecotree->Jet_mcjet_nmcdtrs ; i++) {
                                                if(mcrecotree->Jet_Dtr_TRUE_KEY[h1_index] == mcrecotree->Jet_mcjet_dtrKeys[i]) {
                                                        key1_sim++;
                                                        break;
                                                }
                                        }

                                        // Verifies if the matched particle is in the matched jet
                                        if (key1_sim == 0) 
                                                key1_match = 0;

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
                                        if (abs(mcrecotree->Jet_Dtr_ID[h2_index]) < 100)
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
                                        int key2_sim   = 0;
                                        if (mcrecotree->Jet_Dtr_TRUE_ETA[h2_index] != -999 && std::abs(mcrecotree->Jet_Dtr_TRUE_ID[h2_index]) > 100) {
                                                key2_match++;

                                                for(int i = 0 ; i < mcrecotree->Jet_mcjet_nmcdtrs ; i++) {
                                                        if(mcrecotree->Jet_Dtr_TRUE_KEY[h2_index] == mcrecotree->Jet_mcjet_dtrKeys[i]) {
                                                                key2_sim++;
                                                                break;
                                                        }
                                                }

                                                // Verifies if the matched particle is in the matched jet
                                                if (key2_sim == 0) 
                                                        key2_match = 0;

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

                        // Fill the mc hadrons
                        for (int h1_index = 0 ; h1_index < mcrecotree->Jet_mcjet_nmcdtrs ; h1_index++) {
                                if (abs(mcrecotree->Jet_mcjet_dtrID[h1_index]) < 100)
                                        continue;

                                h1_4vector->SetPxPyPzE(mcrecotree->Jet_mcjet_dtrPX[h1_index]/1000.,
                                                mcrecotree->Jet_mcjet_dtrPY[h1_index]/1000.,
                                                mcrecotree->Jet_mcjet_dtrPZ[h1_index]/1000.,
                                                mcrecotree->Jet_mcjet_dtrE[h1_index]/1000.);
                                
                                if (!apply_chargedtrack_momentum_cuts(mcrecotree->Jet_mcjet_dtrThreeCharge[h1_index],
                                                                h1_4vector->P(),
                                                                h1_4vector->Pt(),
                                                                h1_4vector->Eta())) 
                                        continue;

                                jet_ndtr_mc++;

                                for (int h2_index = h1_index+1 ; h2_index < mcrecotree->Jet_mcjet_nmcdtrs ; h2_index++) {
                                        if (abs(mcrecotree->Jet_mcjet_dtrID[h2_index]) < 100)
                                                continue;

                                        h2_4vector->SetPxPyPzE(mcrecotree->Jet_mcjet_dtrPX[h2_index]/1000.,
                                                        mcrecotree->Jet_mcjet_dtrPY[h2_index]/1000.,
                                                        mcrecotree->Jet_mcjet_dtrPZ[h2_index]/1000.,
                                                        mcrecotree->Jet_mcjet_dtrE[h2_index]/1000.);

                                        if (!apply_chargedtrack_momentum_cuts(mcrecotree->Jet_mcjet_dtrThreeCharge[h2_index],
                                                                        h2_4vector->P(),
                                                                        h2_4vector->Pt(),
                                                                        h2_4vector->Eta())) 
                                                continue;

                                        vars_mc[0]  = mcrecotree->Jet_mcjet_dtrThreeCharge[h1_index]*mcrecotree->Jet_mcjet_dtrThreeCharge[h2_index];
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
                }
                
                last_eventNum = mcrecotree->eventNumber;

                ntuple_reco_jet->Fill(vars_reco_jet);
        }

        fout->cd();
        ntuple_reco->Write();
        ntuple_mc->Write();
        ntuple_reco_jet->Write();
        fout->Close();

        std::cout<<std::endl;

        return 0;
}
