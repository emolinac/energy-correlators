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
        TFile* fout = new TFile((output_folder + namef_ntuple_jet_efficiency).c_str(),"RECREATE");
        
        // Declare the TTrees to be used to build the ntuples
        TZJetsMC* mctree = new TZJetsMC();

        // Create Ntuples
        TNtuple* ntuple = new TNtuple(name_ntuple_jetefficiency.c_str(),"",ntuple_jetefficiency_vars); 
        ntuple->SetAutoSave(0);
        
        // Create necessary 4vectors
        TLorentzVector* Jet_4vector   = new TLorentzVector();
        TLorentzVector* Z0_4vector    = new TLorentzVector();
        TLorentzVector* mum_4vector   = new TLorentzVector();
        TLorentzVector* mup_4vector   = new TLorentzVector();
        TLorentzVector* h_4vector     = new TLorentzVector();
        
        TLorentzVector* true_Jet_4vector   = new TLorentzVector();
        TLorentzVector* true_Z0_4vector    = new TLorentzVector();
        TLorentzVector* true_mum_4vector   = new TLorentzVector();
        TLorentzVector* true_mup_4vector   = new TLorentzVector();

        int eventNum;
        unsigned long long last_eventNum = 0;
        int events = 0;
        bool maxjetpT_found = false;
        
        // Define array carrying the variables
        float vars[Nvars_jetefficiency];
        
        // Fill the matched jets Ntuple
        
        for (int evt = 0 ; evt < mctree->fChain->GetEntries() ; evt++)
        {
                // Access entry of tree
                mctree->GetEntry(evt);

                if (evt%10000 == 0) {
                        double percentage = 100*evt/mctree->fChain->GetEntries();
                        std::cout<<"\r"<<percentage<<"\% jets processed."<< std::flush;
                }

                if (evt != 0)
                        if (last_eventNum == mctree->eventNumber) 
                                continue;
                
                // Apply PV cut
                if (mctree->nPVs!=1) 
                        continue;

                // Set Jet-associated 4 vectors and apply cuts
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

                int ndtrs_mcjet = 0;

                for (int h_index = 0 ; h_index < mctree->MCJet_Dtr_nmcdtrs ; h_index++) {
                        if (mctree->MCJet_Dtr_IsBaryon[h_index] != 1&&mctree->MCJet_Dtr_IsMeson[h_index] != 1) 
                                continue;

                        h_4vector->SetPxPyPzE(mctree->MCJet_Dtr_PX[h_index]/1000.,
                                              mctree->MCJet_Dtr_PY[h_index]/1000.,
                                              mctree->MCJet_Dtr_PZ[h_index]/1000.,
                                              mctree->MCJet_Dtr_E[h_index]/1000.);

                        if (!apply_chargedtrack_momentum_cuts(mctree->MCJet_Dtr_ThreeCharge[h_index],
                                                              h_4vector->P(),
                                                              h_4vector->Pt(),
                                                              h_4vector->Eta())) 
                                continue;

                        ndtrs_mcjet++;
                }

                if (ndtrs_mcjet < 2)
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

                        int ndtrs_mcrecojet = 0;

                        for (int h_index = 0 ; h_index < mctree->MCJet_recojet_nrecodtrs ; h_index++) {
                                if (mctree->MCJet_recojet_Dtr_IsBaryon[h_index] != 1&&mctree->MCJet_recojet_Dtr_IsMeson[h_index] != 1) 
                                        continue;

                                h_4vector->SetPxPyPzE(mctree->MCJet_recojet_Dtr_PX[h_index]/1000.,
                                                      mctree->MCJet_recojet_Dtr_PY[h_index]/1000.,
                                                      mctree->MCJet_recojet_Dtr_PZ[h_index]/1000.,
                                                      mctree->MCJet_recojet_Dtr_E[h_index]/1000.);

                                if (!apply_chargedtrack_cuts(mctree->MCJet_recojet_Dtr_ThreeCharge[h_index],
                                                             h_4vector->P(),
                                                             h_4vector->Pt(),
                                                             mctree->MCJet_recojet_Dtr_TrackChi2[h_index]/mctree->MCJet_recojet_Dtr_TrackNDF[h_index],
                                                             mctree->MCJet_recojet_Dtr_ProbNNghost[h_index],
                                                             h_4vector->Eta())) 
                                        continue;

                                ndtrs_mcrecojet++;
                        }

                        if (apply_jet_cuts(Jet_4vector->Eta(), Jet_4vector->Pt())&&\
                            apply_muon_cuts(Jet_4vector->DeltaR(*mum_4vector), mum_4vector->Pt(), mum_4vector->Eta())&&\
                            apply_muon_cuts(Jet_4vector->DeltaR(*mup_4vector), mup_4vector->Pt(), mup_4vector->Eta())&&\
                            apply_zboson_cuts(TMath::Abs(Jet_4vector->DeltaPhi(*Z0_4vector)), Z0_4vector->M())&&(ndtrs_mcrecojet > 1)) 
                                reco_passed = true;
                }
                
                vars[0] = true_Jet_4vector->Pt();
                vars[1] = true_Jet_4vector->E();
                vars[2] = mctree->MCJet_Dtr_nmcdtrs;
                vars[3] = (reco_passed) ? Jet_4vector->Pt() : -999;
                vars[4] = (reco_passed) ? Jet_4vector->E()  : -999;
                vars[5] = (reco_passed) ? mctree->MCJet_recojet_nrecodtrs : -999;
                vars[6] = true_Jet_4vector->Eta();

                last_eventNum = mctree->eventNumber;

                ntuple->Fill(vars);
        }

        fout->cd();
        ntuple->Write();
        fout->Close();

        std::cout<<std::endl;
        
        return 0;
}
