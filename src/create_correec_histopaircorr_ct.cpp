#include <iostream>
#include "TZJetsTruth.h"
#include "TZJetsTruth.C"
#include "TZJetsPseudoData.h"
#include "TZJetsPseudoData.C"
#include "TZJets2016Data.h"
#include "TZJets2016Data.C"
#include "TROOT.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TH3.h"
#include "TLatex.h"
#include "analysis-constants.h"
#include "analysis-binning.h"
#include "analysis-cuts.h"
#include "analysis-functions.h"
#include "directories.h"
#include "names.h"
#include "utils-algorithms.h"

int main()
{
        // Create output file
        TFile* fout = new TFile((output_folder + namef_paircorr_histos_ct).c_str(),"RECREATE");
        gROOT->cd();
        
        // Declare the TTrees to be used to build the ntuples
        TZJetsPseudoData* pseudodata = new TZJetsPseudoData();
        TZJetsTruth*      truthdata  = new TZJetsTruth();

        // Open correction files
        TFile* fcorrections_pair = new TFile((output_folder + namef_ntuple_eec_paircorrections_ct).c_str());
        TFile* fpurity_jet       = new TFile((output_folder + namef_ntuple_jet_purity_ct).c_str());
        TFile* fefficiency_jet   = new TFile((output_folder + namef_ntuple_jet_efficiency_ct).c_str());
        
        TFile* fefficiency_muon_2017_id  = new TFile((muons_folder + "IDEff_Data_2017.root").c_str());
        TFile* fefficiency_muon_2017_trk = new TFile((muons_folder + "TRKEff_Data_2017.root").c_str());
        TFile* fefficiency_muon_2017_trg = new TFile((muons_folder + "TRGEff_Data_2017.root").c_str());
        
        // Create Ntuples
        TNtuple* ntuple_purity          = (TNtuple*) fcorrections_pair->Get((name_ntuple_correction_reco.c_str()));
        TNtuple* ntuple_efficiency_mc   = (TNtuple*) fcorrections_pair->Get((name_ntuple_correction_mc.c_str()));
        TNtuple* ntuple_efficiency_reco = (TNtuple*) fcorrections_pair->Get((name_ntuple_correction_reco.c_str()));
        TNtuple* ntuple_purity_jet      = (TNtuple*) fpurity_jet->Get((name_ntuple_jetpurity.c_str()));
        TNtuple* ntuple_efficiency_jet  = (TNtuple*) fefficiency_jet->Get((name_ntuple_jetefficiency.c_str()));
        
        // Muon corrections
        TH2D* h2_muon_2017_ideff_data  = (TH2D*) fefficiency_muon_2017_id->Get("Hist_ALL_2017_ETA_PT_Eff");
        TH2D* h2_muon_2017_trkeff_data = (TH2D*) fefficiency_muon_2017_trk->Get("Hist_ALL_2017_ETA_PT_Eff");
        TH2D* h2_muon_2017_trgeff_data = (TH2D*) fefficiency_muon_2017_trg->Get("Hist_ALL_2017_ETA_PT_Eff");
        
        // Jet corrections
        TH1F* hnum_pur_jet = new TH1F("hnum_pur_jet", "", nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
        TH1F* hden_pur_jet = new TH1F("hden_pur_jet", "", nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
        TH1F* hpurity_jet  = new TH1F("hpurity_jet" , "", nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
        hnum_pur_jet->Sumw2();
        hden_pur_jet->Sumw2();

        TH1F* hnum_eff_jet    = new TH1F("hnum_eff_jet"   , "", nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
        TH1F* hden_eff_jet    = new TH1F("hden_eff_jet"   , "", nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
        TH1F* hefficiency_jet = new TH1F("hefficiency_jet", "", nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
        hnum_eff_jet->Sumw2();
        hden_eff_jet->Sumw2();

        ntuple_purity_jet->Project("hnum_pur_jet", "jet_pt", "jet_pt_truth!=-999");
        ntuple_purity_jet->Project("hden_pur_jet", "jet_pt");
        ntuple_efficiency_jet->Project("hnum_eff_jet", "jet_pt_truth", "jet_pt!=-999");
        ntuple_efficiency_jet->Project("hden_eff_jet", "jet_pt_truth");

        hpurity_jet->Divide(hnum_pur_jet, hden_pur_jet, 1, 1, "B");
        hefficiency_jet->Divide(hnum_eff_jet, hden_eff_jet, 1, 1, "B");

        // Pair corrections
        TH2F* hnum_pur    = new TH2F("hnum_pur"   , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
        TH2F* hden_pur    = new TH2F("hden_pur"   , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
        TH2F* hpurity     = new TH2F("hpurity"    , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
        TH2F* hnum_eff    = new TH2F("hnum_eff"   , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
        TH2F* hden_eff    = new TH2F("hden_eff"   , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
        TH2F* hefficiency = new TH2F("hefficiency", "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
        
        hnum_pur->Sumw2();
        hden_pur->Sumw2();
        hnum_eff->Sumw2();
        hden_eff->Sumw2();
        
        ntuple_purity->Project("hnum_pur", "jet_pt:R_L", pair_matching_cut);
        ntuple_purity->Project("hden_pur", "jet_pt:R_L");
        ntuple_efficiency_reco->Project("hnum_eff", "jet_pt_truth:R_L_truth", pair_matching_cut);
        ntuple_efficiency_mc->Project("hden_eff", "jet_pt:R_L");

        hpurity->Divide(hnum_pur, hden_pur, 1, 1, "B");
        hefficiency->Divide(hnum_eff, hden_eff, 1, 1, "B");

        regularize_correction_factors(hpurity);
        regularize_correction_factors(hefficiency);

        hpurity->Smooth();
        hefficiency->Smooth();

        // Create necessary 4vectors
        TLorentzVector* Jet_4vector = new TLorentzVector();
        TLorentzVector* Z0_4vector  = new TLorentzVector();
        TLorentzVector* mum_4vector = new TLorentzVector();
        TLorentzVector* mup_4vector = new TLorentzVector();
        TLorentzVector* h1_4vector  = new TLorentzVector();
        TLorentzVector* h2_4vector  = new TLorentzVector();
        
        // EEC histo
        TH2F* h_eec       = new TH2F("h_eec"      ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
        TH2F* h_eec_truth = new TH2F("h_eec_truth","",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
        
        TH2F* h_npair       = new TH2F("h_npair"      ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
        TH2F* h_npair_truth = new TH2F("h_npair_truth","",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
        h_npair->Sumw2();
        h_npair_truth->Sumw2();
        
        TH1F* h_njet          = new TH1F("h_njet"         ,"",nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
        TH1F* h_njet_wmuoneff = new TH1F("h_njet_wmuoneff","",nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
        TH1F* h_njet_truth    = new TH1F("h_njet_truth"   ,"",nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
        h_njet->Sumw2();
        h_njet_truth->Sumw2();

        unsigned long long last_eventNum = 0;
        
        std::cout<<"Working with Pseudata."<<std::endl;
        for (int evt = 0 ; evt < pseudodata->fChain->GetEntries() ; evt++) {
                // Access entry of tree
                pseudodata->GetEntry(evt);

                if (evt%10000 == 0) {
                        double percentage = 100.*evt/pseudodata->fChain->GetEntries();
                        std::cout<<"\r"<<percentage<<"\% jets processed."<< std::flush;
                }

                if (evt != 0)
                        if (last_eventNum == pseudodata->eventNumber) 
                                continue;

                // Apply PV cut
                if (pseudodata->nPV != 1)
                        continue;

                // Apply trigger cut
                bool mum_trigger = (pseudodata->mum_L0MuonEWDecision_TOS == 1 && 
                                    pseudodata->mum_Hlt1SingleMuonHighPTDecision_TOS == 1 && 
                                    pseudodata->mum_Hlt2EWSingleMuonVHighPtDecision_TOS == 1);

                bool mup_trigger = (pseudodata->mup_L0MuonEWDecision_TOS == 1 && 
                                    pseudodata->mup_Hlt1SingleMuonHighPTDecision_TOS == 1 && 
                                    pseudodata->mup_Hlt2EWSingleMuonVHighPtDecision_TOS == 1);                                    

                if (!mum_trigger && !mup_trigger)
                        continue;
                
                // Set Jet-associated 4 vectors and apply cuts
                Jet_4vector->SetPxPyPzE(pseudodata->Jet_PX/1000.,
                                        pseudodata->Jet_PY/1000.,
                                        pseudodata->Jet_PZ/1000.,
                                        pseudodata->Jet_PE/1000.);
                if (!apply_jet_cuts(Jet_4vector->Eta(), Jet_4vector->Pt())) 
                        continue;
                
                mum_4vector->SetPxPyPzE(pseudodata->mum_PX/1000.,
                                        pseudodata->mum_PY/1000.,
                                        pseudodata->mum_PZ/1000.,
                                        pseudodata->mum_PE/1000.);

                if (!apply_muon_cuts(Jet_4vector->DeltaR(*mum_4vector), mum_4vector->Pt(), mum_4vector->Eta())) 
                        continue;
                
                mup_4vector->SetPxPyPzE(pseudodata->mup_PX/1000.,
                                        pseudodata->mup_PY/1000.,
                                        pseudodata->mup_PZ/1000.,
                                        pseudodata->mup_PE/1000.);

                if (!apply_muon_cuts(Jet_4vector->DeltaR(*mup_4vector), mup_4vector->Pt(), mup_4vector->Eta())) 
                        continue;
                
                Z0_4vector->SetPxPyPzE(mup_4vector->Px()+mum_4vector->Px(),
                                       mup_4vector->Py()+mum_4vector->Py(),
                                       mup_4vector->Pz()+mum_4vector->Pz(),
                                       mup_4vector->E() +mum_4vector->E());

                if (!apply_zboson_cuts(TMath::Abs(Jet_4vector->DeltaPhi(*Z0_4vector)), Z0_4vector->M())) 
                        continue;

                double mup_pt  = (mup_4vector->Pt() >= 70.) ? 69. : mup_4vector->Pt();
                double mum_pt  = (mum_4vector->Pt() >= 70.) ? 69. : mum_4vector->Pt();
                double mup_eta = mup_4vector->Eta();
                double mum_eta = mum_4vector->Eta();

                double mup_eff_id  = h2_muon_2017_ideff_data->GetBinContent(h2_muon_2017_ideff_data->FindBin(mup_eta, mup_pt));
                double mup_eff_trk = h2_muon_2017_trkeff_data->GetBinContent(h2_muon_2017_trkeff_data->FindBin(mup_eta, mup_pt));
                double mup_eff_trg = h2_muon_2017_trgeff_data->GetBinContent(h2_muon_2017_trgeff_data->FindBin(mup_eta, mup_pt));

                double mum_eff_id  = h2_muon_2017_ideff_data->GetBinContent(h2_muon_2017_ideff_data->FindBin(mum_eta, mum_pt));
                double mum_eff_trk = h2_muon_2017_trkeff_data->GetBinContent(h2_muon_2017_trkeff_data->FindBin(mum_eta, mum_pt));
                double mum_eff_trg = h2_muon_2017_trgeff_data->GetBinContent(h2_muon_2017_trgeff_data->FindBin(mum_eta, mum_pt));

                double jet_efficiency = hefficiency_jet->GetBinContent(hefficiency_jet->FindBin(Jet_4vector->Pt()));
                double jet_purity     = hpurity_jet->GetBinContent(hpurity_jet->FindBin(Jet_4vector->Pt()));
                
                double jet_efficiency_error = hefficiency_jet->GetBinError(hefficiency_jet->FindBin(Jet_4vector->Pt()));
                double jet_purity_error     = hpurity_jet->GetBinError(hpurity_jet->FindBin(Jet_4vector->Pt()));

                double muon_weight = 1./(mum_eff_id*mup_eff_id*mum_eff_trk*mup_eff_trk*(mum_eff_trg+mup_eff_trg-mum_eff_trg*mup_eff_trg));

                double jet_ndtr_nominal = 0;

                for (int h1_index = 0 ; h1_index < pseudodata->Jet_NDtr ; h1_index++) {
                        // Skip non-hadronic particles
                        if (pseudodata->Jet_Dtr_IsMeson[h1_index] != 1 && pseudodata->Jet_Dtr_IsBaryon[h1_index] != 1)
                                continue;

                        h1_4vector->SetPxPyPzE(pseudodata->Jet_Dtr_PX[h1_index]/1000.,
                                               pseudodata->Jet_Dtr_PY[h1_index]/1000.,
                                               pseudodata->Jet_Dtr_PZ[h1_index]/1000.,
                                               pseudodata->Jet_Dtr_E[h1_index]/1000.);

                        if (!apply_chargedtrack_cuts(pseudodata->Jet_Dtr_ThreeCharge[h1_index],
                                                     h1_4vector->P(),
                                                     h1_4vector->Pt(),
                                                     pseudodata->Jet_Dtr_TrackChi2[h1_index]/pseudodata->Jet_Dtr_TrackNDF[h1_index],
                                                     pseudodata->Jet_Dtr_ProbNNghost[h1_index],
                                                     h1_4vector->Eta())) 
                                continue;

                        jet_ndtr_nominal++;
                }

                if (jet_ndtr_nominal < 2)
                        continue;
                
                // Loop over hadron 1
                for (int h1_index = 0 ; h1_index < pseudodata->Jet_NDtr ; h1_index++) {
                        // Skip non-hadronic particles
                        if (pseudodata->Jet_Dtr_IsMeson[h1_index] != 1 && pseudodata->Jet_Dtr_IsBaryon[h1_index] != 1) 
                                continue;

                        h1_4vector->SetPxPyPzE(pseudodata->Jet_Dtr_PX[h1_index]/1000.,
                                               pseudodata->Jet_Dtr_PY[h1_index]/1000.,
                                               pseudodata->Jet_Dtr_PZ[h1_index]/1000.,
                                               pseudodata->Jet_Dtr_E[h1_index]/1000.);

                        if (!apply_chargedtrack_cuts(pseudodata->Jet_Dtr_ThreeCharge[h1_index],
                                                     h1_4vector->P(),
                                                     h1_4vector->Pt(),
                                                     pseudodata->Jet_Dtr_TrackChi2[h1_index]/pseudodata->Jet_Dtr_TrackNDF[h1_index],
                                                     pseudodata->Jet_Dtr_ProbNNghost[h1_index],
                                                     h1_4vector->Eta())) 
                                continue;
                        
                        // Loop over hadron 2
                        for (int h2_index = h1_index+1 ; h2_index < pseudodata->Jet_NDtr ; h2_index++) {
                                // Skip non-hadronic particles
                                if (pseudodata->Jet_Dtr_IsMeson[h2_index] != 1 && pseudodata->Jet_Dtr_IsBaryon[h2_index] != 1) 
                                        continue;

                                h2_4vector->SetPxPyPzE(pseudodata->Jet_Dtr_PX[h2_index]/1000.,
                                                       pseudodata->Jet_Dtr_PY[h2_index]/1000.,
                                                       pseudodata->Jet_Dtr_PZ[h2_index]/1000.,
                                                       pseudodata->Jet_Dtr_E[h2_index]/1000.);

                                if (!apply_chargedtrack_cuts(pseudodata->Jet_Dtr_ThreeCharge[h2_index],
                                                             h2_4vector->P(),
                                                             h2_4vector->Pt(),
                                                             pseudodata->Jet_Dtr_TrackChi2[h2_index]/pseudodata->Jet_Dtr_TrackNDF[h2_index],
                                                             pseudodata->Jet_Dtr_ProbNNghost[h2_index],
                                                             h2_4vector->Eta())) 
                                        continue;

                                double momentum_weight = weight(h1_4vector->Pt(), h2_4vector->Pt(), Jet_4vector->Pt());

                                h_eec->Fill(h1_4vector->DeltaR(*h2_4vector), Jet_4vector->Pt(), momentum_weight);
                                h_npair->Fill(h1_4vector->DeltaR(*h2_4vector), Jet_4vector->Pt());
                        }
                }

                h_njet->Fill(Jet_4vector->Pt());
                h_njet_wmuoneff->Fill(Jet_4vector->Pt(), muon_weight);
                
                last_eventNum = pseudodata->eventNumber;
        }

        std::cout<<std::endl;
        std::cout<<"Working the Truth"<<std::endl;
        for (int evt = 0 ; evt < truthdata->fChain->GetEntries() ; evt++) {
                truthdata->GetEntry(evt);

                if (evt%10000 == 0) {
                        double percentage = 100.*evt/truthdata->fChain->GetEntries();
                        std::cout<<"\r"<<percentage<<"\% jets processed."<< std::flush;
                }

                if (truthdata->MCJet_recojet_PX == -999)
                        continue;

                if (evt != 0)
                        if (last_eventNum == truthdata->eventNumber) 
                                continue;

                // Apply PV cut
                if (truthdata->nPVs != 1)
                        continue;

                // Set Jet-associated 4 vectors and apply cuts
                Jet_4vector->SetPxPyPzE(truthdata->MCJet_PX/1000.,
                                        truthdata->MCJet_PY/1000.,
                                        truthdata->MCJet_PZ/1000.,
                                        truthdata->MCJet_PE/1000.);

                if (!apply_jet_cuts(Jet_4vector->Eta(), Jet_4vector->Pt())) 
                        continue;

                mum_4vector->SetPxPyPzE(truthdata->MCJet_truth_mum_PX/1000.,
                                        truthdata->MCJet_truth_mum_PY/1000.,
                                        truthdata->MCJet_truth_mum_PZ/1000.,
                                        truthdata->MCJet_truth_mum_PE/1000.);

                if (!apply_muon_cuts(Jet_4vector->DeltaR(*mum_4vector), mum_4vector->Pt(), mum_4vector->Eta())) 
                        continue;

                mup_4vector->SetPxPyPzE(truthdata->MCJet_truth_mup_PX/1000.,
                                        truthdata->MCJet_truth_mup_PY/1000.,
                                        truthdata->MCJet_truth_mup_PZ/1000.,
                                        truthdata->MCJet_truth_mup_PE/1000.);

                if (!apply_muon_cuts(Jet_4vector->DeltaR(*mup_4vector), mup_4vector->Pt(), mup_4vector->Eta())) 
                        continue;

                Z0_4vector->SetPxPyPzE(mup_4vector->Px()+mum_4vector->Px(),
                                       mup_4vector->Py()+mum_4vector->Py(),
                                       mup_4vector->Pz()+mum_4vector->Pz(),
                                       mup_4vector->E() +mum_4vector->E());

                if (!apply_zboson_cuts(TMath::Abs(Jet_4vector->DeltaPhi(*Z0_4vector)), Z0_4vector->M())) 
                        continue;

                double jet_ndtr_nominal = 0;

                for (int h1_index = 0 ; h1_index < truthdata->MCJet_Dtr_nmcdtrs ; h1_index++) {
                        if (truthdata->MCJet_Dtr_IsMeson[h1_index] != 1 && truthdata->MCJet_Dtr_IsBaryon[h1_index] != 1) 
                                continue;

                        h1_4vector->SetPxPyPzE(truthdata->MCJet_Dtr_PX[h1_index]/1000.,
                                               truthdata->MCJet_Dtr_PY[h1_index]/1000.,
                                               truthdata->MCJet_Dtr_PZ[h1_index]/1000., 
                                               truthdata->MCJet_Dtr_E[h1_index]/1000.);

                        if (!apply_chargedtrack_momentum_cuts(truthdata->MCJet_Dtr_ThreeCharge[h1_index],
                                                              h1_4vector->P(),
                                                              h1_4vector->Pt(),
                                                              h1_4vector->Eta())) 
                                continue;
                        
                        jet_ndtr_nominal++;
                }

                if (jet_ndtr_nominal < 2)
                        continue;

                for (int h1_index = 0 ; h1_index < truthdata->MCJet_Dtr_nmcdtrs ; h1_index++) {
                        if (truthdata->MCJet_Dtr_IsMeson[h1_index] != 1 && truthdata->MCJet_Dtr_IsBaryon[h1_index] != 1) 
                                continue;

                        h1_4vector->SetPxPyPzE(truthdata->MCJet_Dtr_PX[h1_index]/1000.,
                                               truthdata->MCJet_Dtr_PY[h1_index]/1000.,
                                               truthdata->MCJet_Dtr_PZ[h1_index]/1000., 
                                               truthdata->MCJet_Dtr_E[h1_index]/1000.);

                        if (!apply_chargedtrack_momentum_cuts(truthdata->MCJet_Dtr_ThreeCharge[h1_index],
                                                              h1_4vector->P(),
                                                              h1_4vector->Pt(),
                                                              h1_4vector->Eta())) 
                                continue;

                        for (int h2_index = h1_index+1 ; h2_index < truthdata->MCJet_Dtr_nmcdtrs ; h2_index++) {
                                if (truthdata->MCJet_Dtr_IsMeson[h2_index] != 1 && truthdata->MCJet_Dtr_IsBaryon[h2_index] != 1)
                                        continue;

                                h2_4vector->SetPxPyPzE(truthdata->MCJet_Dtr_PX[h2_index]/1000.,
                                                       truthdata->MCJet_Dtr_PY[h2_index]/1000.,
                                                       truthdata->MCJet_Dtr_PZ[h2_index]/1000., 
                                                       truthdata->MCJet_Dtr_E[h2_index]/1000.);

                                if (!apply_chargedtrack_momentum_cuts(truthdata->MCJet_Dtr_ThreeCharge[h2_index],
                                                                      h2_4vector->P(),
                                                                      h2_4vector->Pt(),
                                                                      h2_4vector->Eta()))
                                        continue;

                                double momentum_weight = weight(h1_4vector->Pt(), h2_4vector->Pt(), Jet_4vector->Pt());

                                h_eec_truth->Fill(h1_4vector->DeltaR(*h2_4vector), Jet_4vector->Pt(), momentum_weight);
                                h_npair_truth->Fill(h1_4vector->DeltaR(*h2_4vector), Jet_4vector->Pt());
                        }   
                }
                h_njet_truth->Fill(Jet_4vector->Pt());
                
                last_eventNum = truthdata->eventNumber;
        }

        fout->cd();
        hpurity_jet->Write();
        hefficiency_jet->Write();
        hpurity->Write();
        hefficiency->Write();
        h_eec->Write();
        h_njet->Write();
        h_njet_wmuoneff->Write();
        h_eec_truth->Write();
        h_njet_truth->Write();
        h_npair->Write();
        h_npair_truth->Write();
        fout->Close();
        
        std::cout<<std::endl;

        return 0;
}
