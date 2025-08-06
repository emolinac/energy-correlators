#include <iostream>
#include "TZJetsMC.h"
#include "TZJetsMC.C"
#include "TZJetsMCReco.h"
#include "TZJetsMCReco.C"
#include "TZJets2016Data.h"
#include "TZJets2016Data.C"
#include "TZJets2017Data.h"
#include "TZJets2017Data.C"
#include "TZJets2018Data.h"
#include "TZJets2018Data.C"
#include "TROOT.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TH3.h"
#include "analysis-constants.h"
#include "analysis-binning.h"
#include "analysis-cuts.h"
#include "analysis-functions.h"
#include "directories.h"
#include "names.h"
#include "utils-algorithms.h"

int main()
{
        // Open correction files
        TFile* fcorrections    = new TFile((output_folder + namef_ntuple_e2c_hadroncorrections).c_str());
        TFile* fpurity_jet     = new TFile((output_folder + namef_ntuple_jet_purity).c_str());
        TFile* fefficiency_jet = new TFile((output_folder + namef_ntuple_jet_efficiency).c_str());
        
        TFile* fefficiency_muon_2016_id  = new TFile((muons_folder + "IDEff_Data_2016.root").c_str());
        TFile* fefficiency_muon_2016_trk = new TFile((muons_folder + "TRKEff_Data_2016.root").c_str());
        TFile* fefficiency_muon_2016_trg = new TFile((muons_folder + "TRGEff_Data_2016.root").c_str());
        
        // Create output file
        TFile* fout = new TFile((output_folder + "ntuple_mcreco_hadronmultiplicity.root").c_str(),"RECREATE");
        
        // Declare the TTrees to be used to build the ntuples
        TZJetsMCReco* mcrecotree = new TZJetsMCReco();
        
        // Create Ntuples
        TNtuple* ntuple_purity          = (TNtuple*) fcorrections->Get((name_ntuple_correction_reco.c_str()));
        TNtuple* ntuple_efficiency_mc   = (TNtuple*) fcorrections->Get((name_ntuple_correction_mc.c_str()));
        TNtuple* ntuple_efficiency_reco = (TNtuple*) fcorrections->Get((name_ntuple_correction_reco.c_str()));
        
        TNtuple* ntuple_purity_jet      = (TNtuple*) fpurity_jet->Get((name_ntuple_jetpurity.c_str()));
        TNtuple* ntuple_efficiency_jet  = (TNtuple*) fefficiency_jet->Get((name_ntuple_jetefficiency.c_str()));

        TNtuple* ntuple_data            = new TNtuple(name_ntuple_data.c_str(),"All Data",ntuple_corrdata_vars); 
        TNtuple* ntuple_corrjet         = new TNtuple(name_ntuple_corrjet.c_str(),"All Data",ntuple_jet_vars); 
        
        ntuple_data->SetAutoSave(0);
        ntuple_corrjet->SetAutoSave(0);

        // Muon corrections
        TH2D* h2_muon_2016_ideff_data  = (TH2D*) fefficiency_muon_2016_id->Get("Hist_ALL_2016_ETA_PT_Eff");
        TH2D* h2_muon_2016_trkeff_data = (TH2D*) fefficiency_muon_2016_trk->Get("Hist_ALL_2016_ETA_PT_Eff");
        TH2D* h2_muon_2016_trgeff_data = (TH2D*) fefficiency_muon_2016_trg->Get("Hist_ALL_2016_ETA_PT_Eff");

        // Jet corrections
        TH1F* hnum_pur_jet = new TH1F("hnum_pur_jet", "", nbin_jet_pt_corrections, jet_pt_corrections_binning);
        TH1F* hden_pur_jet = new TH1F("hden_pur_jet", "", nbin_jet_pt_corrections, jet_pt_corrections_binning);
        TH1F* hpurity_jet  = new TH1F("hpurity_jet" , "", nbin_jet_pt_corrections, jet_pt_corrections_binning);
        hnum_pur_jet->Sumw2();
        hden_pur_jet->Sumw2();

        TH1F* hnum_eff_jet    = new TH1F("hnum_eff_jet"   , "", nbin_jet_pt_corrections, jet_pt_corrections_binning);
        TH1F* hden_eff_jet    = new TH1F("hden_eff_jet"   , "", nbin_jet_pt_corrections, jet_pt_corrections_binning);
        TH1F* hefficiency_jet = new TH1F("hefficiency_jet", "", nbin_jet_pt_corrections, jet_pt_corrections_binning);
        hnum_eff_jet->Sumw2();
        hden_eff_jet->Sumw2();

        ntuple_purity_jet->Project("hnum_pur_jet", "jet_pt", "jet_pt_truth!=-999");
        ntuple_purity_jet->Project("hden_pur_jet", "jet_pt");
        ntuple_efficiency_jet->Project("hnum_eff_jet", "jet_pt_truth", "jet_pt!=-999");
        ntuple_efficiency_jet->Project("hden_eff_jet", "jet_pt_truth");

        hpurity_jet->Divide(hnum_pur_jet, hden_pur_jet, 1, 1, "B");
        hefficiency_jet->Divide(hnum_eff_jet, hden_eff_jet, 1, 1, "B");

        // Hadron corrections
        TH3F* hnum_pur    = new TH3F("hnum_pur"   , "", ic_p_nbins, ic_p_binning, sl_eta_nbins, sl_eta_binning, nbin_jet_pt_corrections, jet_pt_corrections_binning);
        TH3F* hden_pur    = new TH3F("hden_pur"   , "", ic_p_nbins, ic_p_binning, sl_eta_nbins, sl_eta_binning, nbin_jet_pt_corrections, jet_pt_corrections_binning);
        TH3F* hpurity     = new TH3F("hpurity"    , "", ic_p_nbins, ic_p_binning, sl_eta_nbins, sl_eta_binning, nbin_jet_pt_corrections, jet_pt_corrections_binning);
        TH3F* hnum_eff    = new TH3F("hnum_eff"   , "", ic_p_nbins, ic_p_binning, sl_eta_nbins, sl_eta_binning, nbin_jet_pt_corrections, jet_pt_corrections_binning);
        TH3F* hden_eff    = new TH3F("hden_eff"   , "", ic_p_nbins, ic_p_binning, sl_eta_nbins, sl_eta_binning, nbin_jet_pt_corrections, jet_pt_corrections_binning);
        TH3F* hefficiency = new TH3F("hefficiency", "", ic_p_nbins, ic_p_binning, sl_eta_nbins, sl_eta_binning, nbin_jet_pt_corrections, jet_pt_corrections_binning);
        
        hnum_pur->Sumw2();
        hden_pur->Sumw2();
        hnum_eff->Sumw2();
        hden_eff->Sumw2();
        
        ntuple_purity->Project("hnum_pur","jet_pt:h_eta:h_p",single_signal_cut);
        ntuple_purity->Project("hden_pur","jet_pt:h_eta:h_p");
        ntuple_efficiency_reco->Project("hnum_eff","jet_pt_truth:h_eta_truth:h_p_truth",single_signal_cut);
        ntuple_efficiency_mc->Project("hden_eff","jet_pt:h_eta:h_p");

        hpurity->Divide(hnum_pur, hden_pur, 1, 1, "B");
        hefficiency->Divide(hnum_eff, hden_eff, 1, 1, "B");

        regularize_correction_factors(hpurity);
        regularize_correction_factors(hefficiency);
        
        TCanvas* c = new TCanvas("c", "", 1920, 1080);
        c->Draw();

        // Create necessary 4vectors
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
        
        // Define array carrying the variables
        float vars[Nvars_corrdata];
        float vars_jet[Nvars_corrjet];

        // Fill the data TNtuple
        std::cout<<"Working with MCReco."<<std::endl;
        for (int evt = 0 ; evt < mcrecotree->fChain->GetEntries() ; evt++) {
                // Access entry of tree
                mcrecotree->GetEntry(evt);

                if (evt%10000 == 0) {
                        double percentage = 100.*evt/mcrecotree->fChain->GetEntries();
                        std::cout<<"\r"<<percentage<<"\% jets processed."<< std::flush;
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

                double mpt = 0; // Min pt at least one track has to have
                for (int h_index = 0 ; h_index < mcrecotree->Jet_NDtr ; h_index++) {
                        // Skip non-hadronic particles
                        if (mcrecotree->Jet_Dtr_IsMeson[h_index] != 1&&mcrecotree->Jet_Dtr_IsBaryon[h_index] != 1) 
                                continue;

                        h1_4vector->SetPxPyPzE(mcrecotree->Jet_Dtr_PX[h_index]/1000.,
                                              mcrecotree->Jet_Dtr_PY[h_index]/1000.,
                                              mcrecotree->Jet_Dtr_PZ[h_index]/1000.,
                                              mcrecotree->Jet_Dtr_E[h_index]/1000.);

                        if (h1_4vector->Pt() > mpt)
                                mpt = h1_4vector->Pt();
                }

                if (!apply_jet_id_cuts(mpt, mcrecotree->Jet_NTrk, mcrecotree->Jet_CPF, mcrecotree->Jet_MTF))
                        continue;

                double mup_pt  = (mup_4vector->Pt() >= 70.) ? 69. : mup_4vector->Pt();
                double mum_pt  = (mum_4vector->Pt() >= 70.) ? 69. : mum_4vector->Pt();
                double mup_eta = mup_4vector->Eta();
                double mum_eta = mum_4vector->Eta();

                double mup_eff_id  = h2_muon_2016_ideff_data->GetBinContent(h2_muon_2016_ideff_data->FindBin(mup_eta, mup_pt));
                double mup_eff_trk = h2_muon_2016_trkeff_data->GetBinContent(h2_muon_2016_trkeff_data->FindBin(mup_eta, mup_pt));
                double mup_eff_trg = h2_muon_2016_trgeff_data->GetBinContent(h2_muon_2016_trgeff_data->FindBin(mup_eta, mup_pt));

                double mum_eff_id  = h2_muon_2016_ideff_data->GetBinContent(h2_muon_2016_ideff_data->FindBin(mum_eta, mum_pt));
                double mum_eff_trk = h2_muon_2016_trkeff_data->GetBinContent(h2_muon_2016_trkeff_data->FindBin(mum_eta, mum_pt));
                double mum_eff_trg = h2_muon_2016_trgeff_data->GetBinContent(h2_muon_2016_trgeff_data->FindBin(mum_eta, mum_pt));

                double jet_efficiency = hefficiency_jet->GetBinContent(hefficiency_jet->FindBin(Jet_4vector->Pt()));
                double jet_purity     = hpurity_jet->GetBinContent(hpurity_jet->FindBin(Jet_4vector->Pt()));
                
                double jet_efficiency_error = hefficiency_jet->GetBinError(hefficiency_jet->FindBin(Jet_4vector->Pt()));
                double jet_purity_error     = hpurity_jet->GetBinError(hpurity_jet->FindBin(Jet_4vector->Pt()));
                
                double jet_ndtr_nominal = 0;
                for (int h1_index = 0 ; h1_index < mcrecotree->Jet_NDtr ; h1_index++) {
                        // Skip non-hadronic particles
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

                        jet_ndtr_nominal++;
                }
                if (jet_ndtr_nominal < 1)
                        continue;

                // Loop over hadron 1
                for (int h1_index = 0 ; h1_index < mcrecotree->Jet_NDtr ; h1_index++) {
                        // Skip non-hadronic particles
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

                        double h1_purity     = hpurity->GetBinContent(hpurity->FindBin(h1_4vector->P(),h1_4vector->Eta(), Jet_4vector->Pt()));
                        double h1_efficiency = hefficiency->GetBinContent(hefficiency->FindBin(h1_4vector->P(),h1_4vector->Eta(), Jet_4vector->Pt()));
                        if (h1_purity > 1. || h1_efficiency > 1.) 
                                continue;

                        double weight_due_to_jet = jet_purity/jet_efficiency/(mum_eff_id*mup_eff_id*mum_eff_trk*mup_eff_trk*(mum_eff_trg+mup_eff_trg-mum_eff_trg*mup_eff_trg));
                        // double weight_due_to_jet = jet_purity/jet_efficiency;

                        vars[0 ] = weight_due_to_jet;
                        vars[1 ] = h1_efficiency;
                        vars[2 ] = h1_purity;
                        vars[3 ] = -999;
                        vars[4 ] = -999;
                        vars[5 ] = -999;
                        vars[6 ] = h1_4vector->Eta();
                        vars[7 ] = -999;
                        vars[8 ] = h1_4vector->Rapidity();
                        vars[9 ] = -999;
                        vars[10] = h1_4vector->P();
                        vars[11] = -999;
                        vars[12] = h1_4vector->Pt();
                        vars[13] = -999;
                        vars[14] = Jet_4vector->Pt();
                        vars[15] = Jet_4vector->Eta();
                        vars[16] = -999;
                        vars[17] = Jet_4vector->E();
                        vars[18] = h1_4vector->E();
                        vars[19] = -999;
                        vars[20] = -999;
                        vars[21] = -999;
                        vars[22] = -999;
                        vars[23] = -999;
                        vars[24] = -999;
                        vars[25] = -999;
                        vars[26] = -999;
                        vars[27] = -999;
                        vars[28] = -999;
                        
                        // Fill the TNtuple
                        ntuple_data->Fill(vars);
                }

                vars_jet[0]  = Jet_4vector->Pt();
                vars_jet[1]  = Jet_4vector->E();
                vars_jet[2]  = mcrecotree->Jet_NDtr;
                vars_jet[3]  = jet_efficiency;
                vars_jet[4]  = jet_purity;
                vars_jet[5]  = jet_efficiency_error;
                vars_jet[6]  = jet_purity_error;
                vars_jet[7]  = mup_eff_id;
                vars_jet[8]  = mup_eff_trk;
                vars_jet[9]  = mup_eff_trg;
                vars_jet[10] = mum_eff_id;
                vars_jet[11] = mum_eff_trk;
                vars_jet[12] = mum_eff_trg;
                vars_jet[13] = 2016;
                vars_jet[14] = Jet_4vector->Eta();
                vars_jet[15] = Z0_4vector->Pt();
                vars_jet[16] = Z0_4vector->Eta();
                vars_jet[17] = Z0_4vector->Rapidity();
                
                ntuple_corrjet->Fill(vars_jet);

                last_eventNum = mcrecotree->eventNumber;
        }


        fout->cd();
        ntuple_data->Write();
        ntuple_corrjet->Write();
        fout->Close();
        
        std::cout<<std::endl;

        return 0;
}
