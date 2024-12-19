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
    // Open correction files
    TFile* fpurity     = new TFile((output_folder+namef_ntuple_e2c_purity).c_str());
    TFile* fefficiency = new TFile((output_folder+namef_ntuple_e2c_efficiency).c_str());
    
    // Create output file
    TFile* fout = new TFile((output_folder+namef_ntuple_e2c_corr).c_str(),"RECREATE");
    
    // Declare the TTrees to be used to build the ntuples
    TZJetsData* datatree = new TZJetsData();
    
    // Create Ntuples
    TNtuple* ntuple_purity          = (TNtuple*) fpurity->Get((name_ntuple_purity.c_str()));
    TNtuple* ntuple_efficiency_mc   = (TNtuple*) fefficiency->Get((name_ntuple_efficiency_mc.c_str()));
    TNtuple* ntuple_efficiency_reco = (TNtuple*) fefficiency->Get((name_ntuple_efficiency_reco.c_str()));
    TNtuple* ntuple_data            = new TNtuple(name_ntuple_data.c_str(),"All Data",ntuple_corrdata_vars); 
    ntuple_data->SetAutoSave(0);

    // Calculate corrections
    TH2F* hsigp   = new TH2F("hsigp"  ,"",ndim_corr,-TMath::Pi(),TMath::Pi(),ndim_corr,eta_min,eta_max);
    TH2F* hallp   = new TH2F("hallp"  ,"",ndim_corr,-TMath::Pi(),TMath::Pi(),ndim_corr,eta_min,eta_max);
    TH2F* hpurity = new TH2F("hpurity","",ndim_corr,-TMath::Pi(),TMath::Pi(),ndim_corr,eta_min,eta_max);
    hsigp->Sumw2();
    hallp->Sumw2();
    hpurity->Sumw2();

    TH2F* hsigeff     = new TH2F("hsigeff"    ,"",ndim_corr,-TMath::Pi(),TMath::Pi(),ndim_corr,eta_min,eta_max);
    TH2F* halleff     = new TH2F("halleff"    ,"",ndim_corr,-TMath::Pi(),TMath::Pi(),ndim_corr,eta_min,eta_max);
    TH2F* hefficiency = new TH2F("hefficiency","",ndim_corr,-TMath::Pi(),TMath::Pi(),ndim_corr,eta_min,eta_max);
    hsigeff->Sumw2();
    halleff->Sumw2();
    hefficiency->Sumw2();

    ntuple_purity->Project("hsigp","h_eta:h_phi",single_signal_cut);
    ntuple_purity->Project("hallp","h_eta:h_phi",pair_cut         );
    hpurity->Divide(hsigp,hallp,1,1,"B");

    ntuple_efficiency_reco->Project("hsigeff","h_eta:h_phi",single_signal_cut);
    ntuple_efficiency_mc->Project("halleff","h_eta:h_phi",pair_cut         );
    hefficiency->Divide(hsigeff,halleff,1,1,"B");

    // Create necessary 4vectors
    TLorentzVector* Jet_4vector = new TLorentzVector();
    TLorentzVector* Z0_4vector  = new TLorentzVector();
    TLorentzVector* mum_4vector = new TLorentzVector();
    TLorentzVector* mup_4vector = new TLorentzVector();

    // Define array carrying the variables
    float vars[Nvars_corrdata];

    // Fill the data TNtuple
    for(int evt = 0 ; evt < datatree->fChain->GetEntries() ; evt++)
    {
        // Access entry of tree
        datatree->GetEntry(evt);

        // Apply PV cut
        if(datatree->nPV!=1) continue;

        // Apply trigger cut
        bool mum_trigger = (datatree->mum_L0MuonEWDecision_TOS==1&&datatree->mum_Hlt1SingleMuonHighPTDecision_TOS==1&&datatree->mum_Hlt2EWSingleMuonVHighPtDecision_TOS==1);
        bool mup_trigger = (datatree->mup_L0MuonEWDecision_TOS==1&&datatree->mup_Hlt1SingleMuonHighPTDecision_TOS==1&&datatree->mup_Hlt2EWSingleMuonVHighPtDecision_TOS==1);

        if(!mum_trigger&&!mup_trigger) continue;

        // Set Jet-associated 4 vectors and apply cuts
        Jet_4vector->SetPxPyPzE(datatree->Jet_PX/1000.,
                                datatree->Jet_PY/1000., 
                                datatree->Jet_PZ/1000., 
                                datatree->Jet_PE/1000.);

        if(Jet_4vector->Eta()<jet_eta_min||Jet_4vector->Eta()>jet_eta_max) continue;
        if(Jet_4vector->Pt()<jet_pt_min_nom) continue;

        mum_4vector->SetPxPyPzE(datatree->mum_PX/1000., 
                                datatree->mum_PY/1000., 
                                datatree->mum_PZ/1000., 
                                datatree->mum_PE/1000.);

        if(Jet_4vector->DeltaR(*mum_4vector, 0/*use pseudorapidity*/)<deltar_muon_jet_min) continue; 
        if(mum_4vector->Pt()<muon_pt_min) continue;
        if(mum_4vector->Eta()<muon_eta_min||mum_4vector->Eta()>muon_eta_max) continue;
        if(datatree->mum_TRACK_PCHI2<muon_trackprob_min) continue;

        mup_4vector->SetPxPyPzE(datatree->mup_PX/1000., 
                                datatree->mup_PY/1000., 
                                datatree->mup_PZ/1000., 
                                datatree->mup_PE/1000.);

        if(Jet_4vector->DeltaR(*mup_4vector, 0/*use pseudorapidity*/)<deltar_muon_jet_min) continue; 
        if(mup_4vector->Pt()<muon_pt_min) continue;
        if(mup_4vector->Eta()<muon_eta_min||mup_4vector->Eta()>muon_eta_max) continue;
        if(datatree->mup_TRACK_PCHI2<muon_trackprob_min) continue;

        Z0_4vector->SetPxPyPzE(mup_4vector->Px()+mum_4vector->Px(), 
                               mup_4vector->Py()+mum_4vector->Py(), 
                               mup_4vector->Pz()+mum_4vector->Pz(), 
                               mup_4vector->E() +mum_4vector->E());

        if(TMath::Abs(Jet_4vector->DeltaPhi(*Z0_4vector))<deltaphi_z_jet_min) continue;
        if(Z0_4vector->M()<muonmuon_mass_min||Z0_4vector->M()>muonmuon_mass_max) continue;

        // Loop over hadron 1
        for(int h1_index = 0 ; h1_index < datatree->Jet_NDtr ; h1_index++)
        {
            // Skip non-hadronic particles
            if(datatree->Jet_Dtr_IsMeson[h1_index]!=1&&datatree->Jet_Dtr_IsBaryon[h1_index]!=1) continue;

            // Skip neutrals
            if(datatree->Jet_Dtr_ThreeCharge[h1_index]==0) continue;

            // Apply cuts
            if(datatree->Jet_Dtr_P[h1_index]/1000.<track_p_min||datatree->Jet_Dtr_P[h1_index]/1000.>track_p_max) continue;
            if(datatree->Jet_Dtr_PT[h1_index]/1000.<track_pt_min) continue;
            if(datatree->Jet_Dtr_TrackChi2[h1_index]/datatree->Jet_Dtr_TrackNDF[h1_index]>track_chi2ndf_max) continue;
            if(datatree->Jet_Dtr_ProbNNghost[h1_index]>track_probnnghost_max) continue;

            // Loop over hadron 2
            for(int h2_index = h1_index+1 ; h2_index < datatree->Jet_NDtr ; h2_index++)
            {
                // Skip non-hadronic particles
                if(datatree->Jet_Dtr_IsMeson[h2_index]!=1&&datatree->Jet_Dtr_IsBaryon[h2_index]!=1) continue;

                // Skip neutrals
                if(datatree->Jet_Dtr_ThreeCharge[h2_index]==0) continue;

                // Apply cuts
                if(datatree->Jet_Dtr_P[h2_index]/1000.<track_p_min||datatree->Jet_Dtr_P[h2_index]/1000.>track_p_max) continue;
                if(datatree->Jet_Dtr_PT[h2_index]/1000.<track_pt_min) continue;
                if(datatree->Jet_Dtr_TrackChi2[h2_index]/datatree->Jet_Dtr_TrackNDF[h2_index]>track_chi2ndf_max) continue;
                if(datatree->Jet_Dtr_ProbNNghost[h2_index]>track_probnnghost_max) continue;

                double purity_correction = (hpurity->GetBinContent(hpurity->FindBin(datatree->Jet_Dtr_PHI[h1_index],datatree->Jet_Dtr_ETA[h1_index])))*
                                           (hpurity->GetBinContent(hpurity->FindBin(datatree->Jet_Dtr_PHI[h2_index],datatree->Jet_Dtr_ETA[h2_index])));

                double purity_error = (hpurity->GetBinContent(hpurity->FindBin(datatree->Jet_Dtr_PHI[h1_index],datatree->Jet_Dtr_ETA[h1_index])))*
                                      (hpurity->GetBinError(hpurity->FindBin(datatree->Jet_Dtr_PHI[h2_index],datatree->Jet_Dtr_ETA[h2_index]))) +
                                      (hpurity->GetBinError(hpurity->FindBin(datatree->Jet_Dtr_PHI[h1_index],datatree->Jet_Dtr_ETA[h1_index])))*
                                      (hpurity->GetBinContent(hpurity->FindBin(datatree->Jet_Dtr_PHI[h2_index],datatree->Jet_Dtr_ETA[h2_index])));

                double efficiency_correction = (hefficiency->GetBinContent(hefficiency->FindBin(datatree->Jet_Dtr_PHI[h1_index],datatree->Jet_Dtr_ETA[h1_index])))*
                                               (hefficiency->GetBinContent(hefficiency->FindBin(datatree->Jet_Dtr_PHI[h2_index],datatree->Jet_Dtr_ETA[h2_index])));

                double efficiency_error = (hefficiency->GetBinContent(hefficiency->FindBin(datatree->Jet_Dtr_PHI[h1_index],datatree->Jet_Dtr_ETA[h1_index])))*
                                          (hefficiency->GetBinError(hefficiency->FindBin(datatree->Jet_Dtr_PHI[h2_index],datatree->Jet_Dtr_ETA[h2_index]))) +
                                          (hefficiency->GetBinError(hefficiency->FindBin(datatree->Jet_Dtr_PHI[h1_index],datatree->Jet_Dtr_ETA[h1_index])))*
                                          (hefficiency->GetBinContent(hefficiency->FindBin(datatree->Jet_Dtr_PHI[h2_index],datatree->Jet_Dtr_ETA[h2_index])));

                vars[0 ] = weight(datatree->Jet_Dtr_E[h1_index], datatree->Jet_Dtr_E[h2_index], datatree->Jet_PE);
                vars[1 ] = efficiency_correction;
                vars[2 ] = purity_correction;
                vars[3 ] = efficiency_error/efficiency_correction;
                vars[4 ] = purity_error/purity_correction;
                vars[5 ] = R_L(rapidity(datatree->Jet_Dtr_E[h1_index],datatree->Jet_Dtr_PZ[h1_index]), rapidity(datatree->Jet_Dtr_E[h2_index],datatree->Jet_Dtr_PZ[h2_index]),
                                        datatree->Jet_Dtr_PHI[h1_index], datatree->Jet_Dtr_PHI[h2_index]);
                vars[6 ] = datatree->Jet_Dtr_ETA[h1_index];
                vars[7 ] = datatree->Jet_Dtr_ETA[h2_index];
                vars[8 ] = rapidity(datatree->Jet_Dtr_E[h1_index],datatree->Jet_Dtr_PZ[h1_index]);
                vars[9 ] = rapidity(datatree->Jet_Dtr_E[h2_index],datatree->Jet_Dtr_PZ[h2_index]);
                vars[10] = datatree->Jet_Dtr_PHI[h1_index];
                vars[11] = datatree->Jet_Dtr_PHI[h2_index];
                vars[12] = datatree->Jet_Dtr_P[h1_index]/1000.;
                vars[13] = datatree->Jet_Dtr_P[h2_index]/1000.;
                vars[14] = datatree->Jet_Dtr_PT[h1_index]/1000.;
                vars[15] = datatree->Jet_Dtr_PT[h2_index]/1000.;
                vars[16] = datatree->Jet_PT/1000.;
                vars[17] = Jet_4vector->Eta();
                vars[18] = Jet_4vector->Phi();
                vars[19] = delta_phi(Jet_4vector->Phi(),Z0_4vector->Phi());
                vars[20] = R_L(Jet_4vector->Rapidity(),mum_4vector->Rapidity(),Jet_4vector->Phi(),mum_4vector->Phi());
                vars[21] = mum_4vector->Pt();
                vars[22] = mum_4vector->Eta();
                vars[23] = R_L(Jet_4vector->Rapidity(),mup_4vector->Rapidity(),Jet_4vector->Phi(),mup_4vector->Phi());
                vars[24] = mup_4vector->Pt();
                vars[25] = mup_4vector->Eta();
                
                // Fill the TNtuple
                ntuple_data->Fill(vars);
            }
        }
    }

    fout->cd();
    ntuple_data->Write();
    fout->Close();

    return 0;
}
