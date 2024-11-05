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
    TFile* fout = new TFile((output_folder+namef_ntuple_e2c_efficiency).c_str(),"RECREATE");
    
    // Declare the TTrees to be used to build the ntuples
    TZJetsMCReco* mcrecotree = new TZJetsMCReco();

    // Create Ntuples
    TNtuple* ntuple_reco = new TNtuple(name_ntuple_efficiency_reco.c_str(),"",ntuple_efficiency_reco_vars); 
    TNtuple* ntuple_mc   = new TNtuple(name_ntuple_efficiency_mc.c_str()  ,"",ntuple_efficiency_mc_vars); 
    
    ntuple_reco->SetAutoSave(0);
    ntuple_mc->SetAutoSave(0);
    
    // Create necessary 4vectors
    TLorentzVector* Jet_4vector = new TLorentzVector();
    TLorentzVector* Z0_4vector  = new TLorentzVector();
    TLorentzVector* mum_4vector = new TLorentzVector();
    TLorentzVector* mup_4vector = new TLorentzVector();

    // Define array carrying the variables
    float vars_reco[Nvars_efficiency_reco];
    float vars_mc[Nvars_efficiency_mc];

    // Fill the matched jets Ntuple
    for(int evt = 0 ; evt < mcrecotree->fChain->GetEntries() ; evt++)
    {
        // Access entry of tree
        mcrecotree->GetEntry(evt);

        // Apply PV cut
        if(mcrecotree->nPV!=1) continue;

        // -999 means there is not matched jet
        if(mcrecotree->Jet_mcjet_nmcdtrs==-999) continue;

        // Apply trigger cut
        if(mcrecotree->mum_L0MuonEWDecision_TOS!=1||mcrecotree->mum_Hlt1SingleMuonHighPTDecision_TOS!=1||mcrecotree->mum_Hlt2EWSingleMuonVHighPtDecision_TOS!=1) continue;
        if(mcrecotree->mup_L0MuonEWDecision_TOS!=1||mcrecotree->mup_Hlt1SingleMuonHighPTDecision_TOS!=1||mcrecotree->mup_Hlt2EWSingleMuonVHighPtDecision_TOS!=1) continue;

        // Set Jet-associated 4 vectors and apply cuts
        Jet_4vector->SetPxPyPzE(mcrecotree->Jet_PX/1000.,
                                mcrecotree->Jet_PY/1000., 
                                mcrecotree->Jet_PZ/1000., 
                                mcrecotree->Jet_PE/1000.);

        if(Jet_4vector->Eta()<jet_eta_min||Jet_4vector->Eta()>jet_eta_max) continue;
        if(Jet_4vector->Pt()<jet_pt_min) continue;

        mum_4vector->SetPxPyPzE(mcrecotree->mum_PX/1000., 
                                mcrecotree->mum_PY/1000., 
                                mcrecotree->mum_PZ/1000., 
                                mcrecotree->mum_PE/1000.);

        if(Jet_4vector->DeltaR(*mum_4vector, 1/*use rapidity*/)<deltar_muon_jet_min) continue; 
        if(mum_4vector->Pt()<muon_pt_min) continue;
        if(mum_4vector->Eta()<muon_eta_min||mum_4vector->Eta()>muon_eta_max) continue;
        if(mcrecotree->mum_TRACK_PCHI2<muon_trackprob_min) continue;

        mup_4vector->SetPxPyPzE(mcrecotree->mup_PX/1000., 
                                mcrecotree->mup_PY/1000., 
                                mcrecotree->mup_PZ/1000., 
                                mcrecotree->mup_PE/1000.);

        if(Jet_4vector->DeltaR(*mup_4vector, 1/*use rapidity*/)<deltar_muon_jet_min) continue; 
        if(mup_4vector->Pt()<muon_pt_min) continue;
        if(mup_4vector->Eta()<muon_eta_min||mup_4vector->Eta()>muon_eta_max) continue;
        if(mcrecotree->mup_TRACK_PCHI2<muon_trackprob_min) continue;

        Z0_4vector->SetPxPyPzE(mup_4vector->Px()+mum_4vector->Px(), 
                               mup_4vector->Py()+mum_4vector->Py(), 
                               mup_4vector->Pz()+mum_4vector->Pz(), 
                               mup_4vector->E() +mum_4vector->E());

        if(TMath::Abs(Jet_4vector->DeltaPhi(*Z0_4vector))<deltaphi_z_jet_min) continue;
        if(Z0_4vector->M()<muonmuon_mass_min||Z0_4vector->M()>muonmuon_mass_max) continue;

        // Loop over hadron 1
        for(int h1_index = 0 ; h1_index < mcrecotree->Jet_NDtr ; h1_index++)
        {
            // Skip non-hadronic particles
            if(mcrecotree->Jet_Dtr_IsMeson[h1_index]!=1&&mcrecotree->Jet_Dtr_IsBaryon[h1_index]!=1) continue;

            // Skip neutrals
            if(mcrecotree->Jet_Dtr_ThreeCharge[h1_index]==0) continue;

            // Apply cuts
            if(mcrecotree->Jet_Dtr_P[h1_index]/1000.<track_p_min||mcrecotree->Jet_Dtr_P[h1_index]/1000.>track_p_max) continue;
            if(mcrecotree->Jet_Dtr_PT[h1_index]/1000.<track_pt_min) continue;
            if(mcrecotree->Jet_Dtr_TrackChi2[h1_index]/mcrecotree->Jet_Dtr_TrackNDF[h1_index]>track_chi2ndf_max) continue;
            if(mcrecotree->Jet_Dtr_ProbNNghost[h1_index]>track_probnnghost_max) continue;

            // Loop over hadron 2
            for(int h2_index = h1_index+1 ; h2_index < mcrecotree->Jet_NDtr ; h2_index++)
            {
                // Skip non-hadronic particles
                if(mcrecotree->Jet_Dtr_IsMeson[h2_index]!=1&&mcrecotree->Jet_Dtr_IsBaryon[h2_index]!=1) continue;

                // Skip neutrals
                if(mcrecotree->Jet_Dtr_ThreeCharge[h2_index]==0) continue;

                // Apply cuts
                if(mcrecotree->Jet_Dtr_P[h2_index]/1000.<track_p_min||mcrecotree->Jet_Dtr_P[h2_index]/1000.>track_p_max) continue;
                if(mcrecotree->Jet_Dtr_PT[h2_index]/1000.<track_pt_min) continue;
                if(mcrecotree->Jet_Dtr_TrackChi2[h2_index]/mcrecotree->Jet_Dtr_TrackNDF[h2_index]>track_chi2ndf_max) continue;
                if(mcrecotree->Jet_Dtr_ProbNNghost[h2_index]>track_probnnghost_max) continue;

                double h1_y = rapidity(mcrecotree->Jet_Dtr_E[h1_index],mcrecotree->Jet_Dtr_PZ[h1_index]);                
                double h2_y = rapidity(mcrecotree->Jet_Dtr_E[h2_index],mcrecotree->Jet_Dtr_PZ[h2_index]);

                // If all good, fill Ntuple
                vars_reco[0]  = weight(mcrecotree->Jet_Dtr_E[h1_index]/1000., mcrecotree->Jet_Dtr_E[h2_index]/1000., mcrecotree->Jet_PE/1000.);
                vars_reco[1]  = R_L(h1_y, h2_y, mcrecotree->Jet_Dtr_PHI[h1_index], mcrecotree->Jet_Dtr_PHI[h2_index]);
                vars_reco[2]  = mcrecotree->Jet_Dtr_ETA[h1_index];
                vars_reco[3]  = mcrecotree->Jet_Dtr_ETA[h2_index];
                vars_reco[4]  = h1_y;
                vars_reco[5]  = h2_y;
                vars_reco[6]  = mcrecotree->Jet_Dtr_PHI[h1_index];
                vars_reco[7]  = mcrecotree->Jet_Dtr_PHI[h2_index];
                vars_reco[8] = mcrecotree->Jet_Dtr_P[h1_index]/1000.;
                vars_reco[9] = mcrecotree->Jet_Dtr_P[h2_index]/1000.;
                vars_reco[10] = mcrecotree->Jet_Dtr_PT[h1_index]/1000.;
                vars_reco[11] = mcrecotree->Jet_Dtr_PT[h2_index]/1000.;
                vars_reco[12] = mcrecotree->Jet_PT/1000.;
                vars_reco[13] = Jet_4vector->Eta();
                vars_reco[14] = Jet_4vector->DeltaPhi(*Z0_4vector);//Jet_4vector->Phi();
                vars_reco[15] = delta_phi(Jet_4vector->Phi(),Z0_4vector->Phi());
                vars_reco[16] = Jet_4vector->DeltaR(*mum_4vector, 1);
                vars_reco[17] = mum_4vector->Pt();
                vars_reco[18] = mum_4vector->Eta();
                vars_reco[19] = Jet_4vector->DeltaR(*mup_4vector, 1);
                vars_reco[20] = mup_4vector->Pt();
                vars_reco[21] = mup_4vector->Eta();
                vars_reco[22] = mcrecotree->Jet_PE/1000.;
                vars_reco[23] = mcrecotree->Jet_mcjet_PE/1000.;
                vars_reco[24] = mcrecotree->Jet_mcjet_nmcdtrs;

                double matchedmc_y1   = rapidity(mcrecotree->Jet_Dtr_TRUE_E[h1_index],mcrecotree->Jet_Dtr_TRUE_PZ[h1_index]);
                double matchedmc_y2   = rapidity(mcrecotree->Jet_Dtr_TRUE_E[h2_index],mcrecotree->Jet_Dtr_TRUE_PZ[h2_index]);
                double matchedmc_phi1 = mcrecotree->Jet_Dtr_TRUE_PHI[h1_index];
                double matchedmc_phi2 = mcrecotree->Jet_Dtr_TRUE_PHI[h2_index];

                vars_reco[25] = (mcrecotree->Jet_Dtr_TRUE_ETA[h1_index]==-999) ? -999 : R_L(h1_y, matchedmc_y1, mcrecotree->Jet_Dtr_PHI[h1_index], matchedmc_phi1);
                vars_reco[26] = (mcrecotree->Jet_Dtr_TRUE_ETA[h2_index]==-999) ? -999 : R_L(h2_y, matchedmc_y2, mcrecotree->Jet_Dtr_PHI[h2_index], matchedmc_phi2);
                vars_reco[27] = (mcrecotree->Jet_Dtr_TRUE_ETA[h1_index]==-999) ? -999 : matchedmc_y1;
                vars_reco[28] = (mcrecotree->Jet_Dtr_TRUE_ETA[h2_index]==-999) ? -999 : matchedmc_y2;
                vars_reco[29] = (mcrecotree->Jet_Dtr_TRUE_ETA[h1_index]==-999) ? -999 : matchedmc_phi1;
                vars_reco[30] = (mcrecotree->Jet_Dtr_TRUE_ETA[h2_index]==-999) ? -999 : matchedmc_phi2;
                vars_reco[31] = (mcrecotree->Jet_Dtr_TRUE_ETA[h2_index]==-999||mcrecotree->Jet_Dtr_TRUE_ETA[h1_index]==-999) ? 
                            -999 : R_L(matchedmc_y1,matchedmc_y2,matchedmc_phi1,matchedmc_phi2);

                // Fill the TNtuple
                ntuple_reco->Fill(vars_reco);
            }
        }        
    }

    for(int evt = 0 ; evt < mcrecotree->fChain->GetEntries() ; evt++)
    {
        // Access entry of tree
        mcrecotree->GetEntry(evt);

        // Apply PV cut
        if(mcrecotree->nPV!=1) continue;

        // Set Jet-associated 4 vectors and apply cuts
        Jet_4vector->SetPxPyPzE(mcrecotree->Jet_mcjet_PX/1000.,
                                mcrecotree->Jet_mcjet_PY/1000., 
                                mcrecotree->Jet_mcjet_PZ/1000., 
                                mcrecotree->Jet_mcjet_PE/1000.);

        if(Jet_4vector->Eta()<jet_eta_min||Jet_4vector->Eta()>jet_eta_max) continue;
        if(Jet_4vector->Pt()<jet_pt_min) continue;

        mum_4vector->SetPxPyPzE(mcrecotree->Jet_mcjet_mum_PX/1000., 
                                mcrecotree->Jet_mcjet_mum_PY/1000., 
                                mcrecotree->Jet_mcjet_mum_PZ/1000., 
                                mcrecotree->Jet_mcjet_mum_PE/1000.);

        if(Jet_4vector->DeltaR(*mum_4vector, 1/*use rapidity*/)<deltar_muon_jet_min) continue; 
        if(mum_4vector->Pt()<muon_pt_min) continue;
        if(mum_4vector->Eta()<muon_eta_min||mum_4vector->Eta()<muon_eta_max) continue;
        
        mup_4vector->SetPxPyPzE(mcrecotree->Jet_mcjet_mup_PX/1000., 
                                mcrecotree->Jet_mcjet_mup_PY/1000., 
                                mcrecotree->Jet_mcjet_mup_PZ/1000., 
                                mcrecotree->Jet_mcjet_mup_PE/1000.);

        if(Jet_4vector->DeltaR(*mup_4vector, 1/*use rapidity*/)<deltar_muon_jet_min) continue; 
        if(mup_4vector->Pt()<muon_pt_min) continue;
        if(mup_4vector->Eta()<muon_eta_min||mup_4vector->Eta()<muon_eta_max) continue;
        
        Z0_4vector->SetPxPyPzE(mup_4vector->Px()+mum_4vector->Px(), 
                               mup_4vector->Py()+mum_4vector->Py(), 
                               mup_4vector->Pz()+mum_4vector->Pz(), 
                               mup_4vector->E() +mum_4vector->E());

        if(TMath::Abs(Jet_4vector->DeltaPhi(*Z0_4vector))<deltaphi_z_jet_min) continue;
        if(Z0_4vector->M()<muonmuon_mass_min||Z0_4vector->M()>muonmuon_mass_max) continue;

        for(int h1_index = 0 ; h1_index < mcrecotree->Jet_mcjet_nmcdtrs ; h1_index++)
        {
            // Skip non-hadronic particles
            if(mcrecotree->Jet_mcjet_dtrIsMeson[h1_index]!=1&&mcrecotree->Jet_mcjet_dtrIsBaryon[h1_index]!=1) continue;

            // Skip neutrals
            if(mcrecotree->Jet_mcjet_dtrThreeCharge[h1_index]==0) continue;

            // Apply cuts
            if(mcrecotree->Jet_mcjet_dtrP[h1_index]/1000.<track_p_min||mcrecotree->Jet_mcjet_dtrP[h1_index]/1000.>track_p_max) continue;
            if(mcrecotree->Jet_mcjet_dtrPT[h1_index]/1000.<track_pt_min) continue;

            for(int h2_index = h1_index+1 ; h2_index < mcrecotree->Jet_mcjet_nmcdtrs ; h2_index++)
            {
                // Skip non-hadronic particles
                if(mcrecotree->Jet_mcjet_dtrIsMeson[h2_index]!=1&&mcrecotree->Jet_mcjet_dtrIsBaryon[h2_index]!=1) continue;

                // Skip neutrals
                if(mcrecotree->Jet_mcjet_dtrThreeCharge[h2_index]==0) continue;

                // Apply cuts
                if(mcrecotree->Jet_mcjet_dtrP[h2_index]/1000.<track_p_min||mcrecotree->Jet_mcjet_dtrP[h2_index]/1000.>track_p_max) continue;
                if(mcrecotree->Jet_mcjet_dtrPT[h2_index]/1000.<track_pt_min) continue;

                double h1_y = rapidity(mcrecotree->Jet_mcjet_dtrE[h1_index],mcrecotree->Jet_mcjet_dtrPZ[h1_index]); 
                double h2_y = rapidity(mcrecotree->Jet_mcjet_dtrE[h2_index],mcrecotree->Jet_mcjet_dtrPZ[h2_index]);

                vars_mc[0]  = weight(mcrecotree->Jet_mcjet_dtrE[h1_index]/1000.,mcrecotree->Jet_mcjet_dtrE[h2_index]/1000.,mcrecotree->Jet_mcjet_PE/1000.);
                vars_mc[1]  = R_L(h1_y, h2_y, mcrecotree->Jet_mcjet_dtrPHI[h1_index], mcrecotree->Jet_mcjet_dtrPHI[h2_index]);
                vars_mc[2]  = mcrecotree->Jet_mcjet_dtrETA[h1_index];
                vars_mc[3]  = mcrecotree->Jet_mcjet_dtrETA[h2_index];
                vars_mc[4]  = h1_y;
                vars_mc[5]  = h2_y;
                vars_mc[6]  = mcrecotree->Jet_mcjet_dtrPHI[h1_index];
                vars_mc[7]  = mcrecotree->Jet_mcjet_dtrPHI[h2_index];
                vars_mc[8]  = mcrecotree->Jet_mcjet_dtrP[h1_index]/1000.;
                vars_mc[9]  = mcrecotree->Jet_mcjet_dtrP[h2_index]/1000.;
                vars_mc[10] = mcrecotree->Jet_mcjet_dtrPT[h1_index]/1000.;
                vars_mc[11] = mcrecotree->Jet_mcjet_dtrPT[h2_index]/1000.; 
                vars_mc[12] = mcrecotree->Jet_mcjet_PT/1000.;
                vars_mc[13] = mcrecotree->Jet_mcjet_ETA;
                vars_mc[14] = mcrecotree->Jet_mcjet_PHI;
                vars_mc[15] = delta_phi(mcrecotree->Jet_mcjet_PHI,Z0_4vector->Phi());
                vars_mc[16] = R_L(Jet_4vector->Rapidity(),mum_4vector->Rapidity(),mcrecotree->Jet_mcjet_PHI,mcrecotree->Jet_mcjet_mum_PHI);
                vars_mc[17] = mcrecotree->Jet_mcjet_mum_PT/1000.;
                vars_mc[18] = mcrecotree->Jet_mcjet_mum_ETA;
                vars_mc[19] = R_L(Jet_4vector->Rapidity(),mup_4vector->Rapidity(),mcrecotree->Jet_mcjet_PHI,mcrecotree->Jet_mcjet_mup_PHI);
                vars_mc[20] = mcrecotree->Jet_mcjet_mup_PT/1000.;
                vars_mc[21] = mcrecotree->Jet_mcjet_mup_ETA;
                
                // Fill the TNtuple
                ntuple_mc->Fill(vars_mc);        
            }   
        }
    }

    fout->cd();
    ntuple_reco->Write();
    ntuple_mc->Write();
    fout->Close();

    return 0;
}
