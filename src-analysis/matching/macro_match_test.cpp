#include <iostream>
#include "../../include/TZJetsMCRecoNoMatchedJetDtrs.h"
#include "../../include/TZJetsMCRecoNoMatchedJetDtrs.C"
#include "../../include/TZJetsMCRecoMatchedJetDtrs.h"
#include "../../include/TZJetsMCRecoMatchedJetDtrs.C"
#include "../../include/analysis-constants.h"
#include "../../include/analysis-cuts.h"
#include "../../include/analysis-functions.h"
#include "../../include/directories.h"
#include "../../include/names.h"

void macro_match_test()
{
    // Declare the TTrees to be used to build the ntuples
    TZJetsMCRecoNoMatchedJetDtrs* tree_nomatchedjetdtrs = new TZJetsMCRecoNoMatchedJetDtrs();
    TZJetsMCRecoMatchedJetDtrs*   tree_matchedjetdtrs   = new TZJetsMCRecoMatchedJetDtrs();

    // Create necessary 4vectors
    TLorentzVector* Jet_4vector = new TLorentzVector();
    TLorentzVector* Z0_4vector  = new TLorentzVector();
    TLorentzVector* mum_4vector = new TLorentzVector();
    TLorentzVector* mup_4vector = new TLorentzVector();

    // Counters
    double N1_total = 0;
    double N1forced_total = 0;
    double N2_total = 0;

    double N1matched_total = 0;
    double N1forcedmatched_total = 0;
    double N2matched_total = 0;
    
    // Fill the matched jets Ntuple
    for(int evt = 0 ; evt < tree_nomatchedjetdtrs->fChain->GetEntries() ; evt++)
    {
        // Access entry of tree
        tree_nomatchedjetdtrs->GetEntry(evt);

        // Apply PV cut
        if(tree_nomatchedjetdtrs->nPV!=1) continue;

        // Apply trigger cut
        if(tree_nomatchedjetdtrs->mum_L0MuonEWDecision_TOS!=1||tree_nomatchedjetdtrs->mum_Hlt1SingleMuonHighPTDecision_TOS!=1||tree_nomatchedjetdtrs->mum_Hlt2EWSingleMuonVHighPtDecision_TOS!=1) continue;
        if(tree_nomatchedjetdtrs->mup_L0MuonEWDecision_TOS!=1||tree_nomatchedjetdtrs->mup_Hlt1SingleMuonHighPTDecision_TOS!=1||tree_nomatchedjetdtrs->mup_Hlt2EWSingleMuonVHighPtDecision_TOS!=1) continue;

        // -999 means there is not matched jet
        if(tree_nomatchedjetdtrs->Jet_mcjet_nmcdtrs==-999) continue;

        // Set Jet-associated 4 vectors and apply cuts
        Jet_4vector->SetPxPyPzE(tree_nomatchedjetdtrs->Jet_PX/1000.,
                                tree_nomatchedjetdtrs->Jet_PY/1000., 
                                tree_nomatchedjetdtrs->Jet_PZ/1000., 
                                tree_nomatchedjetdtrs->Jet_PE/1000.);

        if(Jet_4vector->Eta()<jet_eta_min||Jet_4vector->Eta()>jet_eta_max) continue;
        if(Jet_4vector->Pt()<jet_pt_min_nom) continue;

        mum_4vector->SetPxPyPzE(tree_nomatchedjetdtrs->mum_PX/1000., 
                                tree_nomatchedjetdtrs->mum_PY/1000., 
                                tree_nomatchedjetdtrs->mum_PZ/1000., 
                                tree_nomatchedjetdtrs->mum_PE/1000.);

        if(Jet_4vector->DeltaR(*mum_4vector, 0/*use pseudorapidity*/)<jet_radius) continue; 
        if(mum_4vector->Pt()<muon_pt_min) continue;
        if(mum_4vector->Eta()<muon_eta_min||mum_4vector->Eta()>muon_eta_max) continue;
        if(tree_nomatchedjetdtrs->mum_TRACK_PCHI2<muon_trackprob_min) continue;

        mup_4vector->SetPxPyPzE(tree_nomatchedjetdtrs->mup_PX/1000., 
                                tree_nomatchedjetdtrs->mup_PY/1000., 
                                tree_nomatchedjetdtrs->mup_PZ/1000., 
                                tree_nomatchedjetdtrs->mup_PE/1000.);

        if(Jet_4vector->DeltaR(*mup_4vector, 0/*use pseudorapidity*/)<jet_radius) continue; 
        if(mup_4vector->Pt()<muon_pt_min) continue;
        if(mup_4vector->Eta()<muon_eta_min||mup_4vector->Eta()>muon_eta_max) continue;
        if(tree_nomatchedjetdtrs->mup_TRACK_PCHI2<muon_trackprob_min) continue;

        Z0_4vector->SetPxPyPzE(mup_4vector->Px()+mum_4vector->Px(), 
                               mup_4vector->Py()+mum_4vector->Py(), 
                               mup_4vector->Pz()+mum_4vector->Pz(), 
                               mup_4vector->E() +mum_4vector->E());

        if(TMath::Abs(Jet_4vector->DeltaPhi(*Z0_4vector))<deltaphi_z_jet_min) continue;
        if(Z0_4vector->M()<muonmuon_mass_min||Z0_4vector->M()>muonmuon_mass_max) continue;

        
        // Loop over hadron 2
        for(int h_index = 0 ; h_index < tree_nomatchedjetdtrs->Jet_NDtr ; h_index++)
        {
            // Skip non-hadronic particles
            if(tree_nomatchedjetdtrs->Jet_Dtr_IsMeson[h_index]!=1&&tree_nomatchedjetdtrs->Jet_Dtr_IsBaryon[h_index]!=1) continue;
            // Skip neutrals
            if(tree_nomatchedjetdtrs->Jet_Dtr_ThreeCharge[h_index]==0) continue;
            // Apply cuts
            if(tree_nomatchedjetdtrs->Jet_Dtr_P[h_index]/1000.<track_p_min||tree_nomatchedjetdtrs->Jet_Dtr_P[h_index]/1000.>track_p_max) continue;
            if(tree_nomatchedjetdtrs->Jet_Dtr_PT[h_index]/1000.<track_pt_min) continue;
            if(tree_nomatchedjetdtrs->Jet_Dtr_TrackChi2[h_index]/tree_nomatchedjetdtrs->Jet_Dtr_TrackNDF[h_index]>track_chi2ndf_max) continue;
            if(tree_nomatchedjetdtrs->Jet_Dtr_ProbNNghost[h_index]>track_probnnghost_max) continue;

            // If all good lets count
            N1_total++;
            if(tree_nomatchedjetdtrs->Jet_Dtr_TRUE_ETA[h_index]!=-999) N1matched_total++;
        }
    }

    for(int evt = 0 ; evt < tree_nomatchedjetdtrs->fChain->GetEntries() ; evt++)
    {
        // Access entry of tree
        tree_nomatchedjetdtrs->GetEntry(evt);

        if(tree_nomatchedjetdtrs->Jet_mcjet_dtrE[0]==-999) continue;

        // Apply PV cut
        if(tree_nomatchedjetdtrs->nPV!=1) continue;

        // Apply trigger cut
        if(tree_nomatchedjetdtrs->mum_L0MuonEWDecision_TOS!=1||tree_nomatchedjetdtrs->mum_Hlt1SingleMuonHighPTDecision_TOS!=1||tree_nomatchedjetdtrs->mum_Hlt2EWSingleMuonVHighPtDecision_TOS!=1) continue;
        if(tree_nomatchedjetdtrs->mup_L0MuonEWDecision_TOS!=1||tree_nomatchedjetdtrs->mup_Hlt1SingleMuonHighPTDecision_TOS!=1||tree_nomatchedjetdtrs->mup_Hlt2EWSingleMuonVHighPtDecision_TOS!=1) continue;

        // -999 means there is not matched jet
        if(tree_nomatchedjetdtrs->Jet_mcjet_nmcdtrs==-999) continue;

        // Set Jet-associated 4 vectors and apply cuts
        Jet_4vector->SetPxPyPzE(tree_nomatchedjetdtrs->Jet_PX/1000.,
                                tree_nomatchedjetdtrs->Jet_PY/1000., 
                                tree_nomatchedjetdtrs->Jet_PZ/1000., 
                                tree_nomatchedjetdtrs->Jet_PE/1000.);

        if(Jet_4vector->Eta()<jet_eta_min||Jet_4vector->Eta()>jet_eta_max) continue;
        if(Jet_4vector->Pt()<jet_pt_min_nom) continue;

        mum_4vector->SetPxPyPzE(tree_nomatchedjetdtrs->mum_PX/1000., 
                                tree_nomatchedjetdtrs->mum_PY/1000., 
                                tree_nomatchedjetdtrs->mum_PZ/1000., 
                                tree_nomatchedjetdtrs->mum_PE/1000.);

        if(Jet_4vector->DeltaR(*mum_4vector, 0/*use pseudorapidity*/)<jet_radius) continue; 
        if(mum_4vector->Pt()<muon_pt_min) continue;
        if(mum_4vector->Eta()<muon_eta_min||mum_4vector->Eta()>muon_eta_max) continue;
        if(tree_nomatchedjetdtrs->mum_TRACK_PCHI2<muon_trackprob_min) continue;

        mup_4vector->SetPxPyPzE(tree_nomatchedjetdtrs->mup_PX/1000., 
                                tree_nomatchedjetdtrs->mup_PY/1000., 
                                tree_nomatchedjetdtrs->mup_PZ/1000., 
                                tree_nomatchedjetdtrs->mup_PE/1000.);

        if(Jet_4vector->DeltaR(*mup_4vector, 0/*use pseudorapidity*/)<jet_radius) continue; 
        if(mup_4vector->Pt()<muon_pt_min) continue;
        if(mup_4vector->Eta()<muon_eta_min||mup_4vector->Eta()>muon_eta_max) continue;
        if(tree_nomatchedjetdtrs->mup_TRACK_PCHI2<muon_trackprob_min) continue;

        Z0_4vector->SetPxPyPzE(mup_4vector->Px()+mum_4vector->Px(), 
                               mup_4vector->Py()+mum_4vector->Py(), 
                               mup_4vector->Pz()+mum_4vector->Pz(), 
                               mup_4vector->E() +mum_4vector->E());

        if(TMath::Abs(Jet_4vector->DeltaPhi(*Z0_4vector))<deltaphi_z_jet_min) continue;
        if(Z0_4vector->M()<muonmuon_mass_min||Z0_4vector->M()>muonmuon_mass_max) continue;

        
        // Loop over hadron 2
        for(int h_index = 0 ; h_index < tree_nomatchedjetdtrs->Jet_NDtr ; h_index++)
        {
            // Skip non-hadronic particles
            if(tree_nomatchedjetdtrs->Jet_Dtr_IsMeson[h_index]!=1&&tree_nomatchedjetdtrs->Jet_Dtr_IsBaryon[h_index]!=1) continue;
            // Skip neutrals
            if(tree_nomatchedjetdtrs->Jet_Dtr_ThreeCharge[h_index]==0) continue;
            // Apply cuts
            if(tree_nomatchedjetdtrs->Jet_Dtr_P[h_index]/1000.<track_p_min||tree_nomatchedjetdtrs->Jet_Dtr_P[h_index]/1000.>track_p_max) continue;
            if(tree_nomatchedjetdtrs->Jet_Dtr_PT[h_index]/1000.<track_pt_min) continue;
            if(tree_nomatchedjetdtrs->Jet_Dtr_TrackChi2[h_index]/tree_nomatchedjetdtrs->Jet_Dtr_TrackNDF[h_index]>track_chi2ndf_max) continue;
            if(tree_nomatchedjetdtrs->Jet_Dtr_ProbNNghost[h_index]>track_probnnghost_max) continue;

            // If all good lets count N1forcedmatched_total++;
            N1forced_total++;
            if(tree_nomatchedjetdtrs->Jet_Dtr_TRUE_ETA[h_index]!=-999)
            {
                for(int h_subindex = 0 ; h_subindex < tree_nomatchedjetdtrs->Jet_NDtr ; h_subindex++)
                {
                    if(tree_nomatchedjetdtrs->Jet_Dtr_TRUE_KEY[h_index]==tree_nomatchedjetdtrs->Jet_mcjet_dtrKeys[h_subindex]) N1forcedmatched_total++;
                }
            } 
        }
    }

    std::cout<<"\% matched particles using relatedMCP(dtr)                   = "<<N1matched_total/N1_total<<std::endl;
    std::cout<<"\% matched particles using relatedMCP(dtr) and key matching  = "<<N1forcedmatched_total/N1forced_total<<std::endl;

    for(int evt = 0 ; evt < tree_matchedjetdtrs->fChain->GetEntries() ; evt++)
    {
        // Access entry of tree
        tree_matchedjetdtrs->GetEntry(evt);

        // Apply PV cut
        if(tree_matchedjetdtrs->nPV!=1) continue;

        // Apply trigger cut
        if(tree_matchedjetdtrs->mum_L0MuonEWDecision_TOS!=1||tree_matchedjetdtrs->mum_Hlt1SingleMuonHighPTDecision_TOS!=1||tree_matchedjetdtrs->mum_Hlt2EWSingleMuonVHighPtDecision_TOS!=1) continue;
        if(tree_matchedjetdtrs->mup_L0MuonEWDecision_TOS!=1||tree_matchedjetdtrs->mup_Hlt1SingleMuonHighPTDecision_TOS!=1||tree_matchedjetdtrs->mup_Hlt2EWSingleMuonVHighPtDecision_TOS!=1) continue;

        // -999 means there is not matched jet
        if(tree_matchedjetdtrs->Jet_mcjet_nmcdtrs==-999) continue;

        // Set Jet-associated 4 vectors and apply cuts
        Jet_4vector->SetPxPyPzE(tree_matchedjetdtrs->Jet_PX/1000.,
                                tree_matchedjetdtrs->Jet_PY/1000., 
                                tree_matchedjetdtrs->Jet_PZ/1000., 
                                tree_matchedjetdtrs->Jet_PE/1000.);

        if(Jet_4vector->Eta()<jet_eta_min||Jet_4vector->Eta()>jet_eta_max) continue;
        if(Jet_4vector->Pt()<jet_pt_min_nom) continue;

        mum_4vector->SetPxPyPzE(tree_matchedjetdtrs->mum_PX/1000., 
                                tree_matchedjetdtrs->mum_PY/1000., 
                                tree_matchedjetdtrs->mum_PZ/1000., 
                                tree_matchedjetdtrs->mum_PE/1000.);

        if(Jet_4vector->DeltaR(*mum_4vector, 0/*use pseudorapidity*/)<jet_radius) continue; 
        if(mum_4vector->Pt()<muon_pt_min) continue;
        if(mum_4vector->Eta()<muon_eta_min||mum_4vector->Eta()>muon_eta_max) continue;
        if(tree_matchedjetdtrs->mum_TRACK_PCHI2<muon_trackprob_min) continue;

        mup_4vector->SetPxPyPzE(tree_matchedjetdtrs->mup_PX/1000., 
                                tree_matchedjetdtrs->mup_PY/1000., 
                                tree_matchedjetdtrs->mup_PZ/1000., 
                                tree_matchedjetdtrs->mup_PE/1000.);

        if(Jet_4vector->DeltaR(*mup_4vector, 0/*use pseudorapidity*/)<jet_radius) continue; 
        if(mup_4vector->Pt()<muon_pt_min) continue;
        if(mup_4vector->Eta()<muon_eta_min||mup_4vector->Eta()>muon_eta_max) continue;
        if(tree_matchedjetdtrs->mup_TRACK_PCHI2<muon_trackprob_min) continue;

        Z0_4vector->SetPxPyPzE(mup_4vector->Px()+mum_4vector->Px(), 
                               mup_4vector->Py()+mum_4vector->Py(), 
                               mup_4vector->Pz()+mum_4vector->Pz(), 
                               mup_4vector->E() +mum_4vector->E());

        if(TMath::Abs(Jet_4vector->DeltaPhi(*Z0_4vector))<deltaphi_z_jet_min) continue;
        if(Z0_4vector->M()<muonmuon_mass_min||Z0_4vector->M()>muonmuon_mass_max) continue;

        
        // Loop over hadron 2
        for(int h_index = 0 ; h_index < tree_matchedjetdtrs->Jet_NDtr ; h_index++)
        {
            // Skip non-hadronic particles
            if(tree_matchedjetdtrs->Jet_Dtr_IsMeson[h_index]!=1&&tree_matchedjetdtrs->Jet_Dtr_IsBaryon[h_index]!=1) continue;
            // Skip neutrals
            if(tree_matchedjetdtrs->Jet_Dtr_ThreeCharge[h_index]==0) continue;
            // Apply cuts
            if(tree_matchedjetdtrs->Jet_Dtr_P[h_index]/1000.<track_p_min||tree_matchedjetdtrs->Jet_Dtr_P[h_index]/1000.>track_p_max) continue;
            if(tree_matchedjetdtrs->Jet_Dtr_PT[h_index]/1000.<track_pt_min) continue;
            if(tree_matchedjetdtrs->Jet_Dtr_TrackChi2[h_index]/tree_matchedjetdtrs->Jet_Dtr_TrackNDF[h_index]>track_chi2ndf_max) continue;
            if(tree_matchedjetdtrs->Jet_Dtr_ProbNNghost[h_index]>track_probnnghost_max) continue;

            // If all good lets count
            N2_total++;
            if(tree_matchedjetdtrs->Jet_Dtr_TRUE_ETA[h_index]!=-999) N2matched_total++;
        }
    }

    std::cout<<"\% matched particles using relatedMCP(dtr, matchedjet_dtrs)  = "<<N2matched_total/N2_total<<std::endl;
}
