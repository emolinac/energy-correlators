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
    TFile* fout = new TFile((output_folder+namef_ntuple_e2c_purity).c_str(),"RECREATE");
    
    // Declare the TTrees to be used to build the ntuples
    TZJetsMCReco* mcrecotree = new TZJetsMCReco();

    // Create Ntuples
    TNtuple* ntuple_jet_match = new TNtuple(name_ntuple_purity.c_str(),"Reco jet&dtr matched 2 MC",ntuple_purity_vars); 
    
    ntuple_jet_match->SetAutoSave(0);
    
    // Create necessary 4vectors
    TLorentzVector* Jet_4vector   = new TLorentzVector();
    TLorentzVector* Z0_4vector    = new TLorentzVector();
    TLorentzVector* mum_4vector   = new TLorentzVector();
    TLorentzVector* mup_4vector   = new TLorentzVector();
    TLorentzVector* htrue_4vector = new TLorentzVector();

    // Define array carrying the variables
    float vars[Nvars_purity];

    // Fill the matched jets Ntuple
    for(int evt = 0 ; evt < mcrecotree->fChain->GetEntries() ; evt++)
    {
        // Access entry of tree
        mcrecotree->GetEntry(evt);

        // -999 means there is not matched jet
        if(mcrecotree->Jet_mcjet_nmcdtrs==-999) continue;

        // Apply PV cut
        if(mcrecotree->nPV!=1) continue;

        // Apply trigger cut
        bool mum_trigger = (mcrecotree->mum_L0MuonEWDecision_TOS==1&&mcrecotree->mum_Hlt1SingleMuonHighPTDecision_TOS==1&&mcrecotree->mum_Hlt2EWSingleMuonVHighPtDecision_TOS==1);
        bool mup_trigger = (mcrecotree->mup_L0MuonEWDecision_TOS==1&&mcrecotree->mup_Hlt1SingleMuonHighPTDecision_TOS==1&&mcrecotree->mup_Hlt2EWSingleMuonVHighPtDecision_TOS==1);

        if(!mum_trigger&&!mup_trigger) continue;
        
        // Set Jet-associated 4 vectors and apply cuts
        Jet_4vector->SetPxPyPzE(mcrecotree->Jet_PX/1000.,
                                mcrecotree->Jet_PY/1000., 
                                mcrecotree->Jet_PZ/1000., 
                                mcrecotree->Jet_PE/1000.);

        if(Jet_4vector->Eta()<jet_eta_min||Jet_4vector->Eta()>jet_eta_max) continue;
        if(Jet_4vector->Pt()<jet_pt_min_nom) continue;

        mum_4vector->SetPxPyPzE(mcrecotree->mum_PX/1000., 
                                mcrecotree->mum_PY/1000., 
                                mcrecotree->mum_PZ/1000., 
                                mcrecotree->mum_PE/1000.);

        if(Jet_4vector->DeltaR(*mum_4vector, 0/*use pseudorapidity*/)<deltar_muon_jet_min) continue; 
        if(mum_4vector->Pt()<muon_pt_min) continue;
        if(mum_4vector->Eta()<muon_eta_min||mum_4vector->Eta()>muon_eta_max) continue;
        if(mcrecotree->mum_TRACK_PCHI2<muon_trackprob_min) continue;

        mup_4vector->SetPxPyPzE(mcrecotree->mup_PX/1000., 
                                mcrecotree->mup_PY/1000., 
                                mcrecotree->mup_PZ/1000., 
                                mcrecotree->mup_PE/1000.);

        if(Jet_4vector->DeltaR(*mup_4vector, 0/*use pseudorapidity*/)<deltar_muon_jet_min) continue; 
        if(mup_4vector->Pt()<muon_pt_min) continue;
        if(mup_4vector->Eta()<muon_eta_min||mup_4vector->Eta()>muon_eta_max) continue;
        if(mcrecotree->mup_TRACK_PCHI2<muon_trackprob_min) continue;

        Z0_4vector->SetPxPyPzE(mup_4vector->Px()+mum_4vector->Px(), 
                               mup_4vector->Py()+mum_4vector->Py(), 
                               mup_4vector->Pz()+mum_4vector->Pz(), 
                               mup_4vector->E() +mum_4vector->E());

        if(TMath::Abs(Jet_4vector->DeltaPhi(*Z0_4vector))<deltaphi_z_jet_min) continue;
        if(Z0_4vector->M()<muonmuon_mass_min||Z0_4vector->M()>muonmuon_mass_max) continue;

        
        // Loop over hadron 2
        for(int h_index = 0 ; h_index < mcrecotree->Jet_NDtr ; h_index++)
        {
            // Skip non-hadronic particles
            if(mcrecotree->Jet_Dtr_IsMeson[h_index]!=1&&mcrecotree->Jet_Dtr_IsBaryon[h_index]!=1) continue;
            
            // Skip neutrals
            if(mcrecotree->Jet_Dtr_ThreeCharge[h_index]==0) continue;

            // Apply cuts
            if(mcrecotree->Jet_Dtr_P[h_index]/1000.<track_p_min||mcrecotree->Jet_Dtr_P[h_index]/1000.>track_p_max) continue;
            if(mcrecotree->Jet_Dtr_PT[h_index]/1000.<track_pt_min) continue;
            if(mcrecotree->Jet_Dtr_TrackChi2[h_index]/mcrecotree->Jet_Dtr_TrackNDF[h_index]>track_chi2ndf_max) continue;
            if(mcrecotree->Jet_Dtr_ProbNNghost[h_index]>track_probnnghost_max) continue;

            int key_match = 0;
            if(mcrecotree->Jet_Dtr_TRUE_ETA[h_index]!=-999)
            {
                for(int h_subindex = 0 ; h_subindex < mcrecotree->Jet_NDtr ; h_subindex++)
                {
                    if(mcrecotree->Jet_Dtr_TRUE_KEY[h_index]==mcrecotree->Jet_mcjet_dtrKeys[h_subindex]) key_match++;
                }
            } 

            htrue_4vector->SetPxPyPzE(mcrecotree->Jet_Dtr_TRUE_PX[h_index]/1000.,
                                      mcrecotree->Jet_Dtr_TRUE_PY[h_index]/1000., 
                                      mcrecotree->Jet_Dtr_TRUE_PZ[h_index]/1000., 
                                      mcrecotree->Jet_Dtr_TRUE_E[h_index]/1000.);

            // If all good, fill Ntuple
            vars[0]  = mcrecotree->Jet_Dtr_ETA[h_index];
            vars[1]  = rapidity(mcrecotree->Jet_Dtr_E[h_index],mcrecotree->Jet_Dtr_PZ[h_index]);
            vars[2]  = mcrecotree->Jet_Dtr_PHI[h_index];
            vars[3]  = mcrecotree->Jet_Dtr_P[h_index]/1000.;
            vars[4]  = mcrecotree->Jet_Dtr_PT[h_index]/1000.;
            vars[5]  = mcrecotree->Jet_PT/1000.;
            vars[6]  = Jet_4vector->Eta();
            vars[7]  = Jet_4vector->DeltaPhi(*Z0_4vector);//Jet_4vector->Phi();
            vars[8]  = delta_phi(Jet_4vector->Phi(),Z0_4vector->Phi());
            vars[9]  = Jet_4vector->DeltaR(*mum_4vector, 1);
            vars[10] = mum_4vector->Pt();
            vars[11] = mum_4vector->Eta();
            vars[12] = Jet_4vector->DeltaR(*mup_4vector, 1);
            vars[13] = mup_4vector->Pt();
            vars[14] = mup_4vector->Eta();
            vars[15] = mcrecotree->Jet_PE/1000.;
            vars[16] = mcrecotree->Jet_mcjet_PE/1000.;
            vars[17] = mcrecotree->Jet_mcjet_nmcdtrs;
            vars[18] = (mcrecotree->Jet_Dtr_TRUE_ETA[h_index]==-999) ? -999 : rapidity(mcrecotree->Jet_Dtr_TRUE_E[h_index],mcrecotree->Jet_Dtr_TRUE_PZ[h_index]);
            vars[19] = (mcrecotree->Jet_Dtr_TRUE_ETA[h_index]==-999) ? -999 : mcrecotree->Jet_Dtr_TRUE_ETA[h_index];
            vars[20] = (mcrecotree->Jet_Dtr_TRUE_ETA[h_index]==-999) ? -999 : mcrecotree->Jet_Dtr_TRUE_PHI[h_index];
            vars[21] = Jet_4vector->DeltaR(*htrue_4vector);            
            vars[22] = key_match;            
            // Fill the TNtuple
            ntuple_jet_match->Fill(vars);
        }
 }

    fout->cd();
    ntuple_jet_match->Write();
    fout->Close();

    return 0;
}
