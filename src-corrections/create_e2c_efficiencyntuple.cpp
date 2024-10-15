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

    // Create necessary 3vectors
    TVector3* delta_momentum_a = new TVector3();
    TVector3* delta_momentum_b = new TVector3();

    // Define array carrying the variables
    float vars_reco[Nvars_efficiency_reco];
    float vars_mc[Nvars_efficiency_mc];

    // Fill the TNtuple associated to reco values
    for(int evt = 0 ; evt < mcrecotree->fChain->GetEntries() ; evt++)
    {
        // Access entry of tree
        mcrecotree->GetEntry(evt);

        // -999 means there is not matched jet
        if(mcrecotree->Jet_mcjet_nmcdtrs==-999) continue;

        // Set Jet-associated 4 vectors
        Jet_4vector->SetPxPyPzE(mcrecotree->Jet_PX/1000.,
                                mcrecotree->Jet_PY/1000., 
                                mcrecotree->Jet_PZ/1000., 
                                mcrecotree->Jet_PE/1000.);

        Z0_4vector->SetPxPyPzE(mcrecotree->Z0_PX/1000., 
                               mcrecotree->Z0_PY/1000., 
                               mcrecotree->Z0_PZ/1000., 
                               mcrecotree->Z0_PE/1000.);
                
        mum_4vector->SetPxPyPzE(mcrecotree->mum_PX/1000., 
                                mcrecotree->mum_PY/1000., 
                                mcrecotree->mum_PZ/1000., 
                                mcrecotree->mum_PE/1000.);

        mup_4vector->SetPxPyPzE(mcrecotree->mup_PX/1000., 
                                mcrecotree->mup_PY/1000., 
                                mcrecotree->mup_PZ/1000., 
                                mcrecotree->mup_PE/1000.);

        // Loop over hadron 1
        for(int h1_index = 0 ; h1_index < mcrecotree->Jet_NDtr ; h1_index++)
        {
            // Skip un-id'ed particles
            if(mcrecotree->Jet_Dtr_ID[h1_index]==-999||mcrecotree->Jet_Dtr_ID[h1_index]==0) continue;

            // Skip non-hadronic particles
            if(mcrecotree->Jet_Dtr_IsMeson[h1_index]!=1&&mcrecotree->Jet_Dtr_IsBaryon[h1_index]!=1) continue;

            // Loop over hadron 2
            for(int h2_index = 0 ; h2_index < mcrecotree->Jet_NDtr ; h2_index++)
            {
                // Skip un-id'ed particles
                if(mcrecotree->Jet_Dtr_ID[h2_index]==-999||mcrecotree->Jet_Dtr_ID[h2_index]==0) continue;

                // Skip non-hadronic particles
                if(mcrecotree->Jet_Dtr_IsMeson[h2_index]!=1&&mcrecotree->Jet_Dtr_IsBaryon[h2_index]!=1) continue;

                double h1_y = rapidity(mcrecotree->Jet_Dtr_E[h1_index],mcrecotree->Jet_Dtr_PZ[h1_index]);                
                double h2_y = rapidity(mcrecotree->Jet_Dtr_E[h2_index],mcrecotree->Jet_Dtr_PZ[h2_index]);

                // If all good, fill Ntuple
                vars_reco[0]  = weight(mcrecotree->Jet_Dtr_E[h1_index]/1000., mcrecotree->Jet_Dtr_E[h2_index]/1000., mcrecotree->Jet_PE/1000.);
                vars_reco[1]  = R_L(h1_y, h2_y, mcrecotree->Jet_Dtr_PHI[h1_index], mcrecotree->Jet_Dtr_PHI[h2_index]);
                vars_reco[2]  = mcrecotree->Jet_Dtr_ThreeCharge[h1_index];
                vars_reco[3]  = mcrecotree->Jet_Dtr_ThreeCharge[h2_index];
                vars_reco[4]  = mcrecotree->Jet_Dtr_ETA[h1_index];
                vars_reco[5]  = mcrecotree->Jet_Dtr_ETA[h2_index];
                vars_reco[6] = mcrecotree->Jet_Dtr_TrackChi2[h1_index];
                vars_reco[7] = mcrecotree->Jet_Dtr_TrackChi2[h2_index];
                vars_reco[8] = mcrecotree->Jet_Dtr_TrackNDF[h1_index];
                vars_reco[9] = mcrecotree->Jet_Dtr_TrackNDF[h2_index]; 
                vars_reco[10] = mcrecotree->Jet_Dtr_ProbNNghost[h1_index];
                vars_reco[11] = mcrecotree->Jet_Dtr_ProbNNghost[h2_index];
                vars_reco[12] = mcrecotree->Jet_Dtr_P[h1_index]/1000.;
                vars_reco[13] = mcrecotree->Jet_Dtr_P[h2_index]/1000.;
                vars_reco[14] = mcrecotree->Jet_Dtr_PT[h1_index]/1000.;
                vars_reco[15] = mcrecotree->Jet_Dtr_PT[h2_index]/1000.;
                vars_reco[16] = mcrecotree->Jet_PT/1000.;
                vars_reco[17] = Jet_4vector->Eta();
                vars_reco[18] = Jet_4vector->Phi();
                vars_reco[19] = delta_phi(Jet_4vector->Phi(),Z0_4vector->Phi());
                vars_reco[20] = delta_phi(Jet_4vector->Phi(),mum_4vector->Phi());
                vars_reco[21] = mum_4vector->Pt();
                vars_reco[22] = mum_4vector->Eta();
                vars_reco[23] = mcrecotree->mum_PX/1000.;
                vars_reco[24] = mcrecotree->mum_PY/1000.;
                vars_reco[25] = mcrecotree->mum_PZ/1000.;
                vars_reco[26] = mcrecotree->mum_PE/1000.;
                vars_reco[27] = mum_4vector->M();//mcrecotree->mum_M;
                vars_reco[28] = mcrecotree->mum_TRACK_PCHI2;
                vars_reco[29] = delta_phi(Jet_4vector->Phi(),mup_4vector->Phi());
                vars_reco[30] = mup_4vector->Pt();
                vars_reco[31] = mup_4vector->Eta();
                vars_reco[32] = mcrecotree->mup_PX/1000.;
                vars_reco[33] = mcrecotree->mup_PY/1000.;
                vars_reco[34] = mcrecotree->mup_PZ/1000.;
                vars_reco[35] = mcrecotree->mup_PE/1000.;
                vars_reco[36] = mup_4vector->M();//mcrecotree->mup_M;
                vars_reco[37] = mcrecotree->mup_TRACK_PCHI2;
                vars_reco[38] = mcrecotree->Jet_PE/1000.;
               
                double matchedmc_y1   = rapidity(mcrecotree->Jet_Dtr_TRUE_E[h1_index],mcrecotree->Jet_Dtr_TRUE_PZ[h1_index]);
                double matchedmc_y2   = rapidity(mcrecotree->Jet_Dtr_TRUE_E[h2_index],mcrecotree->Jet_Dtr_TRUE_PZ[h2_index]);
                double matchedmc_phi1 = mcrecotree->Jet_Dtr_TRUE_PHI[h1_index];
                double matchedmc_phi2 = mcrecotree->Jet_Dtr_TRUE_PHI[h2_index];

                vars_reco[39] = (mcrecotree->Jet_Dtr_TRUE_ETA[h1_index]==-999) ? -999 : R_L(h1_y, matchedmc_y1, mcrecotree->Jet_Dtr_PHI[h1_index], matchedmc_phi1);
                vars_reco[40] = (mcrecotree->Jet_Dtr_TRUE_ETA[h2_index]==-999) ? -999 : R_L(h2_y, matchedmc_y2, mcrecotree->Jet_Dtr_PHI[h2_index], matchedmc_phi2);
                vars_reco[41] = (mcrecotree->Jet_Dtr_TRUE_ETA[h1_index]==-999) ? -999 : matchedmc_y1;
                vars_reco[42] = (mcrecotree->Jet_Dtr_TRUE_ETA[h2_index]==-999) ? -999 : matchedmc_y2;
                vars_reco[43] = (mcrecotree->Jet_Dtr_TRUE_ETA[h1_index]==-999) ? -999 : matchedmc_phi1;
                vars_reco[44] = (mcrecotree->Jet_Dtr_TRUE_ETA[h2_index]==-999) ? -999 : matchedmc_phi2;
                vars_reco[45] = (mcrecotree->Jet_Dtr_TRUE_ETA[h2_index]==-999||mcrecotree->Jet_Dtr_TRUE_ETA[h1_index]==-999) ? 
                                -999 : R_L(matchedmc_y1,matchedmc_y2,matchedmc_phi1,matchedmc_phi2);

                // Fill the TNtuple
                ntuple_reco->Fill(vars_reco);
            }
        }        
    }

    // Fill the MC TNtuple
    for(int evt = 0 ; evt < mcrecotree->fChain->GetEntries() ; evt++)
    {
        // Access entry of tree
        mcrecotree->GetEntry(evt);

        for(int h1_index = 0 ; h1_index < mcrecotree->Jet_mcjet_nmcdtrs ; h1_index++)
        {
            // Skip non-hadronic particles
            if(mcrecotree->Jet_mcjet_dtrIsBaryon[h1_index]!=1&&mcrecotree->Jet_mcjet_dtrIsMeson[h1_index]!=1) continue;

            for(int h2_index = 0 ; h2_index < mcrecotree->Jet_mcjet_nmcdtrs ; h2_index++)
            {
                // Skip non-hadronic particles
                if(mcrecotree->Jet_mcjet_dtrIsBaryon[h2_index]!=1&&mcrecotree->Jet_mcjet_dtrIsMeson[h2_index]!=1) continue;

                double h1_y = rapidity(mcrecotree->Jet_mcjet_dtrE[h1_index],mcrecotree->Jet_mcjet_dtrPZ[h1_index]); 
                double h2_y = rapidity(mcrecotree->Jet_mcjet_dtrE[h2_index],mcrecotree->Jet_mcjet_dtrPZ[h2_index]);

                vars_mc[0]  = weight(mcrecotree->Jet_mcjet_dtrE[h1_index]/1000.,mcrecotree->Jet_mcjet_dtrE[h2_index]/1000.,mcrecotree->Jet_mcjet_PE/1000.);
                vars_mc[1]  = R_L(h1_y, h2_y, mcrecotree->Jet_mcjet_dtrPHI[h1_index], mcrecotree->Jet_mcjet_dtrPHI[h2_index]);
                vars_mc[2]  = mcrecotree->Jet_mcjet_dtrThreeCharge[h1_index];
                vars_mc[3]  = mcrecotree->Jet_mcjet_dtrThreeCharge[h2_index];
                vars_mc[4]  = mcrecotree->Jet_mcjet_dtrETA[h1_index];
                vars_mc[5]  = mcrecotree->Jet_mcjet_dtrETA[h2_index];
                vars_mc[6]  = h1_y;
                vars_mc[7]  = h2_y;
                vars_mc[8]  = mcrecotree->Jet_mcjet_dtrPHI[h1_index];
                vars_mc[9]  = mcrecotree->Jet_mcjet_dtrPHI[h2_index];
                vars_mc[10] = mcrecotree->Jet_mcjet_dtrP[h1_index]/1000.;
                vars_mc[11] = mcrecotree->Jet_mcjet_dtrP[h2_index]/1000.;
                vars_mc[12] = mcrecotree->Jet_mcjet_dtrPT[h1_index]/1000.;
                vars_mc[13] = mcrecotree->Jet_mcjet_dtrPT[h2_index]/1000.; 
                vars_mc[14] = mcrecotree->Jet_mcjet_dtrPZ[h1_index]/1000.;
                vars_mc[15] = mcrecotree->Jet_mcjet_dtrPZ[h2_index]/1000.;
                vars_mc[16] = mcrecotree->Jet_mcjet_PT/1000.;
                vars_mc[17] = mcrecotree->Jet_mcjet_ETA;
                vars_mc[18] = mcrecotree->Jet_mcjet_PHI;

                double mum_px = mcrecotree->Jet_mcjet_mum_PX/1000.;
                double mum_py = mcrecotree->Jet_mcjet_mum_PY/1000.;
                double mum_pz = mcrecotree->Jet_mcjet_mum_PZ/1000.;
                double mum_e  = mcrecotree->Jet_mcjet_mum_PE/1000.;

                double mup_px = mcrecotree->Jet_mcjet_mup_PX/1000.;
                double mup_py = mcrecotree->Jet_mcjet_mup_PY/1000.;
                double mup_pz = mcrecotree->Jet_mcjet_mup_PZ/1000.;
                double mup_e  = mcrecotree->Jet_mcjet_mup_PE/1000.;

                // Muon branches
                Z0_4vector->SetPxPyPzE(mum_px+mup_px, mum_py+mup_py, mum_pz+mup_pz, mum_e+mup_e);

                vars_mc[19] = delta_phi(mcrecotree->Jet_mcjet_PHI,Z0_4vector->Phi());
                vars_mc[20] = delta_phi(mcrecotree->Jet_mcjet_PHI,mcrecotree->Jet_mcjet_mum_PHI);
                vars_mc[21] = mcrecotree->Jet_mcjet_mum_PT/1000.;
                vars_mc[22] = mcrecotree->Jet_mcjet_mum_ETA;
                vars_mc[23] = mum_px;
                vars_mc[24] = mum_py;
                vars_mc[25] = mum_pz;
                vars_mc[26] = mum_e;
                vars_mc[27] = mcrecotree->Jet_mcjet_mum_M/1000.;
                vars_mc[28] = delta_phi(mcrecotree->Jet_mcjet_PHI,mcrecotree->Jet_mcjet_mup_PHI);
                vars_mc[29] = mcrecotree->Jet_mcjet_mup_PT/1000.;
                vars_mc[30] = mcrecotree->Jet_mcjet_mup_ETA;
                vars_mc[31] = mup_px;
                vars_mc[32] = mup_py;
                vars_mc[33] = mup_pz;
                vars_mc[34] = mup_e;
                vars_mc[35] = mcrecotree->Jet_mcjet_mup_M/1000.;

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
