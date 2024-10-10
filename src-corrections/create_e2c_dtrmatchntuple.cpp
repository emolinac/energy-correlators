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
    TFile* fout = new TFile((output_folder+namef_ntuple_e2c_dtrmatch).c_str(),"RECREATE");
    
    // Declare the TTrees to be used to build the ntuples
    TZJetsMCReco* mcrecotree = new TZJetsMCReco();

    // Create Ntuples
    TNtuple* ntuple_jet_match = new TNtuple(name_ntuple_reco2mcdtrmatch.c_str(),"Reco jet&dtr matched 2 MC",ntuple_dtrmatch_vars); 
    
    ntuple_jet_match->SetAutoSave(0);
    
    // Create necessary 4vectors
    TLorentzVector* Jet_4vector = new TLorentzVector();
    TLorentzVector* Z0_4vector  = new TLorentzVector();
    TLorentzVector* mum_4vector = new TLorentzVector();
    TLorentzVector* mup_4vector = new TLorentzVector();

    // Create necessary 3vectors
    TVector3* delta_momentum_a = new TVector3();
    TVector3* delta_momentum_b = new TVector3();

    // Define array carrying the variables
    float vars[Nvars_dtrmatch];

    // Fill the matched jets Ntuple
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

        // NOTATION : An entry of value -999 is an invalid entry
        //            An entry of value -1000 is an entry that has to be recalculated

        // Get the locations of the matching particles
        double mcmatch_locations[(const int)mcrecotree->Jet_NDtr];
        for(int h_index = 0 ; h_index < mcrecotree->Jet_NDtr ; h_index++)
        {
            // Skip un-id'ed particles
            if(mcrecotree->Jet_Dtr_ID[h_index]==-999||mcrecotree->Jet_Dtr_ID[h_index]==0) continue;

            // Skip non-hadronic particles
            if(mcrecotree->Jet_Dtr_IsMeson[h_index]!=1&&mcrecotree->Jet_Dtr_IsBaryon[h_index]!=1) continue;

            
        }

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
                vars[0]  = weight(mcrecotree->Jet_Dtr_E[h1_index]/1000., mcrecotree->Jet_Dtr_E[h2_index]/1000., mcrecotree->Jet_PE/1000.);
                vars[1]  = R_L(h1_y, h2_y, mcrecotree->Jet_Dtr_PHI[h1_index], mcrecotree->Jet_Dtr_PHI[h2_index]);
                vars[2]  = mcrecotree->Jet_Dtr_ID[h1_index];
                vars[3]  = mcrecotree->Jet_Dtr_ID[h2_index];
                vars[4]  = mcrecotree->Jet_Dtr_ETA[h1_index];
                vars[5]  = mcrecotree->Jet_Dtr_ETA[h2_index];
                vars[6]  = h1_y;
                vars[7]  = h2_y;
                vars[8]  = mcrecotree->Jet_Dtr_PHI[h1_index];
                vars[9]  = mcrecotree->Jet_Dtr_PHI[h2_index];
                vars[10] = mcrecotree->Jet_Dtr_TrackChi2[h1_index];
                vars[11] = mcrecotree->Jet_Dtr_TrackChi2[h2_index];
                vars[12] = mcrecotree->Jet_Dtr_TrackNDF[h1_index];
                vars[13] = mcrecotree->Jet_Dtr_TrackNDF[h2_index]; 
                vars[14] = mcrecotree->Jet_Dtr_ProbNNghost[h1_index];
                vars[15] = mcrecotree->Jet_Dtr_ProbNNghost[h2_index];
                vars[16] = mcrecotree->Jet_Dtr_P[h1_index]/1000.;
                vars[17] = mcrecotree->Jet_Dtr_P[h2_index]/1000.;
                vars[18] = mcrecotree->Jet_Dtr_PT[h1_index]/1000.;
                vars[19] = mcrecotree->Jet_Dtr_PT[h2_index]/1000.;
                vars[20] = mcrecotree->Jet_Dtr_PZ[h1_index]/1000.;
                vars[21] = mcrecotree->Jet_Dtr_PZ[h2_index]/1000.;
                vars[22] = mcrecotree->Jet_PT/1000.;
                vars[23] = Jet_4vector->Eta();
                vars[24] = Jet_4vector->Phi();
                vars[25] = Z0_4vector->Phi();
                vars[26] = mum_4vector->Phi();
                vars[27] = mum_4vector->Pt();
                vars[28] = mum_4vector->Eta();
                vars[29] = mcrecotree->mum_PX/1000.;
                vars[30] = mcrecotree->mum_PY/1000.;
                vars[31] = mcrecotree->mum_PZ/1000.;
                vars[32] = mcrecotree->mum_PE/1000.;
                vars[33] = mum_4vector->M();//mcrecotree->mum_M;
                vars[34] = mcrecotree->mum_TRACK_PCHI2;
                vars[35] = mup_4vector->Phi();
                vars[36] = mup_4vector->Pt();
                vars[37] = mup_4vector->Eta();
                vars[38] = mcrecotree->mup_PX/1000.;
                vars[39] = mcrecotree->mup_PY/1000.;
                vars[40] = mcrecotree->mup_PZ/1000.;
                vars[41] = mcrecotree->mup_PE/1000.;
                vars[42] = mup_4vector->M();//mcrecotree->mup_M;
                vars[43] = mcrecotree->mup_TRACK_PCHI2;
                vars[44] = mcrecotree->Jet_PE/1000.;
                vars[45] = mcrecotree->Jet_mcjet_PE/1000.;
                vars[46] = mcrecotree->Jet_mcjet_nmcdtrs;

                double matchedmc_y1   = rapidity(mcrecotree->Jet_Dtr_TRUE_E[h1_index],mcrecotree->Jet_Dtr_TRUE_PZ[h1_index]);
                double matchedmc_y2   = rapidity(mcrecotree->Jet_Dtr_TRUE_E[h2_index],mcrecotree->Jet_Dtr_TRUE_PZ[h2_index]);
                double matchedmc_phi1 = mcrecotree->Jet_Dtr_TRUE_PHI[h1_index];
                double matchedmc_phi2 = mcrecotree->Jet_Dtr_TRUE_PHI[h2_index];

                vars[47] = (mcrecotree->Jet_Dtr_TRUE_ETA[h1_index]==-999) ? -999 : R_L(h1_y, matchedmc_y1, mcrecotree->Jet_Dtr_PHI[h1_index], matchedmc_phi1);
                vars[48] = (mcrecotree->Jet_Dtr_TRUE_ETA[h2_index]==-999) ? -999 : R_L(h2_y, matchedmc_y2, mcrecotree->Jet_Dtr_PHI[h2_index], matchedmc_phi2);
                vars[49] = (mcrecotree->Jet_Dtr_TRUE_ETA[h1_index]==-999) ? -999 : matchedmc_y1;
                vars[50] = (mcrecotree->Jet_Dtr_TRUE_ETA[h2_index]==-999) ? -999 : matchedmc_y2;
                vars[51] = (mcrecotree->Jet_Dtr_TRUE_ETA[h1_index]==-999) ? -999 : matchedmc_phi1;
                vars[52] = (mcrecotree->Jet_Dtr_TRUE_ETA[h2_index]==-999) ? -999 : matchedmc_phi2;
                vars[53] = (mcrecotree->Jet_Dtr_TRUE_ETA[h2_index]==-999||mcrecotree->Jet_Dtr_TRUE_ETA[h1_index]==-999) ? 
                            -999 : R_L(matchedmc_y1,matchedmc_y2,matchedmc_phi1,matchedmc_phi2);

                // Fill the TNtuple
                ntuple_jet_match->Fill(vars);
            }
        }        
    }

    fout->cd();
    ntuple_jet_match->Write();
    fout->Close();

    return 0;
}
