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
    TNtuple* ntuple = new TNtuple(name_ntuple_reco2mcdtrmatch.c_str(),"Purity Ntuple",ntuple_dtrmatch_vars); 
    
    ntuple->SetAutoSave(0);
    
    // Create necessary 4vectors
    TLorentzVector* Jet_4vector = new TLorentzVector();
    TLorentzVector* Z0_4vector  = new TLorentzVector();
    TLorentzVector* mum_4vector = new TLorentzVector();
    TLorentzVector* mup_4vector = new TLorentzVector();

    // Define array carrying the variables
    float vars[Nvars_dtrmatch];
    // Fill the MCReco TNtuple
    for(int evt = 0 ; evt < mcrecotree->fChain->GetEntries() ; evt++)
    {

        // Access entry of tree
        mcrecotree->GetEntry(evt);

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
        int signal_h1 = 0;
        int signal_h2 = 0;
        for(int h1_index = 0 ; h1_index < mcrecotree->Jet_NDtr ; h1_index++)
        {
            // Skip un-id'ed particles
            if(mcrecotree->Jet_Dtr_ID[h1_index]==-999||mcrecotree->Jet_Dtr_ID[h1_index]==0) continue;

            // Skip non-hadronic particles
            if(mcrecotree->Jet_Dtr_IsMeson[h1_index]!=1&&mcrecotree->Jet_Dtr_IsBaryon[h1_index]!=1) continue;

            // Determine if h1 has a matched truth particle
            if(mcrecotree->Jet_Dtr_TRUE_ID[h1_index]!=-999) 
            {
                // Determine if the matched truth particle is within the matched jet
                for(int entry = 0 ; entry < mcrecotree->Jet_MatchedNDtr ; entry++)
                {
                    // -999 happens when there is no matched MC jet
                    if(mcrecotree->Jet_mcjet_dtrID[entry]==-999) {signal_h1 = 0; break;}

                    // It is enough to just match a quantity like the momentum
                    if(mcrecotree->Jet_mcjet_dtrPX[h1_index]==mcrecotree->Jet_mcjet_dtrPX[entry]) {signal_h1 = 1; break;}
                }
            }
            // Loop over hadron 2
            for(int h2_index = 0 ; h2_index < mcrecotree->Jet_NDtr ; h2_index++)
            {
                // Skip un-id'ed particles
                if(mcrecotree->Jet_Dtr_ID[h2_index]==-999||mcrecotree->Jet_Dtr_ID[h2_index]==0) continue;

                // Skip non-hadronic particles
                if(mcrecotree->Jet_Dtr_IsMeson[h2_index]!=1&&mcrecotree->Jet_Dtr_IsBaryon[h2_index]!=1) continue;

                // Determine if h2 has a matched truth particle
                if(mcrecotree->Jet_Dtr_TRUE_ID[h2_index]!=-999) 
                {
                    // Determine if the matched truth particle is within the matched jet
                    for(int entry = 0 ; entry < mcrecotree->Jet_MatchedNDtr ; entry++)
                    {
                        // -999 happens when there is no matched MC jet
                        if(mcrecotree->Jet_mcjet_dtrID[entry]==-999) {signal_h2 = 0; break;}

                        // It is enough to just match a quantity like the momentum
                        if(mcrecotree->Jet_mcjet_dtrPX[h2_index]==mcrecotree->Jet_mcjet_dtrPX[entry]) {signal_h2 = 1; break;}
                    }
                }

                // If all good, fill Ntuple
                vars[0]  = weight(mcrecotree->Jet_Dtr_E[h1_index]/1000., mcrecotree->Jet_Dtr_E[h2_index]/1000., mcrecotree->Jet_PE/1000.);
                vars[1]  = R_L(mcrecotree->Jet_Dtr_ETA[h1_index], mcrecotree->Jet_Dtr_ETA[h2_index], mcrecotree->Jet_Dtr_PHI[h1_index], mcrecotree->Jet_Dtr_PHI[h2_index]);
                vars[2]  = mcrecotree->Jet_Dtr_ID[h1_index];
                vars[3]  = mcrecotree->Jet_Dtr_ID[h2_index];
                vars[4]  = mcrecotree->Jet_Dtr_ETA[h1_index];
                vars[5]  = mcrecotree->Jet_Dtr_ETA[h2_index];
                vars[6]  = mcrecotree->Jet_Dtr_PHI[h1_index];
                vars[7]  = mcrecotree->Jet_Dtr_PHI[h2_index];
                vars[8]  = mcrecotree->Jet_Dtr_TrackChi2[h1_index];
                vars[9]  = mcrecotree->Jet_Dtr_TrackChi2[h2_index];
                vars[10] = mcrecotree->Jet_Dtr_TrackNDF[h1_index];
                vars[11] = mcrecotree->Jet_Dtr_TrackNDF[h2_index]; 
                vars[12] = mcrecotree->Jet_Dtr_ProbNNghost[h1_index];
                vars[13] = mcrecotree->Jet_Dtr_ProbNNghost[h2_index];
                vars[14] = mcrecotree->Jet_Dtr_P[h1_index]/1000.;
                vars[15] = mcrecotree->Jet_Dtr_P[h2_index]/1000.;
                vars[16] = mcrecotree->Jet_Dtr_PT[h1_index]/1000.;
                vars[17] = mcrecotree->Jet_Dtr_PT[h2_index]/1000.;
                vars[18] = mcrecotree->Jet_Dtr_PZ[h1_index]/1000.;
                vars[19] = mcrecotree->Jet_Dtr_PZ[h2_index]/1000.;
                vars[20] = mcrecotree->Jet_PT/1000.;
                vars[21] = Jet_4vector->Eta();
                vars[22] = Jet_4vector->Phi();
                vars[23] = Z0_4vector->Phi();
                vars[24] = mum_4vector->Phi();
                vars[25] = mum_4vector->Pt();
                vars[26] = mum_4vector->Eta();
                vars[27] = mcrecotree->mum_PX/1000.;
                vars[28] = mcrecotree->mum_PY/1000.;
                vars[29] = mcrecotree->mum_PZ/1000.;
                vars[30] = mcrecotree->mum_PE/1000.;
                vars[31] = mum_4vector->M();//mcrecotree->mum_M;
                vars[32] = mcrecotree->mum_TRACK_PCHI2;
                vars[33] = mup_4vector->Phi();
                vars[34] = mup_4vector->Pt();
                vars[35] = mup_4vector->Eta();
                vars[36] = mcrecotree->mup_PX/1000.;
                vars[37] = mcrecotree->mup_PY/1000.;
                vars[38] = mcrecotree->mup_PZ/1000.;
                vars[39] = mcrecotree->mup_PE/1000.;
                vars[40] = mup_4vector->M();//mcrecotree->mup_M;
                vars[41] = mcrecotree->mup_TRACK_PCHI2;
                vars[42] = signal_h1;
                vars[43] = signal_h2;

                // Fill the TNtuple
                ntuple->Fill(vars);
            }

            // Reset signal variable
            signal_h1 = 0;
            signal_h2 = 0;
        }        
    }

    fout->cd();
    ntuple->Write();
    fout->Close();

    return 0;
}
