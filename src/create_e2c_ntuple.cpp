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
#include "analysis-constants.h"
#include "analysis-functions.h"
#include "analysis-functions.cpp"
#include "directories.h"
#include "names.h"

int main()
{
    // Create output file
    TFile* fout = new TFile((output_folder+namef_ntuple_e2c).c_str(),"RECREATE");
    gROOT->cd();

    // Declare the TTrees to be used to build the ntuples
    TZJetsData* datatree = new TZJetsData();
    
    // Create Ntuples
    TNtuple* ntuple_data = new TNtuple(name_ntuple_data.c_str(),"",ntuple_data_vars); 

    // Create necessary 4vectors
    TLorentzVector Jet_4vector;
    TLorentzVector Z0data_4vector;
    TLorentzVector mumdata_4vector;
    TLorentzVector mupdata_4vector;

    // Define array carrying the variables
    float vars[Nvars_data];

    // Fille the data TNtuple
    for(int evt = 0 ; evt < datatree->fChain->GetEntries() ; evt++)
    {
        std::cout<<evt*100./datatree->fChain->GetEntries()<<"\% done."<<std::endl;
        
        // Access entry of tree
        datatree->GetEntry(evt);

        // Set Jet-associated 4 vectors
        Jet_4vector.SetPxPyPzE(datatree->Jet_PX/1000.,
                               datatree->Jet_PY/1000., 
                               datatree->Jet_PZ/1000., 
                               datatree->Jet_PE/1000.);

        Z0data_4vector.SetPxPyPzE(datatree->Z0_PX/1000., 
                                  datatree->Z0_PY/1000., 
                                  datatree->Z0_PZ/1000., 
                                  datatree->Z0_PE/1000.);
                
        mumdata_4vector.SetPxPyPzE(datatree->mum_PX/1000., 
                                   datatree->mum_PY/1000., 
                                   datatree->mum_PZ/1000., 
                                   datatree->mum_PE/1000.);

        mupdata_4vector.SetPxPyPzE(datatree->mup_PX/1000., 
                                   datatree->mup_PY/1000., 
                                   datatree->mup_PZ/1000., 
                                   datatree->mup_PE/1000.);

        // Loop over hadron 1
        for(int h1_index = 0 ; h1_index < datatree->Jet_NDtr ; h1_index++)
        {
            // Skip un-id'ed particles
            if(datatree->Jet_Dtr_ID[h1_index]==-999||datatree->Jet_Dtr_ID[h1_index]==0) continue;

            // Skip non-hadronic particles
            if(datatree->Jet_Dtr_IsMeson[h1_index]!=1&&datatree->Jet_Dtr_IsBaryon[h1_index]!=1) continue;

            // Loop over hadron 2
            for(int h2_index = 0 ; h2_index < datatree->Jet_NDtr ; h2_index++)
            {
                // Skip un-id'ed particles
                if(datatree->Jet_Dtr_ID[h2_index]==-999||datatree->Jet_Dtr_ID[h2_index]==0) continue;

                // Skip non-hadronic particles
                if(datatree->Jet_Dtr_IsMeson[h2_index]!=1&&datatree->Jet_Dtr_IsBaryon[h2_index]!=1) continue;

                // If all good, fille Ntuple
                vars[0]  = weight(datatree->Jet_Dtr_E[h1_index]/1000., datatree->Jet_Dtr_E[h2_index]/1000., datatree->Jet_PE/1000.);
                vars[1]  = X_L(datatree->Jet_Dtr_ETA[h1_index], datatree->Jet_Dtr_ETA[h2_index], datatree->Jet_Dtr_PHI[h1_index], datatree->Jet_Dtr_PHI[h2_index]);
                vars[2]  = datatree->Jet_Dtr_ID[h1_index];
                vars[3]  = datatree->Jet_Dtr_ID[h2_index];
                vars[4]  = datatree->Jet_Dtr_ETA[h1_index];
                vars[5]  = datatree->Jet_Dtr_ETA[h2_index];
                vars[6]  = datatree->Jet_Dtr_PHI[h1_index];
                vars[7]  = datatree->Jet_Dtr_PHI[h2_index];
                vars[8]  = datatree->Jet_Dtr_TrackChi2[h1_index];
                vars[9]  = datatree->Jet_Dtr_TrackChi2[h2_index];
                vars[10] = datatree->Jet_Dtr_TrackNDF[h1_index];
                vars[11] = datatree->Jet_Dtr_TrackNDF[h2_index];
                vars[12] = datatree->Jet_Dtr_ProbNNghost[h1_index];
                vars[13] = datatree->Jet_Dtr_ProbNNghost[h2_index];
                vars[14] = datatree->Jet_Dtr_P[h1_index]/1000.;
                vars[15] = datatree->Jet_Dtr_P[h2_index]/1000.;
                vars[16] = datatree->Jet_Dtr_PT[h1_index]/1000.;
                vars[17] = datatree->Jet_Dtr_PT[h2_index]/1000.;
                vars[18] = datatree->Jet_Dtr_PZ[h1_index]/1000.;
                vars[19] = datatree->Jet_Dtr_PZ[h2_index]/1000.;
                vars[20] = datatree->Jet_PT/1000.;
                vars[21] = Jet_4vector.Eta();
                vars[22] = Jet_4vector.Phi();
                vars[23] = Z0data_4vector.Phi();
                vars[24] = mumdata_4vector.Phi();
                vars[25] = mumdata_4vector.Pt();
                vars[26] = mumdata_4vector.Eta();
                vars[27] = datatree->mum_PX/1000.;
                vars[28] = datatree->mum_PY/1000.;
                vars[29] = datatree->mum_PZ/1000.;
                vars[30] = datatree->mum_PE/1000.;
                vars[31] = mumdata_4vector.M();//datatree->mum_M;
                vars[32] = datatree->mum_TRACK_PCHI2;
                vars[33] = mupdata_4vector.Phi();
                vars[34] = mupdata_4vector.Pt();
                vars[35] = mupdata_4vector.Eta();
                vars[36] = datatree->mup_PX/1000.;
                vars[37] = datatree->mup_PY/1000.;
                vars[38] = datatree->mup_PZ/1000.;
                vars[39] = datatree->mup_PE/1000.;
                vars[40] = mupdata_4vector.M();//datatree->mup_M;
                vars[41] = datatree->mup_TRACK_PCHI2;

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
