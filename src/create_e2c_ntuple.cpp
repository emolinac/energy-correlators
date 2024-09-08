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
    TFile* fout = new TFile((output_folder+namef_ntuple_e2c).c_str(),"RECREATE");
    
    // Declare the TTrees to be used to build the ntuples
    TZJetsData*   datatree   = new TZJetsData();
    TZJetsMCReco* mcrecotree = new TZJetsMCReco();
    TZJetsMC*     mctree     = new TZJetsMC();

    //gROOT->cd();

    // Create Ntuples
    TNtuple* ntuple_data       = new TNtuple(name_ntuple_data.c_str()      ,"Data"                        ,ntuple_data_vars      ); 
    TNtuple* ntuple_mcreco     = new TNtuple(name_ntuple_mcreco.c_str()    ,"Reco Sim"                    ,ntuple_mcreco_vars    );
    TNtuple* ntuple_mcjetmatch = new TNtuple(name_ntuple_mcjetmatch.c_str(),"MC jet matched from reco sim",ntuple_mcjetmatch_vars);
    TNtuple* ntuple_mc         = new TNtuple(name_ntuple_mc.c_str()        ,"MC Sim"                      ,ntuple_mc_vars        );
    
    ntuple_data->SetAutoSave(0);
    ntuple_mcreco->SetAutoSave(0);
    ntuple_mcjetmatch->SetAutoSave(0);
    ntuple_mc->SetAutoSave(0);

    // Create necessary 4vectors
    TLorentzVector* Jet_4vector = new TLorentzVector();
    TLorentzVector* Z0_4vector  = new TLorentzVector();
    TLorentzVector* mum_4vector = new TLorentzVector();
    TLorentzVector* mup_4vector = new TLorentzVector();

    // Define array carrying the variables
    float vars[Nvars_data];
    // Fill the data TNtuple
    for(int evt = 0 ; evt < datatree->fChain->GetEntries() ; evt++)
    {
        // Access entry of tree
        datatree->GetEntry(evt);

        // Set Jet-associated 4 vectors
        Jet_4vector->SetPxPyPzE(datatree->Jet_PX/1000.,
                               datatree->Jet_PY/1000., 
                               datatree->Jet_PZ/1000., 
                               datatree->Jet_PE/1000.);

        Z0_4vector->SetPxPyPzE(datatree->Z0_PX/1000., 
                                  datatree->Z0_PY/1000., 
                                  datatree->Z0_PZ/1000., 
                                  datatree->Z0_PE/1000.);
                
        mum_4vector->SetPxPyPzE(datatree->mum_PX/1000., 
                                   datatree->mum_PY/1000., 
                                   datatree->mum_PZ/1000., 
                                   datatree->mum_PE/1000.);

        mup_4vector->SetPxPyPzE(datatree->mup_PX/1000., 
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
                vars[21] = Jet_4vector->Eta();
                vars[22] = Jet_4vector->Phi();
                vars[23] = Z0_4vector->Phi();
                vars[24] = mum_4vector->Phi();
                vars[25] = mum_4vector->Pt();
                vars[26] = mum_4vector->Eta();
                vars[27] = datatree->mum_PX/1000.;
                vars[28] = datatree->mum_PY/1000.;
                vars[29] = datatree->mum_PZ/1000.;
                vars[30] = datatree->mum_PE/1000.;
                vars[31] = mum_4vector->M();//datatree->mum_M;
                vars[32] = datatree->mum_TRACK_PCHI2;
                vars[33] = mup_4vector->Phi();
                vars[34] = mup_4vector->Pt();
                vars[35] = mup_4vector->Eta();
                vars[36] = datatree->mup_PX/1000.;
                vars[37] = datatree->mup_PY/1000.;
                vars[38] = datatree->mup_PZ/1000.;
                vars[39] = datatree->mup_PE/1000.;
                vars[40] = mup_4vector->M();//datatree->mup_M;
                vars[41] = datatree->mup_TRACK_PCHI2;

                // Fill the TNtuple
                ntuple_data->Fill(vars);
            }
        }
        
    }

    std::cout<<"Data Ntuple done."<<std::endl;

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

                // If all good, fille Ntuple
                vars[0]  = weight(mcrecotree->Jet_Dtr_E[h1_index]/1000., mcrecotree->Jet_Dtr_E[h2_index]/1000., mcrecotree->Jet_PE/1000.);
                vars[1]  = X_L(mcrecotree->Jet_Dtr_ETA[h1_index], mcrecotree->Jet_Dtr_ETA[h2_index], mcrecotree->Jet_Dtr_PHI[h1_index], mcrecotree->Jet_Dtr_PHI[h2_index]);
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

                // Fill the TNtuple
                ntuple_mcreco->Fill(vars);
            }
        }
        
    }

    std::cout<<"MCReco Ntuple done."<<std::endl;

/*    // Fill the MCJetMatch TNtuple
    float vars_mcjetmatch[Nvars_mcjetmatch];
    for(int evt = 0 ; evt < mcrecotree->fChain->GetEntries() ; evt++)
    {
        // Access entry of tree
        mcrecotree->GetEntry(evt);

        // Set Jet-associated 4 vectors
        Jet_4vector->SetPxPyPzE(mcrecotree->Jet_PX/1000.,
                               mcrecotree->Jet_PY/1000., 
                               mcrecotree->Jet_PZ/1000., 
                               mcrecotree->Jet_PE/1000.);

        // Loop over hadron 1
        for(int h1_index = 0 ; h1_index < mcrecotree->Jet_mcjet_nmcdtrs ; h1_index++)
        {
            // Skip un-id'ed particles
            if(mcrecotree->Jet_mcjet_dtrID[h1_index]==-999||mcrecotree->Jet_mcjet_dtrID[h1_index]==0) continue;

            // Skip non-hadronic particles
            if(mcrecotree->Jet_mcjet_dtrIsMeson[h1_index]!=1&&mcrecotree->Jet_mcjet_dtrIsBaryon[h1_index]!=1) continue;

            // Loop over hadron 2
            for(int h2_index = 0 ; h2_index < mcrecotree->Jet_mcjet_nmcdtrs ; h2_index++)
            {
                // Skip un-id'ed particles
                if(mcrecotree->Jet_mcjet_dtrID[h2_index]==-999||mcrecotree->Jet_mcjet_dtrID[h2_index]==0) continue;

                // Skip non-hadronic particles
                if(mcrecotree->Jet_mcjet_dtrIsMeson[h2_index]!=1&&mcrecotree->Jet_mcjet_dtrIsBaryon[h2_index]!=1) continue;

                // If all good, fille Ntuple
                vars_mcjetmatch[0]  = weight(mcrecotree->Jet_mcjet_dtrE[h1_index]/1000., mcrecotree->Jet_mcjet_dtrE[h2_index]/1000., mcrecotree->Jet_mcjet_PE/1000.);
                vars_mcjetmatch[1]  = X_L(mcrecotree->Jet_mcjet_dtrETA[h1_index], mcrecotree->Jet_mcjet_dtrETA[h2_index], mcrecotree->Jet_mcjet_dtrPHI[h1_index], mcrecotree->Jet_mcjet_dtrPHI[h2_index]);
                vars_mcjetmatch[2]  = mcrecotree->Jet_mcjet_dtrID[h1_index];
                vars_mcjetmatch[3]  = mcrecotree->Jet_mcjet_dtrID[h2_index];
                vars_mcjetmatch[4]  = mcrecotree->Jet_mcjet_dtrETA[h1_index];
                vars_mcjetmatch[5]  = mcrecotree->Jet_mcjet_dtrETA[h2_index];
                vars_mcjetmatch[6]  = mcrecotree->Jet_mcjet_dtrPHI[h1_index];
                vars_mcjetmatch[7]  = mcrecotree->Jet_mcjet_dtrPHI[h2_index];
                vars_mcjetmatch[8]  = mcrecotree->Jet_mcjet_MotherID[h1_index];
                vars_mcjetmatch[9]  = mcrecotree->Jet_mcjet_MotherID[h2_index];
                vars_mcjetmatch[10] = mcrecotree->Jet_mcjet_TopMotherID[h1_index];
                vars_mcjetmatch[11] = mcrecotree->Jet_mcjet_TopMotherID[h2_index];
                vars_mcjetmatch[12] = mcrecotree->Jet_mcjet_dtrP[h1_index]/1000.;
                vars_mcjetmatch[13] = mcrecotree->Jet_mcjet_dtrP[h2_index]/1000.;
                vars_mcjetmatch[14] = mcrecotree->Jet_mcjet_dtrPT[h1_index]/1000.;
                vars_mcjetmatch[15] = mcrecotree->Jet_mcjet_dtrPT[h2_index]/1000.;
                vars_mcjetmatch[16] = mcrecotree->Jet_mcjet_dtrPZ[h1_index]/1000.;
                vars_mcjetmatch[17] = mcrecotree->Jet_mcjet_dtrPZ[h2_index]/1000.;
                vars_mcjetmatch[18] = mcrecotree->Jet_PT/1000.;
                vars_mcjetmatch[19] = Jet_4vector->Eta();
                vars_mcjetmatch[20] = Jet_4vector->Phi();

                // Fill the TNtuple
                ntuple_mcjetmatch->Fill(vars_mcjetmatch);
            }
        }
        
    }
    std::cout<<"MCJetMatch Ntuple done."<<std::endl;
*/
    // Fill the MC TNtuple
    float vars_mc[Nvars_mc];
    for(int evt = 0 ; evt < mctree->fChain->GetEntries() ; evt++)
    {
        // Access entry of tree
        mctree->GetEntry(evt);

        for(int h1_index = 0 ; h1_index < mctree->MCJet_Dtr_nmcdtrs ; h1_index++)
        {
            // Skip non-hadronic particles
            if(mctree->MCJet_Dtr_IsBaryon[h1_index]!=1&&mctree->MCJet_Dtr_IsMeson[h1_index]!=1) continue;

            for(int h2_index = 0 ; h2_index < mctree->MCJet_Dtr_nmcdtrs ; h2_index++)
            {
                // Skip non-hadronic particles
                if(mctree->MCJet_Dtr_IsBaryon[h2_index]!=1&&mctree->MCJet_Dtr_IsMeson[h2_index]!=1) continue;

                vars_mc[0]  = weight(mctree->MCJet_Dtr_E[h1_index]/1000.,mctree->MCJet_Dtr_E[h2_index]/1000.,mctree->MCJet_PE/1000.);
                vars_mc[1]  = X_L(mctree->MCJet_Dtr_ETA[h1_index],mctree->MCJet_Dtr_ETA[h2_index],mctree->MCJet_Dtr_PHI[h1_index],mctree->MCJet_Dtr_PHI[h2_index]);
                vars_mc[2]  = mctree->MCJet_Dtr_ID[h1_index];
                vars_mc[3]  = mctree->MCJet_Dtr_ID[h2_index];
                vars_mc[4]  = mctree->MCJet_Dtr_ETA[h1_index];
                vars_mc[5]  = mctree->MCJet_Dtr_ETA[h2_index];
                vars_mc[6]  = mctree->MCJet_Dtr_PHI[h1_index];
                vars_mc[7]  = mctree->MCJet_Dtr_PHI[h2_index];
                vars_mc[8]  = mctree->MCJet_Dtr_P[h1_index]/1000.;
                vars_mc[9]  = mctree->MCJet_Dtr_P[h2_index]/1000.;
                vars_mc[10] = mctree->MCJet_Dtr_PT[h1_index]/1000.;
                vars_mc[11] = mctree->MCJet_Dtr_PT[h2_index]/1000.; 
                vars_mc[12] = mctree->MCJet_Dtr_PZ[h1_index]/1000.;
                vars_mc[13] = mctree->MCJet_Dtr_PZ[h2_index]/1000.;
                vars_mc[14] = mctree->MCJet_PT/1000.;
                vars_mc[15] = mctree->MCJet_ETA;
                vars_mc[16] = mctree->MCJet_PHI;

                double mum_px = mctree->MCJet_truth_mum_PX/1000.;
                double mum_py = mctree->MCJet_truth_mum_PY/1000.;
                double mum_pz = mctree->MCJet_truth_mum_PZ/1000.;
                double mum_e  = mctree->MCJet_truth_mum_PE/1000.;

                double mup_px = mctree->MCJet_truth_mup_PX/1000.;
                double mup_py = mctree->MCJet_truth_mup_PY/1000.;
                double mup_pz = mctree->MCJet_truth_mup_PZ/1000.;
                double mup_e  = mctree->MCJet_truth_mup_PE/1000.;

                // Muon branches
                Z0_4vector->SetPxPyPzE(mum_px+mup_px, mum_py+mup_py, mum_pz+mup_pz, mum_e+mup_e);

                vars_mc[17] = Z0_4vector->Phi();
                vars_mc[18] = mctree->MCJet_truth_mum_PHI;
                vars_mc[19] = mctree->MCJet_truth_mum_PT/1000.;
                vars_mc[20] = mctree->MCJet_truth_mum_ETA;
                vars_mc[21] = mum_px;
                vars_mc[22] = mum_py;
                vars_mc[23] = mum_pz;
                vars_mc[24] = mum_e;
                vars_mc[25] = mctree->MCJet_truth_mum_M/1000.;
                vars_mc[26] = mctree->MCJet_truth_mup_PHI;
                vars_mc[27] = mctree->MCJet_truth_mup_PT/1000.;
                vars_mc[28] = mctree->MCJet_truth_mup_ETA;
                vars_mc[29] = mup_px;
                vars_mc[30] = mup_py;
                vars_mc[31] = mup_pz;
                vars_mc[32] = mup_e;
                vars_mc[33] = mctree->MCJet_truth_mup_M/1000.;

                // Fill the TNtuple
                ntuple_mc->Fill(vars_mc);        
            }   
        }
    }

    std::cout<<"MC Ntuple done."<<std::endl;

    fout->cd();
    ntuple_data->Write();
    ntuple_mcreco->Write();
    ntuple_mcjetmatch->Write();
    ntuple_mc->Write();
    fout->Close();

    return 0;
}
