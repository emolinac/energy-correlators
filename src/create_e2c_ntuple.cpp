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
#include "analysis-functions.cpp"
#include "directories.h"
#include "names.h"

int main()
{
    // Create output file
    TFile* fout = new TFile((output_folder+namef_ntuple_e2c).c_str(),"RECREATE");
    //gROOT->cd();

    // Declare the TTrees to be used to build the ntuples
    TZJetsData*   datatree       = new TZJetsData();
    TZJetsMCReco* mcrecotree     = new TZJetsMCReco();
    TZJetsMCReco* mctree         = new TZJetsMCReco();
    TZJetsMC*     originalmctree = new TZJetsMC();

    gROOT->cd();

    // Create Ntuples
    TNtuple* ntuple_data   = new TNtuple(name_ntuple_data.c_str()  ,"",ntuple_data_vars); 
    TNtuple* ntuple_mcreco = new TNtuple(name_ntuple_mcreco.c_str(),"",ntuple_mcreco_vars);
    TNtuple* ntuple_mc     = new TNtuple(name_ntuple_mc.c_str()    ,"",ntuple_mc_vars);
    TNtuple* ntuple_unfold = new TNtuple(name_ntuple_unfold.c_str()    ,"",ntuple_unfold_vars);

    // Create necessary 4vectors
    TLorentzVector Jet_4vector;
    TLorentzVector Z0data_4vector;
    TLorentzVector mumdata_4vector;
    TLorentzVector mupdata_4vector;
    TVector3 hreco1_vector;
    TVector3 hreco2_vector;

    // Define array carrying the variables
    float vars[Nvars_data];

    // Fill the data TNtuple
    for(int evt = 0 ; evt < datatree->fChain->GetEntries() ; evt++)
    {
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

    std::cout<<"Data Ntuple done."<<std::endl;

    // Fill the MCReco TNtuple
    for(int evt = 0 ; evt < mcrecotree->fChain->GetEntries() ; evt++)
    {
        // Access entry of tree
        mcrecotree->GetEntry(evt);

        // Set Jet-associated 4 vectors
        Jet_4vector.SetPxPyPzE(mcrecotree->Jet_PX/1000.,
                               mcrecotree->Jet_PY/1000., 
                               mcrecotree->Jet_PZ/1000., 
                               mcrecotree->Jet_PE/1000.);

        Z0data_4vector.SetPxPyPzE(mcrecotree->Z0_PX/1000., 
                                  mcrecotree->Z0_PY/1000., 
                                  mcrecotree->Z0_PZ/1000., 
                                  mcrecotree->Z0_PE/1000.);
                
        mumdata_4vector.SetPxPyPzE(mcrecotree->mum_PX/1000., 
                                   mcrecotree->mum_PY/1000., 
                                   mcrecotree->mum_PZ/1000., 
                                   mcrecotree->mum_PE/1000.);

        mupdata_4vector.SetPxPyPzE(mcrecotree->mup_PX/1000., 
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
                vars[21] = Jet_4vector.Eta();
                vars[22] = Jet_4vector.Phi();
                vars[23] = Z0data_4vector.Phi();
                vars[24] = mumdata_4vector.Phi();
                vars[25] = mumdata_4vector.Pt();
                vars[26] = mumdata_4vector.Eta();
                vars[27] = mcrecotree->mum_PX/1000.;
                vars[28] = mcrecotree->mum_PY/1000.;
                vars[29] = mcrecotree->mum_PZ/1000.;
                vars[30] = mcrecotree->mum_PE/1000.;
                vars[31] = mumdata_4vector.M();//mcrecotree->mum_M;
                vars[32] = mcrecotree->mum_TRACK_PCHI2;
                vars[33] = mupdata_4vector.Phi();
                vars[34] = mupdata_4vector.Pt();
                vars[35] = mupdata_4vector.Eta();
                vars[36] = mcrecotree->mup_PX/1000.;
                vars[37] = mcrecotree->mup_PY/1000.;
                vars[38] = mcrecotree->mup_PZ/1000.;
                vars[39] = mcrecotree->mup_PE/1000.;
                vars[40] = mupdata_4vector.M();//mcrecotree->mup_M;
                vars[41] = mcrecotree->mup_TRACK_PCHI2;

                // Fill the TNtuple
                ntuple_mcreco->Fill(vars);
            }
        }
        
    }

    std::cout<<"MCReco Ntuple done."<<std::endl;

    float vars_mc[Nvars_mc];
    // Fill the MC TNtuple
    for(int evt = 0 ; evt < mctree->fChain->GetEntries() ; evt++)
    {
        // Access entry of tree
        mctree->GetEntry(evt);

        // Set Jet-associated 4 vectors
        Jet_4vector.SetPxPyPzE(mctree->Jet_PX/1000.,
                               mctree->Jet_PY/1000., 
                               mctree->Jet_PZ/1000., 
                               mctree->Jet_PE/1000.);

        // Loop over hadron 1
        for(int h1_index = 0 ; h1_index < mctree->Jet_mcjet_nmcdtrs ; h1_index++)
        {
            // Skip un-id'ed particles
            if(mctree->Jet_mcjet_dtrID[h1_index]==-999||mctree->Jet_mcjet_dtrID[h1_index]==0) continue;

            // Skip non-hadronic particles
            if(mctree->Jet_mcjet_dtrIsMeson[h1_index]!=1&&mctree->Jet_mcjet_dtrIsBaryon[h1_index]!=1) continue;

            // Loop over hadron 2
            for(int h2_index = 0 ; h2_index < mctree->Jet_mcjet_nmcdtrs ; h2_index++)
            {
                // Skip un-id'ed particles
                if(mctree->Jet_mcjet_dtrID[h2_index]==-999||mctree->Jet_mcjet_dtrID[h2_index]==0) continue;

                // Skip non-hadronic particles
                if(mctree->Jet_mcjet_dtrIsMeson[h2_index]!=1&&mctree->Jet_mcjet_dtrIsBaryon[h2_index]!=1) continue;

                // If all good, fille Ntuple
                vars_mc[0]  = weight(mctree->Jet_mcjet_dtrE[h1_index]/1000., mctree->Jet_mcjet_dtrE[h2_index]/1000., mctree->Jet_mcjet_PE/1000.);
                vars_mc[1]  = X_L(mctree->Jet_mcjet_dtrETA[h1_index], mctree->Jet_mcjet_dtrETA[h2_index], mctree->Jet_mcjet_dtrPHI[h1_index], mctree->Jet_mcjet_dtrPHI[h2_index]);
                vars_mc[2]  = mctree->Jet_mcjet_dtrID[h1_index];
                vars_mc[3]  = mctree->Jet_mcjet_dtrID[h2_index];
                vars_mc[4]  = mctree->Jet_mcjet_dtrETA[h1_index];
                vars_mc[5]  = mctree->Jet_mcjet_dtrETA[h2_index];
                vars_mc[6]  = mctree->Jet_mcjet_dtrPHI[h1_index];
                vars_mc[7]  = mctree->Jet_mcjet_dtrPHI[h2_index];
                vars_mc[8]  = mctree->Jet_mcjet_MotherID[h1_index];
                vars_mc[9]  = mctree->Jet_mcjet_MotherID[h2_index];
                vars_mc[10] = mctree->Jet_mcjet_TopMotherID[h1_index];
                vars_mc[11] = mctree->Jet_mcjet_TopMotherID[h2_index];
                vars_mc[12] = mctree->Jet_mcjet_dtrP[h1_index]/1000.;
                vars_mc[13] = mctree->Jet_mcjet_dtrP[h2_index]/1000.;
                vars_mc[14] = mctree->Jet_mcjet_dtrPT[h1_index]/1000.;
                vars_mc[15] = mctree->Jet_mcjet_dtrPT[h2_index]/1000.;
                vars_mc[16] = mctree->Jet_mcjet_dtrPZ[h1_index]/1000.;
                vars_mc[17] = mctree->Jet_mcjet_dtrPZ[h2_index]/1000.;
                vars_mc[18] = mctree->Jet_PT/1000.;
                vars_mc[19] = Jet_4vector.Eta();
                vars_mc[20] = Jet_4vector.Phi();

                // Fill the TNtuple
                ntuple_mc->Fill(vars_mc);
            }
        }
        
    }

    std::cout<<"MC Ntuple done."<<std::endl;

    // Fill the Unfold TNtuple
    std::cout<<"Note: To unfold I am using a weight based on total momentum"<<std::endl;
    float vars_unfold[Nvars_unfold];
    for(int evt = 0 ; evt < /*originalmctree->fChain->GetEntries()*/100000 ; evt++)
    {
        // Access entry of tree
        originalmctree->GetEntry(evt);

        // Loop over hadron 1
        for(int h1_index = 0 ; h1_index < 50 ; h1_index++)
        {
            // Skip particles out of acceptance
            if(originalmctree->MCJet_Dtr_ETA[h1_index]<eta_min||originalmctree->MCJet_Dtr_ETA[h1_index]>eta_max) continue;

            // Loop over hadron 2
            for(int h2_index = 0 ; h2_index < 50 ; h2_index++)
            {
                // Skip particles out of acceptance
                if(originalmctree->MCJet_Dtr_ETA[h2_index]<eta_min||originalmctree->MCJet_Dtr_ETA[h2_index]>eta_max) continue;

                hreco1_vector.SetXYZ(originalmctree->MCJet_Dtr_mcrecomatch_px[h1_index],originalmctree->MCJet_Dtr_mcrecomatch_py[h1_index],originalmctree->MCJet_Dtr_mcrecomatch_pz[h1_index]);
                hreco2_vector.SetXYZ(originalmctree->MCJet_Dtr_mcrecomatch_px[h2_index],originalmctree->MCJet_Dtr_mcrecomatch_py[h2_index],originalmctree->MCJet_Dtr_mcrecomatch_pz[h2_index]);

                // If all good, fille Ntuple
                vars_unfold[0]  = weight(originalmctree->MCJet_Dtr_P[h1_index]/1000., originalmctree->MCJet_Dtr_P[h2_index]/1000., originalmctree->MCJet_P/1000.);
                vars_unfold[1]  = (originalmctree->MCJet_Dtr_mcrecomatch_px[h1_index]==-999||originalmctree->MCJet_Dtr_mcrecomatch_px[h1_index]==0||
                                   originalmctree->MCJet_Dtr_mcrecomatch_px[h2_index]==-999||originalmctree->MCJet_Dtr_mcrecomatch_px[h2_index]==0)? 
                                   0 : weight(hreco1_vector.Mag()/1000., hreco1_vector.Mag()/1000., originalmctree->MCJet_P/1000.); // FIX: ENERGY OF THE RECO JET
                vars_unfold[2]  = X_L(originalmctree->MCJet_Dtr_ETA[h1_index], originalmctree->MCJet_Dtr_ETA[h2_index], originalmctree->MCJet_Dtr_PHI[h1_index], originalmctree->MCJet_Dtr_PHI[h2_index]);
                vars_unfold[3]  = (originalmctree->MCJet_Dtr_mcrecomatch_px[h1_index]==-999||originalmctree->MCJet_Dtr_mcrecomatch_px[h1_index]==0||
                                   originalmctree->MCJet_Dtr_mcrecomatch_px[h2_index]==-999||originalmctree->MCJet_Dtr_mcrecomatch_px[h2_index]==0)? 
                                   0 : X_L(hreco1_vector.Eta(), hreco2_vector.Eta(), hreco1_vector.Phi(), hreco2_vector.Phi());
                vars_unfold[4]  = originalmctree->MCJet_Dtr_ETA[h1_index];
                vars_unfold[5]  = originalmctree->MCJet_Dtr_ETA[h2_index];
                vars_unfold[6]  = originalmctree->MCJet_Dtr_PHI[h1_index];
                vars_unfold[7]  = originalmctree->MCJet_Dtr_PHI[h2_index];
                vars_unfold[8]  = originalmctree->MCJet_Dtr_P[h1_index]/1000.;
                vars_unfold[9]  = originalmctree->MCJet_Dtr_P[h2_index]/1000.;
                vars_unfold[10] = (originalmctree->MCJet_Dtr_mcrecomatch_px[h1_index]==-999||originalmctree->MCJet_Dtr_mcrecomatch_px[h1_index]==0)?0:hreco1_vector.Eta();
                vars_unfold[11] = (originalmctree->MCJet_Dtr_mcrecomatch_px[h2_index]==-999||originalmctree->MCJet_Dtr_mcrecomatch_px[h2_index]==0)?0:hreco2_vector.Eta();
                vars_unfold[12] = (originalmctree->MCJet_Dtr_mcrecomatch_px[h1_index]==-999||originalmctree->MCJet_Dtr_mcrecomatch_px[h1_index]==0)?0:hreco1_vector.Phi();
                vars_unfold[13] = (originalmctree->MCJet_Dtr_mcrecomatch_px[h2_index]==-999||originalmctree->MCJet_Dtr_mcrecomatch_px[h2_index]==0)?0:hreco2_vector.Phi();
                vars_unfold[14] = (originalmctree->MCJet_Dtr_mcrecomatch_px[h1_index]==-999||originalmctree->MCJet_Dtr_mcrecomatch_px[h1_index]==0)?0:hreco1_vector.Mag()/1000.;
                vars_unfold[15] = (originalmctree->MCJet_Dtr_mcrecomatch_px[h2_index]==-999||originalmctree->MCJet_Dtr_mcrecomatch_px[h2_index]==0)?0:hreco2_vector.Mag()/1000.;
                
                // Fill the TNtuple
                ntuple_unfold->Fill(vars_unfold);
            }
        }
        
    }

    std::cout<<"Unfold Ntuple done."<<std::endl;

    fout->cd();
    ntuple_data->Write();
    ntuple_mcreco->Write();
    ntuple_mc->Write();
    ntuple_unfold->Write();
    fout->Close();

    return 0;
}
