#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_create_weight_correction()
{
        // Open the necessary files
        TFile* fout = new TFile((output_folder + namef_histos_weight_corr).c_str(),"RECREATE");
        gROOT->cd();

        TFile* fcorr = new TFile((output_folder + namef_ntuple_eec_paircorrections).c_str()); 
        if (fcorr->IsZombie()) 
                return;

        TNtuple* ntuple_matches = (TNtuple*) fcorr->Get((name_ntuple_correction_reco.c_str()));

        float h1_pt, h2_pt, jet_pt, h1_pt_truth, h2_pt_truth, jet_pt_truth;
        ntuple_matches->SetBranchAddress("h1_pt",&h1_pt);
        ntuple_matches->SetBranchAddress("h2_pt",&h2_pt);
        ntuple_matches->SetBranchAddress("jet_pt",&jet_pt);
        ntuple_matches->SetBranchAddress("h1_pt_truth",&h1_pt_truth);
        ntuple_matches->SetBranchAddress("h2_pt_truth",&h2_pt_truth);
        ntuple_matches->SetBranchAddress("jet_pt_truth",&jet_pt_truth);

        TH1F* hweight_ratio[nbin_jet_pt][nbin_h_pt][nbin_h_pt];

        for (int jet_pt_bin = 0 ; jet_pt_bin < nbin_jet_pt ; jet_pt_bin++) {
                for (int h1_pt_bin = 0 ; h1_pt_bin < nbin_h_pt ; h1_pt_bin++) {
                        for (int h2_pt_bin = 0 ; h2_pt_bin < nbin_h_pt ; h2_pt_bin++) {
                                std::cout<<"Working on bin "<<jet_pt_bin<<" , "<<h1_pt_bin<<" , "<<h2_pt_bin<<std::endl;

                                hweight_ratio[jet_pt_bin][h1_pt_bin][h2_pt_bin] = new TH1F(Form("hweight_ratio%i%i%i",jet_pt_bin, h1_pt_bin, h2_pt_bin),"",100,0,2);

                                for (int entry = 0 ; entry < ntuple_matches->GetEntries() ; entry++) {
                                        ntuple_matches->GetEntry(entry);

                                        if (h1_pt < h_pt_binning[h1_pt_bin] || h1_pt > h_pt_binning[h1_pt_bin + 1])
                                                continue;

                                        if (h2_pt < h_pt_binning[h2_pt_bin] || h2_pt > h_pt_binning[h2_pt_bin + 1])
                                                continue;

                                        if (jet_pt < jet_pt_binning[jet_pt_bin] || jet_pt > jet_pt_binning[jet_pt_bin + 1])
                                                continue;

                                        double weight_reco  = h1_pt*h2_pt/jet_pt;
                                        double weight_truth = h1_pt_truth*h2_pt_truth/jet_pt_truth;

                                        hweight_ratio[jet_pt_bin][h1_pt_bin][h2_pt_bin]->Fill(weight_truth/weight_reco);
                                }

                                fout->cd();
                                hweight_ratio[jet_pt_bin][h1_pt_bin][h2_pt_bin]->Write();
                                gROOT->cd();
                        }
                }        
        }
}
