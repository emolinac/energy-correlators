#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"

void macro_print_responsematrices()
{
        TCanvas* c = new TCanvas("","",1920,1080);
        c->Draw();

        TFile* f = new TFile((output_folder + namef_ntuple_e2c_paircorrections).c_str());
        TNtuple* ntuple = (TNtuple*) f->Get(name_ntuple_correction_reco.c_str());

        float jet_pt, jet_pt_truth, R_L, R_L_truth, weight, weight_truth;
        set_unfolding_ntuple_branches(ntuple, &R_L, &R_L_truth, &jet_pt, &jet_pt_truth, &weight_pt, &weight_pt_truth);
        
        TH2F* hresp_rl     = new TH2F("hresp_rl"    ,"",nbin_rl_nominal      ,rl_nominal_binning          ,nbin_rl_nominal      ,rl_nominal_binning          );
        TH2F* hresp_jetpt  = new TH2F("hresp_jetpt" ,"",nbin_jet_pt_unfolding,unfolding_jetpt_binning,nbin_jet_pt_unfolding,unfolding_jetpt_binning);
        TH2F* hresp_weight = new TH2F("hresp_weight","",nbin_weight_unfolding,weight_unfoldingbinning,nbin_weight_unfolding,weight_unfoldingbinning);
        
        for (int evt = 0 ; evt < ntuple->GetEntries() ; evt++) {
                // Access entry of ntuple
                ntuple->GetEntry(evt);

                if (abs(R_L_truth - R_L) < rl_resolution) {
                        hresp_jetpt->Fill(jet_pt, jet_pt_truth);
                        hresp_rl->Fill(R_L, R_L_truth);
                        hresp_weight->Fill(weight, weight_truth);
                }
        }

        gStyle->SetPaintTextFormat("1.f");
        hresp_jetpt->Draw("col text");
        hresp_jetpt->SetTitle("Response matrix of p^{jet}_{t};Reco;Truth");
        gPad->SetLogx(1);
        gPad->SetLogy(1);
        c->Print("./plots/responsematrix_jetpt.pdf");
        
        hresp_weight->Draw("col text");
        hresp_weight->SetTitle("Response matrix of momentum weights;Reco;Truth");
        gPad->SetLogx(1);
        gPad->SetLogy(1);
        c->Print("./plots/responsematrix_weight.pdf");
        
        hresp_rl->GetXaxis()->SetRangeUser(0.01,1);
        hresp_rl->GetYaxis()->SetRangeUser(0.01,1); 
        hresp_rl->Draw("col text");
        hresp_rl->SetTitle("Response matrix of R_{L};Reco;Truth");
        gPad->SetLogx(1);
        gPad->SetLogy(1);
        c->Print("./plots/responsematrix_rl.pdf");
}