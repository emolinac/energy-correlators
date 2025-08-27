#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"

void macro_print_responsematrices()
{
        TCanvas* c = new TCanvas("","",1920,1080);
        c->Draw();

        TFile* f = new TFile((output_folder + namef_ntuple_eec_paircorrections).c_str());
        TNtuple* ntuple = (TNtuple*) f->Get(name_ntuple_correction_reco.c_str());

        float jet_pt, jet_pt_truth, R_L, R_L_truth, weight_pt, weight_pt_truth, h1_pt, h1_pt_truth, h2_pt, h2_pt_truth;
        set_unfolding_ntuple_branches(ntuple, &R_L, &R_L_truth, &jet_pt, &jet_pt_truth, &weight_pt, &weight_pt_truth);
        ntuple->SetBranchAddress("h1_pt",&h1_pt);
        ntuple->SetBranchAddress("h1_pt_truth",&h1_pt_truth);
        ntuple->SetBranchAddress("h2_pt",&h2_pt);
        ntuple->SetBranchAddress("h2_pt_truth",&h2_pt_truth);
        
        TH2F* hresp_rl     = new TH2F("hresp_rl"    ,"",200,0.008,0.8,200,0.008,0.8);
        TH2F* hresp_jet_pt = new TH2F("hresp_jet_pt","",200,10,150,200,10,150);
        TH2F* hresp_weight = new TH2F("hresp_weight","",nbin_weight,weight_binning,nbin_weight,weight_binning);
        TH2F* hresp_ptprod = new TH2F("hresp_ptprod","",nbin_weight,ptprod_binning,nbin_weight,ptprod_binning);

        for (int evt = 0 ; evt < ntuple->GetEntries() ; evt++) {
                // Access entry of ntuple
                ntuple->GetEntry(evt);

                if (abs(R_L_truth - R_L) < rl_resolution) {
                        hresp_jet_pt->Fill(jet_pt, jet_pt_truth);
                        hresp_rl->Fill(R_L, R_L_truth);
                        hresp_weight->Fill(weight_pt, weight_pt_truth);
                        hresp_ptprod->Fill(h1_pt*h2_pt,h1_pt_truth*h2_pt_truth);
                }
        }

        gStyle->SetPaintTextFormat("1.f");
        hresp_jet_pt->Draw("col");
        hresp_jet_pt->SetTitle("Response matrix of p_{T,jet};p_{T,jet}^{Reco}(GeV);p_{T,jet}^{Truth}(GeV)");
        hresp_jet_pt->Smooth();
        // gPad->SetLogx(1);
        // gPad->SetLogy(1);
        c->Print("./plots/responsematrix_jet_pt.pdf");
        
        hresp_weight->Draw("col");
        hresp_weight->SetTitle("Response matrix of momentum weights;w^{Reco};w^{Truth}");
        hresp_weight->Smooth();
        gPad->SetLogx(1);
        gPad->SetLogy(1);
        c->Print(Form("./plots/responsematrix_weight_nbinweight%i.pdf",nbin_weight));
        
        hresp_rl->GetXaxis()->SetRangeUser(0.01,1);
        hresp_rl->GetYaxis()->SetRangeUser(0.01,1); 
        hresp_rl->Draw("col");
        hresp_rl->SetTitle("Response matrix of R_{L};R_{L}^{Reco};R_{L}^{Truth}");
        hresp_rl->Smooth();
        gPad->SetLogx(1);
        gPad->SetLogy(1);
        c->Print("./plots/responsematrix_rl.pdf");

        hresp_ptprod->GetXaxis()->SetRangeUser(0.05,1000);
        hresp_ptprod->GetYaxis()->SetRangeUser(0.05,1000); 
        hresp_ptprod->Draw("col");
        hresp_ptprod->SetTitle("Response matrix of p_{1}p_{2};p_{1}p_{2}^{Reco};p_{1}p_{2}^{Truth}");
        hresp_ptprod->Smooth();
        gPad->SetLogx(1);
        gPad->SetLogy(1);
        c->Print("./plots/responsematrix_ptprod.pdf");
}