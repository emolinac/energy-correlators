#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.cpp"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils.cpp"
#include "../include/utils.h"
#include "../include/utils-visual.cpp"
#include "../include/utils-visual.h"

void macro_print_responsematrices()
{
        TCanvas* c = new TCanvas("","",1080,720);
        c->Draw();

        TFile* f = new TFile((output_folder + namef_ntuple_reco2truth_match).c_str());
        TNtuple* ntuple     = (TNtuple*) f->Get(name_ntuple_correction_reco.c_str());
        TNtuple* ntuple_jet = (TNtuple*) f->Get(name_ntuple_jet_reco2truth_match.c_str());

        float jet_pt, jet_pt_truth, R_L, R_L_truth, weight_pt, weight_pt_truth, h1_pt, h1_pt_truth, h2_pt, h2_pt_truth;
        set_unfolding_ntuple_branches(ntuple, &R_L, &R_L_truth, &jet_pt, &jet_pt_truth, &weight_pt, &weight_pt_truth);
        ntuple->SetBranchAddress("h1_pt",&h1_pt);
        ntuple->SetBranchAddress("h1_pt_truth",&h1_pt_truth);
        ntuple->SetBranchAddress("h2_pt",&h2_pt);
        ntuple->SetBranchAddress("h2_pt_truth",&h2_pt_truth);

        float jet_pt_reco_jetunf, jet_pt_truth_jetunf;
        ntuple_jet->SetBranchAddress("jet_pt",&jet_pt_reco_jetunf);
        ntuple_jet->SetBranchAddress("jet_pt_truth",&jet_pt_truth_jetunf);
        
        TH2F* hresp_rl     = new TH2F("hresp_rl"    ,"",nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning,nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning);
        // TH2F* hresp_jet_pt = new TH2F("hresp_jet_pt","",200,10,150,200,10,150);
        TH2F* hresp_jet_pt = new TH2F("hresp_jet_pt","",nbin_jet_pt_unfolding,unfolding_jet_pt_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        TH2F* hresp_weight = new TH2F("hresp_weight","",nbin_weight,weight_binning,nbin_weight,weight_binning);
        TH2F* hresp_ptprod = new TH2F("hresp_ptprod","",nbin_ptprod,ptprod_binning,nbin_ptprod,ptprod_binning);

        for (int evt = 0 ; evt < ntuple->GetEntries() ; evt++) {
                // Access entry of ntuple
                ntuple->GetEntry(evt);

                if (abs(R_L_truth - R_L) < rl_resolution) {
                        // hresp_jet_pt->Fill(jet_pt, jet_pt_truth);
                        hresp_rl->Fill(R_L, R_L_truth);
                        hresp_weight->Fill(weight_pt, weight_pt_truth);
                        hresp_ptprod->Fill(h1_pt*h2_pt,h1_pt_truth*h2_pt_truth);
                }
        }

        for (int evt = 0 ; evt < ntuple_jet->GetEntries() ; evt++) {
                ntuple_jet->GetEntry(evt);

                hresp_jet_pt->Fill(jet_pt_reco_jetunf, jet_pt_truth_jetunf);
        }

        gStyle->SetPaintTextFormat("1.f");
        hresp_jet_pt->Draw("colz");
        hresp_jet_pt->SetTitle("Response matrix of p_{T,jet};p_{T,jet}^{reco}(GeV);p_{T,jet}^{truth}(GeV)");
        gPad->SetRightMargin(0.15);
        // gPad->SetGridx(1);
        // gPad->SetGridy(1);
        gPad->SetLogx(1);
        gPad->SetLogy(1);
        c->Print("./plots/responsematrix_jet_pt.pdf");
        gStyle->SetOptStat(0);
        
        hresp_weight->Draw("colz");
        hresp_weight->SetTitle("Response matrix of momentum weights;w^{reco};w^{truth}");
        gPad->SetRightMargin(0.1);
        gPad->SetLogx(1);
        gPad->SetLogy(1);
        c->Print(Form("./plots/responsematrix_weight_nbinweight%i.pdf",nbin_weight));
        
        hresp_rl->Draw("colz");
        hresp_rl->GetXaxis()->SetRangeUser(0.01,1);
        hresp_rl->GetYaxis()->SetRangeUser(0.01,1); 
        hresp_rl->SetTitle("Response matrix of R_{L};R_{L}^{reco};R_{L}^{truth}");
        // hresp_rl->Smooth();
        gPad->SetLogx(1);
        gPad->SetLogy(1);
        c->Print("./plots/responsematrix_rl.pdf");

        hresp_ptprod->GetXaxis()->SetRangeUser(0.05,1000);
        hresp_ptprod->GetYaxis()->SetRangeUser(0.05,1000); 
        hresp_ptprod->Draw("col text");
        hresp_ptprod->SetTitle("Response matrix of p_{1}p_{2};p_{1}p_{2}^{reco};p_{1}p_{2}^{truth}");
        // hresp_ptprod->Smooth();
        gPad->SetLogx(1);
        gPad->SetLogy(1);
        c->Print("./plots/responsematrix_ptprod.pdf");
}