#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"

void macro_print_weight_distributions()
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
        
        TH1F* hresp_rl     = new TH1F("hresp_rl"    ,"",200,0.008,0.8);
        TH1F* hresp_jet_pt = new TH1F("hresp_jet_pt","",200,10,150);
        TH1F* hresp_weight = new TH1F("hresp_weight","",nbin_weight,weight_binning);
        TH1F* hresp_ptprod = new TH1F("hresp_ptprod","",nbin_weight,ptprod_binning);

        for (int evt = 0 ; evt < ntuple->GetEntries() ; evt++) {
                // Access entry of ntuple
                ntuple->GetEntry(evt);

                if (abs(R_L_truth - R_L) < rl_resolution) {
                        hresp_jet_pt->Fill(jet_pt_truth);
                        hresp_rl->Fill(R_L_truth);
                        hresp_weight->Fill(weight_pt_truth);
                        hresp_ptprod->Fill(h1_pt_truth*h2_pt_truth);
                }
        }

        hresp_weight->Draw();
        hresp_weight->SetTitle("Distribution of momentum weights;w^{Reco};w^{Truth}");
        hresp_weight->Smooth();
        gPad->SetLogx(1);
        // gPad->SetLogy(1);
        c->Print(Form("./plots/weight_nbinweight%i.pdf",nbin_weight));
        
        hresp_ptprod->Draw();
        hresp_ptprod->SetTitle("Dsitribution of p_{1}p_{2};p_{1}p_{2}^{Reco};p_{1}p_{2}^{Truth}");
        hresp_ptprod->Smooth();
        gPad->SetLogx(1);
        // gPad->SetLogy(1);
        c->Print("./plots/ptprod.pdf");
}