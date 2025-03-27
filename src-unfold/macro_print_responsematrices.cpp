#include "../include/analysis-constants.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"

void macro_print_responsematrices()
{
    TCanvas* c = new TCanvas("","",1920,1080);
    c->Draw();

    TFile* f = new TFile((output_folder+namef_ntuple_e2c_pairpurity).c_str());
    TNtuple* ntuple = (TNtuple*) f->Get(name_ntuple_purity.c_str());

    float jet_pt, jet_pt_truth, R_L, R_L_truth, weight, weight_truth;
    ntuple->SetBranchAddress("jet_pt",&jet_pt);
    ntuple->SetBranchAddress("jet_pt_truth",&jet_pt_truth);
    ntuple->SetBranchAddress("R_L",&R_L);
    ntuple->SetBranchAddress("R_L_truth",&R_L_truth);
    ntuple->SetBranchAddress("weight_pt",&weight);
    ntuple->SetBranchAddress("weight_pt_truth",&weight_truth);
    
    TH2F* hresp_rl     = new TH2F("hresp_rl","",Nbin_R_L+2,unfolding_rl_binning,Nbin_R_L+2,unfolding_rl_binning);
    TH2F* hresp_jetpt  = new TH2F("hresp_jetpt","",Nbin_jet_pt+2,unfolding_jetpt_binning,Nbin_jet_pt+2,unfolding_jetpt_binning);
    TH2F* hresp_weight = new TH2F("hresp_weight","",Nbin_weight,weight_binning,Nbin_weight,weight_binning);

    for(int evt = 0 ; evt < ntuple->GetEntries() ; evt++)
    {
        // Access entry of ntuple
        ntuple->GetEntry(evt);

        if(R_L_truth!=-999) 
        {
            hresp_jetpt->Fill(jet_pt, jet_pt_truth);
            hresp_rl->Fill(R_L, R_L_truth);
            hresp_weight->Fill(weight, weight_truth);
        }
    }

    hresp_jetpt->Draw("col text");
    hresp_jetpt->SetTitle("Response matrix of p^{jet}_{t};Reco;Truth");
    gPad->SetLogx(1);
    gPad->SetLogy(1);
    c->Print("./plots/responsematrix_jetpt.pdf");
    hresp_rl->Draw("col text");
    hresp_rl->SetTitle("Response matrix of R_{L};Reco;Truth");
    gPad->SetLogx(1);
    gPad->SetLogy(1);
    c->Print("./plots/responsematrix_rl.pdf");
    hresp_weight->Draw("col text");
    hresp_weight->SetTitle("Response matrix of momentum weights;Reco;Truth");
    gPad->SetLogx(1);
    gPad->SetLogy(1);
    c->Print("./plots/responsematrix_weight.pdf");
}