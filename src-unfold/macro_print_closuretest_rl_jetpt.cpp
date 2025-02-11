#include "../include/analysis-constants.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_closuretest_rl_jetpt(int Niter = 5, double jet_pt_min_local = unfolding_jetpt_binning[0], double jet_pt_max_local = unfolding_jetpt_binning[5])
{
    TFile* f = new TFile((output_folder+namef_ntuple_e2c_pairpurity).c_str());
    TNtuple* ntuple = (TNtuple*) f->Get(name_ntuple_purity.c_str());

    float R_L, R_L_truth, jet_pt, jet_pt_truth;
    ntuple->SetBranchAddress("jet_pt",&jet_pt);
    ntuple->SetBranchAddress("jet_pt_truth",&jet_pt_truth);
    ntuple->SetBranchAddress("R_L",&R_L);
    ntuple->SetBranchAddress("R_L_truth",&R_L_truth);
    
    // Create histograms with the respective true and matched reco 
    TH2F* hmeas = new TH2F("hmeas","",Nbin_R_L+2,unfolding_rl_binning,Nbin_jet_pt+2,unfolding_jetpt_binning);
    TH2F* htrue = new TH2F("htrue","",Nbin_R_L+2,unfolding_rl_binning,Nbin_jet_pt+2,unfolding_jetpt_binning);

    // Create response matrix object
    RooUnfoldResponse* response = new RooUnfoldResponse(hmeas, htrue, "response");

    TH2F* htrue_ref = new TH2F("htrue_ref","",Nbin_R_L+2,unfolding_rl_binning,Nbin_jet_pt+2,unfolding_jetpt_binning);
    TH2F* h_ct      = new TH2F("h_ct"     ,"",Nbin_R_L+2,unfolding_rl_binning,Nbin_jet_pt+2,unfolding_jetpt_binning);

    TRandom3* rndm = new TRandom3();
    for(int evt = 0 ; evt < ntuple->GetEntries() ; evt++)
    {
        // Access entry of ntuple
        ntuple->GetEntry(evt);
        if(jet_pt<jet_pt_min_local||jet_pt>jet_pt_max_local) continue;
        if(R_L_truth==-999) continue;
        if(rndm->Uniform(1)<=0.5) 
        {
            htrue_ref->Fill(R_L_truth,jet_pt_truth);
            continue;
        }

        response->Fill(R_L, jet_pt, R_L_truth, jet_pt_truth);
        hmeas->Fill(R_L, jet_pt);
    }

    TCanvas* c = new TCanvas("c","",1920,1080);
    c->Draw();

    RooUnfoldBayes unfold(response, hmeas, Niter);
    auto* hunfolded_bayes = unfold.Hunfold();
    

    htrue_ref->Scale(1./htrue_ref->Integral());
    hunfolded_bayes->Scale(1./hunfolded_bayes->Integral());

    h_ct->Divide(hunfolded_bayes,htrue_ref,1,1,"B");
    // hmeas->Divide(htrue);

    set_histogram_style(h_ct, kViolet, std_line_width, std_marker_style, std_marker_size+1);
    // set_histogram_style(hmeas, kCyan, std_line_width, std_marker_style, std_marker_size+1);

    h_ct->Draw("COL TEXT");
    gStyle->SetPaintTextFormat(".3f");
    h_ct->SetMarkerSize(.7);
    // h_ct->Smooth();
    h_ct->SetTitle(";R_{L};p^{jet}_{t}");
    
    gPad->SetLogy(1);
    // c->Print(Form("./plots/closuretest_roounfold_rl_bayes_iter%i_jetptfrom%.0fto%.0f.pdf",Niter,jet_pt_min_local,jet_pt_max_local));    
}