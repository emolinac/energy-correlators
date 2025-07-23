#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_roounfold_closuretest_2d(int Niter = 4)
{
    TFile* f = new TFile((output_folder + namef_ntuple_e2c_paircorrections).c_str());
    TNtuple* ntuple = (TNtuple*) f->Get(name_ntuple_correction_reco.c_str());

    float R_L, R_L_truth, jet_pt, jet_pt_truth;
    ntuple->SetBranchAddress("jet_pt",&jet_pt);
    ntuple->SetBranchAddress("jet_pt_truth",&jet_pt_truth);
    ntuple->SetBranchAddress("R_L",&R_L);
    ntuple->SetBranchAddress("R_L_truth",&R_L_truth);
    
    // Create histograms with the respective true and matched reco 
    TH2D* hmeas = new TH2D("hmeas","",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jetpt_binning);
    TH2D* htrue = new TH2D("htrue","",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jetpt_binning);
    RooUnfoldResponse* response = new RooUnfoldResponse(hmeas, htrue, "response");

    TH2D* htrue_ref = new TH2D("htrue_ref","",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jetpt_binning);
    TH2D* hct       = new TH2D("hct"      ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jetpt_binning);

    TRandom3* rndm = new TRandom3();
    for (int evt = 0 ; evt < ntuple->GetEntries() ; evt++)
    {
        // Access entry of ntuple
        ntuple->GetEntry(evt);
        if (abs(R_L_truth-R_L_reco)>rl_resolution) continue;
        if (rndm->Uniform(1)<=0.5) 
        {
            htrue_ref->Fill(R_L_truth,jet_pt_truth);
            continue;
        }

        response->Fill(R_L, jet_pt, R_L_truth, jet_pt_truth);
        hmeas->Fill(R_L, jet_pt);
    }

    TCanvas* c = new TCanvas("c", "", 1920, 1080);
    TCanvas* c2d = new TCanvas("c2d","",1920,1080);
    c->Draw();
    c->cd();

    RooUnfoldBayes unfold(response, hmeas, Niter);
    TH2D* hunfolded_bayes = (TH2D*) unfold.Hunfold();
    htrue_ref->Scale(1./htrue_ref->Integral());
    hunfolded_bayes->Scale(1./hunfolded_bayes->Integral());

    TH1D* htrue_ref_rl          = new TH1D("htrue_ref_rl"         ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
    TH1D* htrue_ref_jetpt       = new TH1D("htrue_ref_jetpt"      ,"",nbin_jet_pt_unfolding    ,unfolding_jetpt_binning);
    TH1D* hunfolded_bayes_rl    = new TH1D("hunfolded_bayes_rl"   ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
    TH1D* hunfolded_bayes_jetpt = new TH1D("hunfolded_bayes_jetpt","",nbin_jet_pt_unfolding    ,unfolding_jetpt_binning);
    TH1D* hct_rl                = new TH1D("hct_rl"               ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
    TH1D* hct_jetpt             = new TH1D("hct_jetpt"            ,"",nbin_jet_pt_unfolding    ,unfolding_jetpt_binning);
    
    TH2D* hct_rl_jetpt          = new TH2D("hct_rl_jetpt","",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jetpt_binning);
    
    htrue_ref_rl->Sumw2();
    htrue_ref_jetpt->Sumw2();
    hunfolded_bayes_rl->Sumw2();
    hunfolded_bayes_jetpt->Sumw2();
    
    htrue_ref_rl           = htrue_ref->ProjectionX("htrue_ref_rl");
    htrue_ref_jetpt        = htrue_ref->ProjectionY("htrue_ref_jetpt");
    hunfolded_bayes_rl     = hunfolded_bayes->ProjectionX("hunfolded_bayes_rl");
    hunfolded_bayes_jetpt  = hunfolded_bayes->ProjectionY("hunfolded_bayes_jetpt");
    
    hct_rl->Divide(htrue_ref_rl,hunfolded_bayes_rl,1,1);
    hct_jetpt->Divide(htrue_ref_jetpt,hunfolded_bayes_jetpt,1,1);
    
    hct_rl_jetpt->Divide(htrue_ref,hunfolded_bayes,1,1);
    
    set_histogram_style(hct_rl       , std_marker_color_jet_pt[1], std_line_width, std_marker_style, std_marker_size+1);
    set_histogram_style(hct_jetpt    , std_marker_color_jet_pt[1], std_line_width, std_marker_style, std_marker_size+1);
    
    TLegend* lrl = new TLegend();
    lrl->AddEntry(hct_rl,"RooUnfold","lpf");
    THStack* hs_rl = new THStack();
    hs_rl->Add(hct_rl);
    hs_rl->Draw("NOSTACK");
    hs_rl->SetTitle(";R_{L};Truth/Unfolded");
    // hs_rl->GetXaxis()->SetRangeUser(rl_binning[0],rl_binning[nbin_rl]);
    hs_rl->SetMinimum(0.89);
    hs_rl->SetMaximum(1.11);
    gPad->SetLogx(1);
    lrl->Draw("SAME");
    c->Print(Form("./plots/unfolded2d_closuretest_rl_niter%i.pdf",Niter));    
    gPad->SetLogx(0);

    TLegend* ljetpt = new TLegend();
    ljetpt->AddEntry(hct_jetpt,"RooUnfold","lpf");
    THStack* hs_jetpt = new THStack();
    hs_jetpt->Add(hct_jetpt);
    hs_jetpt->Draw("NOSTACK");
    hs_jetpt->SetTitle(";p^{jet}_{t}(GeV);Truth/Unfolded");
    // hs_jetpt->GetXaxis()->SetRangeUser(20,100);
    hs_jetpt->SetMinimum(0.89);
    hs_jetpt->SetMaximum(1.11);
    ljetpt->Draw("SAME");
    c->Print(Form("./plots/unfolded2d_closuretest_jetpt_niter%i.pdf",Niter));    

    c2d->Draw();
    c2d->cd();
    hct_rl_jetpt->Draw("COL TEXT");
    hct_rl_jetpt->SetTitle(";R_{L};p^{jet}_{T} GeV");
    hct_rl_jetpt->GetYaxis()->SetRangeUser(20,100);
    hct_rl_jetpt->GetYaxis()->SetRangeUser(rl_nominal_binning[0],rl_nominal_binning[nbin_rl]);
    gStyle->SetPaintTextFormat(".2f");
    gPad->SetLogx(1);
    gPad->SetLogy(1);
    
    c2d->Print(Form("./plots/unfolded2d_closuretest_niter%i.pdf",Niter));    
}