#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_roounfold_closuretest_2d(int Niter = 4)
{
    TFile* f = new TFile((output_folder+namef_ntuple_e2c_paircorrections).c_str());
    TNtuple* ntuple = (TNtuple*) f->Get(name_ntuple_correction_reco.c_str());

    float R_L, R_L_truth, jet_pt, jet_pt_truth;
    ntuple->SetBranchAddress("jet_pt",&jet_pt);
    ntuple->SetBranchAddress("jet_pt_truth",&jet_pt_truth);
    ntuple->SetBranchAddress("R_L",&R_L);
    ntuple->SetBranchAddress("R_L_truth",&R_L_truth);
    
    // Create histograms with the respective true and matched reco 
    TH2D* hmeas = new TH2D("hmeas","",Nbin_R_L_logbin_unfolding,unfolding_rl_logbinning,Nbin_jet_pt_unfolding,unfolding_jetpt_binning);
    TH2D* htrue = new TH2D("htrue","",Nbin_R_L_logbin_unfolding,unfolding_rl_logbinning,Nbin_jet_pt_unfolding,unfolding_jetpt_binning);
    RooUnfoldResponse* response = new RooUnfoldResponse(hmeas, htrue, "response");

    TH2D* htrue_ref = new TH2D("htrue_ref","",Nbin_R_L_logbin_unfolding,unfolding_rl_logbinning,Nbin_jet_pt_unfolding,unfolding_jetpt_binning);
    TH2D* hct       = new TH2D("hct"      ,"",Nbin_R_L_logbin_unfolding,unfolding_rl_logbinning,Nbin_jet_pt_unfolding,unfolding_jetpt_binning);

    TRandom3* rndm = new TRandom3(0);
    for (int evt = 0 ; evt < ntuple->GetEntries() ; evt++)
    {
        // Access entry of ntuple
        ntuple->GetEntry(evt);
        if (abs(R_L_truth-R_L_reco)>0.015) continue;
        if (rndm->Uniform(1)<=0.25) 
        {
            htrue_ref->Fill(R_L_truth,jet_pt_truth);
            continue;
        }
        else if (rndm->Uniform(1)<=0.5&&rndm->Uniform(1)>=0.25)
        {
            hmeas->Fill(R_L, jet_pt);
            continue;
        }

        response->Fill(R_L, jet_pt, R_L_truth, jet_pt_truth);
        // hmeas->Fill(R_L, jet_pt);
    }

    // TCanvas* c = new TCanvas("c","",1920,1080);
    TCanvas* c2d = new TCanvas("c2d","",1920,1080);
    // c->Draw();
    // c->cd();

    RooUnfoldBayes unfold(response, hmeas, Niter);
    TH2D* hunfolded_bayes = (TH2D*) unfold.Hunfold();
    htrue_ref->Scale(1./htrue_ref->Integral());
    hunfolded_bayes->Scale(1./hunfolded_bayes->Integral());

    TH1D* htrue_ref_rl          = new TH1D("htrue_ref_rl"         ,"",Nbin_R_L_logbin_unfolding,unfolding_rl_logbinning);
    TH1D* htrue_ref_jetpt       = new TH1D("htrue_ref_jetpt"      ,"",Nbin_jet_pt_unfolding    ,unfolding_jetpt_binning);
    TH1D* hunfolded_bayes_rl    = new TH1D("hunfolded_bayes_rl"   ,"",Nbin_R_L_logbin_unfolding,unfolding_rl_logbinning);
    TH1D* hunfolded_bayes_jetpt = new TH1D("hunfolded_bayes_jetpt","",Nbin_jet_pt_unfolding    ,unfolding_jetpt_binning);
    TH1D* hct_rl                = new TH1D("hct_rl"               ,"",Nbin_R_L_logbin_unfolding,unfolding_rl_logbinning);
    TH1D* hct_jetpt             = new TH1D("hct_jetpt"            ,"",Nbin_jet_pt_unfolding    ,unfolding_jetpt_binning);
    
    TH2D* hct_rl_jetpt          = new TH2D("hct_rl_jetpt","",Nbin_R_L_logbin_unfolding,unfolding_rl_logbinning,Nbin_jet_pt_unfolding,unfolding_jetpt_binning);
    
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
    
    set_histogram_style(hct_rl    , std_marker_color_jet_pt[1], std_line_width, std_marker_style, std_marker_size+1);
    set_histogram_style(hct_jetpt , std_marker_color_jet_pt[1], std_line_width, std_marker_style, std_marker_size+1);
    
    c2d->Draw();
    c2d->cd();
    hct_rl_jetpt->Draw("COL");

    // Adding content with errors
    TLatex latex;
    latex.SetTextAlign(22); // center alignment
    latex.SetTextSize(0.015);
    latex.SetTextColor(kBlack);

    for (int i = 1; i <= hct_rl_jetpt->GetNbinsX(); ++i) {
        for (int j = 2; j <= hct_rl_jetpt->GetNbinsY(); ++j) {
            double x = hct_rl_jetpt->GetXaxis()->GetBinCenter(i);
            double y = hct_rl_jetpt->GetYaxis()->GetBinCenter(j);
            double content = hct_rl_jetpt->GetBinContent(i, j);
            double error = hct_rl_jetpt->GetBinError(i, j);
            // Draw content and error in the format "content ± error"
            latex.DrawLatex(x, y, Form("%.2f #pm %.2f", content, error));
        }
    }

    hct_rl_jetpt->SetTitle(";R_{L};p^{jet}_{T} GeV");
    hct_rl_jetpt->GetYaxis()->SetRangeUser(20,100);
    hct_rl_jetpt->GetXaxis()->SetRangeUser(rl_logbinning[0],rl_logbinning[Nbin_R_L]);
    gStyle->SetPaintTextFormat(".2f");
    gPad->SetLogx(1);
    gPad->SetLogy(1);
    
    c2d->Print(Form("./plots/unfolded2d_closuretest_niter%i.pdf",Niter));    
}