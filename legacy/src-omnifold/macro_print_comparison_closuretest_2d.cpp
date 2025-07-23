#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_comparison_closuretest_2d(int Niter = 5, double jet_pt_min_local = unfolding_jet_pt_binning[0], double jet_pt_max_local = unfolding_jet_pt_binning[5])
{
    TFile* fof = new TFile("../output-files/multifold_closuretest_2d.root");
    if (fof->IsZombie()) return;
    TH1D* hct_rl_of    = (TH1D*) fof->Get("hct_rl");
    TH1D* hct_jet_pt_of = (TH1D*) fof->Get("hct_jet_pt");

    TFile* f = new TFile((output_folder + namef_ntuple_e2c_paircorrections).c_str());
    TNtuple* ntuple = (TNtuple*) f->Get(name_ntuple_correction_reco.c_str());

    float R_L, R_L_truth, jet_pt, jet_pt_truth;
    ntuple->SetBranchAddress("jet_pt",&jet_pt);
    ntuple->SetBranchAddress("jet_pt_truth",&jet_pt_truth);
    ntuple->SetBranchAddress("R_L",&R_L);
    ntuple->SetBranchAddress("R_L_truth",&R_L_truth);
    
    // Create histograms with the respective true and matched reco 
    TH2D* hmeas = new TH2D("hmeas","",nbin_rl+2,unfolding_rl_binning,nbin_jet_pt+2,unfolding_jet_pt_binning);
    TH2D* htrue = new TH2D("htrue","",nbin_rl+2,unfolding_rl_binning,nbin_jet_pt+2,unfolding_jet_pt_binning);
    RooUnfoldResponse* response = new RooUnfoldResponse(hmeas, htrue, "response");

    TH2D* htrue_ref = new TH2D("htrue_ref","",nbin_rl+2,unfolding_rl_binning,nbin_jet_pt+2,unfolding_jet_pt_binning);
    TH2D* hct       = new TH2D("hct"     ,"",nbin_rl+2,unfolding_rl_binning,nbin_jet_pt+2,unfolding_jet_pt_binning);

    TRandom3* rndm = new TRandom3();
    for (int evt = 0 ; evt < ntuple->GetEntries() ; evt++)
    {
        // Access entry of ntuple
        ntuple->GetEntry(evt);
        if (jet_pt<jet_pt_min_local||jet_pt>jet_pt_max_local) continue;
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
    c->Draw();

    RooUnfoldBayes unfold(response, hmeas, Niter);
    TH2D* hunfolded_bayes = (TH2D*) unfold.Hunfold();
    htrue_ref->Scale(1./htrue_ref->Integral());
    hunfolded_bayes->Scale(1./hunfolded_bayes->Integral());

    TH1D* htrue_ref_rl          = new TH1D("htrue_ref_rl","",nbin_rl+2,unfolding_rl_binning);
    TH1D* htrue_ref_jet_pt       = new TH1D("htrue_ref_jet_pt","",nbin_jet_pt+2,unfolding_jet_pt_binning);
    TH1D* hunfolded_bayes_rl    = new TH1D("hunfolded_bayes_rl","",nbin_rl+2,unfolding_rl_binning);
    TH1D* hunfolded_bayes_jet_pt = new TH1D("hunfolded_bayes_jet_pt","",nbin_jet_pt+2,unfolding_jet_pt_binning);
    TH1D* hct_rl                = new TH1D("hct_rl","",nbin_rl+2,unfolding_rl_binning);
    TH1D* hct_jet_pt             = new TH1D("hct_jet_pt","",nbin_jet_pt+2,unfolding_jet_pt_binning);
    
    htrue_ref_rl->Sumw2();
    htrue_ref_jet_pt->Sumw2();
    hunfolded_bayes_rl->Sumw2();
    hunfolded_bayes_jet_pt->Sumw2();

    htrue_ref_rl          = htrue_ref->ProjectionX("htrue_ref_rl");
    htrue_ref_jet_pt       = htrue_ref->ProjectionY("htrue_ref_jet_pt");
    hunfolded_bayes_rl    = hunfolded_bayes->ProjectionX("hunfolded_bayes_rl");
    hunfolded_bayes_jet_pt = hunfolded_bayes->ProjectionY("hunfolded_bayes_jet_pt");
    
    hct_rl->Divide(htrue_ref_rl,hunfolded_bayes_rl,1,1);
    hct_jet_pt->Divide(htrue_ref_jet_pt,hunfolded_bayes_jet_pt,1,1);
    
    set_histogram_style(hct_rl_of   , std_marker_color_jet_pt[0], std_line_width, std_marker_style, std_marker_size+1);
    set_histogram_style(hct_jet_pt_of, std_marker_color_jet_pt[0], std_line_width, std_marker_style, std_marker_size+1);
    set_histogram_style(hct_rl      , std_marker_color_jet_pt[1], std_line_width, std_marker_style, std_marker_size+1);
    set_histogram_style(hct_jet_pt   , std_marker_color_jet_pt[1], std_line_width, std_marker_style, std_marker_size+1);
    
    TLegend* lrl = new TLegend();
    lrl->AddEntry(hct_rl,"RooUnfold","lpf");
    lrl->AddEntry(hct_rl_of,"Multifold","lpf");
    THStack* hs_rl = new THStack();
    hs_rl->Add(hct_rl);
    hs_rl->Add(hct_rl_of);
    hs_rl->Draw("NOSTACK");
    hs_rl->SetTitle(";R_{L};Truth/Unfolded");
    hs_rl->GetXaxis()->SetRangeUser(rl_binning[0],rl_binning[nbin_rl]);
    hs_rl->SetMinimum(0.89);
    hs_rl->SetMaximum(1.11);
    gPad->SetLogx(1);
    lrl->Draw("SAME");
    c->Print(Form("./plots/closuretest_2d_comparison_rl_iter%i_jet_ptfrom%.0fto%.0f.pdf",Niter,jet_pt_min_local,jet_pt_max_local));    
    gPad->SetLogx(0);

    TLegend* ljet_pt = new TLegend();
    ljet_pt->AddEntry(hct_jet_pt,"RooUnfold","lpf");
    ljet_pt->AddEntry(hct_jet_pt_of,"Multifold","lpf");
    THStack* hs_jet_pt = new THStack();
    hs_jet_pt->Add(hct_jet_pt);
    hs_jet_pt->Add(hct_jet_pt_of);
    hs_jet_pt->Draw("NOSTACK");
    hs_jet_pt->SetTitle(";p^{jet}_{t}(GeV);Truth/Unfolded");
    hs_jet_pt->GetXaxis()->SetRangeUser(20,100);
    hs_jet_pt->SetMinimum(0.89);
    hs_jet_pt->SetMaximum(1.11);
    ljet_pt->Draw("SAME");
    c->Print(Form("./plots/closuretest_2d_comparison_jet_pt_iter%i_jet_ptfrom%.0fto%.0f.pdf",Niter,jet_pt_min_local,jet_pt_max_local));    

}