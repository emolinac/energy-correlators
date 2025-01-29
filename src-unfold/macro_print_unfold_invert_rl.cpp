#include "../include/analysis-constants.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_unfold_invert_rl()
{
    TFile* f = new TFile((output_folder+namef_ntuple_e2c_pairpurity).c_str());
    TNtuple* ntuple = (TNtuple*) f->Get(name_ntuple_purity.c_str());

    float R_L, R_L_truth, jet_pt;
    ntuple->SetBranchAddress("jet_pt",&jet_pt);
    ntuple->SetBranchAddress("R_L",&R_L);
    ntuple->SetBranchAddress("R_L_truth",&R_L_truth);
    
    // Create histograms with the respective true and matched reco 
    TH1F* hmeas = new TH1F("hmeas","",Nbin_R_L,R_L_min,R_L_max);
    TH1F* htrue = new TH1F("htrue","",Nbin_R_L,R_L_min,R_L_max);
    TH2F* hresp = new TH2F("hresp","",Nbin_R_L,R_L_min,R_L_max,Nbin_R_L,R_L_min,R_L_max);

    ntuple->Project("hmeas","R_L"      ,"jet_pt<30.&&jet_pt>20");
    ntuple->Project("htrue","R_L_truth","R_L_truth!=999&&jet_pt<30.&&jet_pt>20");

    // Create response matrix object
    // RooUnfoldResponse response(Nbin_R_L,R_L_min,R_L_max,Nbin_R_L,R_L_min,R_L_max);
    RooUnfoldResponse response(hmeas, htrue, hresp);

    for(int evt = 0 ; evt < ntuple->GetEntries() ; evt++)
    {
        // Access entry of ntuple
        ntuple->GetEntry(evt);

        if(jet_pt<20||jet_pt>30) continue;

        if(R_L_truth!=-999) response.Fill(R_L, R_L_truth);
        else if(R_L_truth==-999) response.Fake(R_L);
    }

    // Draw response matrix
    TH2F* hresponse = (TH2F*) response.HresponseNoOverflow();
    
    RooUnfoldInvert unfold(&response, hmeas);
    auto* hreco_bayes = unfold.Hunfold();

    set_histogram_style(hmeas      , kViolet+2, std_line_width, 24              , std_marker_size+1);
    set_histogram_style(hreco_bayes, kViolet+2, std_line_width, std_marker_style, std_marker_size+1);
    set_histogram_style(htrue      , kGreen   , std_line_width, std_marker_style, std_marker_size+1);

    TCanvas* c = new TCanvas("c","",1920,1080);
    c->Draw();

    THStack* s = new THStack();
    s->Add(hreco_bayes,"E");
    s->Add(hmeas,"E");
    s->Add(htrue,"E");
    // s->Draw("NOSTACK");
    s->SetTitle(";R_{L};N_{pair}");

    gPad->SetLogx(1);
    gPad->SetLogy(1);

    TRatioPlot* rp = new TRatioPlot(htrue, hreco_bayes);
    rp->Draw();

    TLegend* l = new TLegend(.15,.75,.25,.9);
    l->AddEntry(hreco_bayes,"Unfolded","lpf");
    l->AddEntry(hmeas      ,"Meas."   ,"lpf");
    l->AddEntry(htrue      ,"True"    ,"lpf");
    l->Draw("SAME");

    c->Print("../plots/unfolding/unfolded_rl_invert.pdf");
}