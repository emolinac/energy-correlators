#include "../include/analysis-constants.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_unfold_invert_eta(double jet_pt_min_local = 20, double jet_pt_max_local = 100)
{
    TFile* f = new TFile((output_folder+namef_ntuple_e2c_purity).c_str());
    TNtuple* ntuple = (TNtuple*) f->Get(name_ntuple_purity.c_str());

    float eta, eta_truth, jet_pt;
    ntuple->SetBranchAddress("jet_pt",&jet_pt);
    ntuple->SetBranchAddress("h_eta",&eta);
    ntuple->SetBranchAddress("htruth_eta",&eta_truth);
    
    // Create histograms with the respective true and matched reco 
    TH1F* hmeas = new TH1F("hmeas","",sl_eta_nbins,eta_min,eta_max);
    TH1F* htrue = new TH1F("htrue","",sl_eta_nbins,eta_min,eta_max);
    TH2F* hresp = new TH2F("hresp","",sl_eta_nbins,eta_min,eta_max,sl_eta_nbins,eta_min,eta_max);

    ntuple->Project("hmeas","h_eta"     ,Form("jet_pt<%f&&jet_pt>%f",jet_pt_max_local,jet_pt_min_local));
    ntuple->Project("htrue","htruth_eta",Form("htruth_eta!=999&&jet_pt<%f&&jet_pt>%f",jet_pt_max_local,jet_pt_min_local));

    // Create response matrix object
    RooUnfoldResponse response(sl_eta_nbins,eta_min,eta_max,sl_eta_nbins,eta_min,eta_max);
    // RooUnfoldResponse response(hmeas, htrue, hresp);

    for(int evt = 0 ; evt < ntuple->GetEntries() ; evt++)
    {
        // Access entry of ntuple
        ntuple->GetEntry(evt);

        if(jet_pt<jet_pt_min_local||jet_pt>jet_pt_max_local) continue;
        if(eta_truth!=-999) response.Fill(eta, eta_truth);
        else if(eta_truth==-999) response.Fake(eta);
    }

    // Draw response matrix
    TH2F* hresponse = (TH2F*) response.HresponseNoOverflow();
    
    RooUnfoldInvert unfold(&response, hmeas);
    auto* hunfolded_bayes = unfold.Hunfold();

    set_histogram_style(hmeas          , kViolet+2, std_line_width, 24              , std_marker_size);
    set_histogram_style(hunfolded_bayes, kViolet+2, std_line_width, std_marker_style, std_marker_size);
    set_histogram_style(htrue          , kGreen   , std_line_width, std_marker_style, std_marker_size);

    TCanvas* c = new TCanvas("c","",1920,1080);
    c->Draw();

    THStack* s = new THStack();
    s->Add(hunfolded_bayes,"E");
    s->Add(hmeas,"E");
    s->Add(htrue,"E");
    s->Draw("NOSTACK");

    gPad->SetLogx(1);
    gPad->SetLogy(1);

    TRatioPlot* rp = new TRatioPlot(htrue, hunfolded_bayes);
    rp->Draw();

    TLegend* l = new TLegend(.15,.75,.25,.9);
    l->AddEntry(hunfolded_bayes,"Unfolded","lpf");
    l->AddEntry(hmeas          ,"Meas."   ,"lpf");
    l->AddEntry(htrue          ,"True"    ,"lpf");
    l->Draw("SAME");

    // c->Print(Form("../plots/unfolding/unfolded_eta_bayes_iter%i.pdf",Niter));
}