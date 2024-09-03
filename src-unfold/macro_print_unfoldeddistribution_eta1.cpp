#include "../include/analysis-constants.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"
//#include "RooUnfoldResponse.h"

void macro_print_unfoldeddistribution_eta1()
{
    gStyle->SetOptStat();

    // Open file with the Unfold Ntuple
    TFile* f = new TFile((output_folder+namef_ntuple_e2c).c_str());

    // Get the Ntuple and set the needed branches
    float h1_eta_true, h1_eta_meas;
    TNtuple* ntuple = (TNtuple*) f->Get(name_ntuple_unfold.c_str());
    ntuple->SetBranchAddress("h1_eta_true",&h1_eta_true);
    ntuple->SetBranchAddress("h1_eta_reco",&h1_eta_meas);

    // Create histograms with the respective true and matched reco 

    TH1F* hmeas = new TH1F("hmeas","",20,2.5,4);
    TH1F* htrue = new TH1F("htrue","",20,2.5,4);
    TH2F* hresp = new TH2F("hresp","",20,2.5,4,20,2.5,4);

    // Fill the histos
    ntuple->Project("htrue","h1_eta_true","");
    ntuple->Project("hmeas","h1_eta_reco","");

    // Create response matrix object
    RooUnfoldResponse response(hmeas, htrue, hresp);

    for(int i = 0 ; i < ntuple->GetEntries() ; i++)
    {
        ntuple->GetEntry(i);

        if(h1_eta_meas!=0) response.Fill(h1_eta_meas, h1_eta_true);
        else response.Miss(h1_eta_true);
    }

    // Draw response matrix
    TH2F* hresponse = (TH2F*) response.HresponseNoOverflow();
    
    // Invert the response matrix
    RooUnfoldBayes unfold(&response, hmeas, 1);
    TH1* hreco = unfold.Hunfold();
    htrue->Sumw2();
    hmeas->Sumw2();

    // visual
    set_histogram_style(htrue, kCyan  , std_line_width, std_marker_style, std_marker_size+1);
    set_histogram_style(hreco, kViolet, std_line_width, std_marker_style, std_marker_size);
    set_histogram_style(hmeas, kGreen , std_line_width, std_marker_style, std_marker_size);

    htrue->Scale(1./htrue->Integral());
    hreco->Scale(1./hreco->Integral());
    hmeas->Scale(1./hmeas->Integral());    

    THStack* s = new THStack();
    s->Add(htrue);
    s->Add(hreco);
    s->Add(hmeas);

    s->Draw("NOSTACK");

    s->SetTitle(";#eta_{1};Norm. Distributions");
    gPad->SetLogy(1);

    TLegend* l = new TLegend();
    l->AddEntry(htrue, "Truth", "lpf");
    l->AddEntry(hmeas, "Meas" , "lpf");
    l->AddEntry(hreco, "Unf"  , "lpf");

    l->Draw("same");
}