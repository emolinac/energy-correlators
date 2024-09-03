#include "../include/analysis-constants.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
//#include "RooUnfoldResponse.h"

void macro_print_responsematrix()
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

    TH1F* hmeas = new TH1F("hmeas","",100,2.5,4);
    TH1F* htrue = new TH1F("htrue","",100,2.5,4);
    TH2F* hresp = new TH2F("hresp","",100,2.5,4,100,2.5,4);

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
    hresponse->Draw("colz");
}