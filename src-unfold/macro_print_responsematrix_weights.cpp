#include "../include/analysis-constants.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"

// Macro to create response matrix with custom binning

void macro_print_responsematrix_weights()
{
    TFile* f = new TFile((output_folder+namef_ntuple_e2c_pairpurity).c_str());
    TNtuple* ntuple = (TNtuple*) f->Get(name_ntuple_purity.c_str());

    float weight, weight_truth;
    ntuple->SetBranchAddress("weight",&weight);
    ntuple->SetBranchAddress("weight_truth",&weight_truth);
    
    // Create histograms with the respective true and matched reco 
    TH1F* hmeas = new TH1F("hmeas","",Nbin_weight,weight_binning);
    TH1F* htrue = new TH1F("htrue","",Nbin_weight,weight_binning);
    TH2F* hresp = new TH2F("hresp","",Nbin_weight,weight_binning,Nbin_weight,weight_binning);

    for(int evt = 0 ; evt < ntuple->GetEntries() ; evt++)
    {
        // Access entry of ntuple
        ntuple->GetEntry(evt);

        if(weight_truth!=-999) hresp->Fill(weight, weight_truth);
    }

    // ntuple->Project("hmeas","weight"      ,"");
    // ntuple->Project("htrue","weight_truth","weight_truth!=999");

    // Create response matrix object
    RooUnfoldResponse* response = new RooUnfoldResponse(hmeas,htrue,hresp);

    // Draw response matrix
    TH2F* hresponse = (TH2F*) response->Hresponse();
    hresponse->Draw("col text");
}
