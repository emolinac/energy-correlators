#include "../include/analysis-constants.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/user-algorithms.h"
#include "RooUnfoldResponse.h"

void macro_print_responsematrix()
{
    // Open file with the Unfold Ntuple
    TFile* f = new TFile((mother_folder+namef_ntuple_e2c).c_str()),

    // Get the Ntuple and set the needed branches
    float X_L_true, X_L_meas, X_L_true, X_L_meas,
    TNtuple* ntuple = (TNtuple*) f->Get(name_ntuple_unfold.c_str());

    // Create histograms with the respective truth and matched reco 
    double binning[Nbin_X_L];
    determine_log10binning(Nbin_X_L, X_L_min, X_L_max, binning);

    TH1F* hmeas = new TH1F("hmeas","",Nbin_X_L,binning);
    TH1F* htrue = new TH1F("htrue","",Nbin_X_L,binning);

    // Fill the histos
    ntuple->Project("htrue","X_L_true", "");
    ntuple->Project("hmeas","X_L_reco", "");

    // Create response matrix object
    RooUnfoldResponse response(hmeas, htrue);

    for(int i = 0 ; i < ntuple->GetEntries() ; i++)
    {
        ntuple->GetEntry(i);

        if(X_L_reco!=-999||X_L_reco!=0) response.Fill(X_L_meas, X_L_true);
        else response.Fill(X_L_true);
    }
}