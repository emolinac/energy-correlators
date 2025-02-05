#include "../include/analysis-constants.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"

void macro_print_responsematrix_rl()
{
    TFile* f = new TFile((output_folder+namef_ntuple_e2c_pairpurity).c_str());
    TNtuple* ntuple = (TNtuple*) f->Get(name_ntuple_purity.c_str());

    float R_L, R_L_truth;
    ntuple->SetBranchAddress("R_L",&R_L);
    ntuple->SetBranchAddress("R_L_truth",&R_L_truth);
    
    TH2F* hresp = new TH2F("hresp","",Nbin_R_L+2,em_rlunfolding_binning,Nbin_R_L+2,em_rlunfolding_binning);

    for(int evt = 0 ; evt < ntuple->GetEntries() ; evt++)
    {
        // Access entry of ntuple
        ntuple->GetEntry(evt);

        if(R_L_truth!=-999) hresp->Fill(R_L, R_L_truth);
    }

    hresp->Draw("col text");
    hresp->SetTitle("Response matrix of R_{L};Reco;Truth");
}