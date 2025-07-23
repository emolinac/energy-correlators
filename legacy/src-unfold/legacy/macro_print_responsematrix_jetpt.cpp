#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"

void macro_print_responsematrix_jet_pt()
{
    TFile* f = new TFile((output_folder + namef_ntuple_e2c_paircorrections).c_str());
    TNtuple* ntuple = (TNtuple*) f->Get(name_ntuple_correction_reco.c_str());

    float jet_pt, jet_pt_truth;
    ntuple->SetBranchAddress("jet_pt",&jet_pt);
    ntuple->SetBranchAddress("jet_pt_truth",&jet_pt_truth);
    
    TH2F* hresp = new TH2F("hresp","",nbin_jet_pt+2,unfolding_jet_pt_binning,nbin_jet_pt+2,unfolding_jet_pt_binning);

    for (int evt = 0 ; evt < ntuple->GetEntries() ; evt++)
    {
        // Access entry of ntuple
        ntuple->GetEntry(evt);
        if (jet_pt_truth!=-999) hresp->Fill(jet_pt, jet_pt_truth);
    }

    hresp->Draw("COL TEXT");
    hresp->SetTitle("Response matrix of p^{jet}_{T};Reco;Truth");
    gPad->SetLogy(1);
    gPad->SetLogx(1);
}