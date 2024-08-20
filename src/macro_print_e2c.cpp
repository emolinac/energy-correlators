#include "../include/analysis-constants.h"
#include "../include/directories.h"
#include "../include/names.h"

void macro_print_e2c()
{
    // Open ROOT file with ntuple
    TFile* f = new TFile((output_folder+namef_ntuple_e2c).c_str());

    // Get the ntuple
    TNtuple* ntuple = (TNtuple*) f->Get((name_ntuple_data).c_str());

    TH1F* h = new TH1F("h","",100,0.001,2);
    ntuple->Draw("X_L>>h",e2c_nominal_cut/*"(weight)"*/);

    // Other way around
//    float weight, X_L;
//    ntuple->SetBranchAddress("weight",&weight);
//    ntuple->SetBranchAddress("X_L",&X_L);
//
//    for(int i = 0 ; i <ntuple->GetEntries() ; i++)
//    {
//        ntuple->GetEntry(i);
//        if(X_L==0||weight<0) continue;
//        h->Fill(X_L,weight);
//    }
//
//    h->Draw();

    gPad->SetLogx(1);
    gPad->SetLogy(1);
}