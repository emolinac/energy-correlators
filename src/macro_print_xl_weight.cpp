#include "../include/analysis-constants.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_xl_weight()
{
    // Open ROOT file with ntuple
    TFile* f = new TFile((output_folder+namef_ntuple_e2c).c_str());

    // Get the ntuple
    TNtuple* ntuple = (TNtuple*) f->Get((name_ntuple_data).c_str());
    
    // Determine binning
    double binning[Nbin_X_L+1];
    determine_log10binning(Nbin_X_L, X_L_min, X_L_max, binning);

    TH2F* h = new TH2F("h","",Nbin_X_L, binning, 100, 0, .1);
    
    // Create Canvas and draw in it
    TCanvas* c = new TCanvas("","",800,600);
    c->Draw();
    ntuple->Draw("weight:X_L>>h",e2c_cut,"colz");
    
    h->SetTitle(";X_{L};Weight");

    gPad->SetLogx(1);

    c->Print("../plots/XL_weight.pdf");
}