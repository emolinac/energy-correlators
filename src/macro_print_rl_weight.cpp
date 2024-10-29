#include "../include/analysis-constants.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_rl_weight()
{
    // Open ROOT file with ntuple
    TFile* f = new TFile((output_folder+namef_ntuple_e2c).c_str());

    // Get the ntuple
    TNtuple* ntuple = (TNtuple*) f->Get((name_ntuple_data).c_str());
    
    // Determine binning
    double binning[Nbin_R_L+1];
    determine_log10binning(Nbin_R_L, R_L_min, R_L_max, binning);

    TH2F* h = new TH2F("h","",Nbin_R_L,R_L_min,R_L_max, 100, 0, .1);
    
    // Create Canvas and draw in it
    TCanvas* c = new TCanvas("","",800,600);
    c->Draw();
    ntuple->Draw("weight:R_L>>h",pair_data_noneutrals_cut,"colz");
    
    h->SetTitle(";R_{L};Weight");

    gPad->SetLogx(1);

    c->Print("../plots/phase_space_rl_weight.pdf");
}