#include "../include/analysis-constants.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_jet_pt_weight()
{
    // Open ROOT file with ntuple
    TFile* f = new TFile((output_folder+namef_ntuple_e2c).c_str());

    // Get the ntuple
    TNtuple* ntuple = (TNtuple*) f->Get((name_ntuple_data).c_str());
    
    // Determine binning
    TH2F* h = new TH2F("h","",100,jet_pt_min, jet_pt_max, 100, 0, .2);
    
    // Create Canvas and draw in it
    TCanvas* c = new TCanvas("","",800,600);
    c->Draw();
    ntuple->Draw("weight:jet_pt>>h",pair_data_noneutrals_cut,"colz");
    
    h->SetTitle(";Jet P_{T};Weight");

    c->Print("../plots/phase_space_jetpt_weight.pdf");
}