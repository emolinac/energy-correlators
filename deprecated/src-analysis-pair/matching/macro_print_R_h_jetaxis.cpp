#include "../../include/analysis-constants.h"
#include "../../include/analysis-cuts.h"
#include "../../include/directories.h"
#include "../../include/names.h"
#include "../../include/utils-algorithms.h"
#include "../../include/utils-visual.h"

void macro_print_R_h_jetaxis()
{
    // Open the necessary files
    TFile* fpurity = new TFile((output_folder+namef_ntuple_e2c_purity).c_str());

    // Get the corresponding Ntuples
    TNtuple* ntuple_dtrmatch = (TNtuple*) fpurity->Get((name_ntuple_purity).c_str());

    // Define the necessary histograms to calculate purity
    TH1F* hsig    = new TH1F("hsig"   ,"",100,0,6);
    TH1F* hall    = new TH1F("hall"   ,"",100,0,6);
    hsig->Sumw2();
    hall->Sumw2();
    
    set_histogram_style(hsig, kViolet, std_line_width, std_marker_style, std_marker_size);
    set_histogram_style(hall, kCyan  , std_line_width, std_marker_style, std_marker_size);

    // Project into the histograms
    ntuple_dtrmatch->Project("hsig","R_jet_h",pair_cut+"key_match==1");
    ntuple_dtrmatch->Project("hall","R_jet_h",pair_cut);
    
    TCanvas* c = new TCanvas("c","",800,600);
    c->Draw();

    THStack* s = new THStack();
    s->Add(hsig);
    s->Add(hall);
    s->Draw("NOSTACK");
    s->SetTitle(";R_{L};");
    gPad->SetLogy(1);

    c->Print("./plots/R_h_jet.pdf");
}