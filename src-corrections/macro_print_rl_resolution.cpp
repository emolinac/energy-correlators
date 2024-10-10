#include "../include/analysis-constants.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_rl_resolution()
{
    gStyle->SetOptStat(1100);

    // Open the necessary files
    TFile* fpurity = new TFile((output_folder+namef_ntuple_e2c_dtrmatch).c_str());

    // Get the corresponding Ntuples
    TNtuple* ntuple_dtrmatch = (TNtuple*) fpurity->Get((name_ntuple_reco2mcdtrmatch).c_str());

    // Define the necessary histograms to calculate purity
    TH1F* hres = new TH1F("hres","",250,-.6,.6);
    hres->Sumw2();
    set_histogram_style(hres, kViolet, std_line_width, std_marker_style, std_marker_size);
    
    // Project into the histograms
    ntuple_dtrmatch->Project("hres","R_L_truth-R_L",pair_all_cut+"(R_L_truth!=-999)");
    
    TCanvas* c = new TCanvas("c","",800,600);
    c->Draw();

    // MCRECO PLOTS
    hres->Draw();

    hres->SetTitle(";R_{L}_{truth-reco};");

    gPad->SetLogy(1);

    c->Print("../plots/RL_resolution.pdf");
}