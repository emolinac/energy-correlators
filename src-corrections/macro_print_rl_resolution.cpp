#include "../include/analysis-constants.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_rl_resolution(bool include_neutrals = 0)
{
    gStyle->SetOptStat(1100);

    // Open the necessary files
    TFile* fpurity = new TFile((output_folder+namef_ntuple_e2c_purity).c_str());

    // Get the corresponding Ntuples
    TNtuple* ntuple_dtrmatch = (TNtuple*) fpurity->Get((name_ntuple_purity).c_str());

    // Define the necessary histograms to calculate purity
    TH1F* hres = new TH1F("hres","",250,-.05,.05);
    hres->Sumw2();
    set_histogram_style(hres, kViolet, std_line_width, std_marker_style, std_marker_size);
    
    // Project into the histograms
    if(include_neutrals) ntuple_dtrmatch->Project("hres","R_L_truth-R_L",pair_all_cut+"R_L!=0&&R_L_truth!=0");
    else ntuple_dtrmatch->Project("hres","R_L_truth-R_L",pair_all_noneutrals_cut+"R_L!=0&&R_L_truth!=0");
    
    TCanvas* c = new TCanvas("c","",800,600);
    c->Draw();

    // MCRECO PLOTS
    hres->Draw();

    hres->SetTitle(";R_{L}_{truth-reco};");

    //gPad->SetLogy(1);

    if(include_neutrals) c->Print("../plots/RL_resolution.pdf");
    else c->Print("../plots/RL_resolution_noneutrals.pdf");
}