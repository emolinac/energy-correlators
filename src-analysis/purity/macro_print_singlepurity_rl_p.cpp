#include "../include/analysis-constants.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_singlepurity_rl_p(bool include_neutrals = 0)
{
    // Open the necessary files
    TFile* fpurity = new TFile((output_folder+"ntuple_singlepurity.root").c_str());

    // Get the corresponding Ntuples
    TNtuple* ntuple_dtrmatch = (TNtuple*) fpurity->Get((name_ntuple_purity).c_str());

    // Define the necessary histograms to calculate purity
    TH2F* hsig    = new TH2F("hsig"   ,"",20,4,1000,10,eta_min,eta_max);
    TH2F* hall    = new TH2F("hall"   ,"",20,4,1000,10,eta_min,eta_max);
    TH2F* hpurity = new TH2F("hpurity","",20,4,1000,10,eta_min,eta_max);
    
    // Project into the histograms
    ntuple_dtrmatch->Project("hsig","h2_eta:h2_p",single_signal_cut);
    ntuple_dtrmatch->Project("hall","h2_eta:h2_p",pair_cut);
    
    TCanvas* c = new TCanvas("c","",800,600);
    c->Draw();

    // PURITY PLOTS
    hpurity->Divide(hsig,hall,1,1,"B");
    
    hpurity->Draw("COLZ");
}