#include "../../include/analysis-constants.h"
#include "../../include/analysis-cuts.h"
#include "../../include/directories.h"
#include "../../include/names.h"
#include "../../include/utils-algorithms.h"
#include "../../include/utils-visual.h"

void macro_print_singlepurity_rl()
{
    // Open the necessary files
    TFile* fpurity = new TFile((output_folder+"ntuple_singlepurity.root").c_str());

    // Get the corresponding Ntuples
    TNtuple* ntuple_dtrmatch = (TNtuple*) fpurity->Get((name_ntuple_purity).c_str());

    // Define the necessary histograms to calculate purity
    TH1F* hsig    = new TH1F("hsig"   ,"",10,2,4.5);
    TH1F* hall    = new TH1F("hall"   ,"",10,2,4.5);
    TH1F* hpurity = new TH1F("hpurity","",10,2,4.5);
    hsig->Sumw2();
    hall->Sumw2();
    hpurity->Sumw2();

    set_histogram_style(hsig, kViolet, std_line_width, std_marker_style, std_marker_size);
    set_histogram_style(hall, kCyan  , std_line_width, std_marker_style, std_marker_size);

    // Project into the histograms
    ntuple_dtrmatch->Project("hsig","h2_eta",single_signal_cut);
    ntuple_dtrmatch->Project("hall","h2_eta",pair_cut         );
    
    TCanvas* c = new TCanvas("c","",800,600);
    c->Draw();

    // PURITY PLOTS
    hpurity->Divide(hsig,hall,1,1,"B");
    
    set_histogram_style(hpurity, kViolet, std_line_width, std_marker_style, std_marker_size);
    
    hpurity->Draw();
    hpurity->GetXaxis()->SetRangeUser(eta_min,eta_max);
    hpurity->SetTitle(Form("#Delta R_{L}(truth-reco)<%.3f;R_{L};single Purity",R_L_res));

    //c->Print(Form("../../plots/purity/nsingle_purity_eta_deltarleq%.3f.pdf",R_L_res));
}