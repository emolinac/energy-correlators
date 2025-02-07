#include "../../include/analysis-constants.h"
#include "../../include/analysis-cuts.h"
#include "../../include/directories.h"
#include "../../include/names.h"
#include "../../include/utils-algorithms.h"
#include "../../include/utils-visual.h"

void macro_print_singlepurity_phi()
{
    // Open the necessary files
    TFile* fpurity = new TFile((output_folder+namef_ntuple_e2c_purity).c_str());

    // Get the corresponding Ntuples
    TNtuple* ntuple_dtrmatch = (TNtuple*) fpurity->Get((name_ntuple_purity).c_str());

    // Define the necessary histograms to calculate purity
    TH1F* hsig    = new TH1F("hsig"   ,"",ndim_corr,-TMath::Pi(),TMath::Pi());
    TH1F* hall    = new TH1F("hall"   ,"",ndim_corr,-TMath::Pi(),TMath::Pi());
    TH1F* hpurity = new TH1F("hpurity","",ndim_corr,-TMath::Pi(),TMath::Pi());
    hsig->Sumw2();
    hall->Sumw2();
    hpurity->Sumw2();

    set_histogram_style(hsig, kViolet, std_line_width, std_marker_style, std_marker_size);
    set_histogram_style(hall, kCyan  , std_line_width, std_marker_style, std_marker_size);

    // Project into the histograms
    ntuple_dtrmatch->Project("hsig","h_phi",single_signal_cut);
    ntuple_dtrmatch->Project("hall","h_phi",pair_cut         );
    
    TCanvas* c = new TCanvas("c","",800,600);
    c->Draw();

    // PURITY PLOTS
    hpurity->Divide(hsig,hall,1,1,"B");
    
    set_histogram_style(hpurity, kViolet, std_line_width, std_marker_style, std_marker_size);
    
    hpurity->Draw();
    hpurity->GetYaxis()->SetRangeUser(0,1);
    hpurity->SetTitle(";#phi;Single Hadron Purity");

    c->Print("../../plots/purity/singlehadron_purity_phi.pdf");
}