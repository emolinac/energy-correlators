#include "../../include/analysis-constants.h"
#include "../../include/analysis-cuts.h"
#include "../../include/directories.h"
#include "../../include/names.h"
#include "../../include/utils-algorithms.h"
#include "../../include/utils-visual.h"

void macro_print_singlepurity_eta_momentum()
{
    // Open the necessary files
    TFile* fpurity = new TFile((output_folder+namef_ntuple_e2c_purity).c_str());

    // Get the corresponding Ntuples
    TNtuple* ntuple_dtrmatch = (TNtuple*) fpurity->Get((name_ntuple_purity).c_str());

    // Define the necessary histograms to calculate purity
    TH2F* hsig    = new TH2F("hsig"   ,"",ic_p_nbins,ic_p_binning,12,2,4.5);
    TH2F* hall    = new TH2F("hall"   ,"",ic_p_nbins,ic_p_binning,12,2,4.5);
    TH2F* hpurity = new TH2F("hpurity","",ic_p_nbins,ic_p_binning,12,2,4.5);
    hsig->Sumw2();
    hall->Sumw2();
    hpurity->Sumw2();

    set_histogram_style(hsig, kViolet, std_line_width, std_marker_style, std_marker_size);
    set_histogram_style(hall, kCyan  , std_line_width, std_marker_style, std_marker_size);

    // Project into the histograms
    ntuple_dtrmatch->Project("hsig","h_eta:h_p",single_signal_cut);
    ntuple_dtrmatch->Project("hall","h_eta:h_p",pair_cut         );
    
    TCanvas* c = new TCanvas("c","",800,600);
    c->Draw();

    // PURITY PLOTS
    hpurity->Divide(hsig,hall,1,1,"B");
    
    //set_histogram_style(hpurity, kViolet, std_line_width, std_marker_style, std_marker_size);
    
    hpurity->Draw("colz");
    //hpurity->GetYaxis()->SetRangeUser(0,1);
    hpurity->SetTitle(";p (GeV);#eta");

    gPad->SetLogx(1);

    //c->Print("../../plots/purity/singlehadron_purity_eta_momentum.pdf");
}