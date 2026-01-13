#include "../../include/analysis-constants.h"
#include "../../include/analysis-binning.h"
#include "../../include/analysis-cuts.cpp"
#include "../../include/analysis-cuts.h"
#include "../../include/directories.h"
#include "../../include/names.h"
#include "../../include/utils.cpp"
#include "../../include/utils.h"
#include "../../include/utils-visual.cpp"
#include "../../include/utils-visual.h"

void macro_print_singlepurity_momentum()
{
    // Open the necessary files
    TFile* fpurity = new TFile((output_folder + namef_ntuple_reco2truth_singlehadron_match).c_str());

    // Get the corresponding Ntuples
    TNtuple* ntuple_dtrmatch = (TNtuple*) fpurity->Get((name_ntuple_correction_reco).c_str());

    // Define the necessary histograms to calculate purity
    TH1F* hsig    = new TH1F("hsig"   ,"",10,track_p_min,track_p_max);
    TH1F* hall    = new TH1F("hall"   ,"",10,track_p_min,track_p_max);
    TH1F* hpurity = new TH1F("hpurity","",10,track_p_min,track_p_max);
    hsig->Sumw2();
    hall->Sumw2();
    hpurity->Sumw2();

    set_histogram_style(hsig, kViolet, std_line_width, std_marker_style, std_marker_size);
    set_histogram_style(hall, kCyan  , std_line_width, std_marker_style, std_marker_size);

    // Project into the histograms
    ntuple_dtrmatch->Project("hsig","h_p","key_match==1");
    ntuple_dtrmatch->Project("hall","h_p","");
    
    TCanvas* c = new TCanvas("c","",800,600);
    c->Draw();

    // PURITY PLOTS
    hpurity->Divide(hsig,hall,1,1,"B");
    
    set_histogram_style(hpurity, kViolet, std_line_width, std_marker_style, std_marker_size);
    
    hpurity->Draw();
    hpurity->GetYaxis()->SetRangeUser(0,1.19);
    hpurity->SetTitle(";p(GeV);Single Hadron Purity");

    c->Print("./plots/singlehadron_purity_momentum.pdf");
}