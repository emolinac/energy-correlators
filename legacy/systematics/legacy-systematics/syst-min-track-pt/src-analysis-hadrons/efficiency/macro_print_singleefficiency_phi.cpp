#include "../../include/analysis-constants.h"
#include "../../include/analysis-binning.h"
#include "../../include/analysis-cuts.h"
#include "../../include/directories.h"
#include "../../include/names.h"
#include "../../include/utils-algorithms.h"
#include "../../include/utils-visual.h"

void macro_print_singleefficiency_phi()
{
    // Open the necessary files
    TFile* fefficiency = new TFile((output_folder + namef_ntuple_e2c_hadroncorrections).c_str());

    // Get the corresponding Ntuples
    TNtuple* ntuple_mcreco = (TNtuple*) fefficiency->Get((name_ntuple_correction_reco).c_str());
    TNtuple* ntuple_mc     = (TNtuple*) fefficiency->Get((name_ntuple_correction_mc).c_str());

    // Define the necessary histograms to calculate efficiency
    TH1F* hsig        = new TH1F("hsig"       ,"",ndim_corr,-TMath::Pi(),TMath::Pi());
    TH1F* hall        = new TH1F("hall"       ,"",ndim_corr,-TMath::Pi(),TMath::Pi());
    TH1F* hefficiency = new TH1F("hefficiency","",ndim_corr,-TMath::Pi(),TMath::Pi());
    hsig->Sumw2();
    hall->Sumw2();
    hefficiency->Sumw2();

    set_histogram_style(hsig, kViolet, std_line_width, std_marker_style, std_marker_size);
    set_histogram_style(hall, kCyan  , std_line_width, std_marker_style, std_marker_size);

    // Project into the histograms
    ntuple_mcreco->Project("hsig","h_phi",single_signal_cut);
    ntuple_mc->Project("hall","h_phi",pair_cut         );
    
    TCanvas* c = new TCanvas("c","",800,600);
    c->Draw();

    // efficiency PLOTS
    hefficiency->Divide(hsig,hall,1,1,"B");
    
    set_histogram_style(hefficiency, kViolet, std_line_width, std_marker_style, std_marker_size);
    
    hefficiency->Draw();
    hefficiency->GetYaxis()->SetRangeUser(0,1);
    hefficiency->SetTitle(";#phi;Single Hadron efficiency");

    c->Print("./plots/singlehadron_efficiency_phi.pdf");
}