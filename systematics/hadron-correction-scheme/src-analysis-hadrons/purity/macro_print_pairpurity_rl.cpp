#include "../../include/analysis-constants.h"
#include "../../include/analysis-cuts.h"
#include "../../include/directories.h"
#include "../../include/names.h"
#include "../../include/utils-algorithms.h"
#include "../../include/utils-visual.h"

void macro_print_pairpurity_rl()
{
    // Open the necessary files
    TFile* fpurity = new TFile((output_folder+namef_ntuple_e2c_pairpurity).c_str());

    // Get the corresponding Ntuples
    TNtuple* ntuple_purity = (TNtuple*) fpurity->Get((name_ntuple_purity).c_str());

    // Define the necessary histograms to calculate purity
    TH1F* hsig        = new TH1F("hsig"       ,"",Nbin_R_L,rl_binning);
    TH1F* hall        = new TH1F("hall"       ,"",Nbin_R_L,rl_binning);
    TH1F* hpurity = new TH1F("hpurity","",Nbin_R_L,rl_binning);
    hsig->Sumw2();
    hall->Sumw2();
    
    set_histogram_style(hsig, kViolet, std_line_width, std_marker_style, std_marker_size);
    set_histogram_style(hall, kCyan  , std_line_width, std_marker_style, std_marker_size);

    // Project into the histograms
    ntuple_purity->Project("hsig","R_L",pair_signal_cut);
    ntuple_purity->Project("hall","R_L",pair_cut);
    
    TLatex* tex = new TLatex();
    tex->SetTextColorAlpha(16,0.3);
    tex->SetTextSize(0.1991525);
    tex->SetTextAngle(26.15998);
    tex->SetLineWidth(2);
    
    TCanvas* c = new TCanvas("c","",800,600);
    c->Draw();

    gPad->SetLogx(1);

    // purity PLOTS
    hpurity->Divide(hsig,hall,1,1,"B");
    set_histogram_style(hpurity, kViolet, std_line_width, std_marker_style, std_marker_size);
    
    hpurity->Draw();
    hpurity->GetXaxis()->SetRangeUser(R_L_min,1);
    hpurity->GetYaxis()->SetRangeUser(0,1);
    hpurity->SetTitle(Form("#Delta R_{L}(truth-reco)<%.3f;R_{L};Pair purity",R_L_res));

    tex->DrawLatexNDC(0.3,0.3,"simulations");

    c->Print(Form("./plots/npair_purity_rl_deltarleq%.3f.pdf",R_L_res));
}