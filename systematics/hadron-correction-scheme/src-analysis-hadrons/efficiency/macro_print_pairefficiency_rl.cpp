#include "../../include/analysis-constants.h"
#include "../../include/analysis-cuts.h"
#include "../../include/directories.h"
#include "../../include/names.h"
#include "../../include/utils-algorithms.h"
#include "../../include/utils-visual.h"

void macro_print_pairefficiency_rl()
{
    // Open the necessary files
    TFile* fefficiency = new TFile((output_folder+namef_ntuple_e2c_pairefficiency).c_str());

    // Get the corresponding Ntuples
    TNtuple* ntuple_efficiency_reco = (TNtuple*) fefficiency->Get((name_ntuple_efficiency_reco).c_str());
    TNtuple* ntuple_efficiency_mc   = (TNtuple*) fefficiency->Get((name_ntuple_efficiency_mc).c_str());

    // Define the necessary histograms to calculate efficiency
    TH1F* hsig        = new TH1F("hsig"       ,"",Nbin_R_L,rl_binning);
    TH1F* hall        = new TH1F("hall"       ,"",Nbin_R_L,rl_binning);
    TH1F* hefficiency = new TH1F("hefficiency","",Nbin_R_L,rl_binning);
    hsig->Sumw2();
    hall->Sumw2();
    
    set_histogram_style(hsig, kViolet, std_line_width, std_marker_style, std_marker_size);
    set_histogram_style(hall, kCyan  , std_line_width, std_marker_style, std_marker_size);

    // Project into the histograms
    ntuple_efficiency_reco->Project("hsig","R_L_truth",pair_signal_cut);
    ntuple_efficiency_mc->Project("hall","R_L",pair_cut);
    
    TLatex* tex = new TLatex();
    tex->SetTextColorAlpha(16,0.3);
    tex->SetTextSize(0.1991525);
    tex->SetTextAngle(26.15998);
    tex->SetLineWidth(2);
    
    TCanvas* c = new TCanvas("c","",800,600);
    c->Draw();

    gPad->SetLogx(1);

    // efficiency PLOTS
    hefficiency->Divide(hsig,hall,1,1,"B");
    set_histogram_style(hefficiency, kViolet, std_line_width, std_marker_style, std_marker_size);
    
    hefficiency->Draw();
    hefficiency->GetXaxis()->SetRangeUser(R_L_min,1);
    hefficiency->GetYaxis()->SetRangeUser(0,0.6);
    hefficiency->SetTitle(Form("#Delta R_{L}(truth-reco)<%.3f;R_{L};Pair efficiency",R_L_res));

    tex->DrawLatexNDC(0.3,0.3,"simulations");

    c->Print(Form("./plots/npair_efficiency_rl_deltarleq%.3f.pdf",R_L_res));
}