#include "../../include/analysis-constants.h"
#include "../../include/analysis-binning.h"
#include "../../include/analysis-cuts.h"
#include "../../include/directories.h"
#include "../../include/names.h"
#include "../../include/utils-algorithms.h"
#include "../../include/utils-visual.h"

void macro_print_pairefficiency_rl_jet_pt()
{
    // Open the necessary files
    TFile* fefficiency = new TFile((output_folder+namef_ntuple_e2c_pairefficiency).c_str());

    // Get the corresponding Ntuples
    TNtuple* ntuple_efficiency_reco = (TNtuple*) fefficiency->Get((name_ntuple_efficiency_reco).c_str());
    TNtuple* ntuple_efficiency_mc   = (TNtuple*) fefficiency->Get((name_ntuple_efficiency_mc).c_str());

    // Define the necessary histograms to calculate efficiency
    TH1F* hsig[Nbin_jet_pt];
    TH1F* hall[Nbin_jet_pt];
    TH1F* hefficiency[Nbin_jet_pt];

    for(int jet_pt_bin = 0 ; jet_pt_bin < Nbin_jet_pt ; jet_pt_bin++)
    {
        hsig[jet_pt_bin]        = new TH1F(Form("hsig[%i]",jet_pt_bin)       ,"",Nbin_R_L,rl_binning);
        hall[jet_pt_bin]        = new TH1F(Form("hall[%i]",jet_pt_bin)       ,"",Nbin_R_L,rl_binning);
        hefficiency[jet_pt_bin] = new TH1F(Form("hefficiency[%i]",jet_pt_bin),"",Nbin_R_L,rl_binning);

        hsig[jet_pt_bin]->Sumw2();
        hall[jet_pt_bin]->Sumw2();
        hefficiency[jet_pt_bin]->Sumw2();

        set_histogram_style(hsig[jet_pt_bin], corr_marker_color_jet_pt[jet_pt_bin], std_line_width, corr_marker_style_jet_pt[jet_pt_bin], std_marker_size);
        set_histogram_style(hall[jet_pt_bin], std_marker_color_jet_pt[jet_pt_bin], std_line_width, std_marker_style_jet_pt[jet_pt_bin], std_marker_size);
    }

    // Project into the histograms
    for(int jet_pt_bin = 0 ; jet_pt_bin < Nbin_jet_pt ; jet_pt_bin++)
    {
        ntuple_efficiency_reco->Project(Form("hsig[%i]",jet_pt_bin),"R_L",pair_jetpt_signal_cut[jet_pt_bin]);
        ntuple_efficiency_mc->Project(Form("hall[%i]",jet_pt_bin)  ,"R_L",pair_jetpt_cut[jet_pt_bin]);
    }
    
    TCanvas* c = new TCanvas("c","",800,600);
    c->Draw();

    TLatex* tex = new TLatex();
    tex->SetTextColorAlpha(16,0.3);
    tex->SetTextSize(0.1991525);
    tex->SetTextAngle(26.15998);
    tex->SetLineWidth(2);
    
    gPad->SetLogx(1);

    // efficiency PLOTS
    THStack* s_efficiency = new THStack();
    TLegend* l_efficiency = new TLegend();

    for(int jet_pt_bin = 0 ; jet_pt_bin < Nbin_jet_pt ; jet_pt_bin++)
    {
        hefficiency[jet_pt_bin]->Divide(hsig[jet_pt_bin],hall[jet_pt_bin],1,1,"B");
        set_histogram_style(hefficiency[jet_pt_bin], corr_marker_color_jet_pt[jet_pt_bin], std_line_width, std_marker_style_jet_pt[jet_pt_bin], std_marker_size);

        s_efficiency->Add(hefficiency[jet_pt_bin]);
        l_efficiency->AddEntry(hefficiency[jet_pt_bin],Form("%.1f<p^{jet}_{t}(GeV)<%.1f",jet_pt_binning[jet_pt_bin],jet_pt_binning[jet_pt_bin+1]),"lpf");
    }
    
    
    s_efficiency->Draw("NOSTACK");
    s_efficiency->GetXaxis()->SetRangeUser(R_L_min,1);
    s_efficiency->SetMaximum(0.6);
    s_efficiency->SetTitle(Form("#Delta R_{L}(truth-reco)<%.3f;R_{L};Pair efficiency",R_L_res));
    l_efficiency->Draw("SAME");

    tex->DrawLatexNDC(0.3,0.3,"simulations");

    c->Print(Form("./plots/npair_efficiency_rl_jetpt_deltarleq%.3f.pdf",R_L_res));
}