#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

// std::string systematic = "syst-corr-binning";
// std::string systematic = "syst-corr-paradigm";
// std::string systematic = "syst-probnnghost";
// std::string systematic = "syst-min-track-pt";
std::string systematic = "syst-unfolding-dim";

void macro_print_deviation_from_nominal_logbin(bool normalize = false)
{
    std::string syst_name = (systematic=="syst-corr-paradigm") ? namef_histos_corr_e2c_logbin : namef_histos_paircorr_e2c_logbin;
    TFile* fnominal    = new TFile(("../output-files/"+namef_histos_paircorr_e2c_logbin).c_str());
    TFile* fsystematic = new TFile((systematic+"/output-files/"+syst_name).c_str());

    THStack* s = new THStack();
    TLegend* l = new TLegend();

    TH1F* h_nominal[Nbin_jet_pt];
    TH1F* h_systematic[Nbin_jet_pt];
    TH1F* h_deviations[Nbin_jet_pt];

    for(int jet_pt_bin = 0 ; jet_pt_bin < Nbin_jet_pt ; jet_pt_bin++)
    {
        h_nominal[jet_pt_bin]    = (TH1F*) fnominal->Get(Form("hcorr_e2c%i",jet_pt_bin));
        h_systematic[jet_pt_bin] = (TH1F*) fsystematic->Get(Form("hcorr_e2c%i",jet_pt_bin));
        h_deviations[jet_pt_bin] = new TH1F(Form("h_deviations%i",jet_pt_bin),"",Nbin_R_L_logbin,rl_logbinning);

        if(normalize)
        {
            h_nominal[jet_pt_bin]->Scale(1./h_nominal[jet_pt_bin]->Integral());
            h_systematic[jet_pt_bin]->Scale(1./h_systematic[jet_pt_bin]->Integral());
        }

        h_deviations[jet_pt_bin]->Add(h_nominal[jet_pt_bin],h_systematic[jet_pt_bin],1,-1);
        h_deviations[jet_pt_bin]->Divide(h_nominal[jet_pt_bin]);

        set_histogram_style(h_deviations[jet_pt_bin], corr_marker_color_jet_pt[jet_pt_bin], std_line_width-1, corr_marker_style_jet_pt[jet_pt_bin], std_marker_size+1);

        s->Add(h_deviations[jet_pt_bin],"E1 X0 L");
        l->AddEntry(h_deviations[jet_pt_bin],Form("%.1f<p^{jet}_{t}<%.1f GeV",jet_pt_binning[jet_pt_bin],jet_pt_binning[jet_pt_bin+1]),"lf");
    }

    TCanvas* c = new TCanvas("c","",1920,1080);
    c->Draw();

    s->Draw("NOSTACK");
    s->SetTitle(";R_{L};Frac. Dev. From Nominal");
    gPad->SetLogx(1);
    l->Draw("SAME");

    s->SetMaximum(0.9);
    s->SetMinimum(-0.9);

    c->Print(Form("./plots/dev_from_nominal_%s_norm-%s_logbin.pdf",systematic.c_str(),(normalize)?"yes":"no"));
}