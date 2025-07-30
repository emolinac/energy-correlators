#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

// std::string systematic = "hadron-correction-scheme";
std::string systematic = "probnnghost";

void macro_print_deviation_from_nominal()
{
    TFile* fnominal    = new TFile(("../output-files/"+namef_histos_corr_e2c).c_str());
    TFile* fsystematic = new TFile((systematic+"/output-files/"+namef_histos_corr_e2c).c_str());

    THStack* s = new THStack();
    TLegend* l = new TLegend();

    TH1F* h_nominal[nbin_jet_pt];
    TH1F* h_systematic[nbin_jet_pt];
    TH1F* h_deviations[nbin_jet_pt];

    for (int jet_pt_bin = 0 ; jet_pt_bin < nbin_jet_pt ; jet_pt_bin++)
    {
        h_nominal[jet_pt_bin]    = (TH1F*) fnominal->Get(Form("hcorr_e2c[%i]",jet_pt_bin));
        h_systematic[jet_pt_bin] = (TH1F*) fsystematic->Get(Form("hcorr_e2c[%i]",jet_pt_bin));
        h_deviations[jet_pt_bin] = new TH1F(Form("h_deviations%i",jet_pt_bin),"",nbin_rl_nominal,rl_binning);

        h_deviations[jet_pt_bin]->Add(h_nominal[jet_pt_bin],h_systematic[jet_pt_bin],1,-1);
        h_deviations[jet_pt_bin]->Divide(h_nominal[jet_pt_bin]);

        set_histogram_style(h_deviations[jet_pt_bin], corr_marker_color_jet_pt[jet_pt_bin], std_line_width, corr_marker_style_jet_pt[jet_pt_bin], std_marker_size);

        s->Add(h_deviations[jet_pt_bin],"E");
        l->AddEntry(h_deviations[jet_pt_bin],Form("%.1f<p_{T,jet}<%.1f (GeV)",jet_pt_binning[jet_pt_bin],jet_pt_binning[jet_pt_bin + 1]),"lf");
    }

    TCanvas* c = new TCanvas("c", "", 1920, 1080);
    c->Draw();

    s->Draw("NOSTACK");
    s->SetTitle(";R_{L};Frac. Dev. From Nominal");
    gPad->SetLogx(1);
    l->Draw("SAME");

    s->SetMaximum(0.49);
    s->SetMinimum(-0.49);

    c->Print(Form("./plots/fracdev_from_nominal_%s_relerrorcorr%.1f.pdf",systematic.c_str(),corr_rel_error));
}