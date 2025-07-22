#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_deviation_from_nominal_logbin_shapect(bool normalize = true, bool do_print = true, int index_syst = 5)
{
        if (available_systematics[index_syst]!="shape-ct") {
                std::cout<<"This code is only valid for SHAPE-CT !!!"<<std::endl; 
                return;
        }

        std::string systematic = available_systematics[index_syst];
        
        TFile* fout_dev = new TFile((output_folder+devfromnom_namef[systematic]).c_str());
        
        THStack* s     = new THStack();
        THStack* s_tau = new THStack();
        TLegend* l     = new TLegend(0.2,0.7,0.3,0.9);
        TLegend* l_tau = new TLegend(0.2,0.7,0.3,0.9);
        
        TH1F* h_nominal[Nbin_jet_pt];
        TH1F* h_systematic[Nbin_jet_pt];
        TH1F* h_deviations[Nbin_jet_pt];

        TH1F* h_nominal_tau[Nbin_jet_pt];
        TH1F* h_systematic_tau[Nbin_jet_pt];
        TH1F* h_deviations_tau[Nbin_jet_pt];

        for (int jet_pt_bin = 0 ; jet_pt_bin < Nbin_jet_pt ; jet_pt_bin++) {
                h_deviations[jet_pt_bin]     = (TH1F*) fout_dev->Get(Form("h_deviations%i",jet_pt_bin));
                h_deviations_tau[jet_pt_bin] = (TH1F*) fout_dev->Get(Form("h_deviations_tau%i",jet_pt_bin));

                set_histogram_style(h_deviations[jet_pt_bin], corr_marker_color_jet_pt[jet_pt_bin], std_line_width-1, corr_marker_style_jet_pt[jet_pt_bin], std_marker_size+1);
                set_histogram_style(h_deviations_tau[jet_pt_bin], corr_marker_color_jet_pt[jet_pt_bin], std_line_width-1, corr_marker_style_jet_pt[jet_pt_bin], std_marker_size+1);
        }

        TCanvas* c = new TCanvas("c", "", 1920, 1080);
        c->Draw();

        for (int jet_pt_bin = 0 ; jet_pt_bin < Nbin_jet_pt ; jet_pt_bin++) {
                s->Add(h_deviations[jet_pt_bin],"E1 X0");
                l->AddEntry(h_deviations[jet_pt_bin],Form("%.1f<p^{jet}_{t}<%.1f GeV",jet_pt_binning[jet_pt_bin],jet_pt_binning[jet_pt_bin+1]),"lf");

                s_tau->Add(h_deviations_tau[jet_pt_bin],"E1 X0");
                l_tau->AddEntry(h_deviations_tau[jet_pt_bin],Form("%.1f<p^{jet}_{t}<%.1f GeV",jet_pt_binning[jet_pt_bin],jet_pt_binning[jet_pt_bin+1]),"lf");
        }

        s->Draw("NOSTACK");
        s->SetTitle(";R_{L};Corr. Pseudodata / Truth");
        gPad->SetLogx(1);
        l->Draw("SAME");
        s->SetMaximum(1.5);
        s->SetMinimum(0.5);

        if (do_print) 
                c->Print(Form("./plots/dev_from_nominal_%s_norm-%s_logbin.pdf",systematic.c_str(),(normalize)?"yes":"no"));

        s_tau->Draw("NOSTACK");
        s_tau->SetTitle(";R_{L} #LT p^{jet}_{t} #GT(GeV);Corr. Pseudodata / Truth");
        gPad->SetLogx(1);
        l_tau->Draw("SAME");
        s_tau->SetMaximum(1.5);
        s_tau->SetMinimum(0.5);

        if (do_print) 
                c->Print(Form("./plots/dev_from_nominal_tau_%s_norm-%s_logbin.pdf",systematic.c_str(),(normalize)?"yes":"no"));
}