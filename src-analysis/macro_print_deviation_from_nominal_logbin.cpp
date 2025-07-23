#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_deviation_from_nominal_logbin(bool normalize = true, bool do_print = false, int index_syst = 0)
{
        if (available_systematics[index_syst]=="ct") {
                std::cout<<"Remember CT dev from nominal is calculated inside macro_print_fullcorre2c_paircorr_2dunf_ct_niter.cpp !!!"<<std::endl; 
                return;
        }

        std::string systematic = available_systematics[index_syst];
        TFile* fout = new TFile((output_folder + devfromnom_namef[systematic]).c_str(),"RECREATE");
        gROOT->cd();

        TFile* fnominal    = new TFile((output_folder + namef_histos_paircorr_e2c_logbin).c_str());
        TFile* fsystematic = new TFile((output_folder+systematic_namef[systematic]).c_str());

        THStack* s     = new THStack();
        TLegend* l     = new TLegend();
        THStack* s_tau = new THStack();
        TLegend* l_tau = new TLegend();

        TH1F* h_nominal[Nbin_jet_pt];
        TH1F* h_systematic[Nbin_jet_pt];
        TH1F* h_deviations[Nbin_jet_pt];

        TH1F* h_nominal_tau[Nbin_jet_pt];
        TH1F* h_systematic_tau[Nbin_jet_pt];
        TH1F* h_deviations_tau[Nbin_jet_pt];

        for (int jet_pt_bin = 0 ; jet_pt_bin < Nbin_jet_pt ; jet_pt_bin++) {
                h_nominal[jet_pt_bin]    = (TH1F*) fnominal->Get(Form("hcorr_e2c%i",jet_pt_bin));
                h_systematic[jet_pt_bin] = (TH1F*) fsystematic->Get(Form("hcorr_e2c%i",jet_pt_bin));
                h_deviations[jet_pt_bin] = new TH1F(Form("h_deviations%i",jet_pt_bin),"",Nbin_rl_nominal,rl_nominal_binning);

                h_nominal_tau[jet_pt_bin]    = (TH1F*) fnominal->Get(Form("hcorr_tau%i",jet_pt_bin));
                h_systematic_tau[jet_pt_bin] = (TH1F*) fsystematic->Get(Form("hcorr_tau%i",jet_pt_bin));
                h_deviations_tau[jet_pt_bin] = new TH1F(Form("h_deviations_tau%i",jet_pt_bin),"",Nbin_rl_nominal,tau_nominal_binning);

                if (normalize) {
                        h_nominal[jet_pt_bin]->Scale(1./h_nominal[jet_pt_bin]->Integral());
                        h_systematic[jet_pt_bin]->Scale(1./h_systematic[jet_pt_bin]->Integral());

                        h_nominal_tau[jet_pt_bin]->Scale(1./h_nominal_tau[jet_pt_bin]->Integral());
                        h_systematic_tau[jet_pt_bin]->Scale(1./h_systematic_tau[jet_pt_bin]->Integral());
                }

                // h_deviations[jet_pt_bin]->Add(h_nominal[jet_pt_bin],h_systematic[jet_pt_bin],1,-1);
                // h_deviations[jet_pt_bin]->Divide(h_nominal[jet_pt_bin]);
                h_deviations[jet_pt_bin]->Divide(h_systematic[jet_pt_bin],h_nominal[jet_pt_bin],1,1);
                h_deviations_tau[jet_pt_bin]->Divide(h_systematic_tau[jet_pt_bin],h_nominal_tau[jet_pt_bin],1,1);

                set_histogram_style(h_deviations[jet_pt_bin], corr_marker_color_jet_pt[jet_pt_bin], std_line_width-1, corr_marker_style_jet_pt[jet_pt_bin], std_marker_size+1);
                set_histogram_style(h_deviations_tau[jet_pt_bin], corr_marker_color_jet_pt[jet_pt_bin], std_line_width-1, corr_marker_style_jet_pt[jet_pt_bin], std_marker_size+1);
                
                fout->cd();
                h_deviations[jet_pt_bin]->Write();
                h_deviations_tau[jet_pt_bin]->Write();
                gROOT->cd();
        }

        TCanvas* c = new TCanvas("c", "", 1920, 1080);
        c->Draw();

        for (int jet_pt_bin = 0 ; jet_pt_bin < Nbin_jet_pt ; jet_pt_bin++) {
                s->Add(h_deviations[jet_pt_bin],"E1 X0 L");
                l->AddEntry(h_deviations[jet_pt_bin],Form("%.1f<p^{jet}_{t}<%.1f GeV",jet_pt_binning[jet_pt_bin],jet_pt_binning[jet_pt_bin + 1]),"lf");

                s_tau->Add(h_deviations_tau[jet_pt_bin],"E1 X0 L");
                l_tau->AddEntry(h_deviations_tau[jet_pt_bin],Form("%.1f<p^{jet}_{t}<%.1f GeV",jet_pt_binning[jet_pt_bin],jet_pt_binning[jet_pt_bin + 1]),"lf");
        }

        s->Draw("NOSTACK");
        s->SetTitle(";R_{L};Frac. Dev. From Nominal");
        gPad->SetLogx(1);
        l->Draw("SAME");
        s->SetMaximum(1.5);
        s->SetMinimum(0.5);

        if (do_print) 
                c->Print(Form("./plots/dev_from_nominal_%s_norm-%s_logbin.pdf",systematic.c_str(),(normalize)?"yes":"no"));

        s_tau->Draw("NOSTACK");
        s_tau->SetTitle(";R_{L} #LT p^{jet}_{t} #GT(GeV);Frac. Dev. From Nominal");
        gPad->SetLogx(1);
        l_tau->Draw("SAME");
        s_tau->SetMaximum(1.5);
        s_tau->SetMinimum(0.5);

        if (do_print) 
                c->Print(Form("./plots/dev_from_nominal_tau_%s_norm-%s_logbin.pdf",systematic.c_str(),(normalize)?"yes":"no"));
}