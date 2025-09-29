#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_deviation_from_nominal_ct(bool normalize = true, bool do_print = true, int index_syst = 1)
{
        if (available_systematics[index_syst]!="ct") {
                std::cout<<"This code is only valid for CT !!!"<<std::endl; 
                return;
        }

        std::string systematic = available_systematics[index_syst];
        
        TFile* fout_dev = new TFile((output_folder + devfromnom_namef[systematic]).c_str());
        
        TH1F* h_nominal_eec[nbin_jet_pt];
        TH1F* h_systematic_eec[nbin_jet_pt];
        TH1F* h_deviations_eec[nbin_jet_pt];

        TH1F* h_nominal_tau[nbin_jet_pt];
        TH1F* h_systematic_tau[nbin_jet_pt];
        TH1F* h_deviations_tau[nbin_jet_pt];

        TH1F* h_nominal_eec_eqcharge[nbin_jet_pt];
        TH1F* h_systematic_eec_eqcharge[nbin_jet_pt];
        TH1F* h_deviations_eec_eqcharge[nbin_jet_pt];

        TH1F* h_nominal_eec_neqcharge[nbin_jet_pt];
        TH1F* h_systematic_eec_neqcharge[nbin_jet_pt];
        TH1F* h_deviations_eec_neqcharge[nbin_jet_pt];

        for (int jet_pt_bin = 0 ; jet_pt_bin < nbin_jet_pt ; jet_pt_bin++) {
                h_deviations_eec[jet_pt_bin]           = (TH1F*) fout_dev->Get(Form("h_deviations_eec%i",jet_pt_bin));
                h_deviations_tau[jet_pt_bin]           = (TH1F*) fout_dev->Get(Form("h_deviations_tau%i",jet_pt_bin));
                h_deviations_eec_eqcharge[jet_pt_bin]  = (TH1F*) fout_dev->Get(Form("h_deviations_eec_eqcharge%i",jet_pt_bin));
                h_deviations_eec_neqcharge[jet_pt_bin] = (TH1F*) fout_dev->Get(Form("h_deviations_eec_neqcharge%i",jet_pt_bin));

                set_histogram_style(h_deviations_eec[jet_pt_bin], corr_marker_color_jet_pt[jet_pt_bin], std_line_width-1, corr_marker_style_jet_pt[jet_pt_bin], std_marker_size+1);
                set_histogram_style(h_deviations_tau[jet_pt_bin], corr_marker_color_jet_pt[jet_pt_bin], std_line_width-1, corr_marker_style_jet_pt[jet_pt_bin], std_marker_size+1);
                set_histogram_style(h_deviations_eec_eqcharge[jet_pt_bin], corr_marker_color_jet_pt[jet_pt_bin], std_line_width-1, corr_marker_style_jet_pt[jet_pt_bin], std_marker_size+1);
                set_histogram_style(h_deviations_eec_neqcharge[jet_pt_bin], corr_marker_color_jet_pt[jet_pt_bin], std_line_width-1, corr_marker_style_jet_pt[jet_pt_bin], std_marker_size+1);
        }

        TCanvas* c = new TCanvas("c", "", 1920, 1080);
        c->Draw();

        THStack* s           = new THStack();
        THStack* s_tau       = new THStack();
        THStack* s_eqcharge  = new THStack();
        THStack* s_neqcharge = new THStack();
        TLegend* l           = new TLegend(0.2,0.7,0.3,0.9);
        TLegend* l_tau       = new TLegend(0.2,0.7,0.3,0.9);
        TLegend* l_eqcharge  = new TLegend(0.2,0.7,0.3,0.9);
        TLegend* l_neqcharge = new TLegend(0.2,0.7,0.3,0.9);

        for (int jet_pt_bin = 0 ; jet_pt_bin < nbin_jet_pt ; jet_pt_bin++) {
                s->Add(h_deviations_eec[jet_pt_bin],"HIST");
                l->AddEntry(h_deviations_eec[jet_pt_bin],Form("%.1f<p_{T,jet}<%.1f (GeV)",jet_pt_binning[jet_pt_bin],jet_pt_binning[jet_pt_bin + 1]),"lf");

                s_tau->Add(h_deviations_tau[jet_pt_bin],"HIST");
                l_tau->AddEntry(h_deviations_tau[jet_pt_bin],Form("%.1f<p_{T,jet}<%.1f (GeV)",jet_pt_binning[jet_pt_bin],jet_pt_binning[jet_pt_bin + 1]),"lf");
                
                s_eqcharge->Add(h_deviations_eec_eqcharge[jet_pt_bin],"HIST");
                l_eqcharge->AddEntry(h_deviations_eec_eqcharge[jet_pt_bin],Form("%.1f<p_{T,jet}<%.1f (GeV)",jet_pt_binning[jet_pt_bin],jet_pt_binning[jet_pt_bin + 1]),"lf");
                
                s_neqcharge->Add(h_deviations_eec_neqcharge[jet_pt_bin],"HIST");
                l_neqcharge->AddEntry(h_deviations_eec_neqcharge[jet_pt_bin],Form("%.1f<p_{T,jet}<%.1f (GeV)",jet_pt_binning[jet_pt_bin],jet_pt_binning[jet_pt_bin + 1]),"lf");
        }

        s->Draw("NOSTACK");
        s->SetTitle(";R_{L};#LT Corr. Pseudodata / Truth #GT");
        gPad->SetLogx(1);
        l->Draw("SAME");
        s->SetMaximum(1.5);
        s->SetMinimum(0.5);

        if (do_print) 
                c->Print(Form("./plots/dev_from_nominal_%s_norm-%s.pdf",systematic.c_str(),(normalize)?"yes":"no"));

        s_tau->Draw("NOSTACK");
        s_tau->SetTitle(";R_{L} #LT p_{T,jet} #GT(GeV);#LT Corr. Pseudodata / Truth #GT");
        gPad->SetLogx(1);
        l_tau->Draw("SAME");
        s_tau->SetMaximum(1.5);
        s_tau->SetMinimum(0.5);

        if (do_print) 
                c->Print(Form("./plots/dev_from_nominal_tau_%s_norm-%s.pdf",systematic.c_str(),(normalize)?"yes":"no"));

        s_eqcharge->Draw("NOSTACK");
        s_eqcharge->SetTitle(";R_{L};#LT Corr. Pseudodata / Truth #GT");
        gPad->SetLogx(1);
        l_eqcharge->Draw("SAME");
        s_eqcharge->SetMaximum(1.5);
        s_eqcharge->SetMinimum(0.5);

        if (do_print) 
                c->Print(Form("./plots/dev_from_nominal_eec_eqcharge_%s_norm-%s.pdf",systematic.c_str(),(normalize)?"yes":"no"));

        s_neqcharge->Draw("NOSTACK");
        s_neqcharge->SetTitle(";R_{L};#LT Corr. Pseudodata / Truth #GT");
        gPad->SetLogx(1);
        l_neqcharge->Draw("SAME");
        s_neqcharge->SetMaximum(1.5);
        s_neqcharge->SetMinimum(0.5);

        if (do_print) 
                c->Print(Form("./plots/dev_from_nominal_eec_neqcharge_%s_norm-%s.pdf",systematic.c_str(),(normalize)?"yes":"no"));
}