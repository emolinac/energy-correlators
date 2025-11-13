#include "../../include/analysis-constants.h"
#include "../../include/analysis-binning.h"
#include "../../include/analysis-cuts.h"
#include "../../include/directories.h"
#include "../../include/names.h"
#include "../../include/utils.h"
#include "../../include/utils-visual.h"

void macro_print_pairefficiency_rl_jet_pt()
{
    // Open the necessary files
    TFile* fdata       = new TFile((output_folder + namef_ntuple_eec).c_str());
    TFile* fefficiency = new TFile((output_folder + namef_ntuple_eec_hadroncorrections).c_str());

    // Get the corresponding Ntuples
    TNtuple* ntuple_data            = (TNtuple*) fdata->Get((name_ntuple_data).c_str());
    TNtuple* ntuple_efficiency_reco = (TNtuple*) fefficiency->Get((name_ntuple_correction_reco).c_str());
    TNtuple* ntuple_efficiency_mc   = (TNtuple*) fefficiency->Get((name_ntuple_correction_mc).c_str());

    // Determine log binnning
    double binning[nbin_rl_nominal+1];
    determine_log10binning(nbin_rl_nominal, rl_min, rl_max, binning);

    // Define the necessary histograms to calculate efficiency
    TH1F* hsig[nbin_jet_pt];
    TH1F* hall[nbin_jet_pt];
    TH1F* hefficiency[nbin_jet_pt];

    for (int jet_pt_bin = 0 ; jet_pt_bin < nbin_jet_pt ; jet_pt_bin++)
    {
        hsig[jet_pt_bin]    = new TH1F(Form("hsig[%i]",jet_pt_bin)   ,"",nbin_rl_nominal,rl_min, rl_max);
        hall[jet_pt_bin]    = new TH1F(Form("hall[%i]",jet_pt_bin)   ,"",nbin_rl_nominal,rl_min, rl_max);
        hefficiency[jet_pt_bin] = new TH1F(Form("hefficiency[%i]",jet_pt_bin),"",nbin_rl_nominal,rl_min, rl_max);

        hsig[jet_pt_bin]->Sumw2();
        hall[jet_pt_bin]->Sumw2();
        hefficiency[jet_pt_bin]->Sumw2();

        set_histogram_style(hsig[jet_pt_bin], corr_marker_color_jet_pt[jet_pt_bin], std_line_width, corr_marker_style_jet_pt[jet_pt_bin], std_marker_size);
        set_histogram_style(hall[jet_pt_bin], std_marker_color_jet_pt[jet_pt_bin], std_line_width, std_marker_style_jet_pt[jet_pt_bin], std_marker_size);
    }

    // Define the necessary histograms to show data and corrected data
    TH1F* hcorr_data[nbin_jet_pt];
    TH1F* hall_data[nbin_jet_pt];

    for (int jet_pt_bin = 0 ; jet_pt_bin < nbin_jet_pt ; jet_pt_bin++)
    {
        hcorr_data[jet_pt_bin] = new TH1F(Form("hcorr_data[%i]",jet_pt_bin),"",nbin_rl_nominal,rl_min, rl_max);
        hall_data[jet_pt_bin] = new TH1F(Form("hall_data[%i]",jet_pt_bin),"",nbin_rl_nominal,rl_min, rl_max);

        hcorr_data[jet_pt_bin]->Sumw2();
        hall_data[jet_pt_bin]->Sumw2();

        set_histogram_style(hcorr_data[jet_pt_bin], corr_marker_color_jet_pt[jet_pt_bin], std_line_width, corr_marker_style_jet_pt[jet_pt_bin], std_marker_size);
        set_histogram_style(hall_data[jet_pt_bin], std_marker_color_jet_pt[jet_pt_bin] , std_line_width, std_marker_style_jet_pt[jet_pt_bin], std_marker_size);
    }

    // Project into the histograms
    for (int jet_pt_bin = 0 ; jet_pt_bin < nbin_jet_pt ; jet_pt_bin++)
    {
        ntuple_efficiency_reco->Project(Form("hsig[%i]",jet_pt_bin),"R_L",pair_jet_pt_signal_cut[jet_pt_bin]+"(h1_pt<7&&h2_pt<7)");
        ntuple_efficiency_mc->Project(Form("hall[%i]",jet_pt_bin),"R_L"  ,pair_jet_pt_cut[jet_pt_bin]+"(h1_pt<7&&h2_pt<7)");
        ntuple_data->Project(Form("hcorr_data[%i]",jet_pt_bin),"R_L"     ,pair_jet_pt_cut[jet_pt_bin]+"(h1_pt<7&&h2_pt<7)");
        ntuple_data->Project(Form("hall_data[%i]",jet_pt_bin),"R_L"      ,pair_jet_pt_cut[jet_pt_bin]+"(h1_pt<7&&h2_pt<7)");
    }
    
    TCanvas* c = new TCanvas("c","",800,600);
    c->Draw();

    TLatex* tex = new TLatex();
    set_lhcb_watermark_properties(tex);
    
    // MCRECO PLOTS
    gPad->SetLogx(1);
    gPad->SetLogy(1);

    THStack* s = new THStack();
    TLegend* l = new TLegend();

    for (int jet_pt_bin = 0 ; jet_pt_bin < nbin_jet_pt ; jet_pt_bin++)
    {
        s->Add(hsig[jet_pt_bin]);
        s->Add(hall[jet_pt_bin]);
        l->AddEntry(hsig[jet_pt_bin],Form("Reco: %.1f<p_{T,jet}(GeV)<%.1f",jet_pt_binning[jet_pt_bin],jet_pt_binning[jet_pt_bin + 1]),"lpf");
        l->AddEntry(hall[jet_pt_bin],Form("MC  : %.1f<p_{T,jet}(GeV)<%.1f",jet_pt_binning[jet_pt_bin],jet_pt_binning[jet_pt_bin + 1]),"lpf");
    }

    s->Draw("NOSTACK");
    s->GetXaxis()->SetRangeUser(rl_min,1);
    s->SetTitle(Form("#Delta R_{L}(truth-reco)<%.3f;R_{L};N_{pair}",rl_resolution));
    l->Draw("SAME");

    tex->DrawLatexNDC(0.3,0.3,"simulations");

    c->Print(Form("./plots/npair_rl_recovsmc_jet_pt_deltarleq%.3f.pdf",rl_resolution));
    gPad->SetLogy(0);

    // efficiency PLOTS
    THStack* s_efficiency = new THStack();
    TLegend* l_efficiency = new TLegend();

    for (int jet_pt_bin = 0 ; jet_pt_bin < nbin_jet_pt ; jet_pt_bin++)
    {
        hefficiency[jet_pt_bin]->Divide(hsig[jet_pt_bin],hall[jet_pt_bin],1,1,"B");
        set_histogram_style(hefficiency[jet_pt_bin], corr_marker_color_jet_pt[jet_pt_bin], std_line_width, std_marker_style_jet_pt[jet_pt_bin], std_marker_size);

        s_efficiency->Add(hefficiency[jet_pt_bin]);
        l_efficiency->AddEntry(hefficiency[jet_pt_bin],Form("%.1f<p_{T,jet}(GeV)<%.1f",jet_pt_binning[jet_pt_bin],jet_pt_binning[jet_pt_bin + 1]),"lpf");
    }
    
    
    s_efficiency->Draw("NOSTACK");
    s_efficiency->GetXaxis()->SetRangeUser(rl_min,1);
    s_efficiency->SetMaximum(0.6);
    s_efficiency->SetTitle(Form("#Delta R_{L}(truth-reco)<%.3f;R_{L};Pair efficiency",rl_resolution));
    l_efficiency->Draw("SAME");

    tex->DrawLatexNDC(0.3,0.3,"simulations");

    c->Print(Form("./plots/npair_efficiency_rl_jet_pt_deltarleq%.3f.pdf",rl_resolution));
    
    // DATA PLOTS
    THStack* s_data = new THStack();
    TLegend* l_data = new TLegend();

    for (int jet_pt_bin = 0 ; jet_pt_bin < nbin_jet_pt ; jet_pt_bin++)
    {
        hcorr_data[jet_pt_bin]->Divide(hefficiency[jet_pt_bin]);

        s_data->Add(hcorr_data[jet_pt_bin]);
        s_data->Add(hall_data[jet_pt_bin]);

        l_data->AddEntry(hcorr_data[jet_pt_bin],Form("Data w/ efficiency : %.1f<p_{T,jet}(GeV)<%.1f",jet_pt_binning[jet_pt_bin],jet_pt_binning[jet_pt_bin + 1]),"lpf");
        l_data->AddEntry(hall_data[jet_pt_bin],Form("Data All : %.1f<p_{T,jet}(GeV)<%.1f",jet_pt_binning[jet_pt_bin],jet_pt_binning[jet_pt_bin + 1])      ,"lpf");
    }
    
    s_data->Draw("NOSTACK");
    s_data->GetXaxis()->SetRangeUser(rl_min,1);
    s_data->SetTitle(Form("#Delta R_{L}(truth-reco)<%.3f;R_{L};N_{pair}",rl_resolution));
    l_data->Draw("SAME");

    gPad->SetLogx(1);
    gPad->SetLogy(1);

    tex->DrawLatexNDC(0.3,0.3,"simulations");

    c->Print(Form("./plots/npair_wefficiency_rl_data_jet_pt_deltarleq%.3f.pdf",rl_resolution));
}