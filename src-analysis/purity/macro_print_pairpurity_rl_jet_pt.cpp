#include "../../include/analysis-constants.h"
#include "../../include/analysis-cuts.h"
#include "../../include/directories.h"
#include "../../include/names.h"
#include "../../include/utils-algorithms.h"
#include "../../include/utils-visual.h"

void macro_print_pairpurity_rl_jet_pt()
{
    // Open the necessary files
    TFile* fdata   = new TFile((output_folder+namef_ntuple_e2c).c_str());
    TFile* fpurity = new TFile((output_folder+namef_ntuple_e2c_purity).c_str());

    // Get the corresponding Ntuples
    TNtuple* ntuple_data   = (TNtuple*) fdata->Get((name_ntuple_data).c_str());
    TNtuple* ntuple_purity = (TNtuple*) fpurity->Get((name_ntuple_purity).c_str());

    // Determine log binnning
    double binning[Nbin_R_L+1];
    determine_log10binning(Nbin_R_L, R_L_min, R_L_max, binning);

    // Define the necessary histograms to calculate purity
    TH1F* hsig[Nbin_jet_pt];
    TH1F* hall[Nbin_jet_pt];
    TH1F* hpurity[Nbin_jet_pt];

    for(int jet_pt_bin = 0 ; jet_pt_bin < Nbin_jet_pt ; jet_pt_bin++)
    {
        hsig[jet_pt_bin]    = new TH1F(Form("hsig[%i]",jet_pt_bin)   ,"",Nbin_R_L,R_L_min,R_L_max);
        hall[jet_pt_bin]    = new TH1F(Form("hall[%i]",jet_pt_bin)   ,"",Nbin_R_L,R_L_min,R_L_max);
        hpurity[jet_pt_bin] = new TH1F(Form("hpurity[%i]",jet_pt_bin),"",Nbin_R_L,R_L_min,R_L_max);

        hsig[jet_pt_bin]->Sumw2();
        hall[jet_pt_bin]->Sumw2();
        hpurity[jet_pt_bin]->Sumw2();

        set_histogram_style(hsig[jet_pt_bin], corr_marker_color_jet_pt[jet_pt_bin], std_line_width, corr_marker_style_jet_pt[jet_pt_bin], std_marker_size);
        set_histogram_style(hall[jet_pt_bin], std_marker_color_jet_pt[jet_pt_bin], std_line_width, std_marker_style_jet_pt[jet_pt_bin], std_marker_size);
    }

    // Define the necessary histograms to show data and corrected data
    TH1F* hsig_data[Nbin_jet_pt];
    TH1F* hall_data[Nbin_jet_pt];

    for(int jet_pt_bin = 0 ; jet_pt_bin < Nbin_jet_pt ; jet_pt_bin++)
    {
        hsig_data[jet_pt_bin] = new TH1F(Form("hsig_data[%i]",jet_pt_bin),"",Nbin_R_L,R_L_min,R_L_max);
        hall_data[jet_pt_bin] = new TH1F(Form("hall_data[%i]",jet_pt_bin),"",Nbin_R_L,R_L_min,R_L_max);

        hsig_data[jet_pt_bin]->Sumw2();
        hall_data[jet_pt_bin]->Sumw2();

        set_histogram_style(hsig_data[jet_pt_bin], corr_marker_color_jet_pt[jet_pt_bin], std_line_width, corr_marker_style_jet_pt[jet_pt_bin], std_marker_size);
        set_histogram_style(hall_data[jet_pt_bin], std_marker_color_jet_pt[jet_pt_bin] , std_line_width, std_marker_style_jet_pt[jet_pt_bin], std_marker_size);
    }

    // Project into the histograms
    for(int jet_pt_bin = 0 ; jet_pt_bin < Nbin_jet_pt ; jet_pt_bin++)
    {
        ntuple_purity->Project(Form("hsig[%i]",jet_pt_bin),"R_L",pair_jetpt_signal_cut[jet_pt_bin]);
        ntuple_purity->Project(Form("hall[%i]",jet_pt_bin),"R_L",pair_jetpt_cut[jet_pt_bin]);
        ntuple_data->Project(Form("hsig_data[%i]",jet_pt_bin),"R_L",pair_jetpt_cut[jet_pt_bin]);
        ntuple_data->Project(Form("hall_data[%i]",jet_pt_bin),"R_L",pair_jetpt_cut[jet_pt_bin]);
    }
    
    TCanvas* c = new TCanvas("c","",800,600);
    c->Draw();

    TLatex* tex = new TLatex();
    tex->SetTextColorAlpha(16,0.3);
    tex->SetTextSize(0.1991525);
    tex->SetTextAngle(26.15998);
    tex->SetLineWidth(2);

    // MCRECO PLOTS
    gPad->SetLogx(1);
    gPad->SetLogy(1);

    THStack* s = new THStack();
    TLegend* l = new TLegend();

    for(int jet_pt_bin = 0 ; jet_pt_bin < Nbin_jet_pt ; jet_pt_bin++)
    {
        s->Add(hsig[jet_pt_bin]);
        s->Add(hall[jet_pt_bin]);
        l->AddEntry(hsig[jet_pt_bin],Form("MCReco Signal : %.1f<Jet P_{T}(GeV)<%.1f",jet_pt_binning[jet_pt_bin],jet_pt_binning[jet_pt_bin+1]),"lpf");
        l->AddEntry(hall[jet_pt_bin],Form("MCReco All : %.1f<Jet P_{T}(GeV)<%.1f",jet_pt_binning[jet_pt_bin],jet_pt_binning[jet_pt_bin+1])   ,"lpf");
    }

    s->Draw("NOSTACK");
    s->GetXaxis()->SetRangeUser(R_L_min,1);
    s->SetTitle(Form("#Delta R_{L}(truth-reco)<%.3f;R_{L};N_{pair}",R_L_res));
    l->Draw("SAME");

    tex->DrawLatexNDC(0.3,0.3,"simulations");

    c->Print(Form("../../plots/purity/npair_rl_signalvsall_jetpt_deltarleq%.3f.pdf",R_L_res));
    
    gPad->SetLogy(0);

    // PURITY PLOTS
    THStack* s_purity = new THStack();
    TLegend* l_purity = new TLegend();

    for(int jet_pt_bin = 0 ; jet_pt_bin < Nbin_jet_pt ; jet_pt_bin++)
    {
        hpurity[jet_pt_bin]->Divide(hsig[jet_pt_bin],hall[jet_pt_bin],1,1,"B");
        set_histogram_style(hpurity[jet_pt_bin], corr_marker_color_jet_pt[jet_pt_bin], std_line_width, std_marker_style_jet_pt[jet_pt_bin], std_marker_size);

        s_purity->Add(hpurity[jet_pt_bin]);
        l_purity->AddEntry(hpurity[jet_pt_bin],Form("%.1f<Jet P_{T}(GeV)<%.1f",jet_pt_binning[jet_pt_bin],jet_pt_binning[jet_pt_bin+1]),"lpf");
    }
    
    s_purity->Draw("NOSTACK");
    s_purity->GetXaxis()->SetRangeUser(R_L_min,1);
    s_purity->SetMaximum(1.);
    s_purity->SetTitle(Form("#Delta R_{L}(truth-reco)<%.3f;R_{L};Pair Purity",R_L_res));
    l_purity->Draw("SAME");

    tex->DrawLatexNDC(0.3,0.3,"simulations");

    c->Print(Form("../../plots/purity/npair_purity_rl_jetpt_deltarleq%.3f.pdf",R_L_res));
    
    // DATA PLOTS
    THStack* s_data = new THStack();
    TLegend* l_data = new TLegend();

    for(int jet_pt_bin = 0 ; jet_pt_bin < Nbin_jet_pt ; jet_pt_bin++)
    {
        hsig_data[jet_pt_bin]->Multiply(hpurity[jet_pt_bin]);

        s_data->Add(hsig_data[jet_pt_bin]);
        s_data->Add(hall_data[jet_pt_bin]);

        l_data->AddEntry(hsig_data[jet_pt_bin],Form("Data w/ Purity : %.1f<Jet P_{T}(GeV)<%.1f",jet_pt_binning[jet_pt_bin],jet_pt_binning[jet_pt_bin+1]),"lpf");
        l_data->AddEntry(hall_data[jet_pt_bin],Form("Data All : %.1f<Jet P_{T}(GeV)<%.1f",jet_pt_binning[jet_pt_bin],jet_pt_binning[jet_pt_bin+1])      ,"lpf");
    }
    
    s_data->Draw("NOSTACK");
    s_data->GetXaxis()->SetRangeUser(R_L_min,1);
    s_data->SetTitle(Form("#Delta R_{L}(truth-reco)<%.3f;R_{L};N_{pair}",R_L_res));
    l_data->Draw("SAME");

    gPad->SetLogx(1);
    gPad->SetLogy(1);

    tex->DrawLatexNDC(0.3,0.3,"simulations");

    c->Print(Form("../../plots/purity/npair_wpurity_rl_data_jetpt_deltarleq%.3f.pdf",R_L_res));
}