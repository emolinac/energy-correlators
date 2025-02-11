#include "../../include/analysis-constants.h"
#include "../../include/analysis-cuts.h"
#include "../../include/directories.h"
#include "../../include/names.h"
#include "../../include/utils-algorithms.h"
#include "../../include/utils-visual.h"

void macro_print_pairefficiency_rl()
{
    // Open the necessary files
    TFile* fdata   = new TFile((output_folder+namef_ntuple_e2c).c_str());
    TFile* fefficiency = new TFile((output_folder+namef_ntuple_e2c_efficiency).c_str());

    // Get the corresponding Ntuples
    TNtuple* ntuple_data            = (TNtuple*) fdata->Get((name_ntuple_data).c_str());
    TNtuple* ntuple_efficiency_reco = (TNtuple*) fefficiency->Get((name_ntuple_efficiency_reco).c_str());
    TNtuple* ntuple_efficiency_mc   = (TNtuple*) fefficiency->Get((name_ntuple_efficiency_mc).c_str());

    // Determine log binnning
    double binning[Nbin_R_L+1];
    determine_log10binning(Nbin_R_L, R_L_min, R_L_max, binning);

    // Define the necessary histograms to calculate efficiency
    TH1F* hsig        = new TH1F("hsig"       ,"",Nbin_R_L,R_L_min,R_L_max);
    TH1F* hall        = new TH1F("hall"       ,"",Nbin_R_L,R_L_min,R_L_max);
    TH1F* hefficiency = new TH1F("hefficiency","",Nbin_R_L,R_L_min,R_L_max);
    hsig->Sumw2();
    hall->Sumw2();
    
    set_histogram_style(hsig, kViolet, std_line_width, std_marker_style, std_marker_size);
    set_histogram_style(hall, kCyan  , std_line_width, std_marker_style, std_marker_size);

    // Define the necessary histograms to show data and corrected data
    TH1F* hcorr_data = new TH1F("hcorr_data","",Nbin_R_L,R_L_min,R_L_max);
    TH1F* hall_data  = new TH1F("hall_data" ,"",Nbin_R_L,R_L_min,R_L_max);
    hcorr_data->Sumw2();
    hall_data->Sumw2();

    set_histogram_style(hcorr_data, kViolet, std_line_width, std_marker_style, std_marker_size);
    set_histogram_style(hall_data, kCyan  , std_line_width, std_marker_style, std_marker_size);

    // Project into the histograms
    ntuple_efficiency_reco->Project("hsig","R_L",pair_signal_cut+"(h1_pt<7&&h2_pt<7)");
    ntuple_efficiency_mc->Project("hall","R_L",pair_cut+"(h1_pt<7&&h2_pt<7)");
    ntuple_data->Project("hcorr_data","R_L",pair_cut+"(h1_pt<7&&h2_pt<7)");
    ntuple_data->Project("hall_data","R_L",pair_cut+"(h1_pt<7&&h2_pt<7)");
    
    TLatex* tex = new TLatex();
    tex->SetTextColorAlpha(16,0.3);
    tex->SetTextSize(0.1991525);
    tex->SetTextAngle(26.15998);
    tex->SetLineWidth(2);
    
    TCanvas* c = new TCanvas("c","",800,600);
    c->Draw();

    // MCRECO PLOTS
    THStack* s = new THStack();
    s->Add(hsig);
    s->Add(hall);
    s->Draw("NOSTACK");
    s->GetXaxis()->SetRangeUser(R_L_min,1);

    s->SetTitle(Form("#Delta R_{L}(truth-reco)<%.3f;R_{L};N_{pair}",R_L_res));

    gPad->SetLogx(1);
    gPad->SetLogy(1);

    TLegend* l = new TLegend();
    l->AddEntry(hsig,"Reco","lpf");
    l->AddEntry(hall,"MC"  ,"lpf");
    l->Draw("SAME");

    tex->DrawLatexNDC(0.3,0.3,"simulations");

    c->Print(Form("../../plots/efficiency/npair_rl_recovsmc_deltarleq%.3f.pdf",R_L_res));
    
    gPad->SetLogy(0);

    // efficiency PLOTS
    hefficiency->Divide(hsig,hall,1,1,"B");
    set_histogram_style(hefficiency, kViolet, std_line_width, std_marker_style, std_marker_size);
    
    hefficiency->Draw();
    hefficiency->GetXaxis()->SetRangeUser(R_L_min,1);
    hefficiency->GetYaxis()->SetRangeUser(0,0.6);
    hefficiency->SetTitle(Form("#Delta R_{L}(truth-reco)<%.3f;R_{L};Pair efficiency",R_L_res));

    tex->DrawLatexNDC(0.3,0.3,"simulations");

    c->Print(Form("../../plots/efficiency/npair_efficiency_rl_deltarleq%.3f.pdf",R_L_res));
    
    // DATA PLOTS
    hcorr_data->Divide(hefficiency);

    THStack* s_data = new THStack();
    s_data->Add(hcorr_data);
    s_data->Add(hall_data);
    s_data->Draw("NOSTACK");
    s_data->GetXaxis()->SetRangeUser(R_L_min,1);

    s_data->SetTitle(Form("#Delta R_{L}(truth-reco)<%.3f;R_{L};N_{pair}",R_L_res));

    gPad->SetLogx(1);
    gPad->SetLogy(1);

    TLegend* l_data = new TLegend();
    l_data->AddEntry(hcorr_data,"Data w/ efficiency","lpf");
    l_data->AddEntry(hall_data,"Data All"           ,"lpf");
    l_data->Draw("SAME");

    tex->DrawLatexNDC(0.3,0.3,"simulations");
    
    c->Print(Form("../../plots/efficiency/npair_wefficiency_rl_data_deltarleq%.3f.pdf",R_L_res));
}