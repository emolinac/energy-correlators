#include "../../include/analysis-constants.h"
#include "../../include/analysis-binning.h"
#include "../../include/analysis-cuts.h"
#include "../../include/directories.h"
#include "../../include/names.h"
#include "../../include/utils-algorithms.h"
#include "../../include/utils-visual.h"

void macro_print_pairpurity_rl()
{
    // Open the necessary files
    TFile* fdata   = new TFile((output_folder + namef_ntuple_e2c).c_str());
    TFile* fpurity = new TFile((output_folder + namef_ntuple_e2c_pairpurity).c_str());

    // Get the corresponding Ntuples
    TNtuple* ntuple_data     = (TNtuple*) fdata->Get((name_ntuple_data).c_str());
    TNtuple* ntuple_dtrmatch = (TNtuple*) fpurity->Get((name_ntuple_purity).c_str());

    // Determine log binnning
    double binning[nbin_rl+1];
    determine_log10binning(nbin_rl, rl_min, rl_max, binning);

    // Define the necessary histograms to calculate purity
    TH1F* hsig    = new TH1F("hsig"   ,"",nbin_rl,rl_min, rl_max);
    TH1F* hall    = new TH1F("hall"   ,"",nbin_rl,rl_min, rl_max);
    TH1F* hpurity = new TH1F("hpurity","",nbin_rl,rl_min, rl_max);
    hsig->Sumw2();
    hall->Sumw2();
    hpurity->Sumw2();

    set_histogram_style(hsig, kViolet, std_line_width, std_marker_style, std_marker_size);
    set_histogram_style(hall, kCyan  , std_line_width, std_marker_style, std_marker_size);

    // Define the necessary histograms to show data and corrected data
    TH1F* hsig_data = new TH1F("hsig_data","",nbin_rl,rl_min, rl_max);
    TH1F* hall_data = new TH1F("hall_data","",nbin_rl,rl_min, rl_max);
    hsig_data->Sumw2();
    hall_data->Sumw2();

    set_histogram_style(hsig_data, kViolet, std_line_width, std_marker_style, std_marker_size);
    set_histogram_style(hall_data, kCyan  , std_line_width, std_marker_style, std_marker_size);

    // Project into the histograms
    ntuple_dtrmatch->Project("hsig","R_L",pair_signal_cut);
    ntuple_dtrmatch->Project("hall","R_L",pair_cut);
    ntuple_data->Project("hsig_data","R_L",pair_cut);
    ntuple_data->Project("hall_data","R_L",pair_cut);
    
    TCanvas* c = new TCanvas("c","",800,600);
    c->Draw();

    TLatex* tex = new TLatex();
    set_lhcb_watermark_properties(tex);

    // MCRECO PLOTS
    THStack* s = new THStack();
    s->Add(hsig);
    s->Add(hall);
    s->Draw("NOSTACK");
    s->GetXaxis()->SetRangeUser(rl_min,1);

    s->SetTitle(Form("#Delta R_{L}(truth-reco)<%.3f;R_{L};N_{pair}",rl_resolution));

    gPad->SetLogx(1);
    gPad->SetLogy(1);

    TLegend* l = new TLegend();
    l->AddEntry(hsig,"MCReco Signal","lpf");
    l->AddEntry(hall,"MCReco All"   ,"lpf");
    l->Draw("SAME");

    tex->DrawLatexNDC(0.3,0.3,"simulations");

    c->Print(Form("./plots/npair_rl_signalvsall_deltarleq%.3f.pdf",rl_resolution));
    
    gPad->SetLogy(0);

    // PURITY PLOTS
    hpurity->Divide(hsig,hall,1,1,"B");
    
    set_histogram_style(hpurity, kViolet, std_line_width, std_marker_style, std_marker_size);
    
    hpurity->Draw();
    hpurity->GetXaxis()->SetRangeUser(rl_min,1);
    hpurity->GetYaxis()->SetRangeUser(0,1);
    hpurity->SetTitle(Form("#Delta R_{L}(truth-reco)<%.3f;R_{L};Pair Purity",rl_resolution));

    tex->DrawLatexNDC(0.3,0.3,"simulations");

    c->Print(Form("./plots/npair_purity_rl_deltarleq%.3f.pdf",rl_resolution));
    
    // DATA PLOTS
    hsig_data->Multiply(hpurity);

    THStack* s_data = new THStack();
    s_data->Add(hsig_data);
    s_data->Add(hall_data);
    s_data->Draw("NOSTACK");
    s_data->GetXaxis()->SetRangeUser(rl_min,1);

    s_data->SetTitle(Form("#Delta R_{L}(truth-reco)<%.3f;R_{L};N_{pair}",rl_resolution));

    gPad->SetLogx(1);
    gPad->SetLogy(1);

    TLegend* l_data = new TLegend();
    l_data->AddEntry(hsig_data,"Data w/ Purity","lpf");
    l_data->AddEntry(hall_data,"Data All"      ,"lpf");
    l_data->Draw("SAME");

    tex->DrawLatexNDC(0.3,0.3,"simulations");
    
    c->Print(Form("./plots/npair_wpurity_rl_data_deltarleq%.3f.pdf",rl_resolution));
}