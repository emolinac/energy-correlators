#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_corrnpair_singletrackcorr()
{
    // Open the necessary files
    TFile*   fcorr       = new TFile((output_folder + namef_ntuple_e2c_corr).c_str()); 
    TNtuple* ntuple_data = (TNtuple*) fcorr->Get((name_ntuple_data).c_str());
    
    // Determine log binnning
    double binning[nbin_rl_nominal+1];
    determine_log10binning(nbin_rl_nominal, rl_min, rl_max, binning);

    // Define the necessary histograms to calculate efficiency
    TH1F* hsig_eff    = new TH1F("hsig_eff"   ,"",nbin_rl_nominal,rl_min, rl_max);
    TH1F* hall_eff    = new TH1F("hall_eff"   ,"",nbin_rl_nominal,rl_min, rl_max);
    TH1F* hsig_pur    = new TH1F("hsig_pur"   ,"",nbin_rl_nominal,rl_min, rl_max);
    TH1F* hall_pur    = new TH1F("hall_pur"   ,"",nbin_rl_nominal,rl_min, rl_max);
    TH1F* hefficiency = new TH1F("hefficiency","",nbin_rl_nominal,rl_min, rl_max);
    TH1F* hpurity     = new TH1F("hpurity"    ,"",nbin_rl_nominal,rl_min, rl_max);
    hsig_eff->Sumw2();
    hall_eff->Sumw2();
    hsig_pur->Sumw2();
    hall_pur->Sumw2();
    
    // Define the necessary histograms to show data and corrected data
    TH1F* hcorr_data = new TH1F("hcorr_data","",nbin_rl_nominal,rl_min, rl_max);
    TH1F* hall_data  = new TH1F("hall_data" ,"",nbin_rl_nominal,rl_min, rl_max);
    hcorr_data->Sumw2();
    hall_data->Sumw2();

    set_histogram_style(hcorr_data, kViolet, std_line_width, std_marker_style, std_marker_size);
    set_histogram_style(hall_data, kCyan  , std_line_width, std_marker_style, std_marker_size);

    // Project into the histograms
    ntuple_data->Project("hcorr_data","R_L",full_corr_singletrack);
    ntuple_data->Project("hall_data","R_L",pair_cut);

    hcorr_data->Scale(1./hall_data->Integral());
    hall_data->Scale(1./hall_data->Integral());
    
    TCanvas* c = new TCanvas("c","",800,600);
    c->Draw();

    TLatex* tex = new TLatex();
    tex->SetTextColorAlpha(16,0.4);
    tex->SetTextSize(0.1991525);
    tex->SetTextAngle(26.15998);
    tex->SetLineWidth(2);

    THStack* s_data = new THStack();
    s_data->Add(hcorr_data);
    s_data->Add(hall_data);
    s_data->Draw("NOSTACK");

    s_data->SetTitle(";R_{L};Npair");

    gPad->SetLogx(1);
    gPad->SetLogy(1);

    TLegend* l_data = new TLegend();
    l_data->AddEntry(hcorr_data,"Corr. Data","lpf");
    l_data->AddEntry(hall_data ,"Data"      ,"lpf");
    l_data->Draw("SAME");

    tex->DrawLatexNDC(0.25,0.25,"LHCb Internal");

    c->Print("./plots/corr_npair.pdf");
}