#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_corre2c_rl()
{
    // Open the necessary files
    TFile* fdata       = new TFile((output_folder+namef_ntuple_e2c).c_str());
    TFile* fefficiency = new TFile((output_folder+namef_ntuple_e2c_efficiency).c_str());
    TFile* fpurity     = new TFile((output_folder+namef_ntuple_e2c_purity).c_str());


    // Get the corresponding Ntuples
    TNtuple* ntuple_data            = (TNtuple*) fdata->Get((name_ntuple_data).c_str());
    TNtuple* ntuple_efficiency_reco = (TNtuple*) fefficiency->Get((name_ntuple_efficiency_reco).c_str());
    TNtuple* ntuple_efficiency_mc   = (TNtuple*) fefficiency->Get((name_ntuple_efficiency_mc).c_str());
    TNtuple* ntuple_purity          = (TNtuple*) fpurity->Get((name_ntuple_purity).c_str());

    // Determine log binnning
    double binning[Nbin_R_L+1];
    determine_log10binning(Nbin_R_L, R_L_min, R_L_max, binning);

    // Define the necessary histograms to calculate efficiency
    TH1F* hsig_eff    = new TH1F("hsig_eff"   ,"",Nbin_R_L,R_L_min,R_L_max);
    TH1F* hall_eff    = new TH1F("hall_eff"   ,"",Nbin_R_L,R_L_min,R_L_max);
    TH1F* hsig_pur    = new TH1F("hsig_pur"   ,"",Nbin_R_L,R_L_min,R_L_max);
    TH1F* hall_pur    = new TH1F("hall_pur"   ,"",Nbin_R_L,R_L_min,R_L_max);
    TH1F* hefficiency = new TH1F("hefficiency","",Nbin_R_L,R_L_min,R_L_max);
    TH1F* hpurity     = new TH1F("hpurity"    ,"",Nbin_R_L,R_L_min,R_L_max);
    hsig_eff->Sumw2();
    hall_eff->Sumw2();
    hsig_pur->Sumw2();
    hall_pur->Sumw2();
    
    // Define the necessary histograms to show data and corrected data
    TH1F* hcorr_data = new TH1F("hcorr_data","",Nbin_R_L,R_L_min,R_L_max);
    TH1F* hall_data  = new TH1F("hall_data" ,"",Nbin_R_L,R_L_min,R_L_max);
    hcorr_data->Sumw2();
    hall_data->Sumw2();

    set_histogram_style(hcorr_data, kViolet, std_line_width, std_marker_style, std_marker_size);
    set_histogram_style(hall_data, kCyan  , std_line_width, std_marker_style, std_marker_size);

    // Project into the histograms
    ntuple_efficiency_reco->Project("hsig_eff","R_L",pair_signal_cut);
    ntuple_efficiency_mc->Project("hall_eff","R_L",pair_cut);
    ntuple_purity->Project("hsig_pur","R_L",pair_signal_cut);
    ntuple_purity->Project("hall_pur","R_L",pair_cut);
    ntuple_data->Project("hcorr_data","R_L",e2c_cut);
    ntuple_data->Project("hall_data","R_L",e2c_cut);

    hcorr_data->Scale(1./hcorr_data->Integral());
    hall_data->Scale(1./hall_data->Integral());
    
    TCanvas* c = new TCanvas("c","",800,600);
    c->Draw();

    TLatex* tex = new TLatex();
    tex->SetTextColorAlpha(16,0.3);
    tex->SetTextSize(0.1991525);
    tex->SetTextAngle(26.15998);
    tex->SetLineWidth(2);

    // Corrections
    hefficiency->Divide(hsig_eff,hall_eff,1,1,"B");
    hpurity->Divide(hsig_pur,hall_pur,1,1,"B");
    
    // DATA PLOTS
    hcorr_data->Divide(hefficiency);
    hcorr_data->Multiply(hpurity);

    THStack* s_data = new THStack();
    s_data->Add(hcorr_data);
    s_data->Add(hall_data);
    s_data->Draw("NOSTACK");

    s_data->SetTitle(Form("#Delta R_{L}(truth-reco)<%.3f;R_{L};E2C",R_L_res));

    gPad->SetLogx(1);
    gPad->SetLogy(1);

    TLegend* l_data = new TLegend();
    l_data->AddEntry(hcorr_data,"Corr. Data","lpf");
    l_data->AddEntry(hall_data ,"Data"          ,"lpf");
    l_data->Draw("SAME");

    tex->DrawLatexNDC(0.25,0.25,"LHCb Internal");

    c->Print(Form("./plots/corr_e2c_deltarleq%.3f.pdf",R_L_res));
}