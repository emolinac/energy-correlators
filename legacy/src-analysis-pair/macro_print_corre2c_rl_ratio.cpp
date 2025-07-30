#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_corre2c_rl_ratio()
{
    // Open the necessary files
    TFile* fdata       = new TFile((output_folder + namef_ntuple_e2c).c_str());
    TFile* fefficiency = new TFile((output_folder + namef_ntuple_e2c_hadroncorrections).c_str());
    TFile* fpurity     = new TFile((output_folder + namef_ntuple_e2c_hadroncorrections).c_str());


    // Get the corresponding Ntuples
    TNtuple* ntuple_data            = (TNtuple*) fdata->Get((name_ntuple_data).c_str());
    TNtuple* ntuple_efficiency_reco = (TNtuple*) fefficiency->Get((name_ntuple_correction_reco).c_str());
    TNtuple* ntuple_efficiency_mc   = (TNtuple*) fefficiency->Get((name_ntuple_correction_mc).c_str());
    TNtuple* ntuple_purity          = (TNtuple*) fpurity->Get((name_ntuple_purity).c_str());
    TNtuple* ntuple_mc              = (TNtuple*) fdata->Get((name_ntuple_mc).c_str());

    // Determine log binnning
    double binning[nbin_rl_nominal+1];
    determine_log10binning(nbin_rl_nominal, rl_min, rl_max, binning);

    // Define the necessary histograms to calculate efficiency
    TH1F* h_mc        = new TH1F("h_mc"       ,"",nbin_rl_nominal,rl_min, rl_max);
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
    h_mc->Sumw2();
    
    // Define the necessary histograms to show data and corrected data
    TH1F* hcorr_data = new TH1F("hcorr_data","",nbin_rl_nominal,rl_min, rl_max);
    TH1F* hall_data  = new TH1F("hall_data" ,"",nbin_rl_nominal,rl_min, rl_max);
    hcorr_data->Sumw2();
    hall_data->Sumw2();

    // Project into the histograms
    ntuple_mc->Draw("R_L>>h_mc",e2c_cut,"goff");
    ntuple_efficiency_reco->Project("hsig_eff","R_L",pair_signal_cut);
    ntuple_efficiency_mc->Project("hall_eff","R_L",pair_cut);
    ntuple_purity->Project("hsig_pur","R_L",pair_signal_cut);
    ntuple_purity->Project("hall_pur","R_L",pair_cut);
    ntuple_data->Project("hcorr_data","R_L",e2c_cut);
    ntuple_data->Project("hall_data","R_L",e2c_cut);
    
    TCanvas* c = new TCanvas("c","",800,600);
    c->Draw();

    TLatex* tex = new TLatex();
    set_lhcb_watermark_properties(tex);

    set_histogram_style(hcorr_data, kViolet, std_line_width, std_marker_style, std_marker_size);
    set_histogram_style(h_mc      , kCyan  , std_line_width, std_marker_style, std_marker_size);

    // Corrections
    hefficiency->Divide(hsig_eff,hall_eff,1,1,"B");
    hpurity->Divide(hsig_pur,hall_pur,1,1,"B");
    
    // DATA PLOTS
    hcorr_data->Divide(hefficiency);
    hcorr_data->Multiply(hpurity);
    h_mc->SetTitle(";R_{L};Norm. E2C");

    hcorr_data->Scale(1./hcorr_data->Integral());
    h_mc->Scale(1./h_mc->Integral());

    TRatioPlot* rp = new TRatioPlot(h_mc,hcorr_data);
    rp->Draw();
    rp->GetUpperPad()->SetLogx(1);
    rp->GetUpperPad()->SetLogy(1);
    rp->GetLowerPad()->SetLogx(1);
    
    rp->GetLowerRefYaxis()->SetLabelSize(0.03);
    rp->GetUpperRefYaxis()->SetLabelSize(0.03);
    
    rp->GetLowerRefXaxis()->SetRangeUser(rl_min, 1);
    rp->GetUpperRefXaxis()->SetRangeUser(rl_min, 1);
    
    rp->GetLowerRefGraph()->SetMaximum(1.2);
    rp->GetLowerRefGraph()->SetMinimum(0.8);
    rp->GetLowYaxis()->SetNdivisions(505);

    rp->GetUpperPad()->cd();
    TLegend* l_data = new TLegend();
    l_data->AddEntry(hcorr_data,"Corr. Data","lpf");
    l_data->AddEntry(h_mc      ,"MC"        ,"lpf");
    l_data->Draw("SAME");
    c->Update();

    tex->DrawLatexNDC(0.25,0.25,"LHCb Internal");

    c->Print(Form("./plots/ratio_corr_e2c_deltarleq%.3f.pdf",rl_resolution));
}