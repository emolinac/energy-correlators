#include "../include/analysis-constants.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_e2ccorrections_singletrackcorr()
{
    // Open the necessary files
    TFile*   fcorr       = new TFile((output_folder+namef_ntuple_e2c_corr).c_str()); 
    TNtuple* ntuple_data = (TNtuple*) fcorr->Get((name_ntuple_data).c_str());
    
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
    TH1F* hpur_data  = new TH1F("hpur_data","",Nbin_R_L,R_L_min,R_L_max);
    TH1F* heff_data  = new TH1F("heff_data","",Nbin_R_L,R_L_min,R_L_max);
    TH1F* hall_data  = new TH1F("hall_data" ,"",Nbin_R_L,R_L_min,R_L_max);
    hpur_data->Sumw2();
    heff_data->Sumw2();
    hall_data->Sumw2();

    set_histogram_style(hpurity, kViolet, std_line_width, std_marker_style, std_marker_size);
    set_histogram_style(hefficiency, kCyan  , std_line_width, std_marker_style, std_marker_size);

    // Project into the histograms
    ntuple_data->Project("hpur_data","R_L",e2c_purity_corr_singletrack);
    ntuple_data->Project("heff_data","R_L",e2c_efficiency_corr_singletrack);
    ntuple_data->Project("hall_data","R_L",e2c_cut);

    hpurity->Divide(hpur_data,hall_data,1,1,"B");
    hefficiency->Divide(hall_data,heff_data,1,1,"B");

    TCanvas* c = new TCanvas("c","",800,600);
    c->Draw();

    TLatex* tex = new TLatex();
    tex->SetTextColorAlpha(16,0.3);
    tex->SetTextSize(0.0591525);
    // tex->SetTextAngle(26.15998);
    tex->SetLineWidth(2);

    hpurity->SetMaximum(1);
    hpurity->SetMinimum(0);
    hpurity->Draw();
    hpurity->SetTitle(Form("Purity rel. error < %.2f;R_{L};Purity Corr. on E2C",corr_rel_error));

    gPad->SetLogx(1);
    
    tex->DrawLatexNDC(0.2,0.85,"pair-by-pair corr.");

    c->Print(Form("./plots/corr_e2cpurity_singletrackcorr_relerrorleq%.2f_3dcorr.pdf",corr_rel_error));

    hefficiency->SetMaximum(1);
    hefficiency->SetMinimum(0);
    hefficiency->Draw();
    hefficiency->SetTitle(Form("Efficiency rel. error < %.2f;R_{L}; Efficiency Corr. on E2C",corr_rel_error));

    gPad->SetLogx(1);
    
    tex->DrawLatexNDC(0.2,0.85,"pair-by-pair corr.");

    c->Print(Form("./plots/corr_e2cefficiency_singletrackcorr_relerrorleq%.2f_3dcorr.pdf",corr_rel_error));
}