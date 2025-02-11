#include "../include/analysis-constants.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_energyproduct_weights()
{
    // Open the necessary files
    TFile*   fcorr       = new TFile((output_folder+namef_ntuple_e2c_corr).c_str()); 
    TNtuple* ntuple_data = (TNtuple*) fcorr->Get((name_ntuple_data).c_str());
    
    TH1F* hcorr_data[5]; 
    double weight_binning[5+1];
    determine_log10binning(5,0.00001,0.2,weight_binning);    
    // double energyproduct_binning[5+1];
    // determine_log10binning(5,0.00001,3E4,energyproduct_binning);    
    
    TCanvas* c = new TCanvas("c","",1920,1080);
    c->Draw();

    TLatex* tex = new TLatex();
    tex->SetTextColorAlpha(16,0.3);
    tex->SetTextSize(0.1991525);
    tex->SetTextAngle(26.15998);
    tex->SetLineWidth(2);

    THStack* s_data = new THStack();
    TLegend* l_data = new TLegend();
    
    for(int bin = 0 ; bin < 5 ; bin++)
    {
        hcorr_data[bin] = new TH1F(Form("hcorr_data[%i]",bin),"",100,0,3E4);
        // hcorr_data[bin] = new TH1F(Form("hcorr_data[%i]",bin),"",5,energyproduct_binning);
        set_histogram_style(hcorr_data[bin], corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);
        ntuple_data->Project(Form("hcorr_data[%i]",bin),"h1_e*h2_e",Form("weight>%f&&weight<%f",weight_binning[bin],weight_binning[bin+1]));
        s_data->Add(hcorr_data[bin]);
        l_data->AddEntry(hcorr_data[bin],Form("%.3f<weight<%.3f GeV",weight_binning[bin],weight_binning[bin+1]),"lpf");
    }
    
    s_data->Draw("NOSTACK");
    s_data->SetTitle(";E_{i}*E_{j};");
    gPad->SetLogy(1);
    l_data->Draw("SAME");

    tex->DrawLatexNDC(0.25,0.25,"LHCb Internal");

    c->Print("../plots/energyproduct_weight.pdf");
}