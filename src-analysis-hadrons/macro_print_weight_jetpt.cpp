#include "../include/analysis-constants.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_weight_jetpt()
{
    // Open the necessary files
    TFile*   fcorr       = new TFile((output_folder+namef_ntuple_e2c_corr).c_str()); 
    TNtuple* ntuple_data = (TNtuple*) fcorr->Get((name_ntuple_data).c_str());
    
    TH1F* hcorr_data[Nbin_jet_pt]; 
    // TH1F* hall_data[Nbin_jet_pt];  

    TCanvas* c = new TCanvas("c","",1920,1080);
    c->Draw();

    TLatex* tex = new TLatex();
    tex->SetTextColorAlpha(16,0.3);
    tex->SetTextSize(0.1991525);
    tex->SetTextAngle(26.15998);
    tex->SetLineWidth(2);

    THStack* s_data = new THStack();
    TLegend* l_data = new TLegend();
    
    for(int bin = 1 ; bin <= Nbin_jet_pt ; bin++)
    {
        hcorr_data[bin-1] = new TH1F(Form("hcorr_data[%i]",bin-1),"",100,weight_min,weight_max);
        set_histogram_style(hcorr_data[bin-1], corr_marker_color_jet_pt[bin-1], std_line_width, corr_marker_style_jet_pt[bin-1], std_marker_size+1);
        ntuple_data->Project(Form("hcorr_data[%i]",bin-1),"weight",pair_jetpt_cut[bin-1]);
        s_data->Add(hcorr_data[bin-1]);
        l_data->AddEntry(hcorr_data[bin-1],Form("%.1f<p^{jet}_{t}<%.1f GeV",jet_pt_binning[bin-1],jet_pt_binning[bin]),"lpf");
    }
    
    s_data->Draw("NOSTACK");
    s_data->SetTitle(";weight;");
    gPad->SetLogx(1);
    gPad->SetLogy(1);
    l_data->Draw("SAME");

    tex->DrawLatexNDC(0.25,0.25,"LHCb Internal");

    c->Print("../plots/weight_jetpt.pdf");
}