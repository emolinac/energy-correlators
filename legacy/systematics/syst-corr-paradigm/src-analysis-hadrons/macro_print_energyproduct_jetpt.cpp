#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_energyproduct_jet_pt()
{
    // Open the necessary files
    TFile*   fcorr       = new TFile((output_folder + namef_ntuple_e2c_corr).c_str()); 
    TNtuple* ntuple_data = (TNtuple*) fcorr->Get((name_ntuple_data).c_str());
    
    TH1F* hcorr_data[nbin_jet_pt]; 
    // TH1F* hall_data[nbin_jet_pt];  

    TCanvas* c = new TCanvas("c", "", 1920, 1080);
    c->Draw();

    TLatex* tex = new TLatex();
    set_lhcb_watermark_properties(tex);

    THStack* s_data = new THStack();
    TLegend* l_data = new TLegend();
    
    for (int bin = 1 ; bin <= nbin_jet_pt ; bin++)
    {
        hcorr_data[bin-1] = new TH1F(Form("hcorr_data[%i]",bin-1),"",100,0,3E4);
        set_histogram_style(hcorr_data[bin-1], corr_marker_color_jet_pt[bin-1], std_line_width, corr_marker_style_jet_pt[bin-1], std_marker_size+1);
        ntuple_data->Project(Form("hcorr_data[%i]",bin-1),"h1_e*h2_e",pair_jet_pt_cut[bin-1]);
        s_data->Add(hcorr_data[bin-1]);
        l_data->AddEntry(hcorr_data[bin-1],Form("%.1f<p_{T,jet}<%.1f (GeV)",jet_pt_binning[bin-1],jet_pt_binning[bin]),"lpf");
    }
    
    s_data->Draw("NOSTACK");
    s_data->SetTitle(";E_{i}*E_{j};");
    gPad->SetLogx(1);
    gPad->SetLogy(1);
    l_data->Draw("SAME");

    tex->DrawLatexNDC(0.25,0.25,"LHCb Internal");

    c->Print("./plots/energyproduct_jet_pt.pdf");
}