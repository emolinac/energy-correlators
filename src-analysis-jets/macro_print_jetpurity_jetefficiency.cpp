#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_jetpurity_jetefficiency()
{
    // Open the necessary files
    TFile* fpurity = new TFile((output_folder+namef_ntuple_jet_purity).c_str());
    TFile* fefficiency = new TFile((output_folder+namef_ntuple_jet_efficiency).c_str());

    // Get the corresponding Ntuples
    TNtuple* ntuple_purity = (TNtuple*) fpurity->Get((name_ntuple_jetpurity).c_str());
    TNtuple* ntuple        = (TNtuple*) fefficiency->Get((name_ntuple_jetefficiency).c_str());

    // Define the necessary histograms to calculate purity
    TH1F* hsig_purity = new TH1F("hsig_purity"   ,"",Nbin_jet_pt_unfolding,unfolding_jetpt_binning);
    TH1F* hall_purity = new TH1F("hall_purity"   ,"",Nbin_jet_pt_unfolding,unfolding_jetpt_binning);
    TH1F* hpurity     = new TH1F("hpurity","",Nbin_jet_pt_unfolding,unfolding_jetpt_binning);
    hsig_purity->Sumw2();
    hall_purity->Sumw2();
    hpurity->Sumw2();

    set_histogram_style(hpurity, 797 , std_line_width, std_marker_style, std_marker_size);

    // Define the necessary histograms to calculate efficiency
    TH1F* hsig        = new TH1F("hsig"   ,"",Nbin_jet_pt_unfolding,unfolding_jetpt_binning);
    TH1F* hall        = new TH1F("hall"   ,"",Nbin_jet_pt_unfolding,unfolding_jetpt_binning);
    TH1F* hefficiency = new TH1F("hefficiency","",Nbin_jet_pt_unfolding,unfolding_jetpt_binning);
    hsig->Sumw2();
    hall->Sumw2();
    hefficiency->Sumw2();

    set_histogram_style(hefficiency, 868 , std_line_width, std_marker_style, std_marker_size);
        
    // Project into the histograms
    ntuple_purity->Project("hsig_purity","jet_pt","jet_pt_truth!=-999");
    ntuple_purity->Project("hall_purity","jet_pt");
    ntuple->Project("hsig","jet_pt_truth","jet_pt!=-999");
    ntuple->Project("hall","jet_pt_truth");
    
    TCanvas* c = new TCanvas("c","",800,600);
    c->Draw();

    TLatex* tex = new TLatex();
    set_lhcb_watermark_properties(tex);

    hpurity->Divide(hsig_purity,hall_purity,1,1,"B");
    hefficiency->Divide(hsig,hall,1,1,"B");

    hpurity->GetYaxis()->SetRangeUser(0,1.2);

    THStack* hs = new THStack();
    hs->Add(hefficiency);
    hs->Add(hpurity);
    hs->SetTitle(";p^{jet}_{t}(GeV);");
    hs->Draw("NOSTACK");
    
    // tex->DrawLatexNDC(0.3,0.3,"simulations");
    TLegend* l = new TLegend();
    l->AddEntry(hpurity,"Purity","lpf");
    l->AddEntry(hefficiency,"Efficiency","lpf");
    l->Draw("SAME");

    c->Print(Form("./plots/jet_purity_efficiency_ibrahimcuts.pdf"));

    for (int bin = 1 ; bin <= hpurity->GetNbinsX() ; bin++)
    {
        double jet_total_correction = hpurity->GetBinContent(bin)/hefficiency->GetBinContent(bin);
        std::cout<<"Total correction between "<<unfolding_jetpt_binning[bin-1]<<" and "<<unfolding_jetpt_binning[bin]<<" is "<<jet_total_correction<<std::endl;
    }
}