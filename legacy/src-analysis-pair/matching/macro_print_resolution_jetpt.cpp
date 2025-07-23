#include "../../include/analysis-constants.h"
#include "../../include/analysis-binning.h"
#include "../../include/analysis-cuts.h"
#include "../../include/directories.h"
#include "../../include/names.h"
#include "../../include/utils-algorithms.h"
#include "../../include/utils-visual.h"

void macro_print_resolution_jet_pt()
{
    // Open the necessary files
    TFile* fpurity = new TFile((output_folder + namef_ntuple_jet_purity).c_str());
    gStyle->SetOptStat(1110);

    // Get the corresponding Ntuples
    TNtuple* ntuple = (TNtuple*) fpurity->Get((name_ntuple_jetpurity).c_str());

    // Define the necessary histograms to calculate purity
    TH1F* hres   = new TH1F("hres"   ,"",100,-100,100);
    TH1F* hratio = new TH1F("hratio" ,"",100,0,2);
    hres->Sumw2();
    hratio->Sumw2();
    
    // Project into the histograms
    ntuple->Project("hres"  ,"jet_pt-jet_pt_truth");
    ntuple->Project("hratio","jet_pt/jet_pt_truth");
    
    TCanvas* c = new TCanvas("c","",800,600);
    c->Draw();

    TLine* line1 = new TLine();
    TLine* line2 = new TLine();
    line1->SetLineColorAlpha(3,0.4);
    line1->SetLineStyle(9);
    line2->SetLineColorAlpha(3,0.4);
    line2->SetLineStyle(9);

    set_histogram_style(hres, kCyan, std_line_width, std_marker_style, std_marker_size);
    set_histogram_style(hres, kCyan, std_line_width, std_marker_style, std_marker_size);
    hres->Draw();
    hres->SetTitle(";#Delta p^{jet}_{t}(Reco,Truth)(GeV);");
    line1->DrawLine(hres->GetBinCenter(hres->GetMaximumBin())-get_hwhm(hres),0,hres->GetBinCenter(hres->GetMaximumBin())-get_hwhm(hres),hres->GetMaximum());
    line2->DrawLine(hres->GetBinCenter(hres->GetMaximumBin())+get_hwhm(hres),0,hres->GetBinCenter(hres->GetMaximumBin())+get_hwhm(hres),hres->GetMaximum());

    gPad->SetLogy(1);
    
    c->Print("./plots/resolution_jet_pt.pdf");

    hratio->Draw();
    hratio->SetTitle(";p^{jet}_{t}(Reco)/p^{jet}_{t}(Truth);");
    line1->DrawLine(hratio->GetBinCenter(hratio->GetMaximumBin())-get_hwhm(hratio),0,hratio->GetBinCenter(hratio->GetMaximumBin())-get_hwhm(hratio),hratio->GetMaximum());
    line2->DrawLine(hratio->GetBinCenter(hratio->GetMaximumBin())+get_hwhm(hratio),0,hratio->GetBinCenter(hratio->GetMaximumBin())+get_hwhm(hratio),hratio->GetMaximum());

    gPad->SetLogy(0);
    
    c->Print("./plots/resolution_jet_pt_ratio.pdf");

}