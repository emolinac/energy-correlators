#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_resolution_jetpt()
{
    // Open the necessary files
    TFile* fpurity = new TFile((output_folder+namef_ntuple_jet_purity).c_str());
    gStyle->SetOptStat(1110);

    // Get the corresponding Ntuples
    TNtuple* ntuple = (TNtuple*) fpurity->Get((name_ntuple_jetpurity).c_str());

    // Define the necessary histograms to calculate purity
    TH1F* hres    = new TH1F("hres"    ,"",100,-100,100);
    TH1F* hratio  = new TH1F("hratio"  ,"",100,0,2);
    TH1F* hdirres = new TH1F("hdirres" ,"",200,0,0.5);
    hres->Sumw2();
    hratio->Sumw2();
    hdirres->Sumw2();
    
    // Project into the histograms
    ntuple->Project("hres"   ,"jet_pt-jet_pt_truth");
    ntuple->Project("hratio" ,"jet_pt/jet_pt_truth");
    ntuple->Project("hdirres","deltaR_matchedjets");

    for (int bin = 1 ; bin <= hdirres->GetNbinsX() ; bin++)
    {
        double total = hdirres->Integral();
        double percentage = 100.*hdirres->Integral(1,bin)/total;

        std::cout<<percentage<<"\% at "<<hdirres->GetBinCenter(bin)<<std::endl;
    }

    for (int bin = 1 ; bin <= hratio->GetNbinsX() ; bin++)
    {
        double total = hratio->Integral();
        double percentage = 100.*hratio->Integral(bin,hratio->GetNbinsX()-bin+1)/total;

        std::cout<<percentage<<"\% between "<<hratio->GetBinCenter(bin)<<" and "<<hratio->GetBinCenter(hratio->GetNbinsX()-bin+1)<<std::endl;
    }

    TCanvas* c = new TCanvas("c","",800,600);
    c->Draw();

    TLine* line1 = new TLine();
    TLine* line2 = new TLine();
    line1->SetLineColorAlpha(3,0.4);
    line1->SetLineStyle(9);
    line2->SetLineColorAlpha(3,0.4);
    line2->SetLineStyle(9);

    set_histogram_style(hres, kCyan, std_line_width, std_marker_style, std_marker_size);
    set_histogram_style(hdirres, kCyan, std_line_width, std_marker_style, std_marker_size);
    hres->Draw();
    hres->SetTitle(";#Delta p^{jet}_{t}(Reco,Truth)(GeV);");
    // line1->DrawLine(hres->GetBinCenter(hres->GetMaximumBin())-get_hwhm(hres),0,hres->GetBinCenter(hres->GetMaximumBin())-get_hwhm(hres),hres->GetMaximum());
    // line2->DrawLine(hres->GetBinCenter(hres->GetMaximumBin())+get_hwhm(hres),0,hres->GetBinCenter(hres->GetMaximumBin())+get_hwhm(hres),hres->GetMaximum());

    gPad->SetLogy(1);
    
    c->Print("./plots/resolution_jetpt.pdf");

    gPad->SetLogy(0);
    hdirres->Draw();
    hdirres->SetTitle(";#Delta R(Reco,Truth)(GeV);");
    
    c->Print("./plots/dirresolution_jet.pdf");

    hratio->Draw();
    hratio->SetTitle(";p^{jet}_{t}(Reco)/p^{jet}_{t}(Truth);");
    // line1->DrawLine(hratio->GetBinCenter(hratio->GetMaximumBin())-get_hwhm(hratio),0,hratio->GetBinCenter(hratio->GetMaximumBin())-get_hwhm(hratio),hratio->GetMaximum());
    // line2->DrawLine(hratio->GetBinCenter(hratio->GetMaximumBin())+get_hwhm(hratio),0,hratio->GetBinCenter(hratio->GetMaximumBin())+get_hwhm(hratio),hratio->GetMaximum());
    
    c->Print("./plots/resolution_jetpt_ratio.pdf");

}