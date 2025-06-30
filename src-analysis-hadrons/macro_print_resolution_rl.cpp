#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_resolution_rl()
{
    gStyle->SetOptStat(1110);

    // Open the necessary files
    TFile* fpurity = new TFile((output_folder+namef_ntuple_e2c_pairpurity).c_str());

    // Get the corresponding Ntuples
    TNtuple* ntuple_dtrmatch = (TNtuple*) fpurity->Get((name_ntuple_purity).c_str());

    // Define the necessary histograms to calculate purity
    // TH1F* hres = new TH1F("hres","",200,-.06,.06);
    TH1F* hres   = new TH1F("hres","",1000,-.2,.2);
    TH1F* hratio = new TH1F("hratio","",1000,0.6,1.4);
    hres->Sumw2();
    hratio->Sumw2();
    set_histogram_style(hres, kViolet, std_line_width, std_marker_style, std_marker_size);
    set_histogram_style(hratio, kViolet, std_line_width, std_marker_style, std_marker_size);
    
    // Project into the histograms
    ntuple_dtrmatch->Project("hres"  ,"R_L-R_L_truth",pair_cut+"R_L!=0&&R_L_truth!=-999");
    ntuple_dtrmatch->Project("hratio","R_L/R_L_truth",pair_cut+"R_L!=0&&R_L_truth!=-999");
    
    // Calculate event fractions 
    double total = hres->Integral();
    
    std::cout<<"There are "<<total<<" entries."<<std::endl;
    std::cout<<"For resolution:"<<std::endl;
    for(int bin = 1 ; bin <= hres->GetNbinsX()/2. ; bin++)
    {
        double integral = hres->Integral(bin,hres->GetNbinsX()-bin);
        if (integral/total < 0.95) {std::cout<<"Between "<<hres->GetBinCenter(bin)<<" and "<<hres->GetBinCenter(hres->GetNbinsX()-bin)<<" is the 95%% of the sample"<<std::endl; break;}
    }
    for(int bin = 1 ; bin <= hres->GetNbinsX()/2 ; bin++)
    {
        double integral = hres->Integral(bin,hres->GetNbinsX()-bin);
        if (integral/total < 0.68) {std::cout<<"Between "<<hres->GetBinCenter(bin)<<" and "<<hres->GetBinCenter(hres->GetNbinsX()-bin)<<" is the 68%% of the sample"<<std::endl; break;}
    }
    std::cout<<"For resolutionas a percentage:"<<std::endl;
    for(int bin = 1 ; bin <= hratio->GetNbinsX()/2. ; bin++)
    {
        double integral = hratio->Integral(bin,hratio->GetNbinsX()-bin);
        if (integral/total < 0.95) {std::cout<<"Between "<<hratio->GetBinCenter(bin)<<" and "<<hratio->GetBinCenter(hratio->GetNbinsX()-bin)<<" is the 95%% of the sample"<<std::endl; break;}
    }
    
    // Draw
    TCanvas* c = new TCanvas("c","",800,600);
    c->Draw();

    TLatex* tex = new TLatex();
    // tex->SetTextColorAlpha(16,0.3);
    // tex->SetTextSize(0.1991525);
    // tex->SetTextAngle(26.15998);
    tex->SetLineWidth(2);

    hres->Draw();
    hres->SetTitle(";#Delta R_{L}(Reco,Truth);");

    gPad->SetLogy(1);

    TLine* line1 = new TLine();
    TLine* line2 = new TLine();
    line1->SetLineColorAlpha(3,0.4);
    line1->SetLineStyle(9);
    line2->SetLineColorAlpha(3,0.4);
    line2->SetLineStyle(9);

    // line1->DrawLine(hres->GetBinCenter(hres->GetMaximumBin())-get_hwhm(hres),0,hres->GetBinCenter(hres->GetMaximumBin())-get_hwhm(hres),hres->GetMaximum());
    // line2->DrawLine(hres->GetBinCenter(hres->GetMaximumBin())+get_hwhm(hres),0,hres->GetBinCenter(hres->GetMaximumBin())+get_hwhm(hres),hres->GetMaximum());
    
    // tex->DrawLatexNDC(0.2,0.75,Form("HWHM = %.5f",get_hwhm(hres)));

    c->Print("./plots/resolution_rl.pdf");

    hratio->Draw();
    hratio->SetTitle(";R_{L}(Reco)/R_{L}(Truth);");

    // line1->DrawLine(hratio->GetBinCenter(hratio->GetMaximumBin())-get_hwhm(hratio),0,hratio->GetBinCenter(hratio->GetMaximumBin())-get_hwhm(hratio),hratio->GetMaximum());
    // line2->DrawLine(hratio->GetBinCenter(hratio->GetMaximumBin())+get_hwhm(hratio),0,hratio->GetBinCenter(hratio->GetMaximumBin())+get_hwhm(hratio),hratio->GetMaximum());

    // tex->DrawLatexNDC(0.2,0.75,Form("HWHM = %.5f",get_hwhm(hratio)));

    // gPad->SetLogy(0);
    c->Print("./plots/resolution_rl_ratio.pdf");
}