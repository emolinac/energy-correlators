#include "../include/analysis-constants.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_singletrackcorrections(double corr_rel_error_local = corr_rel_error)
{
    gStyle->SetOptStat(1110);
    // Open the necessary files
    TFile*   fcorr       = new TFile((output_folder+namef_ntuple_e2c_corr).c_str()); 
    TNtuple* ntuple_data = (TNtuple*) fcorr->Get((name_ntuple_data).c_str());
    
    TH2F* h = new TH2F("h","",100,0.2,1,100,0.2,1);
    
    // Project into the histograms
    ntuple_data->Project("h","purity:efficiency",Form("efficiency_relerror<%f&&purity_relerror<%f",corr_rel_error_local,corr_rel_error_local));
    // ntuple_data->Project("h","purity:efficiency","efficiency>0&&purity>0");
    
    TCanvas* c = new TCanvas("c","",800,600);
    c->Draw();

    h->Draw("col");

    TLatex* tex = new TLatex();
    tex->SetTextColorAlpha(16,0.3);
    tex->SetTextSize(0.1991525);
    tex->SetTextAngle(26.15998);
    tex->SetLineWidth(2);

    h->SetTitle(";efficiency;purity");
    h->Smooth();

    // TLegend* l_data = new TLegend();
    // l_data->AddEntry(hcorr_data,"Corr. Data","lpf");
    // l_data->AddEntry(hall_data ,"Data"      ,"lpf");
    // l_data->Draw("SAME");

    // tex->DrawLatexNDC(0.25,0.25,"LHCb Internal");

    c->Print(Form("../plots/eff_purity_distribution_relerrorleq%.2f.pdf",corr_rel_error_local));
}