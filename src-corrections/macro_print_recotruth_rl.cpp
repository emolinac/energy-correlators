#include "../include/analysis-constants.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_recotruth_rl(bool include_neutrals = 0)
{
    // Open the necessary files
    TFile* fpurity = new TFile((output_folder+namef_ntuple_e2c_purity).c_str());

    // Get the corresponding Ntuples
    TNtuple* ntuple_dtrmatch = (TNtuple*) fpurity->Get((name_ntuple_purity).c_str());

    // Define the necessary histograms to calculate purity
    TH1F* hres = new TH1F("hres","",500,0,1);
    TH1F* hsel = new TH1F("hsel","",500,0,1);
    
    // Project into the histograms
    if(include_neutrals) ntuple_dtrmatch->Project("hres","deltaR_h1",pair_all_cut);
    else ntuple_dtrmatch->Project("hres","deltaR_h1",pair_all_noneutrals_cut);

    // Check the number of entries
    double selected;
    double selected_bin;
    for(int i = 1 ; i <= hres->GetNbinsX() ; i++)
    {
        if(hres->GetBinCenter(i)>R_match_max){selected = hres->Integral(1,i-1); selected_bin = i; break;}
    }

    std::cout<<"The selected bin is "<<selected_bin<<std::endl;

    if(include_neutrals) ntuple_dtrmatch->Project("hsel","deltaR_h1",Form(pair_all_cut+"deltaR_h1<%f",R_match_max));
    else ntuple_dtrmatch->Project("hsel","deltaR_h1",Form(pair_all_noneutrals_cut+"deltaR_h1<%f",R_match_max));

    set_histogram_style(hres, kCyan, std_line_width, std_marker_style, std_marker_size);
    set_histogram_style(hsel, kCyan, std_line_width, std_marker_style, std_marker_size);
    hsel->SetFillColor(kViolet);
    
    TCanvas* c = new TCanvas("c","",800,600);
    c->Draw();

    // MCRECO PLOTS
    THStack* hs = new THStack();
    hs->Add(hres);    
    hs->Add(hsel);
    hs->Draw("NOSTACK");

    hs->SetTitle(";R_{L}(h_{truth},h_{reco});");

    gPad->SetLogx(1);
    gPad->SetLogy(1);

    double percentage = 100.*selected/(hres->Integral());

    TLegend* l = new TLegend();
    l->AddEntry(hsel,Form("%.2f %% of total events",percentage),"f");
    l->Draw("SAME");

    if(include_neutrals) c->Print("../plots/RL_hrecotruth.pdf");
    else c->Print("../plots/RL_hrecotruth_noneutrals.pdf");
}