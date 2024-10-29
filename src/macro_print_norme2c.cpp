#include "../include/analysis-constants.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_norme2c(bool include_neutrals = 0)
{
    // Open ROOT file with ntuple
    TFile* f = new TFile((output_folder+namef_ntuple_e2c).c_str());

    // Get the ntuple
    TNtuple* ntuple_data   = (TNtuple*) f->Get((name_ntuple_data).c_str());
    TNtuple* ntuple_mcreco = (TNtuple*) f->Get((name_ntuple_mcreco).c_str());
    TNtuple* ntuple_mc     = (TNtuple*) f->Get((name_ntuple_mc).c_str());

    // Determine binning
    double binning[Nbin_R_L+1];
    determine_log10binning(Nbin_R_L, R_L_min, R_L_max, binning);

    TH1F* h_data   = new TH1F("h_data"  ,"",Nbin_R_L,R_L_min,R_L_max);
    TH1F* h_mcreco = new TH1F("h_mcreco","",Nbin_R_L,R_L_min,R_L_max);
    TH1F* h_mc     = new TH1F("h_mc"    ,"",Nbin_R_L,R_L_min,R_L_max);
    h_data->Sumw2();
    h_mcreco->Sumw2();
    h_mc->Sumw2();

    // Create Canvas and draw in it
    TCanvas* c = new TCanvas("","",800,600);
    c->Draw();

    if(include_neutrals)
    {
        ntuple_data->Draw("R_L>>h_data",e2c_data_cut,"goff");
        ntuple_mcreco->Draw("R_L>>h_mcreco",e2c_data_cut,"goff");
        ntuple_mc->Draw("R_L>>h_mc",e2c_mc_cut,"goff");
    }
    else
    {
        ntuple_data->Draw("R_L>>h_data",e2c_data_noneutrals_cut,"goff");
        ntuple_mcreco->Draw("R_L>>h_mcreco",e2c_data_noneutrals_cut,"goff");
        ntuple_mc->Draw("R_L>>h_mc",e2c_mc_noneutrals_cut,"goff");
    }
    
    h_data->Scale(1./h_data->Integral());
    h_mcreco->Scale(1./h_mcreco->Integral());
    h_mc->Scale(1./h_mc->Integral());

    set_histogram_style(h_data   , kViolet+2 , std_line_width, std_marker_style, std_marker_size);
    set_histogram_style(h_mcreco , kCyan     , std_line_width, std_marker_style, std_marker_size);
    set_histogram_style(h_mc     , kGreen+2  , std_line_width, std_marker_style, std_marker_size);

    THStack* s = new THStack();
    s->Add(h_data);
    s->Add(h_mcreco);
    s->Add(h_mc);
    s->Draw("NOSTACK");
    s->SetTitle(";R_{L};Norm. E2C");

    s->SetMaximum(1);

    gPad->SetLogx(1);
    gPad->SetLogy(1);

    TLegend* l = new TLegend();
    l->AddEntry(h_data  ,"Data"  ,"lpf");
    l->AddEntry(h_mc    ,"MC"    ,"lpf");
    l->AddEntry(h_mcreco,"MCReco","lpf");
    l->Draw("SAME");

    if(include_neutrals) c->Print("../plots/norme2c_rl.pdf");
    else c->Print("../plots/norme2c_rl_noneutrals.pdf");
    
}