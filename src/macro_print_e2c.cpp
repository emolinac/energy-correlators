#include "../include/analysis-constants.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_e2c()
{
    // Open ROOT file with ntuple
    TFile* f = new TFile((output_folder+namef_ntuple_e2c).c_str());

    // Get the ntuple
    TNtuple* ntuple_data   = (TNtuple*) f->Get((name_ntuple_data).c_str());
    TNtuple* ntuple_mcreco = (TNtuple*) f->Get((name_ntuple_mcreco).c_str());
    TNtuple* ntuple_mc     = (TNtuple*) f->Get((name_ntuple_mc).c_str());

    // Determine binning
    double binning[Nbin_X_L+1];
    determine_log10binning(Nbin_X_L, X_L_min, X_L_max, binning);

    TH1F* h_data   = new TH1F("h_data"  ,"",Nbin_X_L, binning);
    TH1F* h_mcreco = new TH1F("h_mcreco","",Nbin_X_L, binning);
    TH1F* h_mc     = new TH1F("h_mc"    ,"",Nbin_X_L, binning);
    h_data->Sumw2();
    h_mcreco->Sumw2();
    h_mc->Sumw2();

    // Create Canvas and draw in it
    TCanvas* c = new TCanvas("","",800,600);
    c->Draw();
    ntuple_data->Draw("X_L>>h_data",e2c_cut,"goff");
    ntuple_mcreco->Draw("X_L>>h_mcreco",e2c_cut,"goff");
    ntuple_mc->Draw("X_L>>h_mc",e2c_mc_cut,"goff");

    set_histogram_style(h_data,   kViolet+2, std_line_width, std_marker_style, std_marker_size);
    set_histogram_style(h_mcreco, kCyan, std_line_width, std_marker_style, std_marker_size);
    set_histogram_style(h_mc, kGreen+2, std_line_width, std_marker_style, std_marker_size);

    THStack* s = new THStack();
    s->Add(h_data);
    s->Add(h_mcreco);
    s->Add(h_mc);
    s->Draw("NOSTACK");
    s->SetTitle(";X_{L};E2C");

    gPad->SetLogx(1);
    gPad->SetLogy(1);

    c->Print("../plots/E2C.pdf");
}