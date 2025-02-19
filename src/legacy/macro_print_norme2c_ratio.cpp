#include "../include/analysis-constants.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_norme2c_ratio()
{
    // Open ROOT file with ntuple
    TFile* f = new TFile((output_folder+namef_ntuple_e2c).c_str());

    // Get the ntuple
    TNtuple* ntuple_data = (TNtuple*) f->Get((name_ntuple_data).c_str());
    TNtuple* ntuple_mc   = (TNtuple*) f->Get((name_ntuple_mc).c_str());

    // Determine binning
    double binning[Nbin_R_L+1];
    determine_log10binning(Nbin_R_L, R_L_min, R_L_max, binning);

    TH1F* h_data   = new TH1F("h_data"  ,"",Nbin_R_L,R_L_min,R_L_max);
    TH1F* h_mc     = new TH1F("h_mc"    ,"",Nbin_R_L,R_L_min,R_L_max);
    h_data->Sumw2();
    h_mc->Sumw2();

    // Create Canvas and draw in it
    TCanvas* c = new TCanvas("","",800,600);
    c->Draw();

    ntuple_data->Draw("R_L>>h_data",e2c_cut,"goff");
    ntuple_mc->Draw("R_L>>h_mc",e2c_cut,"goff");
    
    h_data->Scale(1./h_data->Integral());
    h_mc->Scale(1./h_mc->Integral());

    set_histogram_style(h_data , kViolet+2 , std_line_width, std_marker_style, std_marker_size);
    set_histogram_style(h_mc   , kGreen+2  , std_line_width, std_marker_style, std_marker_size);

    // Create the ratio plot
    h_data->SetTitle(";R_{L};Norm E2C");
    TRatioPlot* rp = new TRatioPlot(h_data,h_mc);
    rp->Draw();
    rp->GetUpperPad()->SetLogx(1);
    rp->GetUpperPad()->SetLogy(1);
    rp->GetLowerPad()->SetLogx(1);
    rp->GetLowerPad()->SetLogy(1);
    
    c->Print("../plots/rationorme2c_rl.pdf");
}