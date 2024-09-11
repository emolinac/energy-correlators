#include "../include/analysis-constants.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_purity_rl()
{
    // Open the necessary files
    TFile* fdata   = new TFile((output_folder+namef_ntuple_e2c).c_str());
    TFile* fpurity = new TFile((output_folder+namef_ntuple_e2c_purity).c_str());

    // Get the corresponding Ntuples
    TNtuple* ntuple_data   = (TNtuple*) fdata->Get((name_ntuple_data).c_str());
    TNtuple* ntuple_purity = (TNtuple*) fpurity->Get((name_ntuple_purity).c_str());

    // Determine log binnning
    double binning[Nbin_R_L+1];
    determine_log10binning(Nbin_R_L, R_L_min, R_L_max, binning);

    // Define the necessary histograms to calculate purity
    TH1F* hsig    = new TH1F("hsig"   ,"",Nbin_R_L,binning);
    TH1F* hall    = new TH1F("hall"   ,"",Nbin_R_L,binning);
    TH1F* hpurity = new TH1F("hpurity","",Nbin_R_L,binning);
    hsig->Sumw2();
    hall->Sumw2();
    hpurity->Sumw2();

    set_histogram_style(hsig, kViolet, std_line_width, std_marker_style, std_marker_size);
    set_histogram_style(hall, kCyan  , std_line_width, std_marker_style, std_marker_size);

    // Define the necessary histograms to show data and corrected data
    TH1F* hsig_data = new TH1F("hsig_data","",Nbin_R_L,binning);
    TH1F* hall_data = new TH1F("hall_data","",Nbin_R_L,binning);
    hsig_data->Sumw2();
    hall_data->Sumw2();

    set_histogram_style(hsig_data, kViolet, std_line_width, std_marker_style, std_marker_size);
    set_histogram_style(hall_data, kCyan  , std_line_width, std_marker_style, std_marker_size);

    // Project into the histograms
    ntuple_purity->Project("hsig","R_L",e2c_signal_cut);
    ntuple_purity->Project("hall","R_L",e2c_cut);
    ntuple_data->Project("hsig_data","R_L",e2c_cut);
    ntuple_data->Project("hall_data","R_L",e2c_cut);

    TCanvas* c = new TCanvas("c","",800,600);
    c->Draw();

    // MCRECO PLOTS
    THStack* s = new THStack();
    s->Add(hsig);
    s->Add(hall);
    s->Draw("NOSTACK");

    s->SetTitle(";X_{L};E2C");

    gPad->SetLogx(1);
    gPad->SetLogy(1);

    TLegend* l = new TLegend();
    l->AddEntry(hsig,"MCReco Signal","lpf");
    l->AddEntry(hall,"MCReco All"   ,"lpf");
    l->Draw("SAME");

    c->Print("../plots/purity_xl_signalvsall.pdf");
    gPad->SetLogy(0);

    // PURITY PLOTS
    hpurity->Divide(hsig,hall,1,1,"B");
    hpurity->Scale(100.);

    set_histogram_style(hpurity, kViolet, std_line_width, std_marker_style, std_marker_size);
    
    hpurity->Draw();
    hpurity->SetTitle(";X_{L};Purity(%)");

    c->Print("../plots/purity_xl.pdf");

    // DATA PLOTS
    hpurity->Scale(1./100.);
    hsig_data->Multiply(hpurity);

    THStack* s_data = new THStack();
    s_data->Add(hsig_data);
    s_data->Add(hall_data);
    s_data->Draw("NOSTACK");

    s_data->SetTitle(";X_{L};E2C");

    gPad->SetLogx(1);
    gPad->SetLogy(1);

    TLegend* l_data = new TLegend();
    l_data->AddEntry(hsig_data,"Data w/ Purity","lpf");
    l_data->AddEntry(hall_data,"Data All"      ,"lpf");
    l_data->Draw("SAME");

    c->Print("../plots/purity_xl_data.pdf");

}