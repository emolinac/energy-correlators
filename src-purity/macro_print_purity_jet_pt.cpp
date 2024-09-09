#include "../include/analysis-constants.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_purity_jet_pt()
{
    // Open the necessary files
    TFile* fdata   = new TFile((output_folder+namef_ntuple_e2c).c_str());
    TFile* fpurity = new TFile((output_folder+namef_ntuple_e2c_purity).c_str());

    // Get the corresponding Ntuples
    TNtuple* ntuple_data   = (TNtuple*) fdata->Get((name_ntuple_data).c_str());
    TNtuple* ntuple_purity = (TNtuple*) fpurity->Get((name_ntuple_purity).c_str());

    // Determine log binnning
    // Define the necessary histograms to calculate purity
    TH1F* hsig    = new TH1F("hsig"   ,"",Nbin_jet_pt,jet_pt_binning);
    TH1F* hall    = new TH1F("hall"   ,"",Nbin_jet_pt,jet_pt_binning);
    TH1F* hpurity = new TH1F("hpurity","",Nbin_jet_pt,jet_pt_binning);
    hsig->Sumw2();
    hall->Sumw2();
    hpurity->Sumw2();

    set_histogram_style(hsig, kViolet, std_line_width, std_marker_style, std_marker_size);
    set_histogram_style(hall, kCyan  , std_line_width, std_marker_style, std_marker_size);

    // Define the necessary histograms to show data and corrected data
    TH1F* hsig_data = new TH1F("hsig_data","",Nbin_jet_pt,jet_pt_binning);
    TH1F* hall_data = new TH1F("hall_data","",Nbin_jet_pt,jet_pt_binning);
    hsig_data->Sumw2();
    hall_data->Sumw2();

    set_histogram_style(hsig_data, kViolet, std_line_width, std_marker_style, std_marker_size);
    set_histogram_style(hall_data, kCyan  , std_line_width, std_marker_style, std_marker_size);

    // Project into the histograms
    ntuple_purity->Project("hsig","jet_pt",e2c_signal_cut);
    ntuple_purity->Project("hall","jet_pt",e2c_cut);
    ntuple_data->Project("hsig_data","jet_pt",e2c_cut);
    ntuple_data->Project("hall_data","jet_pt",e2c_cut);

    TCanvas* c = new TCanvas("c","",800,600);
    c->Draw();

    // MCRECO PLOTS
    gPad->SetLogy(1);
    THStack* s = new THStack();
    s->Add(hsig);
    s->Add(hall);
    s->Draw("NOSTACK");

    s->SetTitle(";Jet p_{T};E2C");

    TLegend* l = new TLegend();
    l->AddEntry(hsig,"MCReco Signal","lpf");
    l->AddEntry(hall,"MCReco All"   ,"lpf");
    l->Draw("SAME");

    c->Print("../plots/purity_jet_pt_signalvsall.pdf");

    // PURITY PLOTS
    gPad->SetLogy(0);
    hpurity->Divide(hsig,hall,1,1,"B");
    hpurity->Scale(100.);

    set_histogram_style(hpurity, kViolet, std_line_width, std_marker_style, std_marker_size);
    
    hpurity->Draw();
    hpurity->SetTitle(";Jet p_{T};Purity(%)");

    c->Print("../plots/purity_jet_pt.pdf");

    // DATA PLOTS
    gPad->SetLogy(1);
    hpurity->Scale(1./100.);
    hsig_data->Multiply(hpurity);

    THStack* s_data = new THStack();
    s_data->Add(hsig_data);
    s_data->Add(hall_data);
    s_data->Draw("NOSTACK");

    s_data->SetTitle(";Jet p_{T};E2C");

    TLegend* l_data = new TLegend();
    l_data->AddEntry(hsig_data,"Data w/ Purity","lpf");
    l_data->AddEntry(hall_data,"Data All"      ,"lpf");
    l_data->Draw("SAME");

    c->Print("../plots/purity_jet_pt_data.pdf");

}