#include "../include/analysis-constants.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_recondtrs()
{
    // Open the necessary files
    TFile* f = new TFile((output_folder+namef_ntuple_jet_purity).c_str());
    
    // Get the corresponding Ntuples
    TNtuple* ntuple = (TNtuple*) f->Get((name_ntuple_jetpurity).c_str());

    TH1F* h[Nbin_jet_pt];
    THStack* hs = new THStack();
    TLegend* l = new TLegend();
    for(int bin = 0 ; bin < Nbin_jet_pt ; bin++)
    {
        h[bin] = new TH1F(Form("h[%i]",bin),"",50,0,50);
        ntuple->Project(Form("h[%i]",bin),"jet_ndtr",pair_jetpt_cut[bin]);
        hs->Add(h[bin]);
        l->AddEntry(h[bin],Form("%.1f<Jet p_{T}<%.1f GeV  <Dtr> = %.0f",jet_pt_binning[bin],jet_pt_binning[bin+1],h[bin]->GetMean()),"lpf");
        set_histogram_style(h[bin], corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);
    }

    TCanvas* c = new TCanvas("c","",1920,1080);
    c->Draw();
    hs->Draw("NOSTACK");
    hs->SetTitle(";Reco Ndtrs;");
    l->Draw("SAME");
    c->Print("../plots/avge_ndtrs.pdf");
}