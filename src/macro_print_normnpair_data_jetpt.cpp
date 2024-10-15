#include "../include/analysis-constants.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_normnpair_data_jetpt(bool include_neutrals = 0)
{
    // Open ROOT file with ntuple
    TFile* f = new TFile((output_folder+namef_ntuple_e2c).c_str());

    // Get the ntuple
    TNtuple* ntuple = (TNtuple*) f->Get((name_ntuple_data).c_str());
    
    // Determine binning
    double binning[Nbin_R_L+1];
    determine_log10binning(Nbin_R_L, R_L_min, R_L_max, binning);

    TH1F* h[Nbin_jet_pt];
    THStack* s = new THStack();
    TLegend* l = new TLegend();
    
    for(int jet_pt_bin = 0 ; jet_pt_bin < Nbin_jet_pt ; jet_pt_bin++)
    {
        h[jet_pt_bin] = new TH1F(Form("h[%i]",jet_pt_bin),"",Nbin_R_L, binning);
        h[jet_pt_bin]->Sumw2();

        if(include_neutrals) ntuple->Draw(Form("R_L>>h[%i]",jet_pt_bin),pair_data_jetpt_cut[jet_pt_bin],"goff");
        else ntuple->Draw(Form("R_L>>h[%i]",jet_pt_bin),pair_data_jetpt_noneutrals_cut[jet_pt_bin],"goff");

        h[jet_pt_bin]->Scale(1./h[jet_pt_bin]->Integral());
        set_histogram_style(h[jet_pt_bin] , corr_marker_color_jet_pt[jet_pt_bin] , std_line_width, corr_marker_style_jet_pt[jet_pt_bin], std_marker_size);

        s->Add(h[jet_pt_bin]);
        l->AddEntry(h[jet_pt_bin],Form("%.1f<Jet p_{T}<%.1f GeV",jet_pt_binning[jet_pt_bin],jet_pt_binning[jet_pt_bin+1]),"lpf");    
    }

    // Create Canvas and draw in it
    TCanvas* c = new TCanvas("","",800,600);
    c->Draw();
    s->Draw("NOSTACK");
    s->SetTitle(";R_{L};Norm. N_{pair}");
    s->SetMaximum(1);

    gPad->SetLogx(1);
    gPad->SetLogy(1);

    l->Draw("SAME");

    if(include_neutrals) c->Print("../plots/normnpair_rl_data_jetpt.pdf");
    else c->Print("../plots/normnpair_noneutrals_rl_data_jetpt.pdf");
}